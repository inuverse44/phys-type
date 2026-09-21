//! Driving a two-field spectrum run and exporting its rows.
//!
//! The run has three stages. A background survey locates the end of inflation
//! and fixes the pivot `k_pivot = a_* H_*`. Each requested wavenumber is then
//! carried forward to its own `x_init = k/(aH)` crossing, initialized in the
//! Bunch-Davies vacuum there, and evolved together with the background. The
//! spectra are projected at one requested output time.
//!
//! Re-integrating the background for every wavenumber is deliberate. It keeps
//! the first implementation obviously correct; interpolating a stored
//! background is an optimization that can be added without changing the rows
//! produced here.

use std::io::{self, Write};

use crate::integrate::{
    AdaptiveSystem, ErrorTolerances, StepSizeLimits, adaptive_step, next_step_size, rk4_step,
};
use crate::spectra::{
    ProjectedSpectra, SpectrumError, SpectrumTolerances, field_power_spectrum, mode_covariance,
    project,
};
use crate::{
    ContravariantField, CovariantMultifieldSystem, CurvedFieldSpaceGeometry, KinematicBasisError,
    M0, Mass, ModeInitializationError, MultifieldState, MultifieldStateDelta, Time,
    TwiceDifferentiableMultifieldPotential, TwoFieldKinematics, TwoFieldPerturbationState,
    TwoFieldPerturbationSystem, TwoFieldPerturbationTolerance,
};

/// Everything the run needs beyond the physical model and its initial state.
#[derive(Clone, Copy, Debug)]
pub struct SpectrumRunSettings {
    /// Step size attempted first in every stage.
    pub initial_dt: Time,
    /// Smallest permitted step size.
    pub min_dt: Time,
    /// Largest permitted step size.
    pub max_dt: Time,
    /// Tolerances for the shared adaptive controller.
    ///
    /// The absolute part carries one tolerance per group of components, each
    /// in that group's own dimension. The background gets a floor in Planck
    /// units; the modes get none, because their norm is a modulus that does
    /// not vanish and any floor there would become their accuracy ceiling.
    pub tolerances: ErrorTolerances<TwoFieldPerturbationTolerance>,
    /// A safety limit on attempted steps in any one stage.
    pub max_steps: usize,
    /// Give up on inflation ending after this many e-folds.
    pub max_e_folds: f64,
    /// How many e-folds before its own horizon exit each mode is initialized.
    ///
    /// This is the readable form of the depth `x_init = k/(aH)`, which is its
    /// exponential. The trade-off is steep and regular: one more e-fold costs
    /// about `e` times as much work, because that many more sub-horizon
    /// oscillations must be resolved, and buys about `e^2` in accuracy,
    /// because the Bunch-Davies form is approached as `1/x_init^2`.
    pub e_folds_before_exit: f64,
    /// The pivot exits the horizon this many e-folds before inflation ends.
    pub pivot_e_folds_before_end: f64,
    /// Spectra are reported this many e-folds before inflation ends.
    pub output_e_folds_before_end: f64,
    /// Smallest wavenumber of the logarithmic grid, in units of `k_pivot`.
    pub k_min_over_pivot: f64,
    /// Largest wavenumber of the logarithmic grid, in units of `k_pivot`.
    pub k_max_over_pivot: f64,
    /// Number of grid points, which must be at least one.
    pub k_count: usize,
    /// One background row every this many accepted steps.
    pub background_sample_every: usize,
    /// One mode-evolution row every this many accepted steps.
    pub mode_sample_every: usize,
    /// Reality and range tolerances applied when projecting.
    pub spectrum_tolerances: SpectrumTolerances,
}

impl Default for SpectrumRunSettings {
    fn default() -> Self {
        Self {
            initial_dt: Time::new(100.0),
            min_dt: Time::new(1.0e-6),
            max_dt: Time::new(1.0e5),
            tolerances: ErrorTolerances::new(
                TwoFieldPerturbationTolerance::background_in_planck_units(1.0e-20),
                1.0e-10,
            ),
            max_steps: 1_000_000,
            max_e_folds: 120.0,
            e_folds_before_exit: 5.0,
            pivot_e_folds_before_end: 55.0,
            output_e_folds_before_end: 5.0,
            k_min_over_pivot: 0.1,
            k_max_over_pivot: 10.0,
            k_count: 9,
            background_sample_every: 20,
            mode_sample_every: 50,
            spectrum_tolerances: SpectrumTolerances::default(),
        }
    }
}

/// Why a spectrum run could not be completed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpectrumRunError {
    /// A setting was outside its permitted range.
    InvalidSettings(&'static str),
    /// `epsilon_H` never reached one within `max_e_folds`.
    InflationDidNotEnd {
        /// The e-folding number actually reached.
        e_folds: f64,
    },
    /// Inflation ended before the requested pivot or output time could exist.
    NotEnoughEFolds {
        /// E-folds between the initial state and the end of inflation.
        available: f64,
        /// E-folds the settings asked for before the end of inflation.
        required: f64,
    },
    /// A stage could not meet its tolerances at the minimum step size.
    MinimumStepReached,
    /// A stage exhausted its step budget.
    StepLimitReached,
    /// A mode could not be placed in the Bunch-Davies vacuum.
    ModeInitialization {
        /// The wavenumber in Planck units.
        k: f64,
        /// Why initialization failed.
        source: ModeInitializationError,
    },
    /// A projection could not produce physical spectra.
    Spectrum {
        /// The wavenumber in Planck units.
        k: f64,
        /// Why the projection failed.
        source: SpectrumError,
    },
    /// The adiabatic/entropy basis was undefined along the background.
    KinematicBasis(KinematicBasisError),
}

/// The pivot convention actually realized by a run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pivot {
    /// `k_pivot = a_* H_*`.
    pub k: Mass,
    /// `H_*`, the Hubble rate at the pivot crossing.
    pub hubble: crate::Hubble,
    /// `epsilon_H` at the pivot crossing.
    pub epsilon_h: f64,
    /// E-folds from the initial state to the end of inflation.
    pub e_folds_at_end: f64,
    /// E-folds from the initial state to the pivot crossing.
    pub e_folds_at_pivot: f64,
}

/// One background row, in Planck units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackgroundRecord {
    /// `N = ln(a / a_initial)`.
    pub e_folds: f64,
    /// `phi^I / M_pl`.
    pub phi: [f64; 2],
    /// `(dphi^I/dt) / M_pl^2`.
    pub pi: [f64; 2],
    /// `H / M_pl`.
    pub hubble: f64,
    /// `epsilon_H`.
    pub epsilon_h: f64,
    /// Adiabatic unit vector components.
    pub e_sigma: [f64; 2],
    /// Entropy unit vector components.
    pub e_s: [f64; 2],
    /// Dimensionless signed turn rate.
    pub eta_perp: f64,
}

impl BackgroundRecord {
    /// The CSV header matching [`BackgroundRecord::to_csv_row`].
    pub const CSV_HEADER: &'static str =
        "N,phi_1,phi_2,pi_1,pi_2,H,epsilon_h,e_sigma_1,e_sigma_2,e_s_1,e_s_2,eta_perp";

    /// Formats the row as comma-separated values.
    pub fn to_csv_row(&self) -> String {
        format!(
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            self.e_folds,
            self.phi[0],
            self.phi[1],
            self.pi[0],
            self.pi[1],
            self.hubble,
            self.epsilon_h,
            self.e_sigma[0],
            self.e_sigma[1],
            self.e_s[0],
            self.e_s[1],
            self.eta_perp
        )
    }
}

/// One spectrum row, in Planck units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectrumRecord {
    /// `k / k_pivot`.
    pub k_over_pivot: f64,
    /// `k / M_pl`.
    pub k: f64,
    /// `P^{11} / M_pl^2`.
    pub p_11: f64,
    /// `P^{22} / M_pl^2`.
    pub p_22: f64,
    /// `Re P^{12} / M_pl^2`.
    pub re_p_12: f64,
    /// `P_{Q_sigma} / M_pl^2`.
    pub p_q_sigma: f64,
    /// `P_{Q_s} / M_pl^2`.
    pub p_q_s: f64,
    /// `C_{Q_sigma Q_s} / M_pl^2`.
    pub c_q_sigma_s: f64,
    /// `P_R`.
    pub p_r: f64,
    /// `P_S`.
    pub p_s: f64,
    /// Signed `C_RS`.
    pub c_rs: f64,
    /// `r_C`.
    pub r_c: f64,
}

impl SpectrumRecord {
    /// The CSV header matching [`SpectrumRecord::to_csv_row`].
    pub const CSV_HEADER: &'static str =
        "k_over_pivot,k,P_11,P_22,Re_P_12,P_Qsigma,P_Qs,C_Qsigma_Qs,P_R,P_S,C_RS,r_C";

    /// Formats the row as comma-separated values.
    pub fn to_csv_row(&self) -> String {
        format!(
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            self.k_over_pivot,
            self.k,
            self.p_11,
            self.p_22,
            self.re_p_12,
            self.p_q_sigma,
            self.p_q_s,
            self.c_q_sigma_s,
            self.p_r,
            self.p_s,
            self.c_rs,
            self.r_c
        )
    }
}

/// One row of a single mode's evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModeEvolutionRecord {
    /// `k / k_pivot`.
    pub k_over_pivot: f64,
    /// `N = ln(a / a_initial)`.
    pub e_folds: f64,
    /// `k / (a H)`.
    pub k_over_ah: f64,
    /// `P_R`.
    pub p_r: f64,
    /// `P_S`.
    pub p_s: f64,
    /// Signed `C_RS`.
    pub c_rs: f64,
    /// Dimensionless signed turn rate.
    pub eta_perp: f64,
}

impl ModeEvolutionRecord {
    /// The CSV header matching [`ModeEvolutionRecord::to_csv_row`].
    pub const CSV_HEADER: &'static str = "k_over_pivot,N,k_over_aH,P_R,P_S,C_RS,eta_perp";

    /// Formats the row as comma-separated values.
    pub fn to_csv_row(&self) -> String {
        format!(
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            self.k_over_pivot,
            self.e_folds,
            self.k_over_ah,
            self.p_r,
            self.p_s,
            self.c_rs,
            self.eta_perp
        )
    }
}

/// The rows produced by one spectrum run.
#[derive(Clone, Debug)]
pub struct SpectrumRun {
    /// The realized pivot convention.
    pub pivot: Pivot,
    /// Background rows from the survey stage.
    pub background: Vec<BackgroundRecord>,
    /// One spectrum row per requested wavenumber.
    pub spectra: Vec<SpectrumRecord>,
    /// Sampled mode evolution for every wavenumber.
    pub mode_evolution: Vec<ModeEvolutionRecord>,
}

impl SpectrumRun {
    /// Writes the background rows as CSV.
    pub fn write_background_csv<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_rows(
            writer,
            BackgroundRecord::CSV_HEADER,
            &self.background,
            |row| row.to_csv_row(),
        )
    }

    /// Writes the spectrum rows as CSV.
    pub fn write_spectra_csv<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_rows(writer, SpectrumRecord::CSV_HEADER, &self.spectra, |row| {
            row.to_csv_row()
        })
    }

    /// Writes the mode-evolution rows as CSV.
    pub fn write_mode_evolution_csv<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_rows(
            writer,
            ModeEvolutionRecord::CSV_HEADER,
            &self.mode_evolution,
            |row| row.to_csv_row(),
        )
    }
}

fn write_rows<W: Write, R>(
    writer: &mut W,
    header: &str,
    rows: &[R],
    format: impl Fn(&R) -> String,
) -> io::Result<()> {
    writeln!(writer, "{header}")?;
    for row in rows {
        writeln!(writer, "{}", format(row))?;
    }
    Ok(())
}

/// Writes the metadata comment block that makes a run reproducible.
///
/// The revision is supplied by the caller, because the library does not read
/// the environment or shell out.
pub fn write_metadata<W: Write>(
    writer: &mut W,
    settings: &SpectrumRunSettings,
    pivot: &Pivot,
    revision: Option<&str>,
) -> io::Result<()> {
    writeln!(
        writer,
        "# generator: phys-type {}",
        env!("CARGO_PKG_VERSION")
    )?;
    if let Some(revision) = revision {
        writeln!(writer, "# revision: {revision}")?;
    }
    writeln!(
        writer,
        "# e_folds_before_exit: {:e}",
        settings.e_folds_before_exit
    )?;
    writeln!(writer, "# x_init: {:e}", start_depth(settings))?;
    // The absolute tolerances are per component group, so the metadata
    // records each of them rather than one number.
    let absolute = &settings.tolerances.absolute;
    writeln!(
        writer,
        "# absolute_tolerance_phi: {:e}",
        absolute.background.phi[0].into_value()
    )?;
    writeln!(
        writer,
        "# absolute_tolerance_pi: {:e}",
        absolute.background.pi[0].into_value()
    )?;
    writeln!(
        writer,
        "# absolute_tolerance_a: {:e}",
        absolute.background.a.into_value()
    )?;
    writeln!(
        writer,
        "# absolute_tolerance_mode_amplitude: {:e}",
        absolute.rescaled_amplitude.into_value()
    )?;
    writeln!(
        writer,
        "# absolute_tolerance_mode_velocity: {:e}",
        absolute.rescaled_velocity.into_value()
    )?;
    writeln!(
        writer,
        "# relative_tolerance: {:e}",
        settings.tolerances.relative
    )?;
    writeln!(
        writer,
        "# pivot: k_pivot = a H at N_end - {:e}",
        settings.pivot_e_folds_before_end
    )?;
    writeln!(writer, "# k_pivot_over_m_pl: {:e}", pivot.k.into_value())?;
    writeln!(
        writer,
        "# hubble_at_pivot_over_m_pl: {:e}",
        pivot.hubble.into_value()
    )?;
    writeln!(writer, "# epsilon_h_at_pivot: {:e}", pivot.epsilon_h)?;
    writeln!(writer, "# e_folds_at_end: {:e}", pivot.e_folds_at_end)?;
    writeln!(writer, "# e_folds_at_pivot: {:e}", pivot.e_folds_at_pivot)?;
    writeln!(
        writer,
        "# output_e_folds_before_end: {:e}",
        settings.output_e_folds_before_end
    )
}

/// Advances `system` until `event` stops being positive.
///
/// The crossing is refined inside the last accepted step, so the returned
/// state sits just on the positive side of the event. `on_accepted` observes
/// every accepted step, including the refined final one.
fn advance_until<S, Event, Observe>(
    system: &S,
    initial_state: S::State,
    initial_time: Time,
    settings: &SpectrumRunSettings,
    tolerances: &ErrorTolerances<S::AbsoluteTolerance>,
    mut event: Event,
    mut on_accepted: Observe,
) -> Result<(S::State, Time), SpectrumRunError>
where
    S: AdaptiveSystem<IndependentVariable = Time>,
    S::State: Clone,
    Event: FnMut(&S::State) -> f64,
    Observe: FnMut(&S::State, Time),
{
    let limits = StepSizeLimits {
        min: settings.min_dt,
        max: settings.max_dt,
    };
    let mut state = initial_state;
    let mut time = initial_time;
    let mut dt = settings.initial_dt;

    if event(&state) <= 0.0 {
        return Ok((state, time));
    }

    for _ in 0..settings.max_steps {
        let attempt = adaptive_step(system, time, &state, dt, tolerances);

        if attempt.accepted {
            if event(&attempt.state) <= 0.0 {
                let (refined, refined_dt) = refine_crossing(system, time, &state, dt, &mut event);
                let refined_time = time + refined_dt;
                on_accepted(&refined, refined_time);

                return Ok((refined, refined_time));
            }

            state = attempt.state;
            time = time + dt;
            on_accepted(&state, time);
        } else if dt <= settings.min_dt {
            return Err(SpectrumRunError::MinimumStepReached);
        }

        dt = next_step_size(dt, attempt.error_ratio, &limits);
    }

    Err(SpectrumRunError::StepLimitReached)
}

/// Bisects inside one accepted step for the largest sub-step that keeps the
/// event positive.
///
/// Taking the positive side means a state handed to Bunch-Davies
/// initialization is never a hair shallower than the requested depth.
fn refine_crossing<S, Event>(
    system: &S,
    time: Time,
    state: &S::State,
    dt: Time,
    event: &mut Event,
) -> (S::State, Time)
where
    S: AdaptiveSystem<IndependentVariable = Time>,
    S::State: Clone,
    Event: FnMut(&S::State) -> f64,
{
    let mut low = 0.0;
    let mut high = 1.0;
    let mut best = state.clone();
    let mut best_fraction = 0.0;

    for _ in 0..60 {
        let middle = 0.5 * (low + high);
        let trial = rk4_step(system, time, state, dt * middle);

        if event(&trial) > 0.0 {
            low = middle;
            best = trial;
            best_fraction = middle;
        } else {
            high = middle;
        }

        if high - low < 1.0e-14 {
            break;
        }
    }

    (best, dt * best_fraction)
}

impl<P, G> CovariantMultifieldSystem<P, G, 2>
where
    P: TwiceDifferentiableMultifieldPotential<2> + Clone,
    G: CurvedFieldSpaceGeometry<2> + Clone,
{
    /// Runs the background survey, fixes the pivot, and evolves one mode
    /// matrix per requested wavenumber.
    pub fn run_two_field_spectra(
        &self,
        initial: MultifieldState<2>,
        settings: &SpectrumRunSettings,
    ) -> Result<SpectrumRun, SpectrumRunError> {
        validate(settings)?;

        let (background, e_folds_at_end) = self.survey_background(initial, settings)?;
        let e_folds_at_pivot = e_folds_at_end - settings.pivot_e_folds_before_end;
        let e_folds_at_output = e_folds_at_end - settings.output_e_folds_before_end;
        if e_folds_at_pivot < 0.0 || e_folds_at_output < 0.0 {
            return Err(SpectrumRunError::NotEnoughEFolds {
                available: e_folds_at_end,
                required: settings
                    .pivot_e_folds_before_end
                    .max(settings.output_e_folds_before_end),
            });
        }

        let pivot_state = self.advance_to_e_folds(initial, settings, e_folds_at_pivot)?;
        let pivot_hubble = self.hubble(&pivot_state);
        let pivot = Pivot {
            k: pivot_state.a * pivot_hubble,
            hubble: pivot_hubble,
            epsilon_h: self.epsilon_h(&pivot_state),
            e_folds_at_end,
            e_folds_at_pivot,
        };

        let mut spectra = Vec::with_capacity(settings.k_count);
        let mut mode_evolution = Vec::new();
        for index in 0..settings.k_count {
            let k_over_pivot = grid_point(settings, index);
            let (record, rows) = self.evolve_one_mode(
                initial,
                settings,
                pivot.k * k_over_pivot,
                k_over_pivot,
                e_folds_at_output,
            )?;

            spectra.push(record);
            mode_evolution.extend(rows);
        }

        Ok(SpectrumRun {
            pivot,
            background,
            spectra,
            mode_evolution,
        })
    }

    /// Integrates the background to the end of inflation, recording rows.
    fn survey_background(
        &self,
        initial: MultifieldState<2>,
        settings: &SpectrumRunSettings,
    ) -> Result<(Vec<BackgroundRecord>, f64), SpectrumRunError> {
        let initial_kinematics = self
            .kinematics(&initial, None)
            .map_err(SpectrumRunError::KinematicBasis)?;
        let mut rows = vec![self.background_record(&initial, &initial, &initial_kinematics)];
        let mut previous_normal = Some(initial_kinematics.basis.normal);
        let mut basis_error = None;
        let mut accepted = 0_usize;

        // Inflation ends where epsilon_H reaches one. The survey also stops if
        // the e-folding budget runs out first, which is reported separately.
        let (final_state, _) = advance_until(
            self,
            initial,
            Time::new(0.0),
            settings,
            &background_tolerances(settings),
            |state| {
                let remaining = settings.max_e_folds - state.e_folds(&initial);
                let inflating = 1.0 - self.epsilon_h(state);

                remaining.min(inflating)
            },
            |state, _| {
                accepted += 1;
                if basis_error.is_some() {
                    return;
                }
                match self.kinematics(state, previous_normal.as_ref()) {
                    Ok(kinematics) => {
                        previous_normal = Some(kinematics.basis.normal);
                        if accepted.is_multiple_of(settings.background_sample_every) {
                            rows.push(self.background_record(&initial, state, &kinematics));
                        }
                    }
                    Err(error) => basis_error = Some(error),
                }
            },
        )?;

        if let Some(error) = basis_error {
            return Err(SpectrumRunError::KinematicBasis(error));
        }

        let e_folds = final_state.e_folds(&initial);
        if self.epsilon_h(&final_state) < 1.0 - 1.0e-9 {
            return Err(SpectrumRunError::InflationDidNotEnd { e_folds });
        }

        let final_kinematics = self
            .kinematics(&final_state, previous_normal.as_ref())
            .map_err(SpectrumRunError::KinematicBasis)?;
        rows.push(self.background_record(&initial, &final_state, &final_kinematics));

        Ok((rows, e_folds))
    }

    /// Integrates the background to a requested e-folding number.
    fn advance_to_e_folds(
        &self,
        initial: MultifieldState<2>,
        settings: &SpectrumRunSettings,
        target: f64,
    ) -> Result<MultifieldState<2>, SpectrumRunError> {
        let (state, _) = advance_until(
            self,
            initial,
            Time::new(0.0),
            settings,
            &background_tolerances(settings),
            |state| target - state.e_folds(&initial),
            |_, _| {},
        )?;

        Ok(state)
    }

    fn background_record(
        &self,
        initial: &MultifieldState<2>,
        state: &MultifieldState<2>,
        kinematics: &TwoFieldKinematics,
    ) -> BackgroundRecord {
        let squared = self.m_pl * self.m_pl;

        BackgroundRecord {
            e_folds: state.e_folds(initial),
            phi: core::array::from_fn(|index| (state.phi[index] / self.m_pl).into_value()),
            pi: core::array::from_fn(|index| (state.pi[index] / squared).into_value()),
            hubble: (self.hubble(state) / self.m_pl).into_value(),
            epsilon_h: self.epsilon_h(state),
            e_sigma: component_pair(&kinematics.basis.tangent),
            e_s: component_pair(&kinematics.basis.normal),
            eta_perp: kinematics.eta_perp.into_value(),
        }
    }

    /// Carries one wavenumber from its `x_init` crossing to the output time.
    fn evolve_one_mode(
        &self,
        initial: MultifieldState<2>,
        settings: &SpectrumRunSettings,
        k: Mass,
        k_over_pivot: f64,
        e_folds_at_output: f64,
    ) -> Result<(SpectrumRecord, Vec<ModeEvolutionRecord>), SpectrumRunError> {
        let system = TwoFieldPerturbationSystem::new(self.clone(), k);

        // Carry the background forward until this mode sits at x_init.
        let (start, start_time) = advance_until(
            self,
            initial,
            Time::new(0.0),
            settings,
            &background_tolerances(settings),
            |state| (k / (state.a * self.hubble(state))).into_value() - start_depth(settings),
            |_, _| {},
        )?;

        let mode_state = system
            .bunch_davies_state(start, start_depth(settings))
            .map_err(|source| SpectrumRunError::ModeInitialization {
                k: k.into_value(),
                source,
            })?;

        let (projected, kinematics) = self.project_state(&system, &mode_state, settings, None)?;
        let mut rows = vec![self.mode_record(
            &initial,
            &mode_state,
            k,
            k_over_pivot,
            &projected,
            &kinematics,
        )];
        let mut previous_normal = Some(kinematics.basis.normal);
        let mut sample_error = None;
        let mut accepted = 0_usize;

        let (final_state, _) = advance_until(
            &system,
            mode_state,
            start_time,
            settings,
            &settings.tolerances,
            |state| e_folds_at_output - state.background.e_folds(&initial),
            |state, _| {
                accepted += 1;
                if sample_error.is_some() || !accepted.is_multiple_of(settings.mode_sample_every) {
                    return;
                }
                match self.project_state(&system, state, settings, previous_normal.as_ref()) {
                    Ok((projected, kinematics)) => {
                        previous_normal = Some(kinematics.basis.normal);
                        rows.push(self.mode_record(
                            &initial,
                            state,
                            k,
                            k_over_pivot,
                            &projected,
                            &kinematics,
                        ));
                    }
                    Err(error) => sample_error = Some(error),
                }
            },
        )?;

        if let Some(error) = sample_error {
            return Err(error);
        }

        let (projected, kinematics) =
            self.project_state(&system, &final_state, settings, previous_normal.as_ref())?;
        rows.push(self.mode_record(
            &initial,
            &final_state,
            k,
            k_over_pivot,
            &projected,
            &kinematics,
        ));

        let matrix = field_power_spectrum(
            k,
            &mode_covariance(&final_state.modes, final_state.background.a),
        );
        let components = matrix.components();
        let squared = (self.m_pl * self.m_pl).into_value();

        let record = SpectrumRecord {
            k_over_pivot,
            k: (k / self.m_pl).into_value(),
            p_11: components[0][0].value().re / squared,
            p_22: components[1][1].value().re / squared,
            re_p_12: components[0][1].value().re / squared,
            p_q_sigma: projected.field_adiabatic.into_value() / squared,
            p_q_s: projected.field_entropy.into_value() / squared,
            c_q_sigma_s: projected.field_cross.into_value() / squared,
            p_r: projected.curvature.into_value(),
            p_s: projected.entropy.into_value(),
            c_rs: projected.cross.into_value(),
            r_c: projected.correlation,
        };

        Ok((record, rows))
    }

    fn mode_record(
        &self,
        initial: &MultifieldState<2>,
        state: &TwoFieldPerturbationState,
        k: Mass,
        k_over_pivot: f64,
        projected: &ProjectedSpectra,
        kinematics: &TwoFieldKinematics,
    ) -> ModeEvolutionRecord {
        let hubble = self.hubble(&state.background);

        ModeEvolutionRecord {
            k_over_pivot,
            e_folds: state.background.e_folds(initial),
            k_over_ah: (k / (state.background.a * hubble)).into_value(),
            p_r: projected.curvature.into_value(),
            p_s: projected.entropy.into_value(),
            c_rs: projected.cross.into_value(),
            eta_perp: kinematics.eta_perp.into_value(),
        }
    }

    /// Projects one evolved state, also returning the kinematics so that the
    /// entropy orientation carries into the next sample.
    fn project_state(
        &self,
        system: &TwoFieldPerturbationSystem<P, G>,
        state: &TwoFieldPerturbationState,
        settings: &SpectrumRunSettings,
        previous_normal: Option<&ContravariantField<M0, f64, 2>>,
    ) -> Result<(ProjectedSpectra, TwoFieldKinematics), SpectrumRunError> {
        let kinematics = self
            .kinematics(&state.background, previous_normal)
            .map_err(SpectrumRunError::KinematicBasis)?;
        let metric = self.geometry.metric(&state.background.phi);
        let hubble = self.hubble(&state.background);
        let matrix =
            field_power_spectrum(system.k, &mode_covariance(&state.modes, state.background.a));
        let projected = project(
            &matrix,
            &metric,
            &kinematics.basis,
            hubble,
            settings.spectrum_tolerances,
        )
        .map_err(|source| SpectrumRunError::Spectrum {
            k: system.k.into_value(),
            source,
        })?;

        Ok((projected, kinematics))
    }
}

fn component_pair(vector: &ContravariantField<M0, f64, 2>) -> [f64; 2] {
    let components = vector.components();

    [components[0].into_value(), components[1].into_value()]
}

/// The background-only tolerances implied by a run's settings.
fn background_tolerances(
    settings: &SpectrumRunSettings,
) -> ErrorTolerances<MultifieldStateDelta<2>> {
    ErrorTolerances::new(
        settings.tolerances.absolute.background,
        settings.tolerances.relative,
    )
}

/// The depth `x_init = k/(aH)` at which modes are initialized.
fn start_depth(settings: &SpectrumRunSettings) -> f64 {
    settings.e_folds_before_exit.exp()
}

/// The logarithmically spaced grid point at `index`.
fn grid_point(settings: &SpectrumRunSettings, index: usize) -> f64 {
    if settings.k_count == 1 {
        return settings.k_min_over_pivot;
    }

    let low = settings.k_min_over_pivot.ln();
    let high = settings.k_max_over_pivot.ln();
    let fraction = index as f64 / (settings.k_count - 1) as f64;

    (low + (high - low) * fraction).exp()
}

fn validate(settings: &SpectrumRunSettings) -> Result<(), SpectrumRunError> {
    if settings.k_count == 0 {
        return Err(SpectrumRunError::InvalidSettings(
            "k_count must be positive",
        ));
    }
    if settings.k_min_over_pivot.is_nan()
        || settings.k_max_over_pivot.is_nan()
        || settings.k_min_over_pivot <= 0.0
        || settings.k_max_over_pivot < settings.k_min_over_pivot
    {
        return Err(SpectrumRunError::InvalidSettings(
            "the k grid must be positive and ordered",
        ));
    }
    // Written without a negated comparison so that a NaN is rejected
    // explicitly rather than by falling through a flipped inequality.
    if settings.e_folds_before_exit.is_nan() || settings.e_folds_before_exit <= 0.0 {
        return Err(SpectrumRunError::InvalidSettings(
            "modes must start before their horizon exit",
        ));
    }
    if settings.min_dt.into_value() <= 0.0 || settings.max_dt < settings.min_dt {
        return Err(SpectrumRunError::InvalidSettings(
            "step size limits must be positive and ordered",
        ));
    }
    if settings.background_sample_every == 0 || settings.mode_sample_every == 0 {
        return Err(SpectrumRunError::InvalidSettings(
            "sampling intervals must be positive",
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;

    use super::*;
    use crate::{
        Field, FieldVelocity, FlatFieldSpace, QuadraticMultifieldPotential, ScaleFactor,
        multifield::MultifieldState,
    };

    type TestSystem = CovariantMultifieldSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace, 2>;

    fn system(masses: [f64; 2]) -> TestSystem {
        CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(masses[0]), Mass::new(masses[1])]),
            FlatFieldSpace,
        )
    }

    fn slow_roll_start(system: &TestSystem, phi: [f64; 2]) -> MultifieldState<2> {
        let resting = MultifieldState::new(
            [Field::new(phi[0]), Field::new(phi[1])],
            [FieldVelocity::new(0.0); 2],
            ScaleFactor::new(1.0),
        );

        MultifieldState {
            pi: system.slow_roll_velocity(&resting),
            ..resting
        }
    }

    /// A short run that still exercises every stage.
    fn quick_settings() -> SpectrumRunSettings {
        SpectrumRunSettings {
            e_folds_before_exit: 2.0,
            pivot_e_folds_before_end: 6.0,
            output_e_folds_before_end: 2.0,
            k_min_over_pivot: 1.0,
            k_max_over_pivot: 1.0,
            k_count: 1,
            background_sample_every: 50,
            mode_sample_every: 200,
            ..SpectrumRunSettings::default()
        }
    }

    #[test]
    fn the_grid_is_logarithmic_and_hits_both_ends() {
        let settings = SpectrumRunSettings {
            k_min_over_pivot: 0.1,
            k_max_over_pivot: 10.0,
            k_count: 3,
            ..SpectrumRunSettings::default()
        };

        assert!((grid_point(&settings, 0) - 0.1).abs() < 1.0e-12);
        assert!((grid_point(&settings, 1) - 1.0).abs() < 1.0e-12);
        assert!((grid_point(&settings, 2) - 10.0).abs() < 1.0e-12);
    }

    #[test]
    fn a_single_point_grid_uses_the_lower_end() {
        let settings = SpectrumRunSettings {
            k_count: 1,
            k_min_over_pivot: 2.0,
            ..SpectrumRunSettings::default()
        };

        assert_eq!(grid_point(&settings, 0), 2.0);
    }

    #[test]
    fn invalid_settings_are_rejected_before_any_integration() {
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let settings = SpectrumRunSettings {
            k_count: 0,
            ..SpectrumRunSettings::default()
        };

        assert!(matches!(
            system.run_two_field_spectra(initial, &settings),
            Err(SpectrumRunError::InvalidSettings(_))
        ));
    }

    #[test]
    fn a_run_without_enough_e_folds_is_reported() {
        let system = system([1.0e-6, 1.5e-6]);
        // Starting below the slow-roll attractor ends inflation almost at once.
        let initial = slow_roll_start(&system, [0.4, 0.3]);
        let settings = SpectrumRunSettings {
            pivot_e_folds_before_end: 55.0,
            ..quick_settings()
        };

        assert!(matches!(
            system.run_two_field_spectra(initial, &settings),
            Err(SpectrumRunError::NotEnoughEFolds { .. })
        ));
    }

    #[test]
    fn a_short_run_produces_rows_and_a_positive_pivot() {
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let run = system
            .run_two_field_spectra(initial, &quick_settings())
            .unwrap();

        assert!(run.pivot.k.into_value() > 0.0);
        assert!(run.pivot.e_folds_at_end > run.pivot.e_folds_at_pivot);
        assert!(!run.background.is_empty());
        assert_eq!(run.spectra.len(), 1);
        assert!(!run.mode_evolution.is_empty());

        let spectrum = run.spectra[0];
        assert!(spectrum.p_r > 0.0);
        assert!(spectrum.p_s > 0.0);
        assert!(spectrum.r_c.abs() <= 1.0);
        assert!(spectrum.p_11 > 0.0 && spectrum.p_22 > 0.0);
    }

    #[test]
    fn the_background_survey_ends_where_inflation_ends() {
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let run = system
            .run_two_field_spectra(initial, &quick_settings())
            .unwrap();

        let last = run.background.last().unwrap();
        assert!((last.epsilon_h - 1.0).abs() < 1.0e-6);
        assert!((last.e_folds - run.pivot.e_folds_at_end).abs() < 1.0e-9);
    }

    #[test]
    fn a_separable_straight_trajectory_has_no_turn_and_a_small_cross_spectrum() {
        // Equal masses make the trajectory a straight line in field space, so
        // the entropy direction never sources curvature.
        let system = system([1.0e-6, 1.0e-6]);
        let initial = slow_roll_start(&system, [10.0, 10.0]);
        let run = system
            .run_two_field_spectra(initial, &quick_settings())
            .unwrap();

        for row in &run.background {
            assert!(row.eta_perp.abs() < 1.0e-8, "eta_perp = {}", row.eta_perp);
        }

        let spectrum = run.spectra[0];
        let scale = spectrum.p_q_sigma.max(spectrum.p_q_s);
        assert!(
            spectrum.c_q_sigma_s.abs() < 1.0e-6 * scale,
            "cross = {}",
            spectrum.c_q_sigma_s
        );
    }

    #[test]
    fn the_single_field_limit_reproduces_the_analytic_curvature_spectrum() {
        // A trajectory confined to the first axis has no turn, so the
        // curvature spectrum must match the single-field slow-roll result
        // P_R = H_*^2 / (8 pi^2 epsilon_* M_pl^2) at horizon crossing.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 0.0]);
        let run = system
            .run_two_field_spectra(
                initial,
                &SpectrumRunSettings {
                    pivot_e_folds_before_end: 25.0,
                    ..quick_settings()
                },
            )
            .unwrap();

        let spectrum = run.spectra[0];
        assert!((spectrum.k_over_pivot - 1.0).abs() < 1.0e-12);

        let hubble = run.pivot.hubble.into_value();
        let expected = hubble * hubble / (8.0 * PI * PI * run.pivot.epsilon_h);
        let relative = (spectrum.p_r - expected).abs() / expected;

        // The formula is leading order in slow roll, so its own error is
        // O(epsilon_*). A fixed threshold would either be vacuous or fail
        // whenever the pivot moves.
        let allowed = 3.0 * run.pivot.epsilon_h;
        assert!(
            relative < allowed,
            "P_R = {}, expected {}, relative error {} exceeds {}",
            spectrum.p_r,
            expected,
            relative,
            allowed
        );
    }

    #[test]
    fn the_covariance_stays_positive_semidefinite() {
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let run = system
            .run_two_field_spectra(initial, &quick_settings())
            .unwrap();

        let spectrum = run.spectra[0];
        let determinant = spectrum.p_11 * spectrum.p_22 - spectrum.re_p_12 * spectrum.re_p_12;

        assert!(spectrum.p_11 > 0.0);
        assert!(spectrum.p_11 + spectrum.p_22 > 0.0);
        assert!(
            determinant > -1.0e-12 * spectrum.p_11 * spectrum.p_22,
            "determinant = {determinant}"
        );
    }

    #[test]
    fn the_spectrum_does_not_depend_on_where_the_modes_start() {
        // The Bunch-Davies form is approached as 1/x_init^2, so the two
        // depths the design recommends must agree far more closely than
        // either agrees with a shallow start.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);

        let shallow = system
            .run_two_field_spectra(
                initial,
                &SpectrumRunSettings {
                    e_folds_before_exit: 4.0,
                    ..quick_settings()
                },
            )
            .unwrap();
        let deep = system
            .run_two_field_spectra(
                initial,
                &SpectrumRunSettings {
                    e_folds_before_exit: 5.0,
                    ..quick_settings()
                },
            )
            .unwrap();

        let relative = (shallow.spectra[0].p_r - deep.spectra[0].p_r).abs() / deep.spectra[0].p_r;

        assert!(relative < 1.0e-3, "relative difference {relative}");
    }

    #[test]
    fn starting_deeper_inside_the_horizon_converges() {
        // The Bunch-Davies form used here is leading order, so its error
        // falls off as 1/x_init^2. Successive doublings must therefore move
        // the answer by ever less. A flat sequence would mean the mode
        // components are not actually being controlled to the requested
        // accuracy, which is exactly what a mis-scaled absolute tolerance
        // produces.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let spectrum = |e_folds_before_exit| {
            system
                .run_two_field_spectra(
                    initial,
                    &SpectrumRunSettings {
                        e_folds_before_exit,
                        ..quick_settings()
                    },
                )
                .unwrap()
                .spectra[0]
                .p_r
        };

        let coarse = spectrum(4.0);
        let medium = spectrum(5.0);
        let fine = spectrum(6.0);
        let finest = spectrum(7.0);

        let first = (medium - coarse).abs() / medium;
        let second = (fine - medium).abs() / fine;
        let third = (finest - fine).abs() / finest;

        assert!(second < 0.5 * first, "{first} then {second}");
        assert!(third < 0.5 * second, "{second} then {third}");
    }

    #[test]
    fn the_answer_is_insensitive_to_the_absolute_tolerance() {
        // The background floor and the mode components are now separate
        // settings, so moving the background floor over thirty orders of
        // magnitude cannot reach the modes at all.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let curvature = |absolute| {
            system
                .run_two_field_spectra(
                    initial,
                    &SpectrumRunSettings {
                        tolerances: ErrorTolerances::new(
                            TwoFieldPerturbationTolerance::background_in_planck_units(absolute),
                            1.0e-10,
                        ),
                        ..quick_settings()
                    },
                )
                .unwrap()
                .spectra[0]
                .p_r
        };

        let reference = curvature(1.0e-45);
        for absolute in [1.0e-30, 1.0e-20, 1.0e-16] {
            let value = curvature(absolute);
            let relative = (value - reference).abs() / reference;

            // Raising the floor across twenty-nine orders of magnitude moves
            // the answer by less than a part in ten million.
            assert!(
                relative < 1.0e-7,
                "absolute {absolute:e} moved P_R by {relative:e}"
            );
        }
    }

    #[test]
    fn the_default_leaves_the_modes_under_pure_relative_control() {
        // A mode error is normalized by a modulus, which does not vanish for
        // a physical mode, so no floor is needed to keep the ratio finite.
        // A floor here would only cap the accuracy of the modes.
        let absolute = SpectrumRunSettings::default().tolerances.absolute;

        assert_eq!(absolute.rescaled_amplitude.into_value(), 0.0);
        assert_eq!(absolute.rescaled_velocity.into_value(), 0.0);
        assert!(absolute.background.phi[0].into_value() > 0.0);
    }

    /// Advances a mode state to a requested e-folding number.
    fn advance_modes<P, G>(
        system: &TwoFieldPerturbationSystem<P, G>,
        settings: &SpectrumRunSettings,
        initial: &MultifieldState<2>,
        state: TwoFieldPerturbationState,
        time: Time,
        target: f64,
    ) -> (TwoFieldPerturbationState, Time)
    where
        P: TwiceDifferentiableMultifieldPotential<2> + Clone,
        G: CurvedFieldSpaceGeometry<2> + Clone,
    {
        advance_until(
            system,
            state,
            time,
            settings,
            &settings.tolerances,
            |state| target - state.background.e_folds(initial),
            |_, _| {},
        )
        .unwrap()
    }

    /// `[[R_1, R_2], [S_1, S_2]]` for one evolved state.
    fn projections<P, G>(
        system: &TwoFieldPerturbationSystem<P, G>,
        state: &TwoFieldPerturbationState,
    ) -> crate::ModeProjections
    where
        P: TwiceDifferentiableMultifieldPotential<2> + Clone,
        G: CurvedFieldSpaceGeometry<2> + Clone,
    {
        let kinematics = system
            .background
            .kinematics(&state.background, None)
            .unwrap();

        crate::project_modes(
            &state.modes,
            state.background.a,
            &system.background.geometry.metric(&state.background.phi),
            &kinematics.basis,
            system.background.hubble(&state.background),
            system.k,
        )
    }

    #[test]
    fn two_light_fields_freeze_at_the_de_sitter_amplitude() {
        // The `(H_*/2 pi)^2` result needs a window: earlier than about two
        // e-folds after exit the exact mode still carries its `1 + x^2`
        // correction, and later the entropy field decays under its own mass
        // while the adiabatic one grows with `2 epsilon`. Both effects are
        // measured in the next test; here the window is used.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let run = system
            .run_two_field_spectra(
                initial,
                &SpectrumRunSettings {
                    pivot_e_folds_before_end: 30.0,
                    output_e_folds_before_end: 28.0,
                    ..quick_settings()
                },
            )
            .unwrap();

        let hubble = run.pivot.hubble.into_value();
        let de_sitter = (hubble / (2.0 * PI)).powi(2);
        let spectrum = run.spectra[0];

        for (name, value) in [("P_Qsigma", spectrum.p_q_sigma), ("P_Qs", spectrum.p_q_s)] {
            let ratio = value / de_sitter;
            assert!(
                (ratio - 1.0).abs() < 0.15,
                "{name} / (H/2pi)^2 = {ratio}, heaviest m^2/H^2 = {}",
                (1.5e-6 / hubble).powi(2)
            );
        }
    }

    #[test]
    fn super_horizon_drift_follows_the_analytic_rates() {
        // Between two and eight e-folds after exit the entropy spectrum decays
        // as a light field of mass m does, and the adiabatic field spectrum
        // grows with 2 epsilon because Q_sigma = (dot sigma / H) R and R is
        // conserved.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let at = |after: f64| {
            system
                .run_two_field_spectra(
                    initial,
                    &SpectrumRunSettings {
                        pivot_e_folds_before_end: 30.0,
                        output_e_folds_before_end: 30.0 - after,
                        ..quick_settings()
                    },
                )
                .unwrap()
        };

        let early = at(2.0);
        let late = at(8.0);
        let elapsed = 6.0;

        // A light entropy field decays as P_S ~ a^(-2 m^2 / 3 H^2).
        let mass_ratio = (1.5e-6 / early.pivot.hubble.into_value()).powi(2);
        let predicted_decay = (-2.0 * mass_ratio * elapsed / 3.0).exp();
        let measured_decay = late.spectra[0].p_q_s / early.spectra[0].p_q_s;
        assert!(
            (measured_decay / predicted_decay - 1.0).abs() < 0.1,
            "entropy decayed by {measured_decay}, predicted {predicted_decay}"
        );

        // The adiabatic field spectrum tracks 2 epsilon, which for a quadratic
        // model is 1/N e-folds before the end.
        let predicted_growth = (30.0 - 2.0) / (30.0 - 8.0);
        let measured_growth = late.spectra[0].p_q_sigma / early.spectra[0].p_q_sigma;
        assert!(
            (measured_growth / predicted_growth - 1.0).abs() < 0.1,
            "adiabatic grew by {measured_growth}, predicted {predicted_growth}"
        );
    }

    #[test]
    fn the_super_horizon_transfer_matrix_is_triangular() {
        // Wands' statement: entropy can source curvature, but a purely
        // adiabatic perturbation creates no entropy. Two independent
        // solutions determine the transfer matrix as T = M_late M_early^-1,
        // whose entropy-from-curvature entry must vanish.
        let background = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&background, [12.0, 8.0]);
        let settings = quick_settings();
        let perturbation = TwoFieldPerturbationSystem::new(background, Mass::new(6.0e-5));

        let (start, start_time) = advance_until(
            &background,
            initial,
            Time::new(0.0),
            &settings,
            &background_tolerances(&settings),
            |state| {
                (perturbation.k / (state.a * background.hubble(state))).into_value()
                    - start_depth(&settings)
            },
            |_, _| {},
        )
        .unwrap();
        let exit = start.e_folds(&initial) + settings.e_folds_before_exit;
        let mode_state = perturbation
            .bunch_davies_state(start, start_depth(&settings))
            .unwrap();

        let (early, early_time) = advance_modes(
            &perturbation,
            &settings,
            &initial,
            mode_state,
            start_time,
            exit + 2.0,
        );
        let (late, _) = advance_modes(
            &perturbation,
            &settings,
            &initial,
            early,
            early_time,
            exit + 6.0,
        );

        let m_early = projections(&perturbation, &early);
        let m_late = projections(&perturbation, &late);
        let determinant = m_early[0][0] * m_early[1][1] - m_early[0][1] * m_early[1][0];
        let inverse = [
            [m_early[1][1] / determinant, -m_early[0][1] / determinant],
            [-m_early[1][0] / determinant, m_early[0][0] / determinant],
        ];
        let transfer: [[_; 2]; 2] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                m_late[row][0] * inverse[0][column] + m_late[row][1] * inverse[1][column]
            })
        });

        let entropy_from_curvature = transfer[1][0].norm();
        let entropy_from_entropy = transfer[1][1].norm();
        let curvature_from_entropy = transfer[0][1].norm();

        // The forbidden entry is not exactly zero, because the early time is
        // two e-folds outside the horizon rather than infinitely far, leaving
        // gradient corrections of order `(k/aH)^2 ~= 0.018`. It comes out far
        // below even that.
        assert!(
            entropy_from_curvature < 1.0e-3 * entropy_from_entropy,
            "T_SR = {entropy_from_curvature}, T_SS = {entropy_from_entropy}"
        );
        // The other off-diagonal entry is the physical sourcing, which a
        // turning trajectory must actually produce, so the triangle has a
        // direction rather than simply being small everywhere.
        assert!(
            curvature_from_entropy > 100.0 * entropy_from_curvature,
            "T_RS = {curvature_from_entropy}, T_SR = {entropy_from_curvature}"
        );
    }

    #[test]
    fn tighter_tolerances_do_not_move_the_answer() {
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let curvature = |relative| {
            system
                .run_two_field_spectra(
                    initial,
                    &SpectrumRunSettings {
                        tolerances: ErrorTolerances::new(
                            TwoFieldPerturbationTolerance::background_in_planck_units(1.0e-20),
                            relative,
                        ),
                        ..quick_settings()
                    },
                )
                .unwrap()
                .spectra[0]
                .p_r
        };

        let loose = curvature(1.0e-9);
        let tight = curvature(1.0e-12);

        assert!((loose - tight).abs() / tight < 1.0e-7);
    }

    #[test]
    fn a_denser_grid_does_not_disturb_the_wavenumbers_it_shares() {
        // Each wavenumber is integrated from the same initial state, so a
        // denser grid must reproduce the coarser one exactly. Anything else
        // would mean state leaking between grid points.
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let grid = |count| {
            system
                .run_two_field_spectra(
                    initial,
                    &SpectrumRunSettings {
                        k_min_over_pivot: 0.5,
                        k_max_over_pivot: 2.0,
                        k_count: count,
                        ..quick_settings()
                    },
                )
                .unwrap()
                .spectra
        };

        let coarse = grid(3);
        let dense = grid(5);

        assert_eq!(coarse[1].k_over_pivot, dense[2].k_over_pivot);
        assert_eq!(coarse[1].p_r, dense[2].p_r);
        assert_eq!(coarse[1].p_s, dense[2].p_s);
    }

    #[test]
    fn csv_rows_match_their_headers() {
        let system = system([1.0e-6, 1.5e-6]);
        let initial = slow_roll_start(&system, [12.0, 8.0]);
        let run = system
            .run_two_field_spectra(initial, &quick_settings())
            .unwrap();

        let mut buffer = Vec::new();
        run.write_spectra_csv(&mut buffer).unwrap();
        let text = String::from_utf8(buffer).unwrap();
        let mut lines = text.lines();

        let header = lines.next().unwrap();
        let columns = header.split(',').count();
        assert_eq!(header, SpectrumRecord::CSV_HEADER);
        for line in lines {
            assert_eq!(line.split(',').count(), columns);
        }
    }

    #[test]
    fn metadata_records_the_pivot_and_tolerances() {
        let settings = quick_settings();
        let pivot = Pivot {
            k: Mass::new(1.5e-5),
            hubble: crate::Hubble::new(2.0e-6),
            epsilon_h: 0.01,
            e_folds_at_end: 12.0,
            e_folds_at_pivot: 6.0,
        };

        let mut buffer = Vec::new();
        write_metadata(&mut buffer, &settings, &pivot, Some("abc123")).unwrap();
        let text = String::from_utf8(buffer).unwrap();

        assert!(text.contains("# revision: abc123"));
        assert!(text.contains("# e_folds_before_exit: 2e0"));
        assert!(text.contains("# k_pivot_over_m_pl: 1.5e-5"));
        assert!(text.lines().all(|line| line.starts_with('#')));
    }
}
