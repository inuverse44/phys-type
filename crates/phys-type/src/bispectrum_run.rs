//! A bispectrum driver with convention-labeled CSV output.
//!
//! [`CovariantMultifieldSystem::run_bispectrum`] surveys the background once,
//! fixes the pivot with the same convention as the spectrum run, evolves the
//! canonical tree-level three-point transport for each requested triangle,
//! and projects the result to `zeta` on a superhorizon output slice. The CSV
//! writers record the field-space correlation before projection, the
//! projected spectra and bispectrum, and metadata naming every convention
//! needed to interpret them.

use core::fmt;
use std::io::{self, Write};

use crate::curvature::{CurvatureProjectionError, curvature_correlators};
use crate::integrate::{ErrorTolerances, StepSizeLimits, adaptive_step, next_step_size};
use crate::spectrum_run::{BackgroundTrajectory, Pivot, SpectrumRunError, validate};
use crate::transport::{
    BackgroundNeighborhood, CanonicalThreePointTransportState, CanonicalThreePointTransportSystem,
    CanonicalThreePointTransportTolerance, MomentumTriangle, MomentumTriangleError,
    ThreePointInitialConditionError, canonical_three_point_transport_initial_state,
};
use crate::{
    CovariantMultifieldSystem, CurvatureDerivativeGeometry, Dimensionless, Mass, MultifieldState,
    SpectrumRunSettings, ThirdOrderFieldSpaceGeometry, ThriceDifferentiableMultifieldPotential,
};

/// A triangle in the shape coordinates of Fergusson and Shellard,
/// `k_1 = k_t (1 + alpha + beta) / 4`, `k_2 = k_t (1 - alpha + beta) / 4`,
/// `k_3 = k_t (1 - beta) / 2`, with `k_t` in units of the pivot wavenumber.
///
/// Equilateral is `(alpha, beta) = (0, 1/3)`; the vertices `(+-1, 0)` and
/// `(0, 1)` are squeezed limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriangleShape {
    /// Perimeter `k_t / k_pivot`.
    pub k_t_over_pivot: f64,
    /// Shape coordinate `alpha`.
    pub alpha: f64,
    /// Shape coordinate `beta`.
    pub beta: f64,
}

impl TriangleShape {
    /// The equilateral triangle with the given perimeter.
    pub const fn equilateral(k_t_over_pivot: f64) -> Self {
        Self {
            k_t_over_pivot,
            alpha: 0.0,
            beta: 1.0 / 3.0,
        }
    }

    /// `(k_1, k_2, k_3) / k_pivot` in the order of the shape coordinates.
    pub fn wavenumbers_over_pivot(&self) -> [f64; 3] {
        let k_t = self.k_t_over_pivot;
        [
            k_t * (1.0 + self.alpha + self.beta) / 4.0,
            k_t * (1.0 - self.alpha + self.beta) / 4.0,
            k_t * (1.0 - self.beta) / 2.0,
        ]
    }
}

/// Settings of a bispectrum run.
#[derive(Clone, Copy, Debug)]
pub struct BispectrumRunSettings {
    /// Background, pivot, and output-slice settings shared with the spectrum
    /// run. Its `e_folds_before_exit` is not used for the three-point start.
    pub background: SpectrumRunSettings,
    /// Each triangle starts when its shortest leg has
    /// `k/(aH) = exp(e_folds_before_exit)`.
    pub e_folds_before_exit: f64,
    /// Relative tolerance of the three-point transport.
    pub relative_tolerance: f64,
    /// The output slice must have `k/(aH)` at most this for every leg.
    pub maximum_output_k_over_ah: f64,
    /// Safety limit on attempted transport steps per triangle.
    pub max_steps: usize,
}

impl Default for BispectrumRunSettings {
    fn default() -> Self {
        Self {
            background: SpectrumRunSettings::default(),
            e_folds_before_exit: 5.0,
            relative_tolerance: 1.0e-9,
            maximum_output_k_over_ah: 0.1,
            max_steps: 5_000_000,
        }
    }
}

/// Why a bispectrum run could not be completed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BispectrumRunError {
    /// The background survey or settings failed.
    Spectrum(SpectrumRunError),
    /// The three-point settings are not usable.
    InvalidSettings(&'static str),
    /// A shape does not form a triangle.
    Triangle {
        /// Index of the shape in the request.
        shape: usize,
        /// Why the triangle was rejected.
        source: MomentumTriangleError,
    },
    /// The start slice of a triangle lies outside the surveyed trajectory.
    StartNotFound {
        /// Index of the shape in the request.
        shape: usize,
    },
    /// The subhorizon initial condition could not be formed.
    InitialCondition {
        /// Index of the shape in the request.
        shape: usize,
        /// The underlying failure.
        source: ThreePointInitialConditionError,
    },
    /// The adaptive transport exceeded `max_steps`.
    StepLimit {
        /// Index of the shape in the request.
        shape: usize,
    },
    /// The projection to `zeta` failed on the output slice.
    Projection {
        /// Index of the shape in the request.
        shape: usize,
        /// The underlying failure.
        source: CurvatureProjectionError,
    },
}

impl fmt::Display for BispectrumRunError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Spectrum(error) => write!(formatter, "{error}"),
            Self::InvalidSettings(reason) => write!(formatter, "invalid settings: {reason}"),
            Self::Triangle { shape, source } => write!(formatter, "shape {shape}: {source}"),
            Self::StartNotFound { shape } => write!(
                formatter,
                "shape {shape}: its start slice lies outside the surveyed background"
            ),
            Self::InitialCondition { shape, source } => {
                write!(formatter, "shape {shape}: {source}")
            }
            Self::StepLimit { shape } => {
                write!(
                    formatter,
                    "shape {shape}: the transport exceeded its step limit"
                )
            }
            Self::Projection { shape, source } => write!(formatter, "shape {shape}: {source}"),
        }
    }
}

impl std::error::Error for BispectrumRunError {}

/// One projected triangle on the output slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BispectrumRecord {
    /// The requested shape.
    pub shape: TriangleShape,
    /// Canonically ordered `k_1 >= k_2 >= k_3`, in units of `k_pivot`.
    pub k_over_pivot: [f64; 3],
    /// Elapsed e-folds of the output slice.
    pub e_folds: f64,
    /// Largest `k/(aH)` on the output slice.
    pub max_k_over_ah: f64,
    /// Elapsed e-folds of the common start slice.
    pub start_e_folds: f64,
    /// `k^3 P_zeta / (2 pi^2)` of each canonical leg.
    pub dimensionless_power: [f64; 3],
    /// `(k_1 k_2 k_3)^2 B_zeta`, independent of the scale-factor normalization.
    pub dimensionless_bispectrum: f64,
    /// Reduced bispectrum in the Maldacena/Planck convention.
    pub f_nl: f64,
    /// Fraction of `B_zeta` from `N_a N_b N_c Alpha^abc`.
    pub intrinsic_fraction: f64,
    /// `|Im B| / |Re B|` of the operator-ordered contraction.
    pub imaginary_ratio: f64,
}

impl BispectrumRecord {
    /// Column names matching [`Self::to_csv_row`].
    pub const CSV_HEADER: &'static str = "k_t_over_pivot,alpha,beta,k1_over_pivot,k2_over_pivot,k3_over_pivot,N_start,N,max_k_over_aH,Pcal_zeta_k1,Pcal_zeta_k2,Pcal_zeta_k3,Bcal_zeta,f_NL,intrinsic_fraction,imaginary_ratio";

    /// One CSV row in full precision.
    pub fn to_csv_row(&self) -> String {
        let [k1, k2, k3] = self.k_over_pivot;
        let [p1, p2, p3] = self.dimensionless_power;
        format!(
            "{:e},{:e},{:e},{k1:e},{k2:e},{k3:e},{:e},{:e},{:e},{p1:e},{p2:e},{p3:e},{:e},{:e},{:e},{:e}",
            self.shape.k_t_over_pivot,
            self.shape.alpha,
            self.shape.beta,
            self.start_e_folds,
            self.e_folds,
            self.max_k_over_ah,
            self.dimensionless_bispectrum,
            self.f_nl,
            self.intrinsic_fraction,
            self.imaginary_ratio,
        )
    }
}

/// One component of the canonical field-space three-point function before
/// projection, on the output slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseSpaceBispectrumRecord {
    /// Index of the shape in the request.
    pub shape: usize,
    /// Block name, `QQQ` through `PPP`, where `P` is the canonical momentum.
    pub block: &'static str,
    /// Field indices on the three canonical legs.
    pub indices: [usize; 3],
    /// `(k_1 k_2 k_3)^2 Alpha` in reduced Planck units with comoving
    /// wavenumbers, real part.
    pub real: f64,
    /// Imaginary part of the same.
    pub imaginary: f64,
}

impl PhaseSpaceBispectrumRecord {
    /// Column names matching [`Self::to_csv_row`].
    pub const CSV_HEADER: &'static str = "shape,block,i,j,k,Re_k123sq_alpha,Im_k123sq_alpha";

    /// One CSV row in full precision.
    pub fn to_csv_row(&self) -> String {
        let [i, j, k] = self.indices;
        format!(
            "{},{},{i},{j},{k},{:e},{:e}",
            self.shape, self.block, self.real, self.imaginary
        )
    }
}

/// The result of a bispectrum run.
#[derive(Clone, Debug, PartialEq)]
pub struct BispectrumRun {
    /// The pivot convention realized by the background survey.
    pub pivot: Pivot,
    /// Elapsed e-folds of the output slice.
    pub e_folds_at_output: f64,
    /// One projected record per requested shape.
    pub records: Vec<BispectrumRecord>,
    /// Every canonical component of `Alpha` before projection.
    pub phase_space: Vec<PhaseSpaceBispectrumRecord>,
}

impl BispectrumRun {
    /// Writes the projected records under their header.
    pub fn write_csv<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "{}", BispectrumRecord::CSV_HEADER)?;
        for record in &self.records {
            writeln!(writer, "{}", record.to_csv_row())?;
        }
        Ok(())
    }

    /// Writes the phase-space components under their header.
    pub fn write_phase_space_csv<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "{}", PhaseSpaceBispectrumRecord::CSV_HEADER)?;
        for record in &self.phase_space {
            writeln!(writer, "{}", record.to_csv_row())?;
        }
        Ok(())
    }
}

/// Writes `#`-prefixed metadata naming every convention of a bispectrum run.
///
/// `model` should identify the potential, geometry, parameters, and initial
/// conditions; it is written verbatim.
pub fn write_bispectrum_metadata<W: Write>(
    writer: &mut W,
    settings: &BispectrumRunSettings,
    run: &BispectrumRun,
    model: &str,
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
    writeln!(writer, "# model: {model}")?;
    writeln!(
        writer,
        "# method: tree-level canonical (Q, Pi) phase-space transport; docs/three-point-transport.md"
    )?;
    writeln!(
        writer,
        "# initial_condition: subhorizon Bunch-Davies in-in with first-order vertex rates; start when the shortest leg has k/(aH) = exp({:e})",
        settings.e_folds_before_exit
    )?;
    writeln!(
        writer,
        "# relative_tolerance: {:e}",
        settings.relative_tolerance
    )?;
    writeln!(
        writer,
        "# gauge: uniform-density zeta from flat-gauge fields, simple representation of arXiv:1410.3491 Eqs. (2.34a)-(2.34b); cross coefficients 1/6 (docs/literature-discrepancies.md, entry 1)"
    )?;
    writeln!(
        writer,
        "# zeta_sign: zeta = -phi'_I Q^I / (2 H M_pl^2 epsilon) = -R at linear order (Maldacena)"
    )?;
    writeln!(
        writer,
        "# normalization: <zeta zeta> = (2 pi)^3 delta P, <zeta zeta zeta> = (2 pi)^3 delta B; Pcal = k^3 P / (2 pi^2); Bcal = (k1 k2 k3)^2 B"
    )?;
    writeln!(
        writer,
        "# f_nl_convention: B = (6/5) f_NL [P(k1) P(k2) + P(k1) P(k3) + P(k2) P(k3)]"
    )?;
    writeln!(
        writer,
        "# shape: Fergusson-Shellard, k1 = kt (1 + alpha + beta) / 4, k2 = kt (1 - alpha + beta) / 4, k3 = kt (1 - beta) / 2; k columns are canonically ordered k1 >= k2 >= k3"
    )?;
    writeln!(
        writer,
        "# pivot: k_pivot = a H at N = {:e}, {:e} e-folds before the end of inflation (epsilon_H = 1 at N = {:e})",
        run.pivot.e_folds_at_pivot,
        settings.background.pivot_e_folds_before_end,
        run.pivot.e_folds_at_end
    )?;
    writeln!(
        writer,
        "# output_slice: N = {:e}, {:e} e-folds before the end of inflation; every leg has k/(aH) <= {:e}",
        run.e_folds_at_output,
        settings.background.output_e_folds_before_end,
        settings.maximum_output_k_over_ah
    )?;
    writeln!(
        writer,
        "# reheating: none; correlators are evaluated during inflation and include no post-inflationary transfer"
    )?;
    writeln!(
        writer,
        "# phase_space_units: reduced Planck units with comoving k in units of the arbitrary initial a = 1; Re/Im_k123sq_alpha = (k1 k2 k3)^2 Alpha, blocks QQQ..PPP with canonical momentum P"
    )
}

/// Evolves the canonical transport state of one triangle from its start slice
/// to `target_e_folds` on a surveyed trajectory.
pub(crate) fn evolve_three_point<P, G>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    trajectory: &BackgroundTrajectory,
    triangle: MomentumTriangle,
    e_folds_before_exit: f64,
    target_e_folds: f64,
    relative_tolerance: f64,
    max_steps: usize,
) -> Result<(f64, CanonicalThreePointTransportState), ThreePointEvolutionError>
where
    P: ThriceDifferentiableMultifieldPotential<2> + Copy,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2> + Copy,
{
    let depth = e_folds_before_exit.exp();
    let shortest = Mass::new(triangle.canonical_wavenumbers()[2]);
    let start = trajectory
        .crossing(|state| (shortest / (state.a * system.hubble(state))).into_value() - depth)
        .ok_or(ThreePointEvolutionError::StartNotFound)?;
    // The crossing is found by root refinement, so allow its rounding.
    let mut state = canonical_three_point_transport_initial_state(
        system,
        &BackgroundNeighborhood::from_trajectory(trajectory, start, 1.0e-3),
        triangle,
        depth * (1.0 - 1.0e-6),
    )
    .map_err(ThreePointEvolutionError::InitialCondition)?;
    let transport = CanonicalThreePointTransportSystem::new(system, trajectory, triangle);
    let tolerances = ErrorTolerances::new(
        CanonicalThreePointTransportTolerance::pure_relative(),
        relative_tolerance,
    );
    let limits = StepSizeLimits {
        min: Dimensionless::new(1.0e-9),
        max: Dimensionless::new(5.0e-2),
    };
    let mut e_folds = start;
    let mut step: Dimensionless = Dimensionless::new(1.0e-3);

    for _ in 0..max_steps {
        let remaining = target_e_folds - e_folds;
        if remaining <= 1.0e-13 {
            return Ok((start, state));
        }
        let attempted = Dimensionless::new(step.into_value().min(remaining));
        let outcome = adaptive_step(
            &transport,
            Dimensionless::new(e_folds),
            &state,
            attempted,
            &tolerances,
        );
        if outcome.accepted {
            state = outcome.state;
            e_folds += attempted.into_value();
        }
        step = next_step_size(attempted, outcome.error_ratio, &limits);
    }

    Err(ThreePointEvolutionError::StepLimit)
}

/// Failure of [`CovariantMultifieldSystem::three_point_state`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThreePointEvolutionError {
    /// The shortest leg never reached the requested subhorizon depth.
    StartNotFound,
    /// The subhorizon initial condition could not be built.
    InitialCondition(ThreePointInitialConditionError),
    /// The transport exhausted its step budget.
    StepLimit,
}

impl fmt::Display for ThreePointEvolutionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StartNotFound => {
                formatter.write_str("the triangle never reached its subhorizon start")
            }
            Self::InitialCondition(source) => {
                write!(
                    formatter,
                    "invalid three-point initial condition: {source:?}"
                )
            }
            Self::StepLimit => formatter.write_str("the three-point transport hit its step limit"),
        }
    }
}

impl std::error::Error for ThreePointEvolutionError {}

const BLOCK_NAMES: [&str; 8] = ["QQQ", "QQP", "QPQ", "PQQ", "QPP", "PQP", "PPQ", "PPP"];

fn phase_space_records(
    shape: usize,
    state: &CanonicalThreePointTransportState,
    wavenumbers: [f64; 3],
) -> Vec<PhaseSpaceBispectrumRecord> {
    let weight = (wavenumbers[0] * wavenumbers[1] * wavenumbers[2]).powi(2);
    let alpha = &state.alpha;
    let blocks: [[[[num_complex::Complex64; 2]; 2]; 2]; 8] = [
        alpha
            .qqq
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .qqp
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .qpq
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .pqq
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .qpp
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .pqp
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .ppq
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
        alpha
            .ppp
            .components()
            .map(|m| m.map(|r| r.map(|c| *c.value()))),
    ];
    blocks
        .iter()
        .zip(BLOCK_NAMES)
        .flat_map(|(block, name)| {
            (0..8).map(move |n| {
                let indices = [n / 4, (n / 2) % 2, n % 2];
                let value = block[indices[0]][indices[1]][indices[2]] * weight;
                PhaseSpaceBispectrumRecord {
                    shape,
                    block: name,
                    indices,
                    real: value.re,
                    imaginary: value.im,
                }
            })
        })
        .collect()
}

impl<P, G> CovariantMultifieldSystem<P, G, 2>
where
    P: ThriceDifferentiableMultifieldPotential<2> + Copy,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2> + Copy,
{
    /// The canonical two- and three-point state of `triangle` at
    /// `e_folds_at_output`, evolved on a surveyed `trajectory` from the
    /// subhorizon initial condition where the shortest leg has
    /// `k/(aH) = exp(e_folds_before_exit)`.
    pub fn three_point_state(
        &self,
        trajectory: &BackgroundTrajectory,
        triangle: MomentumTriangle,
        e_folds_before_exit: f64,
        e_folds_at_output: f64,
        relative_tolerance: f64,
        max_steps: usize,
    ) -> Result<CanonicalThreePointTransportState, ThreePointEvolutionError> {
        evolve_three_point(
            self,
            trajectory,
            triangle,
            e_folds_before_exit,
            e_folds_at_output,
            relative_tolerance,
            max_steps,
        )
        .map(|(_, state)| state)
    }

    /// Surveys the background, fixes `k_pivot = a H` as in the spectrum run,
    /// and projects the tree-level bispectrum of every shape on the output
    /// slice.
    pub fn run_bispectrum(
        &self,
        initial: MultifieldState<2>,
        settings: &BispectrumRunSettings,
        shapes: &[TriangleShape],
    ) -> Result<BispectrumRun, BispectrumRunError> {
        validate(&settings.background).map_err(BispectrumRunError::Spectrum)?;
        if !settings.e_folds_before_exit.is_finite() || settings.e_folds_before_exit <= 0.0 {
            return Err(BispectrumRunError::InvalidSettings(
                "e_folds_before_exit must be positive and finite",
            ));
        }
        if !(settings.relative_tolerance > 0.0 && settings.relative_tolerance < 1.0) {
            return Err(BispectrumRunError::InvalidSettings(
                "relative_tolerance must lie in (0, 1)",
            ));
        }
        let (_, trajectory, e_folds_at_end) = self
            .survey_background(initial, &settings.background)
            .map_err(BispectrumRunError::Spectrum)?;
        let e_folds_at_pivot = e_folds_at_end - settings.background.pivot_e_folds_before_end;
        let e_folds_at_output = e_folds_at_end - settings.background.output_e_folds_before_end;
        if e_folds_at_pivot < 0.0 || e_folds_at_output < 0.0 {
            return Err(BispectrumRunError::Spectrum(
                SpectrumRunError::NotEnoughEFolds {
                    available: e_folds_at_end,
                    required: settings
                        .background
                        .pivot_e_folds_before_end
                        .max(settings.background.output_e_folds_before_end),
                },
            ));
        }
        let pivot_state = trajectory.state_at(e_folds_at_pivot);
        let pivot_hubble = self.hubble(&pivot_state);
        let pivot = Pivot {
            k: pivot_state.a * pivot_hubble,
            hubble: pivot_hubble,
            epsilon_h: self.epsilon_h(&pivot_state),
            e_folds_at_end,
            e_folds_at_pivot,
        };
        let output_state = trajectory.state_at(e_folds_at_output);
        let output_hubble = self.hubble(&output_state);

        let mut records = Vec::with_capacity(shapes.len());
        let mut phase_space = Vec::new();
        for (index, shape) in shapes.iter().enumerate() {
            let [k1, k2, k3] = shape
                .wavenumbers_over_pivot()
                .map(|ratio| ratio * pivot.k.into_value());
            let triangle = MomentumTriangle::try_new(k1, k2, k3).map_err(|source| {
                BispectrumRunError::Triangle {
                    shape: index,
                    source,
                }
            })?;
            let (start_e_folds, state) = evolve_three_point(
                self,
                &trajectory,
                triangle,
                settings.e_folds_before_exit,
                e_folds_at_output,
                settings.relative_tolerance,
                settings.max_steps,
            )
            .map_err(|error| match error {
                ThreePointEvolutionError::StartNotFound => {
                    BispectrumRunError::StartNotFound { shape: index }
                }
                ThreePointEvolutionError::InitialCondition(source) => {
                    BispectrumRunError::InitialCondition {
                        shape: index,
                        source,
                    }
                }
                ThreePointEvolutionError::StepLimit => {
                    BispectrumRunError::StepLimit { shape: index }
                }
            })?;
            let projected = curvature_correlators(
                self,
                &output_state,
                triangle,
                &state,
                settings.maximum_output_k_over_ah,
            )
            .map_err(|source| BispectrumRunError::Projection {
                shape: index,
                source,
            })?;
            let wavenumbers = triangle.canonical_wavenumbers();
            let product: f64 = wavenumbers.iter().product();
            records.push(BispectrumRecord {
                shape: *shape,
                k_over_pivot: wavenumbers.map(|k| k / pivot.k.into_value()),
                e_folds: e_folds_at_output,
                max_k_over_ah: wavenumbers[0] / (output_state.a * output_hubble).into_value(),
                start_e_folds,
                dimensionless_power: projected.dimensionless_power(),
                dimensionless_bispectrum: product * product * projected.bispectrum.into_value(),
                f_nl: projected.reduced_f_nl(),
                intrinsic_fraction: (projected.intrinsic / projected.bispectrum).into_value(),
                imaginary_ratio: projected.bispectrum_imaginary_ratio,
            });
            phase_space.extend(phase_space_records(index, &state, wavenumbers));
        }

        Ok(BispectrumRun {
            pivot,
            e_folds_at_output,
            records,
            phase_space,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Field, FieldVelocity, FlatFieldSpace, QuadraticMultifieldPotential, ScaleFactor};

    fn turning_benchmark() -> (
        CovariantMultifieldSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace, 2>,
        MultifieldState<2>,
        BispectrumRunSettings,
    ) {
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(1.5e-6)]),
            FlatFieldSpace,
        );
        let resting = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(0.0); 2],
            ScaleFactor::new(1.0),
        );
        let initial = MultifieldState {
            pi: system.slow_roll_velocity(&resting),
            ..resting
        };
        let settings = BispectrumRunSettings {
            background: SpectrumRunSettings {
                pivot_e_folds_before_end: 6.0,
                output_e_folds_before_end: 2.0,
                k_min_over_pivot: 1.0,
                k_max_over_pivot: 1.0,
                k_count: 1,
                ..SpectrumRunSettings::default()
            },
            e_folds_before_exit: 3.0,
            relative_tolerance: 1.0e-8,
            ..BispectrumRunSettings::default()
        };
        (system, initial, settings)
    }

    #[test]
    fn shape_coordinates_follow_fergusson_and_shellard() {
        let equilateral = TriangleShape::equilateral(3.0).wavenumbers_over_pivot();
        for k in equilateral {
            assert!((k - 1.0).abs() < 1.0e-15);
        }
        let shape = TriangleShape {
            k_t_over_pivot: 4.0,
            alpha: 0.2,
            beta: 0.1,
        };
        let [k1, k2, k3] = shape.wavenumbers_over_pivot();
        assert!((k1 + k2 + k3 - 4.0).abs() < 1.0e-15);
        assert!((k1 - 1.3).abs() < 1.0e-15);
        assert!((k2 - 0.9).abs() < 1.0e-15);
        assert!((k3 - 1.8).abs() < 1.0e-15);
    }

    #[test]
    fn bispectrum_run_reports_consistent_labeled_columns() {
        let (system, initial, settings) = turning_benchmark();
        let shapes = [TriangleShape::equilateral(3.0)];
        let run = system.run_bispectrum(initial, &settings, &shapes).unwrap();
        let record = run.records[0];

        // Recompute f_NL from the dimensionless CSV columns alone.
        let [k1, k2, k3] = record.k_over_pivot;
        let [p1, p2, p3] = record.dimensionless_power;
        let power = |p: f64, k: f64| 2.0 * core::f64::consts::PI.powi(2) * p / k.powi(3);
        let bispectrum = record.dimensionless_bispectrum / (k1 * k2 * k3).powi(2);
        let pairs = power(p1, k1) * power(p2, k2)
            + power(p1, k1) * power(p3, k3)
            + power(p2, k2) * power(p3, k3);
        let f_nl = 5.0 / 6.0 * bispectrum / pairs;
        assert!((f_nl - record.f_nl).abs() < 1.0e-10 * record.f_nl.abs());
        assert!(record.max_k_over_ah < settings.maximum_output_k_over_ah);
        assert!(record.imaginary_ratio < 1.0e-6);
        assert_eq!(run.phase_space.len(), 64);

        let mut csv = Vec::new();
        run.write_csv(&mut csv).unwrap();
        let csv = String::from_utf8(csv).unwrap();
        let mut lines = csv.lines();
        let header = lines.next().unwrap().split(',').count();
        assert_eq!(lines.next().unwrap().split(',').count(), header);

        let mut metadata = Vec::new();
        write_bispectrum_metadata(&mut metadata, &settings, &run, "double quadratic", None)
            .unwrap();
        let metadata = String::from_utf8(metadata).unwrap();
        for key in [
            "# model:",
            "# gauge:",
            "# zeta_sign:",
            "# normalization:",
            "# f_nl_convention:",
            "# shape:",
            "# pivot:",
            "# output_slice:",
            "# reheating: none",
            "# phase_space_units:",
        ] {
            assert!(metadata.contains(key), "missing {key}");
        }
    }

    #[test]
    fn invalid_settings_and_shapes_are_rejected() {
        let (system, initial, settings) = turning_benchmark();
        let bad = BispectrumRunSettings {
            e_folds_before_exit: 0.0,
            ..settings
        };
        assert!(matches!(
            system.run_bispectrum(initial, &bad, &[TriangleShape::equilateral(3.0)]),
            Err(BispectrumRunError::InvalidSettings(_))
        ));
        let open = TriangleShape {
            k_t_over_pivot: 3.0,
            alpha: 0.0,
            beta: 1.2,
        };
        assert!(matches!(
            system.run_bispectrum(initial, &settings, &[open]),
            Err(BispectrumRunError::Triangle { shape: 0, .. })
        ));
    }

    #[test]
    #[ignore = "expensive; run with `cargo test --release -- --ignored`"]
    fn isosceles_shape_is_resolved_by_the_reported_grid() {
        // A straight trajectory, where f_NL varies smoothly with the shape.
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(1.0e-6)]),
            FlatFieldSpace,
        );
        let resting = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(0.0); 2],
            ScaleFactor::new(1.0),
        );
        let initial = MultifieldState {
            pi: system.slow_roll_velocity(&resting),
            ..resting
        };
        let settings = BispectrumRunSettings {
            background: SpectrumRunSettings {
                pivot_e_folds_before_end: 40.0,
                output_e_folds_before_end: 34.0,
                k_min_over_pivot: 1.0,
                k_max_over_pivot: 1.0,
                k_count: 1,
                ..SpectrumRunSettings::default()
            },
            e_folds_before_exit: 5.0,
            ..BispectrumRunSettings::default()
        };
        // beta from 0.1 to 0.7 at spacing 0.05; alpha = 0 is isosceles.
        let betas: Vec<f64> = (0..13).map(|n| 0.1 + 0.05 * f64::from(n)).collect();
        let shapes: Vec<TriangleShape> = betas
            .iter()
            .map(|&beta| TriangleShape {
                k_t_over_pivot: 3.0,
                alpha: 0.0,
                beta,
            })
            .collect();
        let run = system.run_bispectrum(initial, &settings, &shapes).unwrap();
        let f_nl: Vec<f64> = run.records.iter().map(|record| record.f_nl).collect();
        // Linear interpolation from a grid of stride `s` (in units of 0.05),
        // tested at the points that grid skips.
        let interpolation_error = |stride: usize| {
            (0..f_nl.len() - stride)
                .step_by(stride)
                .map(|left| {
                    let middle = left + stride / 2;
                    let predicted = 0.5 * (f_nl[left] + f_nl[left + stride]);
                    (predicted - f_nl[middle]).abs() / f_nl[middle].abs()
                })
                .fold(0.0, f64::max)
        };
        let coarse = interpolation_error(4);
        let fine = interpolation_error(2);
        eprintln!("f_NL {f_nl:?}; interpolation error {coarse:e} then {fine:e}");

        assert!(fine < 0.35 * coarse, "{coarse:e} then {fine:e}");
        assert!(fine < 1.0e-2, "spacing 0.1 misses the shape by {fine:e}");
    }
}
