//! Driving the background integration and exporting its result.
//!
//! The crate does not plot. It produces rows that external analysis tools can
//! read. Stored values stay dimensionful for as long as possible;
//! dimensions are removed only at the export boundary, by dividing out the
//! appropriate power of the reduced Planck mass.

use std::io::{self, Write};

use crate::integrate::rk4_step;
use crate::potential::Potential;
use crate::scalar_field::ScalarFieldSystem;
use crate::state::State;
use crate::{Field, FieldVelocity, ScaleFactor, Time};

/// How far and how finely the background is integrated.
#[derive(Clone, Copy, Debug)]
pub struct RunSettings {
    /// The fixed Runge-Kutta step size.
    pub dt: Time,
    /// Integration stops once this e-folding number is reached.
    pub max_e_folds: f64,
    /// A safety limit on the number of steps.
    pub max_steps: usize,
    /// One row is recorded every `sample_every` steps.
    pub sample_every: usize,
}

impl Default for RunSettings {
    fn default() -> Self {
        Self {
            dt: Time::new(100.0),
            max_e_folds: 60.0,
            max_steps: 10_000_000,
            sample_every: 500,
        }
    }
}

/// Why the integration stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Termination {
    /// The requested e-folding number was reached.
    EFoldsReached,
    /// Inflation ended, that is `epsilon_H >= 1`.
    InflationEnded,
    /// The step limit was hit before either physical condition.
    StepLimitReached,
}

/// One exported row, in Planck units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Record {
    /// `t M_pl`.
    pub t: f64,
    /// The scale factor `a`.
    pub a: f64,
    /// `N = ln(a / a_initial)`.
    pub e_folds: f64,
    /// `phi / M_pl`.
    pub phi: f64,
    /// `(dphi/dt) / M_pl^2`.
    pub dphi: f64,
    /// `H / M_pl`.
    pub hubble: f64,
    /// `epsilon_H`.
    pub epsilon_h: f64,
}

impl Record {
    /// The CSV header matching [`Record::to_csv_row`].
    pub const CSV_HEADER: &'static str = "t,a,e_folds,phi,dphi,H,epsilon_h";

    /// Formats the row as comma-separated values.
    pub fn to_csv_row(&self) -> String {
        format!(
            "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
            self.t, self.a, self.e_folds, self.phi, self.dphi, self.hubble, self.epsilon_h
        )
    }
}

/// The result of integrating a background.
#[derive(Clone, Debug)]
pub struct BackgroundRun {
    /// The sampled rows, including the first and last state.
    pub records: Vec<Record>,
    /// Why the integration stopped.
    pub termination: Termination,
    /// The state at the final step.
    pub final_state: State,
    /// The cosmic time at the final step.
    pub final_time: Time,
    /// The number of steps actually taken.
    pub steps: usize,
}

impl BackgroundRun {
    /// Writes the header and every row as CSV.
    pub fn write_csv<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "{}", Record::CSV_HEADER)?;
        for record in &self.records {
            writeln!(writer, "{}", record.to_csv_row())?;
        }
        Ok(())
    }
}

impl<P: Potential> ScalarFieldSystem<P> {
    /// Builds an initial state from `phi` and `a`, taking the field velocity
    /// from the slow-roll estimate `dphi/dt ~= -V_phi / (3 H)`.
    pub fn slow_roll_initial_state(&self, phi: Field, a: ScaleFactor) -> State {
        let resting = State::new(phi, FieldVelocity::new(0.0), a);

        State {
            pi: self.slow_roll_velocity(&resting),
            ..resting
        }
    }

    /// Records one exported row for `state` at elapsed time `t`.
    fn record(&self, initial: &State, state: &State, t: Time) -> Record {
        Record {
            t: (t * self.m_pl).into_value(),
            a: state.a.into_value(),
            e_folds: state.e_folds(initial),
            phi: (state.phi / self.m_pl).into_value(),
            dphi: (state.pi / (self.m_pl * self.m_pl)).into_value(),
            hubble: (self.hubble(state) / self.m_pl).into_value(),
            epsilon_h: self.epsilon_h(state),
        }
    }

    /// Integrates the background with fixed-step RK4 until the requested
    /// e-folding number is reached or inflation ends.
    pub fn run(&self, initial: State, settings: RunSettings) -> BackgroundRun {
        let mut state = initial;
        let mut time = Time::new(0.0);
        let mut records = vec![self.record(&initial, &state, time)];
        let mut termination = Termination::StepLimitReached;
        let mut steps = 0;

        for step in 1..=settings.max_steps {
            state = rk4_step(self, time, &state, settings.dt);
            time = time + settings.dt;
            steps = step;

            let e_folds = state.e_folds(&initial);
            let epsilon_h = self.epsilon_h(&state);

            let finished = if epsilon_h >= 1.0 {
                termination = Termination::InflationEnded;
                true
            } else if e_folds >= settings.max_e_folds {
                termination = Termination::EFoldsReached;
                true
            } else {
                false
            };

            if finished || step % settings.sample_every == 0 {
                records.push(self.record(&initial, &state, time));
            }

            if finished {
                break;
            }
        }

        BackgroundRun {
            records,
            termination,
            final_state: state,
            final_time: time,
            steps,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Mass;
    use crate::potential::QuadraticPotential;

    fn quadratic_run() -> BackgroundRun {
        let system = ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(1.0e-6)));
        let initial = system.slow_roll_initial_state(Field::new(16.0), ScaleFactor::new(1.0));

        system.run(initial, RunSettings::default())
    }

    #[test]
    fn a_quadratic_background_reaches_sixty_e_folds() {
        let run = quadratic_run();

        assert_eq!(run.termination, Termination::EFoldsReached);
        assert!(run.records.len() > 2);

        let last = run.records.last().expect("at least one record");
        assert!((last.e_folds - 60.0).abs() < 0.1, "N was {}", last.e_folds);
    }

    #[test]
    fn e_folds_match_the_slow_roll_prediction() {
        let run = quadratic_run();
        let last = run.records.last().expect("at least one record");

        // For V = (1/2) m^2 phi^2 slow roll gives N = (phi_i^2 - phi_f^2) / 4.
        let predicted = (16.0 * 16.0 - last.phi * last.phi) / 4.0;

        assert!(
            (predicted - last.e_folds).abs() < 0.5,
            "slow roll predicted {predicted}, integration gave {}",
            last.e_folds
        );
    }

    #[test]
    fn inflation_ends_when_started_close_to_the_end() {
        let system = ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(1.0e-6)));
        let initial = system.slow_roll_initial_state(Field::new(3.0), ScaleFactor::new(1.0));

        let run = system.run(initial, RunSettings::default());

        assert_eq!(run.termination, Termination::InflationEnded);
        let last = run.records.last().expect("at least one record");
        assert!(last.epsilon_h >= 1.0);
        // Slow roll puts the end of inflation at phi = sqrt(2) M_pl, but the
        // exact condition epsilon_H = 1 is reached slightly later, near
        // phi = 1.0 M_pl.
        assert!((last.phi - 1.0).abs() < 0.1, "phi was {}", last.phi);
    }

    #[test]
    fn csv_output_starts_with_the_documented_header() {
        let run = quadratic_run();
        let mut buffer = Vec::new();

        run.write_csv(&mut buffer).expect("writing to a Vec");

        let text = String::from_utf8(buffer).expect("valid UTF-8");
        let mut lines = text.lines();

        assert_eq!(lines.next(), Some("t,a,e_folds,phi,dphi,H,epsilon_h"));
        assert_eq!(lines.next().expect("a first row").split(',').count(), 7);
    }

    #[test]
    fn final_time_matches_the_number_of_fixed_steps() {
        let run = quadratic_run();
        let expected = RunSettings::default().dt.into_value() * run.steps as f64;

        assert_eq!(run.final_time.into_value(), expected);
    }
}
