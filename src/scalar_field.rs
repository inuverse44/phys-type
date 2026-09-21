//! The single-field background system.
//!
//! The evolved variables are `phi`, `pi = dphi/dt` and the scale factor `a`:
//!
//! ```text
//! dphi/dt = pi
//! dpi/dt  = -3 H pi - V_phi
//! da/dt   = a H
//! ```
//!
//! `H` is not integrated independently. The expanding branch is recomputed
//! from the current state at every right-hand-side evaluation:
//!
//! ```text
//! H = +sqrt(((1/2) pi^2 + V(phi)) / (3 M_pl^2))
//! ```

use crate::integrate::{AdaptiveSystem, ErrorTolerances, System, component_error_ratio};
use crate::potential::Potential;
use crate::state::{State, StateDelta, StateDerivative};
use crate::{EnergyDensity, Hubble, Mass};

/// A scalar field evolving in a spatially flat FLRW background.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScalarFieldSystem<P> {
    /// The potential that drives the field.
    pub potential: P,
    /// The reduced Planck mass, numerically `1.0` in Planck units but still
    /// carrying mass dimension `M^1` in the type system.
    pub m_pl: Mass,
}

impl<P: Potential> ScalarFieldSystem<P> {
    /// Builds a system with the given potential and reduced Planck mass.
    pub const fn new(potential: P, m_pl: Mass) -> Self {
        Self { potential, m_pl }
    }

    /// Builds a system in Planck units, where `M_pl = 1`.
    pub const fn in_planck_units(potential: P) -> Self {
        Self::new(potential, Mass::new(1.0))
    }

    /// The total energy density `(1/2) pi^2 + V(phi)`.
    pub fn energy_density(&self, state: &State) -> EnergyDensity {
        0.5 * state.pi * state.pi + self.potential.value(state.phi)
    }

    /// The Hubble rate on the expanding branch of the Friedmann equation.
    pub fn hubble(&self, state: &State) -> Hubble {
        (self.energy_density(state) / (self.m_pl * self.m_pl * 3.0)).sqrt()
    }

    /// The first Hubble slow-roll parameter
    /// `epsilon_H = pi^2 / (2 M_pl^2 H^2)`.
    ///
    /// Inflation ends when this reaches `1`.
    pub fn epsilon_h(&self, state: &State) -> f64 {
        let hubble = self.hubble(state);
        let numerator = state.pi * state.pi;
        let denominator = self.m_pl * self.m_pl * hubble * hubble * 2.0;

        (numerator / denominator).into_value()
    }

    /// The scalar field equation-of-state parameter
    /// `w = (K - V) / (K + V)`, where `K = (1/2) pi^2`.
    pub fn equation_of_state(&self, state: &State) -> f64 {
        let kinetic = state.pi * state.pi * 0.5;
        let potential = self.potential.value(state.phi);

        ((kinetic - potential) / (kinetic + potential)).into_value()
    }

    /// The slow-roll estimate `dphi/dt ~= -V_phi / (3 H)` used to initialise
    /// the field velocity from `phi` alone.
    pub fn slow_roll_velocity(&self, state: &State) -> crate::FieldVelocity {
        let hubble = self.hubble(state);

        -(self.potential.gradient(state.phi) / (hubble * 3.0))
    }
}

impl<P: Potential> System for ScalarFieldSystem<P> {
    type IndependentVariable = crate::Time;
    type State = State;
    type Delta = StateDelta;
    type Derivative = StateDerivative;

    fn derivative(&self, _time: crate::Time, state: &State) -> StateDerivative {
        let hubble = self.hubble(state);

        StateDerivative {
            dphi_dt: state.pi,
            dpi_dt: state.pi * hubble * (-3.0) - self.potential.gradient(state.phi),
            da_dt: state.a * hubble,
        }
    }
}

impl<P: Potential> AdaptiveSystem for ScalarFieldSystem<P> {
    type AbsoluteTolerance = StateDelta;

    fn error_ratio(
        &self,
        current: &State,
        next: &State,
        error: &StateDelta,
        tolerances: &ErrorTolerances<StateDelta>,
    ) -> f64 {
        let relative = tolerances.relative;
        let absolute = &tolerances.absolute;

        component_error_ratio(
            error.phi.into_value(),
            current.phi.into_value(),
            next.phi.into_value(),
            absolute.phi.into_value(),
            relative,
        )
        .max(component_error_ratio(
            error.pi.into_value(),
            current.pi.into_value(),
            next.pi.into_value(),
            absolute.pi.into_value(),
            relative,
        ))
        .max(component_error_ratio(
            error.a.into_value(),
            current.a.into_value(),
            next.a.into_value(),
            absolute.a.into_value(),
            relative,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::QuadraticPotential;
    use crate::{Field, FieldVelocity, ScaleFactor, Time, rk4_step};

    fn slow_roll_start(phi: f64) -> (ScalarFieldSystem<QuadraticPotential>, State) {
        let system = ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(1.0e-6)));
        let resting = State::new(
            Field::new(phi),
            FieldVelocity::new(0.0),
            ScaleFactor::new(1.0),
        );
        let state = State {
            pi: system.slow_roll_velocity(&resting),
            ..resting
        };

        (system, state)
    }

    #[test]
    fn hubble_matches_the_friedmann_equation() {
        let (system, _) = slow_roll_start(16.0);
        let state = State::new(
            Field::new(16.0),
            FieldVelocity::new(0.0),
            ScaleFactor::new(1.0),
        );

        // With pi = 0 the density is (1/2) m^2 phi^2, so H = m phi / sqrt(6).
        let expected = 1.0e-6 * 16.0 / 6.0_f64.sqrt();

        assert!((system.hubble(&state).into_value() - expected).abs() < 1.0e-18);
    }

    #[test]
    fn slow_roll_start_is_deep_in_the_inflating_regime() {
        let (system, state) = slow_roll_start(16.0);

        let epsilon = system.epsilon_h(&state);

        assert!(state.pi.into_value() < 0.0, "the field rolls downhill");
        assert!(epsilon > 0.0 && epsilon < 0.01, "epsilon_H was {epsilon}");
    }

    #[test]
    fn equation_of_state_distinguishes_potential_and_kinetic_energy() {
        let (system, _) = slow_roll_start(16.0);
        let potential_dominated = State::new(
            Field::new(16.0),
            FieldVelocity::new(0.0),
            ScaleFactor::new(1.0),
        );
        let kinetic_dominated = State::new(
            Field::new(0.0),
            FieldVelocity::new(1.0e-6),
            ScaleFactor::new(1.0),
        );

        assert_eq!(system.equation_of_state(&potential_dominated), -1.0);
        assert_eq!(system.equation_of_state(&kinetic_dominated), 1.0);
    }

    #[test]
    fn a_single_step_expands_the_universe_and_lowers_the_field() {
        let (system, state) = slow_roll_start(16.0);

        let next = rk4_step(&system, Time::new(0.0), &state, Time::new(100.0));

        assert!(next.a.into_value() > state.a.into_value());
        assert!(next.phi.into_value() < state.phi.into_value());
    }

    #[test]
    fn de_sitter_limit_reproduces_exponential_expansion() {
        // A field held at rest by a flat potential expands as a = exp(H t).
        #[derive(Clone, Copy)]
        struct FlatPotential {
            density: EnergyDensity,
        }

        impl<V> crate::PotentialValue<V> for FlatPotential
        where
            V: From<f64>,
        {
            fn value(&self, _phi: Field<V>) -> EnergyDensity<V> {
                EnergyDensity::new(self.density.into_value().into())
            }
        }

        let system = ScalarFieldSystem::in_planck_units(FlatPotential {
            density: EnergyDensity::new(3.0e-12),
        });
        let initial = State::new(
            Field::new(1.0),
            FieldVelocity::new(0.0),
            ScaleFactor::new(1.0),
        );

        let hubble = system.hubble(&initial).into_value();
        let dt = Time::new(1.0e3);
        let steps = 1_000;

        let mut state = initial;
        let mut time = Time::new(0.0);
        for _ in 0..steps {
            state = rk4_step(&system, time, &state, dt);
            time = time + dt;
        }

        let elapsed = 1.0e3 * steps as f64;
        let expected = (hubble * elapsed).exp();
        let relative_error = (state.a.into_value() - expected).abs() / expected;

        assert!(relative_error < 1.0e-10, "relative error {relative_error}");
    }
}
