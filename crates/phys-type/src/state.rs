//! The background state of a single-field cosmology and its time derivative.
//!
//! A state, its time derivative, and a finite change over one step are three
//! physically different things, so each is a distinct type. The only way to
//! turn a derivative into a change is to multiply it by a time, and the only
//! way to advance a state is to add such a change:
//!
//! ```text
//! Time * StateDerivative -> StateDelta
//! State + StateDelta     -> State
//! ```

use core::ops::{Add, Mul};

use crate::{Field, FieldVelocity, M1, M2, M3, Quantity, ScaleFactor, Time};

/// The integrated background variables at one instant of cosmic time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct State {
    /// The scalar field `phi`.
    pub phi: Field,
    /// The field velocity `pi = dphi/dt`.
    pub pi: FieldVelocity,
    /// The scale factor `a`.
    pub a: ScaleFactor,
}

/// The time derivative of a [`State`].
///
/// Every component carries one more mass power than the corresponding state
/// component, because differentiating by time divides by `M^-1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StateDerivative {
    /// `dphi/dt`.
    pub dphi_dt: Quantity<M2, f64>,
    /// `dpi/dt`.
    pub dpi_dt: Quantity<M3, f64>,
    /// `da/dt`.
    pub da_dt: Quantity<M1, f64>,
}

/// A finite change of a [`State`] over one step.
///
/// Each component shares the dimension of the state component it updates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StateDelta {
    /// The change of `phi`.
    pub phi: Field,
    /// The change of `pi`.
    pub pi: FieldVelocity,
    /// The change of `a`.
    pub a: ScaleFactor,
}

impl StateDelta {
    /// Attaches the same raw Planck-unit number to every component.
    ///
    /// This is the convenient default, not the safe one: it is exactly the
    /// assumption that one number suits every dimension, which holds only
    /// while a state's components stay within a few orders of each other.
    pub const fn from_planck_units(value: f64) -> Self {
        Self {
            phi: Field::new(value),
            pi: FieldVelocity::new(value),
            a: ScaleFactor::new(value),
        }
    }
}

impl State {
    /// Builds a state from its three background variables.
    pub const fn new(phi: Field, pi: FieldVelocity, a: ScaleFactor) -> Self {
        Self { phi, pi, a }
    }

    /// The e-folding number `N = ln(a / a_initial)`.
    pub fn e_folds(&self, initial: &Self) -> f64 {
        (self.a / initial.a).into_value().ln()
    }
}

impl Add<StateDelta> for State {
    type Output = State;

    fn add(self, rhs: StateDelta) -> Self::Output {
        State {
            phi: self.phi + rhs.phi,
            pi: self.pi + rhs.pi,
            a: self.a + rhs.a,
        }
    }
}

impl Add for StateDerivative {
    type Output = StateDerivative;

    fn add(self, rhs: StateDerivative) -> Self::Output {
        StateDerivative {
            dphi_dt: self.dphi_dt + rhs.dphi_dt,
            dpi_dt: self.dpi_dt + rhs.dpi_dt,
            da_dt: self.da_dt + rhs.da_dt,
        }
    }
}

/// Scaling a derivative by a dimensionless weight, as Runge-Kutta requires.
impl Mul<f64> for StateDerivative {
    type Output = StateDerivative;

    fn mul(self, rhs: f64) -> Self::Output {
        StateDerivative {
            dphi_dt: self.dphi_dt * rhs,
            dpi_dt: self.dpi_dt * rhs,
            da_dt: self.da_dt * rhs,
        }
    }
}

/// The only bridge from a derivative to a state-sized change.
impl Mul<Time> for StateDerivative {
    type Output = StateDelta;

    fn mul(self, rhs: Time) -> Self::Output {
        StateDelta {
            phi: self.dphi_dt * rhs,
            pi: self.dpi_dt * rhs,
            a: self.da_dt * rhs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_state() -> State {
        State::new(
            Field::new(16.0),
            FieldVelocity::new(-1.0e-6),
            ScaleFactor::new(1.0),
        )
    }

    #[test]
    fn a_derivative_times_a_time_is_a_state_sized_change() {
        let derivative = StateDerivative {
            dphi_dt: Quantity::new(2.0),
            dpi_dt: Quantity::new(4.0),
            da_dt: Quantity::new(8.0),
        };

        let delta: StateDelta = derivative * Time::new(0.5);

        assert_eq!(delta.phi.into_value(), 1.0);
        assert_eq!(delta.pi.into_value(), 2.0);
        assert_eq!(delta.a.into_value(), 4.0);
    }

    #[test]
    fn adding_a_change_advances_the_state() {
        let state = sample_state();
        let delta = StateDelta {
            phi: Field::new(-0.5),
            pi: FieldVelocity::new(1.0e-7),
            a: ScaleFactor::new(0.25),
        };

        let advanced = state + delta;

        assert_eq!(advanced.phi.into_value(), 15.5);
        assert_eq!(advanced.a.into_value(), 1.25);
    }

    #[test]
    fn derivatives_combine_with_runge_kutta_weights() {
        let derivative = StateDerivative {
            dphi_dt: Quantity::new(1.0),
            dpi_dt: Quantity::new(1.0),
            da_dt: Quantity::new(1.0),
        };

        let combined = derivative + derivative * 2.0;

        assert_eq!(combined.dphi_dt.into_value(), 3.0);
    }

    #[test]
    fn e_folds_are_derived_from_the_scale_factor() {
        let initial = sample_state();
        let mut later = initial;
        later.a = ScaleFactor::new(core::f64::consts::E);

        assert!((later.e_folds(&initial) - 1.0).abs() < 1.0e-12);
    }
}
