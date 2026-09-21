//! Dimension-aware fixed and embedded Runge-Kutta integrators.
//!
//! The integrator knows nothing about potentials or Friedmann dynamics. It
//! only requires a system that can report the time derivative of a state, and
//! state types whose arithmetic keeps the physical dimensions consistent.

use core::ops::{Add, Mul};

/// A first-order system `dy/dt = f(y)` with dimension-checked state types.
pub trait System {
    /// The independent variable, such as cosmic time `t` or e-folding `N`.
    type IndependentVariable: Copy
        + Add<Output = Self::IndependentVariable>
        + Mul<f64, Output = Self::IndependentVariable>;

    /// The integrated state `y`.
    type State: Clone + Add<Self::Delta, Output = Self::State>;

    /// A finite change of the state over one step.
    type Delta;

    /// The time derivative `f(y)`.
    type Derivative: Clone
        + Add<Output = Self::Derivative>
        + Mul<f64, Output = Self::Derivative>
        + Mul<Self::IndependentVariable, Output = Self::Delta>;

    /// Evaluates the right-hand side at an independent-variable value and
    /// state.
    fn derivative(
        &self,
        independent: Self::IndependentVariable,
        state: &Self::State,
    ) -> Self::Derivative;
}

/// Advances `state` by one classical fourth-order Runge-Kutta step of size
/// `step` from the independent-variable value `independent`.
pub fn rk4_step<S: System>(
    system: &S,
    independent: S::IndependentVariable,
    state: &S::State,
    step: S::IndependentVariable,
) -> S::State {
    let half = step * 0.5;

    let k1 = system.derivative(independent, state);
    let k2 = system.derivative(independent + half, &(state.clone() + k1.clone() * half));
    let k3 = system.derivative(independent + half, &(state.clone() + k2.clone() * half));
    let k4 = system.derivative(independent + step, &(state.clone() + k3.clone() * step));

    let slope = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (1.0 / 6.0);

    state.clone() + slope * step
}

/// The fifth-order state and embedded fifth-minus-fourth-order error estimate
/// produced by one Dormand-Prince 5(4) step.
#[derive(Clone, Debug, PartialEq)]
pub struct EmbeddedStep<State, Delta> {
    /// The fifth-order solution estimate.
    pub state: State,
    /// The difference between the fifth- and fourth-order estimates.
    pub error: Delta,
}

/// Advances `state` by one Dormand-Prince 5(4) step.
///
/// The returned error is state-sized but deliberately has no built-in norm.
/// The caller chooses scales appropriate for each physical component before
/// deciding whether to accept the step.
pub fn rk45_step<S: System>(
    system: &S,
    independent: S::IndependentVariable,
    state: &S::State,
    step: S::IndependentVariable,
) -> EmbeddedStep<S::State, S::Delta> {
    let k1 = system.derivative(independent, state);
    let k2 = system.derivative(
        independent + step * (1.0 / 5.0),
        &(state.clone() + k1.clone() * (1.0 / 5.0) * step),
    );
    let k3 = system.derivative(
        independent + step * (3.0 / 10.0),
        &(state.clone() + (k1.clone() * (3.0 / 40.0) + k2.clone() * (9.0 / 40.0)) * step),
    );
    let k4 = system.derivative(
        independent + step * (4.0 / 5.0),
        &(state.clone()
            + (k1.clone() * (44.0 / 45.0)
                + k2.clone() * (-56.0 / 15.0)
                + k3.clone() * (32.0 / 9.0))
                * step),
    );
    let k5 = system.derivative(
        independent + step * (8.0 / 9.0),
        &(state.clone()
            + (k1.clone() * (19372.0 / 6561.0)
                + k2.clone() * (-25360.0 / 2187.0)
                + k3.clone() * (64448.0 / 6561.0)
                + k4.clone() * (-212.0 / 729.0))
                * step),
    );
    let k6 = system.derivative(
        independent + step,
        &(state.clone()
            + (k1.clone() * (9017.0 / 3168.0)
                + k2.clone() * (-355.0 / 33.0)
                + k3.clone() * (46732.0 / 5247.0)
                + k4.clone() * (49.0 / 176.0)
                + k5.clone() * (-5103.0 / 18656.0))
                * step),
    );

    let fifth_slope = k1.clone() * (35.0 / 384.0)
        + k3.clone() * (500.0 / 1113.0)
        + k4.clone() * (125.0 / 192.0)
        + k5.clone() * (-2187.0 / 6784.0)
        + k6.clone() * (11.0 / 84.0);
    let fifth = state.clone() + fifth_slope * step;
    let k7 = system.derivative(independent + step, &fifth);

    let error_slope = k1 * (71.0 / 57600.0)
        + k3 * (-71.0 / 16695.0)
        + k4 * (71.0 / 1920.0)
        + k5 * (-17253.0 / 339200.0)
        + k6 * (22.0 / 525.0)
        + k7 * (-1.0 / 40.0);

    EmbeddedStep {
        state: fifth,
        error: error_slope * step,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Time;

    #[derive(Clone)]
    struct NonCopyState(Box<f64>);

    struct NonCopyDelta(Box<f64>);

    #[derive(Clone)]
    struct NonCopyDerivative(Box<f64>);

    impl Add<NonCopyDelta> for NonCopyState {
        type Output = Self;

        fn add(self, rhs: NonCopyDelta) -> Self::Output {
            Self(Box::new(*self.0 + *rhs.0))
        }
    }

    impl Add for NonCopyDerivative {
        type Output = Self;

        fn add(self, rhs: Self) -> Self::Output {
            Self(Box::new(*self.0 + *rhs.0))
        }
    }

    impl Mul<f64> for NonCopyDerivative {
        type Output = Self;

        fn mul(self, rhs: f64) -> Self::Output {
            Self(Box::new(*self.0 * rhs))
        }
    }

    impl Mul<Time> for NonCopyDerivative {
        type Output = NonCopyDelta;

        fn mul(self, rhs: Time) -> Self::Output {
            NonCopyDelta(Box::new(*self.0 * rhs.into_value()))
        }
    }

    struct TimeDependentSystem;

    impl System for TimeDependentSystem {
        type IndependentVariable = Time;
        type State = NonCopyState;
        type Delta = NonCopyDelta;
        type Derivative = NonCopyDerivative;

        fn derivative(&self, time: Time, _state: &Self::State) -> Self::Derivative {
            NonCopyDerivative(Box::new(time.into_value()))
        }
    }

    #[test]
    fn supports_non_copy_states_and_passes_stage_time_to_the_system() {
        // dy/dt = t, y(0) = 0 has the exact solution y(1) = 1/2. RK4 is
        // exact for this linear-in-time derivative.
        let result = rk4_step(
            &TimeDependentSystem,
            Time::new(0.0),
            &NonCopyState(Box::new(0.0)),
            Time::new(1.0),
        );

        assert_eq!(*result.0, 0.5);
    }

    #[test]
    fn embedded_step_reports_a_small_error_for_a_smooth_solution() {
        struct ExponentialSystem;

        impl System for ExponentialSystem {
            type IndependentVariable = Time;
            type State = NonCopyState;
            type Delta = NonCopyDelta;
            type Derivative = NonCopyDerivative;

            fn derivative(&self, _time: Time, state: &Self::State) -> Self::Derivative {
                NonCopyDerivative(state.0.clone())
            }
        }

        let result = rk45_step(
            &ExponentialSystem,
            Time::new(0.0),
            &NonCopyState(Box::new(1.0)),
            Time::new(0.1),
        );

        assert!((*result.state.0 - 0.1_f64.exp()).abs() < 1.0e-8);
        assert!(result.error.0.abs() < 1.0e-7);
    }
}
