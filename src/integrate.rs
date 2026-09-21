//! A small fixed-step Runge-Kutta integrator.
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
}
