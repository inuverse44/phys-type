//! A small fixed-step Runge-Kutta integrator.
//!
//! The integrator knows nothing about potentials or Friedmann dynamics. It
//! only requires a system that can report the time derivative of a state, and
//! state types whose arithmetic keeps the physical dimensions consistent.

use core::ops::{Add, Mul};

use crate::Time;

/// A first-order system `dy/dt = f(y)` with dimension-checked state types.
pub trait System {
    /// The integrated state `y`.
    type State: Copy + Add<Self::Delta, Output = Self::State>;

    /// A finite change of the state over one step.
    type Delta;

    /// The time derivative `f(y)`.
    type Derivative: Copy
        + Add<Output = Self::Derivative>
        + Mul<f64, Output = Self::Derivative>
        + Mul<Time, Output = Self::Delta>;

    /// Evaluates the right-hand side at `state`.
    fn derivative(&self, state: &Self::State) -> Self::Derivative;
}

/// Advances `state` by one classical fourth-order Runge-Kutta step of size
/// `dt`.
pub fn rk4_step<S: System>(system: &S, state: &S::State, dt: Time) -> S::State {
    let half = dt * 0.5;

    let k1 = system.derivative(state);
    let k2 = system.derivative(&(*state + k1 * half));
    let k3 = system.derivative(&(*state + k2 * half));
    let k4 = system.derivative(&(*state + k3 * dt));

    let slope = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (1.0 / 6.0);

    *state + slope * dt
}
