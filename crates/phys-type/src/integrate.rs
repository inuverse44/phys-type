//! Dimension-aware fixed and embedded Runge-Kutta integrators.
//!
//! The integrator knows nothing about potentials or Friedmann dynamics. It
//! only requires a system that can report the time derivative of a state, and
//! state types whose arithmetic keeps the physical dimensions consistent.

use core::ops::{Add, Mul};

use crate::Quantity;

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

/// Tolerances that turn an embedded error estimate into a dimensionless
/// ratio.
///
/// `A` is a per-system container holding one absolute tolerance for each
/// group of components, each carrying that group's own physical dimension.
/// A single number cannot serve a whole state: in a spectrum run the scale
/// factor and the mode amplitudes are tens of orders of magnitude apart, and
/// one floor placed for the sake of the first silently caps the accuracy of
/// the second.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ErrorTolerances<A> {
    /// Absolute tolerances, in the dimensions of the components they bound.
    pub absolute: A,
    /// Relative tolerance, applied to every component alike.
    pub relative: f64,
}

impl<A> ErrorTolerances<A> {
    /// Pairs a set of absolute tolerances with one relative tolerance.
    pub const fn new(absolute: A, relative: f64) -> Self {
        Self { absolute, relative }
    }
}

impl<A: Default> Default for ErrorTolerances<A> {
    fn default() -> Self {
        Self {
            absolute: A::default(),
            relative: 1.0e-8,
        }
    }
}

/// The normalized error `|e| / (a_tol + r_tol max(|y|, |y_next|))` of one
/// real component.
///
/// Every argument but the relative tolerance carries the component's own
/// dimension `D`, so a tolerance meant for one quantity cannot be used to
/// bound another. The returned ratio is dimensionless.
///
/// Both endpoints of the step enter the scale, so a component passing through
/// zero is still bounded by the absolute tolerance. A zero scale means the
/// component is being controlled purely relatively and has vanished; a zero
/// error there is exact, and anything else is unbounded.
pub fn component_error_ratio<D>(
    error: Quantity<D, f64>,
    current: Quantity<D, f64>,
    next: Quantity<D, f64>,
    absolute: Quantity<D, f64>,
    relative: f64,
) -> f64 {
    let error = error.into_value();
    let scale =
        absolute.into_value() + relative * current.into_value().abs().max(next.into_value().abs());
    if scale == 0.0 {
        return if error == 0.0 { 0.0 } else { f64::INFINITY };
    }

    error.abs() / scale
}

/// A system whose embedded step error can be normalized to a single ratio.
///
/// The integrator deliberately has no built-in norm, because only the system
/// knows which components exist and what scale each one carries. An
/// implementation reduces the whole state to the largest normalized component
/// error, so one acceptance decision covers the entire step.
pub trait AdaptiveSystem: System {
    /// One absolute tolerance per group of components, each in that group's
    /// dimension.
    ///
    /// Giving this a type rather than a number is what stops a tolerance
    /// meant for one quantity from bounding another.
    type AbsoluteTolerance;

    /// The largest normalized component error of one embedded step.
    ///
    /// A step is accepted when this does not exceed one.
    fn error_ratio(
        &self,
        current: &Self::State,
        next: &Self::State,
        error: &Self::Delta,
        tolerances: &ErrorTolerances<Self::AbsoluteTolerance>,
    ) -> f64;
}

/// Bounds on the step size an adaptive controller may propose.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepSizeLimits<I> {
    /// Smallest permitted step size.
    pub min: I,
    /// Largest permitted step size.
    pub max: I,
}

/// One attempted adaptive step and its acceptance decision.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveAttempt<State> {
    /// The fifth-order solution estimate, whether or not it was accepted.
    pub state: State,
    /// The largest normalized component error of the attempt.
    pub error_ratio: f64,
    /// Whether the attempt met the requested tolerances.
    pub accepted: bool,
}

/// Takes one embedded step and normalizes its error.
///
/// The caller decides what to do with a rejected attempt; this function only
/// reports it.
pub fn adaptive_step<S: AdaptiveSystem>(
    system: &S,
    independent: S::IndependentVariable,
    state: &S::State,
    step: S::IndependentVariable,
    tolerances: &ErrorTolerances<S::AbsoluteTolerance>,
) -> AdaptiveAttempt<S::State> {
    let attempt = rk45_step(system, independent, state, step);
    let error_ratio = system.error_ratio(state, &attempt.state, &attempt.error, tolerances);

    AdaptiveAttempt {
        state: attempt.state,
        error_ratio,
        accepted: error_ratio <= 1.0,
    }
}

/// The step size suggested by the standard fifth-order controller.
///
/// A zero error grows the step by the largest permitted factor, and a
/// non-finite error shrinks it by the smallest one. The result is kept inside
/// `limits`.
pub fn next_step_size<I>(step: I, error_ratio: f64, limits: &StepSizeLimits<I>) -> I
where
    I: Copy + Mul<f64, Output = I> + PartialOrd,
{
    let factor = if error_ratio == 0.0 {
        5.0
    } else if error_ratio.is_finite() {
        (0.9 * error_ratio.powf(-0.2)).clamp(0.2, 5.0)
    } else {
        0.2
    };
    let proposed = step * factor;

    if proposed < limits.min {
        limits.min
    } else if proposed > limits.max {
        limits.max
    } else {
        proposed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Dimensionless, Time};

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

    impl AdaptiveSystem for ExponentialSystem {
        type AbsoluteTolerance = f64;

        fn error_ratio(
            &self,
            current: &Self::State,
            next: &Self::State,
            error: &Self::Delta,
            tolerances: &ErrorTolerances<f64>,
        ) -> f64 {
            component_error_ratio(
                Dimensionless::new(*error.0),
                Dimensionless::new(*current.0),
                Dimensionless::new(*next.0),
                Dimensionless::new(tolerances.absolute),
                tolerances.relative,
            )
        }
    }

    #[test]
    fn a_small_step_of_a_smooth_system_is_accepted() {
        let attempt = adaptive_step(
            &ExponentialSystem,
            Time::new(0.0),
            &NonCopyState(Box::new(1.0)),
            Time::new(0.01),
            &ErrorTolerances::new(1.0e-12, 1.0e-8),
        );

        assert!(attempt.accepted);
        assert!(attempt.error_ratio <= 1.0);
        assert!((*attempt.state.0 - 0.01_f64.exp()).abs() < 1.0e-12);
    }

    #[test]
    fn a_step_that_is_far_too_large_is_rejected() {
        let attempt = adaptive_step(
            &ExponentialSystem,
            Time::new(0.0),
            &NonCopyState(Box::new(1.0)),
            Time::new(5.0),
            &ErrorTolerances::new(1.0e-12, 1.0e-8),
        );

        assert!(!attempt.accepted);
        assert!(attempt.error_ratio > 1.0);
    }

    #[test]
    fn the_controller_grows_shrinks_and_clamps_the_step() {
        let limits = StepSizeLimits {
            min: Time::new(0.5),
            max: Time::new(4.0),
        };
        let step = Time::new(1.0);

        // A vanishing error takes the largest growth factor, clamped by max.
        assert_eq!(next_step_size(step, 0.0, &limits), Time::new(4.0));
        // A large error shrinks the step, clamped by min.
        assert_eq!(next_step_size(step, 1.0e6, &limits), Time::new(0.5));
        // A non-finite error is treated as the worst case.
        assert_eq!(next_step_size(step, f64::NAN, &limits), Time::new(0.5));
        // An error just at tolerance barely changes the step.
        assert!(next_step_size(step, 1.0, &limits) < step);
    }

    #[test]
    fn a_component_passing_through_zero_is_bounded_by_the_absolute_tolerance() {
        let ratio = |error: f64, absolute: f64| {
            component_error_ratio(
                Dimensionless::new(error),
                Dimensionless::new(0.0),
                Dimensionless::new(0.0),
                Dimensionless::new(absolute),
                1.0e-8,
            )
        };

        assert!(ratio(1.0e-13, 1.0e-12) < 1.0);
        assert!(ratio(1.0e-11, 1.0e-12) > 1.0);
    }

    #[test]
    fn a_vanished_component_with_no_absolute_floor_is_exact_or_unbounded() {
        let ratio = |error: f64| {
            component_error_ratio(
                Dimensionless::new(error),
                Dimensionless::new(0.0),
                Dimensionless::new(0.0),
                Dimensionless::new(0.0),
                1.0e-8,
            )
        };

        assert_eq!(ratio(0.0), 0.0);
        assert_eq!(ratio(1.0e-30), f64::INFINITY);
    }

    #[test]
    fn embedded_step_reports_a_small_error_for_a_smooth_solution() {
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
