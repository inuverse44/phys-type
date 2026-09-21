//! Linear two-field perturbations evolved with their background.

use core::ops::{Add, Mul};

use num_complex::Complex64;

use crate::integrate::{AdaptiveSystem, ErrorTolerances, System};
use crate::{
    CovariantMultifieldSystem, CurvedFieldSpaceGeometry, Dimensionless, KinematicBasisError, M1,
    M3Over2, MHalf, MMinusHalf, Mass, MultifieldState, MultifieldStateDelta,
    MultifieldStateDerivative, Quantity, Time, TwiceDifferentiableMultifieldPotential,
};

/// Rescaled Fourier mode amplitude `v^I_A = a U^I_A`, with dimension fixed by
/// `P_Q = k^3 |U|^2 / (2 pi^2)`.
///
/// The rescaling is dimensionless, so this is also the dimension of the
/// unrescaled `U^I_A`. Only the rescaled variable is ever integrated.
pub type ModeAmplitude<V = Complex64> = Quantity<MMinusHalf, V>;
/// Covariant time derivative `D_t v^I_A`.
pub type ModeVelocity<V = Complex64> = Quantity<MHalf, V>;
/// Covariant second derivative of a rescaled mode amplitude.
pub type ModeAcceleration<V = Complex64> = Quantity<M3Over2, V>;
/// A dimensionless complex `2 x 2` diagnostic matrix.
pub type ComplexMatrix2 = [[Complex64; 2]; 2];

/// Converts a dimensioned complex value to a real magnitude for error norms.
///
/// The caller supplies an absolute tolerance with the same dimension before
/// combining this magnitude into a dimensionless normalized error.
pub trait ComplexMagnitude {
    /// Returns the complex modulus of the stored numerical value.
    fn complex_magnitude(&self) -> f64;
}

impl<D> ComplexMagnitude for Quantity<D, Complex64> {
    fn complex_magnitude(&self) -> f64 {
        self.value().norm()
    }
}

/// The dimensionless normalized error of one complex component,
/// `|dz| / (a_tol + r_tol |z|)`.
///
/// The absolute tolerance carries the component's own dimension, so a
/// tolerance meant for a mode amplitude cannot be applied to a mode velocity.
/// The relative tolerance and the returned error are dimensionless. Both the
/// current and the advanced value enter the scale, matching the background
/// integrator, and the modulus keeps a single acceptance decision for the real
/// and imaginary parts together.
pub fn normalized_mode_error<D>(
    error: Quantity<D, Complex64>,
    current: Quantity<D, Complex64>,
    next: Quantity<D, Complex64>,
    absolute_tolerance: Quantity<D, f64>,
    relative_tolerance: f64,
) -> f64 {
    let scale = current.complex_magnitude().max(next.complex_magnitude());

    error.complex_magnitude() / (absolute_tolerance.into_value() + relative_tolerance * scale)
}

/// The two independent complex two-field mode functions, rescaled by the
/// scale factor.
///
/// The integrated variable is `v^I_A = a U^I_A`. Inside the horizon `U`
/// carries a factor `1 / (a sqrt(2k))` and so decays by many orders of
/// magnitude before a mode exits, which forces the absolute error tolerance
/// far below anything the background needs. The rescaled amplitude instead
/// oscillates at fixed size, and its Hubble drag drops from `3H` to `H`.
///
/// Recovering `U` divides by `a`, which
/// [`crate::spectra::mode_covariance`] does once, at the one place where the
/// spectra are formed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldModes {
    /// `v^I_A = a U^I_A`; the first index is field space and the second
    /// labels the independent vacuum oscillator.
    pub rescaled_amplitude: [[ModeAmplitude; 2]; 2],
    /// `D_t v^I_A` with the same index order.
    pub rescaled_velocity: [[ModeVelocity; 2]; 2],
}

/// Background and mode matrix for one comoving wavenumber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldPerturbationState {
    /// Homogeneous two-field background.
    pub background: MultifieldState<2>,
    /// Complex `2 x 2` mode matrix and its covariant velocity.
    pub modes: TwoFieldModes,
}

/// Time derivative of [`TwoFieldPerturbationState`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldPerturbationDerivative {
    /// Background derivative.
    pub background: MultifieldStateDerivative<2>,
    /// Coordinate time derivative of the mode matrix.
    pub amplitude: [[ModeVelocity; 2]; 2],
    /// Coordinate time derivative of the covariant mode velocity.
    pub velocity: [[ModeAcceleration; 2]; 2],
}

/// Finite change of a background-plus-modes state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldPerturbationDelta {
    /// Background change.
    pub background: MultifieldStateDelta<2>,
    /// Change in mode amplitudes.
    pub amplitude: [[ModeAmplitude; 2]; 2],
    /// Change in covariant mode velocities.
    pub velocity: [[ModeVelocity; 2]; 2],
}

impl Add<TwoFieldPerturbationDelta> for TwoFieldPerturbationState {
    type Output = Self;

    fn add(self, rhs: TwoFieldPerturbationDelta) -> Self::Output {
        Self {
            background: self.background + rhs.background,
            modes: TwoFieldModes {
                rescaled_amplitude: add_matrices(self.modes.rescaled_amplitude, rhs.amplitude),
                rescaled_velocity: add_matrices(self.modes.rescaled_velocity, rhs.velocity),
            },
        }
    }
}

impl Add for TwoFieldPerturbationDerivative {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            background: self.background + rhs.background,
            amplitude: add_matrices(self.amplitude, rhs.amplitude),
            velocity: add_matrices(self.velocity, rhs.velocity),
        }
    }
}

impl Mul<f64> for TwoFieldPerturbationDerivative {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            background: self.background * rhs,
            amplitude: scale_matrix(self.amplitude, rhs),
            velocity: scale_matrix(self.velocity, rhs),
        }
    }
}

impl Mul<Time> for TwoFieldPerturbationDerivative {
    type Output = TwoFieldPerturbationDelta;

    fn mul(self, rhs: Time) -> Self::Output {
        TwoFieldPerturbationDelta {
            background: self.background * rhs,
            amplitude: self.amplitude.map(|row| row.map(|entry| entry * rhs)),
            velocity: self.velocity.map(|row| row.map(|entry| entry * rhs)),
        }
    }
}

fn add_matrices<D: Copy>(
    lhs: [[Quantity<D, Complex64>; 2]; 2],
    rhs: [[Quantity<D, Complex64>; 2]; 2],
) -> [[Quantity<D, Complex64>; 2]; 2] {
    core::array::from_fn(|row| core::array::from_fn(|column| lhs[row][column] + rhs[row][column]))
}

fn scale_matrix<D>(
    matrix: [[Quantity<D, Complex64>; 2]; 2],
    coefficient: f64,
) -> [[Quantity<D, Complex64>; 2]; 2] {
    matrix.map(|row| row.map(|entry| entry * coefficient))
}

/// Failure to construct an adiabatic-vacuum mode state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ModeInitializationError {
    /// The comoving wavenumber must be positive and finite.
    InvalidWavenumber,
    /// The requested mode is not sufficiently far inside the Hubble radius.
    NotDeepInsideHorizon {
        /// Actual `k/(aH)`.
        actual: f64,
        /// Required minimum value.
        required: f64,
    },
    /// The background tangent/normal basis is undefined.
    KinematicBasis(KinematicBasisError),
}

/// A two-field background system plus one comoving perturbation wavenumber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldPerturbationSystem<P, G> {
    /// Background dynamics and geometric data.
    pub background: CovariantMultifieldSystem<P, G, 2>,
    /// Positive comoving wavenumber.
    pub k: Mass,
}

impl<P, G> TwoFieldPerturbationSystem<P, G>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    /// Builds a one-wavenumber perturbation system.
    pub const fn new(background: CovariantMultifieldSystem<P, G, 2>, k: Mass) -> Self {
        Self { background, k }
    }

    /// Constructs a Bunch-Davies state using the instantaneous adiabatic and
    /// entropy frame, with a common phase chosen to be zero.
    pub fn bunch_davies_state(
        &self,
        background: MultifieldState<2>,
        minimum_depth: f64,
    ) -> Result<TwoFieldPerturbationState, ModeInitializationError> {
        if self.k.into_value() <= 0.0 || !self.k.into_value().is_finite() {
            return Err(ModeInitializationError::InvalidWavenumber);
        }
        let hubble = self.background.hubble(&background);
        let depth = (self.k / (background.a * hubble)).into_value();
        if !minimum_depth.is_finite() || minimum_depth <= 0.0 || depth < minimum_depth {
            return Err(ModeInitializationError::NotDeepInsideHorizon {
                actual: depth,
                required: minimum_depth,
            });
        }
        let basis = self
            .background
            .kinematics(&background, None)
            .map_err(ModeInitializationError::KinematicBasis)?
            .basis;
        let normalization = (self.k * 2.0).sqrt();
        // The rescaled vacuum amplitude carries no factor of the scale
        // factor, so it is the same size for every mode: |v| = 1/sqrt(2k).
        let scalar_amplitude: ModeAmplitude<f64> = Dimensionless::new(1.0) / normalization;
        let physical_wavenumber = self.k / background.a;
        // D_t v = a D_t U + a H U = -i (k/a) v: rescaling cancels the Hubble
        // term exactly, leaving a pure oscillation.
        let frequency =
            Quantity::<M1, Complex64>::new(Complex64::new(0.0, -physical_wavenumber.into_value()));
        let frames = [basis.tangent, basis.normal];
        let rescaled_amplitude = core::array::from_fn(|field| {
            core::array::from_fn(|mode| {
                let value = (frames[mode].components()[field] * scalar_amplitude).into_value();
                ModeAmplitude::new(Complex64::new(value, 0.0))
            })
        });
        let rescaled_velocity = rescaled_amplitude.map(|row| row.map(|entry| entry * frequency));

        Ok(TwoFieldPerturbationState {
            background,
            modes: TwoFieldModes {
                rescaled_amplitude,
                rescaled_velocity,
            },
        })
    }

    /// Returns `i a^3 (U^dagger P - P^dagger U)`, whose canonical value is
    /// the identity matrix.
    pub fn wronskian(&self, state: &TwoFieldPerturbationState) -> ComplexMatrix2 {
        let metric = self.background.geometry.metric(&state.background.phi);
        // With v = a U the canonical Wronskian keeps one power of a rather
        // than three, because two are cancelled by the rescaling.
        let scale_factor = state.background.a.into_value();
        core::array::from_fn(|left_mode| {
            core::array::from_fn(|right_mode| {
                let difference = (0..2)
                    .flat_map(|first| (0..2).map(move |second| (first, second)))
                    .fold(Complex64::new(0.0, 0.0), |sum, (first, second)| {
                        let g = *metric.components()[first][second].value();
                        let u_left = state.modes.rescaled_amplitude[first][left_mode].into_value();
                        let p_left = state.modes.rescaled_velocity[first][left_mode].into_value();
                        let u_right =
                            state.modes.rescaled_amplitude[second][right_mode].into_value();
                        let p_right =
                            state.modes.rescaled_velocity[second][right_mode].into_value();
                        sum + g * (u_left.conj() * p_right - p_left.conj() * u_right)
                    });

                Complex64::new(0.0, scale_factor) * difference
            })
        })
    }
}

impl<P, G> System for TwoFieldPerturbationSystem<P, G>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    type IndependentVariable = Time;
    type State = TwoFieldPerturbationState;
    type Delta = TwoFieldPerturbationDelta;
    type Derivative = TwoFieldPerturbationDerivative;

    fn derivative(&self, time: Time, state: &Self::State) -> Self::Derivative {
        let background_derivative = self.background.derivative(time, &state.background);
        let hubble = self.background.hubble(&state.background);
        let connection = self.background.geometry.christoffel(&state.background.phi);
        let mass = self.background.perturbation_mass_matrix(&state.background);
        let physical_wavenumber = self.k / state.background.a;
        let wave_number_squared = physical_wavenumber * physical_wavenumber;
        // Rescaling U by a turns `D_t^2 U + 3H D_t U + (M + k^2/a^2) U = 0`
        // into `D_t^2 v + H D_t v + (M + k^2/a^2 - (2 - epsilon_H) H^2) v = 0`.
        // The subtracted term is `a''/a` written in cosmic time.
        let expansion = hubble * hubble * (2.0 - self.background.epsilon_h(&state.background));
        let amplitude = core::array::from_fn(|upper| {
            core::array::from_fn(|mode| {
                let connection_term = (0..2)
                    .flat_map(|first| (0..2).map(move |second| (first, second)))
                    .fold(
                        ModeVelocity::new(Complex64::new(0.0, 0.0)),
                        |sum, (first, second)| {
                            sum + connection.components()[upper][first][second]
                                * state.background.pi[first]
                                * state.modes.rescaled_amplitude[second][mode]
                        },
                    );
                state.modes.rescaled_velocity[upper][mode] - connection_term
            })
        });
        let velocity = core::array::from_fn(|upper| {
            core::array::from_fn(|mode| {
                let restoring = (0..2).fold(
                    ModeAcceleration::new(Complex64::new(0.0, 0.0)),
                    |sum, lower| {
                        let entry = if upper == lower {
                            mass.components()[upper][lower] + wave_number_squared - expansion
                        } else {
                            mass.components()[upper][lower]
                        };
                        sum + entry * state.modes.rescaled_amplitude[lower][mode]
                    },
                );
                let connection_term = (0..2)
                    .flat_map(|first| (0..2).map(move |second| (first, second)))
                    .fold(
                        ModeAcceleration::new(Complex64::new(0.0, 0.0)),
                        |sum, (first, second)| {
                            sum + connection.components()[upper][first][second]
                                * state.background.pi[first]
                                * state.modes.rescaled_velocity[second][mode]
                        },
                    );

                state.modes.rescaled_velocity[upper][mode] * hubble * -1.0
                    - restoring
                    - connection_term
            })
        });

        TwoFieldPerturbationDerivative {
            background: background_derivative,
            amplitude,
            velocity,
        }
    }
}

impl<P, G> AdaptiveSystem for TwoFieldPerturbationSystem<P, G>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    /// Combines the background norm with the complex mode norms.
    ///
    /// Each mode component contributes the modulus of its whole complex
    /// error, so the real and imaginary parts are never accepted separately.
    fn error_ratio(
        &self,
        current: &Self::State,
        next: &Self::State,
        error: &Self::Delta,
        tolerances: ErrorTolerances,
    ) -> f64 {
        let background = self.background.error_ratio(
            &current.background,
            &next.background,
            &error.background,
            tolerances,
        );

        (0..2)
            .flat_map(|field| (0..2).map(move |mode| (field, mode)))
            .map(|(field, mode)| {
                let amplitude = normalized_mode_error(
                    error.amplitude[field][mode],
                    current.modes.rescaled_amplitude[field][mode],
                    next.modes.rescaled_amplitude[field][mode],
                    ModeAmplitude::<f64>::new(tolerances.absolute),
                    tolerances.relative,
                );
                let velocity = normalized_mode_error(
                    error.velocity[field][mode],
                    current.modes.rescaled_velocity[field][mode],
                    next.modes.rescaled_velocity[field][mode],
                    ModeVelocity::<f64>::new(tolerances.absolute),
                    tolerances.relative,
                );

                amplitude.max(velocity)
            })
            .fold(background, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ErrorTolerances, Field, FieldVelocity, FlatFieldSpace, QuadraticMultifieldPotential,
        ScaleFactor, adaptive_step, rk4_step,
    };

    fn system() -> TwoFieldPerturbationSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace> {
        let background = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(1.5e-6)]),
            FlatFieldSpace,
        );
        TwoFieldPerturbationSystem::new(background, Mass::new(1.0e-2))
    }

    fn background() -> MultifieldState<2> {
        MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(-3.0e-6), FieldVelocity::new(-4.0e-6)],
            ScaleFactor::new(1.0),
        )
    }

    #[test]
    fn complex_quantity_retains_dimension_and_uses_modulus() {
        let value = ModeAmplitude::new(Complex64::new(3.0, 4.0));
        assert_eq!(value.complex_magnitude(), 5.0);
    }

    #[test]
    fn normalized_error_uses_the_modulus_of_the_whole_complex_component() {
        let error = ModeAmplitude::new(Complex64::new(3.0e-9, 4.0e-9));
        let current = ModeAmplitude::new(Complex64::new(1.0, 0.0));
        let next = ModeAmplitude::new(Complex64::new(0.0, 2.0));
        let absolute = ModeAmplitude::<f64>::new(0.0);

        // |error| = 5e-9 and the scale is max(1, 2) = 2.
        let ratio = normalized_mode_error(error, current, next, absolute, 1.0e-9);

        assert!((ratio - 2.5).abs() < 1.0e-12);
    }

    #[test]
    fn a_purely_imaginary_error_is_not_ignored() {
        let current = ModeAmplitude::new(Complex64::new(1.0, 0.0));
        let absolute = ModeAmplitude::<f64>::new(0.0);
        let real_error = ModeAmplitude::new(Complex64::new(1.0e-6, 0.0));
        let imaginary_error = ModeAmplitude::new(Complex64::new(0.0, 1.0e-6));

        let from_real = normalized_mode_error(real_error, current, current, absolute, 1.0e-9);
        let from_imaginary =
            normalized_mode_error(imaginary_error, current, current, absolute, 1.0e-9);

        assert_eq!(from_real, from_imaginary);
        assert!(from_real > 1.0);
    }

    #[test]
    fn a_vanishing_component_falls_back_to_the_absolute_tolerance() {
        let error = ModeAmplitude::new(Complex64::new(1.0e-13, 0.0));
        let zero = ModeAmplitude::new(Complex64::new(0.0, 0.0));
        let absolute = ModeAmplitude::<f64>::new(1.0e-12);

        let ratio = normalized_mode_error(error, zero, zero, absolute, 1.0e-9);

        assert!(ratio < 1.0);
    }

    fn zero_delta() -> TwoFieldPerturbationDelta {
        TwoFieldPerturbationDelta {
            background: MultifieldStateDelta {
                phi: [Field::new(0.0); 2],
                pi: [FieldVelocity::new(0.0); 2],
                a: ScaleFactor::new(0.0),
            },
            amplitude: [[ModeAmplitude::new(Complex64::new(0.0, 0.0)); 2]; 2],
            velocity: [[ModeVelocity::new(Complex64::new(0.0, 0.0)); 2]; 2],
        }
    }

    #[test]
    fn a_purely_imaginary_mode_error_still_drives_the_step_decision() {
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let tolerances = ErrorTolerances::default();

        assert_eq!(
            system.error_ratio(&state, &state, &zero_delta(), tolerances),
            0.0
        );

        let mut imaginary = zero_delta();
        imaginary.amplitude[0][0] = ModeAmplitude::new(Complex64::new(0.0, 1.0));

        assert!(system.error_ratio(&state, &state, &imaginary, tolerances) > 1.0);
    }

    #[test]
    fn an_adaptive_step_accepts_a_small_step_and_rejects_a_huge_one() {
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let tolerances = ErrorTolerances::default();

        let small = adaptive_step(&system, Time::new(0.0), &state, Time::new(1.0), tolerances);
        let huge = adaptive_step(
            &system,
            Time::new(0.0),
            &state,
            Time::new(1.0e6),
            tolerances,
        );

        assert!(small.accepted);
        assert!(!huge.accepted);
        assert!(huge.error_ratio > small.error_ratio);
    }

    #[test]
    fn bunch_davies_modes_have_identity_wronskian() {
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let wronskian = system.wronskian(&state);

        for (row, entries) in wronskian.iter().enumerate() {
            for (column, value) in entries.iter().enumerate() {
                let expected = f64::from(row == column);
                assert!((value.re - expected).abs() < 1.0e-12);
                assert!(value.im.abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn initialization_rejects_a_mode_too_close_to_the_horizon() {
        let system = TwoFieldPerturbationSystem::new(system().background, Mass::new(1.0e-8));

        assert!(matches!(
            system.bunch_davies_state(background(), 100.0),
            Err(ModeInitializationError::NotDeepInsideHorizon { .. })
        ));
    }

    #[test]
    fn generic_integrator_advances_background_and_complex_modes() {
        let system = system();
        let initial = system.bunch_davies_state(background(), 100.0).unwrap();

        let advanced = rk4_step(&system, Time::new(0.0), &initial, Time::new(1.0));

        assert!(advanced.background.a.into_value() > initial.background.a.into_value());
        assert_ne!(
            advanced.modes.rescaled_amplitude,
            initial.modes.rescaled_amplitude
        );
    }

    #[test]
    fn wronskian_is_conserved_by_a_small_rk4_step() {
        let system = system();
        let initial = system.bunch_davies_state(background(), 100.0).unwrap();

        let advanced = rk4_step(&system, Time::new(0.0), &initial, Time::new(1.0e-3));
        let wronskian = system.wronskian(&advanced);

        for (row, entries) in wronskian.iter().enumerate() {
            for (column, value) in entries.iter().enumerate() {
                let expected = f64::from(row == column);
                assert!((value.re - expected).abs() < 1.0e-12);
                assert!(value.im.abs() < 1.0e-12);
            }
        }
    }
}
