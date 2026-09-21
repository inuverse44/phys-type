//! Linear two-field perturbations evolved with their background.

use core::ops::{Add, Mul};

use crate::integrate::System;
use crate::{
    CovariantMultifieldSystem, CurvedFieldSpaceGeometry, Dimensionless, KinematicBasisError,
    M3Over2, MHalf, MMinusHalf, Mass, MultifieldState, MultifieldStateDelta,
    MultifieldStateDerivative, Quantity, Time, TwiceDifferentiableMultifieldPotential,
};

/// Fourier mode amplitude `U^I_A`, with dimension fixed by
/// `P_Q = k^3 |U|^2 / (2 pi^2)`.
pub type ModeAmplitude = Quantity<MMinusHalf, f64>;
/// Covariant time derivative `D_t U^I_A`.
pub type ModeVelocity = Quantity<MHalf, f64>;
/// Covariant second derivative of a mode amplitude.
pub type ModeAcceleration = Quantity<M3Over2, f64>;

/// Real and imaginary parts of the two independent two-field mode functions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldModes {
    /// Real part of `U^I_A`; first index is field space and second labels the
    /// independent vacuum oscillator.
    pub real: [[ModeAmplitude; 2]; 2],
    /// Imaginary part of `U^I_A`.
    pub imaginary: [[ModeAmplitude; 2]; 2],
    /// Real part of `D_t U^I_A`.
    pub velocity_real: [[ModeVelocity; 2]; 2],
    /// Imaginary part of `D_t U^I_A`.
    pub velocity_imaginary: [[ModeVelocity; 2]; 2],
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
    /// Coordinate time derivative of the real mode matrix.
    pub real: [[ModeVelocity; 2]; 2],
    /// Coordinate time derivative of the imaginary mode matrix.
    pub imaginary: [[ModeVelocity; 2]; 2],
    /// Coordinate time derivative of the real covariant mode velocity.
    pub velocity_real: [[ModeAcceleration; 2]; 2],
    /// Coordinate time derivative of the imaginary covariant mode velocity.
    pub velocity_imaginary: [[ModeAcceleration; 2]; 2],
}

/// Finite change of a background-plus-modes state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldPerturbationDelta {
    /// Background change.
    pub background: MultifieldStateDelta<2>,
    /// Change in real mode amplitudes.
    pub real: [[ModeAmplitude; 2]; 2],
    /// Change in imaginary mode amplitudes.
    pub imaginary: [[ModeAmplitude; 2]; 2],
    /// Change in real covariant mode velocities.
    pub velocity_real: [[ModeVelocity; 2]; 2],
    /// Change in imaginary covariant mode velocities.
    pub velocity_imaginary: [[ModeVelocity; 2]; 2],
}

impl Add<TwoFieldPerturbationDelta> for TwoFieldPerturbationState {
    type Output = Self;

    fn add(self, rhs: TwoFieldPerturbationDelta) -> Self::Output {
        Self {
            background: self.background + rhs.background,
            modes: TwoFieldModes {
                real: add_matrices(self.modes.real, rhs.real),
                imaginary: add_matrices(self.modes.imaginary, rhs.imaginary),
                velocity_real: add_matrices(self.modes.velocity_real, rhs.velocity_real),
                velocity_imaginary: add_matrices(
                    self.modes.velocity_imaginary,
                    rhs.velocity_imaginary,
                ),
            },
        }
    }
}

impl Add for TwoFieldPerturbationDerivative {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            background: self.background + rhs.background,
            real: add_matrices(self.real, rhs.real),
            imaginary: add_matrices(self.imaginary, rhs.imaginary),
            velocity_real: add_matrices(self.velocity_real, rhs.velocity_real),
            velocity_imaginary: add_matrices(self.velocity_imaginary, rhs.velocity_imaginary),
        }
    }
}

impl Mul<f64> for TwoFieldPerturbationDerivative {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            background: self.background * rhs,
            real: scale_matrix(self.real, rhs),
            imaginary: scale_matrix(self.imaginary, rhs),
            velocity_real: scale_matrix(self.velocity_real, rhs),
            velocity_imaginary: scale_matrix(self.velocity_imaginary, rhs),
        }
    }
}

impl Mul<Time> for TwoFieldPerturbationDerivative {
    type Output = TwoFieldPerturbationDelta;

    fn mul(self, rhs: Time) -> Self::Output {
        TwoFieldPerturbationDelta {
            background: self.background * rhs,
            real: self.real.map(|row| row.map(|entry| entry * rhs)),
            imaginary: self.imaginary.map(|row| row.map(|entry| entry * rhs)),
            velocity_real: self.velocity_real.map(|row| row.map(|entry| entry * rhs)),
            velocity_imaginary: self
                .velocity_imaginary
                .map(|row| row.map(|entry| entry * rhs)),
        }
    }
}

fn add_matrices<D: Copy>(
    lhs: [[Quantity<D, f64>; 2]; 2],
    rhs: [[Quantity<D, f64>; 2]; 2],
) -> [[Quantity<D, f64>; 2]; 2] {
    core::array::from_fn(|row| core::array::from_fn(|column| lhs[row][column] + rhs[row][column]))
}

fn scale_matrix<D>(
    matrix: [[Quantity<D, f64>; 2]; 2],
    coefficient: f64,
) -> [[Quantity<D, f64>; 2]; 2] {
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
        let amplitude: ModeAmplitude = Dimensionless::new(1.0) / (background.a * normalization);
        let physical_wavenumber = self.k / background.a;
        let frames = [basis.tangent, basis.normal];
        let real = core::array::from_fn(|field| {
            core::array::from_fn(|mode| frames[mode].components()[field] * amplitude)
        });
        let imaginary = [[ModeAmplitude::new(0.0); 2]; 2];
        let velocity_real = real.map(|row| row.map(|entry| entry * hubble * -1.0));
        let velocity_imaginary =
            real.map(|row| row.map(|entry| entry * physical_wavenumber * -1.0));

        Ok(TwoFieldPerturbationState {
            background,
            modes: TwoFieldModes {
                real,
                imaginary,
                velocity_real,
                velocity_imaginary,
            },
        })
    }

    /// Returns `i a^3 (U^dagger P - P^dagger U)`, whose canonical value is
    /// the identity matrix.
    pub fn wronskian(&self, state: &TwoFieldPerturbationState) -> ComplexMatrix2 {
        let metric = self.background.geometry.metric(&state.background.phi);
        let mut result = ComplexMatrix2::default();
        let a_cubed = state.background.a.into_value().powi(3);

        for left_mode in 0..2 {
            for right_mode in 0..2 {
                let mut difference_real = 0.0;
                let mut difference_imaginary = 0.0;
                for first in 0..2 {
                    for second in 0..2 {
                        let g = *metric.components()[first][second].value();
                        let ur_a = state.modes.real[first][left_mode].into_value();
                        let ui_a = state.modes.imaginary[first][left_mode].into_value();
                        let pr_a = state.modes.velocity_real[first][left_mode].into_value();
                        let pi_a = state.modes.velocity_imaginary[first][left_mode].into_value();
                        let ur_b = state.modes.real[second][right_mode].into_value();
                        let ui_b = state.modes.imaginary[second][right_mode].into_value();
                        let pr_b = state.modes.velocity_real[second][right_mode].into_value();
                        let pi_b = state.modes.velocity_imaginary[second][right_mode].into_value();

                        difference_real +=
                            g * (ur_a * pr_b + ui_a * pi_b - pr_a * ur_b - pi_a * ui_b);
                        difference_imaginary +=
                            g * (ur_a * pi_b - ui_a * pr_b - pr_a * ui_b + pi_a * ur_b);
                    }
                }
                result.real[left_mode][right_mode] = -a_cubed * difference_imaginary;
                result.imaginary[left_mode][right_mode] = a_cubed * difference_real;
            }
        }

        result
    }
}

/// A raw complex `2 x 2` diagnostic matrix.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ComplexMatrix2 {
    /// Real components.
    pub real: [[f64; 2]; 2],
    /// Imaginary components.
    pub imaginary: [[f64; 2]; 2],
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

        let mode_derivative = |amplitude: &[[ModeAmplitude; 2]; 2],
                               velocity: &[[ModeVelocity; 2]; 2]| {
            let coordinate_velocity = core::array::from_fn(|upper| {
                core::array::from_fn(|mode| {
                    let connection_term = (0..2)
                        .flat_map(|first| (0..2).map(move |second| (first, second)))
                        .fold(ModeVelocity::new(0.0), |sum, (first, second)| {
                            sum + connection.components()[upper][first][second]
                                * state.background.pi[first]
                                * amplitude[second][mode]
                        });
                    velocity[upper][mode] - connection_term
                })
            });
            let coordinate_acceleration = core::array::from_fn(|upper| {
                core::array::from_fn(|mode| {
                    let restoring = (0..2).fold(ModeAcceleration::new(0.0), |sum, lower| {
                        let entry = if upper == lower {
                            mass.components()[upper][lower] + wave_number_squared
                        } else {
                            mass.components()[upper][lower]
                        };
                        sum + entry * amplitude[lower][mode]
                    });
                    let connection_term = (0..2)
                        .flat_map(|first| (0..2).map(move |second| (first, second)))
                        .fold(ModeAcceleration::new(0.0), |sum, (first, second)| {
                            sum + connection.components()[upper][first][second]
                                * state.background.pi[first]
                                * velocity[second][mode]
                        });

                    velocity[upper][mode] * hubble * -3.0 - restoring - connection_term
                })
            });

            (coordinate_velocity, coordinate_acceleration)
        };
        let (real, velocity_real) = mode_derivative(&state.modes.real, &state.modes.velocity_real);
        let (imaginary, velocity_imaginary) =
            mode_derivative(&state.modes.imaginary, &state.modes.velocity_imaginary);

        TwoFieldPerturbationDerivative {
            background: background_derivative,
            real,
            imaginary,
            velocity_real,
            velocity_imaginary,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Field, FieldVelocity, FlatFieldSpace, QuadraticMultifieldPotential, ScaleFactor, rk4_step,
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
    fn bunch_davies_modes_have_identity_wronskian() {
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let wronskian = system.wronskian(&state);

        for row in 0..2 {
            for column in 0..2 {
                let expected = f64::from(row == column);
                assert!((wronskian.real[row][column] - expected).abs() < 1.0e-12);
                assert!(wronskian.imaginary[row][column].abs() < 1.0e-12);
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
        assert_ne!(advanced.modes.real, initial.modes.real);
        assert_ne!(advanced.modes.imaginary, initial.modes.imaginary);
    }

    #[test]
    fn wronskian_is_conserved_by_a_small_rk4_step() {
        let system = system();
        let initial = system.bunch_davies_state(background(), 100.0).unwrap();

        let advanced = rk4_step(&system, Time::new(0.0), &initial, Time::new(1.0e-3));
        let wronskian = system.wronskian(&advanced);

        for row in 0..2 {
            for column in 0..2 {
                let expected = f64::from(row == column);
                assert!((wronskian.real[row][column] - expected).abs() < 1.0e-12);
                assert!(wronskian.imaginary[row][column].abs() < 1.0e-12);
            }
        }
    }
}
