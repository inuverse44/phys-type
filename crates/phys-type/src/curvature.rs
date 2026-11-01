//! Projection of phase-space correlations to the curvature perturbation.
//!
//! On a slice where every wavenumber is outside the horizon, the uniform-
//! density curvature perturbation is related to the flat-gauge field
//! perturbations by a nonlinear gauge transformation,
//! `zeta = N_a X^a + N_ab X^a X^b / 2`. This module evaluates the coefficients
//! in the "simple" representation of Dias, Elliston, Frazer, Mulryne, and
//! Seery (arXiv:1410.3491), written covariantly for the canonical phase-space
//! variables `X = (Q, Pi)`, and contracts them with the transported two- and
//! three-point functions.
//!
//! The sign convention is that of that reference and of Maldacena:
//! `zeta = -phi'_I Q^I / (2 H M_pl^2 epsilon)` at linear order, which is
//! `-R` for the comoving curvature `R = (H / dot_sigma) Q_sigma` reported by
//! [`crate::spectra::project`].

use core::fmt;
use core::ops::Add;

use num_complex::Complex64;

use crate::transport::{
    CanonicalThreePointTransportState, CanonicalTwoPointPhaseSpaceCovariance, MomentumTriangle,
};
use crate::{
    ContravariantField, CovariantMultifieldSystem, FieldSpaceGeometry, FourierLeg, M0, MMinus1,
    MMinus2, MMinus3, MMinus6, Mass, MultifieldPotential, MultifieldState, Quantity,
};

/// The quadratic gauge coefficient `N_ab(k_out; k_in1, k_in2)` for one output
/// leg, with covariant field indices of the two inputs.
///
/// The momentum-momentum block vanishes in the simple representation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvatureGaugeQuadratic<const N: usize> {
    /// `N_{Q^I Q^J}`, dimension `M^-2`.
    pub field_field: [[Quantity<MMinus2, f64>; N]; N],
    /// `N_{Q^I Pi^J}`: field on the first input, momentum on the second,
    /// dimension `M^-3`.
    pub field_momentum: [[Quantity<MMinus3, f64>; N]; N],
    /// `N_{Pi^I Q^J}`: momentum on the first input, field on the second,
    /// dimension `M^-3`.
    pub momentum_field: [[Quantity<MMinus3, f64>; N]; N],
}

/// First- and second-order gauge coefficients on one superhorizon slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvatureGaugeCoefficients<const N: usize> {
    /// `N_I` multiplying `Q^I`, dimension `M^-1`. There is no momentum term.
    pub linear: [Quantity<MMinus1, f64>; N],
    /// Entry `r` has output leg `r` and inputs in ascending order among the
    /// remaining legs: `(2, 3)`, `(1, 3)`, and `(1, 2)`.
    pub quadratic: [CurvatureGaugeQuadratic<N>; 3],
}

/// Why a curvature projection could not be formed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CurvatureProjectionError {
    /// The reduced Planck mass was not positive and finite.
    InvalidPlanckMass(f64),
    /// The Hubble rate was not positive and finite.
    InvalidHubble(f64),
    /// The background has no kinetic energy, so `epsilon = 0`.
    StaticBackground,
    /// The requested maximum `k/(aH)` was not finite and positive.
    InvalidSuperhorizonRequirement(f64),
    /// A leg is not far enough outside the horizon for the simple
    /// representation, which drops terms that decay only on superhorizon
    /// scales.
    InsideHorizon {
        /// Canonical leg that failed the requirement.
        leg: FourierLeg,
        /// Its `k/(aH)` on the slice.
        k_over_ah: f64,
        /// The requested maximum.
        maximum: f64,
    },
}

impl fmt::Display for CurvatureProjectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPlanckMass(value) => write!(
                formatter,
                "the curvature projection requires a positive finite Planck mass, got {value}"
            ),
            Self::InvalidHubble(value) => write!(
                formatter,
                "the curvature projection requires a positive finite Hubble rate, got {value}"
            ),
            Self::StaticBackground => {
                write!(formatter, "the curvature projection requires epsilon > 0")
            }
            Self::InvalidSuperhorizonRequirement(value) => write!(
                formatter,
                "a superhorizon slice requires a finite positive maximum k/(aH), got {value}"
            ),
            Self::InsideHorizon {
                leg,
                k_over_ah,
                maximum,
            } => write!(
                formatter,
                "leg {leg:?} has k/(aH)={k_over_ah:e}, above the allowed {maximum:e}"
            ),
        }
    }
}

impl std::error::Error for CurvatureProjectionError {}

/// Evaluates the gauge coefficients on a superhorizon background slice.
///
/// With `epsilon = phi'_I phi'^I / (2 H^2 M_pl^2)`, the coefficients are
///
/// ```text
/// N_I = -phi'_I / (2 H M_pl^2 epsilon)
/// N_{QQ,IJ} = phi'_I phi'_J / (3 H^2 M_pl^4 epsilon)
///     * [-3/2 + 9/(2 epsilon) + 3 V_,K phi'^K / (4 epsilon^2 M_pl^2 H^3)]
/// N_{Q Pi,IJ}(k; k_a, k_b) = [3 phi'_I phi'_J / (H epsilon M_pl^2)
///     - 3 H (k_a.k_b + k_a^2) G_IJ / k^2] / (6 H^2 M_pl^2 epsilon)
/// N_{Pi Q,IJ}(k; k_a, k_b) = the same with k_a^2 replaced by k_b^2
/// ```
///
/// The cross terms carry `1/6` rather than the `1/3` printed in Eq. (7.4b) of
/// Dias, Frazer, Mulryne, and Seery (arXiv:1609.00379), because the symmetric
/// quadratic form `N_ab X^a X^b / 2` counts the `Q Pi` and `Pi Q` orderings
/// separately; see `docs/literature-discrepancies.md`, entry 1. Every leg must satisfy `k_r/(aH) <= maximum_k_over_ah`.
pub fn curvature_gauge_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
    maximum_k_over_ah: f64,
) -> Result<CurvatureGaugeCoefficients<N>, CurvatureProjectionError>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    let planck_value = system.m_pl.into_value();
    if !planck_value.is_finite() || planck_value <= 0.0 {
        return Err(CurvatureProjectionError::InvalidPlanckMass(planck_value));
    }
    let hubble = system.hubble(state);
    let hubble_value = hubble.into_value();
    if !hubble_value.is_finite() || hubble_value <= 0.0 {
        return Err(CurvatureProjectionError::InvalidHubble(hubble_value));
    }
    if !maximum_k_over_ah.is_finite() || maximum_k_over_ah <= 0.0 {
        return Err(CurvatureProjectionError::InvalidSuperhorizonRequirement(
            maximum_k_over_ah,
        ));
    }
    let wavenumbers = triangle.canonical_wavenumbers();
    for (index, wavenumber) in wavenumbers.into_iter().enumerate() {
        let k_over_ah = (Mass::new(wavenumber) / (state.a * hubble)).into_value();
        if k_over_ah > maximum_k_over_ah {
            return Err(CurvatureProjectionError::InsideHorizon {
                leg: [FourierLeg::One, FourierLeg::Two, FourierLeg::Three][index],
                k_over_ah,
                maximum: maximum_k_over_ah,
            });
        }
    }

    let metric = system.geometry.metric(&state.phi);
    let velocity = ContravariantField::new(state.pi);
    let lowered = metric.lower(&velocity);
    let velocity_squared = lowered.contract(&velocity);
    if velocity_squared.into_value() <= 0.0 {
        return Err(CurvatureProjectionError::StaticBackground);
    }
    let planck = Mass::new(planck_value);
    let planck_squared = planck * planck;
    let epsilon: Quantity<M0, f64> = velocity_squared / (hubble * hubble * planck_squared * 2.0);
    let epsilon_value = epsilon.into_value();
    let gradient = system.potential.gradient(state.phi);
    let slope = (0..N)
        .map(|i| gradient[i] * state.pi[i])
        .fold(Quantity::new(0.0), Add::add);
    let bracket = -1.5
        + 4.5 / epsilon_value
        + 0.75 / (epsilon_value * epsilon_value)
            * (slope / (planck_squared * hubble * hubble * hubble)).into_value();
    let v = lowered.components();

    let linear =
        core::array::from_fn(|i| -(v[i] / (hubble * planck_squared * 2.0 * epsilon_value)));
    let quadratic = [FourierLeg::One, FourierLeg::Two, FourierLeg::Three].map(|output| {
        let (first, second) = match output {
            FourierLeg::One => (FourierLeg::Two, FourierLeg::Three),
            FourierLeg::Two => (FourierLeg::One, FourierLeg::Three),
            FourierLeg::Three => (FourierLeg::One, FourierLeg::Two),
        };
        let k_squared = triangle.wavevector_dot(output, output);
        let inputs_dot = triangle.wavevector_dot(first, second);
        let first_weight =
            ((inputs_dot + triangle.wavevector_dot(first, first)) / k_squared).into_value();
        let second_weight =
            ((inputs_dot + triangle.wavevector_dot(second, second)) / k_squared).into_value();
        let cross = |i: usize, j: usize, weight: f64| -> Quantity<MMinus3, f64> {
            (v[i] * v[j] * 3.0 / (hubble * planck_squared * epsilon_value)
                - hubble * metric.components()[i][j].into_value() * 3.0 * weight)
                / (hubble * hubble * planck_squared * 6.0 * epsilon_value)
        };
        CurvatureGaugeQuadratic {
            field_field: core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    v[i] * v[j] * bracket
                        / (hubble * hubble * planck_squared * planck_squared * 3.0 * epsilon_value)
                })
            }),
            field_momentum: core::array::from_fn(|i| {
                core::array::from_fn(|j| cross(i, j, first_weight))
            }),
            momentum_field: core::array::from_fn(|i| {
                core::array::from_fn(|j| cross(i, j, second_weight))
            }),
        }
    });

    Ok(CurvatureGaugeCoefficients { linear, quadratic })
}

/// Stripped curvature correlators on one slice.
///
/// `power[r]` is `P_zeta(k_r)` with `<zeta zeta> = (2 pi)^3 delta P`, and
/// `bispectrum` is `B_zeta` with `<zeta zeta zeta> = (2 pi)^3 delta B`, both
/// in comoving Fourier normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvatureCorrelators {
    /// Canonically ordered comoving wavenumbers.
    pub wavenumbers: [Mass; 3],
    /// `P_zeta(k_r)`, dimension `M^-3`.
    pub power: [Quantity<MMinus3, f64>; 3],
    /// Real part of `B_zeta`, dimension `M^-6`.
    pub bispectrum: Quantity<MMinus6, f64>,
    /// `|Im B_zeta| / |Re B_zeta|` from the operator-ordered contraction,
    /// which measures residual commutator contributions on the slice.
    pub bispectrum_imaginary_ratio: f64,
    /// The intrinsic part `N_a N_b N_c Alpha^abc`, dimension `M^-6`.
    pub intrinsic: Quantity<MMinus6, f64>,
}

impl CurvatureCorrelators {
    /// `k^3 P_zeta(k) / (2 pi^2)` for each leg.
    pub fn dimensionless_power(&self) -> [f64; 3] {
        core::array::from_fn(|r| {
            let k = self.wavenumbers[r];
            (k * k * k * self.power[r]).into_value() / (2.0 * core::f64::consts::PI.powi(2))
        })
    }

    /// The reduced bispectrum in the Maldacena/Planck convention,
    /// `B_zeta = (6/5) f_NL [P(k_1) P(k_2) + P(k_1) P(k_3) + P(k_2) P(k_3)]`,
    /// with the sign of `zeta` documented in the module.
    pub fn reduced_f_nl(&self) -> f64 {
        let [p1, p2, p3] = self.power;
        let pairs = p1 * p2 + p1 * p3 + p2 * p3;
        5.0 / 6.0 * (self.bispectrum / pairs).into_value()
    }
}

/// `N_I Sigma` with the external field slot projected on `zeta`, in operator
/// order: an external leg to the left of the quadratic insertion contracts as
/// `Sigma^{a c}`, one to the right as `Sigma^{c a}`.
fn projected_field_internal<const N: usize>(
    linear: &[Quantity<MMinus1, f64>; N],
    covariance: &CanonicalTwoPointPhaseSpaceCovariance<N>,
    external_is_left: bool,
    internal: usize,
) -> Quantity<MMinus2, Complex64> {
    (0..N)
        .map(|i| {
            let block = covariance.qq.components();
            let component = if external_is_left {
                block[i][internal]
            } else {
                block[internal][i]
            };
            linear[i] * component
        })
        .fold(Quantity::new(Complex64::new(0.0, 0.0)), Add::add)
}

/// As [`projected_field_internal`], with a momentum slot on the insertion.
fn projected_momentum_internal<const N: usize>(
    linear: &[Quantity<MMinus1, f64>; N],
    covariance: &CanonicalTwoPointPhaseSpaceCovariance<N>,
    external_is_left: bool,
    internal: usize,
) -> Quantity<MMinus1, Complex64> {
    (0..N)
        .map(|i| {
            let component = if external_is_left {
                covariance.qp.components()[i][internal]
            } else {
                covariance.pq.components()[internal][i]
            };
            linear[i] * component
        })
        .fold(Quantity::new(Complex64::new(0.0, 0.0)), Add::add)
}

/// Projects an evolved canonical transport state to `P_zeta` and `B_zeta`.
///
/// `background` must be the slice on which `state` is evaluated, and every
/// leg must satisfy `k/(aH) <= maximum_k_over_ah`. Because the linear gauge
/// coefficient has no momentum part, the canonical and velocity bases give
/// the same result: `Q` is common to both and the quadratic insertion
/// differs only at cubic order.
pub fn curvature_correlators<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    background: &MultifieldState<N>,
    triangle: MomentumTriangle,
    state: &CanonicalThreePointTransportState<N>,
    maximum_k_over_ah: f64,
) -> Result<CurvatureCorrelators, CurvatureProjectionError>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    let coefficients =
        curvature_gauge_coefficients(system, background, triangle, maximum_k_over_ah)?;
    let n = &coefficients.linear;
    let power = state.covariances.map(|covariance| {
        (0..N * N)
            .map(|m| {
                let [i, j] = [m / N, m % N];
                n[i] * n[j]
                    * Quantity::<MMinus1, f64>::new(covariance.qq.components()[i][j].value().re)
            })
            .fold(Quantity::new(0.0), Add::add)
    });

    let intrinsic: Quantity<MMinus6, Complex64> = (0..N * N * N)
        .map(|m| {
            let [i, j, k] = [m / (N * N), (m / N) % N, m % N];
            n[i] * n[j] * n[k] * state.alpha.qqq.components()[i][j][k]
        })
        .fold(Quantity::new(Complex64::new(0.0, 0.0)), Add::add);
    let gauge: Quantity<MMinus6, Complex64> = (0..3)
        .map(|output| {
            let [first, second] = match output {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            let vertex = &coefficients.quadratic[output];
            let field = |leg: usize, index: usize| {
                projected_field_internal(n, &state.covariances[leg], leg < output, index)
            };
            let momentum = |leg: usize, index: usize| {
                projected_momentum_internal(n, &state.covariances[leg], leg < output, index)
            };
            (0..N * N)
                .map(|m| {
                    let [c, d] = [m / N, m % N];
                    vertex.field_field[c][d] * field(first, c) * field(second, d)
                        + vertex.field_momentum[c][d] * field(first, c) * momentum(second, d)
                        + vertex.momentum_field[c][d] * momentum(first, c) * field(second, d)
                })
                .fold(Quantity::new(Complex64::new(0.0, 0.0)), Add::add)
        })
        .fold(Quantity::new(Complex64::new(0.0, 0.0)), Add::add);
    let total = (intrinsic + gauge).into_value();

    Ok(CurvatureCorrelators {
        wavenumbers: triangle.canonical_wavenumbers().map(Mass::new),
        power,
        bispectrum: Quantity::new(total.re),
        bispectrum_imaginary_ratio: total.im.abs() / total.re.abs().max(f64::MIN_POSITIVE),
        intrinsic: Quantity::new(intrinsic.into_value().re),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Field, FieldVelocity, FlatFieldSpace, QuadraticMultifieldPotential, ScaleFactor};

    fn single_field() -> (
        CovariantMultifieldSystem<QuadraticMultifieldPotential<1>, FlatFieldSpace, 1>,
        MultifieldState<1>,
    ) {
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(15.0)],
            [FieldVelocity::new(-8.0e-8)],
            ScaleFactor::new(1.0),
        );
        (system, state)
    }

    #[test]
    fn single_field_linear_coefficient_is_minus_h_over_phi_dot() {
        let (system, state) = single_field();
        let hubble = system.hubble(&state).into_value();
        let k = 1.0e-3 * hubble;
        let triangle = MomentumTriangle::try_new(k, k, k).unwrap();
        let coefficients = curvature_gauge_coefficients(&system, &state, triangle, 0.1).unwrap();
        let expected = -hubble / state.pi[0].into_value();
        let actual = coefficients.linear[0].into_value();

        assert!((actual - expected).abs() <= 1.0e-12 * expected.abs());
    }

    #[test]
    fn cross_coefficients_attach_the_gradient_term_to_the_field_input() {
        let (system, state) = single_field();
        let hubble = system.hubble(&state).into_value();
        let scale = 1.0e-3 * hubble;
        let triangle = MomentumTriangle::try_new(4.0 * scale, 3.0 * scale, 2.0 * scale).unwrap();
        let coefficients = curvature_gauge_coefficients(&system, &state, triangle, 0.1).unwrap();
        let velocity = state.pi[0].into_value();
        let epsilon = velocity * velocity / (2.0 * hubble * hubble);
        let [k1, k2, k3] = triangle.canonical_wavenumbers();
        let dot12 = (k3 * k3 - k1 * k1 - k2 * k2) / 2.0;
        let expected = |field_leg: f64| {
            (3.0 * velocity * velocity / (hubble * epsilon)
                - 3.0 * hubble * (dot12 + field_leg * field_leg) / (k3 * k3))
                / (6.0 * hubble * hubble * epsilon)
        };
        // Output leg 3 has inputs (1, 2).
        let vertex = &coefficients.quadratic[2];
        let field_momentum = vertex.field_momentum[0][0].into_value();
        let momentum_field = vertex.momentum_field[0][0].into_value();

        assert!((field_momentum - expected(k1)).abs() <= 1.0e-12 * field_momentum.abs());
        assert!((momentum_field - expected(k2)).abs() <= 1.0e-12 * momentum_field.abs());
        assert_ne!(field_momentum, momentum_field);
    }

    #[test]
    fn the_projection_requires_every_leg_outside_the_horizon() {
        let (system, state) = single_field();
        let hubble = system.hubble(&state).into_value();
        let triangle =
            MomentumTriangle::try_new(0.2 * hubble, 0.15 * hubble, 0.1 * hubble).unwrap();

        assert!(matches!(
            curvature_gauge_coefficients(&system, &state, triangle, 0.1),
            Err(CurvatureProjectionError::InsideHorizon {
                leg: FourierLeg::One,
                ..
            })
        ));
        assert!(curvature_gauge_coefficients(&system, &state, triangle, 0.5).is_ok());
        assert!(matches!(
            curvature_gauge_coefficients(&system, &state, triangle, 0.0),
            Err(CurvatureProjectionError::InvalidSuperhorizonRequirement(_))
        ));
    }
}
