//! Einstein-frame geometry and potentials of nonminimally assisted inflation.
//!
//! Hyun, Kim, Kodama, Park, and Takahashi,
//! [*Nonminimally Assisted Inflation: A General Analysis*](https://arxiv.org/abs/2302.05866)
//! (arXiv:2302.05866v2), consider an inflaton `phi` and an assistant field `s`
//! with the nonminimal coupling `f(s) = (M_pl^2 / 2)[1 + xi_2 (s/M_pl)^2]`,
//! Eq. (2.5). After the Weyl rescaling the Einstein-frame kinetic terms are
//! Eqs. (2.6)--(2.8),
//!
//! ```text
//! G = diag(K_1(s), K_2(s)),   K_1 = M_pl^2 / (2 f),
//! K_2 = M_pl^2 / (2 f) + kappa (3 M_pl^2 / 2) (f_,s / f)^2,
//! ```
//!
//! with `kappa = 1` for the metric and `kappa = 0` for the Palatini
//! formulation, and the potential is `V = M_pl^4 V_J(phi) / (4 f^2)`,
//! Eq. (2.9). The coordinates here are `(phi^1, phi^2) = (phi, s)`; no field
//! is canonically normalized, so the geometry supplies every connection and
//! curvature term to the covariant transport.
//!
//! Because the metric depends only on `s` and the field space is
//! two-dimensional, the curvature is fixed by the Gaussian curvature `K(s)`:
//! `R^I_JKL = K (delta^I_K G_JL - delta^I_L G_JK)` in the crate convention, and
//! `nabla_M R^I_JKL = (partial_M K)(delta^I_K G_JL - delta^I_L G_JK)`.

use core::fmt;
use core::ops::{Add, Div, Mul};

use crate::{
    ChristoffelDerivatives, ChristoffelSymbols, CovariantRiemannDerivative,
    CurvatureDerivativeGeometry, CurvedFieldSpaceGeometry, EnergyDensity, Field, FieldMetric,
    FieldSpaceGeometry, InverseFieldMetric, MMinus1, MMinus2, MMinus3, Mass,
    MultifieldPotentialValue, Quantity, RiemannTensor, ThirdOrderFieldSpaceGeometry, ThirdOrderJet,
    Transcendental,
};

/// Invalid parameters of a nonminimally assisted model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NonminimalParameterError {
    /// The reduced Planck mass was not positive and finite.
    InvalidPlanckMass(f64),
    /// The nonminimal coupling `xi_2` was not finite and non-negative.
    InvalidCoupling(f64),
    /// `kappa` was not 0 (Palatini) or 1 (metric).
    InvalidFormulation(f64),
}

impl fmt::Display for NonminimalParameterError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPlanckMass(value) => {
                write!(
                    formatter,
                    "the Planck mass must be positive and finite, got {value}"
                )
            }
            Self::InvalidCoupling(value) => {
                write!(
                    formatter,
                    "xi_2 must be finite and non-negative, got {value}"
                )
            }
            Self::InvalidFormulation(value) => write!(
                formatter,
                "kappa must be 0 (Palatini) or 1 (metric formulation), got {value}"
            ),
        }
    }
}

impl std::error::Error for NonminimalParameterError {}

/// The Einstein-frame field space of nonminimally assisted inflation in the
/// coordinates `(phi, s)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonminimalAssistedFieldSpace {
    coupling: f64,
    kappa: f64,
    planck_mass: Mass,
}

impl NonminimalAssistedFieldSpace {
    /// Builds the geometry for `xi_2 = coupling` and `kappa` in `{0, 1}`.
    pub fn new(
        coupling: f64,
        kappa: f64,
        planck_mass: Mass,
    ) -> Result<Self, NonminimalParameterError> {
        let planck = planck_mass.into_value();
        if !planck.is_finite() || planck <= 0.0 {
            return Err(NonminimalParameterError::InvalidPlanckMass(planck));
        }
        if !coupling.is_finite() || coupling < 0.0 {
            return Err(NonminimalParameterError::InvalidCoupling(coupling));
        }
        if kappa != 0.0 && kappa != 1.0 {
            return Err(NonminimalParameterError::InvalidFormulation(kappa));
        }
        Ok(Self {
            coupling,
            kappa,
            planck_mass,
        })
    }

    /// The nonminimal coupling `xi_2`.
    pub const fn coupling(&self) -> f64 {
        self.coupling
    }

    /// `kappa`: 1 for the metric and 0 for the Palatini formulation.
    pub const fn kappa(&self) -> f64 {
        self.kappa
    }

    /// Values and first three `s` derivatives of `K_1` and `K_2`, in reduced
    /// Planck units (`d/ds` carries `M^-1`).
    fn metric_functions(&self, s: Field) -> ([f64; 4], [f64; 4]) {
        let planck = self.planck_mass.into_value();
        let t = ThirdOrderJet::<1>::variable((s / self.planck_mass).into_value(), 0);
        let one = ThirdOrderJet::<1>::constant(1.0);
        let first = one / (one + t * t * self.coupling);
        let second =
            first + t * t * first * first * (6.0 * self.kappa * self.coupling * self.coupling);
        let derivatives = |jet: ThirdOrderJet<1>| {
            [
                jet.value,
                jet.gradient[0] / planck,
                jet.hessian[0][0] / (planck * planck),
                jet.third[0][0][0] / (planck * planck * planck),
            ]
        };
        (derivatives(first), derivatives(second))
    }

    /// Gaussian curvature `K` and `dK/ds` of `K_1 dphi^2 + K_2 ds^2`.
    fn gaussian_curvature(&self, s: Field) -> (f64, f64) {
        let ([a0, a1, a2, a3], [b0, b1, b2, _]) = self.metric_functions(s);
        let curvature =
            (a0 * a1 * b1 - 2.0 * a0 * a2 * b0 + a1 * a1 * b0) / (4.0 * a0 * a0 * b0 * b0);
        let derivative = (a0 * a0 * a1 * b0 * b2 - 2.0 * a0 * a0 * a1 * b1 * b1
            + 3.0 * a0 * a0 * a2 * b0 * b1
            - 2.0 * a0 * a0 * a3 * b0 * b0
            - 2.0 * a0 * a1 * a1 * b0 * b1
            + 4.0 * a0 * a1 * a2 * b0 * b0
            - 2.0 * a1 * a1 * a1 * b0 * b0)
            / (4.0 * a0.powi(3) * b0.powi(3));
        (curvature, derivative)
    }

    /// `delta^I_K G_JL - delta^I_L G_JK` at `position`.
    fn curvature_structure(
        &self,
        position: &[Field; 2],
        i: usize,
        j: usize,
        k: usize,
        l: usize,
    ) -> f64 {
        let metric = self.metric(position);
        let g = |a: usize, b: usize| metric.components()[a][b].into_value();
        f64::from(i == k) * g(j, l) - f64::from(i == l) * g(j, k)
    }
}

impl FieldSpaceGeometry<2> for NonminimalAssistedFieldSpace {
    fn metric(&self, position: &[Field; 2]) -> FieldMetric<2> {
        let ([a, ..], [b, ..]) = self.metric_functions(position[1]);
        FieldMetric::new([
            [Quantity::new(a), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(b)],
        ])
    }

    fn inverse_metric(&self, position: &[Field; 2]) -> InverseFieldMetric<2> {
        let ([a, ..], [b, ..]) = self.metric_functions(position[1]);
        InverseFieldMetric::new([
            [Quantity::new(1.0 / a), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(1.0 / b)],
        ])
    }

    fn christoffel(&self, position: &[Field; 2]) -> ChristoffelSymbols<2> {
        let ([a0, a1, ..], [b0, b1, ..]) = self.metric_functions(position[1]);
        let mut components = [[[Quantity::<MMinus1, f64>::new(0.0); 2]; 2]; 2];
        // Gamma^phi_{phi s} = K_1' / (2 K_1)
        components[0][0][1] = Quantity::new(a1 / (2.0 * a0));
        components[0][1][0] = Quantity::new(a1 / (2.0 * a0));
        // Gamma^s_{phi phi} = -K_1' / (2 K_2), Gamma^s_{s s} = K_2' / (2 K_2)
        components[1][0][0] = Quantity::new(-a1 / (2.0 * b0));
        components[1][1][1] = Quantity::new(b1 / (2.0 * b0));
        ChristoffelSymbols::new(components)
    }
}

impl CurvedFieldSpaceGeometry<2> for NonminimalAssistedFieldSpace {
    fn riemann(&self, position: &[Field; 2]) -> RiemannTensor<2> {
        let (curvature, _) = self.gaussian_curvature(position[1]);
        RiemannTensor::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|l| {
                        Quantity::<MMinus2, f64>::new(
                            curvature * self.curvature_structure(position, i, j, k, l),
                        )
                    })
                })
            })
        }))
    }
}

impl ThirdOrderFieldSpaceGeometry<2> for NonminimalAssistedFieldSpace {
    fn christoffel_derivative(&self, position: &[Field; 2]) -> ChristoffelDerivatives<2> {
        let ([a0, a1, a2, _], [b0, b1, b2, _]) = self.metric_functions(position[1]);
        let mut components = [[[[Quantity::<MMinus2, f64>::new(0.0); 2]; 2]; 2]; 2];
        // Only the s derivative (last index 1) is nonzero.
        let phi_s = (a2 * a0 - a1 * a1) / (2.0 * a0 * a0);
        components[0][0][1][1] = Quantity::new(phi_s);
        components[0][1][0][1] = Quantity::new(phi_s);
        components[1][0][0][1] = Quantity::new(-(a2 * b0 - a1 * b1) / (2.0 * b0 * b0));
        components[1][1][1][1] = Quantity::new((b2 * b0 - b1 * b1) / (2.0 * b0 * b0));
        ChristoffelDerivatives::new(components)
    }
}

impl CurvatureDerivativeGeometry<2> for NonminimalAssistedFieldSpace {
    fn covariant_riemann_derivative(&self, position: &[Field; 2]) -> CovariantRiemannDerivative<2> {
        // Metric compatibility leaves only the derivative of K.
        let (_, derivative) = self.gaussian_curvature(position[1]);
        CovariantRiemannDerivative::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|l| {
                        core::array::from_fn(|m| {
                            Quantity::<MMinus3, f64>::new(
                                f64::from(m == 1)
                                    * derivative
                                    * self.curvature_structure(position, i, j, k, l),
                            )
                        })
                    })
                })
            })
        }))
    }
}

/// Loop inflation, `V_J = Lambda^4 [1 + alpha ln(phi/M_pl)]` (Eq. (4.1) of the
/// reference), in the Einstein frame of the assisted model,
/// `V = V_J / (1 + xi_2 s^2 / M_pl^2)^2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonminimalLoopInflationPotential {
    /// `Lambda^4`, an overall energy density.
    pub scale: EnergyDensity,
    /// The one-loop coefficient `alpha`.
    pub alpha: f64,
    /// The nonminimal coupling `xi_2`.
    pub coupling: f64,
    /// The reduced Planck mass.
    pub planck_mass: Mass,
}

impl<V> MultifieldPotentialValue<V, 2> for NonminimalLoopInflationPotential
where
    V: From<f64>
        + Clone
        + Add<V, Output = V>
        + Mul<V, Output = V>
        + Mul<f64, Output = V>
        + Div<V, Output = V>
        + Transcendental,
{
    fn value(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        let planck = self.planck_mass.into_value();
        let [phi, s] = fields.map(Quantity::into_value);
        let inflaton = V::from(1.0) + (phi * (1.0 / planck)).ln() * self.alpha;
        let scaled = s * (1.0 / planck);
        let coupling = V::from(1.0) + scaled.clone() * scaled * self.coupling;
        EnergyDensity::new(inflaton * self.scale.into_value() / (coupling.clone() * coupling))
    }
}

/// Power-law inflation, `V_J = Lambda^4 exp(-lambda phi / M_pl)` (Sec. 4.2 of
/// the reference), in the Einstein frame of the assisted model. Its slow-roll
/// parameters are constant, so inflation must be ended by a separate sector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonminimalPowerLawPotential {
    /// `Lambda^4`, an overall energy density.
    pub scale: EnergyDensity,
    /// The slope `lambda`.
    pub lambda: f64,
    /// The nonminimal coupling `xi_2`.
    pub coupling: f64,
    /// The reduced Planck mass.
    pub planck_mass: Mass,
}

impl<V> MultifieldPotentialValue<V, 2> for NonminimalPowerLawPotential
where
    V: From<f64>
        + Clone
        + Add<V, Output = V>
        + Mul<V, Output = V>
        + Mul<f64, Output = V>
        + Div<V, Output = V>
        + Transcendental,
{
    fn value(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        let planck = self.planck_mass.into_value();
        let [phi, s] = fields.map(Quantity::into_value);
        let inflaton = (phi * (-self.lambda / planck)).exp();
        let scaled = s * (1.0 / planck);
        let coupling = V::from(1.0) + scaled.clone() * scaled * self.coupling;
        EnergyDensity::new(inflaton * self.scale.into_value() / (coupling.clone() * coupling))
    }
}

/// The inflationary valley of hybrid inflation,
/// `V_J = Lambda^4 [1 + (phi/mu)^2]` (Sec. 4.3 of the reference), in the
/// Einstein frame of the assisted model. The waterfall at `phi = phi_e` ends
/// inflation and is supplied as a field-value end-of-inflation event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonminimalHybridPotential {
    /// `Lambda^4`, an overall energy density.
    pub scale: EnergyDensity,
    /// The scale `mu`, in reduced Planck units.
    pub mu: f64,
    /// The nonminimal coupling `xi_2`.
    pub coupling: f64,
    /// The reduced Planck mass.
    pub planck_mass: Mass,
}

impl<V> MultifieldPotentialValue<V, 2> for NonminimalHybridPotential
where
    V: From<f64>
        + Clone
        + Add<V, Output = V>
        + Mul<V, Output = V>
        + Mul<f64, Output = V>
        + Div<V, Output = V>
        + Transcendental,
{
    fn value(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        let planck = self.planck_mass.into_value();
        let [phi, s] = fields.map(Quantity::into_value);
        let ratio = phi * (1.0 / (self.mu * planck));
        let inflaton = V::from(1.0) + ratio.clone() * ratio;
        let scaled = s * (1.0 / planck);
        let coupling = V::from(1.0) + scaled.clone() * scaled * self.coupling;
        EnergyDensity::new(inflaton * self.scale.into_value() / (coupling.clone() * coupling))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        HyperbolicFieldSpace, LowerFirst, MultifieldPotential, validate_riemann_symmetries,
    };

    fn geometry() -> NonminimalAssistedFieldSpace {
        NonminimalAssistedFieldSpace::new(0.3, 1.0, Mass::new(1.0)).unwrap()
    }

    fn position() -> [Field; 2] {
        [Field::new(7.0), Field::new(0.8)]
    }

    fn shifted(position: [Field; 2], index: usize, step: f64) -> [Field; 2] {
        let mut moved = position;
        moved[index] = moved[index] + Field::new(step);
        moved
    }

    /// Central difference of a component-valued function along coordinate `m`.
    fn derivative(f: impl Fn(&[Field; 2]) -> Vec<f64>, position: [Field; 2], m: usize) -> Vec<f64> {
        let step = 1.0e-5;
        let after = f(&shifted(position, m, step));
        let before = f(&shifted(position, m, -step));
        after
            .iter()
            .zip(before)
            .map(|(a, b)| (a - b) / (2.0 * step))
            .collect()
    }

    fn metric_components(g: &impl FieldSpaceGeometry<2>, p: &[Field; 2]) -> Vec<f64> {
        g.metric(p)
            .components()
            .iter()
            .flatten()
            .map(|c| c.into_value())
            .collect()
    }

    fn christoffel_components(g: &impl FieldSpaceGeometry<2>, p: &[Field; 2]) -> Vec<f64> {
        g.christoffel(p)
            .components()
            .iter()
            .flatten()
            .flatten()
            .map(|c| c.value().to_owned())
            .collect()
    }

    fn riemann_from_connection(
        g: &impl FieldSpaceGeometry<2>,
        p: [Field; 2],
    ) -> [[[[f64; 2]; 2]; 2]; 2] {
        // [nabla_K, nabla_L] A^I = R^I_JKL A^J gives
        // R^I_JKL = d_K G^I_LJ - d_L G^I_KJ + G^I_KM G^M_LJ - G^I_LM G^M_KJ.
        let gamma = christoffel_components(g, &p);
        let at = |i: usize, j: usize, k: usize| gamma[4 * i + 2 * j + k];
        let d = [0, 1].map(|m| derivative(|q| christoffel_components(g, q), p, m));
        let dat = |m: usize, i: usize, j: usize, k: usize| d[m][4 * i + 2 * j + k];
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|l| {
                        dat(k, i, l, j) - dat(l, i, k, j)
                            + (0..2)
                                .map(|m| at(i, k, m) * at(m, l, j) - at(i, l, m) * at(m, k, j))
                                .sum::<f64>()
                    })
                })
            })
        })
    }

    fn assert_close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * scale,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn metric_matches_the_einstein_frame_kinetic_functions() {
        let g = geometry();
        let s: f64 = 0.8;
        let xi = 0.3;
        let k1 = 1.0 / (1.0 + xi * s * s);
        let k2 = k1 + 6.0 * xi * xi * s * s * k1 * k1;
        let metric = g.metric(&position());
        assert_close(metric.components()[0][0].into_value(), k1, 1.0, 1.0e-15);
        assert_close(metric.components()[1][1].into_value(), k2, 1.0, 1.0e-15);
        assert_eq!(metric.components()[0][1].into_value(), 0.0);
        let inverse = g.inverse_metric(&position());
        assert_close(
            inverse.components()[1][1].into_value() * k2,
            1.0,
            1.0,
            1.0e-15,
        );

        let palatini = NonminimalAssistedFieldSpace::new(0.3, 0.0, Mass::new(1.0)).unwrap();
        let metric = palatini.metric(&position());
        assert_eq!(metric.components()[0][0], metric.components()[1][1]);
    }

    #[test]
    fn connection_matches_a_finite_difference_levi_civita_connection() {
        let g = geometry();
        let p = position();
        let inverse = g.inverse_metric(&p);
        let d = [0, 1].map(|m| derivative(|q| metric_components(&g, q), p, m));
        let dg = |m: usize, a: usize, b: usize| d[m][2 * a + b];
        let actual = christoffel_components(&g, &p);
        for (n, value) in actual.iter().enumerate() {
            let [i, j, k] = [n / 4, (n / 2) % 2, n % 2];
            let expected: f64 = (0..2)
                .map(|l| {
                    0.5 * inverse.components()[i][l].into_value()
                        * (dg(j, l, k) + dg(k, l, j) - dg(l, j, k))
                })
                .sum();
            assert_close(*value, expected, 1.0, 1.0e-9);
        }
        let derivatives = g.christoffel_derivative(&p);
        let numeric = [0, 1].map(|m| derivative(|q| christoffel_components(&g, q), p, m));
        for (m, numeric) in numeric.iter().enumerate() {
            for (n, value) in numeric.iter().enumerate() {
                let [i, j, k] = [n / 4, (n / 2) % 2, n % 2];
                assert_close(
                    derivatives.components()[i][j][k][m].value().to_owned(),
                    *value,
                    1.0,
                    1.0e-8,
                );
            }
        }
    }

    #[test]
    fn curvature_matches_the_connection_in_the_crate_convention() {
        // The finite-difference formula reproduces the built-in hyperbolic
        // geometry, which fixes the sign convention used below.
        let hyperbolic = HyperbolicFieldSpace::new(Field::new(2.0)).unwrap();
        let reference = hyperbolic.riemann(&position());
        let numeric = riemann_from_connection(&hyperbolic, position());
        let g = geometry();
        let actual = g.riemann(&position());
        let expected = riemann_from_connection(&g, position());
        for n in 0..16 {
            let [i, j, k, l] = [n / 8, (n / 4) % 2, (n / 2) % 2, n % 2];
            assert_close(
                numeric[i][j][k][l],
                reference.components()[i][j][k][l].into_value(),
                1.0,
                1.0e-7,
            );
            assert_close(
                actual.components()[i][j][k][l].into_value(),
                expected[i][j][k][l],
                1.0,
                1.0e-7,
            );
        }
        let lowered = g.metric(&position()).lower_first(&actual);
        validate_riemann_symmetries(&lowered, 1.0e-14).unwrap();
        assert!(actual.components()[1][0][1][0].into_value().abs() > 1.0e-2);
    }

    #[test]
    fn curvature_derivative_matches_its_covariant_finite_difference() {
        let g = geometry();
        let p = position();
        let riemann = |q: &[Field; 2]| -> Vec<f64> {
            g.riemann(q)
                .components()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .map(|c| c.into_value())
                .collect()
        };
        let r = riemann(&p);
        let at = |i: usize, j: usize, k: usize, l: usize| r[8 * i + 4 * j + 2 * k + l];
        let gamma = christoffel_components(&g, &p);
        let connection = |i: usize, j: usize, k: usize| gamma[4 * i + 2 * j + k];
        let d = [0, 1].map(|m| derivative(riemann, p, m));
        let actual = g.covariant_riemann_derivative(&p);
        for n in 0..32 {
            let [i, j, k, l, m] = [n / 16, (n / 8) % 2, (n / 4) % 2, (n / 2) % 2, n % 2];
            let expected = d[m][8 * i + 4 * j + 2 * k + l]
                + (0..2)
                    .map(|q| {
                        connection(i, m, q) * at(q, j, k, l)
                            - connection(q, m, j) * at(i, q, k, l)
                            - connection(q, m, k) * at(i, j, q, l)
                            - connection(q, m, l) * at(i, j, k, q)
                    })
                    .sum::<f64>();
            assert_close(
                actual.components()[i][j][k][l][m].into_value(),
                expected,
                1.0,
                1.0e-7,
            );
        }
    }

    #[test]
    fn loop_potential_matches_its_closed_form() {
        let potential = NonminimalLoopInflationPotential {
            scale: EnergyDensity::new(2.0e-10),
            alpha: 0.1,
            coupling: 0.3,
            planck_mass: Mass::new(1.0),
        };
        let [phi, s]: [f64; 2] = [7.0, 0.8];
        let expected = 2.0e-10 * (1.0 + 0.1 * phi.ln()) / (1.0 + 0.3 * s * s).powi(2);
        let value = <NonminimalLoopInflationPotential as MultifieldPotentialValue<f64, 2>>::value(
            &potential,
            position(),
        );
        assert_close(value.into_value(), expected, expected, 1.0e-15);
        let gradient = potential.gradient(position());
        let step = 1.0e-6;
        for (index, component) in gradient.iter().enumerate() {
            let at = |q: [Field; 2]| {
                <NonminimalLoopInflationPotential as MultifieldPotentialValue<f64, 2>>::value(
                    &potential, q,
                )
                .into_value()
            };
            let numeric = (at(shifted(position(), index, step))
                - at(shifted(position(), index, -step)))
                / (2.0 * step);
            assert_close(component.into_value(), numeric, numeric.abs(), 1.0e-7);
        }
    }

    #[test]
    fn class_two_potentials_match_their_closed_forms() {
        let [phi, s]: [f64; 2] = [7.0, 0.8];
        let coupling = (1.0_f64 + 0.3 * s * s).powi(2);
        let power_law = NonminimalPowerLawPotential {
            scale: EnergyDensity::new(2.0e-10),
            lambda: 0.05,
            coupling: 0.3,
            planck_mass: Mass::new(1.0),
        };
        let hybrid = NonminimalHybridPotential {
            scale: EnergyDensity::new(2.0e-10),
            mu: 10.0,
            coupling: 0.3,
            planck_mass: Mass::new(1.0),
        };
        let expected_power = 2.0e-10 * (-0.05 * phi).exp() / coupling;
        let expected_hybrid = 2.0e-10 * (1.0 + (phi / 10.0).powi(2)) / coupling;
        let power_value = <NonminimalPowerLawPotential as MultifieldPotentialValue<f64, 2>>::value(
            &power_law,
            position(),
        );
        let hybrid_value = <NonminimalHybridPotential as MultifieldPotentialValue<f64, 2>>::value(
            &hybrid,
            position(),
        );
        assert_close(
            power_value.into_value(),
            expected_power,
            expected_power,
            1.0e-15,
        );
        assert_close(
            hybrid_value.into_value(),
            expected_hybrid,
            expected_hybrid,
            1.0e-15,
        );
        // Automatic-differentiation gradients against the closed forms.
        let power_gradient = power_law.gradient(position());
        assert_close(
            power_gradient[0].into_value(),
            -0.05 * expected_power,
            expected_power,
            1.0e-13,
        );
        let hybrid_gradient = hybrid.gradient(position());
        assert_close(
            hybrid_gradient[0].into_value(),
            2.0e-10 * 2.0 * phi / 100.0 / coupling,
            expected_hybrid,
            1.0e-13,
        );
        assert_close(
            hybrid_gradient[1].into_value(),
            -4.0 * 0.3 * s / (1.0 + 0.3 * s * s) * expected_hybrid,
            expected_hybrid,
            1.0e-13,
        );
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        assert!(matches!(
            NonminimalAssistedFieldSpace::new(0.1, 0.5, Mass::new(1.0)),
            Err(NonminimalParameterError::InvalidFormulation(_))
        ));
        assert!(matches!(
            NonminimalAssistedFieldSpace::new(-0.1, 1.0, Mass::new(1.0)),
            Err(NonminimalParameterError::InvalidCoupling(_))
        ));
        assert!(matches!(
            NonminimalAssistedFieldSpace::new(0.1, 1.0, Mass::new(0.0)),
            Err(NonminimalParameterError::InvalidPlanckMass(_))
        ));
    }
}
