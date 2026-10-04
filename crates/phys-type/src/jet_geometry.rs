//! Field-space geometry computed from a metric and its derivatives.
//!
//! A [`MetricJet`] supplies the covariant metric `G_IJ` at a position as
//! [`ThirdOrderJet`]s, that is, with its first three partial derivatives.
//! [`JetGeometry`] then derives every quantity the engine needs, by the same
//! formulas for any metric:
//!
//! ```text
//! Gamma^I_JK  = G^IL (partial_J G_LK + partial_K G_LJ - partial_L G_JK) / 2
//! R^I_JKL     = partial_K Gamma^I_LJ - partial_L Gamma^I_KJ
//!               + Gamma^I_KM Gamma^M_LJ - Gamma^I_LM Gamma^M_KJ
//! nabla_P R^I_JKL = partial_P R^I_JKL + Gamma^I_PQ R^Q_JKL - Gamma^Q_PJ R^I_QKL
//!               - Gamma^Q_PK R^I_JQL - Gamma^Q_PL R^I_JKQ
//! ```
//!
//! The Riemann tensor follows the crate convention
//! `[nabla_K, nabla_L] A^I = R^I_JKL A^J`. The connection is carried as a
//! second-order [`Jet`], so its first and second partial derivatives, and
//! with them `partial_P R`, come from automatic differentiation rather than
//! from hand-derived formulas. A geometry written by hand can be checked
//! against this one, and a geometry derived from an action, such as the
//! Einstein frame of a nonminimal theory, needs nothing else.
//!
//! Positions and derivatives use the raw field values, which are in reduced
//! Planck units throughout the crate.

use crate::{
    ChristoffelDerivatives, ChristoffelSymbols, CovariantRiemannDerivative,
    CurvatureDerivativeGeometry, CurvedFieldSpaceGeometry, Field, FieldMetric, FieldSpaceGeometry,
    GeometryDomainError, InverseFieldMetric, Jet, Quantity, RiemannTensor,
    ThirdOrderFieldSpaceGeometry, ThirdOrderJet,
};

/// A covariant field-space metric with its first three partial derivatives.
pub trait MetricJet<const N: usize> {
    /// `G_IJ` and its derivatives with respect to the `N` field coordinates,
    /// evaluated at `position`. The matrix must be symmetric.
    fn metric_jet(&self, position: &[Field; N]) -> [[ThirdOrderJet<N>; N]; N];

    /// Checks conditions under which the metric is defined at `position`,
    /// beyond the positive definiteness that [`JetGeometry`] always checks.
    fn check_domain(&self, _position: &[Field; N]) -> Result<(), GeometryDomainError> {
        Ok(())
    }
}

/// The geometry of a [`MetricJet`]: inverse metric, connection, curvature,
/// and their derivatives, all derived from the metric.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JetGeometry<M>(pub M);

/// The second-order jet of `partial_slot` applied to a third-order jet.
fn derivative_jet<const N: usize>(jet: &ThirdOrderJet<N>, slot: usize) -> Jet<N> {
    Jet {
        value: jet.gradient[slot],
        gradient: jet.hessian[slot],
        hessian: jet.third[slot],
    }
}

/// The second-order truncation of a third-order jet.
fn truncate<const N: usize>(jet: &ThirdOrderJet<N>) -> Jet<N> {
    Jet {
        value: jet.value,
        gradient: jet.gradient,
        hessian: jet.hessian,
    }
}

/// Gauss--Jordan inverse of a symmetric positive-definite jet matrix.
fn invert<const N: usize>(matrix: [[Jet<N>; N]; N]) -> [[Jet<N>; N]; N] {
    let mut left = matrix;
    let mut right: [[Jet<N>; N]; N] =
        core::array::from_fn(|i| core::array::from_fn(|j| Jet::constant(f64::from(i == j))));
    for pivot in 0..N {
        let scale = Jet::constant(1.0) / left[pivot][pivot];
        for column in 0..N {
            left[pivot][column] = left[pivot][column] * scale;
            right[pivot][column] = right[pivot][column] * scale;
        }
        for row in 0..N {
            if row != pivot {
                let factor = left[row][pivot];
                for column in 0..N {
                    left[row][column] = left[row][column] - factor * left[pivot][column];
                    right[row][column] = right[row][column] - factor * right[pivot][column];
                }
            }
        }
    }
    right
}

/// Connection, Riemann tensor, and its partial derivative at one position.
struct Curvature<const N: usize> {
    /// `Gamma^I_JK` as a second-order jet.
    connection: [[[Jet<N>; N]; N]; N],
    /// `R^I_JKL`.
    riemann: [[[[f64; N]; N]; N]; N],
    /// `partial_P R^I_JKL`, derivative index last.
    riemann_derivative: [[[[[f64; N]; N]; N]; N]; N],
}

impl<M> JetGeometry<M> {
    /// The connection as a second-order jet.
    fn connection<const N: usize>(&self, position: &[Field; N]) -> [[[Jet<N>; N]; N]; N]
    where
        M: MetricJet<N>,
    {
        let metric = self.0.metric_jet(position);
        let inverse = invert(core::array::from_fn(|i| {
            core::array::from_fn(|j| truncate(&metric[i][j]))
        }));
        // partial_L G_JK as derivative[L][J][K].
        let derivative: [[[Jet<N>; N]; N]; N] = core::array::from_fn(|l| {
            core::array::from_fn(|j| core::array::from_fn(|k| derivative_jet(&metric[j][k], l)))
        });
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    (0..N).fold(Jet::constant(0.0), |sum, l| {
                        sum + inverse[i][l]
                            * (derivative[j][l][k] + derivative[k][l][j] - derivative[l][j][k])
                            * 0.5
                    })
                })
            })
        })
    }

    fn curvature<const N: usize>(&self, position: &[Field; N]) -> Curvature<N>
    where
        M: MetricJet<N>,
    {
        let connection = self.connection(position);
        let g = |i: usize, j: usize, k: usize| connection[i][j][k].value;
        let dg = |i: usize, j: usize, k: usize, p: usize| connection[i][j][k].gradient[p];
        let ddg =
            |i: usize, j: usize, k: usize, p: usize, q: usize| connection[i][j][k].hessian[p][q];
        let riemann: [[[[f64; N]; N]; N]; N] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|l| {
                        dg(i, l, j, k) - dg(i, k, j, l)
                            + (0..N)
                                .map(|m| g(i, k, m) * g(m, l, j) - g(i, l, m) * g(m, k, j))
                                .sum::<f64>()
                    })
                })
            })
        });
        let riemann_derivative = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|l| {
                        core::array::from_fn(|p| {
                            ddg(i, l, j, k, p) - ddg(i, k, j, l, p)
                                + (0..N)
                                    .map(|m| {
                                        dg(i, k, m, p) * g(m, l, j) + g(i, k, m) * dg(m, l, j, p)
                                            - dg(i, l, m, p) * g(m, k, j)
                                            - g(i, l, m) * dg(m, k, j, p)
                                    })
                                    .sum::<f64>()
                        })
                    })
                })
            })
        });
        Curvature {
            connection,
            riemann,
            riemann_derivative,
        }
    }
}

impl<M: MetricJet<N>, const N: usize> FieldSpaceGeometry<N> for JetGeometry<M> {
    fn metric(&self, position: &[Field; N]) -> FieldMetric<N> {
        let metric = self.0.metric_jet(position);
        FieldMetric::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| Quantity::new(metric[i][j].value))
        }))
    }

    fn inverse_metric(&self, position: &[Field; N]) -> InverseFieldMetric<N> {
        let metric = self.0.metric_jet(position);
        let inverse = invert(core::array::from_fn(|i| {
            core::array::from_fn(|j| truncate(&metric[i][j]))
        }));
        InverseFieldMetric::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| Quantity::new(inverse[i][j].value))
        }))
    }

    fn check_domain(&self, position: &[Field; N]) -> Result<(), GeometryDomainError> {
        self.0.check_domain(position)?;
        // Cholesky factorization: every pivot positive iff G is positive definite.
        let metric = self.0.metric_jet(position);
        let mut factor = [[0.0; N]; N];
        for i in 0..N {
            for j in 0..=i {
                let sum: f64 = (0..j).map(|k| factor[i][k] * factor[j][k]).sum();
                if i == j {
                    let pivot = metric[i][i].value - sum;
                    if !(pivot.is_finite() && pivot > 0.0) {
                        return Err(GeometryDomainError {
                            reason: "the field-space metric is not positive definite",
                            value: pivot,
                        });
                    }
                    factor[i][i] = pivot.sqrt();
                } else {
                    factor[i][j] = (metric[i][j].value - sum) / factor[j][j];
                }
            }
        }
        Ok(())
    }

    fn christoffel(&self, position: &[Field; N]) -> ChristoffelSymbols<N> {
        let connection = self.connection(position);
        ChristoffelSymbols::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| Quantity::new(connection[i][j][k].value))
            })
        }))
    }
}

impl<M: MetricJet<N>, const N: usize> CurvedFieldSpaceGeometry<N> for JetGeometry<M> {
    fn riemann(&self, position: &[Field; N]) -> RiemannTensor<N> {
        let riemann = self.curvature(position).riemann;
        RiemannTensor::new(riemann.map(|a| a.map(|b| b.map(|c| c.map(Quantity::new)))))
    }
}

impl<M: MetricJet<N>, const N: usize> ThirdOrderFieldSpaceGeometry<N> for JetGeometry<M> {
    fn christoffel_derivative(&self, position: &[Field; N]) -> ChristoffelDerivatives<N> {
        let connection = self.connection(position);
        ChristoffelDerivatives::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|m| Quantity::new(connection[i][j][k].gradient[m]))
                })
            })
        }))
    }
}

impl<M: MetricJet<N>, const N: usize> CurvatureDerivativeGeometry<N> for JetGeometry<M> {
    fn covariant_riemann_derivative(&self, position: &[Field; N]) -> CovariantRiemannDerivative<N> {
        let Curvature {
            connection,
            riemann,
            riemann_derivative,
        } = self.curvature(position);
        let g = |i: usize, j: usize, k: usize| connection[i][j][k].value;
        CovariantRiemannDerivative::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    core::array::from_fn(|l| {
                        core::array::from_fn(|p| {
                            let correction: f64 = (0..N)
                                .map(|q| {
                                    g(i, p, q) * riemann[q][j][k][l]
                                        - g(q, p, j) * riemann[i][q][k][l]
                                        - g(q, p, k) * riemann[i][j][q][l]
                                        - g(q, p, l) * riemann[i][j][k][q]
                                })
                                .sum();
                            Quantity::new(riemann_derivative[i][j][k][l][p] + correction)
                        })
                    })
                })
            })
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        HyperbolicFieldSpace, LowerFirst, NonminimalAssistedFieldSpace, PolarFieldSpace,
        Transcendental, validate_riemann_symmetries,
    };

    fn variables<const N: usize>(position: &[Field; N]) -> [ThirdOrderJet<N>; N] {
        core::array::from_fn(|i| ThirdOrderJet::variable(position[i].into_value(), i))
    }

    fn diagonal<const N: usize>(entries: [ThirdOrderJet<N>; N]) -> [[ThirdOrderJet<N>; N]; N] {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    entries[i]
                } else {
                    ThirdOrderJet::constant(0.0)
                }
            })
        })
    }

    struct PolarMetric;

    impl MetricJet<2> for PolarMetric {
        fn metric_jet(&self, position: &[Field; 2]) -> [[ThirdOrderJet<2>; 2]; 2] {
            let [rho, _] = variables(position);
            diagonal([ThirdOrderJet::constant(1.0), rho * rho])
        }
    }

    struct HyperbolicMetric(f64);

    impl MetricJet<2> for HyperbolicMetric {
        fn metric_jet(&self, position: &[Field; 2]) -> [[ThirdOrderJet<2>; 2]; 2] {
            let [r, _] = variables(position);
            let x = r / self.0;
            // sinh^2 x = (cosh 2x - 1)/2 = (e^{2x} + e^{-2x} - 2)/4.
            let sinh2 = ((x * 2.0).exp() + (x * -2.0).exp() - ThirdOrderJet::constant(2.0)) / 4.0;
            diagonal([ThirdOrderJet::constant(1.0), sinh2])
        }
    }

    /// The assisted metric diag(K_1(s), K_2(s)) of `NonminimalAssistedFieldSpace`.
    struct AssistedMetric {
        xi: f64,
        kappa: f64,
    }

    impl MetricJet<2> for AssistedMetric {
        fn metric_jet(&self, position: &[Field; 2]) -> [[ThirdOrderJet<2>; 2]; 2] {
            let [_, s] = variables(position);
            let one = ThirdOrderJet::constant(1.0);
            let k1 = one / (one + s * s * self.xi);
            let k2 = k1 + s * s * k1 * k1 * (6.0 * self.kappa * self.xi * self.xi);
            diagonal([k1, k2])
        }
    }

    /// A non-diagonal three-dimensional metric with no symmetry.
    struct GenericMetric;

    impl MetricJet<3> for GenericMetric {
        fn metric_jet(&self, position: &[Field; 3]) -> [[ThirdOrderJet<3>; 3]; 3] {
            let [x, y, z] = variables(position);
            let one = ThirdOrderJet::constant(1.0);
            let a = one + x * x * 0.3 + (y * 0.5).tanh() * 0.2;
            let b = (z * 0.4).exp() + y * y * 0.1;
            let c = one * 2.0 + (x * y).cos() * 0.3;
            let ab = x * z * 0.1;
            let ac = (x + z * 0.5).tanh() * 0.15;
            let bc = y * 0.05 + z * z * 0.02;
            [[a, ab, ac], [ab, b, bc], [ac, bc, c]]
        }
    }

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs()))
    }

    fn assert_same_geometry<A, B>(jet: &A, reference: &B, position: [Field; 2], tolerance: f64)
    where
        A: CurvatureDerivativeGeometry<2> + ThirdOrderFieldSpaceGeometry<2>,
        B: CurvatureDerivativeGeometry<2> + ThirdOrderFieldSpaceGeometry<2>,
    {
        let p = &position;
        for i in 0..2 {
            for j in 0..2 {
                let (m, n) = (jet.metric(p), reference.metric(p));
                assert!(close(
                    m.components()[i][j].into_value(),
                    n.components()[i][j].into_value(),
                    tolerance
                ));
                let (m, n) = (jet.inverse_metric(p), reference.inverse_metric(p));
                assert!(close(
                    m.components()[i][j].into_value(),
                    n.components()[i][j].into_value(),
                    tolerance
                ));
                for k in 0..2 {
                    let (a, b) = (jet.christoffel(p), reference.christoffel(p));
                    assert!(close(
                        a.components()[i][j][k].into_value(),
                        b.components()[i][j][k].into_value(),
                        tolerance
                    ));
                    for l in 0..2 {
                        let (a, b) = (
                            jet.christoffel_derivative(p),
                            reference.christoffel_derivative(p),
                        );
                        assert!(
                            close(
                                a.components()[i][j][k][l].into_value(),
                                b.components()[i][j][k][l].into_value(),
                                tolerance
                            ),
                            "dGamma {i}{j}{k}{l}"
                        );
                        let (a, b) = (jet.riemann(p), reference.riemann(p));
                        assert!(
                            close(
                                a.components()[i][j][k][l].into_value(),
                                b.components()[i][j][k][l].into_value(),
                                tolerance
                            ),
                            "R {i}{j}{k}{l}"
                        );
                        for m in 0..2 {
                            let (a, b) = (
                                jet.covariant_riemann_derivative(p),
                                reference.covariant_riemann_derivative(p),
                            );
                            assert!(
                                close(
                                    a.components()[i][j][k][l][m].into_value(),
                                    b.components()[i][j][k][l][m].into_value(),
                                    tolerance
                                ),
                                "nabla R {i}{j}{k}{l}{m}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn jet_geometry_reproduces_the_polar_field_space() {
        let position = [Field::new(0.3), Field::new(2.0)];
        assert_same_geometry(
            &JetGeometry(PolarMetric),
            &PolarFieldSpace,
            position,
            1.0e-13,
        );
    }

    #[test]
    fn jet_geometry_reproduces_the_hyperbolic_field_space() {
        let position = [Field::new(3.0), Field::new(1.0)];
        let reference = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        assert_same_geometry(
            &JetGeometry(HyperbolicMetric(7.0)),
            &reference,
            position,
            1.0e-12,
        );
    }

    #[test]
    fn jet_geometry_reproduces_the_assisted_field_space() {
        for (xi, kappa) in [(0.022, 1.0), (0.3, 1.0), (0.3, 0.0)] {
            for s in [0.05, 0.7, 2.4] {
                let position = [Field::new(1.3), Field::new(s)];
                let reference =
                    NonminimalAssistedFieldSpace::new(xi, kappa, crate::Mass::new(1.0)).unwrap();
                assert_same_geometry(
                    &JetGeometry(AssistedMetric { xi, kappa }),
                    &reference,
                    position,
                    1.0e-12,
                );
            }
        }
    }

    #[test]
    fn jet_curvature_satisfies_the_algebraic_and_second_bianchi_identities() {
        let geometry = JetGeometry(GenericMetric);
        let position = [Field::new(0.4), Field::new(-0.7), Field::new(0.9)];
        let lowered = geometry
            .metric(&position)
            .lower_first(&geometry.riemann(&position));
        validate_riemann_symmetries(&lowered, 1.0e-12).unwrap();
        // nabla_M R^I_JKL + nabla_K R^I_JLM + nabla_L R^I_JMK = 0.
        let derivative = geometry.covariant_riemann_derivative(&position);
        let d = |i: usize, j: usize, k: usize, l: usize, m: usize| {
            derivative.components()[i][j][k][l][m].into_value()
        };
        let mut largest = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        for m in 0..3 {
                            largest = largest.max(d(i, j, k, l, m).abs());
                            let cyclic = d(i, j, k, l, m) + d(i, j, l, m, k) + d(i, j, m, k, l);
                            assert!(cyclic.abs() < 1.0e-12, "{i}{j}{k}{l}{m}: {cyclic}");
                        }
                    }
                }
            }
        }
        // The identity is not satisfied trivially.
        assert!(largest > 1.0e-3);
    }
}
