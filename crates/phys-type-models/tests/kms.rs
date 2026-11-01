//! The Einstein frame of the two-field model of arXiv:1210.7487 (KMS)
//! against its closed forms: both fields couple to `R`, so the metric is not
//! diagonal and the curvature varies in both directions.

use phys_type::{
    CurvatureDerivativeGeometry, CurvedFieldSpaceGeometry, Field, FieldSpaceGeometry, Formulation,
    einstein_frame,
};
use phys_type_models::kms::KmsTheory;

const FORMULATIONS: [Formulation; 2] = [Formulation::Metric, Formulation::Palatini];

/// Points on the ridge, off it, and near the end of inflation, in units of
/// `M_pl`; KMS start at `phi = 3.1`.
const POINTS: [(f64, f64); 4] = [(3.1, 1.1e-3), (1.5, 0.4), (0.2, -0.3), (0.05, 0.02)];

fn close(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs()))
}

/// `delta^I_K G_JL - delta^I_L G_JK`, the two-dimensional Riemann tensor at
/// unit Gaussian curvature.
fn unit_riemann(g: &[[f64; 2]; 2], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let delta = |a: usize, b: usize| f64::from(u8::from(a == b));
    delta(i, k) * g[j][l] - delta(i, l) * g[j][k]
}

/// A model with every coupling distinct, so that no symmetry hides an error.
fn generic(formulation: Formulation) -> KmsTheory {
    KmsTheory {
        xi_phi: 3.0,
        xi_chi: 0.7,
        lambda_phi: 1.0e-2,
        lambda_chi: 3.0e-3,
        g: 2.0e-3,
        m_phi_squared: 1.0e-4,
        m_chi_squared: 2.0e-4,
        formulation,
    }
}

fn theories() -> Vec<KmsTheory> {
    FORMULATIONS
        .iter()
        .flat_map(|&f| [KmsTheory::benchmark(1.0e-2, f), generic(f)])
        .collect()
}

#[test]
fn the_einstein_metric_is_kms_eq_110() {
    for theory in theories() {
        let (_, geometry) = einstein_frame(theory);
        for (phi, chi) in POINTS {
            let metric = geometry.metric(&[Field::new(phi), Field::new(chi)]);
            let expected = theory.metric(phi, chi);
            for (i, row) in expected.iter().enumerate() {
                for (j, &value) in row.iter().enumerate() {
                    let actual = metric.components()[i][j].into_value();
                    assert!(
                        close(actual, value, 1.0e-13),
                        "{theory:?} at ({phi}, {chi}), G_{i}{j}: {actual} vs {value}"
                    );
                }
            }
        }
    }
}

#[test]
fn the_riemann_tensor_is_fixed_by_the_closed_form_gaussian_curvature() {
    for theory in theories() {
        let (_, geometry) = einstein_frame(theory);
        for (phi, chi) in POINTS {
            let position = [Field::new(phi), Field::new(chi)];
            let riemann = geometry.riemann(&position);
            let g = theory.metric(phi, chi);
            let k = theory.gaussian_curvature(phi, chi);
            for i in 0..2 {
                for j in 0..2 {
                    for m in 0..2 {
                        for l in 0..2 {
                            let expected = k * unit_riemann(&g, i, j, m, l);
                            let actual = riemann.components()[i][j][m][l].into_value();
                            assert!(
                                close(actual, expected, 1.0e-10),
                                "{theory:?} at ({phi}, {chi}), R^{i}_{j}{m}{l}: {actual} vs {expected}"
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn the_covariant_riemann_derivative_follows_the_gradient_of_the_curvature() {
    // In two dimensions nabla_M R^I_JKL = (partial_M K)(delta^I_K G_JL -
    // delta^I_L G_JK), so a central difference of the closed-form K checks
    // the rank-five tensor that the cubic transport uses.
    let step = 1.0e-5;
    for theory in theories() {
        let (_, geometry) = einstein_frame(theory);
        for (phi, chi) in POINTS {
            let derivative =
                geometry.covariant_riemann_derivative(&[Field::new(phi), Field::new(chi)]);
            let gradient = [
                (theory.gaussian_curvature(phi + step, chi)
                    - theory.gaussian_curvature(phi - step, chi))
                    / (2.0 * step),
                (theory.gaussian_curvature(phi, chi + step)
                    - theory.gaussian_curvature(phi, chi - step))
                    / (2.0 * step),
            ];
            let g = theory.metric(phi, chi);
            for i in 0..2 {
                for j in 0..2 {
                    for m in 0..2 {
                        for l in 0..2 {
                            for (n, dk) in gradient.iter().enumerate() {
                                let expected = dk * unit_riemann(&g, i, j, m, l);
                                let actual = derivative.components()[i][j][m][l][n].into_value();
                                assert!(
                                    close(actual, expected, 1.0e-6),
                                    "{theory:?} at ({phi}, {chi}), nabla_{n} R^{i}_{j}{m}{l}: \
                                     {actual} vs {expected}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn the_metric_and_palatini_curvatures_agree_at_zero_coupling() {
    // At xi = 0 both formulations are flat.
    for formulation in FORMULATIONS {
        let theory = KmsTheory {
            xi_phi: 0.0,
            xi_chi: 0.0,
            ..generic(formulation)
        };
        assert_eq!(theory.gaussian_curvature(1.2, -0.4), 0.0);
    }
}
