//! Validation of the numerical in-in reference, `phys_type_models::in_in`.
//! See `docs/constant-turn-bispectrum-amplitudes.md`. The ignored tests take
//! about a second in release mode and minutes in debug mode; run them with
//! `cargo test --release -- --ignored`.

use phys_type_models::hmps::exact_massive_ratio;
use phys_type_models::in_in::{
    CubicOperator, CubicTerm, InInBispectrum, MixingTheory, ModeSettings, ModeTable,
    VertexSettings, bispectrum, covariant_paper_terms, orbit_terms, paper_terms,
    phi_dot_cubed_bispectrum,
};

const SQUEEZED: [f64; 3] = [1.0, 1.0, 0.01];
const EQUILATERAL: [f64; 3] = [1.0, 1.0, 1.0];

fn table(lambda: f64, mu_squared: f64, z_end: f64) -> ModeTable {
    ModeTable::new(
        MixingTheory { lambda, mu_squared },
        ModeSettings {
            z_end,
            ..ModeSettings::default()
        },
    )
}

fn run(table: &ModeTable, terms: &[CubicTerm], k: [f64; 3]) -> InInBispectrum {
    bispectrum(
        table,
        terms,
        k,
        (-12.0_f64).exp(),
        &VertexSettings::default(),
    )
}

#[test]
fn power_spectrum_matches_the_exact_solution() {
    // Eq. (147) of arXiv:2607.14529. A light sigma decays as z^{3/2 - nu}
    // outside the horizon, so the final slice is taken where that is below
    // 1e-6.
    for (lambda, mu_squared, z_end) in [
        (1.0, 1.25, -30.0),
        (2.0, 1.25, -30.0),
        (0.5, 2.0, -20.0),
        (2.0, 2.0, -20.0),
        (1.0, 4.0, -12.0),
        (2.0, 4.0, -12.0),
    ] {
        let z_end = f64::exp(z_end);
        let ratio = table(lambda, mu_squared, z_end).power_ratio(z_end);
        let exact = exact_massive_ratio(lambda, mu_squared);
        assert!(
            (ratio / exact - 1.0).abs() < 2.0e-6,
            "{lambda} {mu_squared}: {ratio} vs {exact}"
        );
    }
}

#[test]
fn phi_dot_cubed_matches_its_closed_form() {
    // B = 3 g / (k1 k2 k3 K^3) without mixing: tests the vertex integral,
    // the window that replaces i0+, and the sign of the in-in formula.
    let t = table(0.0, 2.0, (-14.0_f64).exp());
    let term = [CubicTerm {
        operator: CubicOperator::PhiDotCubed,
        coefficient: 1.0,
    }];
    for k in [EQUILATERAL, [1.0, 1.0, 1.9], [1.0, 1.0, 0.1], SQUEEZED] {
        let b = bispectrum(&t, &term, k, (-8.0_f64).exp(), &VertexSettings::default());
        let exact = phi_dot_cubed_bispectrum(1.0, k);
        assert!(
            (b.bispectrum / exact - 1.0).abs() < 2.0e-6,
            "{k:?}: {} vs {exact}",
            b.bispectrum
        );
    }
}

#[test]
fn entropy_vertices_need_mixing() {
    let t = table(0.0, 1.25, (-12.0_f64).exp());
    for operator in [
        CubicOperator::SigmaPhiDotSquared,
        CubicOperator::SigmaGradientPhiSquared,
        CubicOperator::SigmaSquaredPhiDot,
        CubicOperator::SigmaCubed,
    ] {
        let b = run(
            &t,
            &[CubicTerm {
                operator,
                coefficient: 1.0,
            }],
            EQUILATERAL,
        );
        assert_eq!(b.bispectrum, 0.0, "{operator:?}");
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn the_covariant_basis_reshuffles_the_couplings() {
    // Eq. (8): alpha = alpha~, beta = beta~ + 2 lambda alpha~,
    // gamma = gamma~ + 3 lambda beta~ + 3 lambda^2 alpha~.
    let (lambda, a, b, g) = (1.3, 0.7, -0.4, 0.9);
    let t = table(lambda, 1.6, (-18.0_f64).exp());
    for k in [EQUILATERAL, SQUEEZED] {
        let tilde = run(&t, &paper_terms(a, b, g), k).bispectrum;
        let covariant = run(
            &t,
            &covariant_paper_terms(
                a,
                b + 2.0 * lambda * a,
                g + 3.0 * lambda * b + 3.0 * lambda * lambda * a,
            ),
            k,
        )
        .bispectrum;
        assert!(
            ((tilde - covariant) / tilde).abs() < 1.0e-9,
            "{tilde} vs {covariant}"
        );
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn the_orbit_is_the_paper_basis_plus_a_gradient_coupling() {
    let lambda = 1.5;
    let t = table(lambda, 1.25, (-12.0_f64).exp());
    let orbit = run(&t, &orbit_terms(lambda), EQUILATERAL).bispectrum;
    let paper = run(
        &t,
        &paper_terms(lambda, -0.5 * lambda * lambda, 0.0),
        EQUILATERAL,
    );
    let gradient = run(
        &t,
        &[CubicTerm {
            operator: CubicOperator::SigmaGradientPhiSquared,
            coefficient: -0.5 * lambda,
        }],
        EQUILATERAL,
    );
    let sum = paper.bispectrum + gradient.bispectrum;
    assert!(((orbit - sum) / orbit).abs() < 1.0e-9, "{orbit} vs {sum}");
    // The gradient coupling is not a small part of the orbit's bispectrum.
    assert!(gradient.bispectrum.abs() > 0.1 * orbit.abs());
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn weak_mixing_is_linear_in_lambda() {
    // At one transfer the alpha bispectrum is proportional to lambda, and
    // the power spectrum changes only at order lambda^2.
    let f = |lambda: f64| {
        let t = table(lambda, 2.0, (-14.0_f64).exp());
        run(&t, &paper_terms(1.0, 0.0, 0.0), EQUILATERAL).reduced_f_nl() / lambda
    };
    // f / lambda is even in lambda, so its first correction is O(lambda^2):
    // doubling lambda quadruples the change.
    let [f1, f2, f4] = [1.0e-3, 2.0e-3, 4.0e-3].map(f);
    assert!(f1.abs() > 1.0e-2);
    assert!(((f2 - f1) / f1).abs() < 1.0e-4, "{f2} vs {f1}");
    let ratio = (f4 - f2) / (f2 - f1);
    assert!((ratio - 4.0).abs() < 0.05, "{ratio}");
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn the_orbit_bispectrum_converges_in_every_numerical_setting() {
    let defaults = VertexSettings::default();
    for (lambda, mu_squared) in [(1.0, 1.25), (2.0, 2.0), (2.0, 0.81)] {
        let base = ModeSettings {
            z_start: 6.0_f64.exp(),
            z_end: (-18.0_f64).exp(),
            step: 1.0e-2,
        };
        let f = |modes: ModeSettings, vertex: VertexSettings, k: [f64; 3]| {
            let t = ModeTable::new(MixingTheory { lambda, mu_squared }, modes);
            bispectrum(&t, &orbit_terms(lambda), k, (-12.0_f64).exp(), &vertex).reduced_f_nl()
        };
        for k in [EQUILATERAL, SQUEEZED] {
            let reference = f(base, defaults, k);
            let variants = [
                f(
                    ModeSettings {
                        z_start: 7.5_f64.exp(),
                        ..base
                    },
                    defaults,
                    k,
                ),
                f(
                    ModeSettings {
                        step: 5.0e-3,
                        ..base
                    },
                    defaults,
                    k,
                ),
                f(
                    base,
                    VertexSettings {
                        oscillation_width: 0.5,
                        logarithmic_width: 0.025,
                        ..defaults
                    },
                    k,
                ),
                f(
                    base,
                    VertexSettings {
                        window_onset: 0.25,
                        ..defaults
                    },
                    k,
                ),
                f(
                    base,
                    VertexSettings {
                        response_below: 0.03,
                        ..defaults
                    },
                    k,
                ),
            ];
            for (index, variant) in variants.into_iter().enumerate() {
                assert!(
                    ((variant - reference) / reference).abs() < 1.0e-6,
                    "{lambda} {mu_squared} {k:?} variant {index}: {variant} vs {reference}"
                );
            }
        }
    }
}
