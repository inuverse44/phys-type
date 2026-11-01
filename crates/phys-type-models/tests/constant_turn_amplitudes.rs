//! The transport bispectrum on the constant-turn orbit against the
//! independent in-in reference. See
//! `docs/constant-turn-bispectrum-amplitudes.md`.

use phys_type::{BispectrumRunSettings, TriangleShape};
use phys_type_models::constant_turn::{self, HUBBLE, radius};
use phys_type_models::in_in::{
    MixingTheory, ModeSettings, ModeTable, VertexSettings, bispectrum, orbit_terms,
};

/// Transport f_NL for isosceles triangles k1 = k2 = k_pivot, k3 = r k_pivot
/// on the slice `after_exit` e-folds after the pivot exit, with the lambda
/// and mu^2 measured at that exit.
fn transport(
    epsilon: f64,
    lambda: f64,
    mu_squared: f64,
    after_exit: f64,
    ratios: &[f64],
    depth: f64,
) -> (Vec<f64>, f64, f64) {
    let (system, initial, settings) = constant_turn::orbit(
        lambda,
        mu_squared,
        radius(lambda, epsilon),
        after_exit,
        16.0,
    );
    let survey = system.survey(initial, &settings).unwrap();
    let exit = survey
        .trajectory
        .state_at(survey.e_folds_at_end - settings.pivot_e_folds_before_end);
    let h = system.hubble(&exit).into_value();
    let eta_perp = system
        .oriented_kinematics(&exit, None)
        .unwrap()
        .eta_perp
        .into_value();
    let measured_lambda = 2.0 * eta_perp.abs();
    let measured_mu_squared = constant_turn::radial_mass_squared(lambda, mu_squared)
        * (HUBBLE / h).powi(2)
        + 3.0 * eta_perp * eta_perp;
    let shapes: Vec<TriangleShape> = ratios
        .iter()
        .map(|&r| TriangleShape {
            k_t_over_pivot: 2.0 + r,
            alpha: 0.0,
            beta: (2.0 - r) / (2.0 + r),
        })
        .collect();
    let records = system
        .run_bispectrum(
            initial,
            &BispectrumRunSettings {
                background: settings,
                e_folds_before_exit: depth,
                ..BispectrumRunSettings::default()
            },
            &shapes,
        )
        .unwrap()
        .records;
    (
        records.iter().map(|record| record.f_nl).collect(),
        measured_lambda,
        measured_mu_squared,
    )
}

fn in_in(lambda: f64, mu_squared: f64, after_exit: f64, r: f64) -> f64 {
    let table = ModeTable::new(
        MixingTheory { lambda, mu_squared },
        ModeSettings {
            z_end: (-(after_exit + 6.0)).exp(),
            ..ModeSettings::default()
        },
    );
    bispectrum(
        &table,
        &orbit_terms(lambda),
        [1.0, 1.0, r],
        (-after_exit).exp(),
        &VertexSettings::default(),
    )
    .reduced_f_nl()
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn constant_turn_transport_matches_the_in_in_amplitudes() {
    // The transport evolves the full theory on the orbit; the in-in
    // reference evaluates its decoupling-limit cubic action with the exact
    // linear modes and shares no code with it. With epsilon = 1e-6 the
    // decoupling corrections are below the transport's own start-depth
    // error, about 1e-5 in f_NL.
    use phys_type_models::in_in::{CubicOperator, CubicTerm};
    let ratios = [1.0, 0.01];
    for (lambda, mu_squared, after_exit) in [(1.0, 1.25, 20.0), (2.0, 2.0, 12.0)] {
        let (f_transport, l, m) = transport(1.0e-6, lambda, mu_squared, after_exit, &ratios, 5.0);
        for (r, f) in ratios.iter().zip(&f_transport) {
            let reference = in_in(l, m, after_exit, *r);
            assert!(
                (f - reference).abs() < 5.0e-5,
                "({lambda}, {mu_squared}) r {r}: transport {f} vs in-in {reference}"
            );
        }
        // Without the gradient coupling, which lies outside the basis of
        // arXiv:2607.14529, the in-in result misses the transport by far.
        let table = ModeTable::new(
            MixingTheory {
                lambda: l,
                mu_squared: m,
            },
            ModeSettings {
                z_end: (-(after_exit + 6.0)).exp(),
                ..ModeSettings::default()
            },
        );
        let without_gradient: Vec<CubicTerm> = orbit_terms(l)
            .into_iter()
            .filter(|term| term.operator != CubicOperator::SigmaGradientPhiSquared)
            .collect();
        let partial = bispectrum(
            &table,
            &without_gradient,
            [1.0; 3],
            (-after_exit).exp(),
            &VertexSettings::default(),
        )
        .reduced_f_nl();
        assert!(
            (partial - f_transport[0]).abs() > 0.1,
            "{partial} vs {}",
            f_transport[0]
        );
    }
}
