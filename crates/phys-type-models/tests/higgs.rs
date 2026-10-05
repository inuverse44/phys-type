//! Validation of single-field Higgs inflation (arXiv:0710.3755, metric;
//! arXiv:0803.2664, Palatini) evolved exactly with one field. The tests are
//! expensive in debug builds; run them with `cargo test --release -- --ignored`.

use phys_type::{
    CovariantMultifieldSystem, Field, FieldVelocity, Formulation, Mass, MultifieldState,
    ScaleFactor, SpectrumRunSettings, einstein_frame,
};
use phys_type_models::higgs::HiggsTheory;

/// Exact `(n_s, r)` with the pivot 60 Einstein-frame e-folds before the end.
fn exact(theory: HiggsTheory) -> (f64, f64) {
    let (potential, geometry) = einstein_frame(theory);
    let system = CovariantMultifieldSystem::in_planck_units(potential, geometry);
    let resting = MultifieldState::new(
        [Field::new(theory.field_at_e_folds(72.0))],
        [FieldVelocity::new(0.0)],
        ScaleFactor::new(1.0),
    );
    let initial = MultifieldState {
        pi: system.slow_roll_velocity(&resting),
        ..resting
    };
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: 60.0,
        output_e_folds_before_end: 1.0,
        e_folds_before_exit: 6.0,
        ..SpectrumRunSettings::default()
    };
    let survey = system.survey(initial, &settings).unwrap();
    let output = survey.e_folds_at_end - 1.0;
    let pivot = survey.trajectory.state_at(survey.e_folds_at_end - 60.0);
    let k = (pivot.a * system.hubble(&pivot)).into_value();
    let curvature = |kk: f64| {
        system
            .mode_spectrum(&settings, &survey.trajectory, Mass::new(kk), output)
            .unwrap()
            .spectra
            .curvature
            .into_value()
    };
    let step = 0.02_f64;
    let n_s = 1.0 + (curvature(k * step.exp()) / curvature(k * (-step).exp())).ln() / (2.0 * step);
    let tensor = system
        .tensor_spectrum(&settings, &survey.trajectory, Mass::new(k), output)
        .unwrap();
    (n_s, tensor / curvature(k))
}

fn theory(xi: f64, formulation: Formulation) -> HiggsTheory {
    HiggsTheory {
        xi,
        lambda: 1.0e-12 * (1.0 + xi * xi),
        formulation,
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn higgs_inflation_agrees_with_slow_roll_in_both_formulations() {
    for formulation in [Formulation::Metric, Formulation::Palatini] {
        for xi in [1.0e-2, 1.0, 1.0e3] {
            let model = theory(xi, formulation);
            let (n_s, r) = exact(model);
            let (n_s_slow_roll, r_slow_roll) = model.slow_roll_predictions(60.0);
            // Differences are second order in slow roll.
            assert!(
                (n_s - n_s_slow_roll).abs() < 5.0e-4,
                "{formulation:?} {xi}: n_s {n_s} vs {n_s_slow_roll}"
            );
            assert!(
                (r / r_slow_roll - 1.0).abs() < 0.03,
                "{formulation:?} {xi}: r {r} vs {r_slow_roll}"
            );
        }
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn higgs_inflation_reaches_its_large_coupling_limits() {
    let e_folds = 60.0_f64;
    let (n_s, r) = exact(theory(1.0e4, Formulation::Metric));
    assert!((n_s - (1.0 - 2.0 / e_folds)).abs() < 2.0e-3, "{n_s}");
    // r sits about 10% below 12/N^2 at N = 60, an O(1/N) correction.
    assert!((r / (12.0 / e_folds.powi(2)) - 0.9).abs() < 0.05, "{r}");
    // Palatini: r = 2/(xi N^2) up to O(1/N), so xi r is independent of xi.
    let (_, r_low) = exact(theory(1.0e3, Formulation::Palatini));
    let (_, r_high) = exact(theory(1.0e4, Formulation::Palatini));
    assert!((r_low * 1.0e3 / (r_high * 1.0e4) - 1.0).abs() < 0.01);
    assert!((r_high / (2.0 / (1.0e4 * e_folds.powi(2))) - 1.0).abs() < 0.05);
}
