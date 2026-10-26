//! The three trajectories of arXiv:1210.7487 (KMS) off the ridge `chi = 0`:
//! their background (Figs. 3 and 7) and the crossing-time delta N of
//! Sec. V (Fig. 11). The tests are expensive in debug builds; run them with
//! `cargo test --release -- --ignored`. See
//! `docs/general-nonminimal-inflation.md`.

use phys_type::{BackgroundSurvey, Formulation, Mass, SpectrumRunSettings};
use phys_type_models::kms::{
    KmsSystem, PHI_START, TRAJECTORIES, at_rest, benchmark_system, crossing_delta_n, turn_peak,
};

fn settings() -> SpectrumRunSettings {
    SpectrumRunSettings {
        pivot_e_folds_before_end: 60.0,
        e_folds_before_exit: 6.0,
        background_sample_every: 1,
        ..SpectrumRunSettings::default()
    }
}

fn survey(system: &KmsSystem, chi0: f64) -> BackgroundSurvey<2> {
    system
        .survey(at_rest(PHI_START, chi0), &settings())
        .unwrap()
}

#[test]
#[ignore = "expensive; run in release mode"]
fn the_trajectories_fall_off_the_ridge_as_in_kms_figs_3_and_7() {
    let system = benchmark_system(Formulation::Metric);
    let surveys = TRAJECTORIES.map(|chi0| survey(&system, chi0));
    // Fig. 3: N_total = 71.5 for trajectory 2.
    let total = surveys[1].e_folds_at_end;
    assert!((total - 71.5).abs() < 0.05, "trajectory 2: {total} e-folds");
    // Fig. 7: the turn rate peaks 34.5 and 4.9 e-folds before the end; the
    // exact evolution gives 34.06 and 5.15, a spacing of 28.9 against 29.6.
    let (peak_1, before_1) = turn_peak(&surveys[0]).unwrap();
    let (peak_2, before_2) = turn_peak(&surveys[1]).unwrap();
    assert!(peak_1 > 0.5 && peak_2 > 0.5, "{peak_1}, {peak_2}");
    assert!(
        (before_1 - 34.5).abs() < 1.0,
        "trajectory 1 peak at {before_1}"
    );
    assert!(
        (before_2 - 4.9).abs() < 0.5,
        "trajectory 2 peak at {before_2}"
    );
    // Trajectory 3 stays on the ridge until inflation ends.
    assert_eq!(turn_peak(&surveys[2]), None);
    // Trajectory 2 nearly stops inflating at the turn: epsilon_H reaches
    // 0.99 there, which is why its f_NL is large and steep in chi(tau_0).
    let turning = surveys[1]
        .samples
        .iter()
        .filter(|sample| surveys[1].e_folds_at_end - sample.e_folds > 1.0)
        .map(|sample| sample.epsilon_h)
        .fold(0.0, f64::max);
    assert!(
        turning > 0.98 && turning < 1.0,
        "epsilon_H at the turn {turning}"
    );
}

#[test]
#[ignore = "expensive; run in release mode"]
fn the_crossing_delta_n_reproduces_kms_fig_11() {
    let system = benchmark_system(Formulation::Metric);
    for (index, chi0) in TRAJECTORIES.into_iter().enumerate() {
        let survey = survey(&system, chi0);
        let delta_n = crossing_delta_n(&system, &survey, &settings(), 60.0).unwrap();
        if index == 1 {
            // Sec. V: |f_NL| = 43.3 for trajectory 2.
            assert!(
                (delta_n.f_nl.abs() - 43.3).abs() < 0.3,
                "trajectory 2: f_NL {}",
                delta_n.f_nl
            );
        } else {
            assert!(
                delta_n.f_nl.abs() < 0.1,
                "trajectory {}: f_NL {}",
                index + 1,
                delta_n.f_nl
            );
        }
        // The linear part against the exact mode solver at the end: the
        // slow-roll velocities at crossing cost about 5%.
        let pivot = survey.trajectory.state_at(survey.e_folds_at_end - 60.0);
        let k = (pivot.a * system.hubble(&pivot)).into_value();
        let exact = system
            .mode_spectrum(
                &settings(),
                &survey.trajectory,
                Mass::new(k),
                survey.e_folds_at_end - 0.01,
            )
            .unwrap()
            .spectra
            .curvature
            .into_value();
        assert!(
            (delta_n.power / exact - 1.0).abs() < 0.06,
            "trajectory {}: delta N {} vs mode solver {exact}",
            index + 1,
            delta_n.power
        );
    }
}
