//! The three trajectories of arXiv:1210.7487 (KMS) off the ridge `chi = 0`:
//! their background (Figs. 3 and 7) and the crossing-time delta N of
//! Sec. V (Fig. 11). The tests are expensive in debug builds; run them with
//! `cargo test --release -- --ignored`. See
//! `docs/general-nonminimal-inflation.md`.

use phys_type::{
    BackgroundSurvey, CovariantMultifieldSystem, EinsteinGeometry, EinsteinPotential, Field,
    FieldSpaceGeometry, FieldVelocity, Formulation, Mass, MultifieldState, ScaleFactor,
    SpectrumRunSettings, einstein_frame,
};
use phys_type_models::kms::KmsTheory;

type System =
    CovariantMultifieldSystem<EinsteinPotential<KmsTheory>, EinsteinGeometry<KmsTheory>, 2>;

/// `chi(tau_0)` of trajectories 1--3, KMS Fig. 2; `phi(tau_0) = 3.1`.
const TRAJECTORIES: [f64; 3] = [1.1e-2, 1.1e-3, 1.1e-4];

fn system() -> System {
    let (potential, geometry) = einstein_frame(KmsTheory::benchmark(1.0e-2, Formulation::Metric));
    CovariantMultifieldSystem::in_planck_units(potential, geometry)
}

fn settings() -> SpectrumRunSettings {
    SpectrumRunSettings {
        pivot_e_folds_before_end: 60.0,
        e_folds_before_exit: 6.0,
        background_sample_every: 1,
        ..SpectrumRunSettings::default()
    }
}

/// KMS do not state the initial velocities; starting from rest reproduces
/// their total of 71.5 e-folds, the slow-roll velocity gives 71.2.
fn at_rest(phi: f64, chi: f64) -> MultifieldState<2> {
    MultifieldState::new(
        [Field::new(phi), Field::new(chi)],
        [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
        ScaleFactor::new(1.0),
    )
}

fn survey(system: &System, chi0: f64) -> BackgroundSurvey<2> {
    system.survey(at_rest(3.1, chi0), &settings()).unwrap()
}

/// The largest `|eta_perp|` during inflation and the e-folds before the end
/// at which it occurs. The turn rate is undefined at the resting start.
fn turn_peak(survey: &BackgroundSurvey<2>) -> (f64, f64) {
    survey
        .samples
        .iter()
        .filter(|sample| sample.eta_perp.is_finite())
        .map(|sample| (sample.eta_perp, survey.e_folds_at_end - sample.e_folds))
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
}

#[test]
#[ignore = "expensive; run in release mode"]
fn the_trajectories_fall_off_the_ridge_as_in_kms_figs_3_and_7() {
    let system = system();
    let surveys = TRAJECTORIES.map(|chi0| survey(&system, chi0));
    // Fig. 3: N_total = 71.5 for trajectory 2.
    let total = surveys[1].e_folds_at_end;
    assert!((total - 71.5).abs() < 0.05, "trajectory 2: {total} e-folds");
    // Fig. 7: the turn rate peaks 34.5 and 4.9 e-folds before the end. The
    // exact evolution gives 32.9 and 4.6; the spacing of the two, set by the
    // tachyonic growth of chi on the ridge, is 28.3 here, 27.7 in slow roll,
    // and 29.6 in KMS.
    let (peak_1, before_1) = turn_peak(&surveys[0]);
    let (peak_2, before_2) = turn_peak(&surveys[1]);
    assert!(peak_1 > 10.0 && peak_2 > 10.0, "{peak_1}, {peak_2}");
    assert!(
        (before_1 - 34.5).abs() < 2.0,
        "trajectory 1 peak at {before_1}"
    );
    assert!(
        (before_2 - 4.9).abs() < 0.5,
        "trajectory 2 peak at {before_2}"
    );
    // Trajectory 3 stays on the ridge until inflation ends.
    let (peak_3, _) = turn_peak(&surveys[2]);
    assert!(peak_3 < 1.0, "trajectory 3 turn rate {peak_3}");
}

/// The crossing-time delta N of KMS Sec. V: N to `epsilon_H = 1` from the
/// field values at the crossing of the pivot, 60 e-folds before the end,
/// with slow-roll velocities. Returns `f_NL` in the convention of this
/// repository (KMS Eq. (92) carries the opposite sign) and
/// `P_zeta = (H_*/2 pi)^2 N_a N^a`.
fn crossing_delta_n(system: &System, survey: &BackgroundSurvey<2>) -> (f64, f64) {
    let star = survey.trajectory.state_at(survey.e_folds_at_end - 60.0);
    let x0 = [star.phi[0].into_value(), star.phi[1].into_value()];
    let quiet = SpectrumRunSettings {
        background_sample_every: usize::MAX,
        ..settings()
    };
    let e_folds = |a: f64, b: f64, h: [f64; 2]| {
        let rest = at_rest(x0[0] + a * h[0], x0[1] + b * h[1]);
        let state = MultifieldState {
            pi: system.slow_roll_velocity(&rest),
            ..rest
        };
        system.survey(state, &quiet).unwrap().e_folds_at_end
    };
    let h = [3.0e-4 * x0[0].abs(), 3.0e-4 * x0[1].abs()];
    let centre = e_folds(0.0, 0.0, h);
    let axis = [
        [e_folds(1.0, 0.0, h), e_folds(-1.0, 0.0, h)],
        [e_folds(0.0, 1.0, h), e_folds(0.0, -1.0, h)],
    ];
    let mixed = (e_folds(1.0, 1.0, h) - e_folds(1.0, -1.0, h) - e_folds(-1.0, 1.0, h)
        + e_folds(-1.0, -1.0, h))
        / (4.0 * h[0] * h[1]);
    let first: [f64; 2] = core::array::from_fn(|a| (axis[a][0] - axis[a][1]) / (2.0 * h[a]));
    let partial: [[f64; 2]; 2] = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            if a == b {
                (axis[a][0] - 2.0 * centre + axis[a][1]) / (h[a] * h[a])
            } else {
                mixed
            }
        })
    });
    let position = [Field::new(x0[0]), Field::new(x0[1])];
    let connection = system.geometry.christoffel(&position);
    let inverse = system.geometry.inverse_metric(&position);
    let covariant: [[f64; 2]; 2] = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            partial[a][b]
                - (0..2)
                    .map(|c| connection.components()[c][a][b].value() * first[c])
                    .sum::<f64>()
        })
    });
    let raised: [f64; 2] = core::array::from_fn(|a| {
        (0..2)
            .map(|b| inverse.components()[a][b].into_value() * first[b])
            .sum()
    });
    let norm: f64 = (0..2).map(|a| raised[a] * first[a]).sum();
    let curvature: f64 = (0..4)
        .map(|m| raised[m / 2] * raised[m % 2] * covariant[m / 2][m % 2])
        .sum();
    let hubble = system.hubble(&star).into_value();
    (
        5.0 / 6.0 * curvature / (norm * norm),
        (hubble / (2.0 * core::f64::consts::PI)).powi(2) * norm,
    )
}

#[test]
#[ignore = "expensive; run in release mode"]
fn the_crossing_delta_n_reproduces_kms_fig_11() {
    let system = system();
    for (index, chi0) in TRAJECTORIES.into_iter().enumerate() {
        let survey = survey(&system, chi0);
        let (f_nl, power) = crossing_delta_n(&system, &survey);
        if index == 1 {
            // Sec. V: |f_NL| = 43.3 for trajectory 2.
            assert!((f_nl.abs() - 43.3).abs() < 0.3, "trajectory 2: f_NL {f_nl}");
        } else {
            assert!(f_nl.abs() < 0.1, "trajectory {}: f_NL {f_nl}", index + 1);
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
            (power / exact - 1.0).abs() < 0.06,
            "trajectory {}: delta N {power} vs mode solver {exact}",
            index + 1
        );
    }
}
