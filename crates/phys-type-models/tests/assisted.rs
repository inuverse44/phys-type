//! Validation of nonminimally assisted inflation against arXiv:2302.05866.
//!
//! The model is `AssistedTheory` on the generalized nonminimal theory layer;
//! the suites compare the exact background, mode, and three-point evolution
//! with the reference's slow-roll prescription (`hkkpt`), with numerical
//! delta N on the horizon-crossing slice, and with phase-space delta N on
//! later slices. See `docs/nonminimal-assisted-inflation.md`. All tests here
//! are expensive; run them with `cargo test --release -- --ignored`.

use std::f64::consts::PI;

use phys_type::{
    BackgroundTrajectory, CovariantMultifieldSystem, CurvatureDerivativeGeometry, EinsteinGeometry,
    EinsteinPotential, ErrorTolerances, Field, FieldSpaceGeometry, FieldVelocity, Mass,
    MomentumTriangle, MultifieldPotential, MultifieldState, ScaleFactor, SpectrumRunSettings,
    ThirdOrderFieldSpaceGeometry, ThriceDifferentiableMultifieldPotential, einstein_frame,
};
use phys_type_models::assisted::{
    AssistedInflaton, AssistedTheory, HybridInflaton, LoopInflaton, PowerLawInflaton, bisect,
};
use phys_type_models::hkkpt::{
    assisted_delta_n_prediction, assisted_slow_roll, paper_prescription,
};

type AssistedSystem<I> = CovariantMultifieldSystem<
    EinsteinPotential<AssistedTheory<I>>,
    EinsteinGeometry<AssistedTheory<I>>,
    2,
>;

/// A surveyed system and a triangle of momenta, for three-point transport.
struct SourcedBenchmark<P, G> {
    system: CovariantMultifieldSystem<P, G, 2>,
    trajectory: BackgroundTrajectory,
    triangle: MomentumTriangle,
}

/// Surveys `system` from `initial`. The triangle is `shape` times the
/// pivot wavenumber, or times `pivot` when one is supplied.
fn sourced_benchmark_with<P, G>(
    system: CovariantMultifieldSystem<P, G, 2>,
    initial: MultifieldState<2>,
    shape: [f64; 3],
    pivot: Option<f64>,
    settings: &SpectrumRunSettings,
) -> SourcedBenchmark<P, G>
where
    P: ThriceDifferentiableMultifieldPotential<2> + Copy,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2> + Copy,
{
    let survey = system.survey(initial, settings).unwrap();
    let pivot_state = survey
        .trajectory
        .state_at(survey.e_folds_at_end - settings.pivot_e_folds_before_end);
    let k = pivot.unwrap_or((pivot_state.a * system.hubble(&pivot_state)).into_value());
    let triangle = MomentumTriangle::try_new(shape[0] * k, shape[1] * k, shape[2] * k).unwrap();
    SourcedBenchmark {
        system,
        trajectory: survey.trajectory,
        triangle,
    }
}

impl<P, G> SourcedBenchmark<P, G>
where
    P: ThriceDifferentiableMultifieldPotential<2> + Copy,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2> + Copy,
{
    fn evolve(
        &self,
        e_folds_before_exit: f64,
        target: f64,
        relative_tolerance: f64,
    ) -> phys_type::CanonicalThreePointTransportState<2> {
        self.system
            .three_point_state(
                &self.trajectory,
                self.triangle,
                e_folds_before_exit,
                target,
                relative_tolerance,
                2_000_000,
            )
            .unwrap()
    }
}

/// The curvature correlators of the benchmark on the slice `slice`.
fn project_benchmark<P, G>(
    bench: &SourcedBenchmark<P, G>,
    depth: f64,
    slice: f64,
) -> phys_type::CurvatureCorrelators
where
    P: ThriceDifferentiableMultifieldPotential<2> + Copy,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2> + Copy,
{
    let state = bench.evolve(depth, slice, 1.0e-9);
    phys_type::curvature_correlators(
        &bench.system,
        &bench.trajectory.state_at(slice),
        bench.triangle,
        &state,
        0.1,
    )
    .unwrap()
}

fn assisted_start<P: MultifieldPotential<2>, G: FieldSpaceGeometry<2>>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    phi: f64,
    s: f64,
) -> MultifieldState<2> {
    let resting = MultifieldState::new(
        [Field::new(phi), Field::new(s)],
        [FieldVelocity::new(0.0); 2],
        ScaleFactor::new(1.0),
    );
    MultifieldState {
        pi: system.slow_roll_velocity(&resting),
        ..resting
    }
}

struct CalibratedAssisted<I: AssistedInflaton> {
    inflaton: I,
    s_star: f64,
    coupling: f64,
    system: AssistedSystem<I>,
    initial: MultifieldState<2>,
    trajectory: BackgroundTrajectory,
    end: f64,
    settings: SpectrumRunSettings,
}

/// The assisted model with `kappa = 1`, started about 72 e-folds of phi
/// slow roll before the end, with the initial s tuned so that s = s_star
/// at the pivot 60 e-folds before the end. The output slice lies
/// `output_before_end` e-folds before the end.
fn calibrated_assisted<I: AssistedInflaton>(
    inflaton: I,
    s_star: f64,
    coupling: f64,
    output_before_end: f64,
) -> CalibratedAssisted<I> {
    let (potential, geometry) = einstein_frame(AssistedTheory::metric(inflaton, coupling));
    let system = CovariantMultifieldSystem::in_planck_units(potential, geometry);
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: 60.0,
        output_e_folds_before_end: output_before_end,
        k_min_over_pivot: 1.0,
        k_max_over_pivot: 1.0,
        k_count: 1,
        end_of_inflation: inflaton.end_of_inflation(),
        ..SpectrumRunSettings::default()
    };
    let (low, high) = inflaton.phi_bracket();
    let phi0 = bisect(&|phi| inflaton.slow_roll_e_folds(phi) - 72.0, low, high);
    let mut s0 = 0.35 * s_star;
    let mut survey = system
        .survey(assisted_start(&system, phi0, s0), &settings)
        .unwrap();
    for _ in 0..6 {
        let pivot = survey.trajectory.state_at(survey.e_folds_at_end - 60.0);
        s0 *= s_star / pivot.phi[1].into_value();
        survey = system
            .survey(assisted_start(&system, phi0, s0), &settings)
            .unwrap();
    }
    let (trajectory, end) = (survey.trajectory, survey.e_folds_at_end);
    CalibratedAssisted {
        inflaton,
        s_star,
        coupling,
        system,
        initial: assisted_start(&system, phi0, s0),
        trajectory,
        end,
        settings,
    }
}

fn calibrated_assisted_loop(
    alpha: f64,
    s_star: f64,
    coupling: f64,
) -> CalibratedAssisted<LoopInflaton> {
    calibrated_assisted(LoopInflaton { alpha }, s_star, coupling, 1.0e-3)
}

fn calibrated_assisted_loop_with_output(
    alpha: f64,
    s_star: f64,
    coupling: f64,
    output_before_end: f64,
) -> CalibratedAssisted<LoopInflaton> {
    calibrated_assisted(LoopInflaton { alpha }, s_star, coupling, output_before_end)
}

impl<I: AssistedInflaton> CalibratedAssisted<I> {
    fn pivot(&self) -> MultifieldState<2> {
        self.trajectory.state_at(self.end - 60.0)
    }

    fn k_pivot(&self) -> f64 {
        let pivot = self.pivot();
        (pivot.a * self.system.hubble(&pivot)).into_value()
    }

    fn output(&self) -> f64 {
        self.end - self.settings.output_e_folds_before_end
    }

    fn mode_solver_tilt(&self, slice: f64) -> f64 {
        let k = self.k_pivot();
        let spectrum = |kk: f64| {
            self.system
                .mode_spectrum(&self.settings, &self.trajectory, Mass::new(kk), slice)
                .unwrap()
                .spectra
                .curvature
                .into_value()
        };
        let step = 0.02_f64;
        (spectrum(k * step.exp()).ln() - spectrum(k * (-step).exp()).ln()) / (2.0 * step)
    }

    /// e-folds from `state` to the uniform-density surface of the output.
    fn e_folds_to_output(&self, state: &MultifieldState<2>, shift: [f64; 2]) -> f64 {
        let precise = SpectrumRunSettings {
            tolerances: ErrorTolerances::new(
                phys_type::ConformalPerturbationTolerance::<2>::background_in_planck_units(1.0e-24),
                1.0e-13,
            ),
            ..self.settings
        };
        let target = self.system.hubble(&self.trajectory.state_at(self.output()));
        let start = assisted_start(
            &self.system,
            state.phi[0].into_value() + shift[0],
            state.phi[1].into_value() + shift[1],
        );
        let perturbed = self.system.survey(start, &precise).unwrap().trajectory;
        perturbed
            .crossing(|st| (self.system.hubble(st) - target).into_value())
            .unwrap()
    }

    /// delta N power and f_NL from the horizon-crossing slice at `n_star`,
    /// with massless crossing statistics, covariant in the field metric.
    fn delta_n_at_crossing(&self, n_star: f64) -> (f64, f64, f64) {
        let state = self.trajectory.state_at(n_star);
        let h = [
            1.0e-4 * state.phi[0].into_value(),
            (1.0e-2 * state.phi[1].into_value()).min(1.0e-5),
        ];
        let value = |i: i32, j: i32| {
            self.e_folds_to_output(&state, [f64::from(i) * h[0], f64::from(j) * h[1]])
        };
        let centre = value(0, 0);
        let first = [
            (value(1, 0) - value(-1, 0)) / (2.0 * h[0]),
            (value(0, 1) - value(0, -1)) / (2.0 * h[1]),
        ];
        let mixed =
            (value(1, 1) - value(1, -1) - value(-1, 1) + value(-1, -1)) / (4.0 * h[0] * h[1]);
        let second = [
            [
                (value(1, 0) - 2.0 * centre + value(-1, 0)) / (h[0] * h[0]),
                mixed,
            ],
            [
                mixed,
                (value(0, 1) - 2.0 * centre + value(0, -1)) / (h[1] * h[1]),
            ],
        ];
        // Covariant second derivative N_;IJ = N_,IJ - Gamma^K_IJ N_,K.
        let connection =
            phys_type::FieldSpaceGeometry::christoffel(&self.system.geometry, &state.phi);
        let covariant: [[f64; 2]; 2] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                second[i][j]
                    - (0..2)
                        .map(|k| connection.components()[k][i][j].value() * first[k])
                        .sum::<f64>()
            })
        });
        let metric = phys_type::FieldSpaceGeometry::metric(&self.system.geometry, &state.phi);
        let inverse = [
            1.0 / metric.components()[0][0].into_value(),
            1.0 / metric.components()[1][1].into_value(),
        ];
        let norm: f64 = (0..2).map(|i| inverse[i] * first[i] * first[i]).sum();
        let contraction: f64 = (0..4)
            .map(|m| {
                let [i, j] = [m / 2, m % 2];
                inverse[i] * inverse[j] * first[i] * first[j] * covariant[i][j]
            })
            .sum();
        let hubble = self.system.hubble(&state).into_value();
        let k_star = (state.a * self.system.hubble(&state)).into_value();
        (
            k_star,
            (hubble / (2.0 * PI)).powi(2) * norm,
            5.0 / 6.0 * contraction / (norm * norm),
        )
    }

    /// Eqs. (3.1) and (3.5) with the numerical pivot and s_e. For Class I
    /// the end value of phi follows the reference, epsilon^(0) = 1; for
    /// Class II it is the numerical end, the field-value event itself.
    fn analytic_prediction(&self) -> (f64, f64) {
        let pivot = self.pivot();
        let final_state = self.trajectory.state_at(self.end);
        let phi_e = if self.inflaton.class_one() {
            self.inflaton.paper_phi_e()
        } else {
            final_state.phi[0].into_value()
        };
        assisted_delta_n_prediction(
            &assisted_slow_roll(
                &self.inflaton,
                self.coupling,
                pivot.phi[0].into_value(),
                pivot.phi[1].into_value(),
            ),
            &assisted_slow_roll(
                &self.inflaton,
                self.coupling,
                phi_e,
                final_state.phi[1].into_value(),
            ),
            self.inflaton.class_one(),
        )
    }

    /// Eqs. (3.1) and (3.5) evaluated entirely with the prescription of
    /// the reference: phi_* from 60 e-folds of phi slow roll, and s_e from
    /// Eq. (3.9) with N_e = 60.
    fn paper_prescription(&self) -> (f64, f64, f64) {
        paper_prescription(&self.inflaton, self.coupling, self.s_star)
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_loop_spectral_index_agrees_with_delta_n_and_eq_3_1() {
    let model = calibrated_assisted_loop(0.1, 1.0e-3, 0.022);
    let pivot = model.end - 60.0;
    let mode_solver = 1.0 + model.mode_solver_tilt(model.output());
    let (k_minus, p_minus, _) = model.delta_n_at_crossing(pivot - 0.5);
    let (k_plus, p_plus, _) = model.delta_n_at_crossing(pivot + 0.5);
    let delta_n = 1.0 + (p_plus / p_minus).ln() / (k_plus / k_minus).ln();
    let (analytic, _) = model.analytic_prediction();

    assert!(
        (mode_solver - delta_n).abs() < 2.0e-3,
        "{mode_solver} vs delta N {delta_n}"
    );
    assert!(
        (mode_solver - analytic).abs() < 2.0e-3,
        "{mode_solver} vs Eq. (3.1) {analytic}"
    );
    // zeta is not conserved before the end: half an e-fold earlier the
    // tilt is several percent lower.
    let earlier = 1.0 + model.mode_solver_tilt(model.end - 0.5);
    assert!(
        mode_solver - earlier > 2.0e-2,
        "{earlier} then {mode_solver}"
    );
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_loop_crossing_delta_n_underestimates_the_converged_f_nl() {
    let model = calibrated_assisted_loop(0.1, 1.0e-3, 0.022);
    let k = model.k_pivot();
    let bench = sourced_benchmark_with(
        model.system,
        model.initial,
        [1.0, 1.0, 1.0],
        Some(k),
        &model.settings,
    );
    let shallow = project_benchmark(&bench, 4.0, model.output()).reduced_f_nl();
    let deep = project_benchmark(&bench, 5.0, model.output()).reduced_f_nl();
    let (_, _, delta_n) = model.delta_n_at_crossing(model.end - 60.0);
    let (_, analytic) = model.analytic_prediction();

    assert!(
        (deep - shallow).abs() < 1.0e-4 * deep.abs(),
        "{shallow} then {deep}"
    );
    // The crossing-time delta N assumes massless, uncorrelated field
    // perturbations and so misses the tachyonic growth of the assistant
    // direction; the reference's Eq. (3.5) shares that assumption. Both
    // agree with each other but fall about 20% short of the converged
    // transport, which the phase-space delta N test confirms.
    assert!(deep < 0.0 && delta_n < 0.0 && analytic < 0.0);
    assert!(
        (analytic - delta_n).abs() < 0.15 * delta_n.abs(),
        "{analytic} vs {delta_n}"
    );
    assert!(
        delta_n.abs() < 0.85 * deep.abs() && analytic.abs() < 0.9 * deep.abs(),
        "{deep} vs crossing {delta_n} and Eq. (3.5) {analytic}"
    );
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_loop_eq_3_1_tracks_the_numerical_tilt_across_xi() {
    // Given the numerical s_e, Eq. (3.1) follows the exact tilt while n_s
    // moves by 0.04 across this range of the coupling.
    for xi in [0.021, 0.023, 0.024] {
        let model = calibrated_assisted_loop(0.1, 1.0e-3, xi);
        let numerical = 1.0 + model.mode_solver_tilt(model.output());
        let (analytic, _) = model.analytic_prediction();
        assert!(
            (numerical - analytic).abs() < 1.0e-3,
            "xi {xi}: {numerical} vs Eq. (3.1) {analytic}"
        );
    }
}

/// Real part of a velocity-basis two-point component with phase-space
/// indices a, b in 0..4 = (Q^phi, Q^s, V^phi, V^s).
fn covariance_entry(
    covariance: &phys_type::CanonicalTwoPointPhaseSpaceCovariance<2>,
    a: usize,
    b: usize,
) -> f64 {
    let (i, j) = (a % 2, b % 2);
    match (a < 2, b < 2) {
        (true, true) => covariance.qq.components()[i][j].value().re,
        (true, false) => covariance.qp.components()[i][j].value().re,
        (false, true) => covariance.pq.components()[i][j].value().re,
        (false, false) => covariance.pp.components()[i][j].value().re,
    }
}

fn three_point_entry(
    alpha: &phys_type::ThreePointPhaseSpaceCorrelation<2>,
    a: usize,
    b: usize,
    c: usize,
) -> f64 {
    let (i, j, k) = (a % 2, b % 2, c % 2);
    match (a < 2, b < 2, c < 2) {
        (true, true, true) => alpha.qqq.components()[i][j][k].value().re,
        (true, true, false) => alpha.qqp.components()[i][j][k].value().re,
        (true, false, true) => alpha.qpq.components()[i][j][k].value().re,
        (false, true, true) => alpha.pqq.components()[i][j][k].value().re,
        (true, false, false) => alpha.qpp.components()[i][j][k].value().re,
        (false, true, false) => alpha.pqp.components()[i][j][k].value().re,
        (false, false, true) => alpha.ppq.components()[i][j][k].value().re,
        (false, false, false) => alpha.ppp.components()[i][j][k].value().re,
    }
}

impl<I: AssistedInflaton> CalibratedAssisted<I> {
    /// e-folds from a phase-space point to the output uniform-density surface.
    fn e_folds_from_phase_space(&self, point: [f64; 4]) -> f64 {
        let precise = SpectrumRunSettings {
            tolerances: ErrorTolerances::new(
                phys_type::ConformalPerturbationTolerance::<2>::background_in_planck_units(1.0e-24),
                1.0e-13,
            ),
            ..self.settings
        };
        let target = self.system.hubble(&self.trajectory.state_at(self.output()));
        let start = MultifieldState::new(
            [Field::new(point[0]), Field::new(point[1])],
            [FieldVelocity::new(point[2]), FieldVelocity::new(point[3])],
            ScaleFactor::new(1.0),
        );
        let perturbed = self.system.survey(start, &precise).unwrap().trajectory;
        perturbed
            .crossing(|st| (self.system.hubble(st) - target).into_value())
            .expect("a perturbed trajectory ended before the target surface; move it earlier")
    }

    /// Phase-space delta N from the slice `n_c` for the triangle `shape`
    /// times `k`: returns (f_NL total, f_NL Gaussian part, P of the first
    /// canonical leg). `velocity_scale` is the relative step in the
    /// velocities.
    fn phase_space_delta_n(
        &self,
        k: f64,
        shape: [f64; 3],
        n_c: f64,
        velocity_scale: f64,
    ) -> (f64, f64, f64) {
        let base = self.trajectory.state_at(n_c);
        let x0 = [
            base.phi[0].into_value(),
            base.phi[1].into_value(),
            base.pi[0].into_value(),
            base.pi[1].into_value(),
        ];
        let h = [
            1.0e-4 * x0[0].abs(),
            (1.0e-2 * x0[1].abs()).min(1.0e-5),
            velocity_scale * x0[2].abs(),
            velocity_scale * x0[3].abs(),
        ];
        // The derivatives are taken in the coordinates (delta phi, delta
        // phidot) to which the linear `map` below relates the transported
        // (Q, V). Two nonlinear corrections turn such a point into the
        // initial state of a separate universe. The covariant Q is the
        // tangent of a geodesic, so delta phi = Q - Gamma Q Q / 2 at second
        // order. And the velocity is the proper-time one, phidot / N_l, with
        // the lapse N_l = H / H_loc fixed by the energy constraint
        // 3 H^2 = N_l^2 V + |phidot|^2 / 2 on the flat slice. Without them,
        // f_NL is off by up to 0.1% where it is large and by 1.4% for the
        // small values of the s_* = 0.1 branches.
        let connection =
            phys_type::FieldSpaceGeometry::christoffel(&self.system.geometry, &base.phi);
        let hubble = self.system.hubble(&base).into_value();
        let at = |shift: [f64; 4]| {
            let mut point: [f64; 4] = core::array::from_fn(|a| x0[a] + shift[a]);
            for (i, field) in point.iter_mut().take(2).enumerate() {
                *field -= 0.5
                    * (0..4)
                        .map(|m| {
                            let [j, kk] = [m / 2, m % 2];
                            connection.components()[i][j][kk].value() * shift[j] * shift[kk]
                        })
                        .sum::<f64>();
            }
            let state = |velocity: [f64; 2]| {
                MultifieldState::new(
                    [Field::new(point[0]), Field::new(point[1])],
                    velocity.map(FieldVelocity::new),
                    ScaleFactor::new(1.0),
                )
            };
            let potential = 3.0 * self.system.hubble(&state([0.0; 2])).into_value().powi(2);
            let kinetic = 3.0
                * self
                    .system
                    .hubble(&state([point[2], point[3]]))
                    .into_value()
                    .powi(2)
                - potential;
            let lapse = ((3.0 * hubble * hubble - kinetic) / potential).sqrt();
            point[2] /= lapse;
            point[3] /= lapse;
            self.e_folds_from_phase_space(point)
        };
        let unit = |a: usize, sign: f64| -> [f64; 4] {
            core::array::from_fn(|b| if a == b { sign * h[a] } else { 0.0 })
        };
        let add = |u: [f64; 4], v: [f64; 4]| -> [f64; 4] { core::array::from_fn(|a| u[a] + v[a]) };
        let centre = at([0.0; 4]);
        let plus: [f64; 4] = core::array::from_fn(|a| at(unit(a, 1.0)));
        let minus: [f64; 4] = core::array::from_fn(|a| at(unit(a, -1.0)));
        let first: [f64; 4] = core::array::from_fn(|a| (plus[a] - minus[a]) / (2.0 * h[a]));
        let mut second = [[0.0; 4]; 4];
        for a in 0..4 {
            second[a][a] = (plus[a] - 2.0 * centre + minus[a]) / (h[a] * h[a]);
            for b in 0..a {
                let value = (at(add(unit(a, 1.0), unit(b, 1.0)))
                    - at(add(unit(a, 1.0), unit(b, -1.0)))
                    - at(add(unit(a, -1.0), unit(b, 1.0)))
                    + at(add(unit(a, -1.0), unit(b, -1.0))))
                    / (4.0 * h[a] * h[b]);
                second[a][b] = value;
                second[b][a] = value;
            }
        }

        // Correlations on the slice from the transport, in the velocity basis.
        let bench =
            sourced_benchmark_with(self.system, self.initial, shape, Some(k), &self.settings);
        let state = bench.evolve(5.0, n_c, 1.0e-10);
        let velocity = phys_type::velocity_three_point_from_transport_state(
            &self.system,
            &base,
            bench.triangle,
            &state,
        )
        .unwrap();
        // Coordinate perturbations at linear order: delta phi = Q and
        // delta phidot = V - Gamma phidot Q. The second-order delta phi is in
        // `at` above; the second-order delta phidot depends on the nonlinear
        // definition of V and is left out (its effect here is below 1e-5).
        let mut map = [[0.0; 4]; 4];
        for (a, row) in map.iter_mut().enumerate() {
            row[a] = 1.0;
        }
        for i in 0..2 {
            for (kk, entry) in map[2 + i].iter_mut().take(2).enumerate() {
                *entry = -(0..2)
                    .map(|j| connection.components()[i][j][kk].value() * x0[2 + j])
                    .sum::<f64>();
            }
        }
        let sigma: [[[f64; 4]; 4]; 3] = core::array::from_fn(|leg| {
            core::array::from_fn(|a| {
                core::array::from_fn(|b| {
                    (0..16)
                        .map(|m| {
                            let [c, d] = [m / 4, m % 4];
                            map[a][c] * map[b][d] * covariance_entry(&state.covariances[leg], c, d)
                        })
                        .sum()
                })
            })
        });
        let alpha = |a: usize, b: usize, c: usize| -> f64 {
            (0..64)
                .map(|m| {
                    let [d, e, f] = [m / 16, (m / 4) % 4, m % 4];
                    map[a][d] * map[b][e] * map[c][f] * three_point_entry(&velocity, d, e, f)
                })
                .sum()
        };
        let power: [f64; 3] = core::array::from_fn(|leg| {
            (0..16)
                .map(|m| first[m / 4] * first[m % 4] * sigma[leg][m / 4][m % 4])
                .sum()
        });
        let projected: [[f64; 4]; 3] = core::array::from_fn(|leg| {
            core::array::from_fn(|c| (0..4).map(|a| first[a] * sigma[leg][a][c]).sum())
        });
        let gaussian: f64 = [(1, 2), (0, 2), (0, 1)]
            .into_iter()
            .map(|(u, v)| {
                (0..16)
                    .map(|m| projected[u][m / 4] * projected[v][m % 4] * second[m / 4][m % 4])
                    .sum::<f64>()
            })
            .sum();
        let intrinsic: f64 = (0..64)
            .map(|m| {
                let [a, b, c] = [m / 16, (m / 4) % 4, m % 4];
                first[a] * first[b] * first[c] * alpha(a, b, c)
            })
            .sum();
        let pairs = power[0] * power[1] + power[0] * power[2] + power[1] * power[2];
        let f = |b: f64| 5.0 / 6.0 * b / pairs;
        (f(gaussian + intrinsic), f(gaussian), power[0])
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_loop_transport_matches_phase_space_delta_n() {
    // A separate-universe check of the full tree-level result: phase-space
    // delta N on a superhorizon slice, with the evolved Gaussian covariance
    // and the transported field three-point function, must reproduce the
    // transport f_NL for any slice and shape.
    let model = calibrated_assisted_loop(0.1, 1.0e-3, 0.022);
    let k = model.k_pivot();
    let pivot = model.end - 60.0;
    for shape in [[1.0, 1.0, 1.0], [1.0, 1.0, 0.1]] {
        let bench =
            sourced_benchmark_with(model.system, model.initial, shape, Some(k), &model.settings);
        let transport = project_benchmark(&bench, 5.0, model.output()).reduced_f_nl();
        for offset in [3.0, 5.0] {
            let (total, gaussian, _) = model.phase_space_delta_n(k, shape, pivot + offset, 1.0e-2);
            assert!(
                (total - transport).abs() < 3.0e-4 * transport.abs(),
                "{shape:?} at pivot + {offset}: delta N {total} vs transport {transport}"
            );
            // The field three-point term is small on these slices.
            assert!((total - gaussian).abs() < 1.0e-2 * transport.abs());
        }
    }
}

/// The four (alpha, s_*) branches of Sec. 4.1 of arXiv:2302.05866, each
/// with the numerical coupling that gives n_s = 0.965 at the end of
/// inflation.
const ASSISTED_BRANCHES: [(f64, f64, f64); 4] = [
    (0.1, 1.0e-3, 0.0228),
    (0.1, 0.1, 0.0137),
    (1.0, 1.0e-3, 0.0214),
    (1.0, 0.1, 0.0124),
];

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_loop_branches_follow_eq_3_1_with_the_exact_s_e() {
    for (alpha, s_star, xi) in ASSISTED_BRANCHES {
        let model = calibrated_assisted_loop(alpha, s_star, xi);
        let numerical = 1.0 + model.mode_solver_tilt(model.output());
        let (analytic, _) = model.analytic_prediction();
        let (paper, _, paper_s_e) = model.paper_prescription();
        let s_e = model.trajectory.state_at(model.end).phi[1].into_value();
        assert!(
            (numerical - 0.965).abs() < 2.0e-3,
            "{alpha} {s_star}: n_s {numerical}"
        );
        assert!(
            (numerical - analytic).abs() < 1.5e-3,
            "{alpha} {s_star}: {numerical} vs Eq. (3.1) {analytic}"
        );
        // The slow-roll relation (3.9) overestimates s_e, which lowers n_s.
        assert!(
            paper_s_e > s_e && paper < numerical,
            "{alpha} {s_star}: {paper} vs {numerical}"
        );
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_loop_branches_transport_matches_phase_space_delta_n() {
    for (alpha, s_star, xi) in ASSISTED_BRANCHES {
        // A slightly earlier target keeps every perturbed trajectory
        // inflating until it reaches the common surface.
        let model = calibrated_assisted_loop_with_output(alpha, s_star, xi, 1.0e-2);
        let k = model.k_pivot();
        let bench = sourced_benchmark_with(
            model.system,
            model.initial,
            [1.0, 1.0, 1.0],
            Some(k),
            &model.settings,
        );
        let transport = project_benchmark(&bench, 5.0, model.output()).reduced_f_nl();
        let (delta_n, _, _) =
            model.phase_space_delta_n(k, [1.0, 1.0, 1.0], model.end - 57.0, 1.0e-2);
        assert!(
            (delta_n - transport).abs() < (3.0e-4 * transport.abs()).max(2.0e-5),
            "{alpha} {s_star}: delta N {delta_n} vs transport {transport}"
        );
    }
}

/// Checks one Class II branch of Sec. 4.2/4.3 of arXiv:2302.05866 at the
/// numerical coupling that gives n_s = 0.965 at the end of inflation.
fn check_class_two_branch<I: AssistedInflaton>(
    label: &str,
    inflaton: I,
    s_star: f64,
    xi: f64,
    paper_range: (f64, f64),
) {
    let model = calibrated_assisted(inflaton, s_star, xi, 1.0e-3);
    let end_state = model.trajectory.state_at(model.end);
    // Inflation ends on the field-value surface, not by slow-roll violation.
    assert!(model.system.epsilon_h(&end_state) < 1.0, "{label} {s_star}");
    let numerical = 1.0 + model.mode_solver_tilt(model.output());
    let (analytic, _) = model.analytic_prediction();
    assert!(
        (numerical - 0.965).abs() < 2.0e-3,
        "{label} {s_star}: n_s {numerical}"
    );
    assert!(
        (numerical - analytic).abs() < 1.5e-3,
        "{label} {s_star}: {numerical} vs Eq. (3.1) {analytic}"
    );
    assert!(
        paper_range.0 <= xi && xi <= paper_range.1,
        "{label} {s_star}: coupling {xi} outside {paper_range:?}"
    );
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_class_two_branches_follow_eq_3_1() {
    let power_law = |lambda| PowerLawInflaton { lambda };
    check_class_two_branch(
        "power-law",
        power_law(0.01),
        1.0e-3,
        0.0110,
        (0.0105, 0.0111),
    );
    check_class_two_branch("power-law", power_law(0.01), 0.1, 0.00516, (0.0045, 0.0056));
    check_class_two_branch(
        "power-law",
        power_law(0.05),
        1.0e-3,
        0.01353,
        (0.0129, 0.0136),
    );
    check_class_two_branch("power-law", power_law(0.05), 0.1, 0.00662, (0.0060, 0.0069));
    let hybrid = |mu, phi_e| HybridInflaton { mu, phi_e };
    check_class_two_branch(
        "hybrid",
        hybrid(10.0, 1.0),
        1.0e-3,
        0.01077,
        (0.0100, 0.0110),
    );
    check_class_two_branch("hybrid", hybrid(10.0, 1.0), 0.1, 0.00486, (0.0043, 0.0053));
    check_class_two_branch(
        "hybrid",
        hybrid(6.0, 0.06),
        1.0e-3,
        0.00558,
        (0.0053, 0.0059),
    );
    check_class_two_branch("hybrid", hybrid(6.0, 0.06), 0.1, 0.00340, (0.0022, 0.0043));
}

/// Transport f_NL of one Class II branch against phase-space delta N one
/// hundredth of an e-fold before the end; returns the end-of-inflation
/// transport value and the crossing-time massless delta N value.
fn check_class_two_f_nl<I: AssistedInflaton>(
    label: &str,
    inflaton: I,
    s_star: f64,
    xi: f64,
) -> (f64, f64) {
    let model = calibrated_assisted(inflaton, s_star, xi, 1.0e-3);
    let k = model.k_pivot();
    let bench = sourced_benchmark_with(
        model.system,
        model.initial,
        [1.0, 1.0, 1.0],
        Some(k),
        &model.settings,
    );
    let transport_end = project_benchmark(&bench, 5.0, model.output()).reduced_f_nl();
    let early = calibrated_assisted(inflaton, s_star, xi, 1.0e-2);
    let (_, _, crossing) = early.delta_n_at_crossing(early.end - 60.0);
    let early_bench = sourced_benchmark_with(
        early.system,
        early.initial,
        [1.0, 1.0, 1.0],
        Some(k),
        &early.settings,
    );
    let transport = project_benchmark(&early_bench, 5.0, early.output()).reduced_f_nl();
    let (delta_n, _, _) = early.phase_space_delta_n(k, [1.0, 1.0, 1.0], early.end - 57.0, 1.0e-2);
    assert!(
        (delta_n - transport).abs() < (3.0e-4 * transport.abs()).max(2.0e-5),
        "{label} {s_star}: delta N {delta_n} vs transport {transport}"
    );
    (transport_end, crossing)
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn assisted_class_two_transport_matches_phase_space_delta_n() {
    let (power_law, power_law_crossing) = check_class_two_f_nl(
        "power-law",
        PowerLawInflaton { lambda: 0.01 },
        1.0e-3,
        0.0110,
    );
    let (hybrid, hybrid_crossing) = check_class_two_f_nl(
        "hybrid",
        HybridInflaton {
            mu: 10.0,
            phi_e: 1.0,
        },
        1.0e-3,
        0.01077,
    );
    // With s_* = 1e-3 the tachyonic assistant makes crossing-time massless
    // delta N underestimate |f_NL| by about 8%.
    for (end, crossing) in [(power_law, power_law_crossing), (hybrid, hybrid_crossing)] {
        assert!(
            end < 0.0 && crossing.abs() < 0.95 * end.abs(),
            "{end} vs {crossing}"
        );
    }
    check_class_two_f_nl("power-law", PowerLawInflaton { lambda: 0.05 }, 0.1, 0.00662);
    check_class_two_f_nl(
        "hybrid",
        HybridInflaton {
            mu: 6.0,
            phi_e: 0.06,
        },
        0.1,
        0.00340,
    );
}
