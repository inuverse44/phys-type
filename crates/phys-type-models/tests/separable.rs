//! Three uncoupled quadratic fields against sum-separable delta N.
//!
//! `V = sum_i m_i^2 (phi^i)^2 / 2` on a flat field space is the case of
//! Battefeld and Easther, arXiv:astro-ph/0610296, Sec. III. Their Eq. (35)
//! gives the slow-roll delta N f_NL with the evolution between horizon
//! crossing and a later uniform-density slice; Eq. (3) adds the intrinsic
//! part. The tree-level transport must agree with phase-space delta N on a
//! superhorizon slice to the finite-difference accuracy, and with the
//! slow-roll formulas to first order in slow roll. See
//! `docs/n-field-perturbations.md`. All tests here are expensive; run them
//! with `cargo test --release -- --ignored`.

use std::f64::consts::PI;

use phys_type::{
    BackgroundTrajectory, CovariantMultifieldSystem, ErrorTolerances, Field, FieldVelocity,
    FlatFieldSpace, Mass, MomentumTriangle, MultifieldState, QuadraticMultifieldPotential,
    ScaleFactor, SpectrumRunSettings,
};

const FIELDS: usize = 3;
const PHASE: usize = 2 * FIELDS;

type SeparableSystem =
    CovariantMultifieldSystem<QuadraticMultifieldPotential<3>, FlatFieldSpace, 3>;

/// Squared masses in units of the lightest, `m_i^2 / m_1^2`. Battefeld and
/// Easther's two-field cases use ratios between 1 and 5.
const MASS_RATIOS: [f64; FIELDS] = [1.0, 2.0, 3.0];
const LIGHTEST_MASS: f64 = 1.0e-6;
/// The pivot leaves the horizon this many e-folds before the end.
const PIVOT_BEFORE_END: f64 = 55.0;

struct Separable {
    system: SeparableSystem,
    trajectory: BackgroundTrajectory<3>,
    end: f64,
    settings: SpectrumRunSettings,
}

/// Equal-energy initial conditions, Eq. (51): `phi^i = mu_i phi^1` with
/// `mu_i = m_1 / m_i`, at rest on the slow-roll attractor.
fn separable() -> Separable {
    let system = CovariantMultifieldSystem::in_planck_units(
        QuadraticMultifieldPotential::new(
            MASS_RATIOS.map(|ratio| Mass::new(LIGHTEST_MASS * ratio.sqrt())),
        ),
        FlatFieldSpace,
    );
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: PIVOT_BEFORE_END,
        k_min_over_pivot: 1.0,
        k_max_over_pivot: 1.0,
        k_count: 1,
        ..SpectrumRunSettings::default()
    };
    // About 65 e-folds of slow roll, Eq. (52): 4 N = sum_i (phi^i)^2.
    let phi1 = (4.0 * 65.0 / MASS_RATIOS.iter().map(|ratio| 1.0 / ratio).sum::<f64>()).sqrt();
    let initial = slow_roll_start(&system, MASS_RATIOS.map(|ratio| phi1 / ratio.sqrt()));
    let survey = system.survey(initial, &settings).unwrap();
    Separable {
        system,
        trajectory: survey.trajectory,
        end: survey.e_folds_at_end,
        settings,
    }
}

fn slow_roll_start(system: &SeparableSystem, phi: [f64; FIELDS]) -> MultifieldState<3> {
    let resting = MultifieldState::new(
        phi.map(Field::new),
        [FieldVelocity::new(0.0); FIELDS],
        ScaleFactor::new(1.0),
    );
    MultifieldState {
        pi: system.slow_roll_velocity(&resting),
        ..resting
    }
}

fn potential_terms(phi: [f64; FIELDS]) -> [f64; FIELDS] {
    core::array::from_fn(|i| 0.5 * MASS_RATIOS[i] * phi[i] * phi[i])
}

/// The slow-roll delta N of Battefeld and Easther from the fields at
/// horizon crossing, `star`, to the uniform-density slice through `slice`:
/// returns `sum_i N_i^2` and the delta N term of their Eq. (3),
/// `N_i N_j N_ij / (N_k N_k)^2`, from Eqs. (29) and (34)--(37).
/// The common mass scale cancels, so the potentials are in units of `m_1^2`.
fn battefeld_easther(star: [f64; FIELDS], slice: [f64; FIELDS]) -> (f64, f64) {
    let v_star = potential_terms(star);
    let v_c = potential_terms(slice);
    let w_star: f64 = v_star.iter().sum();
    let w_c: f64 = v_c.iter().sum();
    // Eq. (11): epsilon_i = V_i'^2 / (2 W^2) and eta_i = V_i'' / W.
    let epsilon = |phi: [f64; FIELDS], w: f64| -> [f64; FIELDS] {
        core::array::from_fn(|i| 0.5 * (MASS_RATIOS[i] * phi[i] / w).powi(2))
    };
    let eta = |w: f64| -> [f64; FIELDS] { core::array::from_fn(|i| MASS_RATIOS[i] / w) };
    let (eps_star, eta_star) = (epsilon(star, w_star), eta(w_star));
    let (eps_c, eta_c) = (epsilon(slice, w_c), eta(w_c));
    let eps_c_total: f64 = eps_c.iter().sum();
    // Eq. (37).
    let u: [f64; FIELDS] = core::array::from_fn(|k| {
        (v_star[k] - v_c[k]) / w_star + w_c / w_star * eps_c[k] / eps_c_total
    });
    // Eqs. (33)--(34).
    let delta = |a: usize, b: usize| if a == b { 1.0 } else { 0.0 };
    let a_matrix: [[f64; FIELDS]; FIELDS] = core::array::from_fn(|l| {
        core::array::from_fn(|k| {
            -(w_c / w_star).powi(2)
                * (0..FIELDS)
                    .map(|j| {
                        eps_c[j]
                            * (eps_c[l] / eps_c_total - delta(l, j))
                            * (eps_c[k] / eps_c_total - delta(k, j))
                            * (1.0 - eta_c[j] / eps_c_total)
                    })
                    .sum::<f64>()
        })
    });
    let norm: f64 = (0..FIELDS).map(|k| u[k] * u[k] / eps_star[k]).sum();
    let diagonal: f64 = (0..FIELDS)
        .map(|k| u[k] * u[k] / eps_star[k] * (1.0 - eta_star[k] * u[k] / (2.0 * eps_star[k])))
        .sum();
    let mixing: f64 = (0..FIELDS * FIELDS)
        .map(|m| {
            let [k, l] = [m / FIELDS, m % FIELDS];
            u[k] * u[l] / (eps_star[k] * eps_star[l]) * a_matrix[l][k]
        })
        .sum();
    // Eq. (29): N_k = u_k / sqrt(2 epsilon_k). The factor 2 of Eq. (35)
    // multiplies the whole numerator, as in its two-field form, Eq. (40).
    (0.5 * norm, 2.0 * (diagonal + mixing) / (norm * norm))
}

/// The point of the slow-roll trajectory through `star` with potential
/// `w_c`. The integrals of motion, Eq. (50), make `(phi^i / phi^i_*)^{1/mu_i^2}`
/// equal for all `i`; with that common value `t`, `phi^i = phi^i_* t^{mu_i^2}`.
fn slow_roll_point(star: [f64; FIELDS], w_c: f64) -> [f64; FIELDS] {
    let at =
        |t: f64| -> [f64; FIELDS] { core::array::from_fn(|i| star[i] * t.powf(MASS_RATIOS[i])) };
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        if potential_terms(at(mid)).iter().sum::<f64>() < w_c {
            low = mid;
        } else {
            high = mid;
        }
    }
    at(0.5 * (low + high))
}

/// The slow-roll e-fold number from `star` to the potential `w_c`, Eq. (52).
fn slow_roll_e_folds(star: [f64; FIELDS], w_c: f64) -> f64 {
    let slice = slow_roll_point(star, w_c);
    (0..FIELDS)
        .map(|i| (star[i] * star[i] - slice[i] * slice[i]) / 4.0)
        .sum()
}

fn slow_roll_delta_n_by_differences(star: [f64; FIELDS], slice: [f64; FIELDS]) -> (f64, f64) {
    let w_c: f64 = potential_terms(slice).iter().sum();
    let h = 1.0e-3;
    let at =
        |shift: [f64; FIELDS]| slow_roll_e_folds(core::array::from_fn(|i| star[i] + shift[i]), w_c);
    let unit = |a: usize, sign: f64| -> [f64; FIELDS] {
        core::array::from_fn(|b| if a == b { sign * h } else { 0.0 })
    };
    let centre = at([0.0; FIELDS]);
    let first: [f64; FIELDS] =
        core::array::from_fn(|a| (at(unit(a, 1.0)) - at(unit(a, -1.0))) / (2.0 * h));
    let second: [[f64; FIELDS]; FIELDS] = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            if a == b {
                (at(unit(a, 1.0)) - 2.0 * centre + at(unit(a, -1.0))) / (h * h)
            } else {
                let corner = |sa: f64, sb: f64| {
                    at(core::array::from_fn(|c| unit(a, sa)[c] + unit(b, sb)[c]))
                };
                (corner(1.0, 1.0) - corner(1.0, -1.0) - corner(-1.0, 1.0) + corner(-1.0, -1.0))
                    / (4.0 * h * h)
            }
        })
    });
    let norm: f64 = first.iter().map(|n| n * n).sum();
    let cubic: f64 = (0..FIELDS * FIELDS)
        .map(|m| {
            let [a, b] = [m / FIELDS, m % FIELDS];
            first[a] * first[b] * second[a][b]
        })
        .sum();
    (norm, cubic / (norm * norm))
}

fn fields(state: &MultifieldState<3>) -> [f64; FIELDS] {
    state.phi.map(|phi| phi.into_value())
}

impl Separable {
    fn k_pivot(&self) -> f64 {
        let pivot = self.trajectory.state_at(self.end - PIVOT_BEFORE_END);
        (pivot.a * self.system.hubble(&pivot)).into_value()
    }

    /// The background where the pivot leaves the horizon, `k = a H`.
    fn crossing(&self) -> MultifieldState<3> {
        let k = self.k_pivot();
        let n = self
            .trajectory
            .crossing(|st| k - (st.a * self.system.hubble(st)).into_value())
            .unwrap();
        self.trajectory.state_at(n)
    }

    /// Transport `P_R` and equilateral `f_NL` of the pivot on the slice `slice`.
    fn transport(&self, slice: f64) -> (f64, f64) {
        let k = self.k_pivot();
        let power = self
            .system
            .mode_spectrum(&self.settings, &self.trajectory, Mass::new(k), slice)
            .unwrap()
            .spectra
            .curvature
            .into_value();
        let triangle = MomentumTriangle::try_new(k, k, k).unwrap();
        let state = self
            .system
            .three_point_state(&self.trajectory, triangle, 5.0, slice, 1.0e-9, 2_000_000)
            .unwrap();
        let f_nl = phys_type::curvature_correlators(
            &self.system,
            &self.trajectory.state_at(slice),
            triangle,
            &state,
            0.1,
        )
        .unwrap()
        .reduced_f_nl();
        (power, f_nl)
    }

    /// The slow-roll prediction for the slice `slice` in the Planck sign:
    /// `P_R = (H_*/2 pi)^2 N_i N_i` and, from Eq. (3) with `r = 8 / N_i N_i`
    /// and the equilateral shape factor `f = 5/6`,
    /// `f_NL = (5/6) [r (1 + f) / 16 + N_i N_j N_ij / (N_k N_k)^2]`.
    /// Returns `(P_R, f_NL, delta N term)`.
    fn slow_roll(&self, slice: f64) -> (f64, f64, f64) {
        let star = self.crossing();
        let (norm, delta_n) =
            battefeld_easther(fields(&star), fields(&self.trajectory.state_at(slice)));
        let hubble = self.system.hubble(&star).into_value();
        let power = (hubble / (2.0 * PI)).powi(2) * norm;
        let r = 8.0 / norm;
        let f_nl = 5.0 / 6.0 * (r * (1.0 + 5.0 / 6.0) / 16.0 + delta_n);
        (power, f_nl, delta_n)
    }

    /// e-folds from a phase-space point to the uniform-density slice through
    /// the background at `slice`.
    fn e_folds_from_phase_space(&self, point: [f64; PHASE], slice: f64) -> f64 {
        let precise = SpectrumRunSettings {
            tolerances: ErrorTolerances::new(
                phys_type::ConformalPerturbationTolerance::<2>::background_in_planck_units(1.0e-24),
                1.0e-13,
            ),
            ..self.settings
        };
        let target = self.system.hubble(&self.trajectory.state_at(slice));
        let start = MultifieldState::new(
            core::array::from_fn(|i| Field::new(point[i])),
            core::array::from_fn(|i| FieldVelocity::new(point[FIELDS + i])),
            ScaleFactor::new(1.0),
        );
        let perturbed = self.system.survey(start, &precise).unwrap().trajectory;
        perturbed
            .crossing(|st| (self.system.hubble(st) - target).into_value())
            .expect("a perturbed trajectory ended before the target surface; move it earlier")
    }

    /// Phase-space delta N for the equilateral pivot triangle from the
    /// superhorizon slice `n_c` to the uniform-density slice through
    /// `slice`, with the transported covariance and field three-point
    /// function on `n_c`: returns (f_NL total, f_NL Gaussian part, P_R).
    fn phase_space_delta_n(&self, n_c: f64, slice: f64) -> (f64, f64, f64) {
        let base = self.trajectory.state_at(n_c);
        let x0: [f64; PHASE] = core::array::from_fn(|a| {
            if a < FIELDS {
                base.phi[a].into_value()
            } else {
                base.pi[a - FIELDS].into_value()
            }
        });
        let h: [f64; PHASE] =
            core::array::from_fn(|a| if a < FIELDS { 1.0e-4 } else { 1.0e-2 } * x0[a].abs());
        // The transport perturbs the coordinate-time velocity on a flat
        // slice; a separate universe starts from the proper-time velocity.
        // On a flat slice the local scale factor is the background one, so
        // the lapse is N_l = H / H_loc, and the energy constraint
        // 3 H^2 / N_l^2 = V + |phidot|^2 / (2 N_l^2) fixes it nonlinearly.
        // Without it, P_R and f_NL are both off by 0.8% here.
        let hubble = self.system.hubble(&base).into_value();
        let at = |shift: [f64; PHASE]| {
            let mut point: [f64; PHASE] = core::array::from_fn(|a| x0[a] + shift[a]);
            let phi: [f64; FIELDS] = core::array::from_fn(|i| point[i]);
            let potential: f64 =
                potential_terms(phi).iter().sum::<f64>() * LIGHTEST_MASS * LIGHTEST_MASS;
            let speed2: f64 = point[FIELDS..].iter().map(|v| v * v).sum();
            let lapse = ((hubble * hubble - speed2 / 6.0) / (potential / 3.0)).sqrt();
            for velocity in &mut point[FIELDS..] {
                *velocity /= lapse;
            }
            self.e_folds_from_phase_space(point, slice)
        };
        let unit = |a: usize, sign: f64| -> [f64; PHASE] {
            core::array::from_fn(|b| if a == b { sign * h[a] } else { 0.0 })
        };
        let add = |u: [f64; PHASE], v: [f64; PHASE]| -> [f64; PHASE] {
            core::array::from_fn(|a| u[a] + v[a])
        };
        let centre = at([0.0; PHASE]);
        let plus: [f64; PHASE] = core::array::from_fn(|a| at(unit(a, 1.0)));
        let minus: [f64; PHASE] = core::array::from_fn(|a| at(unit(a, -1.0)));
        let first: [f64; PHASE] = core::array::from_fn(|a| (plus[a] - minus[a]) / (2.0 * h[a]));
        let mut second = [[0.0; PHASE]; PHASE];
        for a in 0..PHASE {
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

        // Correlations on the slice from the transport, in the velocity
        // basis, which on a flat field space are (delta phi, delta phidot).
        let k = self.k_pivot();
        let triangle = MomentumTriangle::try_new(k, k, k).unwrap();
        let state = self
            .system
            .three_point_state(&self.trajectory, triangle, 5.0, n_c, 1.0e-10, 2_000_000)
            .unwrap();
        let alpha = phys_type::velocity_three_point_from_transport_state(
            &self.system,
            &base,
            triangle,
            &state,
        )
        .unwrap();
        let sigma = |leg: usize, a: usize, b: usize| -> f64 {
            let covariance = &state.covariances[leg];
            let (i, j) = (a % FIELDS, b % FIELDS);
            match (a < FIELDS, b < FIELDS) {
                (true, true) => covariance.qq.components()[i][j].value().re,
                (true, false) => covariance.qp.components()[i][j].value().re,
                (false, true) => covariance.pq.components()[i][j].value().re,
                (false, false) => covariance.pp.components()[i][j].value().re,
            }
        };
        let three_point = |a: usize, b: usize, c: usize| -> f64 {
            let (i, j, k) = (a % FIELDS, b % FIELDS, c % FIELDS);
            match (a < FIELDS, b < FIELDS, c < FIELDS) {
                (true, true, true) => alpha.qqq.components()[i][j][k].value().re,
                (true, true, false) => alpha.qqp.components()[i][j][k].value().re,
                (true, false, true) => alpha.qpq.components()[i][j][k].value().re,
                (false, true, true) => alpha.pqq.components()[i][j][k].value().re,
                (true, false, false) => alpha.qpp.components()[i][j][k].value().re,
                (false, true, false) => alpha.pqp.components()[i][j][k].value().re,
                (false, false, true) => alpha.ppq.components()[i][j][k].value().re,
                (false, false, false) => alpha.ppp.components()[i][j][k].value().re,
            }
        };
        let pairs_of = |n: usize| (0..n * n).map(move |m| [m / n, m % n]);
        let power: [f64; 3] = core::array::from_fn(|leg| {
            pairs_of(PHASE)
                .map(|[a, b]| first[a] * first[b] * sigma(leg, a, b))
                .sum()
        });
        let projected: [[f64; PHASE]; 3] = core::array::from_fn(|leg| {
            core::array::from_fn(|c| (0..PHASE).map(|a| first[a] * sigma(leg, a, c)).sum())
        });
        let gaussian: f64 = [(1, 2), (0, 2), (0, 1)]
            .into_iter()
            .map(|(u, v)| {
                pairs_of(PHASE)
                    .map(|[a, b]| projected[u][a] * projected[v][b] * second[a][b])
                    .sum::<f64>()
            })
            .sum();
        let intrinsic: f64 = (0..PHASE * PHASE * PHASE)
            .map(|m| {
                let [a, b, c] = [m / (PHASE * PHASE), (m / PHASE) % PHASE, m % PHASE];
                first[a] * first[b] * first[c] * three_point(a, b, c)
            })
            .sum();
        let pairs = power[0] * power[1] + power[0] * power[2] + power[1] * power[2];
        let f = |b: f64| 5.0 / 6.0 * b / pairs;
        let dimensionless = k.powi(3) * power[0] / (2.0 * PI * PI);
        (f(gaussian + intrinsic), f(gaussian), dimensionless)
    }
}

#[test]
fn eq_35_matches_differences_of_the_slow_roll_e_folds() {
    // A check of the transcription of Eqs. (29) and (34)--(37): the same
    // derivatives by finite differences of the slow-roll e-fold number.
    let star = [12.0, 9.0, 7.0];
    let slice = [8.0, 4.0, 2.0];
    // Put the slice on the slow-roll trajectory through `star`.
    let on_trajectory = slow_roll_point(star, potential_terms(slice).iter().sum());
    let (norm, delta_n) = battefeld_easther(star, on_trajectory);
    let (norm_fd, delta_n_fd) = slow_roll_delta_n_by_differences(star, on_trajectory);
    assert!((norm / norm_fd - 1.0).abs() < 1.0e-8, "{norm} vs {norm_fd}");
    assert!(
        (delta_n / delta_n_fd - 1.0).abs() < 1.0e-6,
        "{delta_n} vs {delta_n_fd}"
    );
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn separable_transport_follows_battefeld_easther() {
    // Slow roll fixes the amplitudes at horizon crossing to O(epsilon_*),
    // about 1%; the superhorizon evolution, which the A-matrix term of
    // Eq. (35) carries, is tested by the change between slices.
    let model = separable();
    let slices = [40.0, 20.0, 5.0].map(|before_end| model.end - before_end);
    let exact = slices.map(|slice| model.transport(slice));
    let slow_roll = slices.map(|slice| model.slow_roll(slice));
    for ((power, f_nl), (power_sr, f_nl_sr, _)) in exact.into_iter().zip(slow_roll) {
        assert!(
            (power / power_sr - 1.0).abs() < 2.0e-2,
            "P_R {power:e} vs {power_sr:e}"
        );
        assert!(
            (f_nl / f_nl_sr - 1.0).abs() < 4.0e-2,
            "f_NL {f_nl} vs {f_nl_sr}"
        );
    }
    let growth = exact[2].0 / exact[0].0;
    let growth_sr = slow_roll[2].0 / slow_roll[0].0;
    assert!(growth > 1.1, "isocurvature conversion {growth}");
    assert!(
        (growth / growth_sr - 1.0).abs() < 5.0e-3,
        "P_R growth {growth} vs {growth_sr}"
    );
    let change = exact[2].1 - exact[0].1;
    let change_sr = slow_roll[2].1 - slow_roll[0].1;
    assert!(
        (change / change_sr - 1.0).abs() < 5.0e-2,
        "f_NL change {change} vs {change_sr}"
    );
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn separable_transport_matches_phase_space_delta_n() {
    // A separate-universe check of the full tree-level result: the field
    // three-point term is more than half of f_NL here, so the comparison
    // tests the transported Alpha as well as the covariance. Both sides use
    // the same start depth, whose own error (0.6% at five e-folds) cancels.
    let model = separable();
    let pivot = model.end - PIVOT_BEFORE_END;
    let slice = model.end - 20.0;
    let (power, transport) = model.transport(slice);
    for (offset, tolerance) in [(5.0, 1.0e-4), (7.0, 1.0e-5)] {
        // Gradient terms fall as (k/aH)^2 = e^{-2 offset}.
        let (total, gaussian, power_dn) = model.phase_space_delta_n(pivot + offset, slice);
        assert!(
            (total / transport - 1.0).abs() < tolerance,
            "pivot + {offset}: delta N {total} vs transport {transport}"
        );
        assert!(
            (power_dn / power - 1.0).abs() < tolerance,
            "pivot + {offset}: P_R {power_dn:e} vs {power:e}"
        );
        assert!((total - gaussian).abs() > 0.5 * total.abs());
    }
}
