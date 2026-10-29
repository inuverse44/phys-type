//! Probe for issue #67: where does the tree-level transport f_NL of KMS
//! trajectory 2 lose precision, and what does phase-space delta N from an
//! earlier slice give? `<index>` is 0, 1, 2 for trajectories 1--3; slices
//! are given in e-folds before the end of inflation.
//!
//! - `transport <index> <tol,...> [<slice,...>]`: P_zeta and f_NL of the
//!   equilateral pivot triangle on a list of slices (decreasing), once per
//!   relative tolerance.
//! - `delta-n <index> <tol> <target> <n_c,...> <phi,chi,velocity>`:
//!   phase-space delta N from each slice `n_c` to the surface of uniform
//!   `H` at `target`, with the transport (tolerance `tol`) up to `n_c`; the
//!   last argument gives the finite-difference steps relative to each
//!   coordinate. The last row is the transport itself at `target`.
//! - `response <index> <tol> <target> <n_c,...> <step>`: the response of
//!   zeta on the surface `target` to the phase space on a slice, integrated
//!   backwards from `target` with RK4 steps of at most `step` e-folds,
//!   combined with the transport (tolerance `tol`) up to each slice `n_c`.
//!   This is the delta N split without finite differences or the
//!   coordinate map: `n_a` and `n_ab` obey the adjoint of the transport's
//!   own flows `u^a_b` and `u^a_bc`.
//! - `noise <index> <n_c> <target>`: `N` along `chi` from `n_c`, to the `H`
//!   surface and to `epsilon_H = 1`, to measure the noise of the endpoint.
//!
//! `KMS_BACKGROUND_TOLERANCE` sets the relative tolerance of the background
//! survey (default `1e-10`); the response needs `1e-13` or below.
//!
//! ```sh
//! cargo run --release --example kms_transport_probe -- transport 1 1e-10,1e-13
//! cargo run --release --example kms_transport_probe -- delta-n 1 1e-12 0.01 50 1e-4,1e-3,1e-2
//! KMS_BACKGROUND_TOLERANCE=1e-13 \
//!     cargo run --release --example kms_transport_probe -- response 1 1e-12 0.01 10,50 1.25e-4
//! ```

use phys_type::{
    BackgroundNeighborhood, BackgroundTrajectory, CanonicalThreePointTransportState,
    CanonicalThreePointTransportSystem, CanonicalThreePointTransportTolerance,
    ConformalPerturbationTolerance, Dimensionless, ErrorTolerances, Field, FieldSpaceGeometry,
    FieldVelocity, Formulation, MomentumTriangle, MultifieldState, ScaleFactor,
    SpectrumRunSettings, StepSizeLimits, adaptive_step, canonical_cubic_hamiltonian_vertices,
    canonical_three_point_transport_initial_state, curvature_correlators,
    curvature_gauge_coefficients, next_step_size, velocity_three_point_from_transport_state,
};
use phys_type_models::kms::{KmsSystem, PHI_START, TRAJECTORIES, at_rest, benchmark_system};

struct Transport<'a> {
    system: &'a KmsSystem,
    trajectory: &'a BackgroundTrajectory<2>,
    triangle: MomentumTriangle,
    state: CanonicalThreePointTransportState<2>,
    e_folds: f64,
    step: Dimensionless,
    tolerance: f64,
}

impl<'a> Transport<'a> {
    fn new(
        system: &'a KmsSystem,
        trajectory: &'a BackgroundTrajectory<2>,
        k: f64,
        tolerance: f64,
    ) -> Self {
        let triangle = MomentumTriangle::try_new(k, k, k).unwrap();
        let depth = 5.0_f64.exp();
        let start = trajectory
            .crossing(|s| k / (s.a * system.hubble(s)).into_value() - depth)
            .unwrap();
        let state = canonical_three_point_transport_initial_state(
            system,
            &BackgroundNeighborhood::from_trajectory(trajectory, start, 1.0e-3),
            triangle,
            depth * (1.0 - 1.0e-6),
        )
        .unwrap();
        Self {
            system,
            trajectory,
            triangle,
            state,
            e_folds: start,
            step: Dimensionless::new(1.0e-3),
            tolerance,
        }
    }

    fn advance(&mut self, target: f64) {
        let transport =
            CanonicalThreePointTransportSystem::new(self.system, self.trajectory, self.triangle);
        let tolerances = ErrorTolerances::new(
            CanonicalThreePointTransportTolerance::pure_relative(),
            self.tolerance,
        );
        let limits = StepSizeLimits {
            min: Dimensionless::new(1.0e-9),
            max: Dimensionless::new(5.0e-2),
        };
        while target - self.e_folds > 1.0e-13 {
            let attempted = Dimensionless::new(self.step.into_value().min(target - self.e_folds));
            let outcome = adaptive_step(
                &transport,
                Dimensionless::new(self.e_folds),
                &self.state,
                attempted,
                &tolerances,
            );
            if outcome.accepted {
                self.state = outcome.state;
                self.e_folds += attempted.into_value();
            }
            self.step = next_step_size(attempted, outcome.error_ratio, &limits);
        }
    }

    fn f_nl(&self) -> (f64, f64) {
        let background = self.trajectory.state_at(self.e_folds);
        let p = curvature_correlators(self.system, &background, self.triangle, &self.state, 0.1)
            .unwrap();
        (p.reduced_f_nl(), p.dimensionless_power()[0])
    }
}

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

/// Phase-space delta N from the transport state at its current slice to the
/// surface `H = target_hubble`: (f_NL, f_NL Gaussian part, P_zeta).
fn phase_space_delta_n(
    transport: &Transport,
    settings: &SpectrumRunSettings,
    target_hubble: f64,
    scale: [f64; 3],
) -> (f64, f64, f64) {
    let system = transport.system;
    let base = transport.trajectory.state_at(transport.e_folds);
    let x0 = [
        base.phi[0].into_value(),
        base.phi[1].into_value(),
        base.pi[0].into_value(),
        base.pi[1].into_value(),
    ];
    let h = [
        scale[0] * x0[0].abs(),
        scale[1] * x0[1].abs(),
        scale[2] * x0[2].abs(),
        scale[2] * x0[3].abs(),
    ];
    let precise = SpectrumRunSettings {
        tolerances: ErrorTolerances::new(
            ConformalPerturbationTolerance::<2>::background_in_planck_units(1.0e-24),
            1.0e-13,
        ),
        background_sample_every: usize::MAX,
        ..*settings
    };
    let connection = system.geometry.christoffel(&base.phi);
    let hubble = system.hubble(&base).into_value();
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
        let potential = 3.0 * system.hubble(&state([0.0; 2])).into_value().powi(2);
        let kinetic = 3.0
            * system
                .hubble(&state([point[2], point[3]]))
                .into_value()
                .powi(2)
            - potential;
        let lapse = ((3.0 * hubble * hubble - kinetic) / potential).sqrt();
        point[2] /= lapse;
        point[3] /= lapse;
        let start = MultifieldState::new(
            [Field::new(point[0]), Field::new(point[1])],
            [FieldVelocity::new(point[2]), FieldVelocity::new(point[3])],
            ScaleFactor::new(1.0),
        );
        let perturbed = system.survey(start, &precise).unwrap().trajectory;
        perturbed
            .crossing(|st| system.hubble(st).into_value() - target_hubble)
            .expect("a perturbed trajectory ended before the target surface")
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

    let state = &transport.state;
    let velocity =
        velocity_three_point_from_transport_state(system, &base, transport.triangle, state)
            .unwrap();
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
    let k = transport.triangle.canonical_wavenumbers()[0];
    (
        f(gaussian + intrinsic),
        f(gaussian),
        k.powi(3) * power[0] / (2.0 * core::f64::consts::PI.powi(2)),
    )
}

/// `u^a_b` on the canonical phase space `(Q^1, Q^2, P^1, P^2)`.
type LinearFlow = [[f64; 4]; 4];
/// `u^a_bc` for output legs 1, 2, 3, inputs in ascending leg order.
type QuadraticFlows = [[[[f64; 4]; 4]; 4]; 3];

/// The block `[i][j][l]` of a quadratic kernel as plain numbers.
fn block<A, B, C, D, const N: usize>(
    kernel: &phys_type::ThreePointBlock<A, B, C, D, f64, N>,
) -> [[[f64; N]; N]; N] {
    let components = kernel.components();
    core::array::from_fn(|i| {
        core::array::from_fn(|j| core::array::from_fn(|l| *components[i][j][l].value()))
    })
}

/// The linear flow and the three quadratic kernels on the background at
/// `e_folds`, the same ones the canonical transport uses.
fn flows(
    system: &KmsSystem,
    trajectory: &BackgroundTrajectory<2>,
    triangle: MomentumTriangle,
    e_folds: f64,
) -> (LinearFlow, QuadraticFlows) {
    let background = trajectory.state_at(e_folds);
    let hubble = system.hubble(&background);
    let h = hubble.into_value();
    let connection = system.geometry.christoffel(&background.phi);
    let mass = system.perturbation_mass_matrix(&background);
    let k = triangle.canonical_wavenumbers()[0] / background.a.into_value();
    let mut u = [[0.0; 4]; 4];
    for i in 0..2 {
        for j in 0..2 {
            let flow: f64 = (0..2)
                .map(|v| connection.components()[i][v][j].value() * background.pi[v].value() / h)
                .sum();
            let delta = f64::from(i == j);
            u[i][j] = -flow;
            u[i][2 + j] = delta / h;
            u[2 + i][j] = -(mass.components()[i][j].value() + delta * k * k) / h;
            u[2 + i][2 + j] = -3.0 * delta - flow;
        }
    }
    let kernels = canonical_cubic_hamiltonian_vertices(system, &background, triangle)
        .unwrap()
        .e_fold_kernels(hubble)
        .unwrap();
    let quadratic = kernels.map(|kernel| {
        // Blocks in the order QQQ, QQP, QPQ, QPP, PQQ, PQP, PPQ, PPP.
        let blocks = [
            block(&kernel.q_qq),
            block(&kernel.q_qp),
            block(&kernel.q_pq),
            block(&kernel.q_pp),
            block(&kernel.p_qq),
            block(&kernel.p_qp),
            block(&kernel.p_pq),
            block(&kernel.p_pp),
        ];
        core::array::from_fn(|a| {
            core::array::from_fn(|b| {
                core::array::from_fn(|c| {
                    blocks[4 * (a / 2) + 2 * (b / 2) + c / 2][a % 2][b % 2][c % 2]
                })
            })
        })
    });
    (u, quadratic)
}

/// The response of the final zeta to the phase space at a slice:
/// `zeta_e = n_a X^a + (1/2) n^(r)_ab X^a X^b`, one quadratic part per output leg.
#[derive(Clone, Copy)]
struct Response {
    linear: [f64; 4],
    quadratic: [[[f64; 4]; 4]; 3],
}

impl Response {
    fn axpy(&self, scale: f64, other: &Self) -> Self {
        Self {
            linear: core::array::from_fn(|a| self.linear[a] + scale * other.linear[a]),
            quadratic: core::array::from_fn(|r| {
                core::array::from_fn(|a| {
                    core::array::from_fn(|b| {
                        self.quadratic[r][a][b] + scale * other.quadratic[r][a][b]
                    })
                })
            }),
        }
    }

    /// `d n / dN` from the adjoint of the flow `dX = u X + (1/2) u_bc X X`.
    fn derivative(&self, u: &LinearFlow, kernels: &QuadraticFlows) -> Self {
        Self {
            linear: core::array::from_fn(|a| {
                -(0..4).map(|b| self.linear[b] * u[b][a]).sum::<f64>()
            }),
            quadratic: core::array::from_fn(|r| {
                core::array::from_fn(|a| {
                    core::array::from_fn(|b| {
                        -(0..4)
                            .map(|c| {
                                self.linear[c] * kernels[r][c][a][b]
                                    + self.quadratic[r][c][b] * u[c][a]
                                    + self.quadratic[r][a][c] * u[c][b]
                            })
                            .sum::<f64>()
                    })
                })
            }),
        }
    }
}

/// f_NL and the dimensionless P_zeta from the response on a slice and the
/// transport state there.
fn split_f_nl(response: &Response, transport: &Transport) -> (f64, f64) {
    let state = &transport.state;
    let n = &response.linear;
    let sigma: [[[f64; 4]; 4]; 3] = core::array::from_fn(|leg| {
        core::array::from_fn(|a| {
            core::array::from_fn(|b| covariance_entry(&state.covariances[leg], a, b))
        })
    });
    let power: [f64; 3] = core::array::from_fn(|leg| {
        (0..16)
            .map(|m| n[m / 4] * n[m % 4] * sigma[leg][m / 4][m % 4])
            .sum()
    });
    let projected: [[f64; 4]; 3] = core::array::from_fn(|leg| {
        core::array::from_fn(|a| (0..4).map(|c| n[c] * sigma[leg][c][a]).sum())
    });
    let intrinsic: f64 = (0..64)
        .map(|m| {
            let [a, b, c] = [m / 16, (m / 4) % 4, m % 4];
            n[a] * n[b] * n[c] * canonical_entry(&state.alpha, a, b, c)
        })
        .sum();
    let gauge: f64 = [(0, 1, 2), (1, 0, 2), (2, 0, 1)]
        .into_iter()
        .map(|(r, s, t)| {
            (0..16)
                .map(|m| {
                    response.quadratic[r][m / 4][m % 4] * projected[s][m / 4] * projected[t][m % 4]
                })
                .sum::<f64>()
        })
        .sum();
    let pairs = power[0] * power[1] + power[0] * power[2] + power[1] * power[2];
    let k = transport.triangle.canonical_wavenumbers()[0];
    (
        5.0 / 6.0 * (intrinsic + gauge) / pairs,
        k.powi(3) * power[0] / (2.0 * core::f64::consts::PI.powi(2)),
    )
}

fn canonical_entry(
    alpha: &phys_type::CanonicalThreePointPhaseSpaceCorrelation<2>,
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args[1].as_str();
    let index: usize = args[2].parse().unwrap();
    let list = |s: &str| -> Vec<f64> { s.split(',').map(|x| x.parse().unwrap()).collect() };
    let system = benchmark_system(Formulation::Metric);
    // KMS_BACKGROUND_TOLERANCE overrides the relative tolerance of the survey.
    let background_tolerance =
        std::env::var("KMS_BACKGROUND_TOLERANCE").map_or(1.0e-10, |value| value.parse().unwrap());
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: 60.0,
        e_folds_before_exit: 6.0,
        background_sample_every: 1,
        tolerances: ErrorTolerances::new(
            ConformalPerturbationTolerance::<2>::background_in_planck_units(
                1.0e-10 * background_tolerance,
            ),
            background_tolerance,
        ),
        ..SpectrumRunSettings::default()
    };
    let survey = system
        .survey(at_rest(PHI_START, TRAJECTORIES[index]), &settings)
        .unwrap();
    let end = survey.e_folds_at_end;
    let trajectory = &survey.trajectory;
    let pivot = trajectory.state_at(end - 60.0);
    let k = (pivot.a * system.hubble(&pivot)).into_value();
    match mode {
        "transport" => {
            let slices = args.get(4).map_or(
                vec![
                    57.0, 50.0, 40.0, 20.0, 10.0, 8.0, 7.0, 6.5, 6.0, 5.5, 5.2, 5.0, 4.8, 4.5, 4.0,
                    3.0, 2.0, 1.0, 0.5, 0.01,
                ],
                |s| list(s),
            );
            println!("tol,before_end,p_zeta,f_nl");
            for tol in list(&args[3]) {
                let mut transport = Transport::new(&system, trajectory, k, tol);
                for &before in &slices {
                    transport.advance(end - before);
                    let (f_nl, p) = transport.f_nl();
                    println!("{tol:e},{before},{p:.6e},{f_nl:.6e}");
                }
            }
        }
        "delta-n" => {
            let tol: f64 = args[3].parse().unwrap();
            let target_before: f64 = args[4].parse().unwrap();
            let slices = list(&args[5]);
            let scale = list(&args[6]);
            let scale = [scale[0], scale[1], scale[2]];
            let target_hubble = system
                .hubble(&trajectory.state_at(end - target_before))
                .into_value();
            let mut transport = Transport::new(&system, trajectory, k, tol);
            println!("n_c,f_nl,f_nl_gaussian,p_zeta");
            for before in slices {
                transport.advance(end - before);
                let (f_nl, gaussian, p) =
                    phase_space_delta_n(&transport, &settings, target_hubble, scale);
                println!("{before},{f_nl:.8e},{gaussian:.8e},{p:.8e}");
            }
            transport.advance(end - target_before);
            let (f_nl, p) = transport.f_nl();
            println!("transport,{f_nl:.8e},,{p:.8e}");
        }
        "response" => {
            // Integrate the response backwards from the surface `target` with
            // fixed RK4 steps `step`, and combine it with the transport on
            // each slice.
            let tol: f64 = args[3].parse().unwrap();
            let target_before: f64 = args[4].parse().unwrap();
            let slices = list(&args[5]);
            let step: f64 = args[6].parse().unwrap();
            let target = end - target_before;
            let triangle = MomentumTriangle::try_new(k, k, k).unwrap();
            let surface = trajectory.state_at(target);
            let gauge = curvature_gauge_coefficients(&system, &surface, triangle, 0.1).unwrap();
            let mut response = Response {
                linear: core::array::from_fn(|a| {
                    if a < 2 {
                        gauge.linear[a].into_value()
                    } else {
                        0.0
                    }
                }),
                quadratic: gauge.quadratic.map(|q| {
                    core::array::from_fn(|a| {
                        core::array::from_fn(|b| match (a < 2, b < 2) {
                            (true, true) => q.field_field[a][b].into_value(),
                            (true, false) => q.field_momentum[a][b - 2].into_value(),
                            (false, true) => q.momentum_field[a - 2][b].into_value(),
                            (false, false) => 0.0,
                        })
                    })
                }),
            };
            let rate = |n: f64, y: &Response| {
                let (u, kernels) = flows(&system, trajectory, triangle, n);
                y.derivative(&u, &kernels)
            };
            let mut transport = Transport::new(&system, trajectory, k, tol);
            let mut responses = Vec::new();
            let mut e_folds = target;
            let mut sorted = slices.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // Never step across a checkpoint of the background: the cubic
            // Hermite interpolation is only C^1 there, which would cost RK4
            // its order.
            let knots: Vec<f64> = survey.samples.iter().map(|sample| sample.e_folds).collect();
            for before in sorted {
                let stop = end - before;
                let mut nodes: Vec<f64> = knots
                    .iter()
                    .copied()
                    .filter(|&n| n < e_folds - 1.0e-12 && n > stop + 1.0e-12)
                    .collect();
                nodes.push(stop);
                nodes.sort_by(|a, b| b.partial_cmp(a).unwrap());
                for node in nodes {
                    let width = e_folds - node;
                    let count = (width / step).ceil().max(1.0);
                    let h = -width / count;
                    for _ in 0..count as usize {
                        let k1 = rate(e_folds, &response);
                        let k2 = rate(e_folds + 0.5 * h, &response.axpy(0.5 * h, &k1));
                        let k3 = rate(e_folds + 0.5 * h, &response.axpy(0.5 * h, &k2));
                        let k4 = rate(e_folds + h, &response.axpy(h, &k3));
                        response = response
                            .axpy(h / 6.0, &k1)
                            .axpy(h / 3.0, &k2)
                            .axpy(h / 3.0, &k3)
                            .axpy(h / 6.0, &k4);
                        e_folds += h;
                    }
                    e_folds = node;
                }
                responses.push((before, response));
            }
            println!("n_c,f_nl,p_zeta");
            for (before, response) in responses.iter().rev() {
                transport.advance(end - before);
                let (f_nl, p) = split_f_nl(response, &transport);
                println!("{before},{f_nl:.8e},{p:.8e}");
            }
            transport.advance(target);
            let (f_nl, p) = transport.f_nl();
            println!("transport,{f_nl:.8e},{p:.8e}");
        }
        "noise" => {
            // N along chi from the slice `n_c`, to the H surface and to epsilon_H = 1.
            let n_c: f64 = args[3].parse().unwrap();
            let target_before: f64 = args[4].parse().unwrap();
            let base = trajectory.state_at(end - n_c);
            let target_hubble = system
                .hubble(&trajectory.state_at(end - target_before))
                .into_value();
            let precise = SpectrumRunSettings {
                tolerances: ErrorTolerances::new(
                    ConformalPerturbationTolerance::<2>::background_in_planck_units(1.0e-24),
                    1.0e-13,
                ),
                background_sample_every: usize::MAX,
                ..settings
            };
            let h = 1.0e-4 * base.phi[1].into_value().abs();
            println!("t,n_h,n_end");
            for step in -8..=8 {
                let t = f64::from(step) / 8.0;
                let mut start = base;
                start.phi[1] = Field::new(base.phi[1].into_value() + t * h);
                start.a = ScaleFactor::new(1.0);
                let run = system.survey(start, &precise).unwrap();
                let n_h = run
                    .trajectory
                    .crossing(|st| system.hubble(st).into_value() - target_hubble)
                    .unwrap();
                println!("{t},{n_h:.15e},{:.15e}", run.e_folds_at_end);
            }
        }
        _ => panic!("unknown mode"),
    }
}
