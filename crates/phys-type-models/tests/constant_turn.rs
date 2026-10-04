//! Validation of the mode solver and the three-point transport against the
//! exact constant-turn solutions of arXiv:2606.18248 and the squeezed
//! scaling of arXiv:2607.14529. See `docs/constant-turn-exact-solutions.md`.

use std::f64::consts::PI;

use phys_type::{
    CovariantMultifieldSystem, EndOfInflation, Field, FieldVelocity, Mass, MultifieldState,
    ScaleFactor, SpectrumRunSettings,
};

/// `V = V_0 - c theta + m^2 (rho - rho_0)^2 / 2` on [`phys_type::PolarFieldSpace`]:
/// a circular orbit of nearly constant turn rate and entropy mass, the
/// background assumed by the exact solutions of arXiv:2606.18248.
#[derive(Clone, Copy, Debug)]
struct ConstantTurnPotential {
    v0: f64,
    slope: f64,
    mass_squared: f64,
    rho0: f64,
}

impl<V> phys_type::MultifieldPotentialValue<V, 2> for ConstantTurnPotential
where
    V: From<f64>
        + Clone
        + core::ops::Add<V, Output = V>
        + core::ops::Mul<V, Output = V>
        + core::ops::Mul<f64, Output = V>,
{
    fn value(&self, fields: [Field<V>; 2]) -> phys_type::EnergyDensity<V> {
        let [rho, theta] = fields.map(phys_type::Quantity::into_value);
        let radial = rho + V::from(-self.rho0);
        phys_type::EnergyDensity::new(
            V::from(self.v0)
                + theta * (-self.slope)
                + radial.clone() * radial * (0.5 * self.mass_squared),
        )
    }
}

/// The mode-solver `P_zeta / Delta_0` on a constant-turn orbit, with the
/// mixing `lambda = 2 |Omega| / H` and entropy mass
/// `mu^2 / H^2 = V_NN / H^2 + 3 (Omega / H)^2` measured at horizon exit.
struct ConstantTurnRun {
    ratio: f64,
    lambda: f64,
    mu_squared: f64,
    k_over_ah: f64,
}

/// Builds the orbit for target `lambda` and `mu^2 / H^2` with first
/// slow-roll parameter `epsilon = rho^2 lambda^2 / 8`, and evaluates the
/// pivot `after_exit` e-folds after horizon exit. The orbit parameters
/// are corrected until the values measured at exit match the targets.
fn constant_turn_run(
    lambda: f64,
    mu_squared: f64,
    epsilon: f64,
    after_exit: f64,
) -> ConstantTurnRun {
    let radius = (8.0 * epsilon).sqrt() / lambda;
    let (mut lambda_in, mut mu_squared_in) = (lambda, mu_squared);
    let mut run = constant_turn_orbit(lambda_in, mu_squared_in, radius, after_exit);
    for _ in 0..3 {
        lambda_in *= lambda / run.lambda;
        mu_squared_in += mu_squared - run.mu_squared;
        run = constant_turn_orbit(lambda_in, mu_squared_in, radius, after_exit);
    }
    run
}

fn constant_turn_orbit(
    lambda: f64,
    mu_squared: f64,
    radius: f64,
    after_exit: f64,
) -> ConstantTurnRun {
    let hubble = 1.0e-5;
    let omega = 0.5 * lambda * hubble;
    let mass_squared = (mu_squared - 0.75 * lambda * lambda) * hubble * hubble;
    // Radial balance m^2 (rho - rho_0) = rho Omega^2 and the angular
    // slow-roll attractor c = 3 H rho^2 Omega.
    let potential = ConstantTurnPotential {
        v0: 3.0 * hubble * hubble,
        slope: 3.0 * hubble * radius * radius * omega,
        mass_squared,
        rho0: radius - radius * omega * omega / mass_squared,
    };
    let system = CovariantMultifieldSystem::in_planck_units(potential, phys_type::PolarFieldSpace);
    let before_exit = 6.0;
    let pivot_before_end = after_exit + 1.0;
    let total = pivot_before_end + before_exit + 8.0;
    let settings = SpectrumRunSettings {
        end_of_inflation: EndOfInflation::FieldValue {
            field: 1,
            value: Field::new(omega / hubble * total),
        },
        e_folds_before_exit: before_exit,
        pivot_e_folds_before_end: pivot_before_end,
        output_e_folds_before_end: 1.0,
        k_min_over_pivot: 1.0,
        k_max_over_pivot: 1.0,
        k_count: 1,
        ..SpectrumRunSettings::default()
    };
    let initial = MultifieldState::new(
        [Field::new(radius), Field::new(0.0)],
        [FieldVelocity::new(0.0), FieldVelocity::new(omega)],
        ScaleFactor::new(1.0),
    );
    let survey = system.survey(initial, &settings).unwrap();
    let (trajectory, end) = (survey.trajectory, survey.e_folds_at_end);
    let exit = trajectory.state_at(end - pivot_before_end);
    let h = system.hubble(&exit).into_value();
    let k = (exit.a * system.hubble(&exit)).into_value();
    let record = system
        .mode_spectrum(
            initial,
            &settings,
            &trajectory,
            Mass::new(k),
            end - pivot_before_end + after_exit,
        )
        .unwrap();
    let eta_perp = system
        .kinematics(&exit, None)
        .unwrap()
        .eta_perp
        .into_value();
    let delta_0 = h * h / (8.0 * PI * PI * system.epsilon_h(&exit));
    ConstantTurnRun {
        ratio: record.p_r / delta_0,
        lambda: 2.0 * eta_perp.abs(),
        mu_squared: mass_squared / (h * h) + 3.0 * eta_perp * eta_perp,
        k_over_ah: record.k_over_ah,
    }
}

/// Minimal complex arithmetic for the closed forms below.
#[derive(Clone, Copy, Debug)]
struct Complex(f64, f64);

impl Complex {
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0, self.1 + other.1)
    }

    fn mul(self, other: Self) -> Self {
        Self(
            self.0 * other.0 - self.1 * other.1,
            self.0 * other.1 + self.1 * other.0,
        )
    }

    fn inv(self) -> Self {
        let norm = self.0 * self.0 + self.1 * self.1;
        Self(self.0 / norm, -self.1 / norm)
    }

    fn ln(self) -> Self {
        Self(self.0.hypot(self.1).ln(), self.1.atan2(self.0))
    }
}

/// `ln |Gamma(z)|` for `Re z > 0`, by recurrence to `Re z >= 20` and
/// the Stirling series.
fn ln_abs_gamma(z: Complex) -> f64 {
    let mut shifted = z;
    let mut correction = 0.0;
    while shifted.0 < 20.0 {
        correction += shifted.ln().0;
        shifted = shifted.add(Complex(1.0, 0.0));
    }
    let w = shifted;
    let inv = w.inv();
    let inv2 = inv.mul(inv);
    let series = inv.mul(
        Complex(1.0 / 12.0, 0.0)
            .add(inv2.mul(Complex(-1.0 / 360.0, 0.0)))
            .add(inv2.mul(inv2).mul(Complex(1.0 / 1260.0, 0.0))),
    );
    let stirling = w
        .add(Complex(-0.5, 0.0))
        .mul(w.ln())
        .add(Complex(-w.0, -w.1))
        .add(Complex(0.5 * (2.0 * PI).ln(), 0.0))
        .add(series);
    stirling.0 - correction
}

/// `Re psi(i y)`, which equals `Re psi(1 + i y)`.
fn re_digamma_imaginary(y: f64) -> f64 {
    let mut shifted = Complex(1.0, y);
    let mut correction = 0.0;
    while shifted.0 < 20.0 {
        correction += shifted.inv().0;
        shifted = shifted.add(Complex(1.0, 0.0));
    }
    let inv = shifted.inv();
    let inv2 = inv.mul(inv);
    let asymptotic = shifted
        .ln()
        .add(inv.mul(Complex(-0.5, 0.0)))
        .add(inv2.mul(Complex(-1.0 / 12.0, 0.0)))
        .add(inv2.mul(inv2).mul(Complex(1.0 / 120.0, 0.0)))
        .add(inv2.mul(inv2).mul(inv2).mul(Complex(-1.0 / 252.0, 0.0)));
    asymptotic.0 - correction
}

/// `Delta_zeta / Delta_0` of arXiv:2606.18248 for a massive entropy
/// field: Eq. (43) for `mu^2 <= 9/4` and Eq. (45) above it.
fn exact_massive_ratio(lambda: f64, mu_squared: f64) -> f64 {
    let g = |re: f64, im: f64| ln_abs_gamma(Complex(re, im));
    let log_ratio = if mu_squared <= 2.25 {
        let nu = (2.25 - mu_squared).sqrt();
        let (a, b) = (0.75 - 0.5 * nu, 0.75 + 0.5 * nu);
        2.0 * (g(a, 0.0) + g(b, 0.0) - g(a, 0.5 * lambda) - g(b, 0.5 * lambda))
    } else {
        let rho = (mu_squared - 2.25).sqrt();
        4.0 * g(0.75, 0.5 * rho)
            - 2.0 * g(0.75, 0.5 * (lambda + rho))
            - 2.0 * g(0.75, 0.5 * (lambda - rho))
    };
    log_ratio.exp()
}

/// Eq. (39) of arXiv:2606.18248 for a massless entropy field, with the
/// logarithm `ln(argument / z_end)`. The printed `argument = 2`
/// disagrees with the equations of motion; `argument = 1/2` agrees.
fn exact_massless_ratio(lambda: f64, z_end: f64, argument: f64) -> f64 {
    let euler = 0.577_215_664_901_532_9;
    let growth = 2.0 + (argument / z_end).ln() - 2.0 * euler - re_digamma_imaginary(0.5 * lambda);
    lambda * (PI * lambda).sinh() / (PI * (1.0 + lambda * lambda)) * growth * growth
        + PI * lambda / (2.0 * (1.0 + lambda * lambda)) / (0.5 * PI * lambda).tanh()
}

#[test]
fn closed_form_gamma_and_digamma_match_known_values() {
    // |Gamma(1/2 + i y)|^2 = pi / cosh(pi y) and Gamma(5) = 24.
    let y = 0.7_f64;
    let reflection = (PI / (PI * y).cosh()).ln();
    assert!((2.0 * ln_abs_gamma(Complex(0.5, y)) - reflection).abs() < 1.0e-12);
    assert!((ln_abs_gamma(Complex(5.0, 0.0)) - 24.0_f64.ln()).abs() < 1.0e-12);
    // Re psi(i y) = -gamma_E + y^2 sum 1 / (n (n^2 + y^2)).
    let series = -0.577_215_664_901_532_9
        + (1..200_000)
            .map(|n| {
                let n = f64::from(n);
                y * y / (n * (n * n + y * y))
            })
            .sum::<f64>();
    assert!((re_digamma_imaginary(y) - series).abs() < 1.0e-9);
    // Eq. (43) at nu = 1/2 is Eq. (44), and Eqs. (43) and (45) meet at
    // mu^2 = 9/4.
    let conformal = (PI * 1.3).sinh() / (PI * 1.3);
    assert!((exact_massive_ratio(1.3, 2.0) / conformal - 1.0).abs() < 1.0e-12);
    assert!(
        (exact_massive_ratio(1.3, 2.25) / exact_massive_ratio(1.3, 2.25 + 1.0e-10) - 1.0).abs()
            < 1.0e-6
    );
}

#[test]
fn constant_turn_spectrum_matches_the_exact_massive_solutions() {
    // (lambda, mu^2 / H^2, e-folds after exit): light, conformal, heavy.
    let cases = [
        (1.0, 1.25, 20.0),
        (2.0, 1.25, 20.0),
        (0.5, 2.0, 12.0),
        (1.0, 2.0, 12.0),
        (2.0, 2.0, 12.0),
        (1.0, 4.0, 12.0),
        (2.0, 4.0, 12.0),
        (1.0, 9.0, 12.0),
    ];
    for (lambda, mu_squared, after_exit) in cases {
        // Compared at the measured parameters, so no orbit correction.
        let radius = (8.0e-5_f64).sqrt() / lambda;
        let run = constant_turn_orbit(lambda, mu_squared, radius, after_exit);
        let exact = exact_massive_ratio(run.lambda, run.mu_squared);
        assert!(
            (run.ratio / exact - 1.0).abs() < 1.0e-3,
            "lambda {lambda}, mu^2 {mu_squared}: {} vs {exact}",
            run.ratio
        );
    }
}

#[test]
fn constant_turn_spectrum_fixes_the_logarithm_of_the_massless_solution() {
    // The slow-roll correction accumulates on superhorizon scales, so
    // strong mixing needs a smaller epsilon.
    for (lambda, epsilon) in [(0.5, 3.0e-6), (1.0, 3.0e-6), (2.0, 3.0e-7)] {
        let run = constant_turn_run(lambda, 0.0, epsilon, 10.0);
        let corrected = exact_massless_ratio(run.lambda, run.k_over_ah, 0.5);
        let printed = exact_massless_ratio(run.lambda, run.k_over_ah, 2.0);
        assert!(
            (run.ratio / corrected - 1.0).abs() < 1.0e-3,
            "lambda {lambda}: {} vs {corrected}",
            run.ratio
        );
        assert!(printed / run.ratio - 1.0 > 0.2, "lambda {lambda}");
    }
}

/// Squeezed `f_NL` on the constant-turn orbit for isosceles triangles
/// `k_1 = k_2 = k_pivot`, `k_3 = r k_pivot`.
fn constant_turn_squeezed_f_nl(
    lambda: f64,
    mu_squared: f64,
    after_exit: f64,
    ratios: &[f64],
) -> Vec<f64> {
    let hubble = 1.0e-5;
    let radius = (8.0e-5_f64).sqrt() / lambda;
    let omega = 0.5 * lambda * hubble;
    let mass_squared = (mu_squared - 0.75 * lambda * lambda) * hubble * hubble;
    let potential = ConstantTurnPotential {
        v0: 3.0 * hubble * hubble,
        slope: 3.0 * hubble * radius * radius * omega,
        mass_squared,
        rho0: radius - radius * omega * omega / mass_squared,
    };
    let system = CovariantMultifieldSystem::in_planck_units(potential, phys_type::PolarFieldSpace);
    let pivot_before_end = after_exit + 1.0;
    let total = pivot_before_end + 16.0;
    let settings = phys_type::BispectrumRunSettings {
        background: SpectrumRunSettings {
            end_of_inflation: EndOfInflation::FieldValue {
                field: 1,
                value: Field::new(omega / hubble * total),
            },
            pivot_e_folds_before_end: pivot_before_end,
            output_e_folds_before_end: 1.0,
            k_min_over_pivot: 1.0,
            k_max_over_pivot: 1.0,
            k_count: 1,
            ..SpectrumRunSettings::default()
        },
        ..phys_type::BispectrumRunSettings::default()
    };
    let initial = MultifieldState::new(
        [Field::new(radius), Field::new(0.0)],
        [FieldVelocity::new(0.0), FieldVelocity::new(omega)],
        ScaleFactor::new(1.0),
    );
    let shapes: Vec<phys_type::TriangleShape> = ratios
        .iter()
        .map(|&r| phys_type::TriangleShape {
            k_t_over_pivot: 2.0 + r,
            alpha: 0.0,
            beta: (2.0 - r) / (2.0 + r),
        })
        .collect();
    system
        .run_bispectrum(initial, &settings, &shapes)
        .unwrap()
        .records
        .iter()
        .map(|record| record.f_nl)
        .collect()
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn constant_turn_squeezed_bispectrum_scales_with_mu_not_the_effective_mass() {
    // arXiv:2607.14529: the leading squeezed shape is r^(1/2 - nu), so the
    // reduced f_NL scales as r^(3/2 - nu) with nu = sqrt(9/4 - mu^2/H^2),
    // not with the effective mass m_eff^2 = mu^2 - lambda^2 H^2.
    for (lambda, mu_squared, after_exit) in
        [(1.0, 1.25, 20.0), (2.0, 1.25, 20.0), (2.0, 0.81, 30.0)]
    {
        let f_nl = constant_turn_squeezed_f_nl(lambda, mu_squared, after_exit, &[0.03, 0.01]);
        let slope = (f_nl[0] / f_nl[1]).ln() / 3.0_f64.ln();
        let exponent = |m2: f64| 1.5 - (2.25 - m2).sqrt();
        let predicted = exponent(mu_squared);
        let effective = exponent(mu_squared - lambda * lambda);
        assert!(
            (slope - predicted).abs() < 0.03,
            "lambda {lambda}, mu^2 {mu_squared}: slope {slope} vs {predicted}"
        );
        assert!((slope - effective).abs() > 0.3, "{slope} vs {effective}");
    }
}
