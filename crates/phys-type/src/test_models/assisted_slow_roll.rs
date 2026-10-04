//! Slow-roll delta N predictions of nonminimally assisted inflation.
//!
//! Hyun, Kim, Kodama, Park, and Takahashi,
//! [*Nonminimally Assisted Inflation: A General Analysis*](https://arxiv.org/abs/2302.05866)
//! (arXiv:2302.05866v2, HKKPT), predict `n_s` and `f_NL` of the assisted
//! model of [`super::nonminimal`] from Eqs. (3.1), (3.5), and (3.9). The
//! repository compares them with the exact background, mode, and transport
//! evolution; see `docs/nonminimal-assisted-inflation.md`.
//!
//! These are the formulas of one reference for one subclass: a single
//! assistant field with `f(s) = (M_pl^2/2)(1 + xi_2 s^2)` and a minimally
//! coupled inflaton. They are not predictions for nonminimal inflation in
//! general, such as an inflaton coupled to curvature or a general `f(phi, s)`.
//! The loop, power-law, and hybrid inflatons are the examples of Secs. 4.1--4.3.
//! The metric formulation, `kappa = 1`, is assumed throughout.

use crate::{EndOfInflation, Field, Mass};

/// An inflaton of the assisted model with the slow-roll data of
/// HKKPT Sec. 4 that its predictions need.
pub trait AssistedInflaton: Copy + core::fmt::Debug {
    /// The Einstein-frame potential, with `Lambda^4 = 1e-10` and `kappa = 1`.
    type Potential: crate::ThriceDifferentiableMultifieldPotential<2> + Copy;
    /// The Einstein-frame potential for the coupling `xi_2`.
    fn potential(&self, coupling: f64) -> Self::Potential;
    /// `epsilon^(0)` and `eta^(0)` of `V_J` at `phi`, Eq. (3.14).
    fn slow_roll(&self, phi: f64) -> (f64, f64);
    /// E-folds of `phi` slow roll from `phi` to the reference's end value.
    fn slow_roll_e_folds(&self, phi: f64) -> f64;
    /// The end-of-inflation value of `phi` in the reference's prescription.
    fn paper_phi_e(&self) -> f64;
    /// A bracket for `phi` in which `slow_roll_e_folds` is monotonic.
    fn phi_bracket(&self) -> (f64, f64);
    /// The end-of-inflation event of the numerical background.
    fn end_of_inflation(&self) -> EndOfInflation;
    /// Class I ends by slow-roll violation, epsilon^phi_e = 1 - epsilon^sigma_e.
    fn class_one(&self) -> bool;
}

/// Bisection root of `f` on `[low, high]`, assuming one sign change.
pub fn bisect(f: &dyn Fn(f64) -> f64, mut low: f64, mut high: f64) -> f64 {
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        if f(low) * f(middle) <= 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    0.5 * (low + high)
}

/// Loop inflation, Sec. 4.1: Class I.
#[derive(Clone, Copy, Debug)]
pub struct LoopInflaton {
    /// The loop coefficient `alpha` of `V_J = Lambda^4 (1 + alpha ln phi)`.
    pub alpha: f64,
}

impl LoopInflaton {
    fn e_folds_from_zero(&self, phi: f64) -> f64 {
        let alpha = self.alpha;
        (phi * phi / 2.0 + alpha * (phi * phi * phi.ln() / 2.0 - phi * phi / 4.0)) / alpha
    }
}

impl AssistedInflaton for LoopInflaton {
    type Potential = crate::test_models::NonminimalLoopInflationPotential;

    fn potential(&self, coupling: f64) -> Self::Potential {
        crate::test_models::NonminimalLoopInflationPotential {
            scale: crate::EnergyDensity::new(1.0e-10),
            alpha: self.alpha,
            coupling,
            planck_mass: Mass::new(1.0),
        }
    }

    fn slow_roll(&self, phi: f64) -> (f64, f64) {
        let factor = 1.0 + self.alpha * phi.ln();
        // Corrected Eq. (4.2): phi^2, not phi, in epsilon^(0).
        (
            self.alpha * self.alpha / (2.0 * phi * phi * factor * factor),
            -self.alpha / (phi * phi * factor),
        )
    }

    fn slow_roll_e_folds(&self, phi: f64) -> f64 {
        self.e_folds_from_zero(phi) - self.e_folds_from_zero(self.paper_phi_e())
    }

    fn paper_phi_e(&self) -> f64 {
        // phi_e = [sqrt(2) W_0(e^{1/alpha} / sqrt(2))]^{-1}; W_0 by Newton.
        let argument = (1.0 / self.alpha).exp() / 2.0_f64.sqrt();
        let mut w = argument.ln().max(0.5);
        for _ in 0..50 {
            w -= (w * w.exp() - argument) / (w.exp() * (w + 1.0));
        }
        1.0 / (2.0_f64.sqrt() * w)
    }

    fn phi_bracket(&self) -> (f64, f64) {
        (self.paper_phi_e() * 1.01, 50.0)
    }

    fn end_of_inflation(&self) -> EndOfInflation {
        EndOfInflation::SlowRollViolation
    }

    fn class_one(&self) -> bool {
        true
    }
}

/// Power-law inflation, Sec. 4.2: Class II with constant slow-roll
/// parameters. Inflation starts at phi = 0 and a separate sector ends it
/// after about 72 e-folds of phi slow roll.
#[derive(Clone, Copy, Debug)]
pub struct PowerLawInflaton {
    /// The slope `lambda` of `V_J = Lambda^4 e^{-lambda phi}`.
    pub lambda: f64,
}

impl AssistedInflaton for PowerLawInflaton {
    type Potential = crate::test_models::NonminimalPowerLawPotential;

    fn potential(&self, coupling: f64) -> Self::Potential {
        crate::test_models::NonminimalPowerLawPotential {
            scale: crate::EnergyDensity::new(1.0e-10),
            lambda: self.lambda,
            coupling,
            planck_mass: Mass::new(1.0),
        }
    }

    fn slow_roll(&self, _phi: f64) -> (f64, f64) {
        (self.lambda * self.lambda / 2.0, self.lambda * self.lambda)
    }

    fn slow_roll_e_folds(&self, phi: f64) -> f64 {
        (self.paper_phi_e() - phi) / self.lambda
    }

    fn paper_phi_e(&self) -> f64 {
        72.0 * self.lambda
    }

    fn phi_bracket(&self) -> (f64, f64) {
        (-100.0, self.paper_phi_e())
    }

    fn end_of_inflation(&self) -> EndOfInflation {
        EndOfInflation::FieldValue {
            field: 0,
            value: Field::new(self.paper_phi_e()),
        }
    }

    fn class_one(&self) -> bool {
        false
    }
}

/// Hybrid inflation, Sec. 4.3: Class II, ended by the waterfall at
/// phi = phi_e.
#[derive(Clone, Copy, Debug)]
pub struct HybridInflaton {
    /// The scale `mu` of `V_J = Lambda^4 [1 + (phi/mu)^2]`.
    pub mu: f64,
    /// The waterfall value `phi_e` that ends inflation.
    pub phi_e: f64,
}

impl AssistedInflaton for HybridInflaton {
    type Potential = crate::test_models::NonminimalHybridPotential;

    fn potential(&self, coupling: f64) -> Self::Potential {
        crate::test_models::NonminimalHybridPotential {
            scale: crate::EnergyDensity::new(1.0e-10),
            mu: self.mu,
            coupling,
            planck_mass: Mass::new(1.0),
        }
    }

    fn slow_roll(&self, phi: f64) -> (f64, f64) {
        let ratio = phi / self.mu;
        let factor = 1.0 + ratio * ratio;
        (
            2.0 * ratio * ratio / (self.mu * self.mu * factor * factor),
            2.0 / (self.mu * self.mu * factor),
        )
    }

    fn slow_roll_e_folds(&self, phi: f64) -> f64 {
        // N = int V_J / V_J' dphi = (mu^2 / 2) ln(phi / phi_e) + (phi^2 - phi_e^2) / 4.
        self.mu * self.mu / 2.0 * (phi / self.phi_e).ln()
            + (phi * phi - self.phi_e * self.phi_e) / 4.0
    }

    fn paper_phi_e(&self) -> f64 {
        self.phi_e
    }

    fn phi_bracket(&self) -> (f64, f64) {
        (self.phi_e * 1.0001, 100.0 * self.mu)
    }

    fn end_of_inflation(&self) -> EndOfInflation {
        EndOfInflation::FieldValue {
            field: 0,
            value: Field::new(self.phi_e),
        }
    }

    fn class_one(&self) -> bool {
        false
    }
}

/// Slow-roll quantities of Eqs. (3.4) and (3.8) of arXiv:2302.05866.
pub struct AssistedSlowRoll {
    /// `epsilon^sigma`, Eq. (3.4).
    pub epsilon_sigma: f64,
    /// `eta^{sigma sigma}`, Eq. (3.4).
    pub eta_sigma_sigma: f64,
    /// `epsilon^b`, Eq. (3.8).
    pub epsilon_b: f64,
    /// `eta^b`, Eq. (3.8).
    pub eta_b: f64,
    /// `epsilon^phi = e^{2b} epsilon^(0)`.
    pub epsilon_phi: f64,
    /// `eta^{phi phi} = e^{2b} eta^(0)`.
    pub eta_phi_phi: f64,
    /// `2b = -ln(1 + xi_2 s^2)`.
    pub two_b: f64,
}

/// Evaluates the slow-roll parameters at `(phi, s)` for coupling `xi`.
pub fn assisted_slow_roll<I: AssistedInflaton>(
    inflaton: &I,
    xi: f64,
    phi: f64,
    s: f64,
) -> AssistedSlowRoll {
    let kappa = 1.0;
    let denominator = 1.0 + xi * s * s * (1.0 + 6.0 * kappa * xi);
    let epsilon_sigma = 8.0 * xi * xi * s * s / denominator;
    let eta_sigma_sigma =
        4.0 * xi * (-1.0 + 3.0 * xi * s * s + 4.0 * xi * xi * s.powi(4) * (1.0 + 6.0 * kappa * xi))
            / (denominator * denominator);
    let eta_b = -16.0 * xi * (1.0 + xi * s * s) / (denominator * denominator);
    let (epsilon_zero, eta_zero) = inflaton.slow_roll(phi);
    let two_b = -(1.0 + xi * s * s).ln();
    AssistedSlowRoll {
        epsilon_sigma,
        eta_sigma_sigma,
        epsilon_b: epsilon_sigma,
        eta_b,
        epsilon_phi: epsilon_zero * (-two_b).exp(),
        eta_phi_phi: eta_zero * (-two_b).exp(),
        two_b,
    }
}

/// Eqs. (3.1) and (3.5) of arXiv:2302.05866. With `class_one`, the Class I
/// end condition epsilon^phi_e = 1 - epsilon^sigma_e replaces the end
/// value of epsilon^phi; Class II uses it directly. Returns `(n_s, f_NL)`
/// with `f_NL` converted to the Maldacena/Planck sign used here: the
/// printed Eq. (3.5) gives `-(6/5) f_NL`, and its single-field limit is
/// `-(5/6)(2 epsilon - eta)`, the opposite of the Planck-sign delta N
/// result `(5/6)(2 epsilon - eta)`.
pub fn assisted_delta_n_prediction(
    star: &AssistedSlowRoll,
    end: &AssistedSlowRoll,
    class_one: bool,
) -> (f64, f64) {
    let signs = 1.0; // s^b s^sigma: both derivatives are negative for s > 0
    let epsilon_sigma_e = end.epsilon_sigma;
    let epsilon_phi_e = if class_one {
        1.0 - epsilon_sigma_e
    } else {
        end.epsilon_phi
    };
    let epsilon_e = epsilon_sigma_e + epsilon_phi_e;
    let u = epsilon_sigma_e / epsilon_e;
    let v = epsilon_phi_e / epsilon_e;
    let x = end.two_b - star.two_b;
    let a = (-x).exp() * (1.0 + epsilon_phi_e / epsilon_sigma_e * (1.0 - x.exp()));
    let (es, ep) = (star.epsilon_sigma, star.epsilon_phi);
    let d = u * u * a * a / es + v * v / ep;
    let n_s = 1.0
        - 2.0 * es
        - 2.0 * ep
        - 4.0 * (-2.0 * x).exp() / d
        - (star.eta_b + 2.0 * star.epsilon_b) / (12.0 * d)
            * (u * a * (ep / es).sqrt() - v * (es / ep).sqrt()).powi(2)
        + 2.0 / d
            * (u * u * a * a * star.eta_sigma_sigma / es
                + v * v * star.eta_phi_phi / ep
                + 4.0 * u * v * a
                + 0.5 * signs * (star.epsilon_b * es).sqrt() * v * (v / ep - 2.0 * u * a / es));
    let c = epsilon_sigma_e * epsilon_phi_e / (epsilon_e * epsilon_e)
        * ((epsilon_sigma_e * end.eta_phi_phi + epsilon_phi_e * end.eta_sigma_sigma) / epsilon_e
            - 4.0 * epsilon_phi_e * epsilon_sigma_e / epsilon_e
            - 0.5 * signs * (end.epsilon_b / epsilon_sigma_e).sqrt() * epsilon_phi_e.powi(2)
                / epsilon_e);
    let bracket = (1.0 - star.eta_sigma_sigma / (2.0 * es)) * (u * a).powi(3) / es
        + (1.0 - star.eta_phi_phi / (2.0 * ep)) * v.powi(3) / ep
        + 0.5 * signs * u * u * v * a * a / es * (star.epsilon_b / es).sqrt()
        + (u * a / es - v / ep).powi(2) * x.exp() * c;
    let minus_six_fifths_f_nl = 2.0 * (-x).exp() / (d * d) * bracket;
    (n_s, 5.0 / 6.0 * minus_six_fifths_f_nl)
}

/// Eqs. (3.1) and (3.5) evaluated entirely with the prescription of
/// HKKPT: `phi_*` from 60 e-folds of `phi` slow roll, and `s_e` from
/// Eq. (3.9) with `N_e = 60`. Returns `(n_s, f_NL, s_e)`, with `f_NL` in the
/// Maldacena/Planck sign.
pub fn paper_prescription<I: AssistedInflaton>(
    inflaton: &I,
    xi: f64,
    s_star: f64,
) -> (f64, f64, f64) {
    let (low, high) = inflaton.phi_bracket();
    let phi_star = bisect(&|phi| inflaton.slow_roll_e_folds(phi) - 60.0, low, high);
    let e_folds = |s_e: f64| {
        ((s_e / s_star).ln()
            + 3.0 * xi * ((1.0 + xi * s_e * s_e) / (1.0 + xi * s_star * s_star)).ln())
            / (4.0 * xi)
            - 60.0
    };
    let s_e = bisect(&e_folds, s_star, 100.0);
    let (n_s, f_nl) = assisted_delta_n_prediction(
        &assisted_slow_roll(inflaton, xi, phi_star, s_star),
        &assisted_slow_roll(inflaton, xi, inflaton.paper_phi_e(), s_e),
        inflaton.class_one(),
    );
    (n_s, f_nl, s_e)
}
