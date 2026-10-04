//! The slow-roll delta N prescription of arXiv:2302.05866 (HKKPT).
//!
//! Eqs. (3.1), (3.5), and (3.9) of the reference predict `n_s` and `f_NL` of
//! the assisted model; the repository compares them with the exact
//! evolution in `docs/nonminimal-assisted-inflation.md`. These are the
//! formulas of one paper for one subclass, in the metric formulation.

use crate::assisted::{AssistedInflaton, bisect};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assisted::{HybridInflaton, LoopInflaton};

    #[test]
    fn the_prescription_reproduces_the_recorded_values() {
        // The values in golden/assisted_predictions.csv, computed before the
        // prescription moved out of the core.
        let (n_s, f_nl, _) = paper_prescription(&LoopInflaton { alpha: 0.1 }, 0.02282, 1.0e-3);
        assert!((n_s - 0.9541).abs() < 1.0e-4, "{n_s}");
        assert!((f_nl + 1.3647).abs() < 1.0e-3, "{f_nl}");
        let hybrid = HybridInflaton {
            mu: 10.0,
            phi_e: 1.0,
        };
        let (n_s, f_nl, _) = paper_prescription(&hybrid, 0.01077, 1.0e-3);
        assert!((n_s - 0.9621).abs() < 1.0e-4, "{n_s}");
        assert!((f_nl + 13.1717).abs() < 1.0e-3, "{f_nl}");
    }
}
