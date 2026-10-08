//! Closed forms of the exact constant-mixing solutions.
//!
//! Huenupi, Muñoz, Palma, and Sypsas solve the two-field system of constant
//! mixing `lambda` and entropy mass `mu` exactly (arXiv:2606.18248 and
//! 2607.14529). These are their power spectra, with the special functions
//! they need. See `docs/constant-turn-exact-solutions.md`.

use std::f64::consts::PI;

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

/// `ln |Gamma(re + i im)|` for `re > 0`, by recurrence to `Re z >= 20` and
/// the Stirling series.
pub fn ln_abs_gamma(re: f64, im: f64) -> f64 {
    let z = Complex(re, im);
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
pub fn re_digamma_imaginary(y: f64) -> f64 {
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
pub fn exact_massive_ratio(lambda: f64, mu_squared: f64) -> f64 {
    let g = ln_abs_gamma;
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
pub fn exact_massless_ratio(lambda: f64, z_end: f64, argument: f64) -> f64 {
    let euler = 0.577_215_664_901_532_9;
    let growth = 2.0 + (argument / z_end).ln() - 2.0 * euler - re_digamma_imaginary(0.5 * lambda);
    lambda * (PI * lambda).sinh() / (PI * (1.0 + lambda * lambda)) * growth * growth
        + PI * lambda / (2.0 * (1.0 + lambda * lambda)) / (0.5 * PI * lambda).tanh()
}
