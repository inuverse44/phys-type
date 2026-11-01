//! Nonminimally assisted inflation as a generalized nonminimal theory.
//!
//! Hyun, Kim, Kodama, Park, and Takahashi,
//! [*Nonminimally Assisted Inflation: A General Analysis*](https://arxiv.org/abs/2302.05866)
//! (arXiv:2302.05866v2, HKKPT): an inflaton `phi` with a Jordan-frame
//! potential `V_J(phi)` and an assistant field `s` coupled to curvature
//! through `f(s) = (M_pl^2/2)(1 + xi_2 s^2)`. In the theory layer of
//! `phys-type` this is `Omega^2 = 1 + xi_2 s^2` with canonical Jordan-frame
//! kinetic terms; `phys_type::einstein_frame` turns it into the Einstein-frame
//! model that the engine evolves.
//!
//! This is one small class of nonminimal inflation: one assistant field,
//! a minimally coupled inflaton, and the examples of the reference's
//! Secs. 4.1--4.3 (loop, power-law, and hybrid inflatons).

use phys_type::{
    Dimensionless, EndOfInflation, EnergyDensity, Field, Formulation, GeneralizedNonminimal,
    Quantity, TheoryScalar,
};

/// `Lambda^4` of every example potential, in reduced Planck units.
pub const SCALE: f64 = 1.0e-10;

/// An inflaton of the assisted model with the slow-roll data of
/// HKKPT Sec. 4 that its predictions need.
pub trait AssistedInflaton: Copy + core::fmt::Debug {
    /// The Jordan-frame potential `V_J(phi)` in reduced Planck units, with
    /// `Lambda^4 = 1e-10`.
    fn jordan_potential<V: TheoryScalar>(&self, phi: V) -> V;
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
    fn jordan_potential<V: TheoryScalar>(&self, phi: V) -> V {
        (V::from(1.0) + phi.ln() * self.alpha) * SCALE
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
    fn jordan_potential<V: TheoryScalar>(&self, phi: V) -> V {
        (phi * -self.lambda).exp() * SCALE
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
    fn jordan_potential<V: TheoryScalar>(&self, phi: V) -> V {
        (V::from(1.0) + phi * phi * (1.0 / (self.mu * self.mu))) * SCALE
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

/// The assisted theory: an inflaton, the coupling `xi_2`, and the
/// formulation of gravity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AssistedTheory<I> {
    /// The inflaton and its Jordan-frame potential.
    pub inflaton: I,
    /// The nonminimal coupling `xi_2` of the assistant field.
    pub coupling: f64,
    /// The metric or Palatini formulation.
    pub formulation: Formulation,
}

impl<I: AssistedInflaton> AssistedTheory<I> {
    /// The metric formulation, used throughout the reference.
    pub const fn metric(inflaton: I, coupling: f64) -> Self {
        Self {
            inflaton,
            coupling,
            formulation: Formulation::Metric,
        }
    }
}

impl<I: AssistedInflaton> GeneralizedNonminimal<2> for AssistedTheory<I> {
    fn formulation(&self) -> Formulation {
        self.formulation
    }

    fn conformal_factor<V: TheoryScalar>(&self, fields: [Field<V>; 2]) -> Dimensionless<V> {
        let [_, s] = fields.map(Quantity::into_value);
        Quantity::new(V::from(1.0) + s * s * self.coupling)
    }

    fn kinetic_matrix<V: TheoryScalar>(
        &self,
        _fields: [Field<V>; 2],
    ) -> [[Dimensionless<V>; 2]; 2] {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| Quantity::new(V::from(f64::from(i == j))))
        })
    }

    fn jordan_potential<V: TheoryScalar>(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        EnergyDensity::new(self.inflaton.jordan_potential(fields[0].into_value()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use phys_type::{MultifieldPotentialValue, einstein_frame};

    #[test]
    fn the_einstein_potential_is_the_jordan_potential_over_omega_to_the_fourth() {
        let cases: [(f64, f64); 3] = [(5.0, 1.0e-3), (2.0, 0.4), (0.8, 2.3)];
        let xi = 0.022;
        let check = |jordan: &dyn Fn(f64) -> f64, einstein: &dyn Fn(f64, f64) -> f64| {
            for (phi, s) in cases {
                let expected = jordan(phi) / (1.0 + xi * s * s).powi(2);
                let actual = einstein(phi, s);
                assert!((actual - expected).abs() <= 1.0e-14 * expected.abs());
            }
        };
        let loop_inflaton = LoopInflaton { alpha: 0.1 };
        let (potential, _) = einstein_frame(AssistedTheory::metric(loop_inflaton, xi));
        check(&|phi| SCALE * (1.0 + 0.1 * phi.ln()), &|phi, s| {
            potential
                .value([Field::new(phi), Field::new(s)])
                .into_value()
        });
        let power_law = PowerLawInflaton { lambda: 0.05 };
        let (potential, _) = einstein_frame(AssistedTheory::metric(power_law, xi));
        check(&|phi| SCALE * (-0.05 * phi).exp(), &|phi, s| {
            potential
                .value([Field::new(phi), Field::new(s)])
                .into_value()
        });
        let hybrid = HybridInflaton {
            mu: 10.0,
            phi_e: 1.0,
        };
        let (potential, _) = einstein_frame(AssistedTheory::metric(hybrid, xi));
        check(&|phi| SCALE * (1.0 + phi * phi / 100.0), &|phi, s| {
            potential
                .value([Field::new(phi), Field::new(s)])
                .into_value()
        });
    }
}
