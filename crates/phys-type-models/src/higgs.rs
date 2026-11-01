//! Single-field Higgs inflation as a generalized nonminimal theory.
//!
//! Bezrukov and Shaposhnikov,
//! [*The Standard Model Higgs boson as the inflaton*](https://arxiv.org/abs/0710.3755):
//! the Higgs field `phi` (in unitary gauge) with
//! `Omega^2 = 1 + xi phi^2 / M_pl^2`, a canonical Jordan-frame kinetic term,
//! and `V_J = lambda phi^4 / 4`. In the metric formulation the large-`xi`
//! predictions are `n_s = 1 - 2/N` and `r = 12/N^2`; in the Palatini
//! formulation ([Bauer and Demir, arXiv:0803.2664](https://arxiv.org/abs/0803.2664))
//! `n_s = 1 - 2/N` and `r = 2/(xi N^2)`. `phys_type::einstein_frame` turns the
//! theory into the Einstein-frame model that the engine evolves; the
//! slow-roll functions here use the same Einstein-frame potential and metric.

use phys_type::{
    Dimensionless, EnergyDensity, Field, Formulation, GeneralizedNonminimal, Quantity, TheoryScalar,
};

/// Higgs inflation with coupling `xi`, quartic `lambda`, and a formulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HiggsTheory {
    /// The nonminimal coupling `xi`.
    pub xi: f64,
    /// The quartic coupling `lambda`.
    pub lambda: f64,
    /// The metric or Palatini formulation.
    pub formulation: Formulation,
}

impl GeneralizedNonminimal<1> for HiggsTheory {
    fn formulation(&self) -> Formulation {
        self.formulation
    }

    fn conformal_factor<V: TheoryScalar>(&self, fields: [Field<V>; 1]) -> Dimensionless<V> {
        let phi = fields[0].into_value();
        Quantity::new(V::from(1.0) + phi * phi * self.xi)
    }

    fn kinetic_matrix<V: TheoryScalar>(
        &self,
        _fields: [Field<V>; 1],
    ) -> [[Dimensionless<V>; 1]; 1] {
        [[Quantity::new(V::from(1.0))]]
    }

    fn jordan_potential<V: TheoryScalar>(&self, fields: [Field<V>; 1]) -> EnergyDensity<V> {
        let phi = fields[0].into_value();
        EnergyDensity::new(phi * phi * phi * phi * (0.25 * self.lambda))
    }
}

impl HiggsTheory {
    /// `kappa` of the Einstein-frame metric.
    fn kappa(&self) -> f64 {
        self.formulation.kappa()
    }

    /// The Einstein-frame metric `G(phi) = (Omega^2 + 6 kappa xi^2 phi^2) / Omega^4`.
    pub fn metric(&self, phi: f64) -> f64 {
        let omega2 = 1.0 + self.xi * phi * phi;
        (omega2 + 6.0 * self.kappa() * self.xi * self.xi * phi * phi) / (omega2 * omega2)
    }

    /// `d ln V / d phi` of the Einstein-frame potential `V = lambda phi^4 / (4 Omega^4)`.
    fn log_slope(&self, phi: f64) -> f64 {
        4.0 / phi - 4.0 * self.xi * phi / (1.0 + self.xi * phi * phi)
    }

    /// First-order slow-roll `epsilon_V = (d ln V / d chi)^2 / 2` with
    /// `d chi = sqrt(G) d phi`.
    pub fn epsilon(&self, phi: f64) -> f64 {
        let slope = self.log_slope(phi);
        0.5 * slope * slope / self.metric(phi)
    }

    /// First-order slow-roll `eta_V = V_{,chi chi} / V`, from the analytic
    /// derivatives of `ln V` and `G` (the Einstein-frame potential is nearly
    /// flat at large `xi`, so finite differences would lose it to round-off).
    pub fn eta(&self, phi: f64) -> f64 {
        let xi = self.xi;
        let omega2 = 1.0 + xi * phi * phi;
        let slope = self.log_slope(phi);
        // d/dphi of d ln V / dphi.
        let curvature = -4.0 / (phi * phi) - 4.0 * xi * (1.0 - xi * phi * phi) / (omega2 * omega2);
        let numerator = omega2 + 6.0 * self.kappa() * xi * xi * phi * phi;
        let metric = numerator / (omega2 * omega2);
        let metric_slope = (2.0 * xi * phi + 12.0 * self.kappa() * xi * xi * phi)
            / (omega2 * omega2)
            - 4.0 * xi * phi * numerator / (omega2 * omega2 * omega2);
        // V_chichi / V = (V_phiphi/V - G'/(2G) V_phi/V) / G.
        (curvature + slope * slope - 0.5 * metric_slope / metric * slope) / metric
    }

    /// The field value where `epsilon_V = 1`.
    pub fn end_of_slow_roll(&self) -> f64 {
        let (mut low, mut high): (f64, f64) = (1.0e-8, 100.0);
        for _ in 0..200 {
            let middle = (low * high).sqrt();
            if self.epsilon(middle) > 1.0 {
                low = middle;
            } else {
                high = middle;
            }
        }
        (low * high).sqrt()
    }

    /// Slow-roll e-folds `N = int G / (d ln V / d phi) d phi` from `phi` to
    /// the end of slow roll.
    pub fn slow_roll_e_folds(&self, phi: f64) -> f64 {
        let end = self.end_of_slow_roll();
        let steps = 20_000;
        let width = (phi - end) / f64::from(steps);
        (0..steps)
            .map(|i| {
                let p = end + (f64::from(i) + 0.5) * width;
                self.metric(p) / self.log_slope(p) * width
            })
            .sum()
    }

    /// The field value `N` slow-roll e-folds before the end.
    pub fn field_at_e_folds(&self, e_folds: f64) -> f64 {
        let end = self.end_of_slow_roll();
        let (mut low, mut high) = (end, end * 1.0e4 + 100.0);
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if self.slow_roll_e_folds(middle) > e_folds {
                high = middle;
            } else {
                low = middle;
            }
        }
        0.5 * (low + high)
    }

    /// First-order slow-roll `(n_s, r)` at the field value `N` e-folds
    /// before the end of slow roll.
    pub fn slow_roll_predictions(&self, e_folds: f64) -> (f64, f64) {
        let phi = self.field_at_e_folds(e_folds);
        let (epsilon, eta) = (self.epsilon(phi), self.eta(phi));
        (1.0 - 6.0 * epsilon + 2.0 * eta, 16.0 * epsilon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slow_roll_reduces_to_quartic_inflation_without_coupling() {
        let quartic = HiggsTheory {
            xi: 0.0,
            lambda: 1.0,
            formulation: Formulation::Metric,
        };
        for phi in [3.0, 10.0, 22.0] {
            assert!((quartic.epsilon(phi) - 8.0 / (phi * phi)).abs() < 1.0e-12);
            assert!((quartic.eta(phi) - 12.0 / (phi * phi)).abs() < 1.0e-12);
        }
    }

    #[test]
    fn slow_roll_reaches_the_large_xi_limits() {
        // Metric: n_s = 1 - 2/N, r = 12/N^2; Palatini: r = 2/(xi N^2). Both
        // carry O(1/N) corrections (effectively a shift of N by a few), so at
        // N = 60 r lies about 10% below the asymptote.
        let e_folds = 60.0;
        let metric = HiggsTheory {
            xi: 1.0e4,
            lambda: 1.0,
            formulation: Formulation::Metric,
        };
        let (n_s, r) = metric.slow_roll_predictions(e_folds);
        assert!((n_s - (1.0 - 2.0 / e_folds)).abs() < 2.0e-3, "{n_s}");
        assert!((r / (12.0 / e_folds / e_folds) - 1.0).abs() < 0.15, "{r}");
        let palatini = HiggsTheory {
            formulation: Formulation::Palatini,
            ..metric
        };
        let (n_s, r) = palatini.slow_roll_predictions(e_folds);
        assert!((n_s - (1.0 - 2.0 / e_folds)).abs() < 2.0e-3, "{n_s}");
        assert!(
            (r / (2.0 / (1.0e4 * e_folds * e_folds)) - 1.0).abs() < 0.1,
            "{r}"
        );
    }
}
