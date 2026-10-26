//! Two nonminimally coupled fields with quartic and mass terms.
//!
//! Kaiser, Mazenc, and Sfakianakis,
//! [*Primordial bispectrum from multifield inflation with nonminimal
//! couplings*](https://arxiv.org/abs/1210.7487) (arXiv:1210.7487, KMS): two
//! fields `(phi, chi)` with canonical Jordan-frame kinetic terms,
//! `f = (M_pl^2 + xi_phi phi^2 + xi_chi chi^2) / 2`, Eq. (17), and
//!
//! ```text
//! V_J = m_phi^2 phi^2 / 2 + m_chi^2 chi^2 / 2 + g phi^2 chi^2 / 2
//!     + lambda_phi phi^4 / 4 + lambda_chi chi^4 / 4,
//! ```
//!
//! Eq. (16). In the theory layer this is
//! `Omega^2 = 1 + xi_phi phi^2 + xi_chi chi^2` in reduced Planck units. Both
//! fields enter `Omega^2`, so the Einstein-frame metric is not diagonal and
//! its curvature varies in both directions, unlike the assisted model.
//!
//! The closed forms below are test oracles for `phys_type::einstein_frame`:
//! the metric, KMS Eq. (110), with the Palatini `kappa = 0` added, and the
//! Gaussian curvature `K = R / 2`, which is half of KMS Eq. (115) as printed
//! (`docs/literature-discrepancies.md` entry 14).
//! See `docs/general-nonminimal-inflation.md`.

use phys_type::{
    Dimensionless, EnergyDensity, Field, Formulation, GeneralizedNonminimal, Quantity, TheoryScalar,
};

/// The two-field model of KMS Eqs. (16)--(17), in reduced Planck units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KmsTheory {
    /// The nonminimal coupling `xi_phi`.
    pub xi_phi: f64,
    /// The nonminimal coupling `xi_chi`.
    pub xi_chi: f64,
    /// The quartic coupling `lambda_phi`.
    pub lambda_phi: f64,
    /// The quartic coupling `lambda_chi`.
    pub lambda_chi: f64,
    /// The cross coupling `g` of `g phi^2 chi^2 / 2`.
    pub g: f64,
    /// The squared mass `m_phi^2`.
    pub m_phi_squared: f64,
    /// The squared mass `m_chi^2`.
    pub m_chi_squared: f64,
    /// The metric or Palatini formulation.
    pub formulation: Formulation,
}

impl KmsTheory {
    /// The benchmark of KMS Figs. 2--3 and 6--11: `xi_phi = 10`,
    /// `xi_chi = 10.02`, `lambda_chi = lambda_phi / 2`, `g = lambda_phi`,
    /// and no mass terms. `lambda_phi` sets only the amplitude.
    pub fn benchmark(lambda_phi: f64, formulation: Formulation) -> Self {
        Self {
            xi_phi: 10.0,
            xi_chi: 10.02,
            lambda_phi,
            lambda_chi: 0.5 * lambda_phi,
            g: lambda_phi,
            m_phi_squared: 0.0,
            m_chi_squared: 0.0,
            formulation,
        }
    }

    /// `Omega^2 = 1 + xi_phi phi^2 + xi_chi chi^2`.
    fn omega2(&self, phi: f64, chi: f64) -> f64 {
        1.0 + self.xi_phi * phi * phi + self.xi_chi * chi * chi
    }

    /// The Einstein-frame metric `G_IJ = delta_IJ / Omega^2
    /// + 6 kappa xi_I xi_J phi^I phi^J / Omega^4`, KMS Eq. (110) at
    /// `kappa = 1`.
    pub fn metric(&self, phi: f64, chi: f64) -> [[f64; 2]; 2] {
        let omega2 = self.omega2(phi, chi);
        let weight = 6.0 * self.formulation.kappa() / (omega2 * omega2);
        let (a, b) = (self.xi_phi * phi, self.xi_chi * chi);
        [
            [1.0 / omega2 + weight * a * a, weight * a * b],
            [weight * a * b, 1.0 / omega2 + weight * b * b],
        ]
    }

    /// The Gaussian curvature `K`, with `R^I_JKL = K (delta^I_K G_JL -
    /// delta^I_L G_JK)`.
    ///
    /// Metric: `K = [(1 + 6 xi_phi)(1 + 6 xi_chi) Omega^4 - C^2] / (6 C^2)`
    /// with `C = Omega^2 + 6 xi_phi^2 phi^2 + 6 xi_chi^2 chi^2`, KMS
    /// Eq. (112); KMS Eq. (115) prints twice this. Palatini: `G = delta / Omega^2` is conformally
    /// flat, so `K = (Omega^2 / 2) Laplacian ln Omega^2
    /// = xi_phi + xi_chi - 2 (xi_phi^2 phi^2 + xi_chi^2 chi^2) / Omega^2`.
    pub fn gaussian_curvature(&self, phi: f64, chi: f64) -> f64 {
        let omega2 = self.omega2(phi, chi);
        let gradient =
            self.xi_phi * self.xi_phi * phi * phi + self.xi_chi * self.xi_chi * chi * chi;
        match self.formulation {
            Formulation::Metric => {
                let c = omega2 + 6.0 * gradient;
                let product = (1.0 + 6.0 * self.xi_phi) * (1.0 + 6.0 * self.xi_chi);
                (product * omega2 * omega2 - c * c) / (6.0 * c * c)
            }
            Formulation::Palatini => self.xi_phi + self.xi_chi - 2.0 * gradient / omega2,
        }
    }
}

impl GeneralizedNonminimal<2> for KmsTheory {
    fn formulation(&self) -> Formulation {
        self.formulation
    }

    fn conformal_factor<V: TheoryScalar>(&self, fields: [Field<V>; 2]) -> Dimensionless<V> {
        let [phi, chi] = fields.map(Quantity::into_value);
        Quantity::new(V::from(1.0) + phi * phi * self.xi_phi + chi * chi * self.xi_chi)
    }

    fn kinetic_matrix<V: TheoryScalar>(
        &self,
        _fields: [Field<V>; 2],
    ) -> [[Dimensionless<V>; 2]; 2] {
        let (one, zero) = (Quantity::new(V::from(1.0)), Quantity::new(V::from(0.0)));
        [[one, zero], [zero, one]]
    }

    fn jordan_potential<V: TheoryScalar>(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        let [phi, chi] = fields.map(Quantity::into_value);
        let (phi2, chi2) = (phi * phi, chi * chi);
        EnergyDensity::new(
            phi2 * (0.5 * self.m_phi_squared)
                + chi2 * (0.5 * self.m_chi_squared)
                + phi2 * chi2 * (0.5 * self.g)
                + phi2 * phi2 * (0.25 * self.lambda_phi)
                + chi2 * chi2 * (0.25 * self.lambda_chi),
        )
    }
}
