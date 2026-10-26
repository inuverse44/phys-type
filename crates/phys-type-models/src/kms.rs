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
    BackgroundSurvey, CovariantMultifieldSystem, Dimensionless, EinsteinGeometry,
    EinsteinPotential, EnergyDensity, Field, FieldSpaceGeometry, FieldVelocity, Formulation,
    GeneralizedNonminimal, MultifieldState, Quantity, ScaleFactor, SpectrumRunError,
    SpectrumRunSettings, TheoryScalar, einstein_frame,
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

/// The Einstein-frame system of a [`KmsTheory`].
pub type KmsSystem =
    CovariantMultifieldSystem<EinsteinPotential<KmsTheory>, EinsteinGeometry<KmsTheory>, 2>;

/// `phi(tau_0) / M_pl` of the KMS trajectories, Fig. 2.
pub const PHI_START: f64 = 3.1;

/// `chi(tau_0) / M_pl` of trajectories 1--3, KMS Fig. 2.
pub const TRAJECTORIES: [f64; 3] = [1.1e-2, 1.1e-3, 1.1e-4];

/// The benchmark system with `lambda_phi = 10^-2`, which sets only the
/// amplitude of the spectra.
pub fn benchmark_system(formulation: Formulation) -> KmsSystem {
    let (potential, geometry) = einstein_frame(KmsTheory::benchmark(1.0e-2, formulation));
    CovariantMultifieldSystem::in_planck_units(potential, geometry)
}

/// The fields `(phi, chi)` at rest. KMS do not state the initial velocities
/// of trajectories 1--3; starting from rest reproduces their 71.5 e-folds of
/// Fig. 3, where the slow-roll velocity gives 71.2.
pub fn at_rest(phi: f64, chi: f64) -> MultifieldState<2> {
    MultifieldState::new(
        [Field::new(phi), Field::new(chi)],
        [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
        ScaleFactor::new(1.0),
    )
}

/// The crossing-time delta N of KMS Sec. V.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrossingDeltaN {
    /// `f_NL` in the convention of this repository; KMS Eq. (92) carries the
    /// opposite sign.
    pub f_nl: f64,
    /// `P_zeta = (H_* / 2 pi)^2 N_a N^a`.
    pub power: f64,
}

/// The delta N of KMS Sec. V: `N` to the end of inflation from the field
/// values at the crossing of the pivot, `pivot_before_end` e-folds before
/// the end of `survey`, with slow-roll velocities. The derivatives are
/// central differences with steps `3e-4` of each field, and the second
/// derivative is covariant, `N_;ab = N_,ab - Gamma^c_ab N_,c`.
pub fn crossing_delta_n(
    system: &KmsSystem,
    survey: &BackgroundSurvey<2>,
    settings: &SpectrumRunSettings,
    pivot_before_end: f64,
) -> Result<CrossingDeltaN, SpectrumRunError> {
    let star = survey
        .trajectory
        .state_at(survey.e_folds_at_end - pivot_before_end);
    let x0 = [star.phi[0].into_value(), star.phi[1].into_value()];
    let quiet = SpectrumRunSettings {
        background_sample_every: usize::MAX,
        ..*settings
    };
    let h = [3.0e-4 * x0[0].abs(), 3.0e-4 * x0[1].abs()];
    let e_folds = |a: f64, b: f64| {
        let rest = at_rest(x0[0] + a * h[0], x0[1] + b * h[1]);
        let state = MultifieldState {
            pi: system.slow_roll_velocity(&rest),
            ..rest
        };
        system.survey(state, &quiet).map(|run| run.e_folds_at_end)
    };
    let centre = e_folds(0.0, 0.0)?;
    let axis = [
        [e_folds(1.0, 0.0)?, e_folds(-1.0, 0.0)?],
        [e_folds(0.0, 1.0)?, e_folds(0.0, -1.0)?],
    ];
    let mixed = (e_folds(1.0, 1.0)? - e_folds(1.0, -1.0)? - e_folds(-1.0, 1.0)?
        + e_folds(-1.0, -1.0)?)
        / (4.0 * h[0] * h[1]);
    let first: [f64; 2] = core::array::from_fn(|a| (axis[a][0] - axis[a][1]) / (2.0 * h[a]));
    let partial: [[f64; 2]; 2] = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            if a == b {
                (axis[a][0] - 2.0 * centre + axis[a][1]) / (h[a] * h[a])
            } else {
                mixed
            }
        })
    });
    let position = [Field::new(x0[0]), Field::new(x0[1])];
    let connection = system.geometry.christoffel(&position);
    let inverse = system.geometry.inverse_metric(&position);
    let covariant: [[f64; 2]; 2] = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            partial[a][b]
                - (0..2)
                    .map(|c| connection.components()[c][a][b].value() * first[c])
                    .sum::<f64>()
        })
    });
    let raised: [f64; 2] = core::array::from_fn(|a| {
        (0..2)
            .map(|b| inverse.components()[a][b].into_value() * first[b])
            .sum()
    });
    let norm: f64 = (0..2).map(|a| raised[a] * first[a]).sum();
    let curvature: f64 = (0..4)
        .map(|m| raised[m / 2] * raised[m % 2] * covariant[m / 2][m % 2])
        .sum();
    let hubble = system.hubble(&star).into_value();
    Ok(CrossingDeltaN {
        f_nl: 5.0 / 6.0 * curvature / (norm * norm),
        power: (hubble / (2.0 * core::f64::consts::PI)).powi(2) * norm,
    })
}

/// The turn as the fields fall off the ridge: the first local maximum of
/// `|eta_perp|` above `0.1`, with the e-folds before the end at which it
/// occurs, or `None` if the fields never turn that sharply. The heavy
/// oscillation of `phi` in the valley afterwards drives `|eta_perp|` far
/// higher, but in spikes; the first maximum is the peak of KMS Fig. 7.
pub fn turn_peak(survey: &BackgroundSurvey<2>) -> Option<(f64, f64)> {
    let turns: Vec<(f64, f64)> = survey
        .samples
        .iter()
        .filter(|sample| sample.eta_perp.is_finite())
        .map(|sample| (sample.eta_perp, survey.e_folds_at_end - sample.e_folds))
        .collect();
    turns
        .windows(3)
        .find(|w| w[1].0 > 0.1 && w[1].0 >= w[0].0 && w[1].0 >= w[2].0)
        .map(|w| w[1])
}
