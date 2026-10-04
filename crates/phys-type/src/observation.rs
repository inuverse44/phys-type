//! Observational provenance at the boundary of the primordial solver.
//!
//! This module names an observational reference; it does not evaluate a CMB
//! likelihood.  In particular, a wavenumber in `Mpc^-1` is kept distinct from
//! the solver's normalization-dependent internal Planck-unit wavenumber.

use std::fmt;
use std::io::{self, Write};

/// A positive comoving wavenumber expressed in inverse megaparsecs.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct ComovingWavenumberMpcInverse(f64);

impl ComovingWavenumberMpcInverse {
    /// Builds an observational comoving wavenumber.
    pub fn try_new(value: f64) -> Result<Self, InvalidObservationalWavenumber> {
        if value.is_finite() && value > 0.0 {
            Ok(Self(value))
        } else {
            Err(InvalidObservationalWavenumber { value })
        }
    }

    const fn new_unchecked(value: f64) -> Self {
        Self(value)
    }

    /// Returns the numerical value in `Mpc^-1`.
    pub const fn into_value(self) -> f64 {
        self.0
    }
}

/// Error returned for a non-positive or non-finite observational wavenumber.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InvalidObservationalWavenumber {
    value: f64,
}

impl fmt::Display for InvalidObservationalWavenumber {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "observational wavenumber must be positive and finite, got {} Mpc^-1",
            self.value
        )
    }
}

impl std::error::Error for InvalidObservationalWavenumber {}

/// A named observational data release.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObservationDataRelease {
    /// The final 2018 Planck legacy release (PR3).
    Planck2018LegacyPr3,
}

impl ObservationDataRelease {
    /// Stable machine-readable identifier.
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::Planck2018LegacyPr3 => "planck_2018_legacy_pr3",
        }
    }
}

/// One CMB likelihood component named by its physical role.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LikelihoodComponent {
    /// High-multipole Plik TT, TE, and EE spectra.
    PlikHighLTTTEEE,
    /// Low-multipole Commander temperature likelihood.
    CommanderLowLTT,
    /// Low-multipole SimAll EE likelihood (`lowE`).
    SimAllLowLEE,
    /// Planck lensing-reconstruction likelihood, when explicitly requested.
    Lensing,
}

impl LikelihoodComponent {
    /// Stable machine-readable identifier.
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::PlikHighLTTTEEE => "plik_high_l_tt_te_ee",
            Self::CommanderLowLTT => "commander_low_l_tt",
            Self::SimAllLowLEE => "simall_low_l_ee",
            Self::Lensing => "lensing_reconstruction",
        }
    }
}

/// How an amplitude appearing in an output relates to the named likelihood.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AmplitudeRole {
    /// A fixed numerical target used to calibrate a benchmark, not a fit.
    BenchmarkTargetNotLikelihoodFit,
    /// An amplitude inferred by an external likelihood analysis.
    LikelihoodFit,
}

impl AmplitudeRole {
    /// Stable machine-readable identifier.
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::BenchmarkTargetNotLikelihoodFit => "benchmark_target_not_likelihood_fit",
            Self::LikelihoodFit => "likelihood_fit",
        }
    }
}

/// Complete provenance needed to identify an observational reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservationalReference {
    /// Stable name for this precise reference combination.
    pub identifier: &'static str,
    /// The observational data release.
    pub data_release: ObservationDataRelease,
    /// The released likelihood package/version.
    pub likelihood_release: &'static str,
    /// Scalar-spectrum pivot in observational comoving units.
    pub scalar_pivot: ComovingWavenumberMpcInverse,
    /// Likelihood components included in the reference combination.
    pub components: &'static [LikelihoodComponent],
}

const PLANCK_2018_BASELINE_COMPONENTS: &[LikelihoodComponent] = &[
    LikelihoodComponent::PlikHighLTTTEEE,
    LikelihoodComponent::CommanderLowLTT,
    LikelihoodComponent::SimAllLowLEE,
];

/// Planck 2018 PR3 baseline `TT,TE,EE+lowE` reference.
///
/// Lensing reconstruction is deliberately absent: callers must name it as an
/// additional likelihood component rather than silently treating it as part
/// of this baseline.
pub const PLANCK_2018_TT_TE_EE_LOW_E: ObservationalReference = ObservationalReference {
    identifier: "planck_2018_tt_te_ee_low_e",
    data_release: ObservationDataRelease::Planck2018LegacyPr3,
    likelihood_release: "COM_Likelihood_Data-baseline_R3.00",
    scalar_pivot: ComovingWavenumberMpcInverse::new_unchecked(0.05),
    components: PLANCK_2018_BASELINE_COMPONENTS,
};

/// Writes machine-readable observational provenance as CSV comments.
pub fn write_observational_metadata<W: Write>(
    writer: &mut W,
    reference: &ObservationalReference,
    amplitude_role: AmplitudeRole,
) -> io::Result<()> {
    writeln!(
        writer,
        "# observational_reference: {}",
        reference.identifier
    )?;
    writeln!(
        writer,
        "# observational_data_release: {}",
        reference.data_release.identifier()
    )?;
    writeln!(
        writer,
        "# likelihood_release: {}",
        reference.likelihood_release
    )?;
    writeln!(
        writer,
        "# likelihood_components: {}",
        reference
            .components
            .iter()
            .map(|component| component.identifier())
            .collect::<Vec<_>>()
            .join(",")
    )?;
    writeln!(
        writer,
        "# scalar_pivot_mpc_inverse: {:e}",
        reference.scalar_pivot.into_value()
    )?;
    writeln!(writer, "# amplitude_role: {}", amplitude_role.identifier())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn observational_wavenumber_requires_a_positive_finite_value() {
        assert!(ComovingWavenumberMpcInverse::try_new(0.05).is_ok());
        assert!(ComovingWavenumberMpcInverse::try_new(0.0).is_err());
        assert!(ComovingWavenumberMpcInverse::try_new(-0.05).is_err());
        assert!(ComovingWavenumberMpcInverse::try_new(f64::NAN).is_err());
    }

    #[test]
    fn planck_baseline_is_explicit_and_excludes_lensing() {
        let reference = PLANCK_2018_TT_TE_EE_LOW_E;
        assert_eq!(reference.scalar_pivot.into_value(), 0.05);
        assert_eq!(reference.components, PLANCK_2018_BASELINE_COMPONENTS);
        assert!(!reference.components.contains(&LikelihoodComponent::Lensing));
    }

    #[test]
    fn metadata_states_release_components_pivot_and_amplitude_role() {
        let mut output = Vec::new();
        write_observational_metadata(
            &mut output,
            &PLANCK_2018_TT_TE_EE_LOW_E,
            AmplitudeRole::BenchmarkTargetNotLikelihoodFit,
        )
        .unwrap();
        let text = String::from_utf8(output).unwrap();

        assert!(text.contains("# observational_data_release: planck_2018_legacy_pr3"));
        assert!(text.contains("plik_high_l_tt_te_ee"));
        assert!(text.contains("commander_low_l_tt"));
        assert!(text.contains("simall_low_l_ee"));
        assert!(text.contains("# scalar_pivot_mpc_inverse: 5e-2"));
        assert!(text.contains("# amplitude_role: benchmark_target_not_likelihood_fit"));
    }
}
