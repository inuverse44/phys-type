//! Adiabatic and entropy power spectra built from a two-field mode matrix.
//!
//! The mode covariance is complex and carries two upper field-space indices.
//! Projecting it onto the adiabatic and entropy directions requires lowered
//! basis covectors, so an accidental same-variance contraction does not
//! compile. A projected spectrum becomes a real number only after its
//! imaginary residual has been checked against an explicit tolerance.

use core::f64::consts::PI;

use num_complex::Complex64;

use crate::{
    ComplexMagnitude, ContravariantField, CovariantField, CurvedFieldSpaceGeometry, Dimensionless,
    FieldMetric, FieldTensor2, Hubble, KinematicBasis2, KinematicBasisError, M0, M2, MMinus1,
    MMinusHalf, Mass, ModeAmplitude, Quantity, ScaleFactor, TwiceDifferentiableMultifieldPotential,
    TwoFieldModes, TwoFieldPerturbationState, TwoFieldPerturbationSystem, Upper,
};

/// Mode covariance `Sigma^{IJ} = sum_A U^I_A U^{J*}_A`.
pub type ModeCovariance = FieldTensor2<Upper, Upper, MMinus1, Complex64, 2>;

/// Field power spectrum matrix `P_Q^{IJ} = k^3 Sigma^{IJ} / (2 pi^2)`.
pub type FieldPowerSpectrum = FieldTensor2<Upper, Upper, M2, Complex64, 2>;

/// A dimensionless spectrum of the curvature or entropy perturbation.
pub type DimensionlessSpectrum = Dimensionless;

/// Tolerances for the reality and positivity checks that guard a projection.
///
/// Every tolerance is relative to the largest modulus in the spectrum matrix,
/// so the same settings apply at any overall normalization.
///
/// The three are not equally tight, because they measure different things.
/// Hermiticity holds by construction and fails only through rounding. The
/// imaginary residual instead measures accumulated integration error: the
/// antisymmetric part of `Sigma` is the equal-time field commutator, which
/// vanishes at Bunch-Davies initialization and is conserved by the mode
/// equations, so whatever survives to the output time was produced by the
/// integrator. Choose it together with the integration tolerances rather than
/// in isolation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectrumTolerances {
    /// Largest accepted relative violation of `P^{IJ} = P^{JI*}`.
    pub hermiticity: f64,
    /// Largest accepted relative imaginary part of a projected spectrum.
    pub imaginary_residual: f64,
    /// Largest accepted overshoot of `|r_C|` beyond one.
    pub correlation_overshoot: f64,
}

impl Default for SpectrumTolerances {
    fn default() -> Self {
        Self {
            hermiticity: 1.0e-8,
            imaginary_residual: 1.0e-6,
            correlation_overshoot: 1.0e-8,
        }
    }
}

/// Why a mode covariance could not be turned into physical spectra.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpectrumError {
    /// The spectrum matrix failed the Hermiticity check.
    NotHermitian {
        /// The observed relative violation.
        residual: f64,
        /// The tolerance it exceeded.
        tolerance: f64,
    },
    /// A projected spectrum kept an imaginary part too large to discard.
    ImaginaryResidual {
        /// The observed relative imaginary part.
        residual: f64,
        /// The tolerance it exceeded.
        tolerance: f64,
    },
    /// An autocorrelation spectrum was not positive, so no correlation
    /// coefficient is defined.
    NonPositiveSpectrum {
        /// `P_R`, which must be positive.
        curvature: f64,
        /// `P_S`, which must be positive.
        entropy: f64,
    },
    /// The correlation coefficient left `[-1, 1]` by more than the tolerance.
    CorrelationOutOfRange {
        /// The unclamped coefficient.
        value: f64,
        /// The tolerance it exceeded.
        tolerance: f64,
    },
    /// The adiabatic and entropy directions could not be constructed.
    KinematicBasis(KinematicBasisError),
}

/// Spectra projected onto the adiabatic and entropy directions.
///
/// The field spectra keep the dimension `M^2`. The curvature and entropy
/// spectra are dimensionless because each carries two factors of
/// `H / dot_sigma`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProjectedSpectra {
    /// `P_{Q_sigma}`, the adiabatic field spectrum.
    pub field_adiabatic: Quantity<M2, f64>,
    /// `P_{Q_s}`, the entropy field spectrum.
    pub field_entropy: Quantity<M2, f64>,
    /// `C_{Q_sigma Q_s}`, the signed field cross spectrum.
    pub field_cross: Quantity<M2, f64>,
    /// `P_R`, the curvature perturbation spectrum.
    pub curvature: DimensionlessSpectrum,
    /// `P_S`, the entropy perturbation spectrum.
    pub entropy: DimensionlessSpectrum,
    /// `C_RS`, the signed curvature-entropy cross spectrum.
    pub cross: DimensionlessSpectrum,
    /// `r_C = C_RS / sqrt(P_R P_S)`, clamped to `[-1, 1]`.
    pub correlation: f64,
}

/// Builds `Sigma^{IJ} = sum_A U^I_A U^{J*}_A` from a perturbation state.
///
/// The modes are integrated as `v^I_A = a U^I_A`, so this is the one place
/// where the rescaling is undone: the sum carries `1 / a^2`. It runs over the
/// independent vacuum oscillators, so the result is Hermitian by construction
/// up to floating-point rounding. Taking the complete state prevents modes
/// from accidentally being paired with a scale factor from another time.
pub fn mode_covariance(state: &TwoFieldPerturbationState) -> ModeCovariance {
    mode_covariance_from_modes(&state.modes, state.background.a)
}

fn mode_covariance_from_modes(modes: &TwoFieldModes, scale_factor: ScaleFactor) -> ModeCovariance {
    let squared = scale_factor * scale_factor;

    ModeCovariance::new(core::array::from_fn(|row| {
        core::array::from_fn(|column| {
            let rescaled = (0..2).fold(
                Quantity::<MMinus1, Complex64>::new(Complex64::new(0.0, 0.0)),
                |sum, mode| {
                    let conjugate =
                        ModeAmplitude::new(modes.rescaled_amplitude[column][mode].value().conj());

                    sum + modes.rescaled_amplitude[row][mode] * conjugate
                },
            );

            rescaled / squared
        })
    }))
}

/// Applies the `k^3 / (2 pi^2)` normalization to a mode covariance.
pub fn field_power_spectrum(k: Mass, covariance: &ModeCovariance) -> FieldPowerSpectrum {
    let normalization = k * k * k / (2.0 * PI * PI);

    FieldPowerSpectrum::new(
        covariance
            .components()
            .map(|row| row.map(|entry| normalization * entry)),
    )
}

/// The largest modulus in a spectrum matrix, used as the scale for every
/// relative tolerance.
pub(crate) fn spectrum_scale(spectrum: &FieldPowerSpectrum) -> f64 {
    spectrum
        .components()
        .iter()
        .flatten()
        .map(ComplexMagnitude::complex_magnitude)
        .fold(0.0, f64::max)
}

/// The relative violation of `P^{IJ} = P^{JI*}`.
///
/// A matrix whose entries all vanish is Hermitian, and reports zero rather
/// than dividing by an empty scale.
pub fn hermiticity_residual(spectrum: &FieldPowerSpectrum) -> f64 {
    let scale = spectrum_scale(spectrum);
    if scale == 0.0 {
        return 0.0;
    }

    let components = spectrum.components();
    let mut residual: f64 = 0.0;
    for (row, entries) in components.iter().enumerate() {
        for (column, entry) in entries.iter().enumerate() {
            let difference = entry.value() - components[column][row].value().conj();
            residual = residual.max(difference.norm() / scale);
        }
    }

    residual
}

/// Contracts a spectrum matrix with two lowered basis covectors.
///
/// Both indices of `P^{IJ}` are upper, so only covectors can be supplied.
fn project_complex(
    spectrum: &FieldPowerSpectrum,
    left: &CovariantField<M0, f64, 2>,
    right: &CovariantField<M0, f64, 2>,
) -> Quantity<M2, Complex64> {
    let components = spectrum.components();
    let left = left.components();
    let right = right.components();

    (0..2)
        .flat_map(|row| (0..2).map(move |column| (row, column)))
        .fold(
            Quantity::<M2, Complex64>::new(Complex64::new(0.0, 0.0)),
            |sum, (row, column)| {
                let weight = left[row].value() * right[column].value();

                sum + components[row][column] * weight
            },
        )
}

/// Discards the imaginary part of a projected spectrum after checking that it
/// is negligible relative to `scale`.
pub(crate) fn discard_imaginary(
    value: Quantity<M2, Complex64>,
    scale: f64,
    tolerance: f64,
) -> Result<Quantity<M2, f64>, SpectrumError> {
    let value = value.into_value();
    if scale > 0.0 {
        let residual = value.im.abs() / scale;
        if residual > tolerance {
            return Err(SpectrumError::ImaginaryResidual {
                residual,
                tolerance,
            });
        }
    }

    Ok(Quantity::new(value.re))
}

/// The dimensionless curvature and entropy amplitude of each independent
/// mode, as `[[R_1, R_2], [S_1, S_2]]`.
///
/// Summing `|R_A|^2` over the independent modes reproduces `P_R`, and
/// likewise for the entropy row, so this is the same content as
/// [`ProjectedSpectra`] before the sum over vacuum oscillators discards which
/// solution each contribution came from. The super-horizon transfer matrix
/// needs that information, because it relates two independent solutions
/// rather than their sum.
pub type ModeProjections = [[Complex64; 2]; 2];

/// Projects each independent mode onto the curvature and entropy directions.
fn project_modes_in_basis(
    modes: &TwoFieldModes,
    scale_factor: ScaleFactor,
    metric: &FieldMetric<2>,
    basis: &KinematicBasis2,
    hubble: Hubble,
    k: Mass,
) -> ModeProjections {
    let adiabatic = metric.lower(&basis.tangent);
    let entropy = metric.lower(&basis.normal);
    let directions = [adiabatic, entropy];
    // R = (H / dot_sigma) Q_sigma, and P_R = (k^3 / 2 pi^2) |R|^2 summed over
    // modes, so each amplitude carries k^(3/2) / (pi sqrt 2).
    let normalization = (k * k * k).sqrt() / (PI * 2.0_f64.sqrt());
    let transfer = hubble / basis.speed;

    core::array::from_fn(|direction| {
        core::array::from_fn(|mode| {
            let projected = (0..2).fold(
                Quantity::<MMinusHalf, Complex64>::new(Complex64::new(0.0, 0.0)),
                |sum, field| {
                    let weight = *directions[direction].components()[field].value();
                    let amplitude = modes.rescaled_amplitude[field][mode] / scale_factor;

                    sum + amplitude * weight
                },
            );

            (normalization * transfer * projected).into_value()
        })
    })
}

/// Projects a spectrum matrix onto an adiabatic and entropy basis.
///
/// `metric` lowers the basis vectors, and `hubble` together with the basis
/// speed `dot_sigma` supplies the `H / dot_sigma` normalization that makes the
/// curvature and entropy spectra dimensionless.
pub fn project(
    spectrum: &FieldPowerSpectrum,
    metric: &FieldMetric<2>,
    basis: &KinematicBasis2,
    hubble: Hubble,
    tolerances: SpectrumTolerances,
) -> Result<ProjectedSpectra, SpectrumError> {
    let residual = hermiticity_residual(spectrum);
    if residual > tolerances.hermiticity {
        return Err(SpectrumError::NotHermitian {
            residual,
            tolerance: tolerances.hermiticity,
        });
    }

    let adiabatic = metric.lower(&basis.tangent);
    let entropy = metric.lower(&basis.normal);
    let scale = spectrum_scale(spectrum);

    let field_adiabatic = discard_imaginary(
        project_complex(spectrum, &adiabatic, &adiabatic),
        scale,
        tolerances.imaginary_residual,
    )?;
    let field_entropy = discard_imaginary(
        project_complex(spectrum, &entropy, &entropy),
        scale,
        tolerances.imaginary_residual,
    )?;
    let field_cross = discard_imaginary(
        project_complex(spectrum, &adiabatic, &entropy),
        scale,
        tolerances.imaginary_residual,
    )?;

    // R = (H / dot_sigma) Q_sigma, so each spectrum picks up two factors.
    let transfer = hubble / basis.speed;
    let normalization = transfer * transfer;
    let curvature = normalization * field_adiabatic;
    let entropy_spectrum = normalization * field_entropy;
    let cross = normalization * field_cross;

    let correlation = correlation_coefficient(
        curvature,
        entropy_spectrum,
        cross,
        tolerances.correlation_overshoot,
    )?;

    Ok(ProjectedSpectra {
        field_adiabatic,
        field_entropy,
        field_cross,
        curvature,
        entropy: entropy_spectrum,
        cross,
        correlation,
    })
}

/// Forms `r_C = C_RS / sqrt(P_R P_S)` and clamps a tolerable overshoot.
fn correlation_coefficient(
    curvature: DimensionlessSpectrum,
    entropy: DimensionlessSpectrum,
    cross: DimensionlessSpectrum,
    tolerance: f64,
) -> Result<f64, SpectrumError> {
    let curvature = curvature.into_value();
    let entropy = entropy.into_value();
    if !(curvature > 0.0 && entropy > 0.0) {
        return Err(SpectrumError::NonPositiveSpectrum { curvature, entropy });
    }

    let value = cross.into_value() / (curvature * entropy).sqrt();
    if value.abs() > 1.0 + tolerance {
        return Err(SpectrumError::CorrelationOutOfRange { value, tolerance });
    }

    Ok(value.clamp(-1.0, 1.0))
}

impl<P, G> TwoFieldPerturbationSystem<P, G>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    /// Projects every independent solution onto the instantaneous curvature
    /// and entropy directions.
    ///
    /// The system and complete state supply `k`, `a`, the metric, the basis,
    /// and `H` together, so values belonging to different evolution times
    /// cannot be combined accidentally.
    pub fn mode_projections(
        &self,
        state: &TwoFieldPerturbationState,
        previous_normal: Option<&ContravariantField<M0, f64, 2>>,
    ) -> Result<ModeProjections, SpectrumError> {
        let kinematics = self
            .background
            .kinematics(&state.background, previous_normal)
            .map_err(SpectrumError::KinematicBasis)?;
        let metric = self.background.geometry.metric(&state.background.phi);
        let hubble = self.background.hubble(&state.background);

        Ok(project_modes_in_basis(
            &state.modes,
            state.background.a,
            &metric,
            &kinematics.basis,
            hubble,
            self.k,
        ))
    }

    /// Projects the spectra of one evolved state onto its instantaneous
    /// adiabatic and entropy directions.
    ///
    /// `previous_normal` carries the entropy orientation forward between
    /// output times, exactly as for the background kinematics.
    pub fn spectra(
        &self,
        state: &TwoFieldPerturbationState,
        previous_normal: Option<&ContravariantField<M0, f64, 2>>,
        tolerances: SpectrumTolerances,
    ) -> Result<ProjectedSpectra, SpectrumError> {
        let kinematics = self
            .background
            .kinematics(&state.background, previous_normal)
            .map_err(SpectrumError::KinematicBasis)?;
        let metric = self.background.geometry.metric(&state.background.phi);
        let hubble = self.background.hubble(&state.background);
        let spectrum = field_power_spectrum(self.k, &mode_covariance(state));

        project(&spectrum, &metric, &kinematics.basis, hubble, tolerances)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        CovariantMultifieldSystem, Field, FieldVelocity, FlatFieldSpace, ModeVelocity,
        MultifieldState, QuadraticMultifieldPotential, ScaleFactor, Time, rk4_step,
    };

    fn system() -> TwoFieldPerturbationSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace> {
        let background = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(1.5e-6)]),
            FlatFieldSpace,
        );
        TwoFieldPerturbationSystem::new(background, Mass::new(1.0e-2))
    }

    fn background() -> MultifieldState<2> {
        MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(-3.0e-6), FieldVelocity::new(-4.0e-6)],
            ScaleFactor::new(1.0),
        )
    }

    fn spectrum_from(entries: [[Complex64; 2]; 2]) -> FieldPowerSpectrum {
        FieldPowerSpectrum::new(entries.map(|row| row.map(Quantity::new)))
    }

    fn euclidean_basis(speed: FieldVelocity) -> KinematicBasis2 {
        KinematicBasis2 {
            tangent: ContravariantField::new([Quantity::new(1.0), Quantity::new(0.0)]),
            normal: ContravariantField::new([Quantity::new(0.0), Quantity::new(1.0)]),
            speed,
        }
    }

    #[test]
    fn mode_covariance_is_hermitian() {
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let spectrum = field_power_spectrum(system.k, &mode_covariance(&state));

        assert!(hermiticity_residual(&spectrum) < 1.0e-14);
    }

    #[test]
    fn bunch_davies_spectrum_matches_the_sub_horizon_vacuum() {
        let system = system();
        let initial = background();
        let state = system.bunch_davies_state(initial, 100.0).unwrap();
        let spectrum = field_power_spectrum(system.k, &mode_covariance(&state));

        // With U^I_A = e_A^I / (a sqrt(2k)), the covariance is the inverse
        // metric over 2 k a^2, so P^{IJ} = (k/a)^2 G^{IJ} / (4 pi^2).
        let physical_wavenumber = (system.k / state.background.a).into_value();
        let expected = physical_wavenumber.powi(2) / (4.0 * PI * PI);
        let components = spectrum.components();

        for (row, entries) in components.iter().enumerate() {
            for (column, entry) in entries.iter().enumerate() {
                let value = entry.value();
                let diagonal = expected * f64::from(row == column);
                assert!((value.re - diagonal).abs() < 1.0e-12 * expected);
                assert!(value.im.abs() < 1.0e-12 * expected);
            }
        }
    }

    #[test]
    fn projection_normalizes_with_the_squared_transfer_factor() {
        let spectrum = spectrum_from([
            [Complex64::new(4.0, 0.0), Complex64::new(1.0, 0.0)],
            [Complex64::new(1.0, 0.0), Complex64::new(9.0, 0.0)],
        ]);
        let hubble = Hubble::new(2.0);
        let speed = FieldVelocity::new(8.0);

        let projected = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(speed),
            hubble,
            SpectrumTolerances::default(),
        )
        .unwrap();

        let transfer = (hubble.into_value() / speed.into_value()).powi(2);
        assert!((projected.field_adiabatic.into_value() - 4.0).abs() < 1.0e-12);
        assert!((projected.field_entropy.into_value() - 9.0).abs() < 1.0e-12);
        assert!((projected.field_cross.into_value() - 1.0).abs() < 1.0e-12);
        assert!((projected.curvature.into_value() - 4.0 * transfer).abs() < 1.0e-15);
        assert!((projected.entropy.into_value() - 9.0 * transfer).abs() < 1.0e-15);
        assert!((projected.cross.into_value() - transfer).abs() < 1.0e-15);
        assert!((projected.correlation - 1.0 / 6.0).abs() < 1.0e-12);
    }

    #[test]
    fn a_negative_cross_spectrum_keeps_its_sign() {
        let spectrum = spectrum_from([
            [Complex64::new(4.0, 0.0), Complex64::new(-1.0, 0.0)],
            [Complex64::new(-1.0, 0.0), Complex64::new(9.0, 0.0)],
        ]);

        let projected = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(FieldVelocity::new(1.0)),
            Hubble::new(1.0),
            SpectrumTolerances::default(),
        )
        .unwrap();

        assert!(projected.field_cross.into_value() < 0.0);
        assert!(projected.cross.into_value() < 0.0);
        assert!((projected.correlation + 1.0 / 6.0).abs() < 1.0e-12);
    }

    #[test]
    fn a_non_hermitian_matrix_is_rejected_before_projection() {
        let spectrum = spectrum_from([
            [Complex64::new(4.0, 0.0), Complex64::new(1.0, 0.0)],
            [Complex64::new(2.0, 0.0), Complex64::new(9.0, 0.0)],
        ]);

        let error = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(FieldVelocity::new(1.0)),
            Hubble::new(1.0),
            SpectrumTolerances::default(),
        )
        .unwrap_err();

        assert!(matches!(error, SpectrumError::NotHermitian { .. }));
    }

    #[test]
    fn a_hermitian_imaginary_cross_term_is_reported_rather_than_silently_dropped() {
        // This matrix is Hermitian, so it passes the first check, but the
        // adiabatic-entropy projection keeps a large imaginary part.
        let spectrum = spectrum_from([
            [Complex64::new(4.0, 0.0), Complex64::new(0.0, 1.0)],
            [Complex64::new(0.0, -1.0), Complex64::new(9.0, 0.0)],
        ]);

        let error = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(FieldVelocity::new(1.0)),
            Hubble::new(1.0),
            SpectrumTolerances::default(),
        )
        .unwrap_err();

        assert!(matches!(error, SpectrumError::ImaginaryResidual { .. }));
    }

    #[test]
    fn a_correlation_beyond_one_is_rejected() {
        let spectrum = spectrum_from([
            [Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)],
            [Complex64::new(2.0, 0.0), Complex64::new(1.0, 0.0)],
        ]);

        let error = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(FieldVelocity::new(1.0)),
            Hubble::new(1.0),
            SpectrumTolerances::default(),
        )
        .unwrap_err();

        assert!(matches!(error, SpectrumError::CorrelationOutOfRange { .. }));
    }

    #[test]
    fn a_perfectly_correlated_spectrum_is_clamped_to_one() {
        let spectrum = spectrum_from([
            [Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0)],
            [Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0)],
        ]);

        let projected = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(FieldVelocity::new(1.0)),
            Hubble::new(1.0),
            SpectrumTolerances::default(),
        )
        .unwrap();

        assert_eq!(projected.correlation, 1.0);
    }

    #[test]
    fn an_empty_spectrum_has_no_correlation_coefficient() {
        let spectrum = spectrum_from([[Complex64::new(0.0, 0.0); 2]; 2]);

        let error = project(
            &spectrum,
            &FieldMetric::euclidean(),
            &euclidean_basis(FieldVelocity::new(1.0)),
            Hubble::new(1.0),
            SpectrumTolerances::default(),
        )
        .unwrap_err();

        assert!(matches!(error, SpectrumError::NonPositiveSpectrum { .. }));
    }

    #[test]
    fn the_imaginary_part_of_the_covariance_starts_at_zero_and_stays_small() {
        // The antisymmetric part of Sigma is the equal-time field commutator.
        // It vanishes exactly in the Bunch-Davies vacuum and is conserved, so
        // any growth measures integration error rather than physics.
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let initial = mode_covariance(&state);

        for row in initial.components() {
            for entry in row {
                assert_eq!(entry.value().im, 0.0);
            }
        }

        let advanced = rk4_step(&system, Time::new(0.0), &state, Time::new(1.0));
        let evolved = mode_covariance(&advanced);
        let scale = evolved
            .components()
            .iter()
            .flatten()
            .map(ComplexMagnitude::complex_magnitude)
            .fold(0.0, f64::max);

        for row in evolved.components() {
            for entry in row {
                assert!(entry.value().im.abs() < 1.0e-12 * scale);
            }
        }
    }

    #[test]
    fn evolved_modes_still_produce_projectable_spectra() {
        let system = system();
        let state = system.bunch_davies_state(background(), 100.0).unwrap();
        let advanced = rk4_step(&system, Time::new(0.0), &state, Time::new(1.0));

        let projected = system
            .spectra(&advanced, None, SpectrumTolerances::default())
            .unwrap();

        assert!(projected.curvature.into_value() > 0.0);
        assert!(projected.entropy.into_value() > 0.0);
        assert!(projected.correlation.abs() <= 1.0);
    }

    #[test]
    fn mode_velocities_do_not_enter_the_covariance() {
        let system = system();
        let mut state = system.bunch_davies_state(background(), 100.0).unwrap();
        let before = mode_covariance(&state);

        state.modes.rescaled_velocity = [[ModeVelocity::new(Complex64::new(7.0, -3.0)); 2]; 2];

        assert_eq!(mode_covariance(&state), before);
    }
}
