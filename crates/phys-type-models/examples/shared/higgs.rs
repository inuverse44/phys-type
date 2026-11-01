//! The exact Higgs-inflation prediction shared by the Higgs examples.

use phys_type::{
    CovariantMultifieldSystem, EinsteinFrame, Field, FieldVelocity, Formulation, Mass,
    MultifieldState, ScaleFactor, SpectrumRunSettings, einstein_frame,
};
use phys_type_models::higgs::HiggsTheory;

/// E-folds before `epsilon_H = 1` at which the spectra are read.
pub const OUTPUT_BEFORE_END: f64 = 1.0;
/// Logarithmic step in `k` of the central difference for `n_s`.
const TILT_STEP: f64 = 0.02;

/// The exact and first-order slow-roll predictions at one pivot.
pub struct Prediction {
    pub n_s: f64,
    pub r: f64,
    pub n_s_slow_roll: f64,
    pub r_slow_roll: f64,
    /// Jordan-frame e-folds over the pivot-to-end interval.
    pub jordan_e_folds: f64,
    pub phi_star: f64,
}

/// Evolves Higgs inflation with coupling `xi` exactly and returns `n_s` and
/// `r` for the pivot `pivot_before_end` Einstein-frame e-folds before the end.
pub fn predict(
    xi: f64,
    formulation: Formulation,
    pivot_before_end: f64,
) -> Result<Prediction, String> {
    // lambda only sets the amplitude; choose it so that H stays small.
    let theory = HiggsTheory {
        xi,
        lambda: 1.0e-12 * (1.0 + xi * xi),
        formulation,
    };
    let (potential, geometry) = einstein_frame(theory);
    let system = CovariantMultifieldSystem::in_planck_units(potential, geometry);
    let phi0 = theory.field_at_e_folds(pivot_before_end + 12.0);
    let resting = MultifieldState::new(
        [Field::new(phi0)],
        [FieldVelocity::new(0.0)],
        ScaleFactor::new(1.0),
    );
    let initial = MultifieldState {
        pi: system.slow_roll_velocity(&resting),
        ..resting
    };
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: pivot_before_end,
        output_e_folds_before_end: OUTPUT_BEFORE_END,
        e_folds_before_exit: 6.0,
        ..SpectrumRunSettings::default()
    };
    let survey = system
        .survey(initial, &settings)
        .map_err(|error| error.to_string())?;
    let crossing = survey.e_folds_at_end - pivot_before_end;
    let output = survey.e_folds_at_end - OUTPUT_BEFORE_END;
    let pivot = survey.trajectory.state_at(crossing);
    let k = (pivot.a * system.hubble(&pivot)).into_value();
    let curvature = |kk: f64| {
        system
            .mode_spectrum(&settings, &survey.trajectory, Mass::new(kk), output)
            .map(|spectrum| spectrum.spectra.curvature.into_value())
    };
    let (Ok(low), Ok(central), Ok(high)) = (
        curvature(k * (-TILT_STEP).exp()),
        curvature(k),
        curvature(k * TILT_STEP.exp()),
    ) else {
        return Err("a mode failed".to_string());
    };
    let n_s = 1.0 + (high / low).ln() / (2.0 * TILT_STEP);
    let tensor = system
        .tensor_spectrum(&settings, &survey.trajectory, Mass::new(k), output)
        .map_err(|error| error.to_string())?;
    let (n_s_slow_roll, r_slow_roll) = theory.slow_roll_predictions(pivot_before_end);
    let frame = EinsteinFrame(theory);
    let end = survey.trajectory.state_at(survey.e_folds_at_end);
    let jordan = |state: &MultifieldState<1>| {
        frame
            .correspondence(state, system.hubble(state))
            .expect("Omega^2 > 0")
            .jordan_e_folds
    };
    Ok(Prediction {
        n_s,
        r: tensor / central,
        n_s_slow_roll,
        r_slow_roll,
        jordan_e_folds: jordan(&end) - jordan(&pivot),
        phi_star: pivot.phi[0].into_value(),
    })
}
