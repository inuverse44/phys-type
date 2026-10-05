//! Higgs inflation (arXiv:0710.3755), metric and Palatini, evolved exactly.
//!
//! For a range of nonminimal couplings `xi`, the theory
//! `Omega^2 = 1 + xi phi^2`, `V_J = lambda phi^4 / 4` is mapped to the
//! Einstein frame by `phys_type::einstein_frame` and evolved with one field.
//! The pivot leaves the horizon 60 Einstein-frame e-folds before `epsilon_H = 1`;
//! `n_s` comes from the mode solver and `r = P_t / P_R` from the tensor and
//! curvature spectra one e-fold before the end. Writes `higgs_inflation.csv`
//! into the directory given as the first argument (default `diagnostics/`),
//! next to the first-order slow-roll values from the same Einstein-frame
//! potential and the number of Jordan-frame e-folds over the same interval.

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::{
    CovariantMultifieldSystem, EinsteinFrame, Field, FieldVelocity, Formulation, Mass,
    MultifieldState, ScaleFactor, SpectrumRunSettings, einstein_frame,
};
use phys_type_models::higgs::HiggsTheory;

#[path = "shared/provenance.rs"]
mod provenance;

use provenance::revision;

const PIVOT_BEFORE_END: f64 = 60.0;
const OUTPUT_BEFORE_END: f64 = 1.0;
const TILT_STEP: f64 = 0.02;

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;
    let path = directory.join("higgs_inflation.csv");
    let mut file = fs::File::create(&path)?;
    writeln!(
        file,
        "# generator: phys-type-models {}",
        env!("CARGO_PKG_VERSION")
    )?;
    if let Some(revision) = revision() {
        writeln!(file, "# revision: {revision}")?;
    }
    writeln!(
        file,
        "# reference: arXiv:0710.3755 (metric), arXiv:0803.2664 (Palatini)"
    )?;
    writeln!(
        file,
        "# model: Omega^2 = 1 + xi phi^2, V_J = lambda phi^4 / 4, one field"
    )?;
    writeln!(
        file,
        "# pivot_e_folds_before_end: {PIVOT_BEFORE_END} (Einstein frame, to epsilon_H = 1)"
    )?;
    writeln!(file, "# output_e_folds_before_end: {OUTPUT_BEFORE_END}")?;
    writeln!(
        file,
        "formulation,xi,n_s,r,n_s_slow_roll,r_slow_roll,jordan_e_folds,phi_star"
    )?;

    let couplings = [1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4];
    for formulation in [Formulation::Metric, Formulation::Palatini] {
        for xi in couplings {
            // lambda only sets the amplitude; choose it so that H stays small.
            let theory = HiggsTheory {
                xi,
                lambda: 1.0e-12 * (1.0 + xi * xi),
                formulation,
            };
            let (potential, geometry) = einstein_frame(theory);
            let system = CovariantMultifieldSystem::in_planck_units(potential, geometry);
            let phi0 = theory.field_at_e_folds(PIVOT_BEFORE_END + 12.0);
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
                pivot_e_folds_before_end: PIVOT_BEFORE_END,
                output_e_folds_before_end: OUTPUT_BEFORE_END,
                e_folds_before_exit: 6.0,
                ..SpectrumRunSettings::default()
            };
            let survey = match system.survey(initial, &settings) {
                Ok(survey) => survey,
                Err(error) => {
                    eprintln!("{formulation:?} xi {xi}: {error}");
                    continue;
                }
            };
            let crossing = survey.e_folds_at_end - PIVOT_BEFORE_END;
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
                eprintln!("{formulation:?} xi {xi}: a mode failed");
                continue;
            };
            let n_s = 1.0 + (high / low).ln() / (2.0 * TILT_STEP);
            let tensor = system
                .tensor_spectrum(&settings, &survey.trajectory, Mass::new(k), output)
                .expect("the tensor mode follows the scalar modes");
            let r = tensor / central;
            let (n_s_slow_roll, r_slow_roll) = theory.slow_roll_predictions(PIVOT_BEFORE_END);
            let frame = EinsteinFrame(theory);
            let hubble = |state: &MultifieldState<1>| system.hubble(state);
            let end = survey.trajectory.state_at(survey.e_folds_at_end);
            let jordan = |state: &MultifieldState<1>| {
                frame
                    .correspondence(state, hubble(state))
                    .expect("Omega^2 > 0")
                    .jordan_e_folds
            };
            let jordan_e_folds = jordan(&end) - jordan(&pivot);
            writeln!(
                file,
                "{formulation:?},{xi:e},{n_s:.8},{r:.8e},{n_s_slow_roll:.8},{r_slow_roll:.8e},{jordan_e_folds:.6},{:.8e}",
                pivot.phi[0].into_value()
            )?;
            eprintln!(
                "{formulation:?} xi {xi:e}: n_s {n_s:.5} (slow roll {n_s_slow_roll:.5}), r {r:.4e} ({r_slow_roll:.4e}), N_J {jordan_e_folds:.3}"
            );
        }
    }
    println!("wrote {}", path.display());
    Ok(())
}
