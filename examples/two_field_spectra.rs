//! A complete two-field adiabatic and entropy spectrum run.
//!
//! Writes four CSV files into the directory given as the first argument
//! (default `diagnostics/`):
//!
//! - `two_field_background.csv` the background trajectory and its
//!   adiabatic/entropy frame
//! - `two_field_spectra.csv`    one row per wavenumber
//! - `two_field_horizon_exit_spectra.csv` one row per wavenumber at `k = aH`
//! - `two_field_modes.csv`      sampled evolution of every mode
//!
//! Each file opens with a metadata comment block recording the model, the
//! tolerances, `x_init`, the pivot convention, and the Git revision when one
//! is available, so a row can be traced back to the run that produced it.

use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::process::Command;

use phys_type::{
    AmplitudeRole, CovariantMultifieldSystem, Field, FieldVelocity, FlatFieldSpace, Mass,
    MultifieldState, PLANCK_2018_TT_TE_EE_LOW_E, QuadraticMultifieldPotential, ScaleFactor,
    SpectrumRun, SpectrumRunSettings, write_metadata, write_observational_metadata,
};

/// Shape parameters for the observationally normalized double-quadratic
/// benchmark. The hierarchy makes the pivot mode cross before, and then
/// experience, the heavy-to-light transition for these initial fields.
const MASS_RATIO: f64 = 7.0;
const PHI_START: [f64; 2] = [14.0, 12.0];
const REFERENCE_LIGHT_MASS: f64 = 1.0e-6;
const TARGET_SCALAR_AMPLITUDE: f64 = 2.1e-9;

type BenchmarkSystem =
    CovariantMultifieldSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace, 2>;

fn model(light_mass: f64) -> BenchmarkSystem {
    CovariantMultifieldSystem::in_planck_units(
        QuadraticMultifieldPotential::new([
            Mass::new(light_mass),
            Mass::new(MASS_RATIO * light_mass),
        ]),
        FlatFieldSpace,
    )
}

/// The slow-roll initial state at `PHI_START`.
fn initial_state(system: &BenchmarkSystem) -> MultifieldState<2> {
    let resting = MultifieldState::new(
        [Field::new(PHI_START[0]), Field::new(PHI_START[1])],
        [FieldVelocity::new(0.0); 2],
        ScaleFactor::new(1.0),
    );

    MultifieldState {
        pi: system.slow_roll_velocity(&resting),
        ..resting
    }
}

/// The current Git revision, when the example is run inside a checkout.
fn revision() -> Option<String> {
    let output = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()?;

    let revision = output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())?;
    let status = Command::new("git")
        .args(["status", "--porcelain", "--untracked-files=no"])
        .output()
        .ok()?;
    let dirty = status.status.success() && !status.stdout.is_empty();

    Some(if dirty {
        format!("{revision}-dirty")
    } else {
        revision
    })
}

fn write_file(
    path: &Path,
    settings: &SpectrumRunSettings,
    run: &SpectrumRun,
    calibration: &Calibration,
    revision: Option<&str>,
    epoch: &str,
    rows: impl Fn(&SpectrumRun, &mut fs::File) -> io::Result<()>,
) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "# benchmark: canonical_double_quadratic")?;
    writeln!(file, "# mass_ratio_m2_over_m1: {:e}", MASS_RATIO)?;
    writeln!(
        file,
        "# masses_over_m_pl: {:e},{:e}",
        calibration.light_mass,
        MASS_RATIO * calibration.light_mass
    )?;
    writeln!(
        file,
        "# phi_start_over_m_pl: {:e},{:e}",
        PHI_START[0], PHI_START[1]
    )?;
    writeln!(
        file,
        "# target_scalar_amplitude: {:e}",
        TARGET_SCALAR_AMPLITUDE
    )?;
    writeln!(file, "# pilot_p_r_at_pivot: {:e}", calibration.pilot_p_r)?;
    writeln!(
        file,
        "# common_mass_rescaling: {:e}",
        calibration.mass_rescaling
    )?;
    writeln!(file, "# evaluation_epoch: {epoch}")?;
    write_observational_metadata(
        &mut file,
        &PLANCK_2018_TT_TE_EE_LOW_E,
        AmplitudeRole::BenchmarkTargetNotLikelihoodFit,
    )?;
    write_metadata(&mut file, settings, &run.pivot, revision)?;

    rows(run, &mut file)
}

#[derive(Clone, Copy, Debug)]
struct Calibration {
    pilot_p_r: f64,
    mass_rescaling: f64,
    light_mass: f64,
}

fn run_or_exit(
    system: &BenchmarkSystem,
    settings: &SpectrumRunSettings,
    stage: &str,
) -> SpectrumRun {
    match system.run_two_field_spectra(initial_state(system), settings) {
        Ok(run) => run,
        Err(error) => {
            eprintln!("the {stage} run could not be completed: {error}");
            std::process::exit(1);
        }
    }
}

fn calibrate(settings: &SpectrumRunSettings) -> Calibration {
    let pilot_settings = SpectrumRunSettings {
        k_min_over_pivot: 1.0,
        k_max_over_pivot: 1.0,
        k_count: 1,
        background_sample_every: settings.background_sample_every,
        mode_sample_every: settings.mode_sample_every,
        ..*settings
    };
    let system = model(REFERENCE_LIGHT_MASS);
    let pilot = run_or_exit(&system, &pilot_settings, "amplitude-calibration pilot");
    let pilot_p_r = pilot.spectra[0].p_r;
    if !pilot_p_r.is_finite() || pilot_p_r <= 0.0 {
        eprintln!("the pilot returned an invalid pivot spectrum: {pilot_p_r:e}");
        std::process::exit(1);
    }
    let mass_rescaling = (TARGET_SCALAR_AMPLITUDE / pilot_p_r).sqrt();

    Calibration {
        pilot_p_r,
        mass_rescaling,
        light_mass: REFERENCE_LIGHT_MASS * mass_rescaling,
    }
}

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;

    // The example reports the final inflationary spectra, after the entropy
    // mode has had the full inflationary interval in which to source the
    // curvature mode.
    let settings = SpectrumRunSettings {
        output_e_folds_before_end: 0.0,
        k_min_over_pivot: 0.01,
        k_max_over_pivot: 100.0,
        k_count: 41,
        // The turn is sharp in e-fold time. Keep enough accepted states to
        // draw its entropy trough and curvature transfer without angular
        // line segments; this changes CSV sampling, not integration accuracy.
        mode_sample_every: 5,
        ..SpectrumRunSettings::default()
    };
    let calibration = calibrate(&settings);
    let system = model(calibration.light_mass);
    let run = run_or_exit(&system, &settings, "calibrated spectrum");
    let revision = revision();
    let revision = revision.as_deref();

    write_file(
        &directory.join("two_field_background.csv"),
        &settings,
        &run,
        &calibration,
        revision,
        "background_survey",
        |run, file| run.write_background_csv(file),
    )?;
    write_file(
        &directory.join("two_field_spectra.csv"),
        &settings,
        &run,
        &calibration,
        revision,
        "requested_output",
        |run, file| run.write_spectra_csv(file),
    )?;
    write_file(
        &directory.join("two_field_horizon_exit_spectra.csv"),
        &settings,
        &run,
        &calibration,
        revision,
        "horizon_exit_k_equals_aH",
        |run, file| run.write_horizon_exit_spectra_csv(file),
    )?;
    write_file(
        &directory.join("two_field_modes.csv"),
        &settings,
        &run,
        &calibration,
        revision,
        "mode_evolution",
        |run, file| run.write_mode_evolution_csv(file),
    )?;

    eprintln!("e-folds of inflation   {:.4}", run.pivot.e_folds_at_end);
    eprintln!("pivot at N             {:.4}", run.pivot.e_folds_at_pivot);
    eprintln!("k_pivot (internal)     {:e} M_pl", run.pivot.k.into_value());
    eprintln!(
        "H at pivot             {:e} M_pl",
        run.pivot.hubble.into_value()
    );
    eprintln!("epsilon_H at pivot     {:e}", run.pivot.epsilon_h);
    eprintln!("target P_R(k_pivot)    {:e}", TARGET_SCALAR_AMPLITUDE);
    eprintln!("pilot P_R(k_pivot)     {:e}", calibration.pilot_p_r);
    eprintln!("common mass rescaling {:.8}", calibration.mass_rescaling);
    eprintln!(
        "calibrated masses      [{:e}, {:e}] M_pl",
        calibration.light_mass,
        MASS_RATIO * calibration.light_mass
    );
    eprintln!("background rows        {}", run.background.len());
    eprintln!("mode evolution rows    {}", run.mode_evolution.len());
    eprintln!("horizon-exit rows     {}", run.horizon_exit_spectra.len());
    eprintln!();
    eprintln!("k/k_pivot        P_R           P_S           r_C");
    for record in &run.spectra {
        eprintln!(
            "{:<12.4} {:>13.6e} {:>13.6e} {:>+9.4}",
            record.k_over_pivot, record.p_r, record.p_s, record.r_c
        );
    }

    Ok(())
}
