//! A complete two-field adiabatic and entropy spectrum run.
//!
//! Writes three CSV files into the directory given as the first argument
//! (default `diagnostics/`):
//!
//! - `two_field_background.csv` the background trajectory and its
//!   adiabatic/entropy frame
//! - `two_field_spectra.csv`    one row per wavenumber
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
    CovariantMultifieldSystem, Field, FieldVelocity, FlatFieldSpace, Mass, MultifieldState,
    QuadraticMultifieldPotential, ScaleFactor, SpectrumRun, SpectrumRunSettings, write_metadata,
};

/// Unequal masses bend the trajectory in field space, so the entropy
/// direction has something to source.
const MASSES: [f64; 2] = [1.0e-6, 2.0e-6];
const PHI_START: [f64; 2] = [16.0, 12.0];

fn model() -> CovariantMultifieldSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace, 2> {
    CovariantMultifieldSystem::in_planck_units(
        QuadraticMultifieldPotential::new([Mass::new(MASSES[0]), Mass::new(MASSES[1])]),
        FlatFieldSpace,
    )
}

/// The slow-roll initial state at `PHI_START`.
fn initial_state(
    system: &CovariantMultifieldSystem<QuadraticMultifieldPotential<2>, FlatFieldSpace, 2>,
) -> MultifieldState<2> {
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
    revision: Option<&str>,
    rows: impl Fn(&SpectrumRun, &mut fs::File) -> io::Result<()>,
) -> io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "# masses_over_m_pl: {:e},{:e}", MASSES[0], MASSES[1])?;
    writeln!(
        file,
        "# phi_start_over_m_pl: {:e},{:e}",
        PHI_START[0], PHI_START[1]
    )?;
    write_metadata(&mut file, settings, &run.pivot, revision)?;

    rows(run, &mut file)
}

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;

    let system = model();
    // The example reports the final inflationary spectra, after the entropy
    // mode has had the full inflationary interval in which to source the
    // curvature mode.
    let settings = SpectrumRunSettings {
        output_e_folds_before_end: 0.0,
        ..SpectrumRunSettings::default()
    };
    let run = match system.run_two_field_spectra(initial_state(&system), &settings) {
        Ok(run) => run,
        Err(error) => {
            eprintln!("the run could not be completed: {error}");
            std::process::exit(1);
        }
    };
    let revision = revision();
    let revision = revision.as_deref();

    write_file(
        &directory.join("two_field_background.csv"),
        &settings,
        &run,
        revision,
        |run, file| run.write_background_csv(file),
    )?;
    write_file(
        &directory.join("two_field_spectra.csv"),
        &settings,
        &run,
        revision,
        |run, file| run.write_spectra_csv(file),
    )?;
    write_file(
        &directory.join("two_field_modes.csv"),
        &settings,
        &run,
        revision,
        |run, file| run.write_mode_evolution_csv(file),
    )?;

    eprintln!("e-folds of inflation   {:.4}", run.pivot.e_folds_at_end);
    eprintln!("pivot at N             {:.4}", run.pivot.e_folds_at_pivot);
    eprintln!("k_pivot                {:e} M_pl", run.pivot.k.into_value());
    eprintln!(
        "H at pivot             {:e} M_pl",
        run.pivot.hubble.into_value()
    );
    eprintln!("epsilon_H at pivot     {:e}", run.pivot.epsilon_h);
    eprintln!("background rows        {}", run.background.len());
    eprintln!("mode evolution rows    {}", run.mode_evolution.len());
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
