//! The two-field nonminimal model of arXiv:1210.7487 (KMS) off its ridge.
//!
//! Evolves trajectories 1--3 of KMS Fig. 2 from rest at `phi = 3.1` and
//! writes, into the directory given as the first argument (default
//! `diagnostics/`),
//!
//! - `kms_trajectories.csv`: the background of each trajectory at every
//!   accepted step, with `|eta_perp|` (the counterpart of Figs. 2--3 and 7);
//! - `kms_f_nl.csv`: the crossing-time delta N `f_NL` of KMS Sec. V against
//!   `chi(tau_0)`, the counterpart of Fig. 11, with the pivot 60 e-folds
//!   before the end. Where the fields leave the ridge close to the end,
//!   `epsilon_H` reaches one during the turn and the survey ends there, at
//!   the first crossing; `e_folds_total` marks those points.

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::{Formulation, SpectrumRunSettings};
use phys_type_models::kms::{PHI_START, TRAJECTORIES, at_rest, benchmark_system, crossing_delta_n};

#[path = "shared/provenance.rs"]
mod provenance;

use provenance::revision;

const PIVOT_BEFORE_END: f64 = 60.0;

fn header(file: &mut fs::File) -> io::Result<()> {
    writeln!(
        file,
        "# generator: phys-type-models {}",
        env!("CARGO_PKG_VERSION")
    )?;
    if let Some(revision) = revision() {
        writeln!(file, "# revision: {revision}")?;
    }
    writeln!(file, "# reference: arXiv:1210.7487")?;
    writeln!(
        file,
        "# model: Omega^2 = 1 + 10 phi^2 + 10.02 chi^2, V_J = lambda (phi^4 + 2 phi^2 chi^2 + chi^4 / 2) / 4, metric"
    )?;
    writeln!(
        file,
        "# start: phi = {PHI_START}, chi = chi0, at rest; end at epsilon_H = 1"
    )
}

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;
    let system = benchmark_system(Formulation::Metric);
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: PIVOT_BEFORE_END,
        e_folds_before_exit: 6.0,
        background_sample_every: 1,
        ..SpectrumRunSettings::default()
    };

    let mut file = fs::File::create(directory.join("kms_trajectories.csv"))?;
    header(&mut file)?;
    writeln!(
        file,
        "trajectory,chi0,e_folds,e_folds_before_end,phi,chi,epsilon_h,eta_perp"
    )?;
    for (index, chi0) in TRAJECTORIES.into_iter().enumerate() {
        let survey = system
            .survey(at_rest(PHI_START, chi0), &settings)
            .map_err(|error| io::Error::other(format!("trajectory {}: {error:?}", index + 1)))?;
        for sample in &survey.samples {
            writeln!(
                file,
                "{},{chi0:e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                index + 1,
                sample.e_folds,
                survey.e_folds_at_end - sample.e_folds,
                sample.phi[0],
                sample.phi[1],
                sample.epsilon_h,
                sample.eta_perp,
            )?;
        }
        eprintln!(
            "trajectory {}: {:.3} e-folds",
            index + 1,
            survey.e_folds_at_end
        );
    }

    let mut file = fs::File::create(directory.join("kms_f_nl.csv"))?;
    header(&mut file)?;
    writeln!(file, "# pivot_e_folds_before_end: {PIVOT_BEFORE_END}")?;
    writeln!(file, "chi0,e_folds_total,f_nl,power")?;
    // Twelve points per decade, and the narrow window where f_NL is large
    // resolved with 81 more.
    let coarse = (0..=36).map(|i| 1.0e-5 * 10f64.powf(f64::from(i) / 12.0));
    let window = (0..=80).map(|i| 4.0e-4 * (2.5e-3 / 4.0e-4_f64).powf(f64::from(i) / 80.0));
    let mut scan: Vec<f64> = coarse.chain(window).chain(TRAJECTORIES).collect();
    scan.sort_by(f64::total_cmp);
    for chi0 in scan {
        let result = system
            .survey(at_rest(PHI_START, chi0), &settings)
            .and_then(|survey| {
                crossing_delta_n(&system, &survey, &settings, PIVOT_BEFORE_END)
                    .map(|delta_n| (survey.e_folds_at_end, delta_n))
            });
        match result {
            Ok((total, delta_n)) => writeln!(
                file,
                "{chi0:.17e},{total:.17e},{:.17e},{:.17e}",
                delta_n.f_nl, delta_n.power
            )?,
            Err(error) => eprintln!("chi0 {chi0:e}: {error:?}"),
        }
    }
    Ok(())
}
