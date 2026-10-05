//! Higgs inflation (arXiv:0710.3755), metric and Palatini, evolved exactly.
//!
//! For a range of nonminimal couplings `xi`, the theory
//! `Omega^2 = 1 + xi phi^2`, `V_J = lambda phi^4 / 4` is mapped to the
//! Einstein frame by `phys_type::einstein_frame` and evolved with one field.
//! The pivot leaves the horizon 50 or 60 Einstein-frame e-folds before
//! `epsilon_H = 1`, the range usually quoted for the reheating uncertainty;
//! `n_s` comes from the mode solver and `r = P_t / P_R` from the tensor and
//! curvature spectra one e-fold before the end. Writes `higgs_inflation.csv`
//! into the directory given as the first argument (default `diagnostics/`),
//! next to the first-order slow-roll values from the same Einstein-frame
//! potential and the number of Jordan-frame e-folds over the same interval.

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::Formulation;

#[path = "shared/higgs.rs"]
mod higgs;
#[path = "shared/provenance.rs"]
mod provenance;

use higgs::{OUTPUT_BEFORE_END, Prediction, predict};
use provenance::revision;

const PIVOTS_BEFORE_END: [f64; 2] = [50.0, 60.0];

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
        "# pivot_e_folds_before_end: column pivot_e_folds (Einstein frame, to epsilon_H = 1)"
    )?;
    writeln!(file, "# output_e_folds_before_end: {OUTPUT_BEFORE_END}")?;
    writeln!(
        file,
        "pivot_e_folds,formulation,xi,n_s,r,n_s_slow_roll,r_slow_roll,jordan_e_folds,phi_star"
    )?;

    let couplings = [1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4];
    for pivot_before_end in PIVOTS_BEFORE_END {
        for formulation in [Formulation::Metric, Formulation::Palatini] {
            for xi in couplings {
                let Prediction {
                    n_s,
                    r,
                    n_s_slow_roll,
                    r_slow_roll,
                    jordan_e_folds,
                    phi_star,
                } = match predict(xi, formulation, pivot_before_end) {
                    Ok(prediction) => prediction,
                    Err(error) => {
                        eprintln!("{formulation:?} xi {xi}: {error}");
                        continue;
                    }
                };
                writeln!(
                    file,
                    "{pivot_before_end},{formulation:?},{xi:e},{n_s:.12},{r:.12e},{n_s_slow_roll:.12},{r_slow_roll:.12e},{jordan_e_folds:.10},{phi_star:.12e}"
                )?;
                eprintln!(
                    "N {pivot_before_end} {formulation:?} xi {xi:e}: n_s {n_s:.5} (slow roll {n_s_slow_roll:.5}), r {r:.4e} ({r_slow_roll:.4e}), N_J {jordan_e_folds:.3}"
                );
            }
        }
    }
    println!("wrote {}", path.display());
    Ok(())
}
