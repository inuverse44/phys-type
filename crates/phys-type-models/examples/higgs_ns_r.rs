//! Higgs inflation on a fine grid of `xi` and `N_*`, for the `n_s`-`r` plane.
//!
//! The same exact evolution as `higgs_inflation`, for `log10 xi` from -5 to 4
//! in steps of 1/4 and pivots 50 to 60 Einstein-frame e-folds before the end
//! in steps of 2, in the metric and Palatini formulations; the layout follows
//! Fig. 4 of Kodama and Takahashi, arXiv:2112.05283. Writes `higgs_ns_r.csv`
//! into the directory given as the first argument (default `diagnostics/`).

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::Formulation;

#[path = "shared/higgs.rs"]
mod higgs;
#[path = "shared/provenance.rs"]
mod provenance;

use higgs::{OUTPUT_BEFORE_END, predict};
use provenance::revision;

const LOG_XI_RANGE: (i32, i32) = (-20, 16);
const LOG_XI_STEPS_PER_DECADE: f64 = 4.0;
const PIVOTS_BEFORE_END: [f64; 6] = [50.0, 52.0, 54.0, 56.0, 58.0, 60.0];

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;
    let path = directory.join("higgs_ns_r.csv");
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
        "# reference: arXiv:0710.3755 (metric), arXiv:0803.2664 (Palatini); layout of arXiv:2112.05283 Fig. 4"
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
        "pivot_e_folds,formulation,log10_xi,n_s,r,n_s_slow_roll,r_slow_roll,jordan_e_folds,phi_star"
    )?;
    for formulation in [Formulation::Metric, Formulation::Palatini] {
        for step in LOG_XI_RANGE.0..=LOG_XI_RANGE.1 {
            let log10_xi = f64::from(step) / LOG_XI_STEPS_PER_DECADE;
            for pivot_before_end in PIVOTS_BEFORE_END {
                match predict(10f64.powf(log10_xi), formulation, pivot_before_end) {
                    Ok(p) => writeln!(
                        file,
                        "{pivot_before_end},{formulation:?},{log10_xi},{:.12},{:.12e},{:.12},{:.12e},{:.10},{:.12e}",
                        p.n_s, p.r, p.n_s_slow_roll, p.r_slow_roll, p.jordan_e_folds, p.phi_star
                    )?,
                    Err(error) => {
                        eprintln!(
                            "{formulation:?} log10 xi {log10_xi} N {pivot_before_end}: {error}"
                        );
                    }
                }
            }
        }
    }
    println!("wrote {}", path.display());
    Ok(())
}
