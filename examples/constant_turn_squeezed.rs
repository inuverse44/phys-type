//! Squeezed bispectra on a constant-turn orbit, for comparison with the
//! nonanalytic squeezed scaling of arXiv:2607.14529.
//!
//! For isosceles triangles `k_1 = k_2 = k_pivot`, `k_3 = r k_pivot`, the
//! leading squeezed shape is `S ~ r^(1/2 - nu)`, so the reduced
//! `f_NL ~ r^(3/2 - nu)` with `nu = sqrt(9/4 - mu^2/H^2)` set by the entropy
//! mass `mu^2 = m_eff^2 + lambda^2 H^2`, not by `m_eff`. Writes
//! `constant_turn_squeezed.csv` into the directory given as the first
//! argument (default `diagnostics/`). The run takes several minutes.

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::{BispectrumRunSettings, TriangleShape};

#[path = "shared/constant_turn.rs"]
mod constant_turn;

use constant_turn::{HUBBLE, radius, revision};

const EPSILON: f64 = 1.0e-5;

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;
    let path = directory.join("constant_turn_squeezed.csv");
    let mut file = fs::File::create(&path)?;
    writeln!(file, "# generator: phys-type {}", env!("CARGO_PKG_VERSION"))?;
    if let Some(revision) = revision() {
        writeln!(file, "# revision: {revision}")?;
    }
    writeln!(file, "# reference: arXiv:2607.14529, Sec. X")?;
    writeln!(
        file,
        "# model: V = V0 - c theta + m^2 (rho - rho0)^2 / 2 on PolarFieldSpace"
    )?;
    writeln!(file, "# epsilon: {EPSILON:e}")?;
    writeln!(file, "# triangle: k1 = k2 = k_pivot, k3 = r k_pivot")?;
    writeln!(file, "# f_nl: tree-level transport, Maldacena/Planck sign")?;
    writeln!(file, "# lambda_and_mu2: measured at the exit of k_pivot")?;
    writeln!(file, "lambda_target,mu2_target,lambda,mu2,r,f_nl")?;

    let ratios = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01];
    let shapes: Vec<TriangleShape> = ratios
        .iter()
        .map(|&r| TriangleShape {
            k_t_over_pivot: 2.0 + r,
            alpha: 0.0,
            beta: (2.0 - r) / (2.0 + r),
        })
        .collect();
    // nu = 1, 1, and 1.2: exponents 1/2, 1/2, and 3/10 apart from the
    // analytic r^1 term, while m_eff^2 < 0 for lambda = 2.
    for (lambda, mu_squared, after_exit) in
        [(1.0, 1.25, 20.0), (2.0, 1.25, 20.0), (2.0, 0.81, 30.0)]
    {
        let (system, initial, settings) = constant_turn::orbit(
            lambda,
            mu_squared,
            radius(lambda, EPSILON),
            after_exit,
            16.0,
        );
        let spectrum = system
            .run_two_field_spectra(initial, &settings)
            .unwrap_or_else(|error| {
                eprintln!("lambda {lambda}, mu^2 {mu_squared}: {error}");
                std::process::exit(1);
            });
        let exit = spectrum
            .background
            .iter()
            .min_by(|a, b| {
                let distance = |row: &&phys_type::BackgroundRecord| {
                    (row.e_folds - spectrum.pivot.e_folds_at_pivot).abs()
                };
                distance(a).total_cmp(&distance(b))
            })
            .expect("the survey records the background");
        let h = spectrum.pivot.hubble.into_value();
        let measured_lambda = 2.0 * exit.eta_perp.abs();
        let measured_mu_squared = constant_turn::radial_mass_squared(lambda, mu_squared)
            * (HUBBLE / h).powi(2)
            + 3.0 * exit.eta_perp * exit.eta_perp;
        let bispectrum = system
            .run_bispectrum(
                initial,
                &BispectrumRunSettings {
                    background: settings,
                    ..BispectrumRunSettings::default()
                },
                &shapes,
            )
            .unwrap_or_else(|error| {
                eprintln!("lambda {lambda}, mu^2 {mu_squared}: {error}");
                std::process::exit(1);
            });
        for (r, record) in ratios.iter().zip(&bispectrum.records) {
            writeln!(
                file,
                "{lambda},{mu_squared},{measured_lambda:.9},{measured_mu_squared:.9},{r},{:.10e}",
                record.f_nl
            )?;
        }
        eprintln!("lambda {lambda}, mu^2 {mu_squared}: done");
    }
    println!("wrote {}", path.display());
    Ok(())
}
