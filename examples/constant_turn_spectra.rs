//! Superhorizon power spectra on a constant-turn orbit, for comparison with
//! the exact solutions of arXiv:2606.18248.
//!
//! The background is the circular orbit of
//! `V = V_0 - c theta + m^2 (rho - rho_0)^2 / 2` on `PolarFieldSpace`, with
//! turn rate `Omega`, mixing `lambda = 2 |Omega| / H`, and entropy mass
//! `mu^2 = m^2 + 3 Omega^2`. Writes two CSV files into the directory given as
//! the first argument (default `diagnostics/`):
//!
//! - `constant_turn_massive.csv`  `P_zeta / Delta_0` for massive entropy
//!   fields, `Delta_0 = H^2 / (8 pi^2 epsilon)` at horizon exit;
//! - `constant_turn_massless.csv` the same for a massless entropy field at
//!   `z_end = k / (aH) = e^-10`.
//!
//! `lambda` and `mu2` are measured at horizon exit; the closed forms are
//! evaluated at those values by `plots/constant_turn.py`.

use std::fs;
use std::io::{self, Write};
use std::ops::{Add, Mul};
use std::path::Path;
use std::process::Command;

use phys_type::{
    CovariantMultifieldSystem, EndOfInflation, EnergyDensity, Field, FieldVelocity,
    MultifieldPotentialValue, MultifieldState, PolarFieldSpace, Quantity, ScaleFactor,
    SpectrumRunSettings,
};

const HUBBLE: f64 = 1.0e-5;
const BEFORE_EXIT: f64 = 6.0;

#[derive(Clone, Copy, Debug)]
struct ConstantTurnPotential {
    v0: f64,
    slope: f64,
    mass_squared: f64,
    rho0: f64,
}

impl<V> MultifieldPotentialValue<V, 2> for ConstantTurnPotential
where
    V: From<f64> + Clone + Add<V, Output = V> + Mul<V, Output = V> + Mul<f64, Output = V>,
{
    fn value(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        let [rho, theta] = fields.map(Quantity::into_value);
        let radial = rho + V::from(-self.rho0);
        EnergyDensity::new(
            V::from(self.v0)
                + theta * (-self.slope)
                + radial.clone() * radial * (0.5 * self.mass_squared),
        )
    }
}

/// One orbit: measured `lambda`, `mu^2 / H^2`, `P_zeta / Delta_0`, and `z`
/// at the output.
#[derive(Clone, Copy, Debug)]
struct Orbit {
    lambda: f64,
    mu_squared: f64,
    ratio: f64,
    z_end: f64,
}

fn orbit(lambda: f64, mu_squared: f64, radius: f64, after_exit: f64) -> Orbit {
    let omega = 0.5 * lambda * HUBBLE;
    let mass_squared = (mu_squared - 0.75 * lambda * lambda) * HUBBLE * HUBBLE;
    // Radial balance m^2 (rho - rho_0) = rho Omega^2 and the angular
    // slow-roll attractor c = 3 H rho^2 Omega.
    let potential = ConstantTurnPotential {
        v0: 3.0 * HUBBLE * HUBBLE,
        slope: 3.0 * HUBBLE * radius * radius * omega,
        mass_squared,
        rho0: radius - radius * omega * omega / mass_squared,
    };
    let system = CovariantMultifieldSystem::in_planck_units(potential, PolarFieldSpace);
    let pivot_before_end = after_exit + 1.0;
    let total = pivot_before_end + BEFORE_EXIT + 8.0;
    let settings = SpectrumRunSettings {
        end_of_inflation: EndOfInflation::FieldValue {
            field: 1,
            value: Field::new(omega / HUBBLE * total),
        },
        e_folds_before_exit: BEFORE_EXIT,
        pivot_e_folds_before_end: pivot_before_end,
        output_e_folds_before_end: 1.0,
        k_min_over_pivot: 1.0,
        k_max_over_pivot: 1.0,
        k_count: 1,
        background_sample_every: 1,
        ..SpectrumRunSettings::default()
    };
    let initial = MultifieldState::new(
        [Field::new(radius), Field::new(0.0)],
        [FieldVelocity::new(0.0), FieldVelocity::new(omega)],
        ScaleFactor::new(1.0),
    );
    let run = system
        .run_two_field_spectra(initial, &settings)
        .unwrap_or_else(|error| {
            eprintln!("lambda {lambda}, mu^2 {mu_squared}: {error}");
            std::process::exit(1);
        });
    let exit = run
        .background
        .iter()
        .min_by(|a, b| {
            let distance = |row: &&phys_type::BackgroundRecord| {
                (row.e_folds - run.pivot.e_folds_at_pivot).abs()
            };
            distance(a).total_cmp(&distance(b))
        })
        .expect("the survey records the background");
    let h = run.pivot.hubble.into_value();
    let delta_0 = h * h / (8.0 * std::f64::consts::PI.powi(2) * run.pivot.epsilon_h);
    let spectrum = run.spectra[0];
    Orbit {
        lambda: 2.0 * exit.eta_perp.abs(),
        mu_squared: mass_squared / (h * h) + 3.0 * exit.eta_perp * exit.eta_perp,
        ratio: spectrum.p_r / delta_0,
        z_end: spectrum.k_over_ah,
    }
}

/// An orbit whose exit values of `lambda` and `mu^2` match the targets.
fn calibrated_orbit(lambda: f64, mu_squared: f64, radius: f64, after_exit: f64) -> Orbit {
    let (mut lambda_in, mut mu_squared_in) = (lambda, mu_squared);
    let mut result = orbit(lambda_in, mu_squared_in, radius, after_exit);
    for _ in 0..3 {
        lambda_in *= lambda / result.lambda;
        mu_squared_in += mu_squared - result.mu_squared;
        result = orbit(lambda_in, mu_squared_in, radius, after_exit);
    }
    result
}

/// `epsilon = rho^2 lambda^2 / 8` fixes the orbit radius.
fn radius(lambda: f64, epsilon: f64) -> f64 {
    (8.0 * epsilon).sqrt() / lambda
}

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

fn header(file: &mut fs::File, entropy_field: &str) -> io::Result<()> {
    writeln!(file, "# generator: phys-type {}", env!("CARGO_PKG_VERSION"))?;
    if let Some(revision) = revision() {
        writeln!(file, "# revision: {revision}")?;
    }
    writeln!(file, "# reference: arXiv:2606.18248")?;
    writeln!(
        file,
        "# model: V = V0 - c theta + m^2 (rho - rho0)^2 / 2 on PolarFieldSpace"
    )?;
    writeln!(file, "# entropy_field: {entropy_field}")?;
    writeln!(file, "# hubble_over_m_pl: {HUBBLE:e}")?;
    writeln!(file, "# e_folds_before_exit: {BEFORE_EXIT}")?;
    writeln!(
        file,
        "# normalization: ratio = P_zeta / (H^2 / (8 pi^2 epsilon)) at exit"
    )?;
    writeln!(file, "# lambda_and_mu2: measured at horizon exit")
}

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;

    let epsilon = 1.0e-5;
    let mut massive = fs::File::create(directory.join("constant_turn_massive.csv"))?;
    header(&mut massive, "massive")?;
    writeln!(massive, "# epsilon: {epsilon:e}")?;
    writeln!(
        massive,
        "lambda_target,mu2_target,lambda,mu2,after_exit,ratio"
    )?;
    // mu^2 = 3 lambda^2 / 4 would leave no radial restoring force (m = 0).
    let masses = [
        1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.25, 4.0, 5.0, 6.5, 9.0, 12.0, 16.0,
    ];
    for lambda in [0.5_f64, 1.0, 2.0] {
        for mu_squared in masses {
            let mu_squared: f64 = mu_squared;
            // A light sigma decays as z^(3/2 - nu); wait until that tail is
            // negligible.
            let decay = 1.5 - (2.25 - mu_squared).max(0.0).sqrt();
            let after_exit = (14.0 / decay).max(12.0);
            let run = orbit(lambda, mu_squared, radius(lambda, epsilon), after_exit);
            writeln!(
                massive,
                "{lambda},{mu_squared},{:.9},{:.9},{after_exit},{:.10e}",
                run.lambda, run.mu_squared, run.ratio
            )?;
        }
    }

    let mut massless = fs::File::create(directory.join("constant_turn_massless.csv"))?;
    header(&mut massless, "massless")?;
    writeln!(
        massless,
        "# orbit: corrected until the exit lambda and mu2 match the targets"
    )?;
    writeln!(massless, "lambda_target,lambda,mu2,epsilon,z_end,ratio")?;
    for step in 1..=12 {
        let lambda = 0.2 * f64::from(step);
        // The slow-roll correction accumulates during the superhorizon
        // growth, so the massless orbit needs a small epsilon.
        let epsilon = 3.0e-7;
        let run = calibrated_orbit(lambda, 0.0, radius(lambda, epsilon), 10.0);
        writeln!(
            massless,
            "{lambda:.1},{:.9},{:.3e},{epsilon:e},{:.9e},{:.10e}",
            run.lambda, run.mu_squared, run.z_end, run.ratio
        )?;
    }

    println!(
        "wrote {} and {}",
        directory.join("constant_turn_massive.csv").display(),
        directory.join("constant_turn_massless.csv").display()
    );
    Ok(())
}
