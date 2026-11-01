//! Reproduction data for the nonminimally assisted inflation examples of
//! arXiv:2302.05866 (Secs. 4.1--4.3): loop (Class I), power-law, and hybrid
//! (Class II) inflatons with an assistant field `s`.
//!
//! This is one subclass of nonminimal inflation: a single assistant with
//! `f(s) = (M_pl^2/2)(1 + xi_2 s^2)` and a minimally coupled inflaton, in the
//! metric formulation. For each branch `(model, s_*)` and a grid of `xi_2`
//! around the coupling that gives `n_s = 0.965`, the example writes into the
//! directory given as the first argument (default `diagnostics/`):
//!
//! - `assisted_predictions.csv` the exact `n_s` (mode solver) and `f_NL`
//!   (tree-level transport, equilateral) 1e-3 e-folds before the end, next to
//!   Eqs. (3.1) and (3.5) of the reference evaluated with the numerical
//!   `phi_*`, `s_e`, and end surface, and with the reference's own
//!   slow-roll prescription; `f_NL` is in the Maldacena/Planck sign;
//! - `assisted_trajectories.csv` the background `(phi, s)` at the coupling
//!   closest to `n_s = 0.965`.
//!
//! The full grid takes several minutes in release mode.

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::{
    BackgroundRecord, BispectrumRunSettings, CovariantMultifieldSystem, EinsteinGeometry,
    EinsteinPotential, Field, FieldVelocity, MultifieldState, ScaleFactor, SpectrumRunSettings,
    TriangleShape, einstein_frame,
};
use phys_type_models::assisted::{
    AssistedInflaton, AssistedTheory, HybridInflaton, LoopInflaton, PowerLawInflaton, bisect,
};
use phys_type_models::hkkpt::{
    assisted_delta_n_prediction, assisted_slow_roll, paper_prescription,
};

#[path = "shared/provenance.rs"]
mod provenance;

use provenance::revision;

type AssistedSystem<I> = CovariantMultifieldSystem<
    EinsteinPotential<AssistedTheory<I>>,
    EinsteinGeometry<AssistedTheory<I>>,
    2,
>;

const PIVOT_BEFORE_END: f64 = 60.0;
const OUTPUT_BEFORE_END: f64 = 1.0e-3;
const START_E_FOLDS: f64 = 72.0;
const TILT_STEP: f64 = 0.02;

/// One numerical point of a branch.
struct Point {
    xi: f64,
    n_s: f64,
    f_nl: f64,
    phi_star: f64,
    s_star: f64,
    phi_end: f64,
    s_end: f64,
    background: Vec<BackgroundRecord>,
}

fn start<I: AssistedInflaton>(system: &AssistedSystem<I>, phi: f64, s: f64) -> MultifieldState<2> {
    let resting = MultifieldState::new(
        [Field::new(phi), Field::new(s)],
        [FieldVelocity::new(0.0); 2],
        ScaleFactor::new(1.0),
    );
    MultifieldState {
        pi: system.slow_roll_velocity(&resting),
        ..resting
    }
}

/// Linear interpolation of the background fields at `e_folds`.
fn fields_at(rows: &[BackgroundRecord], e_folds: f64) -> [f64; 2] {
    let index = rows
        .partition_point(|row| row.e_folds < e_folds)
        .clamp(1, rows.len() - 1);
    let (a, b) = (&rows[index - 1], &rows[index]);
    let t = (e_folds - a.e_folds) / (b.e_folds - a.e_folds);
    [0, 1].map(|i| a.phi[i] + t * (b.phi[i] - a.phi[i]))
}

/// The branch at coupling `xi`, with the initial `s` tuned so that `s = s_*`
/// at the pivot, `PIVOT_BEFORE_END` e-folds before the end.
fn point<I: AssistedInflaton>(inflaton: I, s_star: f64, xi: f64) -> Result<Point, String> {
    let (potential, geometry) = einstein_frame(AssistedTheory::metric(inflaton, xi));
    let system = CovariantMultifieldSystem::in_planck_units(potential, geometry);
    let settings = SpectrumRunSettings {
        pivot_e_folds_before_end: PIVOT_BEFORE_END,
        output_e_folds_before_end: OUTPUT_BEFORE_END,
        k_min_over_pivot: (-TILT_STEP).exp(),
        k_max_over_pivot: TILT_STEP.exp(),
        k_count: 2,
        background_sample_every: 1,
        end_of_inflation: inflaton.end_of_inflation(),
        ..SpectrumRunSettings::default()
    };
    let (low, high) = inflaton.phi_bracket();
    let phi0 = bisect(
        &|phi| inflaton.slow_roll_e_folds(phi) - START_E_FOLDS,
        low,
        high,
    );
    let mut s0 = 0.35 * s_star;
    let mut run = None;
    for _ in 0..7 {
        let result = system
            .run_two_field_spectra(start::<I>(&system, phi0, s0), &settings)
            .map_err(|e| e.to_string())?;
        let pivot = fields_at(&result.background, result.pivot.e_folds_at_pivot);
        s0 *= s_star / pivot[1];
        run = Some(result);
    }
    let run = run.expect("the calibration runs at least once");
    let star = fields_at(&run.background, run.pivot.e_folds_at_pivot);
    let end = fields_at(&run.background, run.pivot.e_folds_at_end);
    let [low_k, high_k] = [run.spectra[0], run.spectra[1]];
    let n_s = 1.0 + (high_k.p_r / low_k.p_r).ln() / (2.0 * TILT_STEP);

    let bispectrum = system
        .run_bispectrum(
            start::<I>(&system, phi0, s0),
            &BispectrumRunSettings {
                background: SpectrumRunSettings {
                    k_min_over_pivot: 1.0,
                    k_max_over_pivot: 1.0,
                    k_count: 1,
                    ..settings
                },
                ..BispectrumRunSettings::default()
            },
            &[TriangleShape::equilateral(3.0)],
        )
        .map_err(|e| e.to_string())?;

    Ok(Point {
        xi,
        n_s,
        f_nl: bispectrum.records[0].f_nl,
        phi_star: star[0],
        s_star: star[1],
        phi_end: end[0],
        s_end: end[1],
        background: run.background,
    })
}

fn header(file: &mut fs::File) -> io::Result<()> {
    writeln!(
        file,
        "# generator: phys-type-models {}",
        env!("CARGO_PKG_VERSION")
    )?;
    if let Some(revision) = revision() {
        writeln!(file, "# revision: {revision}")?;
    }
    writeln!(file, "# reference: arXiv:2302.05866v2, Secs. 4.1-4.3")?;
    writeln!(
        file,
        "# subclass: nonminimally assisted inflation, f(s) = (1 + xi_2 s^2)/2, kappa = 1"
    )?;
    writeln!(file, "# pivot_e_folds_before_end: {PIVOT_BEFORE_END}")?;
    writeln!(file, "# output_e_folds_before_end: {OUTPUT_BEFORE_END:e}")?;
    writeln!(
        file,
        "# f_nl_sign: Maldacena/Planck (the reference uses the opposite)"
    )
}

/// Writes every point of one branch and returns the one closest to
/// `n_s = 0.965`.
fn branch<I: AssistedInflaton>(
    predictions: &mut fs::File,
    label: &str,
    inflaton: I,
    s_star: f64,
    couplings: &[f64],
) -> io::Result<Option<Point>> {
    let mut best: Option<Point> = None;
    for &xi in couplings {
        let point = match point(inflaton, s_star, xi) {
            Ok(point) => point,
            Err(error) => {
                eprintln!("{label} s* {s_star} xi {xi}: skipped ({error})");
                continue;
            }
        };
        let phi_e = if inflaton.class_one() {
            inflaton.paper_phi_e()
        } else {
            point.phi_end
        };
        let (eq_n_s, eq_f_nl) = assisted_delta_n_prediction(
            &assisted_slow_roll(&inflaton, xi, point.phi_star, point.s_star),
            &assisted_slow_roll(&inflaton, xi, phi_e, point.s_end),
            inflaton.class_one(),
        );
        let (paper_n_s, paper_f_nl, paper_s_e) = paper_prescription(&inflaton, xi, s_star);
        writeln!(
            predictions,
            "{label},{s_star},{xi},{:.8},{:.8},{:.8},{:.8e},{:.8e},{:.8e},{:.8},{:.8},{:.8}",
            point.n_s,
            eq_n_s,
            paper_n_s,
            point.f_nl,
            eq_f_nl,
            paper_f_nl,
            point.phi_star,
            point.s_end,
            paper_s_e
        )?;
        eprintln!(
            "{label} s* {s_star} xi {xi}: n_s {:.4} f_NL {:.4}",
            point.n_s, point.f_nl
        );
        if best
            .as_ref()
            .is_none_or(|b| (point.n_s - 0.965).abs() < (b.n_s - 0.965).abs())
        {
            best = Some(point);
        }
    }
    Ok(best)
}

/// Writes one branch on couplings around `center`, the coupling that gives
/// `n_s = 0.965`, and the trajectory of the point closest to it.
fn record<I: AssistedInflaton>(
    predictions: &mut fs::File,
    trajectories: &mut fs::File,
    label: &str,
    inflaton: I,
    s_star: f64,
    center: f64,
) -> io::Result<()> {
    let couplings = [0.8, 0.9, 0.95, 1.0, 1.05, 1.1].map(|factor| factor * center);
    if let Some(best) = branch(predictions, label, inflaton, s_star, &couplings)? {
        for row in &best.background {
            writeln!(
                trajectories,
                "{label},{s_star},{},{:.6},{:.6},{:.8e},{:.8e}",
                best.xi, best.n_s, row.e_folds, row.phi[0], row.phi[1]
            )?;
        }
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let directory = std::env::args().nth(1).unwrap_or("diagnostics".to_string());
    let directory = Path::new(&directory);
    fs::create_dir_all(directory)?;

    let mut predictions = fs::File::create(directory.join("assisted_predictions.csv"))?;
    header(&mut predictions)?;
    writeln!(
        predictions,
        "model,s_star,xi,n_s,n_s_eq_3_1,n_s_paper,f_nl,f_nl_eq_3_5,f_nl_paper,phi_star,s_e,s_e_paper"
    )?;
    let mut trajectories = fs::File::create(directory.join("assisted_trajectories.csv"))?;
    header(&mut trajectories)?;
    writeln!(trajectories, "model,s_star,xi,n_s,N,phi,s")?;

    for (alpha, s_star, center) in [
        (0.1, 1.0e-3, 0.02282),
        (0.1, 0.1, 0.01367),
        (1.0, 1.0e-3, 0.02138),
        (1.0, 0.1, 0.01235),
    ] {
        let label = format!("loop alpha={alpha}");
        let inflaton = LoopInflaton { alpha };
        record(
            &mut predictions,
            &mut trajectories,
            &label,
            inflaton,
            s_star,
            center,
        )?;
    }
    for (lambda, s_star, center) in [
        (0.01, 1.0e-3, 0.0110),
        (0.01, 0.1, 0.00516),
        (0.05, 1.0e-3, 0.01353),
        (0.05, 0.1, 0.00662),
    ] {
        let label = format!("power-law lambda={lambda}");
        let inflaton = PowerLawInflaton { lambda };
        record(
            &mut predictions,
            &mut trajectories,
            &label,
            inflaton,
            s_star,
            center,
        )?;
    }
    for (mu, phi_e, s_star, center) in [
        (10.0, 1.0, 1.0e-3, 0.01077),
        (10.0, 1.0, 0.1, 0.00486),
        (6.0, 0.06, 1.0e-3, 0.00558),
        (6.0, 0.06, 0.1, 0.00340),
    ] {
        let label = format!("hybrid mu={mu} phi_e={phi_e}");
        let inflaton = HybridInflaton { mu, phi_e };
        record(
            &mut predictions,
            &mut trajectories,
            &label,
            inflaton,
            s_star,
            center,
        )?;
    }
    println!(
        "wrote assisted_predictions.csv and assisted_trajectories.csv in {}",
        directory.display()
    );
    Ok(())
}
