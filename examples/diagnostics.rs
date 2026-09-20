//! Diagnostic runs used to inspect the integrator's behaviour.
//!
//! Writes three CSV files into the directory given as the first argument
//! (default `diagnostics/`):
//!
//! - `trajectory.csv`  the reference quadratic background
//! - `attractor.csv`   several initial velocities converging on slow roll
//! - `convergence.csv` the Runge-Kutta error as a function of step size

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::{
    Field, FieldVelocity, Mass, Potential, QuadraticPotential, RunSettings, ScalarFieldSystem,
    ScaleFactor, State, Time, rk4_step,
};

const MASS: f64 = 1.0e-6;
const PHI_START: f64 = 16.0;

fn system() -> ScalarFieldSystem<QuadraticPotential> {
    ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(MASS)))
}

/// The reference background, sampled along the e-folding number.
fn trajectory(dir: &Path) -> io::Result<()> {
    let system = system();
    let initial = system.slow_roll_initial_state(Field::new(PHI_START), ScaleFactor::new(1.0));
    let run = system.run(initial, RunSettings::default());

    let mut file = fs::File::create(dir.join("trajectory.csv"))?;
    run.write_csv(&mut file)
}

/// Phase-space convergence onto the slow-roll attractor.
///
/// Each run starts at the same `phi` but with a very different `dphi/dt`,
/// spanning kinetic-dominated and uphill initial conditions. Hubble friction
/// should pull all of them onto the same slow-roll curve within a few Hubble
/// times.
fn attractor(dir: &Path) -> io::Result<()> {
    let system = system();
    let reference = system.slow_roll_initial_state(Field::new(PHI_START), ScaleFactor::new(1.0));
    let slow_roll = reference.pi.into_value();

    let starts = [
        ("slow roll", slow_roll),
        ("at rest", 0.0),
        ("10x downhill", slow_roll * 10.0),
        ("100x downhill", slow_roll * 100.0),
        ("1000x downhill", slow_roll * 1000.0),
        ("100x uphill", -slow_roll * 100.0),
    ];

    let mut file = fs::File::create(dir.join("attractor.csv"))?;
    writeln!(
        file,
        "label,step,t,e_folds,phi,dphi,H,epsilon_h,slow_roll_dphi"
    )?;

    // Sample densely at early times, where the fast transient happens, and
    // sparsely later.
    let sampled = |step: usize| match step {
        0..=200 => step.is_multiple_of(5),
        201..=2_000 => step.is_multiple_of(50),
        2_001..=20_000 => step.is_multiple_of(500),
        _ => step.is_multiple_of(2_000),
    };

    let dt = Time::new(20.0);
    for (label, pi) in starts {
        let initial = State::new(
            Field::new(PHI_START),
            FieldVelocity::new(pi),
            ScaleFactor::new(1.0),
        );
        let mut state = initial;

        for step in 0..=100_000 {
            if sampled(step) {
                writeln!(
                    file,
                    "{label},{step},{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
                    20.0 * step as f64,
                    state.e_folds(&initial),
                    state.phi.into_value(),
                    state.pi.into_value(),
                    system.hubble(&state).into_value(),
                    system.epsilon_h(&state),
                    system.slow_roll_velocity(&state).into_value(),
                )?;
            }
            state = rk4_step(&system, &state, dt);
        }
    }

    Ok(())
}

/// The error of one fixed-step run against a much finer reference run.
///
/// For classical RK4 the global error should fall as `dt^4`. The step sizes
/// have to start coarse enough to be resolution-limited rather than
/// roundoff-limited: `1/H` is about `1.5e5` here, so `dt = 1e5` is a genuinely
/// coarse step.
fn convergence(dir: &Path) -> io::Result<()> {
    let system = system();
    let initial = system.slow_roll_initial_state(Field::new(PHI_START), ScaleFactor::new(1.0));
    let horizon = 4.0e6;

    let integrate = |dt: f64| -> State {
        let steps = (horizon / dt).round() as usize;
        let dt = Time::new(dt);
        let mut state = initial;
        for _ in 0..steps {
            state = rk4_step(&system, &state, dt);
        }
        state
    };

    let reference = integrate(25.0);
    let reference_phi = reference.phi.into_value();

    let mut file = fs::File::create(dir.join("convergence.csv"))?;
    writeln!(file, "dt,h_dt,steps,phi,abs_error,rel_error,order")?;

    let hubble = system.hubble(&initial).into_value();
    let mut previous: Option<(f64, f64)> = None;

    for exponent in 0..10 {
        let dt = 2.0e5 / 2.0_f64.powi(exponent);
        let state = integrate(dt);
        let error = (state.phi.into_value() - reference_phi).abs();

        // The empirical convergence order between successive halvings.
        let order = match previous {
            Some((prev_dt, prev_error)) if error > 0.0 => {
                (prev_error / error).log2() / (prev_dt / dt).log2()
            }
            _ => f64::NAN,
        };
        previous = Some((dt, error));

        writeln!(
            file,
            "{dt:.6},{:.6e},{},{:.15e},{error:.6e},{:.6e},{order:.4}",
            hubble * dt,
            (horizon / dt).round() as usize,
            state.phi.into_value(),
            error / reference_phi.abs(),
        )?;
    }

    Ok(())
}

fn main() -> io::Result<()> {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "diagnostics".into());
    let dir = Path::new(&dir);
    fs::create_dir_all(dir)?;

    trajectory(dir)?;
    attractor(dir)?;
    convergence(dir)?;

    // Touch the trait import so the example also exercises it directly.
    let potential = QuadraticPotential::new(Mass::new(MASS));
    eprintln!(
        "V(16 M_pl) = {:.6e} M_pl^4",
        potential.value(Field::new(PHI_START)).into_value()
    );
    eprintln!(
        "wrote trajectory.csv, attractor.csv and convergence.csv to {}",
        dir.display()
    );

    Ok(())
}
