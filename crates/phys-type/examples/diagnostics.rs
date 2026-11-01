//! Diagnostic runs used to inspect the integrator's behaviour.
//!
//! Writes four CSV files into the directory given as the first argument
//! (default `diagnostics/`):
//!
//! - `trajectory.csv`  the reference quadratic background
//! - `attractor.csv`   several initial velocities converging on slow roll
//! - `convergence.csv` the Runge-Kutta error as a function of step size
//! - `oscillation.csv` ten post-inflation oscillations of the scalar field

use std::fs;
use std::io::{self, Write};
use std::path::Path;

use phys_type::{
    Field, FieldVelocity, Mass, PotentialValue, QuadraticPotential, RunSettings, ScalarFieldSystem,
    ScaleFactor, State, Time, rk4_step,
};

const MASS: f64 = 1.0e-6;
const PHI_START: f64 = 16.0;
const STEPS_PER_SHORTEST_SCALE: f64 = 100.0;
const OSCILLATION_ZERO_CROSSINGS: usize = 20;
const POST_INFLATION_MAX_STEPS: usize = 1_000_000;

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
        let mut time = Time::new(0.0);

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
            state = rk4_step(&system, time, &state, dt);
            time = time + dt;
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
        let mut time = Time::new(0.0);
        for _ in 0..steps {
            state = rk4_step(&system, time, &state, dt);
            time = time + dt;
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

/// Chooses a step from the shorter of the Hubble and quadratic-oscillation
/// time scales: `dt = 1 / (100 max(H, m))`.
fn physical_step(
    system: &ScalarFieldSystem<QuadraticPotential>,
    state: &State,
    steps_per_scale: f64,
) -> Time {
    let frequency = system.hubble(state).into_value().max(MASS);
    Time::new(1.0 / (steps_per_scale * frequency))
}

/// Continues a state until the exact Hubble slow-roll condition reaches one.
fn advance_to_inflation_end(
    system: &ScalarFieldSystem<QuadraticPotential>,
    mut state: State,
    mut time: Time,
) -> io::Result<(State, Time)> {
    for _ in 0..POST_INFLATION_MAX_STEPS {
        if system.epsilon_h(&state) >= 1.0 {
            return Ok((state, time));
        }

        let dt = physical_step(system, &state, STEPS_PER_SHORTEST_SCALE);
        state = rk4_step(system, time, &state, dt);
        time = time + dt;
    }

    Err(io::Error::other(
        "inflation did not end before the post-inflation step limit",
    ))
}

#[derive(Clone, Copy, Debug)]
struct OscillationSummary {
    elapsed: f64,
    a_ratio: f64,
    mean_w: f64,
    steps: usize,
}

/// Evolves through the requested zero crossings and optionally observes every
/// row. Keeping the evolution in one function makes the CSV and convergence
/// calculation follow exactly the same path.
fn evolve_oscillations<F>(
    system: &ScalarFieldSystem<QuadraticPotential>,
    initial_state: State,
    initial_time: Time,
    steps_per_scale: f64,
    mut observe: F,
) -> io::Result<OscillationSummary>
where
    F: FnMut(usize, usize, f64, &State, Time) -> io::Result<()>,
{
    let mut state = initial_state;
    let mut time = initial_time;
    let initial_time = initial_time.into_value();
    let initial_scale_factor = initial_state.a.into_value();
    let mut previous_phi = state.phi.into_value();
    let mut zero_crossings = 0;
    let mut weighted_w = 0.0;
    let mut accumulated_time = 0.0;

    for step in 0..=POST_INFLATION_MAX_STEPS {
        let elapsed = time.into_value() - initial_time;
        let dt = physical_step(system, &state, steps_per_scale);
        observe(step, zero_crossings, elapsed, &state, dt)?;

        if zero_crossings >= OSCILLATION_ZERO_CROSSINGS {
            return Ok(OscillationSummary {
                elapsed,
                a_ratio: state.a.into_value() / initial_scale_factor,
                mean_w: weighted_w / accumulated_time,
                steps: step,
            });
        }

        let dt_value = dt.into_value();
        weighted_w += system.equation_of_state(&state) * dt_value;
        accumulated_time += dt_value;

        let next = rk4_step(system, time, &state, dt);
        let next_phi = next.phi.into_value();
        if previous_phi != 0.0
            && next_phi != 0.0
            && previous_phi.is_sign_positive() != next_phi.is_sign_positive()
        {
            zero_crossings += 1;
        }

        previous_phi = next_phi;
        state = next;
        time = time + dt;
    }

    Err(io::Error::other(
        "the requested oscillations did not finish before the step limit",
    ))
}

fn relative_difference(lhs: f64, rhs: f64) -> f64 {
    (lhs - rhs).abs() / rhs.abs()
}

/// Continues the reference background through ten coherent oscillations.
///
/// This models a non-decaying scalar in a quadratic potential. It tests the
/// matter-like averaged behaviour, not particle production or reheating.
fn oscillation(dir: &Path) -> io::Result<()> {
    let system = system();
    let initial = system.slow_roll_initial_state(Field::new(PHI_START), ScaleFactor::new(1.0));
    let inflation = system.run(initial, RunSettings::default());
    let (end_state, end_time) =
        advance_to_inflation_end(&system, inflation.final_state, inflation.final_time)?;

    let end_time_value = end_time.into_value();
    let end_scale_factor = end_state.a.into_value();
    let end_hubble = system.hubble(&end_state).into_value();

    let mut file = fs::File::create(dir.join("oscillation.csv"))?;
    writeln!(
        file,
        "step,zero_crossings,t_since_end,m_t_since_end,phi,dphi,a_over_a_end,H,epsilon_h,w,matter_a_over_a_end,dt"
    )?;

    let standard = evolve_oscillations(
        &system,
        end_state,
        end_time,
        STEPS_PER_SHORTEST_SCALE,
        |step, zero_crossings, elapsed, state, dt| {
            let a_ratio = state.a.into_value() / end_scale_factor;
            let matter_a_ratio = (1.0 + 1.5 * end_hubble * elapsed).powf(2.0 / 3.0);

            writeln!(
                file,
                "{step},{zero_crossings},{elapsed:.9e},{:.9e},{:.9e},{:.9e},{a_ratio:.9e},{:.9e},{:.9e},{:.9e},{matter_a_ratio:.9e},{:.9e}",
                MASS * elapsed,
                state.phi.into_value(),
                state.pi.into_value(),
                system.hubble(state).into_value(),
                system.epsilon_h(state),
                system.equation_of_state(state),
                dt.into_value(),
            )
        },
    )?;

    let refined = evolve_oscillations(
        &system,
        end_state,
        end_time,
        STEPS_PER_SHORTEST_SCALE * 2.0,
        |_, _, _, _, _| Ok(()),
    )?;

    let elapsed_difference = relative_difference(standard.elapsed, refined.elapsed);
    let scale_factor_difference = relative_difference(standard.a_ratio, refined.a_ratio);
    let mean_w_difference = (standard.mean_w - refined.mean_w).abs();
    let tolerance = 1.0e-3;

    if elapsed_difference > tolerance
        || scale_factor_difference > tolerance
        || mean_w_difference > tolerance
    {
        return Err(io::Error::other(format!(
            "post-inflation step convergence failed: elapsed={elapsed_difference:.3e}, a={scale_factor_difference:.3e}, mean_w={mean_w_difference:.3e}"
        )));
    }

    eprintln!(
        "post-inflation: {} zero crossings in {:.6} oscillation time units; <w>={:.6}; half-step differences: time={:.3e}, a={:.3e}, <w>={:.3e}; t_end={end_time_value:.6e}",
        OSCILLATION_ZERO_CROSSINGS,
        MASS * standard.elapsed,
        standard.mean_w,
        elapsed_difference,
        scale_factor_difference,
        mean_w_difference,
    );
    debug_assert!(standard.steps > 0);

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
    oscillation(dir)?;

    // Touch the trait import so the example also exercises it directly.
    let potential = QuadraticPotential::new(Mass::new(MASS));
    eprintln!(
        "V(16 M_pl) = {:.6e} M_pl^4",
        potential.value(Field::new(PHI_START)).into_value()
    );
    eprintln!(
        "wrote trajectory.csv, attractor.csv, convergence.csv and oscillation.csv to {}",
        dir.display()
    );

    Ok(())
}
