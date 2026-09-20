//! Integrates a quadratic single-field inflationary background and writes the
//! result as CSV.
//!
//! ```sh
//! cargo run --release --example quadratic_inflation > inflation.csv
//! ```
//!
//! The CSV goes to standard output and a short summary goes to standard error,
//! so the two can be separated by redirection. Plot it with, for example:
//!
//! ```text
//! gnuplot -e "set datafile separator ','; \
//!             set key autotitle columnhead; \
//!             plot 'inflation.csv' using 3:4 with lines; pause -1"
//! ```

use std::io::{self, Write};

use phys_type::{
    Field, Mass, QuadraticPotential, RunSettings, ScalarFieldSystem, ScaleFactor, Termination,
};

fn main() -> io::Result<()> {
    // V(phi) = (1/2) m^2 phi^2 with m = 10^-6 M_pl, in Planck units.
    let system = ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(1.0e-6)));

    // The field velocity comes from the slow-roll estimate, so only phi and a
    // have to be chosen.
    let initial = system.slow_roll_initial_state(Field::new(16.0), ScaleFactor::new(1.0));

    let run = system.run(initial, RunSettings::default());

    let stdout = io::stdout();
    let mut out = stdout.lock();
    run.write_csv(&mut out)?;
    out.flush()?;

    let last = run
        .records
        .last()
        .expect("the initial row is always recorded");
    let reason = match run.termination {
        Termination::EFoldsReached => "reached the requested e-folding number",
        Termination::InflationEnded => "inflation ended (epsilon_H >= 1)",
        Termination::StepLimitReached => "hit the step limit",
    };

    eprintln!("steps        {}", run.steps);
    eprintln!("rows         {}", run.records.len());
    eprintln!("stopped      {reason}");
    eprintln!("final N      {:.4}", last.e_folds);
    eprintln!("final phi    {:.6} M_pl", last.phi);
    eprintln!("final H      {:.6e} M_pl", last.hubble);
    eprintln!("final eps_H  {:.6e}", last.epsilon_h);

    Ok(())
}
