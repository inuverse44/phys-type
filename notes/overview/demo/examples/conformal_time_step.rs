//! A conformal-time step passed to code written in cosmic time. Both have
//! dimension M^-1, so a dimension check alone cannot tell them apart, yet
//! d eta = dt / a differs by e^60 near the end of inflation. ConformalTime is
//! a distinct type, so this does not compile.

use phys_type::{ConformalTime, Time};

/// Advances a background by a cosmic-time step dt.
fn advance_in_cosmic_time(dt: Time) -> f64 {
    dt.into_value()
}

fn main() {
    let d_eta = ConformalTime::new(1.0e3);

    let elapsed = advance_in_cosmic_time(d_eta);

    println!("{elapsed}");
}
