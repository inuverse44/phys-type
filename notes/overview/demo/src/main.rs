//! Demonstrations of phys-type for the overview note. This binary is the
//! Friedmann equation H = sqrt((phi_dot^2 / 2 + V) / (3 M_pl^2)), with every
//! intermediate dimension checked by the compiler.
//!
//!     cargo run                            # compiles and prints H
//!     cargo build --example forgot_square  # rejected: M^(3/2) is not M
//!     cargo build --example add_velocity   # rejected: M^2 + M^4
//!     cargo run --example slow_roll        # epsilon_H in Planck units and GeV
//!     cargo build --example forgot_planck_mass   # rejected: M^2 is not M^0
//!     cargo build --example conformal_time_step  # rejected: eta is not t
//!     cargo build --example pivot_in_mpc         # rejected: Mpc^-1 is not k/M_pl
//!     cargo run --example jet              # one formula, f64 and Jet
//!     cargo run --example tensor           # field-space tensors; see its
//!                                          # header for the broken ones

use phys_type::{EnergyDensity, FieldVelocity, Hubble, M4, Mass, Quantity};

fn main() {
    let phi_dot = FieldVelocity::new(1.0e-6); // d phi / dt, M^2
    let potential = EnergyDensity::new(1.0e-10); // V, M^4
    let m_pl = Mass::new(1.0); // M

    // M^2 * M^2 + M^4 = M^4
    let rho: Quantity<M4, f64> = phi_dot * phi_dot * 0.5 + potential;
    // sqrt(M^4 / M^2) = M
    let hubble: Hubble = (rho / (m_pl * m_pl * 3.0)).sqrt();

    println!("H = {:.6e} M_pl", hubble.into_value());
}
