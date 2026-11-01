//! phi_dot instead of phi_dot^2 / 2: M^2 + M^4 has no addition, so this does
//! not compile.

use phys_type::{EnergyDensity, FieldVelocity};

fn main() {
    let phi_dot = FieldVelocity::new(1.0e-6);
    let potential = EnergyDensity::new(1.0e-10);

    let rho = phi_dot + potential;

    println!("rho = {:e}", rho.into_value());
}
