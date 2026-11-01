//! epsilon_H written as phi_dot^2 / (2 H^2), the M_pl^2 dropped as one does
//! in Planck units. The result has dimension M^2, not M^0, so this does not
//! compile.

use phys_type::{Dimensionless, FieldVelocity, Hubble};

fn main() {
    let hubble = Hubble::new(1.0e-5);
    let phi_dot = FieldVelocity::new(1.0e-6);

    let epsilon: Dimensionless = phi_dot * phi_dot / (hubble * hubble * 2.0);

    println!("epsilon_H = {}", epsilon.into_value());
}
