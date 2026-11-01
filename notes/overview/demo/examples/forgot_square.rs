//! M_pl instead of M_pl^2 in the denominator: sqrt(M^4 / M) = M^(3/2), which
//! is not the dimension M of a Hubble rate, so this does not compile.

use phys_type::{EnergyDensity, FieldVelocity, Hubble, M4, Mass, Quantity};

fn main() {
    let phi_dot = FieldVelocity::new(1.0e-6);
    let potential = EnergyDensity::new(1.0e-10);
    let m_pl = Mass::new(1.0);

    let rho: Quantity<M4, f64> = phi_dot * phi_dot * 0.5 + potential;
    let hubble: Hubble = (rho / (m_pl * 3.0)).sqrt();

    println!("H = {:e}", hubble.into_value());
}
