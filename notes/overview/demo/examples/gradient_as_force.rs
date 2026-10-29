//! The equation of motion phi_ddot^I = -3 H phi_dot^I - G^IJ V_,J written with
//! V_,I in place of G^IJ V_,J. Right in Cartesian coordinates, wrong in polar
//! ones; an upper and a lower index cannot be added, so this does not compile.

use phys_type::{ContravariantField, CovariantField, M3, Quantity};

fn main() {
    // 3 H phi_dot^I, already a vector of dimension M^3.
    let friction = ContravariantField::<M3, f64, 2>::new([Quantity::new(3.0e-4), Quantity::new(9.0e-4)]);
    let gradient = CovariantField::<M3, f64, 2>::new([Quantity::new(1.0e-3), Quantity::new(4.0e-3)]);

    let acceleration = -(friction + gradient);

    println!("{:?}", acceleration);
}
