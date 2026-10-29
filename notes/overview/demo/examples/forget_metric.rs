//! sigma_dot^2 written as phi_dot^I phi_dot^I, without the metric. Right in
//! Cartesian coordinates, wrong in polar ones; two upper indices cannot be
//! contracted, so this does not compile.

use phys_type::{ContravariantField, M2, Quantity};

fn main() {
    let velocity = ContravariantField::<M2, f64, 2>::new([Quantity::new(0.1), Quantity::new(0.3)]);

    let kinetic = velocity.contract(&velocity);

    println!("sigma_dot^2 = {}", kinetic.into_value());
}
