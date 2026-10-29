//! One potential formula, written once over the number type `V`, evaluated
//! with `f64` (the value) and with `Jet<1>` (the value, V' and V'').
//!
//! V(phi) = lambda phi^4 / 4 at phi = 2 with lambda = 1:
//! V = 4, V' = lambda phi^3 = 8, V'' = 3 lambda phi^2 = 12.

use core::ops::Mul;

use phys_type::{EnergyDensity, Field, Jet};

/// lambda phi^4 / 4 with lambda = 1, generic over the number type.
fn quartic<V>(phi: Field<V>) -> EnergyDensity<V>
where
    V: Copy + Mul<Output = V> + Mul<f64, Output = V>,
{
    // M * M * M * M = M^4, whatever V is.
    phi * phi * phi * phi * 0.25
}

fn main() {
    let value = quartic(Field::new(2.0_f64));
    println!("f64:    V = {}", value.into_value());

    // The active variable phi = 2 with seed dphi/dphi = 1.
    let jet = quartic(Field::new(Jet::<1>::variable(2.0, 0))).into_value();
    println!(
        "Jet<1>: V = {}, V' = {}, V'' = {}",
        jet.value, jet.gradient[0], jet.hessian[0][0]
    );
}
