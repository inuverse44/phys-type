//! `V(phi)` is an energy density `M^4`, not a field value `M^1`.

use phys_type::{Field, Mass, Potential, QuadraticPotential};

fn main() {
    let potential = QuadraticPotential::new(Mass::new(1.0e-6));

    let _invalid: Field = potential.value(Field::new(16.0));
}
