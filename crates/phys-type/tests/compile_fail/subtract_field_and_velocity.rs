//! `phi` is `M^1` and `dphi/dt` is `M^2`, so their difference is meaningless.

use phys_type::{Field, FieldVelocity};

fn main() {
    let phi = Field::new(1.0);
    let pi = FieldVelocity::new(1.0);

    let _invalid = phi - pi;
}
