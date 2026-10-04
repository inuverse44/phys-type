//! `sqrt(M^(1/2))` would be `M^(1/4)`, which is off the supported grid.

use phys_type::{MHalf, Quantity};

fn main() {
    let root: Quantity<MHalf, f64> = Quantity::new(4.0);

    let _invalid = root.sqrt();
}
