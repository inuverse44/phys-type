//! Dimensions must be attached explicitly; a raw number is not a quantity.

use phys_type::Field;

fn main() {
    let _invalid: Field = 0.17;
}
