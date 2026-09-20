//! A field and a time have different mass dimensions, so they cannot be added.

use phys_type::{Field, Time};

fn main() {
    let phi = Field::new(1.0);
    let t = Time::new(1.0);

    let _invalid = phi + t;
}
