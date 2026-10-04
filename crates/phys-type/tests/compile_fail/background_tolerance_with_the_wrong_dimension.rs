use phys_type::{Field, FieldVelocity, component_error_ratio};

fn main() {
    let value = Field::new(1.0);

    // A tolerance for field velocities cannot bound a field value.
    let _invalid = component_error_ratio(value, value, value, FieldVelocity::new(1.0e-12), 1.0e-8);
}
