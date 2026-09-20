//! A derivative must be multiplied by a time before it can advance a state.

use phys_type::{Field, FieldVelocity, Quantity, ScaleFactor, State, StateDerivative};

fn main() {
    let state = State::new(
        Field::new(16.0),
        FieldVelocity::new(0.0),
        ScaleFactor::new(1.0),
    );
    let derivative = StateDerivative {
        dphi_dt: Quantity::new(0.0),
        dpi_dt: Quantity::new(0.0),
        da_dt: Quantity::new(0.0),
    };

    let _invalid = state + derivative;
}
