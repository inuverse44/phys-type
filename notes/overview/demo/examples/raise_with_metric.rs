//! G_IJ mistaken for G^IJ when raising the index of V_,I. The metric has two
//! lower indices and only lowers, so this does not compile.

use phys_type::{CovariantField, Field, FieldSpaceGeometry, M3, PolarFieldSpace, Quantity};

fn main() {
    let position = [Field::new(2.0), Field::new(0.5)];
    let gradient = CovariantField::<M3, f64, 2>::new([Quantity::new(1.0e-3), Quantity::new(4.0e-3)]);

    let metric = PolarFieldSpace.metric(&position);
    let force = metric.raise(&gradient);

    println!("{:?}", force);
}
