use phys_type::{ContravariantField, CovariantField, M1, Quantity};

fn main() {
    let upper = ContravariantField::<M1, f64, 2>::new([
        Quantity::new(1.0),
        Quantity::new(2.0),
    ]);
    let lower = CovariantField::<M1, f64, 2>::new([
        Quantity::new(1.0),
        Quantity::new(2.0),
    ]);

    let _invalid = upper + lower;
}
