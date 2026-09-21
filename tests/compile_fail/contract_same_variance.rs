use phys_type::{ContravariantField, M1, Quantity};

fn main() {
    let lhs = ContravariantField::<M1, f64, 2>::new([
        Quantity::new(1.0),
        Quantity::new(2.0),
    ]);
    let rhs = ContravariantField::<M1, f64, 2>::new([
        Quantity::new(3.0),
        Quantity::new(4.0),
    ]);

    let _invalid = lhs.contract(&rhs);
}
