use phys_type::{
    CanonicalCubicHamiltonianCoefficients, FieldPerturbation, M3Over2, Quantity,
    ThreePointBlock,
};

fn main() {
    let wrong_a = ThreePointBlock::<
        FieldPerturbation,
        FieldPerturbation,
        FieldPerturbation,
        M3Over2,
        f64,
        1,
    >::new([[[Quantity::new(0.0)]]]);

    let _ = CanonicalCubicHamiltonianCoefficients::<1> {
        a_first_raised: wrong_a,
        b_first_raised: ThreePointBlock::new([[[Quantity::new(0.0)]]]),
        b_third_raised: ThreePointBlock::new([[[Quantity::new(0.0)]]]),
        c_first_raised: ThreePointBlock::new([[[Quantity::new(0.0)]]]),
        c_third_raised: ThreePointBlock::new([[[Quantity::new(0.0)]]]),
    };
}
