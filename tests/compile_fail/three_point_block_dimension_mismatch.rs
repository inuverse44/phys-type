use phys_type::{
    Complex64, FieldPerturbation, MMinus3Over2, MMinusHalf, MomentumPerturbation, Quantity,
    ThreePointBlock,
};

fn main() {
    let qqp = ThreePointBlock::<
        FieldPerturbation,
        FieldPerturbation,
        MomentumPerturbation,
        MMinusHalf,
        Complex64,
        2,
    >::new([[[Quantity::new(Complex64::new(0.0, 0.0)); 2]; 2]; 2]);

    let _: ThreePointBlock<
        FieldPerturbation,
        FieldPerturbation,
        MomentumPerturbation,
        MMinus3Over2,
        Complex64,
        2,
    > = qqp;
}
