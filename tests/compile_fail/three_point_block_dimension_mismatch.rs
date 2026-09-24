use phys_type::{
    Complex64, FieldPerturbation, MMinus2, MMinus3, MomentumPerturbation, Quantity,
    ThreePointBlock,
};

fn main() {
    let qqp = ThreePointBlock::<
        FieldPerturbation,
        FieldPerturbation,
        MomentumPerturbation,
        MMinus2,
        Complex64,
        2,
    >::new([[[Quantity::new(Complex64::new(0.0, 0.0)); 2]; 2]; 2]);

    let _: ThreePointBlock<
        FieldPerturbation,
        FieldPerturbation,
        MomentumPerturbation,
        MMinus3,
        Complex64,
        2,
    > = qqp;
}
