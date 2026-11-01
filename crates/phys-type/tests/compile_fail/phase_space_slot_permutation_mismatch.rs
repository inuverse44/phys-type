use phys_type::{
    FieldPerturbation, M0, MomentumPerturbation, Quantity, TwoPointBlock,
};

fn main() {
    let block = TwoPointBlock::<FieldPerturbation, MomentumPerturbation, M0, f64, 1>::new([[
        Quantity::new(1.0),
    ]]);

    let _: TwoPointBlock<FieldPerturbation, FieldPerturbation, M0, f64, 1> = block;
}
