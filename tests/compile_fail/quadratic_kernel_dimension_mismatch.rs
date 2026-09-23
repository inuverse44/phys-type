use phys_type::{
    FieldPerturbation, MHalf, MomentumPerturbation, QuadraticPhaseSpaceKernel, Quantity,
    ThreePointBlock,
};

fn main() {
    let mut kernel = QuadraticPhaseSpaceKernel::<1>::zero();
    let wrong = ThreePointBlock::<
        MomentumPerturbation,
        FieldPerturbation,
        FieldPerturbation,
        MHalf,
        f64,
        1,
    >::new([[[Quantity::new(1.0)]]]);

    kernel.p_qq = wrong;
}
