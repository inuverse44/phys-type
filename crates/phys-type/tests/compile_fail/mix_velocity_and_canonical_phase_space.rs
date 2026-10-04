use num_complex::Complex64;
use phys_type::{
    CanonicalQuadraticPhaseSpaceKernel, Quantity, TwoPointBlock,
    TwoPointPhaseSpaceCovariance, three_point_source_from_kernels,
};

fn main() {
    let kernel = CanonicalQuadraticPhaseSpaceKernel::<1>::zero();
    let velocity_covariance = TwoPointPhaseSpaceCovariance::<1> {
        qq: TwoPointBlock::new([[Quantity::new(Complex64::new(1.0, 0.0))]]),
        qp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
        pq: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
        pp: TwoPointBlock::new([[Quantity::new(Complex64::new(1.0, 0.0))]]),
    };

    let _ = three_point_source_from_kernels(
        [&kernel, &kernel, &kernel],
        [
            &velocity_covariance,
            &velocity_covariance,
            &velocity_covariance,
        ],
    );
}
