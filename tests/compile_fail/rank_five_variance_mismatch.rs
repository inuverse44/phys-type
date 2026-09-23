use phys_type::{
    CovariantRiemannDerivative, FieldTensor5, Lower, MMinus3, Quantity,
};

fn main() {
    let all_lower = FieldTensor5::<Lower, Lower, Lower, Lower, Lower, MMinus3, f64, 2>::new(
        [[[[[Quantity::new(0.0); 2]; 2]; 2]; 2]; 2],
    );

    let _invalid: CovariantRiemannDerivative<2> = all_lower;
}
