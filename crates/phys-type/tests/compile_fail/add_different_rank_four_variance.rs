use phys_type::{FieldTensor4, Lower, MMinus2, Quantity, Upper};

fn main() {
    let upper = FieldTensor4::<Upper, Lower, Lower, Lower, MMinus2, f64, 2>::new(
        [[[[Quantity::new(0.0); 2]; 2]; 2]; 2],
    );
    let lower = FieldTensor4::<Lower, Lower, Lower, Lower, MMinus2, f64, 2>::new(
        [[[[Quantity::new(0.0); 2]; 2]; 2]; 2],
    );

    let _invalid = upper + lower;
}
