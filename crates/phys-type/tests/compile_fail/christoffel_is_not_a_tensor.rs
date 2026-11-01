use phys_type::{ChristoffelSymbols, FieldTensor3, Lower, MMinus1, Upper};

fn main() {
    let connection = ChristoffelSymbols::<2>::zero();

    let _invalid: FieldTensor3<Upper, Lower, Lower, MMinus1, f64, 2> = connection;
}
