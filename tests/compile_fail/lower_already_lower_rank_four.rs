use phys_type::{FieldMetric, FieldTensor4, Lower, LowerFirst, MMinus2, Quantity};

fn main() {
    let metric = FieldMetric::<2>::euclidean();
    let tensor = FieldTensor4::<Lower, Lower, Lower, Lower, MMinus2, f64, 2>::new(
        [[[[Quantity::new(0.0); 2]; 2]; 2]; 2],
    );

    let _invalid = metric.lower_first(&tensor);
}
