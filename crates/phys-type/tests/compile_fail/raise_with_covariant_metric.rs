use phys_type::{CovariantField, FieldMetric, M1, Quantity};

fn main() {
    let metric = FieldMetric::<2>::euclidean();
    let covector = CovariantField::<M1, f64, 2>::new([
        Quantity::new(1.0),
        Quantity::new(2.0),
    ]);

    let _invalid = metric.raise(&covector);
}
