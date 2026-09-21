use phys_type::{Complex64, ModeAmplitude, ModeVelocity};

fn main() {
    let amplitude = ModeAmplitude::new(Complex64::new(1.0, 2.0));
    let velocity = ModeVelocity::new(Complex64::new(3.0, 4.0));

    let _invalid = amplitude + velocity;
}
