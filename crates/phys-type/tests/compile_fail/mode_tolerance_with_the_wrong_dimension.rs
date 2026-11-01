use phys_type::{ModeVelocity, PerturbationTolerance};

fn main() {
    let mut tolerance = PerturbationTolerance::default();

    // A tolerance for mode velocities cannot bound mode amplitudes.
    tolerance.rescaled_amplitude = ModeVelocity::<f64>::new(1.0e-20);
}
