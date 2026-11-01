use phys_type::{
    Complex64, Field, FieldVelocity, Mass, ModeAmplitude, ModeVelocity, MultifieldState,
    ScaleFactor, Modes, PerturbationState, field_power_spectrum, mode_covariance,
};

fn main() {
    let modes = Modes {
        rescaled_amplitude: [[ModeAmplitude::new(Complex64::new(1.0, 0.0)); 2]; 2],
        rescaled_velocity: [[ModeVelocity::new(Complex64::new(0.0, 0.0)); 2]; 2],
    };
    let state = PerturbationState {
        background: MultifieldState::new(
            [Field::new(1.0); 2],
            [FieldVelocity::new(0.0); 2],
            ScaleFactor::new(1.0),
        ),
        modes,
    };
    let k = Mass::new(1.0);
    let spectrum = field_power_spectrum(k, &mode_covariance(&state));

    let _invalid = field_power_spectrum(k, &spectrum);
}
