use phys_type::{
    Complex64, Mass, ModeAmplitude, ModeVelocity, ScaleFactor, TwoFieldModes,
    field_power_spectrum, mode_covariance,
};

fn main() {
    let modes = TwoFieldModes {
        rescaled_amplitude: [[ModeAmplitude::new(Complex64::new(1.0, 0.0)); 2]; 2],
        rescaled_velocity: [[ModeVelocity::new(Complex64::new(0.0, 0.0)); 2]; 2],
    };
    let k = Mass::new(1.0);
    let spectrum = field_power_spectrum(k, &mode_covariance(&modes, ScaleFactor::new(1.0)));

    let _invalid = field_power_spectrum(k, &spectrum);
}
