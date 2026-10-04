//! A circular orbit of constant turn rate and constant entropy mass.
//!
//! `V = V_0 - c theta + m^2 (rho - rho_0)^2 / 2` on `PolarFieldSpace` has a
//! circular orbit of radius `rho` with turn rate `Omega`, mixing
//! `lambda = 2 |Omega| / H`, and entropy mass `mu^2 = m^2 + 3 Omega^2`, the
//! background of the exact solutions of arXiv:2606.18248 and 2607.14529.

use std::ops::{Add, Mul};

use phys_type::{
    CovariantMultifieldSystem, EndOfInflation, EnergyDensity, Field, FieldVelocity,
    MultifieldPotentialValue, MultifieldState, PolarFieldSpace, Quantity, ScaleFactor,
    SpectrumRunSettings,
};

/// The Hubble rate of the orbits, in reduced Planck units.
pub const HUBBLE: f64 = 1.0e-5;

/// `V = V_0 - c theta + m^2 (rho - rho_0)^2 / 2` on `PolarFieldSpace`.
#[derive(Clone, Copy, Debug)]
pub struct ConstantTurnPotential {
    v0: f64,
    slope: f64,
    mass_squared: f64,
    rho0: f64,
}

impl<V> MultifieldPotentialValue<V, 2> for ConstantTurnPotential
where
    V: From<f64> + Clone + Add<V, Output = V> + Mul<V, Output = V> + Mul<f64, Output = V>,
{
    fn value(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
        let [rho, theta] = fields.map(Quantity::into_value);
        let radial = rho + V::from(-self.rho0);
        EnergyDensity::new(
            V::from(self.v0)
                + theta * (-self.slope)
                + radial.clone() * radial * (0.5 * self.mass_squared),
        )
    }
}

/// The engine system of a constant-turn orbit.
pub type ConstantTurnSystem = CovariantMultifieldSystem<ConstantTurnPotential, PolarFieldSpace, 2>;

/// The orbit for design values `lambda` and `mu^2 / H^2` at radius
/// `radius`, its initial state, and spectrum settings whose pivot leaves the
/// horizon `after_exit + 1` e-folds before the end, `before_pivot` e-folds
/// after the start. `mu^2 = 3 lambda^2 / 4` leaves no radial restoring force.
pub fn orbit(
    lambda: f64,
    mu_squared: f64,
    radius: f64,
    after_exit: f64,
    before_pivot: f64,
) -> (ConstantTurnSystem, MultifieldState<2>, SpectrumRunSettings) {
    let omega = 0.5 * lambda * HUBBLE;
    let mass_squared = (mu_squared - 0.75 * lambda * lambda) * HUBBLE * HUBBLE;
    // Radial balance m^2 (rho - rho_0) = rho Omega^2 and the angular
    // slow-roll attractor c = 3 H rho^2 Omega.
    let potential = ConstantTurnPotential {
        v0: 3.0 * HUBBLE * HUBBLE,
        slope: 3.0 * HUBBLE * radius * radius * omega,
        mass_squared,
        rho0: radius - radius * omega * omega / mass_squared,
    };
    let pivot_before_end = after_exit + 1.0;
    let total = pivot_before_end + before_pivot;
    let settings = SpectrumRunSettings {
        end_of_inflation: EndOfInflation::FieldValue {
            field: 1,
            value: Field::new(omega / HUBBLE * total),
        },
        pivot_e_folds_before_end: pivot_before_end,
        output_e_folds_before_end: 1.0,
        k_min_over_pivot: 1.0,
        k_max_over_pivot: 1.0,
        k_count: 1,
        background_sample_every: 1,
        ..SpectrumRunSettings::default()
    };
    let initial = MultifieldState::new(
        [Field::new(radius), Field::new(0.0)],
        [FieldVelocity::new(0.0), FieldVelocity::new(omega)],
        ScaleFactor::new(1.0),
    );
    (
        CovariantMultifieldSystem::in_planck_units(potential, PolarFieldSpace),
        initial,
        settings,
    )
}

/// `V_NN / H_0^2 = m^2 / H_0^2` for design values, with `H_0 = HUBBLE`.
pub fn radial_mass_squared(lambda: f64, mu_squared: f64) -> f64 {
    mu_squared - 0.75 * lambda * lambda
}

/// `epsilon = rho^2 lambda^2 / 8` fixes the orbit radius.
pub fn radius(lambda: f64, epsilon: f64) -> f64 {
    (8.0 * epsilon).sqrt() / lambda
}
