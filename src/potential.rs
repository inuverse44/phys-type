//! Scalar field potentials.

use crate::{EnergyDensity, Field, Mass, PotentialGradient};

/// A scalar field potential and its first derivative.
///
/// The gradient is supplied explicitly here. Automatic differentiation can
/// later be added by using an automatic-differentiation number as the value
/// type of `Quantity` and adapting potentials to it, without changing this
/// trait's meaning.
pub trait Potential {
    /// The potential energy density `V(phi)`.
    fn value(&self, phi: Field) -> EnergyDensity;

    /// The gradient `dV/dphi`.
    fn gradient(&self, phi: Field) -> PotentialGradient;
}

/// The quadratic potential `V(phi) = (1/2) m^2 phi^2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticPotential {
    /// The field mass `m`.
    pub mass: Mass,
}

impl QuadraticPotential {
    /// Builds a quadratic potential with field mass `mass`.
    pub const fn new(mass: Mass) -> Self {
        Self { mass }
    }
}

impl Potential for QuadraticPotential {
    fn value(&self, phi: Field) -> EnergyDensity {
        self.mass * self.mass * phi * phi * 0.5
    }

    fn gradient(&self, phi: Field) -> PotentialGradient {
        self.mass * self.mass * phi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_value_matches_the_closed_form() {
        let potential = QuadraticPotential::new(Mass::new(2.0));

        let value = potential.value(Field::new(3.0));

        assert_eq!(value.into_value(), 18.0);
    }

    #[test]
    fn quadratic_gradient_matches_the_closed_form() {
        let potential = QuadraticPotential::new(Mass::new(2.0));

        let gradient = potential.gradient(Field::new(3.0));

        assert_eq!(gradient.into_value(), 12.0);
    }

    #[test]
    fn gradient_agrees_with_a_finite_difference() {
        let potential = QuadraticPotential::new(Mass::new(1.0e-6));
        let phi = 16.0;
        let step = 1.0e-4;

        let forward = potential.value(Field::new(phi + step)).into_value();
        let backward = potential.value(Field::new(phi - step)).into_value();
        let numerical = (forward - backward) / (2.0 * step);
        let analytic = potential.gradient(Field::new(phi)).into_value();

        assert!((numerical - analytic).abs() < 1.0e-18);
    }
}
