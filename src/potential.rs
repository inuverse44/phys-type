//! Scalar field potentials.

use core::ops::Mul;

use crate::{EnergyDensity, Field, Mass, PotentialGradient};

/// A scalar field potential formula evaluated with numerical value type `V`.
///
/// A potential can implement this once for every compatible value type, such
/// as `f64` today and an automatic-differentiation number in the future.
pub trait PotentialValue<V> {
    /// The potential energy density `V(phi)` evaluated with value type `V`.
    fn value(&self, phi: Field<V>) -> EnergyDensity<V>;
}

/// The temporary, manually differentiated interface used by the background
/// MVP.
///
/// Automatic differentiation will eventually replace the handwritten
/// gradient. Requiring [`PotentialValue<f64>`] keeps the current background
/// calculation available during that transition.
pub trait Potential: PotentialValue<f64> {
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

impl<V> PotentialValue<V> for QuadraticPotential
where
    V: Clone + Mul<V, Output = V> + Mul<f64, Output = V>,
    f64: Mul<V, Output = V>,
{
    fn value(&self, phi: Field<V>) -> EnergyDensity<V> {
        self.mass * self.mass * phi.clone() * phi * 0.5
    }
}

impl Potential for QuadraticPotential {
    fn gradient(&self, phi: Field) -> PotentialGradient {
        self.mass * self.mass * phi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Copy, Debug, PartialEq)]
    struct TestNumber(f64);

    impl Mul for TestNumber {
        type Output = Self;

        fn mul(self, rhs: Self) -> Self::Output {
            Self(self.0 * rhs.0)
        }
    }

    impl Mul<f64> for TestNumber {
        type Output = Self;

        fn mul(self, rhs: f64) -> Self::Output {
            Self(self.0 * rhs)
        }
    }

    impl Mul<TestNumber> for f64 {
        type Output = TestNumber;

        fn mul(self, rhs: TestNumber) -> Self::Output {
            TestNumber(self * rhs.0)
        }
    }

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
    fn quadratic_formula_is_not_fixed_to_f64() {
        let potential = QuadraticPotential::new(Mass::new(2.0));

        let value = potential.value(Field::new(TestNumber(3.0)));

        assert_eq!(value.into_value(), TestNumber(18.0));
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
