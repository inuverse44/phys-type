//! Scalar field potentials.

use core::ops::{Add, Mul};

use crate::{Dual, EnergyDensity, Field, Jet, Mass, PotentialGradient, PotentialHessian};

/// A scalar field potential formula evaluated with numerical value type `V`.
///
/// A potential can implement this once for every compatible value type, such
/// as `f64` today and an automatic-differentiation number in the future.
pub trait PotentialValue<V> {
    /// The potential energy density `V(phi)` evaluated with value type `V`.
    fn value(&self, phi: Field<V>) -> EnergyDensity<V>;
}

/// A scalar potential that can be evaluated normally and with automatic
/// differentiation.
pub trait Potential: PotentialValue<f64> + PotentialValue<Dual> {
    /// Computes `dV/dphi` by evaluating the potential with an active field.
    fn gradient(&self, phi: Field) -> PotentialGradient {
        let active = Field::new(Dual::variable(phi.into_value()));
        let result = <Self as PotentialValue<Dual>>::value(self, active).into_value();

        PotentialGradient::new(result.derivative)
    }
}

impl<T> Potential for T where T: PotentialValue<f64> + PotentialValue<Dual> {}

/// A multifield potential formula evaluated with numerical value type `V`.
pub trait MultifieldPotentialValue<V, const N: usize> {
    /// Evaluates `V(phi^1, ..., phi^N)`.
    fn value(&self, phi: [Field<V>; N]) -> EnergyDensity<V>;
}

/// A multifield potential with a gradient obtained by seeding one active
/// field at a time.
pub trait MultifieldPotential<const N: usize>:
    MultifieldPotentialValue<f64, N> + MultifieldPotentialValue<Dual, N>
{
    /// Computes all components `dV/dphi^I` by forward-mode autodiff.
    fn gradient(&self, phi: [Field; N]) -> [PotentialGradient; N] {
        core::array::from_fn(|active| {
            let fields = core::array::from_fn(|index| {
                let value = phi[index].into_value();
                Field::new(if index == active {
                    Dual::variable(value)
                } else {
                    Dual::constant(value)
                })
            });
            let result = <Self as MultifieldPotentialValue<Dual, N>>::value(self, fields);

            PotentialGradient::new(result.into_value().derivative)
        })
    }
}

impl<T, const N: usize> MultifieldPotential<N> for T where
    T: MultifieldPotentialValue<f64, N> + MultifieldPotentialValue<Dual, N>
{
}

/// Value, gradient, and coordinate Hessian from one potential evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PotentialDerivatives<const N: usize> {
    /// Potential energy density.
    pub value: EnergyDensity,
    /// Coordinate gradient `V_,I`.
    pub gradient: [PotentialGradient; N],
    /// Coordinate Hessian `V_,IJ`.
    pub hessian: [[PotentialHessian; N]; N],
}

/// A multifield potential that supports second-order automatic
/// differentiation.
pub trait TwiceDifferentiableMultifieldPotential<const N: usize>:
    MultifieldPotential<N> + MultifieldPotentialValue<Jet<N>, N>
{
    /// Evaluates the value, every gradient component, and the full Hessian in
    /// one second-order forward-mode pass.
    fn derivatives(&self, phi: [Field; N]) -> PotentialDerivatives<N> {
        let fields =
            core::array::from_fn(|index| Field::new(Jet::variable(phi[index].into_value(), index)));
        let result =
            <Self as MultifieldPotentialValue<Jet<N>, N>>::value(self, fields).into_value();

        PotentialDerivatives {
            value: EnergyDensity::new(result.value),
            gradient: result.gradient.map(PotentialGradient::new),
            hessian: result.hessian.map(|row| row.map(PotentialHessian::new)),
        }
    }

    /// Computes the coordinate Hessian `V_,IJ`.
    fn hessian(&self, phi: [Field; N]) -> [[PotentialHessian; N]; N] {
        self.derivatives(phi).hessian
    }
}

impl<T, const N: usize> TwiceDifferentiableMultifieldPotential<N> for T where
    T: MultifieldPotential<N> + MultifieldPotentialValue<Jet<N>, N>
{
}

/// The quadratic potential `V(phi) = (1/2) m^2 phi^2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticPotential {
    /// The field mass `m`.
    pub mass: Mass,
}

/// A separable quadratic potential
/// `V(phi) = sum_I (1/2) m_I^2 (phi^I)^2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticMultifieldPotential<const N: usize> {
    /// One mass for each scalar field.
    pub masses: [Mass; N],
}

impl<const N: usize> QuadraticMultifieldPotential<N> {
    /// Builds a separable quadratic potential from its field masses.
    pub const fn new(masses: [Mass; N]) -> Self {
        Self { masses }
    }
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
{
    fn value(&self, phi: Field<V>) -> EnergyDensity<V> {
        phi.clone() * phi * self.mass * self.mass * 0.5
    }
}

impl<V, const N: usize> MultifieldPotentialValue<V, N> for QuadraticMultifieldPotential<N>
where
    V: From<f64> + Clone + Add<V, Output = V> + Mul<V, Output = V> + Mul<f64, Output = V>,
{
    fn value(&self, phi: [Field<V>; N]) -> EnergyDensity<V> {
        phi.into_iter().zip(self.masses).fold(
            EnergyDensity::new(V::from(0.0)),
            |total, (field, mass)| {
                let term: EnergyDensity<V> = field.clone() * field * mass * mass * 0.5;
                total + term
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use core::ops::Sub;

    use super::*;

    #[derive(Clone, Copy, Debug, PartialEq)]
    struct TestNumber(f64);

    struct StarobinskyTestPotential {
        scale: f64,
        alpha: f64,
    }

    impl<V> MultifieldPotentialValue<V, 1> for StarobinskyTestPotential
    where
        V: From<f64>
            + Clone
            + Sub<V, Output = V>
            + Mul<V, Output = V>
            + Mul<f64, Output = V>
            + crate::Transcendental,
    {
        fn value(&self, phi: [Field<V>; 1]) -> EnergyDensity<V> {
            let exponential = (phi[0].clone().into_value() * -self.alpha).exp();
            let plateau = V::from(1.0) - exponential;

            EnergyDensity::new(plateau.clone() * plateau * self.scale)
        }
    }

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
    fn gradient_comes_from_the_generic_potential_formula() {
        struct CubicPotential;

        impl<V> PotentialValue<V> for CubicPotential
        where
            V: Clone + Mul<V, Output = V>,
        {
            fn value(&self, phi: Field<V>) -> EnergyDensity<V> {
                // The numerical formula is phi^3. This synthetic test focuses
                // on differentiation; its declared output dimension matches
                // the PotentialValue interface.
                let value = phi.into_value();
                EnergyDensity::new(value.clone() * value.clone() * value)
            }
        }

        let gradient = CubicPotential.gradient(Field::new(2.0));

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

    #[test]
    fn multifield_quadratic_value_and_gradient_match_closed_forms() {
        let potential = QuadraticMultifieldPotential::new([Mass::new(2.0), Mass::new(3.0)]);
        let fields = [Field::new(3.0), Field::new(4.0)];

        let value = potential.value(fields).into_value();
        let gradient = potential.gradient(fields).map(|entry| entry.into_value());

        assert_eq!(value, 90.0);
        assert_eq!(gradient, [12.0, 36.0]);
    }

    #[test]
    fn multifield_quadratic_hessian_matches_the_mass_matrix() {
        let potential = QuadraticMultifieldPotential::new([Mass::new(2.0), Mass::new(3.0)]);
        let fields = [Field::new(3.0), Field::new(4.0)];

        let derivatives = potential.derivatives(fields);
        let hessian = derivatives
            .hessian
            .map(|row| row.map(|entry| entry.into_value()));

        assert_eq!(derivatives.value.into_value(), 90.0);
        assert_eq!(
            derivatives.gradient.map(|entry| entry.into_value()),
            [12.0, 36.0]
        );
        assert_eq!(hessian, [[4.0, 0.0], [0.0, 9.0]]);
    }

    #[test]
    fn transcendental_potential_gets_gradient_and_hessian_from_one_formula() {
        let potential = StarobinskyTestPotential {
            scale: 2.0,
            alpha: 0.7,
        };
        let phi = 1.3;
        let exponential = (-potential.alpha * phi).exp();
        let derivatives = potential.derivatives([Field::new(phi)]);
        let expected_value = potential.scale * (1.0 - exponential).powi(2);
        let expected_gradient =
            2.0 * potential.scale * potential.alpha * exponential * (1.0 - exponential);
        let expected_hessian = 2.0
            * potential.scale
            * potential.alpha.powi(2)
            * (-exponential + 2.0 * exponential * exponential);

        assert!((derivatives.value.into_value() - expected_value).abs() < 1.0e-15);
        assert!((derivatives.gradient[0].into_value() - expected_gradient).abs() < 1.0e-15);
        assert!((derivatives.hessian[0][0].into_value() - expected_hessian).abs() < 1.0e-15);
    }
}
