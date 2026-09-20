//! The core value type: a number that carries a compile-time dimension.

use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::dimension::{DimDiv, DimMul, DimSqrt};

/// A physical quantity whose dimension is checked at compile time.
///
/// `D` describes the physical dimension and `V` stores the numerical value.
/// The dimension marker occupies no space at runtime.
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct Quantity<D, V> {
    value: V,
    dimension: PhantomData<D>,
}

impl<D, V> Quantity<D, V> {
    /// Attaches the dimension `D` to a raw numerical value.
    pub const fn new(value: V) -> Self {
        Self {
            value,
            dimension: PhantomData,
        }
    }

    /// Borrows the stored numerical value without removing its dimension.
    pub const fn value(&self) -> &V {
        &self.value
    }

    /// Removes the dimension marker and returns the stored numerical value.
    pub fn into_value(self) -> V {
        self.value
    }
}

impl<D, Lhs, Rhs> Add<Quantity<D, Rhs>> for Quantity<D, Lhs>
where
    Lhs: Add<Rhs>,
{
    type Output = Quantity<D, <Lhs as Add<Rhs>>::Output>;

    fn add(self, rhs: Quantity<D, Rhs>) -> Self::Output {
        Quantity::new(self.value + rhs.value)
    }
}

impl<D, Lhs, Rhs> Sub<Quantity<D, Rhs>> for Quantity<D, Lhs>
where
    Lhs: Sub<Rhs>,
{
    type Output = Quantity<D, <Lhs as Sub<Rhs>>::Output>;

    fn sub(self, rhs: Quantity<D, Rhs>) -> Self::Output {
        Quantity::new(self.value - rhs.value)
    }
}

impl<D, V> Neg for Quantity<D, V>
where
    V: Neg,
{
    type Output = Quantity<D, <V as Neg>::Output>;

    fn neg(self) -> Self::Output {
        Quantity::new(-self.value)
    }
}

impl<Dl, Dr, Vl, Vr> Mul<Quantity<Dr, Vr>> for Quantity<Dl, Vl>
where
    Dl: DimMul<Dr>,
    Vl: Mul<Vr>,
{
    type Output = Quantity<<Dl as DimMul<Dr>>::Output, <Vl as Mul<Vr>>::Output>;

    fn mul(self, rhs: Quantity<Dr, Vr>) -> Self::Output {
        Quantity::new(self.value * rhs.value)
    }
}

impl<Dl, Dr, Vl, Vr> Div<Quantity<Dr, Vr>> for Quantity<Dl, Vl>
where
    Dl: DimDiv<Dr>,
    Vl: Div<Vr>,
{
    type Output = Quantity<<Dl as DimDiv<Dr>>::Output, <Vl as Div<Vr>>::Output>;

    fn div(self, rhs: Quantity<Dr, Vr>) -> Self::Output {
        Quantity::new(self.value / rhs.value)
    }
}

/// Scaling by a raw number, which is treated as a dimensionless coefficient.
impl<D, V> Mul<f64> for Quantity<D, V>
where
    V: Mul<f64>,
{
    type Output = Quantity<D, <V as Mul<f64>>::Output>;

    fn mul(self, rhs: f64) -> Self::Output {
        Quantity::new(self.value * rhs)
    }
}

impl<D, V> Div<f64> for Quantity<D, V>
where
    V: Div<f64>,
{
    type Output = Quantity<D, <V as Div<f64>>::Output>;

    fn div(self, rhs: f64) -> Self::Output {
        Quantity::new(self.value / rhs)
    }
}

impl<D> Mul<Quantity<D, f64>> for f64 {
    type Output = Quantity<D, f64>;

    fn mul(self, rhs: Quantity<D, f64>) -> Self::Output {
        Quantity::new(self * rhs.into_value())
    }
}

impl<D, V> Quantity<D, V>
where
    D: DimSqrt,
{
    /// The positive square root, which halves the mass exponent.
    ///
    /// Available only when the resulting exponent stays on the supported
    /// half-integer grid.
    pub fn sqrt(self) -> Quantity<<D as DimSqrt>::Output, V>
    where
        V: Sqrt,
    {
        Quantity::new(self.value.sqrt_value())
    }
}

/// The numerical part of a square root, kept separate so that later value
/// types such as automatic-differentiation numbers can opt in.
pub trait Sqrt {
    /// Returns the positive square root of the stored number.
    fn sqrt_value(self) -> Self;
}

impl Sqrt for f64 {
    fn sqrt_value(self) -> Self {
        self.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        Dimensionless, EnergyDensity, Field, FieldVelocity, Hubble, M1, M2, M3, M4, MHalf, Mass,
        Quantity, Time,
    };

    #[test]
    fn creates_and_extracts_a_value() {
        let phi = Field::new(0.17);

        assert_eq!(*phi.value(), 0.17);
        assert_eq!(phi.into_value(), 0.17);
    }

    #[test]
    fn adds_quantities_with_the_same_dimension() {
        let lhs = Field::new(1.25);
        let rhs = Field::new(2.75);

        assert_eq!((lhs + rhs).into_value(), 4.0);
    }

    #[test]
    fn subtracts_quantities_with_the_same_dimension() {
        let lhs = FieldVelocity::new(3.5);
        let rhs = FieldVelocity::new(1.25);

        assert_eq!((lhs - rhs).into_value(), 2.25);
    }

    #[test]
    fn negates_a_quantity() {
        let gradient: Quantity<M3, f64> = Quantity::new(2.5);

        assert_eq!((-gradient).into_value(), -2.5);
    }

    #[test]
    fn dimension_marker_has_no_runtime_size() {
        assert_eq!(core::mem::size_of::<Field>(), core::mem::size_of::<f64>());
    }

    #[test]
    fn multiplication_adds_mass_exponents() {
        let phi = Field::new(2.0);
        let mass = Mass::new(3.0);

        let product: Quantity<M2, f64> = phi * mass;

        assert_eq!(product.into_value(), 6.0);
    }

    #[test]
    fn multiplying_two_half_powers_gives_one_mass_power() {
        let root: Quantity<MHalf, f64> = Quantity::new(4.0);

        let product: Quantity<M1, f64> = root * root;

        assert_eq!(product.into_value(), 16.0);
    }

    #[test]
    fn division_subtracts_mass_exponents() {
        let density = EnergyDensity::new(12.0);
        let squared_mass: Quantity<M2, f64> = Quantity::new(4.0);

        let quotient: Quantity<M2, f64> = density / squared_mass;

        assert_eq!(quotient.into_value(), 3.0);
    }

    #[test]
    fn dividing_equal_dimensions_gives_a_dimensionless_quantity() {
        let lhs = Field::new(6.0);
        let rhs = Field::new(4.0);

        let ratio: Dimensionless = lhs / rhs;

        assert_eq!(ratio.into_value(), 1.5);
    }

    #[test]
    fn field_velocity_over_time_has_the_expected_dimension() {
        let pi = FieldVelocity::new(9.0);
        let dt = Time::new(3.0);

        let rate: Quantity<M3, f64> = pi / dt;

        assert_eq!(rate.into_value(), 3.0);
    }

    #[test]
    fn scales_by_a_raw_number_from_either_side() {
        let phi = Field::new(2.0);

        assert_eq!((phi * 3.0).into_value(), 6.0);
        assert_eq!((3.0 * phi).into_value(), 6.0);
        assert_eq!((phi / 4.0).into_value(), 0.5);
    }

    #[test]
    fn square_root_halves_the_mass_exponent() {
        let squared: Quantity<M2, f64> = Quantity::new(9.0);

        let root: Hubble = squared.sqrt();

        assert_eq!(root.into_value(), 3.0);
    }

    #[test]
    fn square_root_of_an_energy_density_is_a_squared_mass() {
        let density = EnergyDensity::new(16.0);

        let root: Quantity<M2, f64> = density.sqrt();

        assert_eq!(root.into_value(), 4.0);
    }

    #[test]
    fn square_root_reaches_the_half_integer_grid() {
        let field = Field::new(4.0);

        let root: Quantity<MHalf, f64> = field.sqrt();

        assert_eq!(root.into_value(), 2.0);
    }

    #[test]
    fn friedmann_style_expression_type_checks() {
        let pi = FieldVelocity::new(0.5);
        let potential = EnergyDensity::new(1.0);
        let m_pl = Mass::new(1.0);

        let total: Quantity<M4, f64> = pi * pi * 0.5 + potential;
        let hubble: Hubble = (total / (m_pl * m_pl * 3.0)).sqrt();

        assert!(hubble.into_value() > 0.0);
    }
}
