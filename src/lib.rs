//! Dimension-aware physical quantities for numerical physics.

use core::marker::PhantomData;
use core::ops::{Add, Sub};

use typenum::{N2, P1, P2, P3, P4, P6, P8, Z0};

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

/// A mass dimension whose exponent is encoded in half-power steps.
///
/// For example, `MassPower<P1>` is `M^(1/2)` and `MassPower<P2>` is `M^1`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MassPower<E>(PhantomData<E>);

pub type M0 = MassPower<Z0>;
pub type MHalf = MassPower<P1>;
pub type M1 = MassPower<P2>;
pub type M3Over2 = MassPower<P3>;
pub type M2 = MassPower<P4>;
pub type M3 = MassPower<P6>;
pub type M4 = MassPower<P8>;
pub type MMinus1 = MassPower<N2>;

pub type Dimensionless<V = f64> = Quantity<M0, V>;
pub type Field<V = f64> = Quantity<M1, V>;
pub type FieldVelocity<V = f64> = Quantity<M2, V>;
pub type Hubble<V = f64> = Quantity<M1, V>;
pub type EnergyDensity<V = f64> = Quantity<M4, V>;
pub type Time<V = f64> = Quantity<MMinus1, V>;

#[cfg(test)]
mod tests {
    use super::*;

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
    fn dimension_marker_has_no_runtime_size() {
        assert_eq!(core::mem::size_of::<Field>(), core::mem::size_of::<f64>());
    }
}
