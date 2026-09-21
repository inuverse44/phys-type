//! Fixed-size field-space tensors with compile-time index variance.
//!
//! Tensor variance is independent of physical mass dimension. A vector with
//! an upper field-space index and one with a lower index are therefore
//! different Rust types even when their components have the same dimension.

use core::marker::PhantomData;
use core::ops::{Add, Mul, Neg, Sub};

use crate::dimension::DimMul;
use crate::{M0, Quantity};

/// Marker for a contravariant (upper) field-space index.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Upper;

/// Marker for a covariant (lower) field-space index.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Lower;

/// A fixed-size field-space vector with variance `I` and physical dimension
/// `D`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldVector<I, D, V, const N: usize> {
    components: [Quantity<D, V>; N],
    variance: PhantomData<I>,
}

impl<I, D, V, const N: usize> FieldVector<I, D, V, N> {
    /// Builds a vector from dimension-carrying components.
    pub const fn new(components: [Quantity<D, V>; N]) -> Self {
        Self {
            components,
            variance: PhantomData,
        }
    }

    /// Borrows all components in coordinate order.
    pub const fn components(&self) -> &[Quantity<D, V>; N] {
        &self.components
    }

    /// Removes the variance wrapper and returns all components.
    pub fn into_components(self) -> [Quantity<D, V>; N] {
        self.components
    }
}

impl<I, D, Lhs, Rhs, const N: usize> Add<FieldVector<I, D, Rhs, N>> for FieldVector<I, D, Lhs, N>
where
    Lhs: Add<Rhs>,
{
    type Output = FieldVector<I, D, <Lhs as Add<Rhs>>::Output, N>;

    fn add(self, rhs: FieldVector<I, D, Rhs, N>) -> Self::Output {
        let mut lhs = self.components.into_iter();
        let mut rhs = rhs.components.into_iter();

        FieldVector::new(core::array::from_fn(|_| {
            lhs.next().expect("array length is fixed") + rhs.next().expect("array length is fixed")
        }))
    }
}

impl<I, D, Lhs, Rhs, const N: usize> Sub<FieldVector<I, D, Rhs, N>> for FieldVector<I, D, Lhs, N>
where
    Lhs: Sub<Rhs>,
{
    type Output = FieldVector<I, D, <Lhs as Sub<Rhs>>::Output, N>;

    fn sub(self, rhs: FieldVector<I, D, Rhs, N>) -> Self::Output {
        let mut lhs = self.components.into_iter();
        let mut rhs = rhs.components.into_iter();

        FieldVector::new(core::array::from_fn(|_| {
            lhs.next().expect("array length is fixed") - rhs.next().expect("array length is fixed")
        }))
    }
}

impl<I, D, V, const N: usize> Neg for FieldVector<I, D, V, N>
where
    V: Neg,
{
    type Output = FieldVector<I, D, <V as Neg>::Output, N>;

    fn neg(self) -> Self::Output {
        FieldVector::new(self.components.map(Neg::neg))
    }
}

impl<I, D, V, const N: usize> Mul<f64> for FieldVector<I, D, V, N>
where
    V: Mul<f64>,
{
    type Output = FieldVector<I, D, <V as Mul<f64>>::Output, N>;

    fn mul(self, rhs: f64) -> Self::Output {
        FieldVector::new(self.components.map(|component| component * rhs))
    }
}

impl<Dl, const N: usize> FieldVector<Lower, Dl, f64, N> {
    /// Contracts this covector with a contravariant vector.
    ///
    /// There is deliberately no same-variance overload.
    pub fn contract<Dr>(
        &self,
        rhs: &FieldVector<Upper, Dr, f64, N>,
    ) -> Quantity<<Dl as DimMul<Dr>>::Output, f64>
    where
        Dl: DimMul<Dr>,
    {
        let value = self
            .components
            .iter()
            .zip(rhs.components.iter())
            .map(|(lhs, rhs)| lhs.value() * rhs.value())
            .sum();

        Quantity::new(value)
    }
}

/// A fixed-size rank-two field-space tensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldTensor2<I, J, D, V, const N: usize> {
    components: [[Quantity<D, V>; N]; N],
    variance: PhantomData<(I, J)>,
}

impl<I, J, D, V, const N: usize> FieldTensor2<I, J, D, V, N> {
    /// Builds a rank-two tensor from row-major components.
    pub const fn new(components: [[Quantity<D, V>; N]; N]) -> Self {
        Self {
            components,
            variance: PhantomData,
        }
    }

    /// Borrows the row-major tensor components.
    pub const fn components(&self) -> &[[Quantity<D, V>; N]; N] {
        &self.components
    }

    /// Removes the variance wrapper and returns the components.
    pub fn into_components(self) -> [[Quantity<D, V>; N]; N] {
        self.components
    }
}

/// A metric `G_IJ` with two lower, dimensionless field-space indices.
pub type FieldMetric<const N: usize> = FieldTensor2<Lower, Lower, M0, f64, N>;

/// An inverse metric `G^IJ` with two upper, dimensionless field-space indices.
pub type InverseFieldMetric<const N: usize> = FieldTensor2<Upper, Upper, M0, f64, N>;

impl<const N: usize> FieldTensor2<Lower, Lower, M0, f64, N> {
    /// Returns the Euclidean metric in `N` Cartesian field coordinates.
    pub fn euclidean() -> Self {
        Self::new(core::array::from_fn(|row| {
            core::array::from_fn(|column| Quantity::new(f64::from(row == column)))
        }))
    }

    /// Lowers the index of a contravariant field-space vector.
    pub fn lower<D>(
        &self,
        vector: &FieldVector<Upper, D, f64, N>,
    ) -> FieldVector<Lower, D, f64, N> {
        FieldVector::new(core::array::from_fn(|row| {
            Quantity::new(
                self.components[row]
                    .iter()
                    .zip(vector.components.iter())
                    .map(|(metric, component)| metric.value() * component.value())
                    .sum(),
            )
        }))
    }
}

impl<const N: usize> FieldTensor2<Upper, Upper, M0, f64, N> {
    /// Returns the inverse Euclidean metric in `N` Cartesian field
    /// coordinates.
    pub fn euclidean() -> Self {
        Self::new(core::array::from_fn(|row| {
            core::array::from_fn(|column| Quantity::new(f64::from(row == column)))
        }))
    }

    /// Raises the index of a covariant field-space vector.
    pub fn raise<D>(
        &self,
        covector: &FieldVector<Lower, D, f64, N>,
    ) -> FieldVector<Upper, D, f64, N> {
        FieldVector::new(core::array::from_fn(|row| {
            Quantity::new(
                self.components[row]
                    .iter()
                    .zip(covector.components.iter())
                    .map(|(metric, component)| metric.value() * component.value())
                    .sum(),
            )
        }))
    }
}

/// A field-space vector with one upper index.
pub type ContravariantField<D, V, const N: usize> = FieldVector<Upper, D, V, N>;

/// A field-space vector with one lower index.
pub type CovariantField<D, V, const N: usize> = FieldVector<Lower, D, V, N>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{M1, M2, M3};

    #[test]
    fn euclidean_metric_lowers_and_raises_an_index() {
        let metric = FieldMetric::<2>::euclidean();
        let inverse = InverseFieldMetric::<2>::euclidean();
        let vector =
            ContravariantField::<M1, f64, 2>::new([Quantity::new(3.0), Quantity::new(-4.0)]);

        let covector = metric.lower(&vector);
        let round_trip = inverse.raise(&covector);

        assert_eq!(
            covector.into_components().map(Quantity::into_value),
            [3.0, -4.0]
        );
        assert_eq!(round_trip, vector);
    }

    #[test]
    fn nontrivial_metric_performs_matrix_vector_multiplication() {
        let metric = FieldMetric::<2>::new([
            [Quantity::new(2.0), Quantity::new(0.5)],
            [Quantity::new(0.5), Quantity::new(3.0)],
        ]);
        let vector =
            ContravariantField::<M1, f64, 2>::new([Quantity::new(4.0), Quantity::new(2.0)]);

        let lowered = metric.lower(&vector).into_components();

        assert_eq!(lowered.map(Quantity::into_value), [9.0, 8.0]);
    }

    #[test]
    fn contraction_combines_physical_dimensions() {
        let covector = CovariantField::<M1, f64, 2>::new([Quantity::new(2.0), Quantity::new(3.0)]);
        let vector =
            ContravariantField::<M2, f64, 2>::new([Quantity::new(5.0), Quantity::new(7.0)]);

        let contraction: Quantity<M3, f64> = covector.contract(&vector);

        assert_eq!(contraction.into_value(), 31.0);
    }

    #[test]
    fn vector_arithmetic_preserves_variance() {
        let lhs = ContravariantField::<M1, f64, 2>::new([Quantity::new(1.0), Quantity::new(2.0)]);
        let rhs = ContravariantField::<M1, f64, 2>::new([Quantity::new(3.0), Quantity::new(4.0)]);

        let result: ContravariantField<M1, f64, 2> = (lhs + rhs) * 0.5;

        assert_eq!(
            result.into_components().map(Quantity::into_value),
            [2.0, 3.0]
        );
    }
}
