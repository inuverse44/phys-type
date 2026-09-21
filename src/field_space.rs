//! Fixed-size field-space tensors with compile-time index variance.
//!
//! Tensor variance is independent of physical mass dimension. A vector with
//! an upper field-space index and one with a lower index are therefore
//! different Rust types even when their components have the same dimension.

use core::marker::PhantomData;
use core::ops::{Add, Mul, Neg, Sub};

use crate::dimension::DimMul;
use crate::{Field, FieldVelocity, M0, M2, M3, MMinus1, Quantity};

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

/// Christoffel symbols `Gamma^I_JK` for scalar fields of mass dimension one.
///
/// A derivative with respect to a field contributes dimension `M^-1`, so the
/// connection components carry that dimension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChristoffelSymbols<const N: usize> {
    components: [[[Quantity<MMinus1, f64>; N]; N]; N],
}

impl<const N: usize> ChristoffelSymbols<N> {
    /// Builds connection coefficients in `upper, lower, lower` index order.
    pub const fn new(components: [[[Quantity<MMinus1, f64>; N]; N]; N]) -> Self {
        Self { components }
    }

    /// Returns a vanishing connection.
    pub fn zero() -> Self {
        Self::new([[[Quantity::new(0.0); N]; N]; N])
    }

    /// Borrows the connection coefficients.
    pub const fn components(&self) -> &[[[Quantity<MMinus1, f64>; N]; N]; N] {
        &self.components
    }

    /// Computes `Gamma^I_JK velocity^J velocity^K`.
    pub fn quadratic_velocity(
        &self,
        velocity: &ContravariantField<M2, f64, N>,
    ) -> ContravariantField<M3, f64, N> {
        ContravariantField::new(core::array::from_fn(|upper| {
            Quantity::new(
                (0..N)
                    .flat_map(|first| (0..N).map(move |second| (first, second)))
                    .map(|(first, second)| {
                        self.components[upper][first][second].value()
                            * velocity.components[first].value()
                            * velocity.components[second].value()
                    })
                    .sum(),
            )
        }))
    }
}

/// Supplies field-space geometry at a background position.
pub trait FieldSpaceGeometry<const N: usize> {
    /// The covariant metric `G_IJ`.
    fn metric(&self, position: &[Field; N]) -> FieldMetric<N>;

    /// The inverse metric `G^IJ`.
    fn inverse_metric(&self, position: &[Field; N]) -> InverseFieldMetric<N>;

    /// The Levi-Civita connection `Gamma^I_JK`.
    fn christoffel(&self, position: &[Field; N]) -> ChristoffelSymbols<N>;
}

/// Euclidean field space in Cartesian coordinates.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FlatFieldSpace;

impl<const N: usize> FieldSpaceGeometry<N> for FlatFieldSpace {
    fn metric(&self, _position: &[Field; N]) -> FieldMetric<N> {
        FieldMetric::euclidean()
    }

    fn inverse_metric(&self, _position: &[Field; N]) -> InverseFieldMetric<N> {
        InverseFieldMetric::euclidean()
    }

    fn christoffel(&self, _position: &[Field; N]) -> ChristoffelSymbols<N> {
        ChristoffelSymbols::zero()
    }
}

/// Why an orthonormal two-field kinematic basis could not be constructed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KinematicBasisError {
    /// The background is instantaneously stationary, so its tangent is
    /// undefined.
    ZeroSpeed,
    /// The supplied metric is not symmetric to floating-point tolerance.
    NonSymmetricMetric,
    /// The supplied two-dimensional metric is not positive definite.
    NonPositiveMetric,
}

/// Instantaneous tangent and normal directions for a two-field trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KinematicBasis2 {
    /// Unit vector tangent to the background velocity.
    pub tangent: ContravariantField<M0, f64, 2>,
    /// Unit vector normal to the background velocity.
    pub normal: ContravariantField<M0, f64, 2>,
    /// Covariant magnitude `dot sigma` of the background velocity.
    pub speed: FieldVelocity,
}

impl KinematicBasis2 {
    /// Constructs an oriented orthonormal basis.
    ///
    /// If `previous_normal` is supplied, the sign is selected to give a
    /// positive overlap with it. Otherwise the positive coordinate
    /// orientation is used.
    pub fn from_velocity(
        metric: &FieldMetric<2>,
        velocity: &ContravariantField<M2, f64, 2>,
        previous_normal: Option<&ContravariantField<M0, f64, 2>>,
    ) -> Result<Self, KinematicBasisError> {
        let components = metric.components();
        let g00 = *components[0][0].value();
        let g01 = *components[0][1].value();
        let g10 = *components[1][0].value();
        let g11 = *components[1][1].value();
        let scale = g01.abs().max(g10.abs()).max(1.0);

        if (g01 - g10).abs() > 1.0e-12 * scale {
            return Err(KinematicBasisError::NonSymmetricMetric);
        }

        let determinant = g00 * g11 - g01 * g10;
        if g00 <= 0.0 || determinant <= 0.0 {
            return Err(KinematicBasisError::NonPositiveMetric);
        }

        let lowered_velocity = metric.lower(velocity);
        let speed_squared = lowered_velocity.contract(velocity).into_value();
        if speed_squared <= 0.0 || !speed_squared.is_finite() {
            return Err(KinematicBasisError::ZeroSpeed);
        }

        let speed = speed_squared.sqrt();
        let tangent = ContravariantField::new([
            Quantity::new(velocity.components[0].value() / speed),
            Quantity::new(velocity.components[1].value() / speed),
        ]);
        let tangent_lower = metric.lower(&tangent);
        let root_determinant = determinant.sqrt();
        let mut normal = ContravariantField::new([
            Quantity::new(-tangent_lower.components[1].value() / root_determinant),
            Quantity::new(tangent_lower.components[0].value() / root_determinant),
        ]);

        if let Some(previous) = previous_normal
            && metric.lower(&normal).contract(previous).into_value() < 0.0
        {
            normal = normal * -1.0;
        }

        Ok(Self {
            tangent,
            normal,
            speed: FieldVelocity::new(speed),
        })
    }
}

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

    #[test]
    fn christoffel_symbols_contract_both_lower_indices() {
        let mut components = [[[Quantity::new(0.0); 2]; 2]; 2];
        components[0][0][0] = Quantity::new(0.5);
        let connection = ChristoffelSymbols::new(components);
        let velocity =
            ContravariantField::<M2, f64, 2>::new([Quantity::new(2.0), Quantity::new(3.0)]);

        let result = connection
            .quadratic_velocity(&velocity)
            .into_components()
            .map(Quantity::into_value);

        assert_eq!(result, [2.0, 0.0]);
    }

    #[test]
    fn two_field_basis_is_oriented_and_orthonormal() {
        let metric = FieldMetric::<2>::euclidean();
        let velocity =
            ContravariantField::<M2, f64, 2>::new([Quantity::new(3.0), Quantity::new(4.0)]);

        let basis = KinematicBasis2::from_velocity(&metric, &velocity, None).unwrap();
        let tangent_norm = metric.lower(&basis.tangent).contract(&basis.tangent);
        let normal_norm = metric.lower(&basis.normal).contract(&basis.normal);
        let overlap = metric.lower(&basis.tangent).contract(&basis.normal);

        assert_eq!(basis.speed.into_value(), 5.0);
        assert!((tangent_norm.into_value() - 1.0).abs() < 1.0e-15);
        assert!((normal_norm.into_value() - 1.0).abs() < 1.0e-15);
        assert!(overlap.into_value().abs() < 1.0e-15);
        assert_eq!(
            basis.normal.into_components().map(Quantity::into_value),
            [-0.8, 0.6]
        );
    }

    #[test]
    fn previous_normal_prevents_an_orientation_flip() {
        let metric = FieldMetric::<2>::euclidean();
        let velocity =
            ContravariantField::<M2, f64, 2>::new([Quantity::new(3.0), Quantity::new(4.0)]);
        let previous =
            ContravariantField::<M0, f64, 2>::new([Quantity::new(0.8), Quantity::new(-0.6)]);

        let basis = KinematicBasis2::from_velocity(&metric, &velocity, Some(&previous)).unwrap();

        assert_eq!(basis.normal, previous);
    }
}
