//! Fixed-size field-space tensors with compile-time index variance.
//!
//! Tensor variance is independent of physical mass dimension. A vector with
//! an upper field-space index and one with a lower index are therefore
//! different Rust types even when their components have the same dimension.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Mul, Neg, Sub};

use crate::dimension::DimMul;
use crate::{Field, FieldVelocity, M0, M1, M2, M3, MMinus1, MMinus2, MMinus3, Quantity};

/// Dense components of a rank-three field-space tensor.
pub type Tensor3Components<D, V, const N: usize> = [[[Quantity<D, V>; N]; N]; N];

/// Dense components of a rank-four field-space tensor.
pub type Tensor4Components<D, V, const N: usize> = [[[[Quantity<D, V>; N]; N]; N]; N];

/// Dense components of a rank-five field-space tensor.
pub type Tensor5Components<D, V, const N: usize> = [[[[[Quantity<D, V>; N]; N]; N]; N]; N];

/// Coordinate components of `partial_M Gamma^I_JK`.
pub type ChristoffelDerivativeComponents<const N: usize> = Tensor4Components<MMinus2, f64, N>;

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

/// A fixed-size rank-three field-space tensor.
///
/// This type is only for genuine tensors. Coordinate-dependent objects such
/// as [`ChristoffelSymbols`] intentionally use a separate type.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldTensor3<I, J, K, D, V, const N: usize> {
    components: Tensor3Components<D, V, N>,
    variance: PhantomData<(I, J, K)>,
}

impl<I, J, K, D, V, const N: usize> FieldTensor3<I, J, K, D, V, N> {
    /// Builds a rank-three tensor in first-second-third index order.
    pub const fn new(components: Tensor3Components<D, V, N>) -> Self {
        Self {
            components,
            variance: PhantomData,
        }
    }

    /// Borrows all components.
    pub const fn components(&self) -> &Tensor3Components<D, V, N> {
        &self.components
    }

    /// Removes the variance wrapper and returns all components.
    pub fn into_components(self) -> Tensor3Components<D, V, N> {
        self.components
    }
}

impl<I, J, K, D: Copy, const N: usize> Add for FieldTensor3<I, J, K, D, f64, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    self.components[first][second][third] + rhs.components[first][second][third]
                })
            })
        }))
    }
}

impl<I, J, K, D: Copy, const N: usize> Sub for FieldTensor3<I, J, K, D, f64, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    self.components[first][second][third] - rhs.components[first][second][third]
                })
            })
        }))
    }
}

impl<I, J, K, D, const N: usize> Neg for FieldTensor3<I, J, K, D, f64, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(
            self.components
                .map(|matrix| matrix.map(|row| row.map(Neg::neg))),
        )
    }
}

impl<I, J, K, D, const N: usize> Mul<f64> for FieldTensor3<I, J, K, D, f64, N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(
            self.components
                .map(|matrix| matrix.map(|row| row.map(|component| component * rhs))),
        )
    }
}

/// A fixed-size rank-four field-space tensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldTensor4<I, J, K, L, D, V, const N: usize> {
    components: Tensor4Components<D, V, N>,
    variance: PhantomData<(I, J, K, L)>,
}

/// A fixed-size rank-five field-space tensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldTensor5<I, J, K, L, M, D, V, const N: usize> {
    components: Tensor5Components<D, V, N>,
    variance: PhantomData<(I, J, K, L, M)>,
}

impl<I, J, K, L, M, D, V, const N: usize> FieldTensor5<I, J, K, L, M, D, V, N> {
    /// Builds a rank-five tensor in first-through-fifth index order.
    pub const fn new(components: Tensor5Components<D, V, N>) -> Self {
        Self {
            components,
            variance: PhantomData,
        }
    }

    /// Borrows all components.
    pub const fn components(&self) -> &Tensor5Components<D, V, N> {
        &self.components
    }

    /// Removes the variance wrapper and returns all components.
    pub fn into_components(self) -> Tensor5Components<D, V, N> {
        self.components
    }
}

impl<I, J, K, L, D, V, const N: usize> FieldTensor4<I, J, K, L, D, V, N> {
    /// Builds a rank-four tensor in first-second-third-fourth index order.
    pub const fn new(components: Tensor4Components<D, V, N>) -> Self {
        Self {
            components,
            variance: PhantomData,
        }
    }

    /// Borrows all components.
    pub const fn components(&self) -> &Tensor4Components<D, V, N> {
        &self.components
    }

    /// Removes the variance wrapper and returns all components.
    pub fn into_components(self) -> Tensor4Components<D, V, N> {
        self.components
    }
}

impl<I, J, K, L, D: Copy, const N: usize> Add for FieldTensor4<I, J, K, L, D, f64, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    core::array::from_fn(|fourth| {
                        self.components[first][second][third][fourth]
                            + rhs.components[first][second][third][fourth]
                    })
                })
            })
        }))
    }
}

impl<I, J, K, L, D: Copy, const N: usize> Sub for FieldTensor4<I, J, K, L, D, f64, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    core::array::from_fn(|fourth| {
                        self.components[first][second][third][fourth]
                            - rhs.components[first][second][third][fourth]
                    })
                })
            })
        }))
    }
}

impl<I, J, K, L, D, const N: usize> Neg for FieldTensor4<I, J, K, L, D, f64, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(
            self.components
                .map(|cube| cube.map(|matrix| matrix.map(|row| row.map(Neg::neg)))),
        )
    }
}

impl<I, J, K, L, D, const N: usize> Mul<f64> for FieldTensor4<I, J, K, L, D, f64, N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(
            self.components.map(|cube| {
                cube.map(|matrix| matrix.map(|row| row.map(|component| component * rhs)))
            }),
        )
    }
}

/// Lowers the first upper index of a tensor using a metric.
pub trait LowerFirst<Tensor> {
    /// Tensor type after lowering the first index.
    type Output;

    /// Applies `G_AI` to the first index.
    fn lower_first(&self, tensor: &Tensor) -> Self::Output;
}

/// Raises the first lower index of a tensor using an inverse metric.
pub trait RaiseFirst<Tensor> {
    /// Tensor type after raising the first index.
    type Output;

    /// Applies `G^AI` to the first index.
    fn raise_first(&self, tensor: &Tensor) -> Self::Output;
}

impl<J, K, D, const N: usize> LowerFirst<FieldTensor3<Upper, J, K, D, f64, N>> for FieldMetric<N> {
    type Output = FieldTensor3<Lower, J, K, D, f64, N>;

    fn lower_first(&self, tensor: &FieldTensor3<Upper, J, K, D, f64, N>) -> Self::Output {
        FieldTensor3::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    Quantity::new(
                        (0..N)
                            .map(|contracted| {
                                self.components[first][contracted].value()
                                    * tensor.components[contracted][second][third].value()
                            })
                            .sum(),
                    )
                })
            })
        }))
    }
}

impl<J, K, D, const N: usize> RaiseFirst<FieldTensor3<Lower, J, K, D, f64, N>>
    for InverseFieldMetric<N>
{
    type Output = FieldTensor3<Upper, J, K, D, f64, N>;

    fn raise_first(&self, tensor: &FieldTensor3<Lower, J, K, D, f64, N>) -> Self::Output {
        FieldTensor3::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    Quantity::new(
                        (0..N)
                            .map(|contracted| {
                                self.components[first][contracted].value()
                                    * tensor.components[contracted][second][third].value()
                            })
                            .sum(),
                    )
                })
            })
        }))
    }
}

impl<J, K, L, D, const N: usize> LowerFirst<FieldTensor4<Upper, J, K, L, D, f64, N>>
    for FieldMetric<N>
{
    type Output = FieldTensor4<Lower, J, K, L, D, f64, N>;

    fn lower_first(&self, tensor: &FieldTensor4<Upper, J, K, L, D, f64, N>) -> Self::Output {
        FieldTensor4::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    core::array::from_fn(|fourth| {
                        Quantity::new(
                            (0..N)
                                .map(|contracted| {
                                    self.components[first][contracted].value()
                                        * tensor.components[contracted][second][third][fourth]
                                            .value()
                                })
                                .sum(),
                        )
                    })
                })
            })
        }))
    }
}

impl<J, K, L, D, const N: usize> RaiseFirst<FieldTensor4<Lower, J, K, L, D, f64, N>>
    for InverseFieldMetric<N>
{
    type Output = FieldTensor4<Upper, J, K, L, D, f64, N>;

    fn raise_first(&self, tensor: &FieldTensor4<Lower, J, K, L, D, f64, N>) -> Self::Output {
        FieldTensor4::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    core::array::from_fn(|fourth| {
                        Quantity::new(
                            (0..N)
                                .map(|contracted| {
                                    self.components[first][contracted].value()
                                        * tensor.components[contracted][second][third][fourth]
                                            .value()
                                })
                                .sum(),
                        )
                    })
                })
            })
        }))
    }
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

/// A covariant potential Hessian `V_;IJ`.
pub type CovariantPotentialHessian<const N: usize> = FieldTensor2<Lower, Lower, M2, f64, N>;

/// The fully symmetrized covariant third potential derivative `V_;(IJK)`.
pub type CovariantPotentialThirdDerivative<const N: usize> =
    FieldTensor3<Lower, Lower, Lower, M1, f64, N>;

/// Field-space Riemann curvature in the convention `R^I_JKL`.
pub type RiemannTensor<const N: usize> = FieldTensor4<Upper, Lower, Lower, Lower, MMinus2, f64, N>;

/// Field-space Riemann curvature with all four indices lowered.
pub type LoweredRiemannTensor<const N: usize> =
    FieldTensor4<Lower, Lower, Lower, Lower, MMinus2, f64, N>;

/// The covariant curvature derivative `nabla_M R^I_JKL`.
pub type CovariantRiemannDerivative<const N: usize> =
    FieldTensor5<Upper, Lower, Lower, Lower, Lower, MMinus3, f64, N>;

impl<const N: usize> FieldTensor5<Upper, Lower, Lower, Lower, Lower, MMinus3, f64, N> {
    /// Returns an identically zero covariant curvature derivative.
    pub fn zero() -> Self {
        Self::new([[[[[Quantity::new(0.0); N]; N]; N]; N]; N])
    }
}

impl<const N: usize> FieldTensor4<Upper, Lower, Lower, Lower, MMinus2, f64, N> {
    /// Returns identically zero curvature.
    pub fn zero() -> Self {
        Self::new([[[[Quantity::new(0.0); N]; N]; N]; N])
    }

    /// Forms `R^I_LMJ velocity^L velocity^M`.
    ///
    /// The result has mixed indices `(I, J)` and mass dimension `M^2`, ready
    /// for the linear perturbation mass matrix.
    pub fn contract_velocities(
        &self,
        velocity: &ContravariantField<M2, f64, N>,
    ) -> FieldTensor2<Upper, Lower, M2, f64, N> {
        FieldTensor2::new(core::array::from_fn(|upper| {
            core::array::from_fn(|lower| {
                Quantity::new(
                    (0..N)
                        .flat_map(|first| (0..N).map(move |second| (first, second)))
                        .map(|(first, second)| {
                            self.components[upper][first][second][lower].value()
                                * velocity.components[first].value()
                                * velocity.components[second].value()
                        })
                        .sum(),
                )
            })
        }))
    }
}

/// Christoffel symbols `Gamma^I_JK` for scalar fields of mass dimension one.
///
/// A derivative with respect to a field contributes dimension `M^-1`, so the
/// connection components carry that dimension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChristoffelSymbols<const N: usize> {
    components: [[[Quantity<MMinus1, f64>; N]; N]; N],
}

/// Coordinate derivatives `partial_M Gamma^I_JK` of Christoffel symbols.
///
/// Like [`ChristoffelSymbols`], this is intentionally not represented as a
/// tensor. Components use upper, first lower, second lower, derivative-index
/// order and carry mass dimension `M^-2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChristoffelDerivatives<const N: usize> {
    components: ChristoffelDerivativeComponents<N>,
}

impl<const N: usize> ChristoffelDerivatives<N> {
    /// Builds coordinate connection derivatives.
    pub const fn new(components: ChristoffelDerivativeComponents<N>) -> Self {
        Self { components }
    }

    /// Returns identically zero coordinate connection derivatives.
    pub fn zero() -> Self {
        Self::new([[[[Quantity::new(0.0); N]; N]; N]; N])
    }

    /// Borrows components in upper, lower, lower, derivative order.
    pub const fn components(&self) -> &ChristoffelDerivativeComponents<N> {
        &self.components
    }
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

/// A position where a field-space geometry, or the representation of a
/// theory it encodes, is not valid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeometryDomainError {
    /// The condition that failed.
    pub reason: &'static str,
    /// The offending value, such as `Omega^2` or the smallest Cholesky pivot.
    pub value: f64,
}

impl fmt::Display for GeometryDomainError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{} (value {:e})", self.reason, self.value)
    }
}

impl std::error::Error for GeometryDomainError {}

/// Supplies field-space geometry at a background position.
pub trait FieldSpaceGeometry<const N: usize> {
    /// The covariant metric `G_IJ`.
    fn metric(&self, position: &[Field; N]) -> FieldMetric<N>;

    /// The inverse metric `G^IJ`.
    fn inverse_metric(&self, position: &[Field; N]) -> InverseFieldMetric<N>;

    /// The Levi-Civita connection `Gamma^I_JK`.
    fn christoffel(&self, position: &[Field; N]) -> ChristoffelSymbols<N>;

    /// Checks that the geometry is valid at `position`. The background
    /// survey calls it on every accepted step. The default accepts every
    /// position; geometries derived from a theory check their validity
    /// conditions, such as a positive conformal factor.
    fn check_domain(&self, _position: &[Field; N]) -> Result<(), GeometryDomainError> {
        Ok(())
    }
}

/// A field-space geometry that supplies analytic Riemann curvature.
///
/// Curvature is a separate capability so background-only geometries do not
/// need to manufacture it. Implementations must use the crate convention
/// `[nabla_K, nabla_L] A^I = R^I_JKL A^J`.
pub trait CurvedFieldSpaceGeometry<const N: usize>: FieldSpaceGeometry<N> {
    /// Returns `R^I_JKL` at the supplied background position.
    fn riemann(&self, position: &[Field; N]) -> RiemannTensor<N>;
}

/// A geometry supplying the coordinate connection derivatives needed to
/// covariantize third derivatives of scalar potentials.
pub trait ThirdOrderFieldSpaceGeometry<const N: usize>: FieldSpaceGeometry<N> {
    /// Returns `partial_M Gamma^I_JK` at the supplied position.
    fn christoffel_derivative(&self, position: &[Field; N]) -> ChristoffelDerivatives<N>;
}

/// A curved geometry supplying the rank-five derivative needed by the cubic
/// transport kernel.
pub trait CurvatureDerivativeGeometry<const N: usize>: CurvedFieldSpaceGeometry<N> {
    /// Returns `nabla_M R^I_JKL` at the supplied position.
    fn covariant_riemann_derivative(&self, position: &[Field; N]) -> CovariantRiemannDerivative<N>;
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

impl<const N: usize> CurvedFieldSpaceGeometry<N> for FlatFieldSpace {
    fn riemann(&self, _position: &[Field; N]) -> RiemannTensor<N> {
        RiemannTensor::zero()
    }
}

impl<const N: usize> ThirdOrderFieldSpaceGeometry<N> for FlatFieldSpace {
    fn christoffel_derivative(&self, _position: &[Field; N]) -> ChristoffelDerivatives<N> {
        ChristoffelDerivatives::zero()
    }
}

impl<const N: usize> CurvatureDerivativeGeometry<N> for FlatFieldSpace {
    fn covariant_riemann_derivative(
        &self,
        _position: &[Field; N],
    ) -> CovariantRiemannDerivative<N> {
        CovariantRiemannDerivative::zero()
    }
}

/// Euclidean field space in polar coordinates `(phi^1, phi^2) = (rho, theta)`.
///
/// The line element is
///
/// ```text
/// ds^2 = (drho)^2 + rho^2 (dtheta)^2,
/// ```
///
/// so the Riemann tensor vanishes while the connection does not. Circular
/// orbits at fixed `rho` are the simplest trajectories of constant turn rate,
/// with `rho` normal and `theta` tangent to the motion. The angle `theta` is
/// measured in units of the reduced Planck mass, like every field coordinate.
///
/// The coordinates degenerate at `rho = 0`, where `G_22` vanishes and the
/// inverse metric diverges. Keep trajectories away from the origin.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PolarFieldSpace;

impl FieldSpaceGeometry<2> for PolarFieldSpace {
    fn metric(&self, position: &[Field; 2]) -> FieldMetric<2> {
        let rho = position[0].into_value();

        FieldMetric::new([
            [Quantity::new(1.0), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(rho * rho)],
        ])
    }

    fn inverse_metric(&self, position: &[Field; 2]) -> InverseFieldMetric<2> {
        let rho = position[0].into_value();

        InverseFieldMetric::new([
            [Quantity::new(1.0), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(1.0 / (rho * rho))],
        ])
    }

    fn christoffel(&self, position: &[Field; 2]) -> ChristoffelSymbols<2> {
        let rho = position[0].into_value();
        let mut components = [[[Quantity::<MMinus1, f64>::new(0.0); 2]; 2]; 2];

        // Gamma^1_22 = -rho
        components[0][1][1] = Quantity::new(-rho);
        // Gamma^2_12 = Gamma^2_21 = 1 / rho
        let connection = Quantity::new(1.0 / rho);
        components[1][0][1] = connection;
        components[1][1][0] = connection;

        ChristoffelSymbols::new(components)
    }
}

impl CurvedFieldSpaceGeometry<2> for PolarFieldSpace {
    fn riemann(&self, _position: &[Field; 2]) -> RiemannTensor<2> {
        RiemannTensor::zero()
    }
}

impl ThirdOrderFieldSpaceGeometry<2> for PolarFieldSpace {
    fn christoffel_derivative(&self, position: &[Field; 2]) -> ChristoffelDerivatives<2> {
        let rho = position[0].into_value();
        let mut components = [[[[Quantity::<MMinus2, f64>::new(0.0); 2]; 2]; 2]; 2];

        // partial_1 Gamma^1_22 = -1
        components[0][1][1][0] = Quantity::new(-1.0);
        // partial_1 Gamma^2_12 = partial_1 Gamma^2_21 = -1 / rho^2
        let angular = Quantity::new(-1.0 / (rho * rho));
        components[1][0][1][0] = angular;
        components[1][1][0][0] = angular;

        ChristoffelDerivatives::new(components)
    }
}

impl CurvatureDerivativeGeometry<2> for PolarFieldSpace {
    fn covariant_riemann_derivative(
        &self,
        _position: &[Field; 2],
    ) -> CovariantRiemannDerivative<2> {
        CovariantRiemannDerivative::zero()
    }
}

/// A maximally symmetric field space of constant negative curvature.
///
/// In coordinates `(phi^1, phi^2)` the line element is
///
/// ```text
/// ds^2 = (dphi^1)^2 + sinh^2(phi^1 / L) (dphi^2)^2,
/// ```
///
/// whose Gaussian curvature is `-1/L^2`. This is the two-field geometry that
/// appears in alpha-attractor and hyperbolic-inflation models, and it is the
/// smallest example that exercises a non-vanishing connection and a
/// non-vanishing Riemann tensor at once.
///
/// The coordinate `phi^1` is a geodesic radius, so the metric degenerates at
/// `phi^1 = 0` where `G_22` vanishes and the inverse metric diverges. Keep
/// trajectories away from the origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HyperbolicFieldSpace {
    /// The curvature radius `L`; the Gaussian curvature is `-1/L^2`.
    curvature_radius: Field,
}

/// An invalid parameter for [`HyperbolicFieldSpace`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InvalidCurvatureRadius;

impl fmt::Display for InvalidCurvatureRadius {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("the hyperbolic curvature radius must be positive and finite")
    }
}

impl std::error::Error for InvalidCurvatureRadius {}

impl HyperbolicFieldSpace {
    /// Builds a field space of positive, finite curvature radius `L`.
    pub fn new(curvature_radius: Field) -> Result<Self, InvalidCurvatureRadius> {
        let value = curvature_radius.into_value();
        if !value.is_finite() || value <= 0.0 {
            return Err(InvalidCurvatureRadius);
        }

        Ok(Self { curvature_radius })
    }

    /// Returns the curvature radius `L`.
    pub const fn curvature_radius(&self) -> Field {
        self.curvature_radius
    }

    /// `sinh(phi^1 / L)` and `cosh(phi^1 / L)` at a position.
    fn hyperbolics(&self, position: &[Field; 2]) -> (f64, f64) {
        let scaled = (position[0] / self.curvature_radius).into_value();

        (scaled.sinh(), scaled.cosh())
    }
}

impl FieldSpaceGeometry<2> for HyperbolicFieldSpace {
    fn metric(&self, position: &[Field; 2]) -> FieldMetric<2> {
        let (sinh, _) = self.hyperbolics(position);

        FieldMetric::new([
            [Quantity::new(1.0), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(sinh * sinh)],
        ])
    }

    fn inverse_metric(&self, position: &[Field; 2]) -> InverseFieldMetric<2> {
        let (sinh, _) = self.hyperbolics(position);

        InverseFieldMetric::new([
            [Quantity::new(1.0), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(1.0 / (sinh * sinh))],
        ])
    }

    fn christoffel(&self, position: &[Field; 2]) -> ChristoffelSymbols<2> {
        let (sinh, cosh) = self.hyperbolics(position);
        let radius = self.curvature_radius.into_value();
        let mut components = [[[Quantity::<MMinus1, f64>::new(0.0); 2]; 2]; 2];

        // Gamma^1_22 = -sinh cosh / L
        components[0][1][1] = Quantity::new(-sinh * cosh / radius);
        // Gamma^2_12 = Gamma^2_21 = coth / L
        let connection = Quantity::new(cosh / (sinh * radius));
        components[1][0][1] = connection;
        components[1][1][0] = connection;

        ChristoffelSymbols::new(components)
    }
}

impl CurvedFieldSpaceGeometry<2> for HyperbolicFieldSpace {
    /// `R^I_JKL = K (delta^I_K G_JL - delta^I_L G_JK)` with `K = -1/L^2`,
    /// which is the maximally symmetric form in the crate convention
    /// `[nabla_K, nabla_L] A^I = R^I_JKL A^J`.
    fn riemann(&self, position: &[Field; 2]) -> RiemannTensor<2> {
        let metric = self.metric(position);
        let radius = self.curvature_radius.into_value();
        let curvature = -1.0 / (radius * radius);
        let components = metric.components();

        RiemannTensor::new(core::array::from_fn(|upper| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    core::array::from_fn(|fourth| {
                        let first =
                            f64::from(upper == third) * components[second][fourth].into_value();
                        let second_term =
                            f64::from(upper == fourth) * components[second][third].into_value();

                        Quantity::new(curvature * (first - second_term))
                    })
                })
            })
        }))
    }
}

impl ThirdOrderFieldSpaceGeometry<2> for HyperbolicFieldSpace {
    fn christoffel_derivative(&self, position: &[Field; 2]) -> ChristoffelDerivatives<2> {
        let (sinh, cosh) = self.hyperbolics(position);
        let radius_squared = self.curvature_radius.into_value().powi(2);
        let mut components = [[[[Quantity::<MMinus2, f64>::new(0.0); 2]; 2]; 2]; 2];

        // partial_1 Gamma^1_22 = -(cosh^2 + sinh^2) / L^2
        components[0][1][1][0] = Quantity::new(-(cosh * cosh + sinh * sinh) / radius_squared);
        // partial_1 Gamma^2_12 = partial_1 Gamma^2_21 = -csch^2 / L^2
        let angular = Quantity::new(-1.0 / (sinh * sinh * radius_squared));
        components[1][0][1][0] = angular;
        components[1][1][0][0] = angular;

        ChristoffelDerivatives::new(components)
    }
}

impl CurvatureDerivativeGeometry<2> for HyperbolicFieldSpace {
    fn covariant_riemann_derivative(
        &self,
        _position: &[Field; 2],
    ) -> CovariantRiemannDerivative<2> {
        // The hyperbolic metric is maximally symmetric with constant
        // sectional curvature, so metric compatibility gives nabla R = 0.
        CovariantRiemannDerivative::zero()
    }
}

/// Riemann identity checked by [`validate_riemann_symmetries`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RiemannIdentity {
    /// `R_IJKL = -R_JIKL`.
    FirstPairAntisymmetry,
    /// `R_IJKL = -R_IJLK`.
    SecondPairAntisymmetry,
    /// `R_IJKL = R_KLIJ`.
    PairExchange,
    /// `R_IJKL + R_IKLJ + R_ILJK = 0`.
    FirstBianchi,
    /// Validation requires a finite, non-negative absolute tolerance.
    InvalidTolerance,
}

/// A component-level failure of a field-space curvature identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeometryValidationError {
    /// The identity that failed.
    pub identity: RiemannIdentity,
    /// Component indices `(I, J, K, L)` at the first failure.
    pub indices: [usize; 4],
    /// Absolute residual of the identity.
    pub residual: f64,
    /// Requested absolute tolerance.
    pub tolerance: f64,
}

/// Validates the Levi-Civita symmetries of an all-lower Riemann tensor.
pub fn validate_riemann_symmetries<const N: usize>(
    riemann: &LoweredRiemannTensor<N>,
    tolerance: f64,
) -> Result<(), GeometryValidationError> {
    if tolerance < 0.0 || !tolerance.is_finite() {
        return Err(GeometryValidationError {
            identity: RiemannIdentity::InvalidTolerance,
            indices: [0; 4],
            residual: tolerance,
            tolerance,
        });
    }

    let component = |first: usize, second: usize, third: usize, fourth: usize| {
        *riemann.components[first][second][third][fourth].value()
    };

    for first in 0..N {
        for second in 0..N {
            for third in 0..N {
                for fourth in 0..N {
                    let indices = [first, second, third, fourth];
                    let checks = [
                        (
                            RiemannIdentity::FirstPairAntisymmetry,
                            component(first, second, third, fourth)
                                + component(second, first, third, fourth),
                        ),
                        (
                            RiemannIdentity::SecondPairAntisymmetry,
                            component(first, second, third, fourth)
                                + component(first, second, fourth, third),
                        ),
                        (
                            RiemannIdentity::PairExchange,
                            component(first, second, third, fourth)
                                - component(third, fourth, first, second),
                        ),
                        (
                            RiemannIdentity::FirstBianchi,
                            component(first, second, third, fourth)
                                + component(first, third, fourth, second)
                                + component(first, fourth, second, third),
                        ),
                    ];

                    for (identity, residual) in checks {
                        if residual.abs() > tolerance {
                            return Err(GeometryValidationError {
                                identity,
                                indices,
                                residual: residual.abs(),
                                tolerance,
                            });
                        }
                    }
                }
            }
        }
    }

    Ok(())
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
    use crate::{M1, M2, M3, MMinus2};

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

    #[test]
    fn rank_three_arithmetic_preserves_every_index() {
        let mut lhs_components = [[[Quantity::new(0.0); 2]; 2]; 2];
        let mut rhs_components = [[[Quantity::new(0.0); 2]; 2]; 2];
        lhs_components[1][0][1] = Quantity::new(2.0);
        rhs_components[1][0][1] = Quantity::new(3.0);
        let lhs = FieldTensor3::<Upper, Lower, Lower, M1, f64, 2>::new(lhs_components);
        let rhs = FieldTensor3::<Upper, Lower, Lower, M1, f64, 2>::new(rhs_components);

        let result = (lhs + rhs) * 2.0;

        assert_eq!(result.components()[1][0][1].into_value(), 10.0);
    }

    #[test]
    fn rank_four_first_index_lowering_and_raising_round_trip() {
        let metric = FieldMetric::<2>::new([
            [Quantity::new(2.0), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(3.0)],
        ]);
        let inverse = InverseFieldMetric::<2>::new([
            [Quantity::new(0.5), Quantity::new(0.0)],
            [Quantity::new(0.0), Quantity::new(1.0 / 3.0)],
        ]);
        let mut components = [[[[Quantity::new(0.0); 2]; 2]; 2]; 2];
        components[0][0][1][1] = Quantity::new(5.0);
        components[1][0][1][1] = Quantity::new(7.0);
        let tensor = FieldTensor4::<Upper, Lower, Lower, Lower, MMinus2, f64, 2>::new(components);

        let lowered = metric.lower_first(&tensor);
        let restored = inverse.raise_first(&lowered);

        assert_eq!(lowered.components()[0][0][1][1].into_value(), 10.0);
        assert_eq!(lowered.components()[1][0][1][1].into_value(), 21.0);
        assert_eq!(restored, tensor);
    }

    fn constant_curvature_two_field_space(curvature: f64) -> LoweredRiemannTensor<2> {
        LoweredRiemannTensor::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    core::array::from_fn(|fourth| {
                        let first_term = f64::from(first == third && second == fourth);
                        let second_term = f64::from(first == fourth && second == third);
                        Quantity::new(curvature * (first_term - second_term))
                    })
                })
            })
        }))
    }

    #[test]
    fn constant_curvature_tensor_satisfies_riemann_identities() {
        let lowered = constant_curvature_two_field_space(2.0);

        assert_eq!(validate_riemann_symmetries(&lowered, 1.0e-14), Ok(()));
    }

    #[test]
    fn riemann_velocity_contraction_has_expected_indices_and_dimension() {
        let lowered = constant_curvature_two_field_space(2.0);
        let riemann: RiemannTensor<2> = InverseFieldMetric::<2>::euclidean().raise_first(&lowered);
        let velocity =
            ContravariantField::<M2, f64, 2>::new([Quantity::new(3.0), Quantity::new(4.0)]);

        let contracted: FieldTensor2<Upper, Lower, M2, f64, 2> =
            riemann.contract_velocities(&velocity);
        let values = contracted
            .into_components()
            .map(|row| row.map(Quantity::into_value));

        assert_eq!(values, [[-32.0, 24.0], [24.0, -18.0]]);
    }

    #[test]
    fn validation_reports_the_failed_identity_and_component() {
        let mut components = constant_curvature_two_field_space(1.0).into_components();
        components[0][1][0][1] = Quantity::new(2.0);
        let invalid = LoweredRiemannTensor::new(components);

        let error = validate_riemann_symmetries(&invalid, 1.0e-14).unwrap_err();

        assert_eq!(error.identity, RiemannIdentity::FirstBianchi);
        assert_eq!(error.indices, [0, 0, 1, 1]);
        assert_eq!(error.residual, 1.0);
    }

    #[test]
    fn flat_geometry_supplies_exactly_zero_curvature() {
        let position = [Field::new(1.0), Field::new(2.0)];
        let riemann =
            <FlatFieldSpace as CurvedFieldSpaceGeometry<2>>::riemann(&FlatFieldSpace, &position);

        assert!(
            riemann
                .components()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .all(|component| component.into_value() == 0.0)
        );
    }

    #[test]
    fn flat_geometry_supplies_zero_third_order_geometry() {
        let position = [Field::new(1.0), Field::new(2.0)];
        let connection_derivative = FlatFieldSpace.christoffel_derivative(&position);
        let curvature_derivative = FlatFieldSpace.covariant_riemann_derivative(&position);

        assert!(
            connection_derivative
                .components()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .all(|component| component.into_value() == 0.0)
        );
        assert!(
            curvature_derivative
                .components()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .flatten()
                .all(|component| component.into_value() == 0.0)
        );
    }

    #[test]
    fn the_polar_metric_and_its_inverse_are_reciprocal() {
        let position = [Field::new(0.3), Field::new(2.0)];
        let metric = PolarFieldSpace.metric(&position);
        let inverse = PolarFieldSpace.inverse_metric(&position);

        for row in 0..2 {
            for column in 0..2 {
                let product: f64 = (0..2)
                    .map(|index| {
                        metric.components()[row][index].into_value()
                            * inverse.components()[index][column].into_value()
                    })
                    .sum();

                assert!((product - f64::from(row == column)).abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn the_polar_connection_is_the_levi_civita_connection_of_its_metric() {
        // Gamma^I_JK = G^IL (partial_J G_LK + partial_K G_LJ - partial_L G_JK) / 2
        let position = [Field::new(0.3), Field::new(2.0)];
        let step = 1.0e-6;
        let metric_derivative = |derivative: usize| {
            let mut ahead = position;
            let mut behind = position;
            ahead[derivative] = Field::new(position[derivative].into_value() + step);
            behind[derivative] = Field::new(position[derivative].into_value() - step);
            let forward = PolarFieldSpace.metric(&ahead);
            let backward = PolarFieldSpace.metric(&behind);
            core::array::from_fn::<_, 2, _>(|first| {
                core::array::from_fn::<_, 2, _>(|second| {
                    (forward.components()[first][second].into_value()
                        - backward.components()[first][second].into_value())
                        / (2.0 * step)
                })
            })
        };
        let derivatives = [metric_derivative(0), metric_derivative(1)];
        let partial =
            |derivative: usize, first: usize, second: usize| derivatives[derivative][first][second];
        let inverse = PolarFieldSpace.inverse_metric(&position);
        let connection = PolarFieldSpace.christoffel(&position);

        for upper in 0..2 {
            for first in 0..2 {
                for second in 0..2 {
                    let expected: f64 = (0..2)
                        .map(|index| {
                            0.5 * inverse.components()[upper][index].into_value()
                                * (partial(first, index, second) + partial(second, index, first)
                                    - partial(index, first, second))
                        })
                        .sum();
                    let analytic = connection.components()[upper][first][second].into_value();

                    assert!((analytic - expected).abs() < 1.0e-8);
                }
            }
        }
    }

    #[test]
    fn polar_connection_derivative_matches_finite_differences() {
        let position = [Field::new(0.3), Field::new(2.0)];
        let analytic = PolarFieldSpace.christoffel_derivative(&position);
        let step = 1.0e-6;

        for derivative in 0..2 {
            let mut forward = position;
            let mut backward = position;
            forward[derivative] = Field::new(position[derivative].value() + step);
            backward[derivative] = Field::new(position[derivative].value() - step);
            let forward_connection = PolarFieldSpace.christoffel(&forward);
            let backward_connection = PolarFieldSpace.christoffel(&backward);

            for upper in 0..2 {
                for first in 0..2 {
                    for second in 0..2 {
                        let numerical = (forward_connection.components()[upper][first][second]
                            .value()
                            - backward_connection.components()[upper][first][second].value())
                            / (2.0 * step);
                        let expected =
                            analytic.components()[upper][first][second][derivative].value();

                        assert!((numerical - expected).abs() < 1.0e-7);
                    }
                }
            }
        }
    }

    #[test]
    fn the_polar_connection_has_vanishing_curvature() {
        // R^I_JKL = partial_K Gamma^I_LJ - partial_L Gamma^I_KJ
        //         + Gamma^I_KM Gamma^M_LJ - Gamma^I_LM Gamma^M_KJ
        let position = [Field::new(0.3), Field::new(2.0)];
        let gamma = PolarFieldSpace.christoffel(&position);
        let partial = PolarFieldSpace.christoffel_derivative(&position);
        let g = |upper: usize, first: usize, second: usize| {
            gamma.components()[upper][first][second].into_value()
        };
        let dg = |upper: usize, first: usize, second: usize, derivative: usize| {
            partial.components()[upper][first][second][derivative].into_value()
        };

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        let quadratic: f64 = (0..2)
                            .map(|m| g(i, k, m) * g(m, l, j) - g(i, l, m) * g(m, k, j))
                            .sum();
                        let riemann = dg(i, l, j, k) - dg(i, k, j, l) + quadratic;

                        assert!(riemann.abs() < 1.0e-12, "R^{i}_{j}{k}{l} = {riemann}");
                    }
                }
            }
        }
    }

    #[test]
    fn the_hyperbolic_geometry_satisfies_the_riemann_identities() {
        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        let position = [Field::new(3.0), Field::new(1.0)];
        let lowered = geometry
            .metric(&position)
            .lower_first(&geometry.riemann(&position));

        validate_riemann_symmetries(&lowered, 1.0e-12).unwrap();
    }

    #[test]
    fn hyperbolic_connection_derivative_matches_finite_differences() {
        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        let position = [Field::new(3.0), Field::new(1.0)];
        let analytic = geometry.christoffel_derivative(&position);
        let step = 1.0e-5;

        for derivative in 0..2 {
            let mut forward = position;
            let mut backward = position;
            forward[derivative] = Field::new(position[derivative].value() + step);
            backward[derivative] = Field::new(position[derivative].value() - step);
            let forward_connection = geometry.christoffel(&forward);
            let backward_connection = geometry.christoffel(&backward);

            for upper in 0..2 {
                for first in 0..2 {
                    for second in 0..2 {
                        let numerical = (forward_connection.components()[upper][first][second]
                            .value()
                            - backward_connection.components()[upper][first][second].value())
                            / (2.0 * step);
                        let expected =
                            analytic.components()[upper][first][second][derivative].value();

                        assert!((numerical - expected).abs() < 1.0e-10);
                    }
                }
            }
        }
    }

    #[test]
    fn hyperbolic_covariant_curvature_derivative_is_zero() {
        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        let derivative = geometry.covariant_riemann_derivative(&[Field::new(3.0), Field::new(1.0)]);

        assert!(
            derivative
                .components()
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .flatten()
                .all(|component| component.into_value() == 0.0)
        );
    }

    #[test]
    fn hyperbolic_geometry_rejects_an_invalid_curvature_radius() {
        for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                HyperbolicFieldSpace::new(Field::new(radius)),
                Err(InvalidCurvatureRadius)
            );
        }

        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        assert_eq!(geometry.curvature_radius(), Field::new(7.0));
    }

    #[test]
    fn the_hyperbolic_metric_and_its_inverse_are_reciprocal() {
        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        let position = [Field::new(3.0), Field::new(1.0)];
        let metric = geometry.metric(&position);
        let inverse = geometry.inverse_metric(&position);

        for row in 0..2 {
            for column in 0..2 {
                let product: f64 = (0..2)
                    .map(|index| {
                        metric.components()[row][index].into_value()
                            * inverse.components()[index][column].into_value()
                    })
                    .sum();

                assert!((product - f64::from(row == column)).abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn the_hyperbolic_curvature_scalar_is_minus_two_over_the_radius_squared() {
        // In two dimensions R = 2K, and K = -1/L^2 for this geometry.
        let radius = 7.0;
        let geometry = HyperbolicFieldSpace::new(Field::new(radius)).unwrap();
        let position = [Field::new(3.0), Field::new(1.0)];
        let riemann = geometry.riemann(&position);
        let inverse = geometry.inverse_metric(&position);

        // R = G^JL R^I_JIL
        let scalar: f64 = (0..2)
            .flat_map(|second| (0..2).map(move |fourth| (second, fourth)))
            .map(|(second, fourth)| {
                let trace: f64 = (0..2)
                    .map(|index| riemann.components()[index][second][index][fourth].into_value())
                    .sum();

                inverse.components()[second][fourth].into_value() * trace
            })
            .sum();

        assert!(
            (scalar + 2.0 / (radius * radius)).abs() < 1.0e-12,
            "R = {scalar}"
        );
    }

    #[test]
    fn the_hyperbolic_connection_is_metric_compatible() {
        // nabla_K G_IJ = 0, checked against a finite-difference derivative of
        // the metric.
        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        let position = [Field::new(3.0), Field::new(1.0)];
        let step = 1.0e-6;
        let connection = geometry.christoffel(&position);

        for lower in 0..2 {
            let mut ahead = position;
            let mut behind = position;
            ahead[lower] = Field::new(position[lower].into_value() + step);
            behind[lower] = Field::new(position[lower].into_value() - step);
            let forward = geometry.metric(&ahead);
            let backward = geometry.metric(&behind);

            for first in 0..2 {
                for second in 0..2 {
                    let derivative = (forward.components()[first][second].into_value()
                        - backward.components()[first][second].into_value())
                        / (2.0 * step);
                    let corrections: f64 = (0..2)
                        .map(|index| {
                            connection.components()[index][lower][first].into_value()
                                * geometry.metric(&position).components()[index][second]
                                    .into_value()
                                + connection.components()[index][lower][second].into_value()
                                    * geometry.metric(&position).components()[first][index]
                                        .into_value()
                        })
                        .sum();

                    assert!(
                        (derivative - corrections).abs() < 1.0e-6,
                        "nabla_{lower} G_{first}{second} = {}",
                        derivative - corrections
                    );
                }
            }
        }
    }
}
