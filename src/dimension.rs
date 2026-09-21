//! Compile-time physical dimensions expressed as powers of mass.
//!
//! Natural units with `c = hbar = 1` leave a single independent dimension, so
//! every physical dimension is a power of mass. The exponent is stored in
//! half-power steps so that `sqrt(M^1) = M^(1/2)` stays representable:
//!
//! ```text
//! internal 1 -> M^(1/2)
//! internal 2 -> M^1
//! internal 4 -> M^2
//! ```

use core::marker::PhantomData;
use core::ops::{Add, Sub};

use typenum::{Diff, N1, N2, N4, P1, P2, P3, P4, P5, P6, P8, PartialDiv, PartialQuot, Sum, Z0};

/// A mass dimension whose exponent is encoded in half-power steps.
///
/// For example, `MassPower<P1>` is `M^(1/2)` and `MassPower<P2>` is `M^1`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MassPower<E>(PhantomData<E>);

/// Type-level multiplication of dimensions, which adds mass exponents.
pub trait DimMul<Rhs> {
    /// The dimension of the product.
    type Output;
}

impl<A, B> DimMul<MassPower<B>> for MassPower<A>
where
    A: Add<B>,
{
    type Output = MassPower<Sum<A, B>>;
}

/// Type-level division of dimensions, which subtracts mass exponents.
pub trait DimDiv<Rhs> {
    /// The dimension of the quotient.
    type Output;
}

impl<A, B> DimDiv<MassPower<B>> for MassPower<A>
where
    A: Sub<B>,
{
    type Output = MassPower<Diff<A, B>>;
}

/// Type-level square root, which halves the mass exponent.
///
/// This is implemented only when the internal half-power exponent is even, so
/// that the result stays on the supported half-integer grid.
pub trait DimSqrt {
    /// The dimension of the square root.
    type Output;
}

impl<A> DimSqrt for MassPower<A>
where
    A: PartialDiv<P2>,
{
    type Output = MassPower<PartialQuot<A, P2>>;
}

/// `M^0`, a dimensionless quantity.
pub type M0 = MassPower<Z0>;
/// `M^(1/2)`.
pub type MHalf = MassPower<P1>;
/// `M^1`.
pub type M1 = MassPower<P2>;
/// `M^(3/2)`.
pub type M3Over2 = MassPower<P3>;
/// `M^2`.
pub type M2 = MassPower<P4>;
/// `M^(5/2)`.
pub type M5Over2 = MassPower<P5>;
/// `M^3`.
pub type M3 = MassPower<P6>;
/// `M^4`.
pub type M4 = MassPower<P8>;
/// `M^-1`.
pub type MMinus1 = MassPower<N2>;
/// `M^(-1/2)`.
pub type MMinusHalf = MassPower<N1>;
/// `M^-2`.
pub type MMinus2 = MassPower<N4>;
