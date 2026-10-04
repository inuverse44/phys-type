//! Forward-mode automatic differentiation for scalar functions.

use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::Sqrt;

/// Elementary transcendental functions supported by potential value types.
///
/// A potential written generically over this trait can be evaluated with
/// ordinary `f64`, first-order [`Dual`], second-order [`Jet`], or
/// [`ThirdOrderJet`] values without maintaining separate analytic derivative
/// formulas.
pub trait Transcendental: Sized {
    /// Natural exponential.
    fn exp(self) -> Self;
    /// Natural logarithm.
    fn ln(self) -> Self;
    /// Cosine in radians.
    fn cos(self) -> Self;
    /// Hyperbolic tangent.
    fn tanh(self) -> Self;
}

impl Transcendental for f64 {
    fn exp(self) -> Self {
        self.exp()
    }

    fn ln(self) -> Self {
        self.ln()
    }

    fn cos(self) -> Self {
        self.cos()
    }

    fn tanh(self) -> Self {
        self.tanh()
    }
}

/// A value and its derivative with respect to one active variable.
///
/// Constants have derivative zero and an active variable has derivative one.
/// Arithmetic propagates both components using the usual differentiation
/// rules.
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct Dual {
    /// The ordinary numerical value.
    pub value: f64,
    /// The derivative with respect to the active variable.
    pub derivative: f64,
}

impl Dual {
    /// Creates a constant, whose derivative is zero.
    pub const fn constant(value: f64) -> Self {
        Self {
            value,
            derivative: 0.0,
        }
    }

    /// Creates the active differentiation variable, seeded with derivative one.
    pub const fn variable(value: f64) -> Self {
        Self {
            value,
            derivative: 1.0,
        }
    }

    fn map_unary(self, value: f64, first_derivative: f64) -> Self {
        Self {
            value,
            derivative: first_derivative * self.derivative,
        }
    }
}

impl Transcendental for Dual {
    fn exp(self) -> Self {
        let value = self.value.exp();
        self.map_unary(value, value)
    }

    fn ln(self) -> Self {
        self.map_unary(self.value.ln(), 1.0 / self.value)
    }

    fn cos(self) -> Self {
        self.map_unary(self.value.cos(), -self.value.sin())
    }

    fn tanh(self) -> Self {
        let value = self.value.tanh();
        self.map_unary(value, 1.0 - value * value)
    }
}

impl From<f64> for Dual {
    fn from(value: f64) -> Self {
        Self::constant(value)
    }
}

impl Add for Dual {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            derivative: self.derivative + rhs.derivative,
        }
    }
}

impl Sub for Dual {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value - rhs.value,
            derivative: self.derivative - rhs.derivative,
        }
    }
}

impl Neg for Dual {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            value: -self.value,
            derivative: -self.derivative,
        }
    }
}

impl Mul for Dual {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value * rhs.value,
            derivative: self.derivative * rhs.value + self.value * rhs.derivative,
        }
    }
}

impl Div for Dual {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value / rhs.value,
            derivative: (self.derivative * rhs.value - self.value * rhs.derivative)
                / (rhs.value * rhs.value),
        }
    }
}

impl Mul<f64> for Dual {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            value: self.value * rhs,
            derivative: self.derivative * rhs,
        }
    }
}

impl Mul<Dual> for f64 {
    type Output = Dual;

    fn mul(self, rhs: Dual) -> Self::Output {
        rhs * self
    }
}

impl Div<f64> for Dual {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            value: self.value / rhs,
            derivative: self.derivative / rhs,
        }
    }
}

impl Sqrt for Dual {
    fn sqrt_value(self) -> Self {
        let root = self.value.sqrt();
        Self {
            value: root,
            derivative: self.derivative / (2.0 * root),
        }
    }
}

/// A value, gradient, and Hessian with respect to `N` active variables.
///
/// `Jet` is a second-order forward-mode automatic-differentiation value. All
/// first and second partial derivatives propagate through ordinary arithmetic,
/// so a scalar formula needs to be evaluated only once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Jet<const N: usize> {
    /// The ordinary numerical value.
    pub value: f64,
    /// First derivatives in active-variable order.
    pub gradient: [f64; N],
    /// Symmetric second derivatives in row-major order.
    pub hessian: [[f64; N]; N],
}

impl<const N: usize> Jet<N> {
    /// Creates a constant with vanishing derivatives.
    pub const fn constant(value: f64) -> Self {
        Self {
            value,
            gradient: [0.0; N],
            hessian: [[0.0; N]; N],
        }
    }

    /// Creates active variable `index` with a unit gradient seed.
    pub fn variable(value: f64, index: usize) -> Self {
        assert!(index < N, "active variable index is out of bounds");
        let mut result = Self::constant(value);
        result.gradient[index] = 1.0;
        result
    }

    fn map_unary(self, value: f64, first_derivative: f64, second_derivative: f64) -> Self {
        Self {
            value,
            gradient: self.gradient.map(|entry| first_derivative * entry),
            hessian: core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    second_derivative * self.gradient[row] * self.gradient[column]
                        + first_derivative * self.hessian[row][column]
                })
            }),
        }
    }
}

impl<const N: usize> Transcendental for Jet<N> {
    fn exp(self) -> Self {
        let value = self.value.exp();
        self.map_unary(value, value, value)
    }

    fn ln(self) -> Self {
        let first = 1.0 / self.value;
        self.map_unary(self.value.ln(), first, -first * first)
    }

    fn cos(self) -> Self {
        self.map_unary(self.value.cos(), -self.value.sin(), -self.value.cos())
    }

    fn tanh(self) -> Self {
        let value = self.value.tanh();
        let first = 1.0 - value * value;
        self.map_unary(value, first, -2.0 * value * first)
    }
}

impl<const N: usize> Default for Jet<N> {
    fn default() -> Self {
        Self::constant(0.0)
    }
}

impl<const N: usize> From<f64> for Jet<N> {
    fn from(value: f64) -> Self {
        Self::constant(value)
    }
}

impl<const N: usize> Add for Jet<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            gradient: core::array::from_fn(|index| self.gradient[index] + rhs.gradient[index]),
            hessian: core::array::from_fn(|row| {
                core::array::from_fn(|column| self.hessian[row][column] + rhs.hessian[row][column])
            }),
        }
    }
}

impl<const N: usize> Sub for Jet<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value - rhs.value,
            gradient: core::array::from_fn(|index| self.gradient[index] - rhs.gradient[index]),
            hessian: core::array::from_fn(|row| {
                core::array::from_fn(|column| self.hessian[row][column] - rhs.hessian[row][column])
            }),
        }
    }
}

impl<const N: usize> Neg for Jet<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            value: -self.value,
            gradient: self.gradient.map(|entry| -entry),
            hessian: self.hessian.map(|row| row.map(|entry| -entry)),
        }
    }
}

// The additions are the first- and second-order product rules, not an
// accidental replacement of multiplication.
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const N: usize> Mul for Jet<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value * rhs.value,
            gradient: core::array::from_fn(|index| {
                self.gradient[index] * rhs.value + self.value * rhs.gradient[index]
            }),
            hessian: core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    self.hessian[row][column] * rhs.value
                        + self.gradient[row] * rhs.gradient[column]
                        + self.gradient[column] * rhs.gradient[row]
                        + self.value * rhs.hessian[row][column]
                })
            }),
        }
    }
}

impl<const N: usize> Div for Jet<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let inverse_value = 1.0 / rhs.value;
        let inverse = Self {
            value: inverse_value,
            gradient: core::array::from_fn(|index| {
                -rhs.gradient[index] * inverse_value * inverse_value
            }),
            hessian: core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    2.0 * rhs.gradient[row] * rhs.gradient[column] * inverse_value.powi(3)
                        - rhs.hessian[row][column] * inverse_value * inverse_value
                })
            }),
        };

        self * inverse
    }
}

impl<const N: usize> Mul<f64> for Jet<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            value: self.value * rhs,
            gradient: self.gradient.map(|entry| entry * rhs),
            hessian: self.hessian.map(|row| row.map(|entry| entry * rhs)),
        }
    }
}

impl<const N: usize> Mul<Jet<N>> for f64 {
    type Output = Jet<N>;

    fn mul(self, rhs: Jet<N>) -> Self::Output {
        rhs * self
    }
}

impl<const N: usize> Div<f64> for Jet<N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self * (1.0 / rhs)
    }
}

impl<const N: usize> Sqrt for Jet<N> {
    fn sqrt_value(self) -> Self {
        let root = self.value.sqrt();
        let first = 1.0 / (2.0 * root);
        let second = -1.0 / (4.0 * self.value * root);

        Self {
            value: root,
            gradient: self.gradient.map(|entry| first * entry),
            hessian: core::array::from_fn(|row| {
                core::array::from_fn(|column| {
                    second * self.gradient[row] * self.gradient[column]
                        + first * self.hessian[row][column]
                })
            }),
        }
    }
}

/// A value and every derivative through third order for `N` active variables.
///
/// This is separate from [`Jet`] so adding cubic transport support does not
/// enlarge or change the established value/gradient/Hessian representation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThirdOrderJet<const N: usize> {
    /// Ordinary numerical value.
    pub value: f64,
    /// First partial derivatives.
    pub gradient: [f64; N],
    /// Symmetric second partial derivatives.
    pub hessian: [[f64; N]; N],
    /// Fully symmetric third partial derivatives.
    pub third: [[[f64; N]; N]; N],
}

impl<const N: usize> ThirdOrderJet<N> {
    /// Creates a constant with every derivative set to zero.
    pub const fn constant(value: f64) -> Self {
        Self {
            value,
            gradient: [0.0; N],
            hessian: [[0.0; N]; N],
            third: [[[0.0; N]; N]; N],
        }
    }

    /// Creates active variable `index` with a unit gradient seed.
    pub fn variable(value: f64, index: usize) -> Self {
        assert!(index < N, "active variable index is out of bounds");
        let mut result = Self::constant(value);
        result.gradient[index] = 1.0;
        result
    }

    fn map_unary(
        self,
        value: f64,
        first_derivative: f64,
        second_derivative: f64,
        third_derivative: f64,
    ) -> Self {
        Self {
            value,
            gradient: self.gradient.map(|entry| first_derivative * entry),
            hessian: core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    second_derivative * self.gradient[first] * self.gradient[second]
                        + first_derivative * self.hessian[first][second]
                })
            }),
            third: core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    core::array::from_fn(|third| {
                        third_derivative
                            * self.gradient[first]
                            * self.gradient[second]
                            * self.gradient[third]
                            + second_derivative
                                * (self.hessian[first][second] * self.gradient[third]
                                    + self.hessian[first][third] * self.gradient[second]
                                    + self.hessian[second][third] * self.gradient[first])
                            + first_derivative * self.third[first][second][third]
                    })
                })
            }),
        }
    }

    fn reciprocal(self) -> Self {
        let inverse = 1.0 / self.value;
        self.map_unary(
            inverse,
            -inverse.powi(2),
            2.0 * inverse.powi(3),
            -6.0 * inverse.powi(4),
        )
    }
}

impl<const N: usize> Default for ThirdOrderJet<N> {
    fn default() -> Self {
        Self::constant(0.0)
    }
}

impl<const N: usize> From<f64> for ThirdOrderJet<N> {
    fn from(value: f64) -> Self {
        Self::constant(value)
    }
}

impl<const N: usize> Transcendental for ThirdOrderJet<N> {
    fn exp(self) -> Self {
        let value = self.value.exp();
        self.map_unary(value, value, value, value)
    }

    fn ln(self) -> Self {
        let inverse = 1.0 / self.value;
        self.map_unary(
            self.value.ln(),
            inverse,
            -inverse.powi(2),
            2.0 * inverse.powi(3),
        )
    }

    fn cos(self) -> Self {
        self.map_unary(
            self.value.cos(),
            -self.value.sin(),
            -self.value.cos(),
            self.value.sin(),
        )
    }

    fn tanh(self) -> Self {
        let value = self.value.tanh();
        let first = 1.0 - value * value;
        let second = -2.0 * value * first;
        let third = -2.0 * first * first - 2.0 * value * second;
        self.map_unary(value, first, second, third)
    }
}

impl<const N: usize> Add for ThirdOrderJet<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            gradient: core::array::from_fn(|index| self.gradient[index] + rhs.gradient[index]),
            hessian: core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    self.hessian[first][second] + rhs.hessian[first][second]
                })
            }),
            third: core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    core::array::from_fn(|third| {
                        self.third[first][second][third] + rhs.third[first][second][third]
                    })
                })
            }),
        }
    }
}

impl<const N: usize> Sub for ThirdOrderJet<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<const N: usize> Neg for ThirdOrderJet<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            value: -self.value,
            gradient: self.gradient.map(|entry| -entry),
            hessian: self.hessian.map(|row| row.map(|entry| -entry)),
            third: self
                .third
                .map(|matrix| matrix.map(|row| row.map(|entry| -entry))),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const N: usize> Mul for ThirdOrderJet<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value * rhs.value,
            gradient: core::array::from_fn(|index| {
                self.gradient[index] * rhs.value + self.value * rhs.gradient[index]
            }),
            hessian: core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    self.hessian[first][second] * rhs.value
                        + self.gradient[first] * rhs.gradient[second]
                        + self.gradient[second] * rhs.gradient[first]
                        + self.value * rhs.hessian[first][second]
                })
            }),
            third: core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    core::array::from_fn(|third| {
                        self.third[first][second][third] * rhs.value
                            + self.hessian[first][second] * rhs.gradient[third]
                            + self.hessian[first][third] * rhs.gradient[second]
                            + self.gradient[first] * rhs.hessian[second][third]
                            + self.hessian[second][third] * rhs.gradient[first]
                            + self.gradient[second] * rhs.hessian[first][third]
                            + self.gradient[third] * rhs.hessian[first][second]
                            + self.value * rhs.third[first][second][third]
                    })
                })
            }),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const N: usize> Div for ThirdOrderJet<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.reciprocal()
    }
}

impl<const N: usize> Mul<f64> for ThirdOrderJet<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            value: self.value * rhs,
            gradient: self.gradient.map(|entry| entry * rhs),
            hessian: self.hessian.map(|row| row.map(|entry| entry * rhs)),
            third: self
                .third
                .map(|matrix| matrix.map(|row| row.map(|entry| entry * rhs))),
        }
    }
}

impl<const N: usize> Mul<ThirdOrderJet<N>> for f64 {
    type Output = ThirdOrderJet<N>;

    fn mul(self, rhs: ThirdOrderJet<N>) -> Self::Output {
        rhs * self
    }
}

impl<const N: usize> Div<f64> for ThirdOrderJet<N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self * (1.0 / rhs)
    }
}

impl<const N: usize> Sqrt for ThirdOrderJet<N> {
    fn sqrt_value(self) -> Self {
        let root = self.value.sqrt();
        self.map_unary(
            root,
            1.0 / (2.0 * root),
            -1.0 / (4.0 * root.powi(3)),
            3.0 / (8.0 * root.powi(5)),
        )
    }
}

/// A value and its first four partial derivatives with respect to `N`
/// variables.
///
/// The Einstein-frame metric of a nonminimal theory contains first
/// derivatives of the conformal factor, and the cubic transport needs the
/// covariant derivative of the Riemann tensor, a third derivative of the
/// metric. The conformal factor is therefore differentiated four times.
///
/// Products follow the Leibniz rule as a sum over subsets of the derivative
/// indices, and composition with a univariate function follows Faà di Bruno's
/// formula as a sum over set partitions of the indices, so every component
/// is computed by the same rule rather than by order-specific formulas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FourthOrderJet<const N: usize> {
    /// Ordinary numerical value.
    pub value: f64,
    /// First partial derivatives.
    pub gradient: [f64; N],
    /// Symmetric second partial derivatives.
    pub hessian: [[f64; N]; N],
    /// Fully symmetric third partial derivatives.
    pub third: [[[f64; N]; N]; N],
    /// Fully symmetric fourth partial derivatives.
    pub fourth: [[[[f64; N]; N]; N]; N],
}

/// Set partitions of the positions `0..k` for `k <= 4`, each block a bit
/// mask over the positions.
fn set_partitions(k: usize) -> Vec<Vec<u8>> {
    fn extend(position: usize, k: usize, blocks: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        if position == k {
            out.push(blocks.clone());
            return;
        }
        for index in 0..blocks.len() {
            blocks[index] |= 1 << position;
            extend(position + 1, k, blocks, out);
            blocks[index] &= !(1 << position);
        }
        blocks.push(1 << position);
        extend(position + 1, k, blocks, out);
        blocks.pop();
    }
    let mut out = Vec::new();
    extend(0, k, &mut Vec::new(), &mut out);
    out
}

impl<const N: usize> FourthOrderJet<N> {
    /// Creates a constant with every derivative set to zero.
    pub const fn constant(value: f64) -> Self {
        Self {
            value,
            gradient: [0.0; N],
            hessian: [[0.0; N]; N],
            third: [[[0.0; N]; N]; N],
            fourth: [[[[0.0; N]; N]; N]; N],
        }
    }

    /// Creates active variable `index` with a unit gradient seed.
    pub fn variable(value: f64, index: usize) -> Self {
        assert!(index < N, "active variable index is out of bounds");
        let mut result = Self::constant(value);
        result.gradient[index] = 1.0;
        result
    }

    /// The derivative with respect to the variables in `indices`, of order
    /// `indices.len() <= 4`; the empty slice gives the value.
    pub fn derivative(&self, indices: &[usize]) -> f64 {
        match *indices {
            [] => self.value,
            [i] => self.gradient[i],
            [i, j] => self.hessian[i][j],
            [i, j, k] => self.third[i][j][k],
            [i, j, k, l] => self.fourth[i][j][k][l],
            _ => panic!("FourthOrderJet carries derivatives up to fourth order"),
        }
    }

    /// The jet whose every component is `component(indices)`.
    fn from_components(component: impl Fn(&[usize]) -> f64) -> Self {
        Self {
            value: component(&[]),
            gradient: core::array::from_fn(|i| component(&[i])),
            hessian: core::array::from_fn(|i| core::array::from_fn(|j| component(&[i, j]))),
            third: core::array::from_fn(|i| {
                core::array::from_fn(|j| core::array::from_fn(|k| component(&[i, j, k])))
            }),
            fourth: core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| core::array::from_fn(|l| component(&[i, j, k, l])))
                })
            }),
        }
    }

    /// `indices` restricted to the positions in `mask`.
    fn select(indices: &[usize], mask: u8) -> ([usize; 4], usize) {
        let mut selected = [0; 4];
        let mut count = 0;
        for (position, &index) in indices.iter().enumerate() {
            if mask & (1 << position) != 0 {
                selected[count] = index;
                count += 1;
            }
        }
        (selected, count)
    }

    /// `h(self)` given `h` and its first four derivatives at `self.value`.
    fn map_unary(self, derivatives: [f64; 5]) -> Self {
        let partitions: [Vec<Vec<u8>>; 5] = core::array::from_fn(set_partitions);
        Self::from_components(|indices| {
            if indices.is_empty() {
                return derivatives[0];
            }
            partitions[indices.len()]
                .iter()
                .map(|blocks| {
                    derivatives[blocks.len()]
                        * blocks
                            .iter()
                            .map(|&mask| {
                                let (selected, count) = Self::select(indices, mask);
                                self.derivative(&selected[..count])
                            })
                            .product::<f64>()
                })
                .sum()
        })
    }

    fn reciprocal(self) -> Self {
        let inverse = 1.0 / self.value;
        self.map_unary([
            inverse,
            -inverse.powi(2),
            2.0 * inverse.powi(3),
            -6.0 * inverse.powi(4),
            24.0 * inverse.powi(5),
        ])
    }

    fn zip(self, rhs: Self, combine: impl Fn(f64, f64) -> f64) -> Self {
        Self::from_components(|indices| combine(self.derivative(indices), rhs.derivative(indices)))
    }
}

impl<const N: usize> Default for FourthOrderJet<N> {
    fn default() -> Self {
        Self::constant(0.0)
    }
}

impl<const N: usize> From<f64> for FourthOrderJet<N> {
    fn from(value: f64) -> Self {
        Self::constant(value)
    }
}

impl<const N: usize> Transcendental for FourthOrderJet<N> {
    fn exp(self) -> Self {
        let value = self.value.exp();
        self.map_unary([value; 5])
    }

    fn ln(self) -> Self {
        let inverse = 1.0 / self.value;
        self.map_unary([
            self.value.ln(),
            inverse,
            -inverse.powi(2),
            2.0 * inverse.powi(3),
            -6.0 * inverse.powi(4),
        ])
    }

    fn cos(self) -> Self {
        let (sin, cos) = self.value.sin_cos();
        self.map_unary([cos, -sin, -cos, sin, cos])
    }

    fn tanh(self) -> Self {
        let t = self.value.tanh();
        let first = 1.0 - t * t;
        // d/dx of a polynomial p(t) is p'(t) (1 - t^2).
        let second = -2.0 * t * first;
        let third = (-2.0 + 6.0 * t * t) * first;
        let fourth = (16.0 * t - 24.0 * t * t * t) * first;
        self.map_unary([t, first, second, third, fourth])
    }
}

impl<const N: usize> Add for FourthOrderJet<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.zip(rhs, |a, b| a + b)
    }
}

impl<const N: usize> Sub for FourthOrderJet<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.zip(rhs, |a, b| a - b)
    }
}

impl<const N: usize> Neg for FourthOrderJet<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self * -1.0
    }
}

impl<const N: usize> Mul for FourthOrderJet<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::from_components(|indices| {
            let full = (1u8 << indices.len()) - 1;
            (0..=full)
                .map(|mask| {
                    let (left, left_count) = Self::select(indices, mask);
                    let (right, right_count) = Self::select(indices, full & !mask);
                    self.derivative(&left[..left_count]) * rhs.derivative(&right[..right_count])
                })
                .sum()
        })
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const N: usize> Div for FourthOrderJet<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.reciprocal()
    }
}

impl<const N: usize> Mul<f64> for FourthOrderJet<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::from_components(|indices| self.derivative(indices) * rhs)
    }
}

impl<const N: usize> Mul<FourthOrderJet<N>> for f64 {
    type Output = FourthOrderJet<N>;

    fn mul(self, rhs: FourthOrderJet<N>) -> Self::Output {
        rhs * self
    }
}

impl<const N: usize> Div<f64> for FourthOrderJet<N> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self * (1.0 / rhs)
    }
}

impl<const N: usize> Sqrt for FourthOrderJet<N> {
    fn sqrt_value(self) -> Self {
        let root = self.value.sqrt();
        self.map_unary([
            root,
            0.5 / root,
            -0.25 / root.powi(3),
            0.375 / root.powi(5),
            -0.9375 / root.powi(7),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn transcendental_formula<V>(x: V, y: V) -> V
    where
        V: Clone + Add<Output = V> + Sub<Output = V> + Mul<Output = V> + Transcendental,
    {
        x.clone().exp() + y.clone().ln() + (x.clone() * y.clone()).cos() + (x - y).tanh()
    }

    #[test]
    fn differentiates_a_polynomial() {
        let x = Dual::variable(3.0);
        let result = x * x * x + Dual::constant(2.0) * x;

        assert_eq!(result.value, 33.0);
        assert_eq!(result.derivative, 29.0);
    }

    #[test]
    fn differentiates_a_quotient_and_square_root() {
        let x = Dual::variable(4.0);
        let result = (x / Dual::constant(2.0)).sqrt_value();

        assert!((result.value - 2.0_f64.sqrt()).abs() < 1.0e-15);
        assert!((result.derivative - 1.0 / (4.0 * 2.0_f64.sqrt())).abs() < 1.0e-15);
    }

    #[test]
    fn jet_computes_a_full_gradient_and_hessian() {
        let x = Jet::<2>::variable(2.0, 0);
        let y = Jet::<2>::variable(3.0, 1);
        let result = x * x * y;

        assert_eq!(result.value, 12.0);
        assert_eq!(result.gradient, [12.0, 4.0]);
        assert_eq!(result.hessian, [[6.0, 4.0], [4.0, 0.0]]);
    }

    #[test]
    fn jet_propagates_quotients_and_square_roots() {
        let x = Jet::<2>::variable(4.0, 0);
        let y = Jet::<2>::variable(2.0, 1);
        let quotient = x / y;
        let root = x.sqrt_value();

        assert_eq!(quotient.gradient, [0.5, -1.0]);
        assert_eq!(quotient.hessian, [[0.0, -0.25], [-0.25, 1.0]]);
        assert_eq!(root.gradient, [0.25, 0.0]);
        assert_eq!(root.hessian, [[-0.03125, 0.0], [0.0, 0.0]]);
    }

    #[test]
    fn dual_differentiates_transcendental_functions() {
        let x = Dual::variable(0.4);

        let exponential = x.exp();
        let logarithm = x.ln();
        let cosine = x.cos();
        let hyperbolic_tangent = x.tanh();

        assert!((exponential.derivative - 0.4_f64.exp()).abs() < 1.0e-15);
        assert!((logarithm.derivative - 2.5).abs() < 1.0e-15);
        assert!((cosine.derivative + 0.4_f64.sin()).abs() < 1.0e-15);
        assert!((hyperbolic_tangent.derivative - (1.0 - 0.4_f64.tanh().powi(2))).abs() < 1.0e-15);
    }

    #[test]
    fn jet_transcendental_second_derivatives_match_closed_forms() {
        let x = Jet::<1>::variable(0.4, 0);
        let exponential = x.exp();
        let logarithm = x.ln();
        let cosine = x.cos();
        let hyperbolic_tangent = x.tanh();
        let tangent = 0.4_f64.tanh();

        assert!((exponential.hessian[0][0] - 0.4_f64.exp()).abs() < 1.0e-15);
        assert!((logarithm.hessian[0][0] + 1.0 / 0.4_f64.powi(2)).abs() < 1.0e-15);
        assert!((cosine.hessian[0][0] + 0.4_f64.cos()).abs() < 1.0e-15);
        assert!(
            (hyperbolic_tangent.hessian[0][0] + 2.0 * tangent * (1.0 - tangent * tangent)).abs()
                < 1.0e-15
        );
    }

    #[test]
    fn jet_transcendental_gradient_and_hessian_match_finite_differences() {
        let point = [0.4, 1.3];
        let x = Jet::<2>::variable(point[0], 0);
        let y = Jet::<2>::variable(point[1], 1);
        let result = transcendental_formula(x, y);
        let value = |x: f64, y: f64| transcendental_formula(x, y);
        let step = 1.0e-4;

        let gradient = [
            (value(point[0] + step, point[1]) - value(point[0] - step, point[1])) / (2.0 * step),
            (value(point[0], point[1] + step) - value(point[0], point[1] - step)) / (2.0 * step),
        ];
        let diagonal = [
            (value(point[0] + step, point[1]) - 2.0 * value(point[0], point[1])
                + value(point[0] - step, point[1]))
                / (step * step),
            (value(point[0], point[1] + step) - 2.0 * value(point[0], point[1])
                + value(point[0], point[1] - step))
                / (step * step),
        ];
        let mixed = (value(point[0] + step, point[1] + step)
            - value(point[0] + step, point[1] - step)
            - value(point[0] - step, point[1] + step)
            + value(point[0] - step, point[1] - step))
            / (4.0 * step * step);

        for index in 0..2 {
            assert!((result.gradient[index] - gradient[index]).abs() < 1.0e-7);
            assert!((result.hessian[index][index] - diagonal[index]).abs() < 1.0e-6);
        }
        assert!((result.hessian[0][1] - mixed).abs() < 1.0e-7);
        assert_eq!(result.hessian[0][1], result.hessian[1][0]);
    }

    #[test]
    fn third_order_jet_computes_all_mixed_permutations() {
        let x = ThirdOrderJet::<3>::variable(2.0, 0);
        let y = ThirdOrderJet::<3>::variable(3.0, 1);
        let z = ThirdOrderJet::<3>::variable(5.0, 2);
        let result = x * x * y + x * y * z;

        assert_eq!(result.value, 42.0);
        assert_eq!(result.gradient, [27.0, 14.0, 6.0]);
        assert_eq!(
            result.hessian,
            [[6.0, 9.0, 3.0], [9.0, 0.0, 2.0], [3.0, 2.0, 0.0]]
        );

        for permutation in [[0, 0, 1], [0, 1, 0], [1, 0, 0]] {
            assert_eq!(
                result.third[permutation[0]][permutation[1]][permutation[2]],
                2.0
            );
        }
        for permutation in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            assert_eq!(
                result.third[permutation[0]][permutation[1]][permutation[2]],
                1.0
            );
        }
    }

    #[test]
    fn third_order_jet_propagates_quotients_and_square_roots() {
        let x = ThirdOrderJet::<1>::variable(4.0, 0);
        let reciprocal = ThirdOrderJet::constant(1.0) / x;
        let root = x.sqrt_value();

        assert_eq!(reciprocal.third[0][0][0], -6.0 / 4.0_f64.powi(4));
        assert_eq!(root.third[0][0][0], 3.0 / (8.0 * 2.0_f64.powi(5)));
    }

    #[test]
    fn third_order_jet_transcendental_derivatives_match_closed_forms() {
        let point = 0.4_f64;
        let x = ThirdOrderJet::<1>::variable(point, 0);
        let exponential = x.exp();
        let logarithm = x.ln();
        let cosine = x.cos();
        let hyperbolic_tangent = x.tanh();
        let tangent = point.tanh();
        let first = 1.0 - tangent * tangent;
        let second = -2.0 * tangent * first;
        let expected_tanh_third = -2.0 * first * first - 2.0 * tangent * second;

        assert!((exponential.third[0][0][0] - point.exp()).abs() < 1.0e-15);
        assert!((logarithm.third[0][0][0] - 2.0 / point.powi(3)).abs() < 1.0e-14);
        assert!((cosine.third[0][0][0] - point.sin()).abs() < 1.0e-15);
        assert!((hyperbolic_tangent.third[0][0][0] - expected_tanh_third).abs() < 1.0e-15);
    }

    fn fourth_order_sample<V>(x: V, y: V) -> V
    where
        V: Clone
            + Add<Output = V>
            + Mul<Output = V>
            + Div<Output = V>
            + Mul<f64, Output = V>
            + From<f64>
            + Transcendental
            + Sqrt,
    {
        let product = x.clone() * y.clone();
        let denominator = V::from(1.0) + x.clone() * x.clone() * 0.3;
        (product.clone() * 0.7).exp() / denominator
            + (V::from(2.0) + y.clone()).ln() * x.clone().cos()
            + (y * 0.5).tanh() * (x + V::from(3.0)).sqrt_value()
    }

    #[test]
    fn fourth_order_jet_agrees_with_the_third_order_jet() {
        let (x0, y0) = (0.4, -0.3);
        let fourth = fourth_order_sample(
            FourthOrderJet::<2>::variable(x0, 0),
            FourthOrderJet::<2>::variable(y0, 1),
        );
        let third = fourth_order_sample(
            ThirdOrderJet::<2>::variable(x0, 0),
            ThirdOrderJet::<2>::variable(y0, 1),
        );
        assert!((fourth.value - third.value).abs() < 1.0e-14);
        for i in 0..2 {
            assert!((fourth.gradient[i] - third.gradient[i]).abs() < 1.0e-13);
            for j in 0..2 {
                assert!((fourth.hessian[i][j] - third.hessian[i][j]).abs() < 1.0e-12);
                for k in 0..2 {
                    assert!((fourth.third[i][j][k] - third.third[i][j][k]).abs() < 1.0e-11);
                }
            }
        }
    }

    #[test]
    fn fourth_order_jet_matches_finite_differences_of_third_derivatives() {
        let (x0, y0) = (0.4, -0.3);
        let fourth = fourth_order_sample(
            FourthOrderJet::<2>::variable(x0, 0),
            FourthOrderJet::<2>::variable(y0, 1),
        );
        let step = 1.0e-4;
        for l in 0..2 {
            let shifted = |sign: f64| {
                let (dx, dy) = if l == 0 {
                    (sign * step, 0.0)
                } else {
                    (0.0, sign * step)
                };
                fourth_order_sample(
                    ThirdOrderJet::<2>::variable(x0 + dx, 0),
                    ThirdOrderJet::<2>::variable(y0 + dy, 1),
                )
            };
            let (plus, minus) = (shifted(1.0), shifted(-1.0));
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let numerical = (plus.third[i][j][k] - minus.third[i][j][k]) / (2.0 * step);
                        let analytic = fourth.fourth[i][j][k][l];
                        assert!(
                            (numerical - analytic).abs() < 1.0e-6 * (1.0 + analytic.abs()),
                            "{i}{j}{k}{l}: {numerical} vs {analytic}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn fourth_order_jet_reproduces_closed_form_fourth_derivatives() {
        // d^4/dx^4 of exp(2x), ln x, cos x, tanh x, and sqrt x.
        let x0 = 0.7_f64;
        let x = FourthOrderJet::<1>::variable(x0, 0);
        let fourth = |jet: FourthOrderJet<1>| jet.fourth[0][0][0][0];
        let t = x0.tanh();
        let sech2 = 1.0 - t * t;
        let checks = [
            (fourth((x * 2.0).exp()), 16.0 * (2.0 * x0).exp()),
            (fourth(x.ln()), -6.0 / x0.powi(4)),
            (fourth(x.cos()), x0.cos()),
            (fourth(x.tanh()), sech2 * (16.0 * t - 24.0 * t.powi(3))),
            (fourth(x.sqrt_value()), -15.0 / 16.0 * x0.powf(-3.5)),
            (fourth(FourthOrderJet::from(1.0) / x), 24.0 / x0.powi(5)),
        ];
        for (index, (analytic, expected)) in checks.into_iter().enumerate() {
            assert!(
                (analytic - expected).abs() < 1.0e-12 * (1.0 + expected.abs()),
                "check {index}: {analytic} vs {expected}"
            );
        }
    }
}
