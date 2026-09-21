//! Forward-mode automatic differentiation for scalar functions.

use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::Sqrt;

/// Elementary transcendental functions supported by potential value types.
///
/// A potential written generically over this trait can be evaluated with
/// ordinary `f64`, first-order [`Dual`], or second-order [`Jet`] values without
/// maintaining separate analytic derivative formulas.
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
}
