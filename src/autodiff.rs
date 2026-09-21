//! Forward-mode automatic differentiation for scalar functions.

use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::Sqrt;

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

#[cfg(test)]
mod tests {
    use super::*;

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
}
