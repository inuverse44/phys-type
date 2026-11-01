//! Semantically distinct time coordinates.

use core::ops::{Add, Mul};

use crate::Time;

/// Conformal time `d eta = dt / a`.
///
/// Cosmic and conformal time have the same physical dimension, but they are
/// different coordinates.  Keeping a separate wrapper prevents an integrator
/// step in one coordinate from being passed to a system written for the other.
#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct ConformalTime(Time);

impl ConformalTime {
    /// Attaches the conformal-time meaning to a value in Planck units.
    pub const fn new(value: f64) -> Self {
        Self(Time::new(value))
    }

    /// Returns the raw value in Planck units.
    pub fn into_value(self) -> f64 {
        self.0.into_value()
    }

    /// Returns the dimensioned interval for quantity arithmetic.
    pub(crate) const fn as_time(self) -> Time {
        self.0
    }
}

impl Add for ConformalTime {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Mul<f64> for ConformalTime {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conformal_time_has_time_arithmetic_without_becoming_cosmic_time() {
        let interval = ConformalTime::new(2.0) + ConformalTime::new(3.0) * 2.0;

        assert_eq!(interval.into_value(), 8.0);
    }
}
