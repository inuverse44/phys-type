//! Homogeneous backgrounds containing a fixed number of scalar fields.

use core::ops::{Add, Mul};

use crate::integrate::System;
use crate::potential::MultifieldPotential;
use crate::{
    EnergyDensity, Field, FieldVelocity, Hubble, M1, M2, M3, Mass, PotentialGradient, Quantity,
    ScaleFactor, Time,
};

/// The values, velocities, and common scale factor of `N` scalar fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultifieldState<const N: usize> {
    /// Scalar field values `phi^I`.
    pub phi: [Field; N],
    /// Field velocities `dphi^I/dt`.
    pub pi: [FieldVelocity; N],
    /// The shared FLRW scale factor.
    pub a: ScaleFactor,
}

impl<const N: usize> MultifieldState<N> {
    /// Builds a multifield state.
    pub const fn new(phi: [Field; N], pi: [FieldVelocity; N], a: ScaleFactor) -> Self {
        Self { phi, pi, a }
    }

    /// The e-folding number `N = ln(a / a_initial)`.
    pub fn e_folds(&self, initial: &Self) -> f64 {
        (self.a / initial.a).into_value().ln()
    }
}

/// The time derivative of a [`MultifieldState`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultifieldStateDerivative<const N: usize> {
    /// `dphi^I/dt`.
    pub dphi_dt: [Quantity<M2, f64>; N],
    /// `d pi^I/dt`.
    pub dpi_dt: [Quantity<M3, f64>; N],
    /// `da/dt`.
    pub da_dt: Quantity<M1, f64>,
}

/// A finite change of a [`MultifieldState`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultifieldStateDelta<const N: usize> {
    /// Changes in the field values.
    pub phi: [Field; N],
    /// Changes in the field velocities.
    pub pi: [FieldVelocity; N],
    /// Change in the scale factor.
    pub a: ScaleFactor,
}

impl<const N: usize> Add<MultifieldStateDelta<N>> for MultifieldState<N> {
    type Output = Self;

    fn add(self, rhs: MultifieldStateDelta<N>) -> Self::Output {
        Self {
            phi: core::array::from_fn(|index| self.phi[index] + rhs.phi[index]),
            pi: core::array::from_fn(|index| self.pi[index] + rhs.pi[index]),
            a: self.a + rhs.a,
        }
    }
}

impl<const N: usize> Add for MultifieldStateDerivative<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            dphi_dt: core::array::from_fn(|index| self.dphi_dt[index] + rhs.dphi_dt[index]),
            dpi_dt: core::array::from_fn(|index| self.dpi_dt[index] + rhs.dpi_dt[index]),
            da_dt: self.da_dt + rhs.da_dt,
        }
    }
}

impl<const N: usize> Mul<f64> for MultifieldStateDerivative<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            dphi_dt: self.dphi_dt.map(|value| value * rhs),
            dpi_dt: self.dpi_dt.map(|value| value * rhs),
            da_dt: self.da_dt * rhs,
        }
    }
}

impl<const N: usize> Mul<Time> for MultifieldStateDerivative<N> {
    type Output = MultifieldStateDelta<N>;

    fn mul(self, rhs: Time) -> Self::Output {
        MultifieldStateDelta {
            phi: self.dphi_dt.map(|value| value * rhs),
            pi: self.dpi_dt.map(|value| value * rhs),
            a: self.da_dt * rhs,
        }
    }
}

/// `N` canonical scalar fields in a spatially flat FLRW background.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultifieldSystem<P, const N: usize> {
    /// The potential shared by the fields.
    pub potential: P,
    /// The reduced Planck mass.
    pub m_pl: Mass,
}

impl<P, const N: usize> MultifieldSystem<P, N>
where
    P: MultifieldPotential<N>,
{
    /// Builds a system with the given potential and reduced Planck mass.
    pub const fn new(potential: P, m_pl: Mass) -> Self {
        Self { potential, m_pl }
    }

    /// Builds a system in reduced Planck units.
    pub const fn in_planck_units(potential: P) -> Self {
        Self::new(potential, Mass::new(1.0))
    }

    /// Total canonical kinetic plus potential energy density.
    pub fn energy_density(&self, state: &MultifieldState<N>) -> EnergyDensity {
        let kinetic = state
            .pi
            .iter()
            .fold(EnergyDensity::new(0.0), |sum, velocity| {
                sum + *velocity * *velocity * 0.5
            });

        kinetic + self.potential.value(state.phi)
    }

    /// The expanding-branch Hubble rate.
    pub fn hubble(&self, state: &MultifieldState<N>) -> Hubble {
        (self.energy_density(state) / (self.m_pl * self.m_pl * 3.0)).sqrt()
    }

    /// `epsilon_H = sum_I (pi^I)^2 / (2 M_pl^2 H^2)`.
    pub fn epsilon_h(&self, state: &MultifieldState<N>) -> f64 {
        let kinetic_twice = state
            .pi
            .iter()
            .fold(EnergyDensity::new(0.0), |sum, velocity| {
                sum + *velocity * *velocity
            });
        let hubble = self.hubble(state);
        let denominator = self.m_pl * self.m_pl * hubble * hubble * 2.0;

        (kinetic_twice / denominator).into_value()
    }

    /// Slow-roll velocities `dphi^I/dt ~= -V_,I / (3H)`.
    pub fn slow_roll_velocity(&self, state: &MultifieldState<N>) -> [FieldVelocity; N] {
        let hubble = self.hubble(state);
        self.potential
            .gradient(state.phi)
            .map(|gradient: PotentialGradient| -(gradient / (hubble * 3.0)))
    }
}

impl<P, const N: usize> System for MultifieldSystem<P, N>
where
    P: MultifieldPotential<N>,
{
    type IndependentVariable = Time;
    type State = MultifieldState<N>;
    type Delta = MultifieldStateDelta<N>;
    type Derivative = MultifieldStateDerivative<N>;

    fn derivative(&self, _time: Time, state: &Self::State) -> Self::Derivative {
        let hubble = self.hubble(state);
        let gradient = self.potential.gradient(state.phi);

        MultifieldStateDerivative {
            dphi_dt: state.pi,
            dpi_dt: core::array::from_fn(|index| {
                state.pi[index] * hubble * (-3.0) - gradient[index]
            }),
            da_dt: state.a * hubble,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticMultifieldPotential, rk4_step, rk45_step};

    fn two_field_system() -> MultifieldSystem<QuadraticMultifieldPotential<2>, 2> {
        MultifieldSystem::in_planck_units(QuadraticMultifieldPotential::new([
            Mass::new(1.0e-6),
            Mass::new(2.0e-6),
        ]))
    }

    #[test]
    fn energy_density_sums_both_fields() {
        let system = two_field_system();
        let state = MultifieldState::new(
            [Field::new(3.0), Field::new(4.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );

        let expected = 0.5 * (1.0e-12 * 9.0 + 4.0e-12 * 16.0);
        assert!((system.energy_density(&state).into_value() - expected).abs() < 1.0e-25);
    }

    #[test]
    fn generic_integrators_advance_array_backed_states() {
        let system = two_field_system();
        let resting = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );
        let initial = MultifieldState {
            pi: system.slow_roll_velocity(&resting),
            ..resting
        };
        let time = Time::new(0.0);
        let dt = Time::new(100.0);

        let rk4 = rk4_step(&system, time, &initial, dt);
        let rk45 = rk45_step(&system, time, &initial, dt).state;

        assert!(rk4.a.into_value() > initial.a.into_value());
        assert!(rk4.phi[0].into_value() < initial.phi[0].into_value());
        assert!(rk4.phi[1].into_value() < initial.phi[1].into_value());
        assert!((rk4.phi[0].into_value() - rk45.phi[0].into_value()).abs() < 1.0e-12);
    }
}
