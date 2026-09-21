//! Homogeneous backgrounds containing a fixed number of scalar fields.

use core::ops::{Add, Mul};

use crate::integrate::System;
use crate::potential::MultifieldPotential;
use crate::{
    ContravariantField, CovariantPotentialHessian, CurvedFieldSpaceGeometry, EnergyDensity, Field,
    FieldSpaceGeometry, FieldTensor2, FieldVelocity, Hubble, KinematicBasis2, KinematicBasisError,
    Lower, M0, M1, M2, M3, Mass, PotentialDerivatives, PotentialGradient, PotentialHessian,
    Quantity, ScaleFactor, Time, TwiceDifferentiableMultifieldPotential, Upper,
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

/// `N` scalar fields evolving in an explicitly supplied field-space geometry.
///
/// This is separate from [`MultifieldSystem`] while the original canonical API
/// remains stable. `FlatFieldSpace` reproduces that system's dynamics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CovariantMultifieldSystem<P, G, const N: usize> {
    /// The scalar potential.
    pub potential: P,
    /// Metric, inverse metric, and connection provider.
    pub geometry: G,
    /// The reduced Planck mass.
    pub m_pl: Mass,
}

impl<P, G, const N: usize> CovariantMultifieldSystem<P, G, N>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    /// Builds a covariant multifield system.
    pub const fn new(potential: P, geometry: G, m_pl: Mass) -> Self {
        Self {
            potential,
            geometry,
            m_pl,
        }
    }

    /// Builds a covariant system in reduced Planck units.
    pub const fn in_planck_units(potential: P, geometry: G) -> Self {
        Self::new(potential, geometry, Mass::new(1.0))
    }

    /// Total covariant kinetic plus potential energy density.
    pub fn energy_density(&self, state: &MultifieldState<N>) -> EnergyDensity {
        let velocity = ContravariantField::new(state.pi);
        let metric = self.geometry.metric(&state.phi);
        let kinetic = metric.lower(&velocity).contract(&velocity) * 0.5;

        kinetic + self.potential.value(state.phi)
    }

    /// The expanding-branch Hubble rate.
    pub fn hubble(&self, state: &MultifieldState<N>) -> Hubble {
        (self.energy_density(state) / (self.m_pl * self.m_pl * 3.0)).sqrt()
    }

    /// `epsilon_H = G_IJ pi^I pi^J / (2 M_pl^2 H^2)`.
    pub fn epsilon_h(&self, state: &MultifieldState<N>) -> f64 {
        let velocity = ContravariantField::new(state.pi);
        let kinetic_twice = self
            .geometry
            .metric(&state.phi)
            .lower(&velocity)
            .contract(&velocity);
        let hubble = self.hubble(state);
        let denominator = self.m_pl * self.m_pl * hubble * hubble * 2.0;

        (kinetic_twice / denominator).into_value()
    }

    /// Slow-roll velocity `pi^I ~= -G^IJ V_;J / (3H)`.
    pub fn slow_roll_velocity(&self, state: &MultifieldState<N>) -> [FieldVelocity; N] {
        let gradient = crate::CovariantField::new(self.potential.gradient(state.phi));
        let raised = self.geometry.inverse_metric(&state.phi).raise(&gradient);
        let hubble = self.hubble(state);

        raised
            .into_components()
            .map(|component| -(component / (hubble * 3.0)))
    }
}

impl<P, G, const N: usize> CovariantMultifieldSystem<P, G, N>
where
    P: TwiceDifferentiableMultifieldPotential<N>,
    G: CurvedFieldSpaceGeometry<N>,
{
    /// Builds the mixed-index covariant mass matrix `C^I_J` for linear field
    /// perturbations in spatially flat gauge.
    pub fn perturbation_mass_matrix(
        &self,
        state: &MultifieldState<N>,
    ) -> FieldTensor2<Upper, Lower, M2, f64, N> {
        let metric = self.geometry.metric(&state.phi);
        let inverse = self.geometry.inverse_metric(&state.phi);
        let velocity = ContravariantField::new(state.pi);
        let lowered_velocity = metric.lower(&velocity);
        let covariant_hessian = self.covariant_potential_hessian(state);
        let raised_hessian: FieldTensor2<Upper, Lower, M2, f64, N> =
            FieldTensor2::new(core::array::from_fn(|upper| {
                core::array::from_fn(|lower| {
                    Quantity::new(
                        (0..N)
                            .map(|contracted| {
                                inverse.components()[upper][contracted].value()
                                    * covariant_hessian.components()[contracted][lower].value()
                            })
                            .sum(),
                    )
                })
            }));
        let curvature = self
            .geometry
            .riemann(&state.phi)
            .contract_velocities(&velocity);
        let gradient = crate::CovariantField::new(self.potential.gradient(state.phi));
        let raised_gradient = inverse.raise(&gradient).into_components();
        let hubble = self.hubble(state);
        let covariant_acceleration = ContravariantField::new(core::array::from_fn(|index| {
            state.pi[index] * hubble * -3.0 - raised_gradient[index]
        }));
        let lowered_acceleration = metric.lower(&covariant_acceleration);
        let epsilon = self.epsilon_h(state);
        let planck_squared = self.m_pl * self.m_pl;

        FieldTensor2::new(core::array::from_fn(|upper| {
            core::array::from_fn(|lower| {
                let gravitational =
                    (state.pi[upper] * lowered_velocity.components()[lower] * (3.0 + epsilon)
                        + (covariant_acceleration.components()[upper]
                            * lowered_velocity.components()[lower]
                            + state.pi[upper] * lowered_acceleration.components()[lower])
                            / hubble)
                        / planck_squared;

                raised_hessian.components()[upper][lower]
                    - curvature.components()[upper][lower]
                    - gravitational
            })
        }))
    }
}

impl<P, G, const N: usize> CovariantMultifieldSystem<P, G, N>
where
    P: TwiceDifferentiableMultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    /// Computes potential value, coordinate gradient, and coordinate Hessian.
    pub fn potential_derivatives(&self, state: &MultifieldState<N>) -> PotentialDerivatives<N> {
        self.potential.derivatives(state.phi)
    }

    /// Computes `V_;IJ = V_,IJ - Gamma^K_IJ V_,K`.
    pub fn covariant_potential_hessian(
        &self,
        state: &MultifieldState<N>,
    ) -> CovariantPotentialHessian<N> {
        let derivatives = self.potential.derivatives(state.phi);
        let connection = self.geometry.christoffel(&state.phi);

        CovariantPotentialHessian::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                let correction: f64 = (0..N)
                    .map(|upper| {
                        connection.components()[upper][first][second].value()
                            * derivatives.gradient[upper].value()
                    })
                    .sum();

                PotentialHessian::new(derivatives.hessian[first][second].into_value() - correction)
            })
        }))
    }
}

/// Two-field kinematics attached to a background state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoFieldKinematics {
    /// Instantaneous adiabatic and entropy basis.
    pub basis: KinematicBasis2,
    /// Signed angular velocity `e_sI D_t e_sigma^I`.
    pub turn_rate: Hubble,
    /// Dimensionless signed turn rate `turn_rate / H`.
    pub eta_perp: Quantity<M0, f64>,
}

impl<P, G> CovariantMultifieldSystem<P, G, 2>
where
    P: MultifieldPotential<2>,
    G: FieldSpaceGeometry<2>,
{
    /// Computes the oriented adiabatic/entropy basis and signed turn rate.
    pub fn kinematics(
        &self,
        state: &MultifieldState<2>,
        previous_normal: Option<&ContravariantField<M0, f64, 2>>,
    ) -> Result<TwoFieldKinematics, KinematicBasisError> {
        let metric = self.geometry.metric(&state.phi);
        let inverse = self.geometry.inverse_metric(&state.phi);
        let velocity = ContravariantField::new(state.pi);
        let basis = KinematicBasis2::from_velocity(&metric, &velocity, previous_normal)?;
        let gradient = crate::CovariantField::new(self.potential.gradient(state.phi));
        let raised_gradient = inverse.raise(&gradient).into_components();
        let hubble = self.hubble(state);
        let covariant_acceleration = ContravariantField::new(core::array::from_fn(|index| {
            state.pi[index] * hubble * -3.0 - raised_gradient[index]
        }));
        let normal_lower = metric.lower(&basis.normal);
        let turn_rate = normal_lower.contract(&covariant_acceleration) / basis.speed;

        Ok(TwoFieldKinematics {
            basis,
            turn_rate,
            eta_perp: turn_rate / hubble,
        })
    }
}

impl<P, G, const N: usize> System for CovariantMultifieldSystem<P, G, N>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    type IndependentVariable = Time;
    type State = MultifieldState<N>;
    type Delta = MultifieldStateDelta<N>;
    type Derivative = MultifieldStateDerivative<N>;

    fn derivative(&self, _time: Time, state: &Self::State) -> Self::Derivative {
        let hubble = self.hubble(state);
        let gradient = crate::CovariantField::new(self.potential.gradient(state.phi));
        let raised_gradient = self
            .geometry
            .inverse_metric(&state.phi)
            .raise(&gradient)
            .into_components();
        let velocity = ContravariantField::new(state.pi);
        let connection = self
            .geometry
            .christoffel(&state.phi)
            .quadratic_velocity(&velocity)
            .into_components();

        MultifieldStateDerivative {
            dphi_dt: state.pi,
            dpi_dt: core::array::from_fn(|index| {
                state.pi[index] * hubble * -3.0 - raised_gradient[index] - connection[index]
            }),
            da_dt: state.a * hubble,
        }
    }
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
    use crate::{
        ChristoffelSymbols, FieldMetric, FlatFieldSpace, InverseFieldMetric, MMinus1,
        QuadraticMultifieldPotential, rk4_step, rk45_step,
    };

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

    #[test]
    fn flat_covariant_system_reproduces_the_canonical_system() {
        let canonical = two_field_system();
        let covariant =
            CovariantMultifieldSystem::in_planck_units(canonical.potential, FlatFieldSpace);
        let state = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(-1.0e-6), FieldVelocity::new(-2.0e-6)],
            ScaleFactor::new(1.0),
        );

        assert_eq!(
            covariant.energy_density(&state),
            canonical.energy_density(&state)
        );
        assert_eq!(covariant.hubble(&state), canonical.hubble(&state));
        assert_eq!(
            covariant.derivative(Time::new(0.0), &state),
            canonical.derivative(Time::new(0.0), &state)
        );
        assert!((covariant.epsilon_h(&state) - canonical.epsilon_h(&state)).abs() < 1.0e-15);
    }

    #[derive(Clone, Copy)]
    struct ConnectedFlatSpace;

    impl FieldSpaceGeometry<2> for ConnectedFlatSpace {
        fn metric(&self, _position: &[Field; 2]) -> FieldMetric<2> {
            FieldMetric::euclidean()
        }

        fn inverse_metric(&self, _position: &[Field; 2]) -> InverseFieldMetric<2> {
            InverseFieldMetric::euclidean()
        }

        fn christoffel(&self, _position: &[Field; 2]) -> ChristoffelSymbols<2> {
            let mut components = [[[Quantity::<MMinus1, f64>::new(0.0); 2]; 2]; 2];
            components[0][1][1] = Quantity::new(2.0);
            ChristoffelSymbols::new(components)
        }
    }

    #[test]
    fn covariant_equation_includes_the_connection_force() {
        let potential = QuadraticMultifieldPotential::new([Mass::new(0.0), Mass::new(0.0)]);
        let system = CovariantMultifieldSystem::in_planck_units(potential, ConnectedFlatSpace);
        let state = MultifieldState::new(
            [Field::new(0.0), Field::new(0.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(1.0)],
            ScaleFactor::new(1.0),
        );

        let derivative = system.derivative(Time::new(0.0), &state);

        assert_eq!(derivative.dpi_dt[0].into_value(), -2.0);
    }

    #[test]
    fn covariant_hessian_includes_the_connection_correction() {
        let potential = QuadraticMultifieldPotential::new([Mass::new(2.0), Mass::new(3.0)]);
        let system = CovariantMultifieldSystem::in_planck_units(potential, ConnectedFlatSpace);
        let state = MultifieldState::new(
            [Field::new(3.0), Field::new(4.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );

        let hessian = system
            .covariant_potential_hessian(&state)
            .into_components()
            .map(|row| row.map(Quantity::into_value));

        assert_eq!(hessian, [[4.0, 0.0], [0.0, -15.0]]);
    }

    #[test]
    fn perturbation_mass_matrix_reduces_to_the_potential_hessian_at_rest() {
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(2.0), Mass::new(3.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(3.0), Field::new(4.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );

        let mass = system
            .perturbation_mass_matrix(&state)
            .into_components()
            .map(|row| row.map(Quantity::into_value));

        assert_eq!(mass, [[4.0, 0.0], [0.0, 9.0]]);
    }

    #[test]
    fn two_field_kinematics_reports_an_orthogonal_signed_turn() {
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(2.0e-6)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(-3.0e-6), FieldVelocity::new(-4.0e-6)],
            ScaleFactor::new(1.0),
        );

        let kinematics = system.kinematics(&state, None).unwrap();
        let metric = FieldMetric::<2>::euclidean();
        let overlap = metric
            .lower(&kinematics.basis.tangent)
            .contract(&kinematics.basis.normal)
            .into_value();

        assert!(overlap.abs() < 1.0e-15);
        assert!(kinematics.turn_rate.into_value().is_finite());
        assert!(kinematics.eta_perp.into_value().is_finite());
    }
}
