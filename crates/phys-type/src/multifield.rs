//! Homogeneous backgrounds containing a fixed number of scalar fields.

use core::ops::{Add, Mul};

use crate::integrate::{AdaptiveSystem, ErrorTolerances, System, component_error_ratio};
use crate::potential::MultifieldPotential;
use crate::{
    ContravariantField, CovariantPotentialHessian, CovariantPotentialThirdDerivative,
    CurvedFieldSpaceGeometry, EnergyDensity, Field, FieldSpaceGeometry, FieldTensor2,
    FieldVelocity, Hubble, KinematicBasis2, KinematicBasisError, Lower, M0, M1, M2, M3, Mass,
    PotentialDerivatives, PotentialDerivativesThroughThird, PotentialGradient, PotentialHessian,
    PotentialThirdDerivative, Quantity, ScaleFactor, ThirdOrderFieldSpaceGeometry,
    ThriceDifferentiableMultifieldPotential, Time, TwiceDifferentiableMultifieldPotential, Upper,
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

// A const-generic array is not `Default`, so this is written out.
impl<const N: usize> Default for MultifieldStateDelta<N> {
    fn default() -> Self {
        Self::from_planck_units(0.0)
    }
}

impl<const N: usize> MultifieldStateDelta<N> {
    /// Attaches the same raw Planck-unit number to every component.
    ///
    /// Convenient rather than safe: it assumes one number suits every
    /// dimension, which holds only while a state's components stay within a
    /// few orders of each other.
    pub fn from_planck_units(value: f64) -> Self {
        Self {
            phi: [Field::new(value); N],
            pi: [FieldVelocity::new(value); N],
            a: ScaleFactor::new(value),
        }
    }
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
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N>,
{
    /// Computes potential value and coordinate derivatives through third
    /// order.
    pub fn potential_derivatives_through_third(
        &self,
        state: &MultifieldState<N>,
    ) -> PotentialDerivativesThroughThird<N> {
        self.potential.derivatives_through_third(state.phi)
    }

    /// Computes the fully symmetrized covariant third derivative `V_;(IJK)`.
    ///
    /// The implementation first forms `nabla_K (V_;IJ)`. Covariant derivatives
    /// do not generally commute on the gradient in curved field space, so the
    /// three index positions are then explicitly symmetrized as required by
    /// the cubic interaction.
    pub fn covariant_potential_third_derivative(
        &self,
        state: &MultifieldState<N>,
    ) -> CovariantPotentialThirdDerivative<N> {
        let derivatives = self.potential.derivatives_through_third(state.phi);
        let connection = self.geometry.christoffel(&state.phi);
        let connection_derivative = self.geometry.christoffel_derivative(&state.phi);
        let covariant_hessian: CovariantPotentialHessian<N> =
            CovariantPotentialHessian::new(core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    let correction: f64 = (0..N)
                        .map(|upper| {
                            connection.components()[upper][first][second].value()
                                * derivatives.gradient[upper].value()
                        })
                        .sum();

                    PotentialHessian::new(derivatives.hessian[first][second].value() - correction)
                })
            }));

        let ordered: [[[PotentialThirdDerivative; N]; N]; N] = core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|derivative| {
                    let differentiated_connection: f64 = (0..N)
                        .map(|upper| {
                            connection_derivative.components()[upper][first][second][derivative]
                                .value()
                                * derivatives.gradient[upper].value()
                                + connection.components()[upper][first][second].value()
                                    * derivatives.hessian[upper][derivative].value()
                        })
                        .sum();
                    let tensor_connections: f64 = (0..N)
                        .map(|upper| {
                            connection.components()[upper][derivative][first].value()
                                * covariant_hessian.components()[upper][second].value()
                                + connection.components()[upper][derivative][second].value()
                                    * covariant_hessian.components()[first][upper].value()
                        })
                        .sum();

                    PotentialThirdDerivative::new(
                        derivatives.third[first][second][derivative].value()
                            - differentiated_connection
                            - tensor_connections,
                    )
                })
            })
        });

        CovariantPotentialThirdDerivative::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    PotentialThirdDerivative::new(
                        (ordered[first][second][third].value()
                            + ordered[first][third][second].value()
                            + ordered[second][third][first].value())
                            / 3.0,
                    )
                })
            })
        }))
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

/// Kinematics of a background state for any number of fields.
///
/// The tangent `e_sigma^I = phidot^I / sigmadot` is unique for every `N`. The
/// entropy directions form its orthogonal complement, of dimension `N - 1`,
/// which has no preferred basis for `N >= 3`; it is described by the
/// projector `perp^IJ = G^IJ - e_sigma^I e_sigma^J`. The turn vector is
/// `omega^I = D_t e_sigma^I = perp^I_J D_t phidot^J / sigmadot`, orthogonal to
/// the tangent, with `D_t phidot^I = -3 H phidot^I - G^IJ V_J`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Kinematics<const N: usize> {
    /// Unit tangent `e_sigma^I`.
    pub tangent: ContravariantField<M0, f64, N>,
    /// Covariant speed `sigmadot`.
    pub speed: FieldVelocity,
    /// Entropy projector `perp^IJ = G^IJ - e_sigma^I e_sigma^J`.
    pub entropy_projector: FieldTensor2<Upper, Upper, M0, f64, N>,
    /// Turn vector `omega^I = D_t e_sigma^I`.
    pub turn: ContravariantField<M1, f64, N>,
    /// Turn rate `|omega| = sqrt(G_IJ omega^I omega^J)`.
    pub turn_rate: Hubble,
    /// Dimensionless turn rate `|omega| / H`.
    pub eta_perp: Quantity<M0, f64>,
}

impl<P, G, const N: usize> CovariantMultifieldSystem<P, G, N>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    /// The tangent, entropy projector, and turn vector at `state`.
    pub fn kinematics(
        &self,
        state: &MultifieldState<N>,
    ) -> Result<Kinematics<N>, KinematicBasisError> {
        let metric = self.geometry.metric(&state.phi);
        let inverse = self.geometry.inverse_metric(&state.phi);
        let velocity = ContravariantField::new(state.pi);
        let speed_squared = metric.lower(&velocity).contract(&velocity).into_value();
        if speed_squared <= 0.0 || !speed_squared.is_finite() {
            return Err(KinematicBasisError::ZeroSpeed);
        }
        let speed = speed_squared.sqrt();
        let unit: [f64; N] = core::array::from_fn(|i| state.pi[i].into_value() / speed);
        let tangent: ContravariantField<M0, f64, N> =
            ContravariantField::new(unit.map(Quantity::new));
        let projector: [[f64; N]; N] = core::array::from_fn(|i| {
            core::array::from_fn(|j| inverse.components()[i][j].into_value() - unit[i] * unit[j])
        });
        let gradient = crate::CovariantField::new(self.potential.gradient(state.phi));
        let raised_gradient = inverse.raise(&gradient).into_components();
        let hubble = self.hubble(state).into_value();
        let acceleration: [f64; N] = core::array::from_fn(|i| {
            -3.0 * hubble * state.pi[i].into_value() - raised_gradient[i].into_value()
        });
        // perp^I_J a^J = perp^IK G_KJ a^J.
        let lowered: [f64; N] = core::array::from_fn(|k| {
            (0..N)
                .map(|j| metric.components()[k][j].into_value() * acceleration[j])
                .sum()
        });
        let turn: [f64; N] = core::array::from_fn(|i| {
            (0..N).map(|k| projector[i][k] * lowered[k]).sum::<f64>() / speed
        });
        let turn_rate = (0..N)
            .flat_map(|i| (0..N).map(move |j| (i, j)))
            .map(|(i, j)| metric.components()[i][j].into_value() * turn[i] * turn[j])
            .sum::<f64>()
            .max(0.0)
            .sqrt();
        Ok(Kinematics {
            tangent,
            speed: FieldVelocity::new(speed),
            entropy_projector: FieldTensor2::new(projector.map(|row| row.map(Quantity::new))),
            turn: ContravariantField::new(turn.map(Quantity::new)),
            turn_rate: Hubble::new(turn_rate),
            eta_perp: Quantity::new(turn_rate / hubble),
        })
    }
}

/// Oriented two-field kinematics attached to a background state: the
/// general [`Kinematics`] together with an entropy normal whose sign is
/// continued along the trajectory, and the signed turn rate.
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
    pub fn oriented_kinematics(
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

/// The largest normalized component error of a multifield background step.
///
/// Both background systems evolve the same state, so they share this norm.
pub(crate) fn multifield_error_ratio<const N: usize>(
    current: &MultifieldState<N>,
    next: &MultifieldState<N>,
    error: &MultifieldStateDelta<N>,
    tolerances: &ErrorTolerances<MultifieldStateDelta<N>>,
) -> f64 {
    let relative = tolerances.relative;
    let absolute = &tolerances.absolute;
    let fields = (0..N).map(|index| {
        component_error_ratio(
            error.phi[index],
            current.phi[index],
            next.phi[index],
            absolute.phi[index],
            relative,
        )
    });
    let velocities = (0..N).map(|index| {
        component_error_ratio(
            error.pi[index],
            current.pi[index],
            next.pi[index],
            absolute.pi[index],
            relative,
        )
    });
    let scale_factor = component_error_ratio(error.a, current.a, next.a, absolute.a, relative);

    fields.chain(velocities).fold(scale_factor, f64::max)
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

impl<P, G, const N: usize> AdaptiveSystem for CovariantMultifieldSystem<P, G, N>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    type AbsoluteTolerance = MultifieldStateDelta<N>;

    fn error_ratio(
        &self,
        current: &Self::State,
        next: &Self::State,
        error: &Self::Delta,
        tolerances: &ErrorTolerances<MultifieldStateDelta<N>>,
    ) -> f64 {
        multifield_error_ratio(current, next, error, tolerances)
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

impl<P, const N: usize> AdaptiveSystem for MultifieldSystem<P, N>
where
    P: MultifieldPotential<N>,
{
    type AbsoluteTolerance = MultifieldStateDelta<N>;

    fn error_ratio(
        &self,
        current: &Self::State,
        next: &Self::State,
        error: &Self::Delta,
        tolerances: &ErrorTolerances<MultifieldStateDelta<N>>,
    ) -> f64 {
        multifield_error_ratio(current, next, error, tolerances)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ChristoffelSymbols, FieldMetric, FlatFieldSpace, HyperbolicFieldSpace, InverseFieldMetric,
        MMinus1, QuadraticMultifieldPotential, rk4_step, rk45_step,
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
    fn flat_covariant_third_derivative_matches_coordinate_derivative() {
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(2.0), Mass::new(3.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(3.0), Field::new(4.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );

        let third = system
            .covariant_potential_third_derivative(&state)
            .into_components();

        assert_eq!(
            third.map(|matrix| matrix.map(|row| row.map(Quantity::into_value))),
            [[[0.0; 2]; 2]; 2]
        );
    }

    #[test]
    fn hyperbolic_covariant_third_derivative_matches_hessian_finite_difference() {
        let geometry = HyperbolicFieldSpace::new(Field::new(7.0)).unwrap();
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(2.0), Mass::new(3.0)]),
            geometry,
        );
        let state = MultifieldState::new(
            [Field::new(3.0), Field::new(1.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );
        let analytic = system.covariant_potential_third_derivative(&state);
        let central_hessian = system.covariant_potential_hessian(&state);
        let connection = geometry.christoffel(&state.phi);
        let step = 1.0e-5;

        let ordered: [[[f64; 2]; 2]; 2] = core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|derivative| {
                    let mut forward = state;
                    let mut backward = state;
                    forward.phi[derivative] = Field::new(state.phi[derivative].value() + step);
                    backward.phi[derivative] = Field::new(state.phi[derivative].value() - step);
                    let forward_hessian = system.covariant_potential_hessian(&forward);
                    let backward_hessian = system.covariant_potential_hessian(&backward);
                    let coordinate_derivative = (forward_hessian.components()[first][second]
                        .value()
                        - backward_hessian.components()[first][second].value())
                        / (2.0 * step);
                    let correction: f64 = (0..2)
                        .map(|upper| {
                            connection.components()[upper][derivative][first].value()
                                * central_hessian.components()[upper][second].value()
                                + connection.components()[upper][derivative][second].value()
                                    * central_hessian.components()[first][upper].value()
                        })
                        .sum();

                    coordinate_derivative - correction
                })
            })
        });

        for (first, matrix) in ordered.iter().enumerate() {
            for (second, row) in matrix.iter().enumerate() {
                for (third, _) in row.iter().enumerate() {
                    let numerical = (ordered[first][second][third]
                        + ordered[first][third][second]
                        + ordered[second][third][first])
                        / 3.0;
                    let expected = analytic.components()[first][second][third].value();

                    assert!((numerical - expected).abs() < 1.0e-8);
                    assert_eq!(
                        analytic.components()[first][second][third],
                        analytic.components()[third][second][first]
                    );
                    assert_eq!(
                        analytic.components()[first][second][third],
                        analytic.components()[first][third][second]
                    );
                }
            }
        }
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

        let kinematics = system.oriented_kinematics(&state, None).unwrap();
        let metric = FieldMetric::<2>::euclidean();
        let overlap = metric
            .lower(&kinematics.basis.tangent)
            .contract(&kinematics.basis.normal)
            .into_value();

        assert!(overlap.abs() < 1.0e-15);
        assert!(kinematics.turn_rate.into_value().is_finite());
        assert!(kinematics.eta_perp.into_value().is_finite());
    }
    /// The flat three-field potential `sum m_i^2 phi_i^2 / 2`.
    fn three_field_system()
    -> CovariantMultifieldSystem<QuadraticMultifieldPotential<3>, FlatFieldSpace, 3> {
        CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([
                Mass::new(1.0e-6),
                Mass::new(2.3e-6),
                Mass::new(4.1e-6),
            ]),
            FlatFieldSpace,
        )
    }

    fn slow_roll_state<P, G, const N: usize>(
        system: &CovariantMultifieldSystem<P, G, N>,
        phi: [f64; N],
    ) -> MultifieldState<N>
    where
        P: MultifieldPotential<N>,
        G: FieldSpaceGeometry<N>,
    {
        let resting = MultifieldState::new(
            phi.map(Field::new),
            [FieldVelocity::new(0.0); N],
            ScaleFactor::new(1.0),
        );
        // Tilt the slow-roll velocity off the gradient, so that the state has
        // a generic acceleration and a nonzero turn.
        let slow_roll = system.slow_roll_velocity(&resting);
        let scale = slow_roll
            .iter()
            .map(|v| v.into_value().abs())
            .fold(0.0, f64::max);
        MultifieldState {
            pi: core::array::from_fn(|i| {
                FieldVelocity::new(slow_roll[i].into_value() + 0.3 * scale * (i as f64 - 0.7))
            }),
            ..resting
        }
    }

    #[test]
    fn n_field_kinematics_reduce_to_the_oriented_two_field_basis() {
        let system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(7.0e-6)]),
            HyperbolicFieldSpace::new(Field::new(5.0)).unwrap(),
        );
        let state = slow_roll_state(&system, [12.0, 8.0]);
        let general = system.kinematics(&state).unwrap();
        let oriented = system.oriented_kinematics(&state, None).unwrap();
        let close = |a: f64, b: f64| (a - b).abs() <= 1.0e-12 * (1.0 + a.abs().max(b.abs()));
        assert!(close(
            general.speed.into_value(),
            oriented.basis.speed.into_value()
        ));
        assert!(close(
            general.turn_rate.into_value(),
            oriented.turn_rate.into_value().abs()
        ));
        for i in 0..2 {
            assert!(close(
                general.tangent.components()[i].into_value(),
                oriented.basis.tangent.components()[i].into_value()
            ));
            // omega = (signed turn rate) e_s.
            assert!(close(
                general.turn.components()[i].into_value(),
                oriented.turn_rate.into_value()
                    * oriented.basis.normal.components()[i].into_value()
            ));
            for j in 0..2 {
                // perp^IJ = e_s^I e_s^J in two dimensions.
                let normal = oriented.basis.normal.components();
                assert!(close(
                    general.entropy_projector.components()[i][j].into_value(),
                    normal[i].into_value() * normal[j].into_value()
                ));
            }
        }
    }

    #[test]
    fn the_entropy_projector_is_an_orthogonal_projector_of_rank_n_minus_one() {
        let system = three_field_system();
        let state = slow_roll_state(&system, [10.0, 7.0, 4.0]);
        let kinematics = system.kinematics(&state).unwrap();
        let p = |i: usize, j: usize| kinematics.entropy_projector.components()[i][j].into_value();
        let e = |i: usize| kinematics.tangent.components()[i].into_value();
        let w = |i: usize| kinematics.turn.components()[i].into_value();
        // Flat metric: perp perp = perp, perp e = 0, trace = N - 1, omega . e = 0.
        for i in 0..3 {
            assert!((0..3).map(|k| p(i, k) * e(k)).sum::<f64>().abs() < 1.0e-14);
            for j in 0..3 {
                let square: f64 = (0..3).map(|k| p(i, k) * p(k, j)).sum();
                assert!((square - p(i, j)).abs() < 1.0e-14);
            }
        }
        assert!(((0..3).map(|i| p(i, i)).sum::<f64>() - 2.0).abs() < 1.0e-14);
        let overlap: f64 = (0..3).map(|i| w(i) * e(i)).sum();
        assert!(overlap.abs() < 1.0e-12 * kinematics.turn_rate.into_value());
        assert!(kinematics.turn_rate.into_value() > 0.0);
    }

    /// A curved, non-diagonal three-field metric.
    #[derive(Clone, Copy, Debug)]
    struct CurvedThreeField;

    impl crate::MetricJet<3> for CurvedThreeField {
        fn metric_jet(&self, position: &[Field; 3]) -> [[crate::ThirdOrderJet<3>; 3]; 3] {
            use crate::{ThirdOrderJet, Transcendental};
            let [x, y, z]: [ThirdOrderJet<3>; 3] =
                core::array::from_fn(|i| ThirdOrderJet::variable(position[i].into_value(), i));
            let one = ThirdOrderJet::constant(1.0);
            let a = one + x * x * 0.01;
            let b = (y * 0.05).exp();
            let c = one + z * z * 0.02;
            let ab = (x * y * 0.001).tanh();
            let ac = z * 0.01;
            let bc = one * 0.0;
            [[a, ab, ac], [ab, b, bc], [ac, bc, c]]
        }
    }

    /// The numerical `D_t e_sigma^I` from a short integration both ways.
    fn assert_turn_matches_finite_difference<P, G>(
        system: &CovariantMultifieldSystem<P, G, 3>,
        state: MultifieldState<3>,
    ) where
        P: TwiceDifferentiableMultifieldPotential<3>,
        G: CurvedFieldSpaceGeometry<3>,
    {
        let kinematics = system.kinematics(&state).unwrap();
        let step = 2.0e2;
        let advance = |dt: f64| rk4_step(system, Time::new(0.0), &state, Time::new(dt));
        let tangent_at = |s: &MultifieldState<3>| system.kinematics(s).unwrap().tangent;
        let (plus, minus) = (tangent_at(&advance(step)), tangent_at(&advance(-step)));
        let connection = system.geometry.christoffel(&state.phi);
        for i in 0..3 {
            let derivative = (plus.components()[i].into_value()
                - minus.components()[i].into_value())
                / (2.0 * step);
            // D_t e^I = d e^I / dt + Gamma^I_JK e^J phidot^K.
            let correction: f64 = (0..9)
                .map(|m| {
                    let [j, k] = [m / 3, m % 3];
                    connection.components()[i][j][k].into_value()
                        * kinematics.tangent.components()[j].into_value()
                        * state.pi[k].into_value()
                })
                .sum();
            let analytic = kinematics.turn.components()[i].into_value();
            assert!(
                (derivative + correction - analytic).abs()
                    < 1.0e-4 * kinematics.turn_rate.into_value(),
                "{i}: {} vs {analytic}",
                derivative + correction
            );
        }
    }

    #[test]
    fn the_turn_vector_is_the_covariant_time_derivative_of_the_tangent() {
        let flat = three_field_system();
        assert_turn_matches_finite_difference(&flat, slow_roll_state(&flat, [10.0, 7.0, 4.0]));
        let curved = CovariantMultifieldSystem::in_planck_units(
            flat.potential,
            crate::JetGeometry(CurvedThreeField),
        );
        assert_turn_matches_finite_difference(&curved, slow_roll_state(&curved, [10.0, 7.0, 4.0]));
    }
}
