//! Theories stated by their action, and their computational representations.
//!
//! A [`GeneralizedNonminimal`] theory is given by its Jordan-frame action
//!
//! ```text
//! S = int d^4x sqrt(-g) [ (M_pl^2 / 2) Omega^2(phi) R
//!                         - (1/2) k_IJ(phi) g^{mu nu} d_mu phi^I d_nu phi^J
//!                         - V_J(phi) ],
//! ```
//!
//! in the metric or the Palatini formulation. The engine evolves the
//! Einstein-frame sigma model, so [`einstein_frame`] maps the theory exactly
//! to the Einstein-frame field-space metric and potential,
//!
//! ```text
//! G_IJ = k_IJ / Omega^2 + kappa (3/2) M_pl^2 (d_I Omega^2)(d_J Omega^2) / Omega^4,
//! V    = V_J / Omega^4,
//! ```
//!
//! with `kappa = 1` for the metric and `kappa = 0` for the Palatini
//! formulation, under `g^E = Omega^2 g^J`. The Einstein frame is a
//! computational representation, not a statement about which frame is
//! physical; [`EinsteinFrame::correspondence`] returns the Jordan-frame scale
//! factor, e-fold number, Hubble rate, and field velocities along an
//! Einstein-frame state. The derivation and the frame dependence of the
//! curvature perturbation are recorded in `docs/frames.md`.
//!
//! The geometry is a [`JetGeometry`] of [`EinsteinMetric`]: the conformal
//! factor is differentiated four times and the kinetic matrix three times by
//! automatic differentiation, so a theory only states its three functions.

use core::ops::{Add, Div, Mul, Sub};

use crate::{
    Dimensionless, EnergyDensity, Field, FourthOrderJet, GeometryDomainError, Jet, JetGeometry,
    Mass, MetricJet, MultifieldPotentialValue, MultifieldState, Quantity, Sqrt, ThirdOrderJet,
    Transcendental,
};

/// A scalar type with which a theory's functions can be evaluated: `f64` and
/// every forward-mode jet of the crate.
pub trait TheoryScalar:
    Copy
    + From<f64>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Mul<f64, Output = Self>
    + Transcendental
    + Sqrt
{
}

impl<T> TheoryScalar for T where
    T: Copy
        + From<f64>
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Mul<f64, Output = T>
        + Transcendental
        + Sqrt
{
}

/// The variational formulation of the gravitational action.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Formulation {
    /// The metric (second-order) formulation, `kappa = 1`.
    Metric,
    /// The Palatini (first-order) formulation, `kappa = 0`.
    Palatini,
}

impl Formulation {
    /// `kappa` in the Einstein-frame metric.
    pub const fn kappa(self) -> f64 {
        match self {
            Self::Metric => 1.0,
            Self::Palatini => 0.0,
        }
    }
}

/// A scalar-tensor theory with nonminimal couplings `Omega^2(phi) R`,
/// a field-dependent kinetic matrix, and a potential, in the Jordan frame.
///
/// Every function is written once, generically over [`TheoryScalar`], and is
/// evaluated with `f64` and with jets. Fields are in reduced Planck units.
pub trait GeneralizedNonminimal<const N: usize> {
    /// The formulation of the gravitational action.
    fn formulation(&self) -> Formulation;

    /// The reduced Planck mass.
    fn planck_mass(&self) -> Mass {
        Mass::new(1.0)
    }

    /// `Omega^2(phi) = 2 f(phi) / M_pl^2`, positive where the theory is
    /// defined.
    fn conformal_factor<V: TheoryScalar>(&self, fields: [Field<V>; N]) -> Dimensionless<V>;

    /// The symmetric kinetic matrix `k_IJ(phi)`.
    fn kinetic_matrix<V: TheoryScalar>(&self, fields: [Field<V>; N]) -> [[Dimensionless<V>; N]; N];

    /// The Jordan-frame potential `V_J(phi)`.
    fn jordan_potential<V: TheoryScalar>(&self, fields: [Field<V>; N]) -> EnergyDensity<V>;
}

/// The Einstein-frame field-space metric of a theory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EinsteinMetric<T>(pub T);

impl<T> EinsteinMetric<T> {
    /// `kappa (3/2) M_pl^2`, the weight of `d_I Omega^2 d_J Omega^2 / Omega^4`.
    fn gradient_weight<const N: usize>(&self) -> f64
    where
        T: GeneralizedNonminimal<N>,
    {
        let planck = self.0.planck_mass().into_value();
        self.0.formulation().kappa() * 1.5 * planck * planck
    }
}

/// The Einstein-frame potential `V_J / Omega^4` of a theory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EinsteinPotential<T>(pub T);

/// The Einstein-frame geometry of a theory.
pub type EinsteinGeometry<T> = JetGeometry<EinsteinMetric<T>>;

/// The Einstein-frame potential and geometry of `theory`, ready for
/// `CovariantMultifieldSystem::in_planck_units`.
pub fn einstein_frame<T: Clone>(theory: T) -> (EinsteinPotential<T>, EinsteinGeometry<T>) {
    (
        EinsteinPotential(theory.clone()),
        JetGeometry(EinsteinMetric(theory)),
    )
}

fn active<V, const N: usize>(
    position: &[Field; N],
    variable: impl Fn(f64, usize) -> V,
) -> [Field<V>; N] {
    core::array::from_fn(|i| Quantity::new(variable(position[i].into_value(), i)))
}

fn check_conformal_factor<T: GeneralizedNonminimal<N>, const N: usize>(
    theory: &T,
    position: &[Field; N],
) -> Result<(), GeometryDomainError> {
    let omega_squared = theory.conformal_factor(*position).into_value();
    if omega_squared.is_finite() && omega_squared > 0.0 {
        Ok(())
    } else {
        Err(GeometryDomainError {
            reason: "the conformal factor Omega^2 is not positive and finite",
            value: omega_squared,
        })
    }
}

impl<T: GeneralizedNonminimal<N>, const N: usize> MetricJet<N> for EinsteinMetric<T> {
    fn check_domain(&self, position: &[Field; N]) -> Result<(), GeometryDomainError> {
        check_conformal_factor(&self.0, position)
    }

    fn metric_values(&self, position: &[Field; N]) -> [[f64; N]; N] {
        // The values need Omega^2 through first order only.
        let theory = &self.0;
        let omega_squared = theory
            .conformal_factor(active(position, Jet::<N>::variable))
            .into_value();
        let kinetic = theory.kinetic_matrix(*position);
        let weight = self.gradient_weight();
        let inverse = 1.0 / omega_squared.value;
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                kinetic[i][j].into_value() * inverse
                    + omega_squared.gradient[i]
                        * omega_squared.gradient[j]
                        * inverse
                        * inverse
                        * weight
            })
        })
    }

    fn metric_second(&self, position: &[Field; N]) -> [[Jet<N>; N]; N] {
        // Second derivatives of G need Omega^2 through third order.
        let theory = &self.0;
        let omega_squared = theory
            .conformal_factor(active(position, ThirdOrderJet::<N>::variable))
            .into_value();
        let value = Jet {
            value: omega_squared.value,
            gradient: omega_squared.gradient,
            hessian: omega_squared.hessian,
        };
        let slope: [Jet<N>; N] = core::array::from_fn(|i| Jet {
            value: omega_squared.gradient[i],
            gradient: omega_squared.hessian[i],
            hessian: omega_squared.third[i],
        });
        let kinetic = theory.kinetic_matrix(active(position, Jet::<N>::variable));
        let weight = self.gradient_weight();
        let inverse = Jet::constant(1.0) / value;
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                kinetic[i][j].into_value() * inverse
                    + slope[i] * slope[j] * inverse * inverse * weight
            })
        })
    }

    fn metric_jet(&self, position: &[Field; N]) -> [[ThirdOrderJet<N>; N]; N] {
        let theory = &self.0;
        let omega_squared = theory
            .conformal_factor(active(position, FourthOrderJet::<N>::variable))
            .into_value();
        // Omega^2 and d_I Omega^2 as third-order jets.
        let value = ThirdOrderJet {
            value: omega_squared.value,
            gradient: omega_squared.gradient,
            hessian: omega_squared.hessian,
            third: omega_squared.third,
        };
        let slope: [ThirdOrderJet<N>; N] = core::array::from_fn(|i| ThirdOrderJet {
            value: omega_squared.gradient[i],
            gradient: omega_squared.hessian[i],
            hessian: omega_squared.third[i],
            third: omega_squared.fourth[i],
        });
        let kinetic = theory.kinetic_matrix(active(position, ThirdOrderJet::<N>::variable));
        let gradient_weight = self.gradient_weight();
        let inverse = ThirdOrderJet::constant(1.0) / value;
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                kinetic[i][j].into_value() * inverse
                    + slope[i] * slope[j] * inverse * inverse * gradient_weight
            })
        })
    }
}

impl<T, V, const N: usize> MultifieldPotentialValue<V, N> for EinsteinPotential<T>
where
    T: GeneralizedNonminimal<N>,
    V: TheoryScalar,
{
    fn value(&self, fields: [Field<V>; N]) -> EnergyDensity<V> {
        let omega_squared = self.0.conformal_factor(fields).into_value();
        let jordan = self.0.jordan_potential(fields).into_value();
        EnergyDensity::new(jordan / (omega_squared * omega_squared))
    }
}

/// Jordan-frame quantities along an Einstein-frame state, under
/// `g^E = Omega^2 g^J`: `a_E = Omega a_J` and `dt_E = Omega dt_J`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameCorrespondence<const N: usize> {
    /// The conformal factor `Omega`.
    pub omega: f64,
    /// The Jordan-frame scale factor `a_J = a_E / Omega`.
    pub jordan_scale_factor: f64,
    /// `N_J = ln a_J = N_E - ln Omega`, with `N_E = ln a_E`.
    pub jordan_e_folds: f64,
    /// `H_J = Omega (H_E - phi'^I d_I ln Omega)`, the dot an Einstein-time
    /// derivative.
    pub jordan_hubble: f64,
    /// `d phi^I / dt_J = Omega d phi^I / dt_E`.
    pub jordan_field_velocity: [f64; N],
}

/// A theory together with its Einstein-frame representation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EinsteinFrame<T>(pub T);

impl<T> EinsteinFrame<T> {
    /// Checks that the map to the Einstein frame is defined at `position`.
    pub fn check<const N: usize>(&self, position: &[Field; N]) -> Result<(), GeometryDomainError>
    where
        T: GeneralizedNonminimal<N>,
    {
        check_conformal_factor(&self.0, position)
    }

    /// Jordan-frame quantities at an Einstein-frame state with Hubble rate
    /// `hubble` (Einstein time).
    pub fn correspondence<const N: usize>(
        &self,
        state: &MultifieldState<N>,
        hubble: crate::Hubble,
    ) -> Result<FrameCorrespondence<N>, GeometryDomainError>
    where
        T: GeneralizedNonminimal<N>,
    {
        self.check(&state.phi)?;
        let jet = self
            .0
            .conformal_factor(active(&state.phi, FourthOrderJet::<N>::variable))
            .into_value();
        let omega = jet.value.sqrt();
        // d_I ln Omega = d_I Omega^2 / (2 Omega^2).
        let rate: f64 = (0..N)
            .map(|i| jet.gradient[i] / (2.0 * jet.value) * state.pi[i].into_value())
            .sum();
        let a_e = state.a.into_value();
        Ok(FrameCorrespondence {
            omega,
            jordan_scale_factor: a_e / omega,
            jordan_e_folds: a_e.ln() - omega.ln(),
            jordan_hubble: omega * (hubble.into_value() - rate),
            jordan_field_velocity: state.pi.map(|v| omega * v.into_value()),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_models::{NonminimalAssistedFieldSpace, NonminimalLoopInflationPotential};
    use crate::{
        CurvatureDerivativeGeometry, CurvedFieldSpaceGeometry, FieldSpaceGeometry, ScaleFactor,
        ThirdOrderFieldSpaceGeometry, ThriceDifferentiableMultifieldPotential,
    };

    /// The nonminimally assisted model of arXiv:2302.05866 as a theory:
    /// `Omega^2 = 1 + xi s^2`, canonical Jordan kinetic terms, loop inflation.
    #[derive(Clone, Copy, Debug)]
    struct AssistedLoop {
        xi: f64,
        alpha: f64,
        formulation: Formulation,
    }

    impl GeneralizedNonminimal<2> for AssistedLoop {
        fn formulation(&self) -> Formulation {
            self.formulation
        }

        fn conformal_factor<V: TheoryScalar>(&self, fields: [Field<V>; 2]) -> Dimensionless<V> {
            let [_, s] = fields.map(Quantity::into_value);
            Quantity::new(V::from(1.0) + s * s * self.xi)
        }

        fn kinetic_matrix<V: TheoryScalar>(
            &self,
            _fields: [Field<V>; 2],
        ) -> [[Dimensionless<V>; 2]; 2] {
            core::array::from_fn(|i| {
                core::array::from_fn(|j| Quantity::new(V::from(f64::from(i == j))))
            })
        }

        fn jordan_potential<V: TheoryScalar>(&self, fields: [Field<V>; 2]) -> EnergyDensity<V> {
            let [phi, _] = fields.map(Quantity::into_value);
            EnergyDensity::new((V::from(1.0) + phi.ln() * self.alpha) * 1.0e-10)
        }
    }

    /// Single-field Higgs-like inflation, `Omega^2 = 1 + xi phi^2`.
    #[derive(Clone, Copy, Debug)]
    struct Higgs {
        xi: f64,
        formulation: Formulation,
    }

    impl GeneralizedNonminimal<1> for Higgs {
        fn formulation(&self) -> Formulation {
            self.formulation
        }

        fn conformal_factor<V: TheoryScalar>(&self, fields: [Field<V>; 1]) -> Dimensionless<V> {
            let phi = fields[0].into_value();
            Quantity::new(V::from(1.0) + phi * phi * self.xi)
        }

        fn kinetic_matrix<V: TheoryScalar>(
            &self,
            _fields: [Field<V>; 1],
        ) -> [[Dimensionless<V>; 1]; 1] {
            [[Quantity::new(V::from(1.0))]]
        }

        fn jordan_potential<V: TheoryScalar>(&self, fields: [Field<V>; 1]) -> EnergyDensity<V> {
            let phi = fields[0].into_value();
            EnergyDensity::new(phi * phi * phi * phi * 0.25e-12)
        }
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1.0e-12 * (1.0 + a.abs().max(b.abs()))
    }

    #[test]
    fn the_einstein_frame_reproduces_the_hand_written_assisted_model() {
        for formulation in [Formulation::Metric, Formulation::Palatini] {
            let theory = AssistedLoop {
                xi: 0.022,
                alpha: 0.1,
                formulation,
            };
            let (potential, geometry) = einstein_frame(theory);
            let reference_geometry =
                NonminimalAssistedFieldSpace::new(theory.xi, formulation.kappa(), Mass::new(1.0))
                    .unwrap();
            let reference_potential = NonminimalLoopInflationPotential {
                scale: EnergyDensity::new(1.0e-10),
                alpha: theory.alpha,
                coupling: theory.xi,
                planck_mass: Mass::new(1.0),
            };
            for (phi, s) in [(5.0, 1.0e-3), (2.0, 0.4), (0.8, 2.3)] {
                let position = [Field::new(phi), Field::new(s)];
                let a = potential.derivatives_through_third(position);
                let b = reference_potential.derivatives_through_third(position);
                assert!(close(a.value.into_value(), b.value.into_value()));
                for i in 0..2 {
                    assert!(close(
                        a.gradient[i].into_value(),
                        b.gradient[i].into_value()
                    ));
                    for j in 0..2 {
                        assert!(close(
                            a.hessian[i][j].into_value(),
                            b.hessian[i][j].into_value()
                        ));
                        for k in 0..2 {
                            assert!(close(
                                a.third[i][j][k].into_value(),
                                b.third[i][j][k].into_value()
                            ));
                        }
                    }
                }
                let (g, r) = (
                    geometry.metric(&position),
                    reference_geometry.metric(&position),
                );
                let (dg, dr) = (
                    geometry.christoffel_derivative(&position),
                    reference_geometry.christoffel_derivative(&position),
                );
                let (rg, rr) = (
                    geometry.riemann(&position),
                    reference_geometry.riemann(&position),
                );
                let (ng, nr) = (
                    geometry.covariant_riemann_derivative(&position),
                    reference_geometry.covariant_riemann_derivative(&position),
                );
                for i in 0..2 {
                    for j in 0..2 {
                        assert!(close(
                            g.components()[i][j].into_value(),
                            r.components()[i][j].into_value()
                        ));
                        for k in 0..2 {
                            for l in 0..2 {
                                assert!(close(
                                    dg.components()[i][j][k][l].into_value(),
                                    dr.components()[i][j][k][l].into_value()
                                ));
                                assert!(close(
                                    rg.components()[i][j][k][l].into_value(),
                                    rr.components()[i][j][k][l].into_value()
                                ));
                                for m in 0..2 {
                                    assert!(close(
                                        ng.components()[i][j][k][l][m].into_value(),
                                        nr.components()[i][j][k][l][m].into_value()
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn the_single_field_einstein_metric_has_its_closed_form() {
        // G = (1 + xi phi^2 + 6 kappa xi^2 phi^2) / (1 + xi phi^2)^2.
        let (xi, phi) = (1.0e3, 0.02);
        for formulation in [Formulation::Metric, Formulation::Palatini] {
            let (_, geometry) = einstein_frame(Higgs { xi, formulation });
            let metric = geometry.metric(&[Field::new(phi)]).components()[0][0].into_value();
            let omega2 = 1.0 + xi * phi * phi;
            let expected =
                (omega2 + 6.0 * formulation.kappa() * xi * xi * phi * phi) / (omega2 * omega2);
            assert!(
                close(metric, expected),
                "{formulation:?}: {metric} vs {expected}"
            );
        }
    }

    #[test]
    fn the_einstein_frame_rejects_a_vanishing_conformal_factor() {
        let theory = Higgs {
            xi: -1.0,
            formulation: Formulation::Metric,
        };
        let (_, geometry) = einstein_frame(theory);
        assert!(geometry.check_domain(&[Field::new(0.5)]).is_ok());
        let error = geometry.check_domain(&[Field::new(1.5)]).unwrap_err();
        assert!(error.value < 0.0, "{error}");
        assert!(EinsteinFrame(theory).check(&[Field::new(1.5)]).is_err());
    }

    #[test]
    fn the_frame_correspondence_is_consistent_along_a_path() {
        // The relations are kinematic, so any path phi(t_E), a_E(t_E) tests
        // them: d N_J / d t_J must equal the reported H_J.
        let frame = EinsteinFrame(AssistedLoop {
            xi: 0.3,
            alpha: 0.1,
            formulation: Formulation::Metric,
        });
        let hubble = 2.0e-5;
        let velocity = [-3.0e-6, 4.0e-6];
        let state_at = |t: f64| {
            MultifieldState::new(
                [
                    Field::new(2.0 + velocity[0] * t),
                    Field::new(0.7 + velocity[1] * t),
                ],
                velocity.map(crate::FieldVelocity::new),
                ScaleFactor::new((hubble * t).exp()),
            )
        };
        let (t, step) = (3.0e4, 10.0);
        let at = |t: f64| {
            frame
                .correspondence(&state_at(t), crate::Hubble::new(hubble))
                .unwrap()
        };
        let middle = at(t);
        let derivative_in_einstein_time =
            (at(t + step).jordan_e_folds - at(t - step).jordan_e_folds) / (2.0 * step);
        // dt_J = dt_E / Omega.
        let numerical = derivative_in_einstein_time * middle.omega;
        assert!(
            (numerical - middle.jordan_hubble).abs() < 1.0e-8 * middle.jordan_hubble.abs(),
            "{numerical} vs {}",
            middle.jordan_hubble
        );
        let state = state_at(t);
        assert!(close(
            middle.jordan_scale_factor * middle.omega,
            state.a.into_value()
        ));
        assert!(close(
            middle.jordan_field_velocity[1],
            middle.omega * velocity[1]
        ));
    }
}
