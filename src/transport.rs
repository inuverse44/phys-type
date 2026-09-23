//! Typed kinematics and correlation containers for phase-space transport.
//!
//! This module separates Fourier legs, phase-space slots, field coordinates,
//! and physical dimensions. It provides the linear mode solver's two-point
//! oracle and direct two-point covariance dynamics; nonlinear three-point
//! evolution is not implemented yet.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Mul};

use num_complex::Complex64;

use crate::integrate::System;
use crate::{
    CovariantMultifieldSystem, CurvedFieldSpaceGeometry, Dimensionless, Hubble, M0, M1, MMinus1,
    Mass, ModeAmplitude, ModeVelocity, MultifieldState, Quantity,
    TwiceDifferentiableMultifieldPotential, TwoFieldPerturbationState,
};

/// Default relative tolerance for numerical triangle closure.
pub const DEFAULT_TRIANGLE_CLOSURE_TOLERANCE: f64 = 64.0 * f64::EPSILON;

/// A Fourier leg of a three-point function.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum FourierLeg {
    /// The leg carrying `k_1`.
    One,
    /// The leg carrying `k_2`.
    Two,
    /// The leg carrying `k_3`.
    Three,
}

impl FourierLeg {
    /// Zero-based position of this leg.
    pub const fn index(self) -> usize {
        match self {
            Self::One => 0,
            Self::Two => 1,
            Self::Three => 2,
        }
    }

    const fn from_index(index: usize) -> Self {
        match index {
            0 => Self::One,
            1 => Self::Two,
            2 => Self::Three,
            _ => panic!("a Fourier leg index is always below three"),
        }
    }
}

/// Why three wavenumber magnitudes do not form a usable momentum triangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MomentumTriangleError {
    /// A leg was zero, negative, or non-finite.
    InvalidWavenumber {
        /// The offending input leg.
        leg: FourierLeg,
        /// The supplied normalized wavenumber.
        value: f64,
    },
    /// The requested relative closure tolerance was invalid.
    InvalidClosureTolerance {
        /// The supplied tolerance.
        value: f64,
    },
    /// The longest leg exceeded the sum of the other two beyond tolerance.
    DoesNotClose {
        /// Longest supplied magnitude.
        longest: f64,
        /// Sum of the two shorter magnitudes.
        shorter_sum: f64,
        /// Relative tolerance used in the comparison.
        relative_tolerance: f64,
    },
}

impl fmt::Display for MomentumTriangleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidWavenumber { leg, value } => write!(
                formatter,
                "momentum-triangle leg {} must be positive and finite, got {value:e}",
                leg.index() + 1
            ),
            Self::InvalidClosureTolerance { value } => write!(
                formatter,
                "triangle closure tolerance must be finite and non-negative, got {value:e}"
            ),
            Self::DoesNotClose {
                longest,
                shorter_sum,
                relative_tolerance,
            } => write!(
                formatter,
                "momentum triangle does not close: {longest:e} exceeds {shorter_sum:e} with relative tolerance {relative_tolerance:e}"
            ),
        }
    }
}

impl std::error::Error for MomentumTriangleError {}

/// Three positive normalized wavenumbers in canonical descending order.
///
/// The original caller order is retained as a permutation. This allows one
/// canonical triangle integration to serve all permutations without
/// confusing Fourier legs with field-space tensor indices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MomentumTriangle {
    perimeter: f64,
    second_over_first: f64,
    third_over_first: f64,
    canonical_to_input: [usize; 3],
    input_to_canonical: [usize; 3],
    closure_tolerance: f64,
}

impl MomentumTriangle {
    /// Builds a triangle using [`DEFAULT_TRIANGLE_CLOSURE_TOLERANCE`].
    pub fn try_new(k1: f64, k2: f64, k3: f64) -> Result<Self, MomentumTriangleError> {
        Self::try_new_with_tolerance(k1, k2, k3, DEFAULT_TRIANGLE_CLOSURE_TOLERANCE)
    }

    /// Builds a triangle with an explicit relative closure tolerance.
    pub fn try_new_with_tolerance(
        k1: f64,
        k2: f64,
        k3: f64,
        closure_tolerance: f64,
    ) -> Result<Self, MomentumTriangleError> {
        let input = [k1, k2, k3];
        for (index, value) in input.into_iter().enumerate() {
            if !value.is_finite() || value <= 0.0 {
                return Err(MomentumTriangleError::InvalidWavenumber {
                    leg: FourierLeg::from_index(index),
                    value,
                });
            }
        }
        if !closure_tolerance.is_finite() || closure_tolerance < 0.0 {
            return Err(MomentumTriangleError::InvalidClosureTolerance {
                value: closure_tolerance,
            });
        }

        let mut indexed = [(k1, 0_usize), (k2, 1_usize), (k3, 2_usize)];
        indexed.sort_by(|left, right| right.0.total_cmp(&left.0));
        let canonical = indexed.map(|(value, _)| value);
        let canonical_to_input = indexed.map(|(_, index)| index);
        let mut input_to_canonical = [0_usize; 3];
        for (canonical_index, input_index) in canonical_to_input.into_iter().enumerate() {
            input_to_canonical[input_index] = canonical_index;
        }

        let shorter_sum = canonical[1] + canonical[2];
        if canonical[0] > shorter_sum * (1.0 + closure_tolerance) {
            return Err(MomentumTriangleError::DoesNotClose {
                longest: canonical[0],
                shorter_sum,
                relative_tolerance: closure_tolerance,
            });
        }

        Ok(Self {
            perimeter: canonical.into_iter().sum(),
            second_over_first: canonical[1] / canonical[0],
            third_over_first: canonical[2] / canonical[0],
            canonical_to_input,
            input_to_canonical,
            closure_tolerance,
        })
    }

    /// Canonically ordered `[k_1, k_2, k_3]`, with `k_1 >= k_2 >= k_3`.
    pub fn canonical_wavenumbers(self) -> [f64; 3] {
        let first = self.perimeter / (1.0 + self.second_over_first + self.third_over_first);
        [
            first,
            first * self.second_over_first,
            first * self.third_over_first,
        ]
    }

    /// The perimeter `k_t = k_1 + k_2 + k_3`.
    pub const fn perimeter(self) -> f64 {
        self.perimeter
    }

    /// Shape coordinates `(k_2/k_1, k_3/k_1)`.
    pub const fn shape_ratios(self) -> (f64, f64) {
        (self.second_over_first, self.third_over_first)
    }

    /// Relative tolerance used to validate closure.
    pub const fn closure_tolerance(self) -> f64 {
        self.closure_tolerance
    }

    /// For each canonical leg, returns its position in the caller's input.
    pub const fn canonical_to_input(self) -> [usize; 3] {
        self.canonical_to_input
    }

    /// For each input leg, returns its position in canonical order.
    pub const fn input_to_canonical(self) -> [usize; 3] {
        self.input_to_canonical
    }

    /// Reorders leg-associated values from caller order to canonical order.
    pub fn canonicalize<T: Copy>(self, input: [T; 3]) -> [T; 3] {
        self.canonical_to_input.map(|index| input[index])
    }

    /// Restores leg-associated values from canonical order to caller order.
    pub fn restore_input_order<T: Copy>(self, canonical: [T; 3]) -> [T; 3] {
        self.input_to_canonical.map(|index| canonical[index])
    }
}

mod sealed {
    pub trait Sealed {}
}

/// Runtime phase-space slot tag used when iterating heterogeneous slots.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhaseSpaceSlotKind {
    /// A `Q^I` slot.
    Field,
    /// A `P^I` slot.
    Momentum,
}

/// Marker implemented by the two phase-space slot kinds.
pub trait PhaseSpaceSlot: sealed::Sealed + Copy + fmt::Debug + 'static {
    /// Stable machine-readable name.
    const NAME: &'static str;
    /// Runtime tag used when the compile-time marker is erased.
    const KIND: PhaseSpaceSlotKind;
}

/// The field-perturbation slot `Q^I`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FieldPerturbation;

/// The momentum-perturbation slot `P^I`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MomentumPerturbation;

impl sealed::Sealed for FieldPerturbation {}
impl sealed::Sealed for MomentumPerturbation {}

impl PhaseSpaceSlot for FieldPerturbation {
    const NAME: &'static str = "field";
    const KIND: PhaseSpaceSlotKind = PhaseSpaceSlotKind::Field;
}

impl PhaseSpaceSlot for MomentumPerturbation {
    const NAME: &'static str = "momentum";
    const KIND: PhaseSpaceSlotKind = PhaseSpaceSlotKind::Momentum;
}

/// An out-of-range field coordinate inside a phase-space slot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhaseSpaceIndexError {
    /// Requested field coordinate.
    pub field_index: usize,
    /// Number of available field coordinates.
    pub field_count: usize,
}

impl fmt::Display for PhaseSpaceIndexError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "phase-space field index {} is outside 0..{}",
            self.field_index, self.field_count
        )
    }
}

impl std::error::Error for PhaseSpaceIndexError {}

/// A field coordinate tagged with its phase-space slot at compile time.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhaseSpaceIndex<S, const N: usize> {
    field_index: usize,
    slot: PhantomData<S>,
}

impl<S: PhaseSpaceSlot, const N: usize> PhaseSpaceIndex<S, N> {
    /// Builds an index when `field_index < N`.
    pub const fn try_new(field_index: usize) -> Result<Self, PhaseSpaceIndexError> {
        if field_index < N {
            Ok(Self {
                field_index,
                slot: PhantomData,
            })
        } else {
            Err(PhaseSpaceIndexError {
                field_index,
                field_count: N,
            })
        }
    }

    /// Returns the field coordinate inside this slot.
    pub const fn field_index(self) -> usize {
        self.field_index
    }

    /// Erases the compile-time slot while retaining an explicit runtime tag.
    pub const fn erase(self) -> AnyPhaseSpaceIndex<N> {
        AnyPhaseSpaceIndex {
            slot: S::KIND,
            field_index: self.field_index,
        }
    }
}

/// A runtime-tagged phase-space index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AnyPhaseSpaceIndex<const N: usize> {
    /// Field or momentum slot.
    slot: PhaseSpaceSlotKind,
    /// Coordinate within the field space.
    field_index: usize,
}

impl<const N: usize> AnyPhaseSpaceIndex<N> {
    /// Builds a runtime-tagged index when `field_index < N`.
    pub const fn try_new(
        slot: PhaseSpaceSlotKind,
        field_index: usize,
    ) -> Result<Self, PhaseSpaceIndexError> {
        if field_index < N {
            Ok(Self { slot, field_index })
        } else {
            Err(PhaseSpaceIndexError {
                field_index,
                field_count: N,
            })
        }
    }

    /// Returns the runtime slot tag.
    pub const fn slot(self) -> PhaseSpaceSlotKind {
        self.slot
    }

    /// Returns the coordinate within field space.
    pub const fn field_index(self) -> usize {
        self.field_index
    }

    /// Flattened position in `[Q^0..Q^N, P^0..P^N]` order.
    pub const fn flat_index(self) -> usize {
        match self.slot {
            PhaseSpaceSlotKind::Field => self.field_index,
            PhaseSpaceSlotKind::Momentum => N + self.field_index,
        }
    }
}

/// A dimensioned two-point phase-space block with typed left/right slots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoPointBlock<L, R, D, V, const N: usize> {
    components: [[Quantity<D, V>; N]; N],
    slots: PhantomData<(L, R)>,
}

impl<L, R, D, V, const N: usize> TwoPointBlock<L, R, D, V, N> {
    /// Builds a block in left-slot/right-slot field-coordinate order.
    pub const fn new(components: [[Quantity<D, V>; N]; N]) -> Self {
        Self {
            components,
            slots: PhantomData,
        }
    }

    /// Borrows all block components.
    pub const fn components(&self) -> &[[Quantity<D, V>; N]; N] {
        &self.components
    }

    /// Transposes field coordinates and swaps the typed phase-space slots.
    pub fn transpose(self) -> TwoPointBlock<R, L, D, V, N>
    where
        Quantity<D, V>: Copy,
    {
        TwoPointBlock::new(core::array::from_fn(|left| {
            core::array::from_fn(|right| self.components[right][left])
        }))
    }
}

impl<L, R, D, const N: usize> TwoPointBlock<L, R, D, Complex64, N> {
    /// Conjugate-transposes field coordinates and swaps the typed slots.
    pub fn adjoint(self) -> TwoPointBlock<R, L, D, Complex64, N>
    where
        D: Copy,
    {
        TwoPointBlock::new(core::array::from_fn(|left| {
            core::array::from_fn(|right| Quantity::new(self.components[right][left].value().conj()))
        }))
    }
}

/// The complete equal-time covariance for `X^a = (Q^I, P^I)`.
///
/// The P3c convention is `Q^I = U^I` and `P^I = D_t U^I`, where the existing
/// linear solver integrates `v^I = a U^I`. The four blocks therefore have
/// dimensions `M^-1`, `M^0`, `M^0`, and `M^1`, respectively.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoPointPhaseSpaceCovariance<const N: usize> {
    /// `Sigma_QQ^{IJ} = <Q^I Q^{J*}>`.
    pub qq: TwoPointBlock<FieldPerturbation, FieldPerturbation, MMinus1, Complex64, N>,
    /// `Sigma_QP^{IJ} = <Q^I P^{J*}>`.
    pub qp: TwoPointBlock<FieldPerturbation, MomentumPerturbation, M0, Complex64, N>,
    /// `Sigma_PQ^{IJ} = <P^I Q^{J*}>`.
    pub pq: TwoPointBlock<MomentumPerturbation, FieldPerturbation, M0, Complex64, N>,
    /// `Sigma_PP^{IJ} = <P^I P^{J*}>`.
    pub pp: TwoPointBlock<MomentumPerturbation, MomentumPerturbation, M1, Complex64, N>,
}

impl<const N: usize> Add for TwoPointPhaseSpaceCovariance<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            qq: add_two_point_blocks(self.qq, rhs.qq),
            qp: add_two_point_blocks(self.qp, rhs.qp),
            pq: add_two_point_blocks(self.pq, rhs.pq),
            pp: add_two_point_blocks(self.pp, rhs.pp),
        }
    }
}

impl<const N: usize> Mul<f64> for TwoPointPhaseSpaceCovariance<N> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            qq: scale_two_point_block(self.qq, rhs),
            qp: scale_two_point_block(self.qp, rhs),
            pq: scale_two_point_block(self.pq, rhs),
            pp: scale_two_point_block(self.pp, rhs),
        }
    }
}

impl<const N: usize> Mul<Dimensionless> for TwoPointPhaseSpaceCovariance<N> {
    type Output = Self;

    fn mul(self, rhs: Dimensionless) -> Self::Output {
        self * rhs.into_value()
    }
}

fn add_two_point_blocks<L, R, D: Copy, const N: usize>(
    lhs: TwoPointBlock<L, R, D, Complex64, N>,
    rhs: TwoPointBlock<L, R, D, Complex64, N>,
) -> TwoPointBlock<L, R, D, Complex64, N> {
    TwoPointBlock::new(core::array::from_fn(|left| {
        core::array::from_fn(|right| lhs.components[left][right] + rhs.components[left][right])
    }))
}

fn scale_two_point_block<L, R, D: Copy, const N: usize>(
    block: TwoPointBlock<L, R, D, Complex64, N>,
    scale: f64,
) -> TwoPointBlock<L, R, D, Complex64, N> {
    TwoPointBlock::new(
        block
            .components
            .map(|row| row.map(|component| component * scale)),
    )
}

/// Supplies the immutable P2 background trajectory to transport systems.
pub trait TwoFieldBackgroundProvider {
    /// Interpolates the background at elapsed e-fold `N`.
    fn state_at_e_folds(&self, e_folds: f64) -> MultifieldState<2>;
}

/// Direct two-point phase-space transport on an immutable background.
///
/// The state is the covariance itself; no linear mode matrix is retained or
/// consulted after the initial covariance has been constructed.
pub struct TwoPointTransportSystem<'a, P, G, B> {
    background: &'a CovariantMultifieldSystem<P, G, 2>,
    trajectory: &'a B,
    k: Mass,
}

impl<'a, P, G, B> TwoPointTransportSystem<'a, P, G, B> {
    /// Builds a direct covariance system for one comoving wavenumber.
    pub const fn new(
        background: &'a CovariantMultifieldSystem<P, G, 2>,
        trajectory: &'a B,
        k: Mass,
    ) -> Self {
        Self {
            background,
            trajectory,
            k,
        }
    }
}

impl<P, G, B> System for TwoPointTransportSystem<'_, P, G, B>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
    B: TwoFieldBackgroundProvider,
{
    type IndependentVariable = Dimensionless;
    type State = TwoPointPhaseSpaceCovariance<2>;
    type Delta = TwoPointPhaseSpaceCovariance<2>;
    type Derivative = TwoPointPhaseSpaceCovariance<2>;

    fn derivative(&self, e_folds: Dimensionless, state: &Self::State) -> Self::Derivative {
        let background = self.trajectory.state_at_e_folds(e_folds.into_value());
        let hubble = self.background.hubble(&background).into_value();
        let connection = self.background.geometry.christoffel(&background.phi);
        let mass = self.background.perturbation_mass_matrix(&background);
        let physical_wavenumber = (self.k / background.a).into_value();
        let restoring: [[f64; 2]; 2] = core::array::from_fn(|upper| {
            core::array::from_fn(|lower| {
                mass.components()[upper][lower].value()
                    + f64::from(upper == lower) * physical_wavenumber.powi(2)
            })
        });
        let connection_flow: [[f64; 2]; 2] = core::array::from_fn(|upper| {
            core::array::from_fn(|lower| {
                (0..2)
                    .map(|velocity| {
                        connection.components()[upper][velocity][lower].value()
                            * background.pi[velocity].value()
                            / hubble
                    })
                    .sum()
            })
        });

        let qq = TwoPointBlock::new(core::array::from_fn(|left| {
            core::array::from_fn(|right| {
                let connection_terms: Complex64 = (0..2)
                    .map(|contracted| {
                        connection_flow[left][contracted]
                            * state.qq.components()[contracted][right].value()
                            + connection_flow[right][contracted]
                                * state.qq.components()[left][contracted].value()
                    })
                    .sum();
                Quantity::new(
                    (state.pq.components()[left][right].value()
                        + state.qp.components()[left][right].value())
                        / hubble
                        - connection_terms,
                )
            })
        }));
        let qp = TwoPointBlock::new(core::array::from_fn(|left| {
            core::array::from_fn(|right| {
                let (restoring_term, connection_term) = (0..2).fold(
                    (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)),
                    |(restoring_term, connection_term), index| {
                        (
                            restoring_term
                                + restoring[right][index]
                                    * state.qq.components()[left][index].value(),
                            connection_term
                                + connection_flow[left][index]
                                    * state.qp.components()[index][right].value()
                                + connection_flow[right][index]
                                    * state.qp.components()[left][index].value(),
                        )
                    },
                );
                Quantity::new(
                    state.pp.components()[left][right].value() / hubble
                        - state.qp.components()[left][right].value() * 3.0
                        - restoring_term / hubble
                        - connection_term,
                )
            })
        }));
        let pq = TwoPointBlock::new(core::array::from_fn(|left| {
            core::array::from_fn(|right| {
                let (restoring_term, connection_term) = (0..2).fold(
                    (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)),
                    |(restoring_term, connection_term), index| {
                        (
                            restoring_term
                                + restoring[left][index]
                                    * state.qq.components()[index][right].value(),
                            connection_term
                                + connection_flow[left][index]
                                    * state.pq.components()[index][right].value()
                                + connection_flow[right][index]
                                    * state.pq.components()[left][index].value(),
                        )
                    },
                );
                Quantity::new(
                    state.pp.components()[left][right].value() / hubble
                        - state.pq.components()[left][right].value() * 3.0
                        - restoring_term / hubble
                        - connection_term,
                )
            })
        }));
        let pp = TwoPointBlock::new(core::array::from_fn(|left| {
            core::array::from_fn(|right| {
                let (restoring_term, connection_term) = (0..2).fold(
                    (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)),
                    |(restoring_term, connection_term), index| {
                        (
                            restoring_term
                                + restoring[left][index]
                                    * state.qp.components()[index][right].value()
                                + restoring[right][index]
                                    * state.pq.components()[left][index].value(),
                            connection_term
                                + connection_flow[left][index]
                                    * state.pp.components()[index][right].value()
                                + connection_flow[right][index]
                                    * state.pp.components()[left][index].value(),
                        )
                    },
                );
                Quantity::new(
                    -state.pp.components()[left][right].value() * 6.0
                        - restoring_term / hubble
                        - connection_term,
                )
            })
        }));

        TwoPointPhaseSpaceCovariance { qq, qp, pq, pp }
    }
}

/// Builds every phase-space covariance block from the existing linear mode
/// matrix at one background state.
///
/// The supplied system computes `H` from the same state, preventing a mode
/// matrix from being paired with a Hubble value from a different epoch.
pub fn two_point_covariance_from_mode_state<P, G>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    state: &TwoFieldPerturbationState,
) -> TwoPointPhaseSpaceCovariance<2>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    let hubble = system.hubble(&state.background);
    phase_space_covariance(&state.modes, state.background.a.into_value(), hubble)
}

fn phase_space_covariance(
    modes: &crate::TwoFieldModes,
    scale_factor: f64,
    hubble: Hubble,
) -> TwoPointPhaseSpaceCovariance<2> {
    let q: [[ModeAmplitude; 2]; 2] = modes
        .rescaled_amplitude
        .map(|row| row.map(|entry| entry / scale_factor));
    let p: [[ModeVelocity; 2]; 2] = core::array::from_fn(|field| {
        core::array::from_fn(|mode| {
            ModeVelocity::new(
                (modes.rescaled_velocity[field][mode]
                    - modes.rescaled_amplitude[field][mode] * hubble)
                    .into_value()
                    / scale_factor,
            )
        })
    });

    let qq = TwoPointBlock::new(core::array::from_fn(|left| {
        core::array::from_fn(|right| {
            (0..2).fold(Quantity::new(Complex64::new(0.0, 0.0)), |sum, mode| {
                sum + q[left][mode] * ModeAmplitude::new(q[right][mode].value().conj())
            })
        })
    }));
    let qp = TwoPointBlock::new(core::array::from_fn(|left| {
        core::array::from_fn(|right| {
            (0..2).fold(Quantity::new(Complex64::new(0.0, 0.0)), |sum, mode| {
                sum + q[left][mode] * ModeVelocity::new(p[right][mode].value().conj())
            })
        })
    }));
    let pq = TwoPointBlock::new(core::array::from_fn(|left| {
        core::array::from_fn(|right| {
            (0..2).fold(Quantity::new(Complex64::new(0.0, 0.0)), |sum, mode| {
                sum + p[left][mode] * ModeAmplitude::new(q[right][mode].value().conj())
            })
        })
    }));
    let pp = TwoPointBlock::new(core::array::from_fn(|left| {
        core::array::from_fn(|right| {
            (0..2).fold(Quantity::new(Complex64::new(0.0, 0.0)), |sum, mode| {
                sum + p[left][mode] * ModeVelocity::new(p[right][mode].value().conj())
            })
        })
    }));

    TwoPointPhaseSpaceCovariance { qq, qp, pq, pp }
}

/// A dimensioned three-point block with one typed slot per Fourier leg.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThreePointBlock<A, B, C, D, V, const N: usize> {
    components: [[[Quantity<D, V>; N]; N]; N],
    slots: PhantomData<(A, B, C)>,
}

impl<A, B, C, D, V, const N: usize> ThreePointBlock<A, B, C, D, V, N> {
    /// Builds a block in first-/second-/third-leg field-coordinate order.
    pub const fn new(components: [[[Quantity<D, V>; N]; N]; N]) -> Self {
        Self {
            components,
            slots: PhantomData,
        }
    }

    /// Borrows all block components.
    pub const fn components(&self) -> &[[[Quantity<D, V>; N]; N]; N] {
        &self.components
    }

    /// Exchanges Fourier legs one and two, including their typed slots.
    pub fn swap_first_second(self) -> ThreePointBlock<B, A, C, D, V, N>
    where
        Quantity<D, V>: Copy,
    {
        ThreePointBlock::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| self.components[second][first][third])
            })
        }))
    }

    /// Exchanges Fourier legs two and three, including their typed slots.
    pub fn swap_second_third(self) -> ThreePointBlock<A, C, B, D, V, N>
    where
        Quantity<D, V>: Copy,
    {
        ThreePointBlock::new(core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| self.components[first][third][second])
            })
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Field, FieldVelocity, FlatFieldSpace, M0, Mass, MultifieldState,
        QuadraticMultifieldPotential, ScaleFactor, TwoFieldPerturbationSystem,
    };

    #[test]
    fn triangle_canonicalizes_scale_shape_and_leg_data() {
        let triangle = MomentumTriangle::try_new(2.0, 3.0, 4.0).unwrap();

        assert_eq!(triangle.canonical_wavenumbers(), [4.0, 3.0, 2.0]);
        assert_eq!(triangle.perimeter(), 9.0);
        assert_eq!(triangle.shape_ratios(), (0.75, 0.5));
        assert_eq!(triangle.canonical_to_input(), [2, 1, 0]);
        assert_eq!(triangle.input_to_canonical(), [2, 1, 0]);
        assert_eq!(triangle.canonicalize(['a', 'b', 'c']), ['c', 'b', 'a']);
        assert_eq!(
            triangle.restore_input_order(['c', 'b', 'a']),
            ['a', 'b', 'c']
        );
    }

    #[test]
    fn triangle_rejects_invalid_legs_tolerance_and_nonclosure() {
        assert!(matches!(
            MomentumTriangle::try_new(0.0, 1.0, 1.0),
            Err(MomentumTriangleError::InvalidWavenumber {
                leg: FourierLeg::One,
                ..
            })
        ));
        assert!(matches!(
            MomentumTriangle::try_new(f64::NAN, 1.0, 1.0),
            Err(MomentumTriangleError::InvalidWavenumber { .. })
        ));
        assert!(matches!(
            MomentumTriangle::try_new_with_tolerance(1.0, 1.0, 1.0, -1.0),
            Err(MomentumTriangleError::InvalidClosureTolerance { .. })
        ));
        assert!(matches!(
            MomentumTriangle::try_new(3.0, 1.0, 1.0),
            Err(MomentumTriangleError::DoesNotClose { .. })
        ));
    }

    #[test]
    fn folded_and_squeezed_boundaries_are_valid() {
        assert!(MomentumTriangle::try_new(2.0, 1.0, 1.0).is_ok());
        assert!(MomentumTriangle::try_new(1.0, 0.999, 0.001).is_ok());
        assert!(MomentumTriangle::try_new_with_tolerance(2.0 + 1.0e-14, 1.0, 1.0, 1.0e-14).is_ok());
    }

    #[test]
    fn phase_space_indices_keep_slot_and_field_coordinate_separate() {
        let field = PhaseSpaceIndex::<FieldPerturbation, 2>::try_new(1).unwrap();
        let momentum = PhaseSpaceIndex::<MomentumPerturbation, 2>::try_new(1).unwrap();

        assert_eq!(field.field_index(), 1);
        assert_eq!(field.erase().flat_index(), 1);
        assert_eq!(momentum.erase().flat_index(), 3);
        assert!(PhaseSpaceIndex::<FieldPerturbation, 2>::try_new(2).is_err());
    }

    #[test]
    fn two_point_transpose_swaps_slots_and_coordinates() {
        let block = TwoPointBlock::<FieldPerturbation, MomentumPerturbation, M0, f64, 2>::new([
            [Quantity::new(1.0), Quantity::new(2.0)],
            [Quantity::new(3.0), Quantity::new(4.0)],
        ]);
        let transposed: TwoPointBlock<MomentumPerturbation, FieldPerturbation, M0, f64, 2> =
            block.transpose();

        assert_eq!(transposed.components()[0][1].into_value(), 3.0);
        assert_eq!(transposed.components()[1][0].into_value(), 2.0);
    }

    #[test]
    fn complex_two_point_adjoint_conjugates_and_swaps_slots() {
        let block =
            TwoPointBlock::<FieldPerturbation, MomentumPerturbation, M0, Complex64, 2>::new([
                [
                    Quantity::new(Complex64::new(1.0, 2.0)),
                    Quantity::new(Complex64::new(3.0, 4.0)),
                ],
                [
                    Quantity::new(Complex64::new(5.0, 6.0)),
                    Quantity::new(Complex64::new(7.0, 8.0)),
                ],
            ]);
        let adjoint: TwoPointBlock<MomentumPerturbation, FieldPerturbation, M0, Complex64, 2> =
            block.adjoint();

        assert_eq!(
            adjoint.components()[0][1].into_value(),
            Complex64::new(5.0, -6.0)
        );
        assert_eq!(
            adjoint.components()[1][0].into_value(),
            Complex64::new(3.0, -4.0)
        );
    }

    #[test]
    fn mode_oracle_builds_all_bunch_davies_phase_space_blocks() {
        let background_system = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(2.0e-6)]),
            FlatFieldSpace,
        );
        let k = Mass::new(1.0e-2);
        let mode_system = TwoFieldPerturbationSystem::new(background_system, k);
        let background = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(-3.0e-6), FieldVelocity::new(-4.0e-6)],
            ScaleFactor::new(1.0),
        );
        let state = mode_system.bunch_davies_state(background, 10.0).unwrap();
        let covariance = two_point_covariance_from_mode_state(&mode_system.background, &state);
        let hubble = mode_system.background.hubble(&background).into_value();
        let frequency = (k / background.a).into_value();
        let qq_diagonal = 1.0 / (2.0 * k.into_value() * background.a.into_value().powi(2));

        assert_eq!(covariance.qq, covariance.qq.adjoint());
        assert_eq!(covariance.pp, covariance.pp.adjoint());
        assert_eq!(covariance.pq, covariance.qp.adjoint());

        for field in 0..2 {
            let qq = covariance.qq.components()[field][field].into_value();
            let qp = covariance.qp.components()[field][field].into_value();
            let pp = covariance.pp.components()[field][field].into_value();

            assert!((qq.re - qq_diagonal).abs() < 1.0e-12);
            assert!(qq.im.abs() < 1.0e-12);
            assert!((qp - Complex64::new(-hubble, frequency) * qq).norm() < 1.0e-12);
            assert!(
                (pp.re - (hubble * hubble + frequency * frequency) * qq_diagonal).abs() < 1.0e-12
            );
            assert!(pp.im.abs() < 1.0e-12);
        }
        assert!(covariance.qq.components()[0][1].value().norm() < 1.0e-12);
    }

    #[test]
    fn three_point_permutation_moves_slots_with_fourier_legs() {
        let components = core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    Quantity::new((100 * first + 10 * second + third) as f64)
                })
            })
        });
        let block = ThreePointBlock::<
            FieldPerturbation,
            MomentumPerturbation,
            FieldPerturbation,
            M0,
            f64,
            2,
        >::new(components);
        let swapped: ThreePointBlock<
            MomentumPerturbation,
            FieldPerturbation,
            FieldPerturbation,
            M0,
            f64,
            2,
        > = block.swap_first_second();

        assert_eq!(swapped.components()[1][0][1].into_value(), 11.0);
    }
}
