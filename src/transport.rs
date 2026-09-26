//! Typed kinematics and correlation containers for phase-space transport.
//!
//! This module separates Fourier legs, phase-space slots, field coordinates,
//! and physical dimensions. It provides the linear mode solver's two-point
//! oracle, direct two-point covariance dynamics, and typed three-point
//! transport infrastructure. Model-derived cubic coefficients and sourced
//! end-to-end evolution are not implemented yet.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Mul};

use num_complex::Complex64;

use crate::integrate::{AdaptiveSystem, ErrorTolerances, System};
use crate::{
    ContravariantField, CovariantMultifieldSystem, CurvedFieldSpaceGeometry, Dimensionless,
    FieldSpaceGeometry, Hubble, LowerFirst, M0, M1, M2, M3, MMinus1, MMinus2, MMinus3, Mass,
    ModeAmplitude, ModeVelocity, MultifieldPotential, MultifieldState, Quantity,
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

    /// Returns `k_i . k_j` in the canonical all-incoming convention.
    ///
    /// Momentum conservation is `k_1 + k_2 + k_3 = 0`. Therefore, for
    /// distinct legs, `k_i . k_j = (k_l^2 - k_i^2 - k_j^2) / 2`, where `l`
    /// is the remaining leg. The result has mass dimension two.
    pub fn wavevector_dot(self, left: FourierLeg, right: FourierLeg) -> Quantity<M2, f64> {
        let wavenumbers = self.canonical_wavenumbers();
        let left_index = left.index();
        let right_index = right.index();
        let value = if left_index == right_index {
            wavenumbers[left_index].powi(2)
        } else {
            let remaining_index = 3 - left_index - right_index;
            (wavenumbers[remaining_index].powi(2)
                - wavenumbers[left_index].powi(2)
                - wavenumbers[right_index].powi(2))
                / 2.0
        };
        Quantity::new(value)
    }

    /// Returns the cosine between two canonical all-incoming wavevectors.
    ///
    /// The result is clamped to `[-1, 1]` to remove closure-tolerance
    /// roundoff at folded boundaries.
    pub fn wavevector_cosine(self, left: FourierLeg, right: FourierLeg) -> f64 {
        let wavenumbers = self.canonical_wavenumbers();
        let denominator = wavenumbers[left.index()] * wavenumbers[right.index()];
        (self.wavevector_dot(left, right).into_value() / denominator).clamp(-1.0, 1.0)
    }

    /// Returns `1 - (k_i . k_j)^2 / (k_i^2 k_j^2)`.
    ///
    /// This dimensionless transverse factor occurs repeatedly in the cubic
    /// `A`, `B`, and `C` coefficients.
    pub fn transverse_factor(self, left: FourierLeg, right: FourierLeg) -> f64 {
        let cosine = self.wavevector_cosine(left, right);
        (1.0 - cosine * cosine).max(0.0)
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
    /// The second phase-space slot, either velocity or canonical momentum.
    ///
    /// The semantic distinction is retained by its compile-time marker and is
    /// intentionally erased by this structural tag.
    Momentum,
}

/// Marker implemented by a semantic phase-space slot.
///
/// `KIND` describes the structural position used by the two-slot contraction
/// algebra. Distinct momentum semantics remain distinct through the marker
/// type even when they occupy the same structural position.
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
///
/// In the validated P3c convention this means the covariant velocity
/// `P^I = D_t Q^I`; it is not the nonlinear canonical momentum.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MomentumPerturbation;

/// The rescaled canonical-momentum slot `Pi^I = a^-3 p^I`.
///
/// It agrees with [`MomentumPerturbation`] at linear order but differs at
/// quadratic order in perturbations. Keeping a distinct marker prevents a
/// canonical cubic kernel from being contracted with velocity correlators.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CanonicalMomentumPerturbation;

impl sealed::Sealed for FieldPerturbation {}
impl sealed::Sealed for MomentumPerturbation {}
impl sealed::Sealed for CanonicalMomentumPerturbation {}

impl PhaseSpaceSlot for FieldPerturbation {
    const NAME: &'static str = "field";
    const KIND: PhaseSpaceSlotKind = PhaseSpaceSlotKind::Field;
}

impl PhaseSpaceSlot for MomentumPerturbation {
    const NAME: &'static str = "momentum";
    const KIND: PhaseSpaceSlotKind = PhaseSpaceSlotKind::Momentum;
}

impl PhaseSpaceSlot for CanonicalMomentumPerturbation {
    const NAME: &'static str = "canonical_momentum";
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
/// dimensions `M^-1`, `M^0`, `M^0`, and `M^1`, respectively. The trailing
/// marker can instead be [`CanonicalMomentumPerturbation`] for a canonical
/// phase space; the default preserves the P3c velocity convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoPointPhaseSpaceCovariance<const N: usize, P: PhaseSpaceSlot = MomentumPerturbation> {
    /// `Sigma_QQ^{IJ} = <Q^I Q^{J*}>`.
    pub qq: TwoPointBlock<FieldPerturbation, FieldPerturbation, MMinus1, Complex64, N>,
    /// `Sigma_QP^{IJ} = <Q^I P^{J*}>`.
    pub qp: TwoPointBlock<FieldPerturbation, P, M0, Complex64, N>,
    /// `Sigma_PQ^{IJ} = <P^I Q^{J*}>`.
    pub pq: TwoPointBlock<P, FieldPerturbation, M0, Complex64, N>,
    /// `Sigma_PP^{IJ} = <P^I P^{J*}>`.
    pub pp: TwoPointBlock<P, P, M1, Complex64, N>,
}

/// A covariance whose second phase-space slot is rescaled canonical momentum.
pub type CanonicalTwoPointPhaseSpaceCovariance<const N: usize> =
    TwoPointPhaseSpaceCovariance<N, CanonicalMomentumPerturbation>;

impl<const N: usize, P: PhaseSpaceSlot> Add for TwoPointPhaseSpaceCovariance<N, P> {
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

impl<const N: usize, P: PhaseSpaceSlot> TwoPointPhaseSpaceCovariance<N, P> {
    fn component(&self, slots: [PhaseSpaceSlotKind; 2], indices: [usize; 2]) -> Complex64 {
        let [left, right] = indices;
        match slots {
            [PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Field] => {
                self.qq.components[left][right].into_value()
            }
            [PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Momentum] => {
                self.qp.components[left][right].into_value()
            }
            [PhaseSpaceSlotKind::Momentum, PhaseSpaceSlotKind::Field] => {
                self.pq.components[left][right].into_value()
            }
            [PhaseSpaceSlotKind::Momentum, PhaseSpaceSlotKind::Momentum] => {
                self.pp.components[left][right].into_value()
            }
        }
    }
}

impl<const N: usize, P: PhaseSpaceSlot> Mul<f64> for TwoPointPhaseSpaceCovariance<N, P> {
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

impl<const N: usize, P: PhaseSpaceSlot> Mul<Dimensionless> for TwoPointPhaseSpaceCovariance<N, P> {
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

/// Absolute error floors for direct two-point covariance evolution.
///
/// Each block retains its own physical dimension. Zero floors select pure
/// relative control, which is appropriate while the complete block norm is
/// nonzero.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TwoPointTransportTolerance {
    /// Absolute tolerance for the `QQ` block.
    pub qq: Quantity<MMinus1, f64>,
    /// Absolute tolerance for the `QP` block.
    pub qp: Quantity<M0, f64>,
    /// Absolute tolerance for the `PQ` block.
    pub pq: Quantity<M0, f64>,
    /// Absolute tolerance for the `PP` block.
    pub pp: Quantity<M1, f64>,
}

impl TwoPointTransportTolerance {
    /// Uses pure relative control for every block.
    pub const fn pure_relative() -> Self {
        Self {
            qq: Quantity::new(0.0),
            qp: Quantity::new(0.0),
            pq: Quantity::new(0.0),
            pp: Quantity::new(0.0),
        }
    }
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

impl<P, G, B> AdaptiveSystem for TwoPointTransportSystem<'_, P, G, B>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
    B: TwoFieldBackgroundProvider,
{
    type AbsoluteTolerance = TwoPointTransportTolerance;

    fn error_ratio(
        &self,
        current: &Self::State,
        next: &Self::State,
        error: &Self::Delta,
        tolerances: &ErrorTolerances<Self::AbsoluteTolerance>,
    ) -> f64 {
        covariance_block_error_ratio(
            &current.qq,
            &next.qq,
            &error.qq,
            tolerances.absolute.qq,
            tolerances.relative,
        )
        .max(covariance_block_error_ratio(
            &current.qp,
            &next.qp,
            &error.qp,
            tolerances.absolute.qp,
            tolerances.relative,
        ))
        .max(covariance_block_error_ratio(
            &current.pq,
            &next.pq,
            &error.pq,
            tolerances.absolute.pq,
            tolerances.relative,
        ))
        .max(covariance_block_error_ratio(
            &current.pp,
            &next.pp,
            &error.pp,
            tolerances.absolute.pp,
            tolerances.relative,
        ))
    }
}

fn covariance_block_error_ratio<L, R, D: Copy>(
    current: &TwoPointBlock<L, R, D, Complex64, 2>,
    next: &TwoPointBlock<L, R, D, Complex64, 2>,
    error: &TwoPointBlock<L, R, D, Complex64, 2>,
    absolute: Quantity<D, f64>,
    relative: f64,
) -> f64 {
    let reference = current
        .components()
        .iter()
        .chain(next.components().iter())
        .flatten()
        .map(|component| component.value().norm())
        .fold(0.0, f64::max);
    let scale = absolute.into_value() + relative * reference;

    if scale == 0.0 {
        return if error
            .components()
            .iter()
            .flatten()
            .all(|component| component.value().norm() == 0.0)
        {
            0.0
        } else {
            f64::INFINITY
        };
    }

    error
        .components()
        .iter()
        .flatten()
        .map(|component| component.value().norm() / scale)
        .fold(0.0, f64::max)
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

impl<A, B, C, D: Copy, const N: usize> ThreePointBlock<A, B, C, D, Complex64, N> {
    /// Complex-conjugates every component without moving its Fourier leg.
    pub fn conjugate(self) -> Self {
        Self::new(self.components.map(|matrix| {
            matrix.map(|row| row.map(|component| Quantity::new(component.value().conj())))
        }))
    }
}

/// The complete connected three-point function for `X^a = (Q^I, P^I)`.
///
/// After the single momentum-conserving delta function has been removed, the
/// blocks with zero through three momentum slots have dimensions `M^-3`,
/// `M^-2`, `M^-1`, and `M^0`. These are not obtained by simply cubing the
/// `M^-1/2` linear mode-function dimension, because a stripped three-point
/// function retains a different Fourier-volume dimension than a stripped
/// two-point function. Fourier-leg permutations move both the component axis
/// and its phase-space slot. The trailing marker fixes whether the second slot
/// is velocity or rescaled canonical momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThreePointPhaseSpaceCorrelation<const N: usize, P: PhaseSpaceSlot = MomentumPerturbation>
{
    /// `Alpha_QQQ^{IJK}`.
    pub qqq: ThreePointBlock<
        FieldPerturbation,
        FieldPerturbation,
        FieldPerturbation,
        MMinus3,
        Complex64,
        N,
    >,
    /// `Alpha_QQP^{IJK}`.
    pub qqp: ThreePointBlock<FieldPerturbation, FieldPerturbation, P, MMinus2, Complex64, N>,
    /// `Alpha_QPQ^{IJK}`.
    pub qpq: ThreePointBlock<FieldPerturbation, P, FieldPerturbation, MMinus2, Complex64, N>,
    /// `Alpha_PQQ^{IJK}`.
    pub pqq: ThreePointBlock<P, FieldPerturbation, FieldPerturbation, MMinus2, Complex64, N>,
    /// `Alpha_QPP^{IJK}`.
    pub qpp: ThreePointBlock<FieldPerturbation, P, P, MMinus1, Complex64, N>,
    /// `Alpha_PQP^{IJK}`.
    pub pqp: ThreePointBlock<P, FieldPerturbation, P, MMinus1, Complex64, N>,
    /// `Alpha_PPQ^{IJK}`.
    pub ppq: ThreePointBlock<P, P, FieldPerturbation, MMinus1, Complex64, N>,
    /// `Alpha_PPP^{IJK}`.
    pub ppp: ThreePointBlock<P, P, P, M0, Complex64, N>,
}

/// A connected three-point function in canonical phase-space variables.
pub type CanonicalThreePointPhaseSpaceCorrelation<const N: usize> =
    ThreePointPhaseSpaceCorrelation<N, CanonicalMomentumPerturbation>;

/// The quadratic phase-space flow `u^a_bc` for one ordered momentum vertex.
///
/// The first slot is the output and the last two are ordered inputs. Every
/// block dimension is `dim(output) - dim(input one) - dim(input two)`, so its
/// contraction with two covariance blocks has the dimension of the target
/// three-point block.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticPhaseSpaceKernel<const N: usize, P: PhaseSpaceSlot = MomentumPerturbation> {
    /// `u^Q_QQ`, dimension `M^-1`.
    pub q_qq:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
    /// `u^Q_QP`, dimension `M^-2`.
    pub q_qp: ThreePointBlock<FieldPerturbation, FieldPerturbation, P, MMinus2, f64, N>,
    /// `u^Q_PQ`, dimension `M^-2`.
    pub q_pq: ThreePointBlock<FieldPerturbation, P, FieldPerturbation, MMinus2, f64, N>,
    /// `u^Q_PP`, dimension `M^-3`.
    pub q_pp: ThreePointBlock<FieldPerturbation, P, P, MMinus3, f64, N>,
    /// `u^P_QQ`, dimension `M^0`.
    pub p_qq: ThreePointBlock<P, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `u^P_QP`, dimension `M^-1`.
    pub p_qp: ThreePointBlock<P, FieldPerturbation, P, MMinus1, f64, N>,
    /// `u^P_PQ`, dimension `M^-1`.
    pub p_pq: ThreePointBlock<P, P, FieldPerturbation, MMinus1, f64, N>,
    /// `u^P_PP`, dimension `M^-2`.
    pub p_pp: ThreePointBlock<P, P, P, MMinus2, f64, N>,
}

/// A quadratic flow whose second phase-space slot is canonical momentum.
pub type CanonicalQuadraticPhaseSpaceKernel<const N: usize> =
    QuadraticPhaseSpaceKernel<N, CanonicalMomentumPerturbation>;

/// Raised cubic-Hamiltonian coefficients in output/first-input/second-input order.
///
/// All component arrays use `[I][J][K]`, even when `I` originated as the
/// third covariant index. This makes the subsequent phase-space block mapping
/// explicit and prevents first- and third-index raising from being confused.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalCubicHamiltonianCoefficients<const N: usize> {
    /// `A^I_JK`, stored as `[I][J][K]`, dimension `M^1`.
    pub a_first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M1, f64, N>,
    /// `B^I_JK`, with the first original index raised, dimension `M^0`.
    pub b_first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `B_JK^I`, with the third original index raised, stored as `[I][J][K]`.
    pub b_third_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `C^I_JK`, with the first original index raised, dimension `M^-1`.
    pub c_first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
    /// `C_JK^I`, with the third original index raised, stored as `[I][J][K]`.
    pub c_third_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
}

/// The two raised forms of the symmetrized cubic coefficient `B_IJK`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalCubicBCoefficients<const N: usize> {
    /// `B^I_JK`, stored in output/first-input/second-input order.
    pub first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `B_JK^I`, also stored in output/first-input/second-input order.
    pub third_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
}

/// The two raised forms of the symmetrized cubic coefficient `C_IJK`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalCubicCCoefficients<const N: usize> {
    /// `C^I_JK`, stored in output/first-input/second-input order.
    pub first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
    /// `C_JK^I`, also stored in output/first-input/second-input order.
    pub third_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
}

/// Invalid background data supplied while evaluating a cubic coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CubicCoefficientError {
    /// The reduced Planck mass was not positive and finite.
    InvalidPlanckMass(f64),
    /// The expanding-branch Hubble rate was not positive and finite.
    InvalidHubble(f64),
}

impl fmt::Display for CubicCoefficientError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPlanckMass(value) => write!(
                formatter,
                "cubic coefficients require a positive finite reduced Planck mass, got {value}"
            ),
            Self::InvalidHubble(value) => write!(
                formatter,
                "cubic coefficients require a positive finite Hubble rate, got {value}"
            ),
        }
    }
}

impl std::error::Error for CubicCoefficientError {}

/// Evaluates the symmetrized curved-field-space coefficient `B_IJK` and
/// returns the two raised forms required by the canonical phase-space kernel.
///
/// The background equation of motion is used covariantly to form
/// `Z_I = D_t dot(phi)_I + dot(phi)_I dot(phi)^2 / (2 H M_pl^2)`.
/// Field indices and their attached canonical Fourier legs are symmetrized
/// together.
pub fn canonical_cubic_b_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicBCoefficients<N>, CubicCoefficientError>
where
    P: MultifieldPotential<N>,
    G: CurvedFieldSpaceGeometry<N>,
{
    let planck_value = system.m_pl.into_value();
    if !planck_value.is_finite() || planck_value <= 0.0 {
        return Err(CubicCoefficientError::InvalidPlanckMass(planck_value));
    }
    let hubble = system.hubble(state);
    let hubble_value = hubble.into_value();
    if !hubble_value.is_finite() || hubble_value <= 0.0 {
        return Err(CubicCoefficientError::InvalidHubble(hubble_value));
    }

    let metric = system.geometry.metric(&state.phi);
    let inverse = system.geometry.inverse_metric(&state.phi);
    let lowered_riemann = metric.lower_first(&system.geometry.riemann(&state.phi));
    let velocity = ContravariantField::new(state.pi);
    let lowered_velocity = metric.lower(&velocity);
    let velocity_squared = lowered_velocity.contract(&velocity);
    let gradient = system.potential.gradient(state.phi);
    let planck = Mass::new(planck_value);
    let planck_squared = planck * planck;
    let z: [Quantity<M3, f64>; N] = core::array::from_fn(|i| {
        -(lowered_velocity.components()[i] * hubble * 3.0) - gradient[i]
            + lowered_velocity.components()[i] * velocity_squared / (hubble * planck_squared * 2.0)
    });
    let legs = [FourierLeg::One, FourierLeg::Two, FourierLeg::Three];
    let wavenumbers = triangle.canonical_wavenumbers();

    let unsymmetrized = |i: usize,
                         j: usize,
                         k: usize,
                         first_leg: FourierLeg,
                         second_leg: FourierLeg,
                         third_leg: FourierLeg| {
        let velocity_i = lowered_velocity.components()[i];
        let velocity_j = lowered_velocity.components()[j];
        let velocity_k = lowered_velocity.components()[k];
        let curvature = (0..N)
            .map(|l| {
                (lowered_riemann.components()[k][i][j][l]
                    + lowered_riemann.components()[k][j][i][l])
                    * velocity.components()[l]
                    * (2.0 / 3.0)
            })
            .fold(Quantity::new(0.0), Add::add);
        let first_wavenumber = wavenumbers[first_leg.index()];
        let first_second_ratio =
            triangle.wavevector_dot(first_leg, second_leg).into_value() / first_wavenumber.powi(2);

        curvature
            - velocity_i * z[j] * velocity_k
                / (hubble * hubble * hubble * planck_squared * planck_squared)
                * (triangle.transverse_factor(second_leg, third_leg) / 4.0)
            + velocity_i * velocity_j * velocity_k
                / (hubble * hubble * planck_squared * planck_squared)
                * 0.25
            - z[i] / (hubble * planck_squared)
                * metric.components()[j][k].into_value()
                * first_second_ratio
    };

    let lower: [[[Quantity<M0, f64>; N]; N]; N] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (unsymmetrized(i, j, k, legs[0], legs[1], legs[2])
                    + unsymmetrized(j, i, k, legs[1], legs[0], legs[2]))
                    * 0.5
            })
        })
    });

    let first_raised = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        lower[contracted][j][k] * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));
    let third_raised = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        lower[j][k][contracted] * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));

    Ok(CanonicalCubicBCoefficients {
        first_raised,
        third_raised,
    })
}

/// Evaluates the symmetrized curved-field-space coefficient `C_IJK` and
/// returns the two raised forms required by the canonical phase-space kernel.
///
/// Field indices and their attached canonical Fourier legs are symmetrized
/// together. All arithmetic that carries physical dimension remains inside
/// [`Quantity`]; only metric components and momentum ratios are raw
/// dimensionless factors.
pub fn canonical_cubic_c_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicCCoefficients<N>, CubicCoefficientError>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    let planck_value = system.m_pl.into_value();
    if !planck_value.is_finite() || planck_value <= 0.0 {
        return Err(CubicCoefficientError::InvalidPlanckMass(planck_value));
    }
    let hubble = system.hubble(state);
    let hubble_value = hubble.into_value();
    if !hubble_value.is_finite() || hubble_value <= 0.0 {
        return Err(CubicCoefficientError::InvalidHubble(hubble_value));
    }

    let metric = system.geometry.metric(&state.phi);
    let inverse = system.geometry.inverse_metric(&state.phi);
    let velocity = ContravariantField::new(state.pi);
    let lowered_velocity = metric.lower(&velocity);
    let planck = Mass::new(planck_value);
    let legs = [FourierLeg::One, FourierLeg::Two, FourierLeg::Three];
    let wavenumbers = triangle.canonical_wavenumbers();

    let unsymmetrized = |i: usize,
                         j: usize,
                         k: usize,
                         first_leg: FourierLeg,
                         second_leg: FourierLeg,
                         third_leg: FourierLeg| {
        let velocity_i = lowered_velocity.components()[i];
        let velocity_j = lowered_velocity.components()[j];
        let velocity_k = lowered_velocity.components()[k];
        let planck_squared = planck * planck;
        let first_wavenumber = wavenumbers[first_leg.index()];
        let first_third_ratio =
            triangle.wavevector_dot(first_leg, third_leg).into_value() / first_wavenumber.powi(2);

        -(velocity_k / (hubble * planck_squared) * metric.components()[i][j].into_value() * 0.5)
            + velocity_i * velocity_j * velocity_k
                / (hubble * hubble * hubble * planck_squared * planck_squared)
                * (triangle.transverse_factor(first_leg, second_leg) / 8.0)
            + velocity_i / (hubble * planck_squared)
                * metric.components()[j][k].into_value()
                * first_third_ratio
    };

    let lower: [[[Quantity<MMinus1, f64>; N]; N]; N] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (unsymmetrized(i, j, k, legs[0], legs[1], legs[2])
                    + unsymmetrized(j, i, k, legs[1], legs[0], legs[2]))
                    * 0.5
            })
        })
    });

    let first_raised = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        lower[contracted][j][k] * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));
    let third_raised = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        lower[j][k][contracted] * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));

    Ok(CanonicalCubicCCoefficients {
        first_raised,
        third_raised,
    })
}

/// An invalid Hubble rate supplied while converting a cosmic-time cubic flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalKernelError {
    /// Rejected Hubble value in reduced Planck units.
    pub hubble: f64,
}

impl fmt::Display for CanonicalKernelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "canonical cubic kernel requires a positive finite Hubble rate, got {}",
            self.hubble
        )
    }
}

impl std::error::Error for CanonicalKernelError {}

impl<const N: usize> QuadraticPhaseSpaceKernel<N, CanonicalMomentumPerturbation> {
    /// Maps cosmic-time `A`, `B`, and `C` coefficients to the eight e-fold
    /// canonical phase-space blocks.
    ///
    /// The conversion divides every block by `H` exactly once and applies the
    /// ordered input-index exchanges required by the Hamilton equations.
    pub fn from_cubic_hamiltonian(
        coefficients: CanonicalCubicHamiltonianCoefficients<N>,
        hubble: Hubble,
    ) -> Result<Self, CanonicalKernelError> {
        let hubble_value = hubble.into_value();
        if !hubble_value.is_finite() || hubble_value <= 0.0 {
            return Err(CanonicalKernelError {
                hubble: hubble_value,
            });
        }
        let hubble = Hubble::new(hubble_value);

        Ok(Self {
            q_qq: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        -(coefficients.b_third_raised.components[i][j][k] / hubble)
                    })
                })
            })),
            q_qp: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        -(coefficients.c_first_raised.components[i][k][j] / hubble)
                    })
                })
            })),
            q_pq: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        -(coefficients.c_first_raised.components[i][j][k] / hubble)
                    })
                })
            })),
            q_pp: zero_real_three_point_block(),
            p_qq: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        coefficients.a_first_raised.components[i][j][k] * 3.0 / hubble
                    })
                })
            })),
            p_qp: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        coefficients.b_first_raised.components[i][j][k] / hubble
                    })
                })
            })),
            p_pq: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        coefficients.b_first_raised.components[i][k][j] / hubble
                    })
                })
            })),
            p_pp: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        coefficients.c_third_raised.components[i][k][j] / hubble
                    })
                })
            })),
        })
    }
}

impl<const N: usize, P: PhaseSpaceSlot> QuadraticPhaseSpaceKernel<N, P> {
    /// Returns a kernel whose eight blocks vanish.
    pub fn zero() -> Self {
        Self {
            q_qq: zero_real_three_point_block(),
            q_qp: zero_real_three_point_block(),
            q_pq: zero_real_three_point_block(),
            q_pp: zero_real_three_point_block(),
            p_qq: zero_real_three_point_block(),
            p_qp: zero_real_three_point_block(),
            p_pq: zero_real_three_point_block(),
            p_pp: zero_real_three_point_block(),
        }
    }

    fn component(&self, slots: [PhaseSpaceSlotKind; 3], indices: [usize; 3]) -> f64 {
        let [output, first, second] = indices;
        match slots {
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
            ] => self.q_qq.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ] => self.q_qp.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ] => self.q_pq.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
            ] => self.q_pp.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
            ] => self.p_qq.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ] => self.p_qp.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ] => self.p_pq.components[output][first][second].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
            ] => self.p_pp.components[output][first][second].into_value(),
        }
    }
}

/// Contracts three ordered quadratic kernels with the three leg covariances.
///
/// `kernels[a]` has output leg `a`. Its two input slots follow the remaining
/// legs in ascending order: `[1, 2]`, `[0, 2]`, and `[0, 1]`. The returned
/// value is exactly the three cyclic `u Sigma Sigma` terms in `D_N Alpha`; it
/// does not include homogeneous propagation or an initial condition.
pub fn three_point_source_from_kernels<const N: usize, P: PhaseSpaceSlot>(
    kernels: [&QuadraticPhaseSpaceKernel<N, P>; 3],
    covariances: [&TwoPointPhaseSpaceCovariance<N, P>; 3],
) -> ThreePointPhaseSpaceCorrelation<N, P> {
    let source_block = |output_slots: [PhaseSpaceSlotKind; 3]| {
        core::array::from_fn(|first| {
            core::array::from_fn(|second| {
                core::array::from_fn(|third| {
                    let output_indices = [first, second, third];
                    (0..3)
                        .map(|output_leg| {
                            let input_legs = match output_leg {
                                0 => [1, 2],
                                1 => [0, 2],
                                2 => [0, 1],
                                _ => unreachable!(),
                            };
                            [PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Momentum]
                                .into_iter()
                                .flat_map(|first_input_slot| {
                                    [PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Momentum]
                                        .into_iter()
                                        .flat_map(move |second_input_slot| {
                                            (0..N).flat_map(move |first_input_index| {
                                                (0..N).map(move |second_input_index| {
                                                    (
                                                        first_input_slot,
                                                        second_input_slot,
                                                        first_input_index,
                                                        second_input_index,
                                                    )
                                                })
                                            })
                                        })
                                })
                                .map(
                                    |(
                                        first_input_slot,
                                        second_input_slot,
                                        first_input_index,
                                        second_input_index,
                                    )| {
                                        kernels[output_leg].component(
                                            [
                                                output_slots[output_leg],
                                                first_input_slot,
                                                second_input_slot,
                                            ],
                                            [
                                                output_indices[output_leg],
                                                first_input_index,
                                                second_input_index,
                                            ],
                                        ) * covariances[input_legs[0]].component(
                                            [first_input_slot, output_slots[input_legs[0]]],
                                            [first_input_index, output_indices[input_legs[0]]],
                                        ) * covariances[input_legs[1]].component(
                                            [second_input_slot, output_slots[input_legs[1]]],
                                            [second_input_index, output_indices[input_legs[1]]],
                                        )
                                    },
                                )
                                .sum::<Complex64>()
                        })
                        .sum()
                })
            })
        })
    };

    ThreePointPhaseSpaceCorrelation {
        qqq: complex_three_point_block(source_block([PhaseSpaceSlotKind::Field; 3])),
        qqp: complex_three_point_block(source_block([
            PhaseSpaceSlotKind::Field,
            PhaseSpaceSlotKind::Field,
            PhaseSpaceSlotKind::Momentum,
        ])),
        qpq: complex_three_point_block(source_block([
            PhaseSpaceSlotKind::Field,
            PhaseSpaceSlotKind::Momentum,
            PhaseSpaceSlotKind::Field,
        ])),
        pqq: complex_three_point_block(source_block([
            PhaseSpaceSlotKind::Momentum,
            PhaseSpaceSlotKind::Field,
            PhaseSpaceSlotKind::Field,
        ])),
        qpp: complex_three_point_block(source_block([
            PhaseSpaceSlotKind::Field,
            PhaseSpaceSlotKind::Momentum,
            PhaseSpaceSlotKind::Momentum,
        ])),
        pqp: complex_three_point_block(source_block([
            PhaseSpaceSlotKind::Momentum,
            PhaseSpaceSlotKind::Field,
            PhaseSpaceSlotKind::Momentum,
        ])),
        ppq: complex_three_point_block(source_block([
            PhaseSpaceSlotKind::Momentum,
            PhaseSpaceSlotKind::Momentum,
            PhaseSpaceSlotKind::Field,
        ])),
        ppp: complex_three_point_block(source_block([PhaseSpaceSlotKind::Momentum; 3])),
    }
}

fn complex_three_point_block<A, B, C, D, const N: usize>(
    components: [[[Complex64; N]; N]; N],
) -> ThreePointBlock<A, B, C, D, Complex64, N> {
    ThreePointBlock::new(components.map(|matrix| matrix.map(|row| row.map(Quantity::new))))
}

impl<const N: usize, P: PhaseSpaceSlot> ThreePointPhaseSpaceCorrelation<N, P> {
    /// Returns the additive identity in every phase-space block.
    pub fn zero() -> Self {
        Self {
            qqq: zero_three_point_block(),
            qqp: zero_three_point_block(),
            qpq: zero_three_point_block(),
            pqq: zero_three_point_block(),
            qpp: zero_three_point_block(),
            pqp: zero_three_point_block(),
            ppq: zero_three_point_block(),
            ppp: zero_three_point_block(),
        }
    }

    /// Exchanges Fourier legs one and two, including their phase-space slots.
    pub fn swap_first_second(self) -> Self {
        Self {
            qqq: self.qqq.swap_first_second(),
            qqp: self.qqp.swap_first_second(),
            qpq: self.pqq.swap_first_second(),
            pqq: self.qpq.swap_first_second(),
            qpp: self.pqp.swap_first_second(),
            pqp: self.qpp.swap_first_second(),
            ppq: self.ppq.swap_first_second(),
            ppp: self.ppp.swap_first_second(),
        }
    }

    /// Exchanges Fourier legs two and three, including their phase-space slots.
    pub fn swap_second_third(self) -> Self {
        Self {
            qqq: self.qqq.swap_second_third(),
            qqp: self.qpq.swap_second_third(),
            qpq: self.qqp.swap_second_third(),
            pqq: self.pqq.swap_second_third(),
            qpp: self.qpp.swap_second_third(),
            pqp: self.ppq.swap_second_third(),
            ppq: self.pqp.swap_second_third(),
            ppp: self.ppp.swap_second_third(),
        }
    }

    /// Complex-conjugates all eight phase-space blocks.
    pub fn conjugate(self) -> Self {
        Self {
            qqq: self.qqq.conjugate(),
            qqp: self.qqp.conjugate(),
            qpq: self.qpq.conjugate(),
            pqq: self.pqq.conjugate(),
            qpp: self.qpp.conjugate(),
            pqp: self.pqp.conjugate(),
            ppq: self.ppq.conjugate(),
            ppp: self.ppp.conjugate(),
        }
    }

    fn component(&self, slots: [PhaseSpaceSlotKind; 3], indices: [usize; 3]) -> Complex64 {
        let [first, second, third] = indices;
        match slots {
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
            ] => self.qqq.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ] => self.qqp.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ] => self.qpq.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
            ] => self.pqq.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
            ] => self.qpp.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ] => self.pqp.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ] => self.ppq.components[first][second][third].into_value(),
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
            ] => self.ppp.components[first][second][third].into_value(),
        }
    }
}

impl<const N: usize, P: PhaseSpaceSlot> Add for ThreePointPhaseSpaceCorrelation<N, P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            qqq: add_three_point_blocks(self.qqq, rhs.qqq),
            qqp: add_three_point_blocks(self.qqp, rhs.qqp),
            qpq: add_three_point_blocks(self.qpq, rhs.qpq),
            pqq: add_three_point_blocks(self.pqq, rhs.pqq),
            qpp: add_three_point_blocks(self.qpp, rhs.qpp),
            pqp: add_three_point_blocks(self.pqp, rhs.pqp),
            ppq: add_three_point_blocks(self.ppq, rhs.ppq),
            ppp: add_three_point_blocks(self.ppp, rhs.ppp),
        }
    }
}

impl<const N: usize, P: PhaseSpaceSlot> Mul<f64> for ThreePointPhaseSpaceCorrelation<N, P> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            qqq: scale_three_point_block(self.qqq, rhs),
            qqp: scale_three_point_block(self.qqp, rhs),
            qpq: scale_three_point_block(self.qpq, rhs),
            pqq: scale_three_point_block(self.pqq, rhs),
            qpp: scale_three_point_block(self.qpp, rhs),
            pqp: scale_three_point_block(self.pqp, rhs),
            ppq: scale_three_point_block(self.ppq, rhs),
            ppp: scale_three_point_block(self.ppp, rhs),
        }
    }
}

impl<const N: usize, P: PhaseSpaceSlot> Mul<Dimensionless>
    for ThreePointPhaseSpaceCorrelation<N, P>
{
    type Output = Self;

    fn mul(self, rhs: Dimensionless) -> Self::Output {
        self * rhs.into_value()
    }
}

fn add_three_point_blocks<A, B, C, D: Copy, const N: usize>(
    lhs: ThreePointBlock<A, B, C, D, Complex64, N>,
    rhs: ThreePointBlock<A, B, C, D, Complex64, N>,
) -> ThreePointBlock<A, B, C, D, Complex64, N> {
    ThreePointBlock::new(core::array::from_fn(|first| {
        core::array::from_fn(|second| {
            core::array::from_fn(|third| {
                lhs.components[first][second][third] + rhs.components[first][second][third]
            })
        })
    }))
}

fn scale_three_point_block<A, B, C, D: Copy, const N: usize>(
    block: ThreePointBlock<A, B, C, D, Complex64, N>,
    scale: f64,
) -> ThreePointBlock<A, B, C, D, Complex64, N> {
    ThreePointBlock::new(
        block
            .components
            .map(|matrix| matrix.map(|row| row.map(|component| component * scale))),
    )
}

fn zero_three_point_block<A, B, C, D: Copy, const N: usize>()
-> ThreePointBlock<A, B, C, D, Complex64, N> {
    ThreePointBlock::new([[[Quantity::new(Complex64::new(0.0, 0.0)); N]; N]; N])
}

fn zero_real_three_point_block<A, B, C, D: Copy, const N: usize>()
-> ThreePointBlock<A, B, C, D, f64, N> {
    ThreePointBlock::new([[[Quantity::new(0.0); N]; N]; N])
}

/// Homogeneous three-point evolution under the linear flow `u^a_b`.
///
/// This system applies the already validated P3c linear dynamics to all three
/// Fourier legs. It intentionally omits the quadratic `u^a_bc Sigma Sigma`
/// source; that source is the next P3d layer and cannot be confused with this
/// homogeneous propagation in the public type name.
pub struct ThreePointHomogeneousTransportSystem<'a, P, G, B> {
    background: &'a CovariantMultifieldSystem<P, G, 2>,
    trajectory: &'a B,
    triangle: MomentumTriangle,
}

impl<'a, P, G, B> ThreePointHomogeneousTransportSystem<'a, P, G, B> {
    /// Builds a system whose state legs use the triangle's canonical order.
    pub const fn new(
        background: &'a CovariantMultifieldSystem<P, G, 2>,
        trajectory: &'a B,
        triangle: MomentumTriangle,
    ) -> Self {
        Self {
            background,
            trajectory,
            triangle,
        }
    }

    /// Returns the canonically ordered normalized wavenumbers used by the ODE.
    pub fn wavenumbers(&self) -> [Mass; 3] {
        self.triangle.canonical_wavenumbers().map(Mass::new)
    }
}

impl<P, G, B> System for ThreePointHomogeneousTransportSystem<'_, P, G, B>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
    B: TwoFieldBackgroundProvider,
{
    type IndependentVariable = Dimensionless;
    type State = ThreePointPhaseSpaceCorrelation<2>;
    type Delta = ThreePointPhaseSpaceCorrelation<2>;
    type Derivative = ThreePointPhaseSpaceCorrelation<2>;

    fn derivative(&self, e_folds: Dimensionless, state: &Self::State) -> Self::Derivative {
        let background = self.trajectory.state_at_e_folds(e_folds.into_value());
        let hubble = self.background.hubble(&background).into_value();
        let connection = self.background.geometry.christoffel(&background.phi);
        let mass = self.background.perturbation_mass_matrix(&background);
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
        let flows = self.wavenumbers().map(|wavenumber| {
            let physical_wavenumber = (wavenumber / background.a).into_value();
            LinearPhaseSpaceFlow {
                qq: core::array::from_fn(|upper| {
                    core::array::from_fn(|lower| -connection_flow[upper][lower])
                }),
                qp: core::array::from_fn(|upper| {
                    core::array::from_fn(|lower| f64::from(upper == lower) / hubble)
                }),
                pq: core::array::from_fn(|upper| {
                    core::array::from_fn(|lower| {
                        -(mass.components()[upper][lower].value()
                            + f64::from(upper == lower) * physical_wavenumber.powi(2))
                            / hubble
                    })
                }),
                pp: core::array::from_fn(|upper| {
                    core::array::from_fn(|lower| {
                        -connection_flow[upper][lower] - 3.0 * f64::from(upper == lower)
                    })
                }),
            }
        });

        ThreePointPhaseSpaceCorrelation {
            qqq: homogeneous_three_point_block(state, [PhaseSpaceSlotKind::Field; 3], &flows),
            qqp: homogeneous_three_point_block(
                state,
                [
                    PhaseSpaceSlotKind::Field,
                    PhaseSpaceSlotKind::Field,
                    PhaseSpaceSlotKind::Momentum,
                ],
                &flows,
            ),
            qpq: homogeneous_three_point_block(
                state,
                [
                    PhaseSpaceSlotKind::Field,
                    PhaseSpaceSlotKind::Momentum,
                    PhaseSpaceSlotKind::Field,
                ],
                &flows,
            ),
            pqq: homogeneous_three_point_block(
                state,
                [
                    PhaseSpaceSlotKind::Momentum,
                    PhaseSpaceSlotKind::Field,
                    PhaseSpaceSlotKind::Field,
                ],
                &flows,
            ),
            qpp: homogeneous_three_point_block(
                state,
                [
                    PhaseSpaceSlotKind::Field,
                    PhaseSpaceSlotKind::Momentum,
                    PhaseSpaceSlotKind::Momentum,
                ],
                &flows,
            ),
            pqp: homogeneous_three_point_block(
                state,
                [
                    PhaseSpaceSlotKind::Momentum,
                    PhaseSpaceSlotKind::Field,
                    PhaseSpaceSlotKind::Momentum,
                ],
                &flows,
            ),
            ppq: homogeneous_three_point_block(
                state,
                [
                    PhaseSpaceSlotKind::Momentum,
                    PhaseSpaceSlotKind::Momentum,
                    PhaseSpaceSlotKind::Field,
                ],
                &flows,
            ),
            ppp: homogeneous_three_point_block(state, [PhaseSpaceSlotKind::Momentum; 3], &flows),
        }
    }
}

#[derive(Clone, Copy)]
struct LinearPhaseSpaceFlow {
    qq: [[f64; 2]; 2],
    qp: [[f64; 2]; 2],
    pq: [[f64; 2]; 2],
    pp: [[f64; 2]; 2],
}

impl LinearPhaseSpaceFlow {
    fn component(
        self,
        output: PhaseSpaceSlotKind,
        input: PhaseSpaceSlotKind,
        upper: usize,
        lower: usize,
    ) -> f64 {
        match (output, input) {
            (PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Field) => self.qq[upper][lower],
            (PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Momentum) => self.qp[upper][lower],
            (PhaseSpaceSlotKind::Momentum, PhaseSpaceSlotKind::Field) => self.pq[upper][lower],
            (PhaseSpaceSlotKind::Momentum, PhaseSpaceSlotKind::Momentum) => self.pp[upper][lower],
        }
    }
}

fn homogeneous_three_point_block<A, B, C, D>(
    state: &ThreePointPhaseSpaceCorrelation<2>,
    output_slots: [PhaseSpaceSlotKind; 3],
    flows: &[LinearPhaseSpaceFlow; 3],
) -> ThreePointBlock<A, B, C, D, Complex64, 2> {
    ThreePointBlock::new(core::array::from_fn(|first| {
        core::array::from_fn(|second| {
            core::array::from_fn(|third| {
                let output_indices = [first, second, third];
                let value = (0..3)
                    .flat_map(|leg| {
                        [PhaseSpaceSlotKind::Field, PhaseSpaceSlotKind::Momentum]
                            .into_iter()
                            .flat_map(move |input_slot| {
                                (0..2).map(move |input_index| (leg, input_slot, input_index))
                            })
                    })
                    .map(|(leg, input_slot, input_index)| {
                        let mut input_slots = output_slots;
                        input_slots[leg] = input_slot;
                        let mut input_indices = output_indices;
                        input_indices[leg] = input_index;
                        flows[leg].component(
                            output_slots[leg],
                            input_slot,
                            output_indices[leg],
                            input_index,
                        ) * state.component(input_slots, input_indices)
                    })
                    .sum();
                Quantity::new(value)
            })
        })
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Field, FieldVelocity, FlatFieldSpace, HyperbolicFieldSpace, M0, Mass, MultifieldState,
        QuadraticMultifieldPotential, ScaleFactor, TwoFieldPerturbationSystem,
    };

    #[derive(Clone, Copy)]
    struct FixedTrajectory(MultifieldState<2>);

    impl TwoFieldBackgroundProvider for FixedTrajectory {
        fn state_at_e_folds(&self, _e_folds: f64) -> MultifieldState<2> {
            self.0
        }
    }

    fn sample_three_point_state() -> ThreePointPhaseSpaceCorrelation<2> {
        fn block<A, B, C, D>(offset: f64) -> ThreePointBlock<A, B, C, D, Complex64, 2> {
            ThreePointBlock::new(core::array::from_fn(|first| {
                core::array::from_fn(|second| {
                    core::array::from_fn(|third| {
                        Quantity::new(Complex64::new(
                            offset + (100 * first + 10 * second + third) as f64,
                            -offset,
                        ))
                    })
                })
            }))
        }

        ThreePointPhaseSpaceCorrelation {
            qqq: block(0.0),
            qqp: block(1_000.0),
            qpq: block(2_000.0),
            pqq: block(3_000.0),
            qpp: block(4_000.0),
            pqp: block(5_000.0),
            ppq: block(6_000.0),
            ppp: block(7_000.0),
        }
    }

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

        assert_eq!(
            triangle
                .wavevector_dot(FourierLeg::One, FourierLeg::One)
                .into_value(),
            16.0
        );
        assert_eq!(
            triangle
                .wavevector_dot(FourierLeg::One, FourierLeg::Two)
                .into_value(),
            -10.5
        );
        assert_eq!(
            triangle
                .wavevector_dot(FourierLeg::One, FourierLeg::Three)
                .into_value(),
            -5.5
        );
        assert_eq!(
            triangle
                .wavevector_dot(FourierLeg::Two, FourierLeg::Three)
                .into_value(),
            1.5
        );
        assert_eq!(
            triangle.wavevector_dot(FourierLeg::Two, FourierLeg::One),
            triangle.wavevector_dot(FourierLeg::One, FourierLeg::Two)
        );
    }

    #[test]
    fn triangle_transverse_factors_cover_equilateral_and_folded_limits() {
        let equilateral = MomentumTriangle::try_new(2.0, 2.0, 2.0).unwrap();
        let folded = MomentumTriangle::try_new(2.0, 1.0, 1.0).unwrap();

        assert_eq!(
            equilateral.wavevector_cosine(FourierLeg::One, FourierLeg::Two),
            -0.5
        );
        assert_eq!(
            equilateral.transverse_factor(FourierLeg::One, FourierLeg::Two),
            0.75
        );
        assert_eq!(
            folded.wavevector_cosine(FourierLeg::Two, FourierLeg::Three),
            1.0
        );
        assert_eq!(
            folded.transverse_factor(FourierLeg::Two, FourierLeg::Three),
            0.0
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
        let canonical = PhaseSpaceIndex::<CanonicalMomentumPerturbation, 2>::try_new(1).unwrap();

        assert_eq!(field.field_index(), 1);
        assert_eq!(field.erase().flat_index(), 1);
        assert_eq!(momentum.erase().flat_index(), 3);
        assert_eq!(canonical.erase().flat_index(), 3);
        assert_eq!(MomentumPerturbation::NAME, "momentum");
        assert_eq!(CanonicalMomentumPerturbation::NAME, "canonical_momentum");
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

    #[test]
    fn complete_three_point_state_permutations_move_blocks_and_axes() {
        let alpha = sample_three_point_state();

        let first_second = alpha.swap_first_second();
        assert_eq!(
            first_second.qpq.components()[1][0][1].into_value(),
            Complex64::new(3_011.0, -3_000.0)
        );
        assert_eq!(
            first_second.pqp.components()[1][0][1].into_value(),
            Complex64::new(4_011.0, -4_000.0)
        );
        assert_eq!(first_second.swap_first_second(), alpha);

        let second_third = alpha.swap_second_third();
        assert_eq!(
            second_third.qqp.components()[1][0][1].into_value(),
            Complex64::new(2_110.0, -2_000.0)
        );
        assert_eq!(
            second_third.pqp.components()[1][0][1].into_value(),
            Complex64::new(6_110.0, -6_000.0)
        );
        assert_eq!(second_third.swap_second_third(), alpha);
        assert_eq!(alpha.conjugate().conjugate(), alpha);
        assert_eq!(alpha + ThreePointPhaseSpaceCorrelation::zero(), alpha);
        assert_eq!(alpha * Dimensionless::new(1.0), alpha);
    }

    #[test]
    fn complete_three_point_zero_has_the_dimensioned_additive_identity() {
        let alpha = ThreePointPhaseSpaceCorrelation::<2>::zero();

        assert_eq!(
            alpha.qqq.components()[0][0][0].into_value(),
            Complex64::new(0.0, 0.0)
        );
        assert_eq!(
            alpha.ppp.components()[1][1][1].into_value(),
            Complex64::new(0.0, 0.0)
        );
    }

    #[test]
    fn homogeneous_three_point_flow_applies_the_linear_system_to_each_leg() {
        let background = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(0.0), FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(2.0e-6)]),
            FlatFieldSpace,
        );
        let trajectory = FixedTrajectory(background);
        let triangle = MomentumTriangle::try_new(3.0e-5, 2.0e-5, 2.0e-5).unwrap();
        let system = ThreePointHomogeneousTransportSystem::new(&model, &trajectory, triangle);
        let alpha = sample_three_point_state();

        let derivative = system.derivative(Dimensionless::new(0.0), &alpha);
        let hubble = model.hubble(&background).into_value();
        let expected = (alpha.pqq.components()[0][0][0].into_value()
            + alpha.qpq.components()[0][0][0].into_value()
            + alpha.qqp.components()[0][0][0].into_value())
            / hubble;
        let actual = derivative.qqq.components()[0][0][0].into_value();

        assert!((actual - expected).norm() / expected.norm() < 1.0e-15);
        assert_eq!(
            system.wavenumbers().map(Quantity::into_value),
            [3.0e-5, 2.0e-5, 2.0e-5]
        );
    }

    #[test]
    fn equilateral_homogeneous_flow_commutes_with_fourier_leg_exchange() {
        let background = MultifieldState::new(
            [Field::new(12.0), Field::new(8.0)],
            [FieldVelocity::new(-3.0e-6), FieldVelocity::new(-4.0e-6)],
            ScaleFactor::new(1.0),
        );
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-6), Mass::new(2.0e-6)]),
            FlatFieldSpace,
        );
        let trajectory = FixedTrajectory(background);
        let triangle = MomentumTriangle::try_new(2.0e-5, 2.0e-5, 2.0e-5).unwrap();
        let system = ThreePointHomogeneousTransportSystem::new(&model, &trajectory, triangle);
        let alpha = sample_three_point_state();

        let permuted_then_evolved =
            system.derivative(Dimensionless::new(0.0), &alpha.swap_first_second());
        let evolved_then_permuted = system
            .derivative(Dimensionless::new(0.0), &alpha)
            .swap_first_second();
        let slots = [
            [PhaseSpaceSlotKind::Field; 3],
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ],
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ],
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
            ],
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
            ],
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ],
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ],
            [PhaseSpaceSlotKind::Momentum; 3],
        ];

        for block_slots in slots {
            for first in 0..2 {
                for second in 0..2 {
                    for third in 0..2 {
                        let indices = [first, second, third];
                        let expected = evolved_then_permuted.component(block_slots, indices);
                        let actual = permuted_then_evolved.component(block_slots, indices);
                        let scale = expected.norm().max(actual.norm()).max(1.0);
                        assert!((actual - expected).norm() / scale < 1.0e-15);
                    }
                }
            }
        }
    }

    #[test]
    fn quadratic_kernel_source_sums_the_three_ordered_uss_terms() {
        fn covariance(qq: Complex64) -> TwoPointPhaseSpaceCovariance<1> {
            TwoPointPhaseSpaceCovariance {
                qq: TwoPointBlock::new([[Quantity::new(qq)]]),
                qp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
                pq: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
                pp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            }
        }
        fn kernel(q_qq: f64) -> QuadraticPhaseSpaceKernel<1> {
            QuadraticPhaseSpaceKernel {
                q_qq: ThreePointBlock::new([[[Quantity::new(q_qq)]]]),
                ..QuadraticPhaseSpaceKernel::zero()
            }
        }

        let kernels = [kernel(2.0), kernel(3.0), kernel(5.0)];
        let covariances = [
            covariance(Complex64::new(7.0, 1.0)),
            covariance(Complex64::new(11.0, -2.0)),
            covariance(Complex64::new(13.0, 3.0)),
        ];
        let source = three_point_source_from_kernels(
            [&kernels[0], &kernels[1], &kernels[2]],
            [&covariances[0], &covariances[1], &covariances[2]],
        );
        let expected = 2.0
            * covariances[1].qq.components()[0][0].into_value()
            * covariances[2].qq.components()[0][0].into_value()
            + 3.0
                * covariances[0].qq.components()[0][0].into_value()
                * covariances[2].qq.components()[0][0].into_value()
            + 5.0
                * covariances[0].qq.components()[0][0].into_value()
                * covariances[1].qq.components()[0][0].into_value();

        assert_eq!(source.qqq.components()[0][0][0].into_value(), expected);
        assert_eq!(
            source.ppp.components()[0][0][0].into_value(),
            Complex64::new(0.0, 0.0)
        );
    }

    #[test]
    fn canonical_kernel_contracts_only_with_canonical_covariances() {
        let kernel = CanonicalQuadraticPhaseSpaceKernel::<1>::zero();
        let covariance = CanonicalTwoPointPhaseSpaceCovariance {
            qq: TwoPointBlock::new([[Quantity::new(Complex64::new(1.0, 0.0))]]),
            qp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            pq: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            pp: TwoPointBlock::new([[Quantity::new(Complex64::new(1.0, 0.0))]]),
        };

        let source: CanonicalThreePointPhaseSpaceCorrelation<1> = three_point_source_from_kernels(
            [&kernel, &kernel, &kernel],
            [&covariance, &covariance, &covariance],
        );

        assert_eq!(source, CanonicalThreePointPhaseSpaceCorrelation::zero());
    }

    #[test]
    fn cubic_hamiltonian_map_preserves_raised_index_order_and_divides_by_hubble() {
        fn block<D: Copy>(
            offset: f64,
        ) -> ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, D, f64, 2>
        {
            ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| Quantity::new(offset + (100 * i + 10 * j + k) as f64))
                })
            }))
        }

        let coefficients = CanonicalCubicHamiltonianCoefficients {
            a_first_raised: block::<M1>(1_000.0),
            b_first_raised: block::<M0>(2_000.0),
            b_third_raised: block::<M0>(3_000.0),
            c_first_raised: block::<MMinus1>(4_000.0),
            c_third_raised: block::<MMinus1>(5_000.0),
        };
        let kernel = CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
            coefficients,
            Hubble::new(2.0),
        )
        .unwrap();

        assert_eq!(*kernel.q_qq.components()[1][0][1].value(), -3_101.0 / 2.0);
        assert_eq!(*kernel.q_qp.components()[1][0][1].value(), -4_110.0 / 2.0);
        assert_eq!(*kernel.q_pq.components()[1][0][1].value(), -4_101.0 / 2.0);
        assert_eq!(*kernel.q_pp.components()[1][0][1].value(), 0.0);
        assert_eq!(
            *kernel.p_qq.components()[1][0][1].value(),
            3.0 * 1_101.0 / 2.0
        );
        assert_eq!(*kernel.p_qp.components()[1][0][1].value(), 2_101.0 / 2.0);
        assert_eq!(*kernel.p_pq.components()[1][0][1].value(), 2_110.0 / 2.0);
        assert_eq!(*kernel.p_pp.components()[1][0][1].value(), 5_110.0 / 2.0);
    }

    #[test]
    fn cubic_hamiltonian_map_rejects_nonphysical_hubble_rates() {
        let coefficients = CanonicalCubicHamiltonianCoefficients::<1> {
            a_first_raised: zero_real_three_point_block(),
            b_first_raised: zero_real_three_point_block(),
            b_third_raised: zero_real_three_point_block(),
            c_first_raised: zero_real_three_point_block(),
            c_third_raised: zero_real_three_point_block(),
        };

        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
                    coefficients,
                    Hubble::new(invalid),
                )
                .is_err()
            );
        }
    }

    #[test]
    fn canonical_cubic_c_matches_the_flat_single_field_limit() {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(2.0)],
            [FieldVelocity::new(0.5)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(1.0, 1.0, 1.0).unwrap();

        let coefficients = canonical_cubic_c_coefficients(&model, &state, triangle).unwrap();
        let hubble = model.hubble(&state).into_value();
        let velocity: f64 = 0.5;
        let expected = -velocity / hubble + 3.0 * velocity.powi(3) / (32.0 * hubble.powi(3));

        for actual in [
            coefficients.first_raised.components()[0][0][0].into_value(),
            coefficients.third_raised.components()[0][0][0].into_value(),
        ] {
            let scale = actual.abs().max(expected.abs()).max(1.0);
            assert!((actual - expected).abs() / scale < 1.0e-14);
        }
    }

    #[test]
    fn canonical_cubic_c_preserves_the_two_raised_index_layouts() {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0), Mass::new(2.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(1.0)],
            [FieldVelocity::new(0.5), FieldVelocity::new(-0.25)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();

        let coefficients = canonical_cubic_c_coefficients(&model, &state, triangle).unwrap();

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let actual = coefficients.third_raised.components()[i][j][k].into_value();
                    let expected = coefficients.first_raised.components()[j][k][i].into_value();
                    let scale = actual.abs().max(expected.abs()).max(1.0);
                    assert!((actual - expected).abs() / scale < 1.0e-14);
                }
            }
        }
    }

    #[test]
    fn canonical_cubic_c_rejects_nonphysical_background_scales() {
        let potential = QuadraticMultifieldPotential::new([Mass::new(1.0)]);
        let state = MultifieldState::new(
            [Field::new(1.0)],
            [FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(1.0, 1.0, 1.0).unwrap();
        let invalid_planck =
            CovariantMultifieldSystem::new(potential, FlatFieldSpace, Mass::new(0.0));

        assert_eq!(
            canonical_cubic_c_coefficients(&invalid_planck, &state, triangle),
            Err(CubicCoefficientError::InvalidPlanckMass(0.0))
        );
        assert_eq!(
            canonical_cubic_b_coefficients(&invalid_planck, &state, triangle),
            Err(CubicCoefficientError::InvalidPlanckMass(0.0))
        );

        let model = CovariantMultifieldSystem::in_planck_units(potential, FlatFieldSpace);
        let zero_energy = MultifieldState::new(
            [Field::new(0.0)],
            [FieldVelocity::new(0.0)],
            ScaleFactor::new(1.0),
        );
        assert_eq!(
            canonical_cubic_c_coefficients(&model, &zero_energy, triangle),
            Err(CubicCoefficientError::InvalidHubble(0.0))
        );
        assert_eq!(
            canonical_cubic_b_coefficients(&model, &zero_energy, triangle),
            Err(CubicCoefficientError::InvalidHubble(0.0))
        );
    }

    #[test]
    fn canonical_cubic_b_matches_the_flat_single_field_limit() {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(2.0)],
            [FieldVelocity::new(0.5)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(1.0, 1.0, 1.0).unwrap();

        let coefficients = canonical_cubic_b_coefficients(&model, &state, triangle).unwrap();
        let hubble = model.hubble(&state).into_value();
        let velocity: f64 = 0.5;
        let gradient = 2.0;
        let z = -3.0 * hubble * velocity - gradient + velocity.powi(3) / (2.0 * hubble);
        let expected = -3.0 * velocity.powi(2) * z / (16.0 * hubble.powi(3))
            + velocity.powi(3) / (4.0 * hubble.powi(2))
            + z / (2.0 * hubble);

        for actual in [
            coefficients.first_raised.components()[0][0][0].into_value(),
            coefficients.third_raised.components()[0][0][0].into_value(),
        ] {
            let scale = actual.abs().max(expected.abs()).max(1.0);
            assert!((actual - expected).abs() / scale < 1.0e-14);
        }
    }

    #[test]
    fn canonical_cubic_b_preserves_the_two_raised_index_layouts() {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0), Mass::new(2.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(1.0)],
            [FieldVelocity::new(0.5), FieldVelocity::new(-0.25)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();

        let coefficients = canonical_cubic_b_coefficients(&model, &state, triangle).unwrap();

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let actual = coefficients.third_raised.components()[i][j][k].into_value();
                    let expected = coefficients.first_raised.components()[j][k][i].into_value();
                    let scale = actual.abs().max(expected.abs()).max(1.0);
                    assert!((actual - expected).abs() / scale < 1.0e-14);
                }
            }
        }
    }

    #[test]
    fn canonical_cubic_b_contains_the_covariant_field_space_curvature_term() {
        #[derive(Clone, Copy)]
        struct WithoutCurvature(HyperbolicFieldSpace);

        impl FieldSpaceGeometry<2> for WithoutCurvature {
            fn metric(&self, position: &[Field; 2]) -> crate::FieldMetric<2> {
                self.0.metric(position)
            }

            fn inverse_metric(&self, position: &[Field; 2]) -> crate::InverseFieldMetric<2> {
                self.0.inverse_metric(position)
            }

            fn christoffel(&self, position: &[Field; 2]) -> crate::ChristoffelSymbols<2> {
                self.0.christoffel(position)
            }
        }

        impl CurvedFieldSpaceGeometry<2> for WithoutCurvature {
            fn riemann(&self, _position: &[Field; 2]) -> crate::RiemannTensor<2> {
                crate::RiemannTensor::zero()
            }
        }

        let geometry = HyperbolicFieldSpace::new(Field::new(3.0)).unwrap();
        let potential = QuadraticMultifieldPotential::new([Mass::new(1.0), Mass::new(0.75)]);
        let curved = CovariantMultifieldSystem::in_planck_units(potential, geometry);
        let suppressed =
            CovariantMultifieldSystem::in_planck_units(potential, WithoutCurvature(geometry));
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(0.5)],
            [FieldVelocity::new(0.4), FieldVelocity::new(-0.3)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();

        let actual = canonical_cubic_b_coefficients(&curved, &state, triangle).unwrap();
        let without = canonical_cubic_b_coefficients(&suppressed, &state, triangle).unwrap();
        let metric = geometry.metric(&state.phi);
        let inverse = geometry.inverse_metric(&state.phi);
        let riemann = metric.lower_first(&geometry.riemann(&state.phi));

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let expected = (0..2)
                        .map(|contracted| {
                            let lower_curvature = (0..2)
                                .map(|l| {
                                    (riemann.components()[k][contracted][j][l]
                                        + riemann.components()[k][j][contracted][l])
                                        .into_value()
                                        * state.pi[l].into_value()
                                        * (2.0 / 3.0)
                                })
                                .sum::<f64>();
                            inverse.components()[i][contracted].into_value() * lower_curvature
                        })
                        .sum::<f64>();
                    let difference = actual.first_raised.components()[i][j][k].into_value()
                        - without.first_raised.components()[i][j][k].into_value();
                    let scale = difference.abs().max(expected.abs()).max(1.0);
                    assert!((difference - expected).abs() / scale < 1.0e-14);
                }
            }
        }
    }
}
