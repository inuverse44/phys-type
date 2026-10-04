//! Typed kinematics and correlation containers for phase-space transport.
//!
//! This module separates Fourier legs, phase-space slots, field coordinates,
//! and physical dimensions. It provides the linear mode solver's two-point
//! oracle, direct two-point covariance dynamics, and typed three-point
//! transport infrastructure. Model-derived cubic coefficients and the
//! subhorizon tree-level three-point initial condition are complete; sourced
//! end-to-end evolution is not implemented yet.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Mul};

use num_complex::Complex64;

use crate::integrate::{AdaptiveSystem, ErrorTolerances, System};
use crate::{
    ContravariantField, CovariantMultifieldSystem, CurvatureDerivativeGeometry,
    CurvedFieldSpaceGeometry, Dimensionless, FieldSpaceGeometry, Hubble, LowerFirst, M0, M1, M2,
    M3, M4, MMinus1, MMinus2, MMinus3, Mass, ModeAmplitude, ModeVelocity, MultifieldPotential,
    MultifieldState, Quantity, ThirdOrderFieldSpaceGeometry,
    ThriceDifferentiableMultifieldPotential, TwiceDifferentiableMultifieldPotential,
    TwoFieldPerturbationState,
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

/// Relabels a velocity covariance as canonical at linear perturbative order.
///
/// Since `Pi^I = V^I + O(X^2)`, all four two-point blocks are numerically
/// identical at the order used to source tree-level three-point transport.
/// The explicit function keeps this perturbative boundary visible in the API.
pub fn canonical_two_point_from_velocity_linear<const N: usize>(
    covariance: TwoPointPhaseSpaceCovariance<N>,
) -> CanonicalTwoPointPhaseSpaceCovariance<N> {
    CanonicalTwoPointPhaseSpaceCovariance {
        qq: covariance.qq,
        qp: TwoPointBlock::new(covariance.qp.components),
        pq: TwoPointBlock::new(covariance.pq.components),
        pp: TwoPointBlock::new(covariance.pp.components),
    }
}

/// Relabels a canonical covariance as velocity at linear perturbative order.
pub fn velocity_two_point_from_canonical_linear<const N: usize>(
    covariance: CanonicalTwoPointPhaseSpaceCovariance<N>,
) -> TwoPointPhaseSpaceCovariance<N> {
    TwoPointPhaseSpaceCovariance {
        qq: covariance.qq,
        qp: TwoPointBlock::new(covariance.qp.components),
        pq: TwoPointBlock::new(covariance.pq.components),
        pp: TwoPointBlock::new(covariance.pp.components),
    }
}

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
        two_point_linear_derivative(self.background, &background, self.k, state)
    }
}

/// Linear covariance flow in elapsed e-folds on one background state.
///
/// At linear order the velocity `V` and canonical momentum `Pi` obey the same
/// equations, so the marker `S` is carried through unchanged.
fn two_point_linear_derivative<P, G, S: PhaseSpaceSlot>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    background: &MultifieldState<2>,
    k: Mass,
    state: &TwoPointPhaseSpaceCovariance<2, S>,
) -> TwoPointPhaseSpaceCovariance<2, S>
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    let hubble = system.hubble(background).into_value();
    let connection = system.geometry.christoffel(&background.phi);
    let mass = system.perturbation_mass_matrix(background);
    let physical_wavenumber = (k / background.a).into_value();
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
                            + restoring[right][index] * state.qq.components()[left][index].value(),
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
                            + restoring[left][index] * state.qq.components()[index][right].value(),
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
                            + restoring[left][index] * state.qp.components()[index][right].value()
                            + restoring[right][index] * state.pq.components()[left][index].value(),
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
    /// `B^I_KJ(k_out; k_in2, k_in1)`, stored as `[I][J][K]`.
    ///
    /// The second covariant slot carries index `K` and the second input
    /// momentum, so this is a separate evaluation rather than a transpose of
    /// `b_first_raised`.
    pub b_first_raised_input_exchanged:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `B_JK^I`, with the third original index raised, stored as `[I][J][K]`.
    pub b_third_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `C^I_JK`, with the first original index raised, dimension `M^-1`.
    pub c_first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
    /// `C^I_KJ(k_out; k_in2, k_in1)`, stored as `[I][J][K]`.
    pub c_first_raised_input_exchanged:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
    /// `C_JK^I`, with the third original index raised, stored as `[I][J][K]`.
    pub c_third_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>,
}

/// Cubic-Hamiltonian coefficients ordered for each possible output leg.
///
/// For output legs `1`, `2`, and `3`, the two input legs are respectively
/// `(2, 3)`, `(1, 3)`, and `(1, 2)`. Keeping this ordering in a dedicated type
/// prevents a single `(1, 2, 3)` vertex from being reused incorrectly in all
/// three cyclic source terms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalCubicHamiltonianVertices<const N: usize> {
    coefficients: [CanonicalCubicHamiltonianCoefficients<N>; 3],
}

impl<const N: usize> CanonicalCubicHamiltonianVertices<N> {
    /// Returns the coefficients whose first slot is the requested output leg.
    pub const fn for_output(
        &self,
        output: FourierLeg,
    ) -> &CanonicalCubicHamiltonianCoefficients<N> {
        &self.coefficients[output.index()]
    }

    /// Converts all three ordered vertices to e-fold quadratic-flow kernels.
    pub fn e_fold_kernels(
        self,
        hubble: Hubble,
    ) -> Result<[CanonicalQuadraticPhaseSpaceKernel<N>; 3], CanonicalKernelError> {
        Ok([
            CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
                self.coefficients[0],
                hubble,
            )?,
            CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
                self.coefficients[1],
                hubble,
            )?,
            CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
                self.coefficients[2],
                hubble,
            )?,
        ])
    }
}

/// The raised form of the fully symmetrized cubic coefficient `A_IJK`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalCubicACoefficients<const N: usize> {
    /// `A^I_JK`, stored in output/first-input/second-input order.
    pub first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M1, f64, N>,
}

/// The two raised forms of the symmetrized cubic coefficient `B_IJK`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalCubicBCoefficients<const N: usize> {
    /// `B^I_JK`, stored in output/first-input/second-input order.
    pub first_raised:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>,
    /// `B^I_KJ` with its two input momenta exchanged, stored as `[I][J][K]`.
    pub first_raised_input_exchanged:
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
    /// `C^I_KJ` with its two input momenta exchanged, stored as `[I][J][K]`.
    pub first_raised_input_exchanged:
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
    /// The scale factor needed by a physical-wavenumber term was not positive
    /// and finite.
    InvalidScaleFactor(f64),
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
            Self::InvalidScaleFactor(value) => write!(
                formatter,
                "cubic coefficients require a positive finite scale factor, got {value}"
            ),
        }
    }
}

impl std::error::Error for CubicCoefficientError {}

/// Evaluates the fully symmetrized curved-field-space coefficient `A_IJK` and
/// returns the raised form required by the canonical phase-space kernel.
///
/// This is the coefficient that exercises all third-order geometry
/// capabilities: the symmetrized covariant potential derivative, the Riemann
/// tensor, and its covariant derivative. Each permutation moves a field index
/// and its attached canonical Fourier leg together.
fn canonical_cubic_a_coefficients_for_legs<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
    legs: [FourierLeg; 3],
) -> Result<CanonicalCubicACoefficients<N>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
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
    let scale_value = state.a.into_value();
    if !scale_value.is_finite() || scale_value <= 0.0 {
        return Err(CubicCoefficientError::InvalidScaleFactor(scale_value));
    }

    let metric = system.geometry.metric(&state.phi);
    let inverse = system.geometry.inverse_metric(&state.phi);
    let lowered_riemann = metric.lower_first(&system.geometry.riemann(&state.phi));
    let riemann_derivative = system.geometry.covariant_riemann_derivative(&state.phi);
    let potential_hessian = system.covariant_potential_hessian(state);
    let potential_third = system.covariant_potential_third_derivative(state);
    let velocity = ContravariantField::new(state.pi);
    let lowered_velocity = metric.lower(&velocity);
    let velocity_squared = lowered_velocity.contract(&velocity);
    let gradient = system.potential.gradient(state.phi);
    let planck = Mass::new(planck_value);
    let planck_squared = planck * planck;
    let planck_fourth = planck_squared * planck_squared;
    let planck_sixth = planck_fourth * planck_squared;
    let scale = crate::ScaleFactor::new(scale_value);
    let scale_squared = scale * scale;
    let z: [Quantity<M3, f64>; N] = core::array::from_fn(|i| {
        -(lowered_velocity.components()[i] * hubble * 3.0) - gradient[i]
            + lowered_velocity.components()[i] * velocity_squared / (hubble * planck_squared * 2.0)
    });
    let unsymmetrized = |i: usize,
                         j: usize,
                         k: usize,
                         _first_leg: FourierLeg,
                         second_leg: FourierLeg,
                         third_leg: FourierLeg| {
        let velocity_i = lowered_velocity.components()[i];
        let velocity_j = lowered_velocity.components()[j];
        let velocity_k = lowered_velocity.components()[k];
        let riemann_term = (0..N)
            .flat_map(|l| (0..N).map(move |m| (l, m)))
            .map(|(l, m)| {
                velocity_k
                    * velocity.components()[l]
                    * velocity.components()[m]
                    * (lowered_riemann.components()[l][i][j][m]
                        + lowered_riemann.components()[l][j][i][m])
                    / (hubble * planck_squared)
                    * -0.25
            })
            .fold(Quantity::new(0.0), Add::add);
        let riemann_derivative_term = (0..N)
            .flat_map(|upper| (0..N).flat_map(move |l| (0..N).map(move |m| (upper, l, m))))
            .map(|(upper, l, m)| {
                riemann_derivative.components()[upper][l][m][j][k]
                    * metric.components()[i][upper].into_value()
                    * velocity.components()[l]
                    * velocity.components()[m]
                    / 3.0
            })
            .fold(Quantity::new(0.0), Add::add);
        let physical_dot = triangle.wavevector_dot(second_leg, third_leg) / scale_squared;

        -(potential_third.components()[i][j][k] / 3.0)
            - velocity_i * potential_hessian.components()[j][k] / (hubble * planck_squared * 2.0)
            + velocity_i * velocity_j * z[k] / (hubble * hubble * planck_fourth * 4.0)
            + velocity_i * z[j] * z[k] / (hubble * hubble * hubble * planck_fourth)
                * (triangle.transverse_factor(second_leg, third_leg) / 8.0)
            + velocity_i * velocity_j * velocity_k / (hubble * hubble * hubble * planck_sixth)
                * (hubble * hubble * planck_squared * 6.0 - velocity_squared)
                / 8.0
            + riemann_term
            + riemann_derivative_term
            + velocity_i * physical_dot * metric.components()[j][k].into_value()
                / (hubble * planck_squared * 2.0)
    };

    let lower: [[[Quantity<M1, f64>; N]; N]; N] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (unsymmetrized(i, j, k, legs[0], legs[1], legs[2])
                    + unsymmetrized(i, k, j, legs[0], legs[2], legs[1])
                    + unsymmetrized(j, i, k, legs[1], legs[0], legs[2])
                    + unsymmetrized(j, k, i, legs[1], legs[2], legs[0])
                    + unsymmetrized(k, i, j, legs[2], legs[0], legs[1])
                    + unsymmetrized(k, j, i, legs[2], legs[1], legs[0]))
                    / 6.0
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

    Ok(CanonicalCubicACoefficients { first_raised })
}

/// Evaluates `A` in canonical Fourier-leg order `(k_1, k_2, k_3)`.
pub fn canonical_cubic_a_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicACoefficients<N>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    canonical_cubic_a_coefficients_for_legs(
        system,
        state,
        triangle,
        [FourierLeg::One, FourierLeg::Two, FourierLeg::Three],
    )
}

/// Evaluates the symmetrized curved-field-space coefficient `B_IJK` and
/// returns the two raised forms required by the canonical phase-space kernel.
///
/// The background equation of motion is used covariantly to form
/// `Z_I = D_t dot(phi)_I + dot(phi)_I dot(phi)^2 / (2 H M_pl^2)`.
/// Field indices and their attached canonical Fourier legs are symmetrized
/// together.
fn canonical_cubic_b_coefficients_for_legs<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
    legs: [FourierLeg; 3],
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

    let lower_for_legs = |coefficient_legs: [FourierLeg; 3]| {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    (unsymmetrized(
                        i,
                        j,
                        k,
                        coefficient_legs[0],
                        coefficient_legs[1],
                        coefficient_legs[2],
                    ) + unsymmetrized(
                        j,
                        i,
                        k,
                        coefficient_legs[1],
                        coefficient_legs[0],
                        coefficient_legs[2],
                    )) * 0.5
                })
            })
        })
    };
    let first_lower: [[[Quantity<M0, f64>; N]; N]; N] = lower_for_legs(legs);
    let exchanged_lower: [[[Quantity<M0, f64>; N]; N]; N] =
        lower_for_legs([legs[0], legs[2], legs[1]]);
    let third_lower: [[[Quantity<M0, f64>; N]; N]; N] = lower_for_legs([legs[1], legs[2], legs[0]]);

    let first_raised = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        first_lower[contracted][j][k]
                            * inverse.components()[i][contracted].into_value()
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
                        third_lower[j][k][contracted]
                            * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));

    let first_raised_input_exchanged = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        exchanged_lower[contracted][k][j]
                            * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));

    Ok(CanonicalCubicBCoefficients {
        first_raised,
        first_raised_input_exchanged,
        third_raised,
    })
}

/// Evaluates `B` in canonical Fourier-leg order `(k_1, k_2, k_3)`.
pub fn canonical_cubic_b_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicBCoefficients<N>, CubicCoefficientError>
where
    P: MultifieldPotential<N>,
    G: CurvedFieldSpaceGeometry<N>,
{
    canonical_cubic_b_coefficients_for_legs(
        system,
        state,
        triangle,
        [FourierLeg::One, FourierLeg::Two, FourierLeg::Three],
    )
}

/// Evaluates the symmetrized curved-field-space coefficient `C_IJK` and
/// returns the two raised forms required by the canonical phase-space kernel.
///
/// Field indices and their attached canonical Fourier legs are symmetrized
/// together. All arithmetic that carries physical dimension remains inside
/// [`Quantity`]; only metric components and momentum ratios are raw
/// dimensionless factors.
fn canonical_cubic_c_coefficients_for_legs<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
    legs: [FourierLeg; 3],
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

    let lower_for_legs = |coefficient_legs: [FourierLeg; 3]| {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    (unsymmetrized(
                        i,
                        j,
                        k,
                        coefficient_legs[0],
                        coefficient_legs[1],
                        coefficient_legs[2],
                    ) + unsymmetrized(
                        j,
                        i,
                        k,
                        coefficient_legs[1],
                        coefficient_legs[0],
                        coefficient_legs[2],
                    )) * 0.5
                })
            })
        })
    };
    let first_lower: [[[Quantity<MMinus1, f64>; N]; N]; N] = lower_for_legs(legs);
    let exchanged_lower: [[[Quantity<MMinus1, f64>; N]; N]; N] =
        lower_for_legs([legs[0], legs[2], legs[1]]);
    let third_lower: [[[Quantity<MMinus1, f64>; N]; N]; N] =
        lower_for_legs([legs[1], legs[2], legs[0]]);

    let first_raised = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        first_lower[contracted][j][k]
                            * inverse.components()[i][contracted].into_value()
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
                        third_lower[j][k][contracted]
                            * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));

    let first_raised_input_exchanged = ThreePointBlock::new(core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .map(|contracted| {
                        exchanged_lower[contracted][k][j]
                            * inverse.components()[i][contracted].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    }));

    Ok(CanonicalCubicCCoefficients {
        first_raised,
        first_raised_input_exchanged,
        third_raised,
    })
}

/// Evaluates `C` in canonical Fourier-leg order `(k_1, k_2, k_3)`.
pub fn canonical_cubic_c_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicCCoefficients<N>, CubicCoefficientError>
where
    P: MultifieldPotential<N>,
    G: FieldSpaceGeometry<N>,
{
    canonical_cubic_c_coefficients_for_legs(
        system,
        state,
        triangle,
        [FourierLeg::One, FourierLeg::Two, FourierLeg::Three],
    )
}

/// Evaluates all model-derived cubic Hamiltonian coefficients for one
/// background state and canonical momentum triangle.
///
/// The returned tensors can be passed directly to
/// [`CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian`].
fn canonical_cubic_hamiltonian_coefficients_for_legs<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
    legs: [FourierLeg; 3],
) -> Result<CanonicalCubicHamiltonianCoefficients<N>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    let a = canonical_cubic_a_coefficients_for_legs(system, state, triangle, legs)?;
    let b = canonical_cubic_b_coefficients_for_legs(system, state, triangle, legs)?;
    let c = canonical_cubic_c_coefficients_for_legs(system, state, triangle, legs)?;

    Ok(CanonicalCubicHamiltonianCoefficients {
        a_first_raised: a.first_raised,
        b_first_raised: b.first_raised,
        b_first_raised_input_exchanged: b.first_raised_input_exchanged,
        b_third_raised: b.third_raised,
        c_first_raised: c.first_raised,
        c_first_raised_input_exchanged: c.first_raised_input_exchanged,
        c_third_raised: c.third_raised,
    })
}

/// Evaluates all cubic coefficients in canonical Fourier-leg order.
pub fn canonical_cubic_hamiltonian_coefficients<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicHamiltonianCoefficients<N>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    canonical_cubic_hamiltonian_coefficients_for_legs(
        system,
        state,
        triangle,
        [FourierLeg::One, FourierLeg::Two, FourierLeg::Three],
    )
}

/// Evaluates the three ordered vertices needed by the cyclic three-point
/// transport source.
pub fn canonical_cubic_hamiltonian_vertices<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalCubicHamiltonianVertices<N>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    Ok(CanonicalCubicHamiltonianVertices {
        coefficients: [
            canonical_cubic_hamiltonian_coefficients_for_legs(
                system,
                state,
                triangle,
                [FourierLeg::One, FourierLeg::Two, FourierLeg::Three],
            )?,
            canonical_cubic_hamiltonian_coefficients_for_legs(
                system,
                state,
                triangle,
                [FourierLeg::Two, FourierLeg::One, FourierLeg::Three],
            )?,
            canonical_cubic_hamiltonian_coefficients_for_legs(
                system,
                state,
                triangle,
                [FourierLeg::Three, FourierLeg::One, FourierLeg::Two],
            )?,
        ],
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
    /// ordered input exchanges required by the Hamilton equations. An input
    /// exchange moves the field index and its Fourier momentum together:
    /// `q_qp` and `p_pq` read the separately evaluated input-exchanged `C` and
    /// `B` tensors, while `p_pp` uses `C_KJ^I(k_2, k_1) = C_JK^I(k_1, k_2)`.
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
                        -(coefficients.c_first_raised_input_exchanged.components[i][j][k] / hubble)
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
                        coefficients.b_first_raised_input_exchanged.components[i][j][k] / hubble
                    })
                })
            })),
            p_pp: ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        coefficients.c_third_raised.components[i][j][k] / hubble
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

/// Selects the equal-time Wick contraction in operator order.
///
/// An external leg standing to the left of the quadratic insertion contracts
/// as `Sigma^{a d}`; one standing to its right contracts as `Sigma^{d a}`.
/// The two differ by the commutator carried in the imaginary part.
fn ordered_covariance<const N: usize, P: PhaseSpaceSlot>(
    covariance: &TwoPointPhaseSpaceCovariance<N, P>,
    external_is_left: bool,
    internal: (PhaseSpaceSlotKind, usize),
    external: (PhaseSpaceSlotKind, usize),
) -> Complex64 {
    if external_is_left {
        covariance.component([external.0, internal.0], [external.1, internal.1])
    } else {
        covariance.component([internal.0, external.0], [internal.1, external.1])
    }
}

/// Contracts three ordered quadratic kernels with the three leg covariances.
///
/// `kernels[a]` has output leg `a`. Its two input slots follow the remaining
/// legs in ascending order: `[1, 2]`, `[0, 2]`, and `[0, 1]`. The returned
/// value is exactly the three cyclic `u Sigma Sigma` terms in `D_N Alpha`; it
/// does not include homogeneous propagation or an initial condition.
///
/// The covariances keep the operator order of `<X^a X^b X^c>`:
/// `u^a_de Sigma^db Sigma^ec + u^b_de Sigma^ad Sigma^ec + u^c_de Sigma^ad
/// Sigma^be`. Transposing a covariance on a leg to the left of the insertion
/// would conjugate its commutator part.
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
                                        ) * ordered_covariance(
                                            covariances[input_legs[0]],
                                            input_legs[0] < output_leg,
                                            (first_input_slot, first_input_index),
                                            (
                                                output_slots[input_legs[0]],
                                                output_indices[input_legs[0]],
                                            ),
                                        ) * ordered_covariance(
                                            covariances[input_legs[1]],
                                            input_legs[1] < output_leg,
                                            (second_input_slot, second_input_index),
                                            (
                                                output_slots[input_legs[1]],
                                                output_indices[input_legs[1]],
                                            ),
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

fn canonical_to_velocity_redefinition_kernel<const N: usize>(
    coefficients: CanonicalCubicHamiltonianCoefficients<N>,
) -> CanonicalQuadraticPhaseSpaceKernel<N> {
    CanonicalQuadraticPhaseSpaceKernel {
        q_qq: zero_real_three_point_block(),
        q_qp: zero_real_three_point_block(),
        q_pq: zero_real_three_point_block(),
        q_pp: zero_real_three_point_block(),
        p_qq: ThreePointBlock::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| -coefficients.b_third_raised.components[i][j][k])
            })
        })),
        p_qp: ThreePointBlock::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    -coefficients.c_first_raised_input_exchanged.components[i][j][k]
                })
            })
        })),
        p_pq: ThreePointBlock::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| -coefficients.c_first_raised.components[i][j][k])
            })
        })),
        p_pp: zero_real_three_point_block(),
    }
}

/// Converts a canonical equal-time three-point function to the velocity basis
/// at tree level.
///
/// The nonlinear relation used here is
/// `V^I = Pi^I - B_JK^I Q^J Q^K / 2 - C^I_JK Pi^J Q^K + O(X^3)`.
/// Each quadratic insertion is contracted with the two linear covariances on
/// the other legs. The ordered vertices and covariances must describe the same
/// background state and canonical momentum triangle.
pub fn velocity_three_point_from_canonical_tree<const N: usize>(
    canonical: CanonicalThreePointPhaseSpaceCorrelation<N>,
    vertices: &CanonicalCubicHamiltonianVertices<N>,
    covariances: [&CanonicalTwoPointPhaseSpaceCovariance<N>; 3],
) -> ThreePointPhaseSpaceCorrelation<N> {
    let kernels = [
        canonical_to_velocity_redefinition_kernel(*vertices.for_output(FourierLeg::One)),
        canonical_to_velocity_redefinition_kernel(*vertices.for_output(FourierLeg::Two)),
        canonical_to_velocity_redefinition_kernel(*vertices.for_output(FourierLeg::Three)),
    ];
    let transformed = canonical
        + three_point_source_from_kernels([&kernels[0], &kernels[1], &kernels[2]], covariances);

    ThreePointPhaseSpaceCorrelation {
        qqq: transformed.qqq,
        qqp: ThreePointBlock::new(transformed.qqp.components),
        qpq: ThreePointBlock::new(transformed.qpq.components),
        pqq: ThreePointBlock::new(transformed.pqq.components),
        qpp: ThreePointBlock::new(transformed.qpp.components),
        pqp: ThreePointBlock::new(transformed.pqp.components),
        ppq: ThreePointBlock::new(transformed.ppq.components),
        ppp: ThreePointBlock::new(transformed.ppp.components),
    }
}

/// Fully raised cubic tensors consumed by the subhorizon three-point initial
/// condition.
///
/// Every array is stored as `[I][J][K]`, where `I`, `J`, and `K` belong to
/// Fourier legs one, two, and three of the correlation being initialized. This
/// is the index order of [`ThreePointPhaseSpaceCorrelation`], not the
/// output/input order of a quadratic kernel. The leg-resolved `B` and `C`
/// tensors are the three differently ordered evaluations required by the
/// symmetrized interaction vertex; they are not index permutations of one
/// array.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalSubhorizonInitialVertices<const N: usize> {
    /// `a^2 A_fast^{IJK}`, the gradient term of `A` with its explicit
    /// `a^{-2}` removed, dimension `M^1`.
    pub fast_a_times_scale_squared:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M1, f64, N>,
    /// `A^{IJK} - A_fast^{IJK}`, every term of `A` without a physical
    /// wavenumber, dimension `M^1`.
    pub slow_a:
        ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M1, f64, N>,
    /// Entry `r` is `B` with its momentum slot on leg `r`, dimension `M^0`.
    ///
    /// For example, entry three is `B^{IJK}(k_1, k_2, k_3)` and entry one is
    /// `B^{JKI}(k_2, k_3, k_1)`, both stored at `[I][J][K]`.
    pub b_momentum_on_leg:
        [ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, M0, f64, N>; 3],
    /// Entry `r` is `C` with its field slot on leg `r`, dimension `M^-1`.
    pub c_field_on_leg:
        [ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, MMinus1, f64, N>;
            3],
}

/// Why a subhorizon three-point initial condition could not be formed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThreePointInitialConditionError {
    /// The required minimum `k/(aH)` was not finite and greater than one.
    InvalidSubhorizonRequirement(f64),
    /// The shortest canonical leg is not deep enough inside the horizon.
    InsufficientSubhorizonDepth {
        /// Canonical leg that failed the requirement.
        leg: FourierLeg,
        /// Its `k/(aH)` on the requested slice.
        k_over_ah: f64,
        /// The requested minimum.
        required: f64,
    },
    /// The scale factor was not positive and finite.
    InvalidScaleFactor(f64),
    /// The Hubble rate was not positive and finite.
    InvalidHubble(f64),
    /// A wavenumber was not positive and finite.
    InvalidWavenumber {
        /// Leg carrying the rejected value.
        leg: FourierLeg,
        /// Rejected value in reduced Planck units.
        value: f64,
    },
    /// Evaluating the model-derived cubic coefficients failed.
    Coefficient(CubicCoefficientError),
    /// A leg could not be placed in the Bunch--Davies vacuum.
    ModeInitialization(crate::ModeInitializationError),
    /// The finite-difference spacing was not positive and finite.
    InvalidEFoldStep(f64),
}

impl fmt::Display for ThreePointInitialConditionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSubhorizonRequirement(value) => write!(
                formatter,
                "a subhorizon initial slice requires a finite minimum k/(aH) above one, got {value}"
            ),
            Self::InsufficientSubhorizonDepth {
                leg,
                k_over_ah,
                required,
            } => write!(
                formatter,
                "leg {leg:?} has k/(aH)={k_over_ah:e}, below the required {required:e}"
            ),
            Self::InvalidScaleFactor(value) => write!(
                formatter,
                "the initial slice requires a positive finite scale factor, got {value}"
            ),
            Self::InvalidHubble(value) => write!(
                formatter,
                "the initial slice requires a positive finite Hubble rate, got {value}"
            ),
            Self::InvalidWavenumber { leg, value } => write!(
                formatter,
                "leg {leg:?} requires a positive finite wavenumber, got {value}"
            ),
            Self::Coefficient(error) => write!(formatter, "{error}"),
            Self::ModeInitialization(error) => write!(formatter, "{error:?}"),
            Self::InvalidEFoldStep(value) => write!(
                formatter,
                "vertex rates require a positive finite e-fold step, got {value}"
            ),
        }
    }
}

impl std::error::Error for ThreePointInitialConditionError {}

impl From<CubicCoefficientError> for ThreePointInitialConditionError {
    fn from(error: CubicCoefficientError) -> Self {
        Self::Coefficient(error)
    }
}

fn raise_second_and_third<D: Copy, const N: usize>(
    components: &[[[Quantity<D, f64>; N]; N]; N],
    inverse: &[[Dimensionless; N]; N],
) -> [[[Quantity<D, f64>; N]; N]; N] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (0..N)
                    .flat_map(|l| (0..N).map(move |m| (l, m)))
                    .map(|(l, m)| {
                        components[i][l][m]
                            * inverse[j][l].into_value()
                            * inverse[k][m].into_value()
                    })
                    .fold(Quantity::new(0.0), Add::add)
            })
        })
    })
}

/// Reorders a vertex stored as `[output][input one][input two]` into the
/// `[leg 1][leg 2][leg 3]` order of a three-point correlation.
fn vertex_to_leg_order<D: Copy, const N: usize>(
    output: FourierLeg,
    components: &[[[Quantity<D, f64>; N]; N]; N],
) -> [[[Quantity<D, f64>; N]; N]; N] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| match output {
                FourierLeg::One => components[i][j][k],
                FourierLeg::Two => components[j][i][k],
                FourierLeg::Three => components[k][i][j],
            })
        })
    })
}

/// Evaluates the fully raised tensors used by the subhorizon initial
/// condition on one background state and canonical momentum triangle.
///
/// The fast part of `A` is the only cubic term proportional to a physical
/// wavenumber squared,
/// `A_fast_IJK = Sym[dot(phi)_I G_JK (k_J . k_K)] / (2 H M_pl^2 a^2)`.
/// It is separated because its explicit `a^{-2}` changes the time integral.
/// The `B` and `C` tensors reuse the output-leg vertices of
/// [`canonical_cubic_hamiltonian_vertices`], so every field index still moves
/// with its Fourier leg.
pub fn canonical_subhorizon_initial_vertices<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    state: &MultifieldState<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalSubhorizonInitialVertices<N>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    let vertices = canonical_cubic_hamiltonian_vertices(system, state, triangle)?;
    let hubble = system.hubble(state);
    let planck = Mass::new(system.m_pl.into_value());
    let scale = crate::ScaleFactor::new(state.a.into_value());
    let inverse = *system.geometry.inverse_metric(&state.phi).components();
    let velocity = state.pi;
    let legs = [FourierLeg::One, FourierLeg::Two, FourierLeg::Three];

    let fast_unsymmetrized =
        |i: usize, j: usize, k: usize, second: FourierLeg, third: FourierLeg| {
            velocity[i] * triangle.wavevector_dot(second, third) * inverse[j][k].into_value()
                / (hubble * planck * planck * 2.0)
        };
    let fast: [[[Quantity<M1, f64>; N]; N]; N] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                (fast_unsymmetrized(i, j, k, legs[1], legs[2])
                    + fast_unsymmetrized(i, k, j, legs[2], legs[1])
                    + fast_unsymmetrized(j, i, k, legs[0], legs[2])
                    + fast_unsymmetrized(j, k, i, legs[2], legs[0])
                    + fast_unsymmetrized(k, i, j, legs[0], legs[1])
                    + fast_unsymmetrized(k, j, i, legs[1], legs[0]))
                    / 6.0
            })
        })
    });
    let full_a = raise_second_and_third(
        vertices
            .for_output(FourierLeg::One)
            .a_first_raised
            .components(),
        &inverse,
    );
    let slow_a = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| full_a[i][j][k] - fast[i][j][k] / (scale * scale))
        })
    });
    let b_momentum_on_leg = legs.map(|leg| {
        ThreePointBlock::new(vertex_to_leg_order(
            leg,
            &raise_second_and_third(
                vertices.for_output(leg).b_third_raised.components(),
                &inverse,
            ),
        ))
    });
    let c_field_on_leg = legs.map(|leg| {
        ThreePointBlock::new(vertex_to_leg_order(
            leg,
            &raise_second_and_third(
                vertices.for_output(leg).c_third_raised.components(),
                &inverse,
            ),
        ))
    });

    Ok(CanonicalSubhorizonInitialVertices {
        fast_a_times_scale_squared: ThreePointBlock::new(fast),
        slow_a: ThreePointBlock::new(slow_a),
        b_momentum_on_leg,
        c_field_on_leg,
    })
}

/// Assembles the tree-level Bunch--Davies three-point function on one early
/// slice from fully raised cubic tensors and their e-fold rates.
///
/// `wavenumbers` are comoving and follow the leg order of `vertices`. The
/// result keeps, separately for each block and for each of the fast-`A`,
/// slow-`A`, `B`, and `C` contributions, the highest power of `k/(aH)` that
/// survives in the in-in integral with massless de Sitter mode functions.
/// Blocks with an odd number of momentum slots vanish at the naive leading
/// order; their first surviving order is retained.
///
/// `rates` holds the covariant e-fold derivative `D_N` of every tensor in
/// `vertices`, in the same layout and dimensions. Each coefficient is expanded
/// to first order about the slice, and the rate terms at the same power of
/// `k/(aH)` as the retained static term of that contribution are kept. They
/// are a relative correction of the order of the background-variation rate.
/// Zero rates give the frozen-coefficient approximation. Every component is
/// real.
pub fn canonical_three_point_from_subhorizon_vertices<const N: usize>(
    vertices: &CanonicalSubhorizonInitialVertices<N>,
    rates: &CanonicalSubhorizonInitialVertices<N>,
    wavenumbers: [Mass; 3],
    scale_factor: crate::ScaleFactor,
    hubble: Hubble,
) -> Result<CanonicalThreePointPhaseSpaceCorrelation<N>, ThreePointInitialConditionError> {
    for (index, wavenumber) in wavenumbers.into_iter().enumerate() {
        let value = wavenumber.into_value();
        if !value.is_finite() || value <= 0.0 {
            return Err(ThreePointInitialConditionError::InvalidWavenumber {
                leg: FourierLeg::from_index(index),
                value,
            });
        }
    }
    let scale_value = scale_factor.into_value();
    if !scale_value.is_finite() || scale_value <= 0.0 {
        return Err(ThreePointInitialConditionError::InvalidScaleFactor(
            scale_value,
        ));
    }
    let hubble_value = hubble.into_value();
    if !hubble_value.is_finite() || hubble_value <= 0.0 {
        return Err(ThreePointInitialConditionError::InvalidHubble(hubble_value));
    }

    let k = wavenumbers;
    let perimeter = k[0] + k[1] + k[2];
    let product = k[0] * k[1] * k[2];
    let pair_sum = k[0] * k[1] + k[0] * k[2] + k[1] * k[2];
    let normalization = Quantity::<M0, f64>::new(0.25) / (product * perimeter);
    let a2 = scale_value * scale_value;
    let a4 = a2 * a2;
    let a6 = a4 * a2;
    let fast = vertices.fast_a_times_scale_squared.components();
    let slow = vertices.slow_a.components();
    let b = vertices.b_momentum_on_leg.map(|block| *block.components());
    let c = vertices.c_field_on_leg.map(|block| *block.components());
    let others = |leg: usize| match leg {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };

    let qqq_bracket = |i: usize, j: usize, l: usize| -> Quantity<M1, f64> {
        fast[i][j][l] * 3.0 / a4 + slow[i][j][l] * 3.0 / a2
            - hubble * (b[0][i][j][l] + b[1][i][j][l] + b[2][i][j][l]) / a2
            - (k[1] * k[2] * c[0][i][j][l]
                + k[0] * k[2] * c[1][i][j][l]
                + k[0] * k[1] * c[2][i][j][l])
                / a4
    };
    // One momentum slot on leg `p`; `u` and `v` are the two field legs.
    let one_momentum_bracket = |p: usize, i: usize, j: usize, l: usize| -> Quantity<M2, f64> {
        let (u, v) = others(p);
        let fast_weight = hubble * (perimeter + k[p]) / perimeter * -3.0;
        let slow_weight = hubble * (perimeter - k[p]) / perimeter * -3.0;
        let b_term = (0..3)
            .map(|s| k[p] * k[s] * b[s][i][j][l])
            .fold(Quantity::new(0.0), Add::add);
        let same_leg_c_weight = -(hubble
            * (k[p] * k[p] * (k[u] + k[v]) + k[p] * (k[u] * k[u] + k[v] * k[v])
                - k[u] * k[v] * (k[u] + k[v]))
            / perimeter);
        let other_leg_c_weight = |s: usize| -(hubble * k[p] * k[p] * (k[p] + k[s]) / perimeter);

        fast[i][j][l] * fast_weight / a4
            + slow[i][j][l] * slow_weight / a2
            + b_term / a4
            + (c[p][i][j][l] * same_leg_c_weight
                + c[u][i][j][l] * other_leg_c_weight(u)
                + c[v][i][j][l] * other_leg_c_weight(v))
                / a4
    };
    // Two momentum slots; `q` is the remaining field leg.
    let two_momentum_bracket = |q: usize, i: usize, j: usize, l: usize| -> Quantity<M3, f64> {
        let (u, v) = others(q);
        let momentum_pair = k[u] * k[v];
        fast[i][j][l] * momentum_pair * -3.0 / a6 + slow[i][j][l] * momentum_pair * -3.0 / a4
            - hubble
                * (k[u] * k[u] * b[u][i][j][l]
                    + k[v] * k[v] * b[v][i][j][l]
                    + (k[q] * (k[u] + k[v]) - momentum_pair) * b[q][i][j][l])
                / a4
            + (product * (k[v] * c[u][i][j][l] + k[u] * c[v][i][j][l])
                + momentum_pair * momentum_pair * c[q][i][j][l])
                / a6
    };
    let ppp_bracket = |i: usize, j: usize, l: usize| -> Quantity<M4, f64> {
        let c_term = (0..3)
            .map(|s| {
                let (u, v) = others(s);
                c[s][i][j][l] * k[u] * k[u] * k[v] * k[v] * (perimeter + k[s]) / perimeter
            })
            .fold(Quantity::new(0.0), Add::add);
        fast[i][j][l] * hubble * (pair_sum + product / perimeter) * 3.0 / a6
            + slow[i][j][l] * hubble * (pair_sum - product / perimeter) * 3.0 / a4
            - product * (k[0] * b[0][i][j][l] + k[1] * b[1][i][j][l] + k[2] * b[2][i][j][l]) / a6
            - hubble * c_term / a6
    };

    let fast_rate = rates.fast_a_times_scale_squared.components();
    let slow_rate = rates.slow_a.components();
    let b_rate = rates.b_momentum_on_leg.map(|block| *block.components());
    let c_rate = rates.c_field_on_leg.map(|block| *block.components());
    let weighted_b_rate = |i: usize, j: usize, l: usize| {
        (0..3)
            .map(|r| k[r] * b_rate[r][i][j][l])
            .fold(Quantity::new(0.0), Add::add)
    };
    let qqq_rate = |i: usize, j: usize, l: usize| -> Quantity<M1, f64> {
        -(hubble * weighted_b_rate(i, j, l) / perimeter / a2)
    };
    let one_momentum_rate = |p: usize, i: usize, j: usize, l: usize| -> Quantity<M2, f64> {
        let c_term = (0..3)
            .map(|r| product / k[r] * c_rate[r][i][j][l])
            .fold(Quantity::new(0.0), Add::add);
        hubble * k[p] / perimeter
            * (fast_rate[i][j][l] * 3.0 / a4 + slow_rate[i][j][l] * 3.0 / a2 - c_term / a4)
    };
    let two_momentum_rate = |q: usize, i: usize, j: usize, l: usize| -> Quantity<M3, f64> {
        let (u, v) = others(q);
        hubble * k[u] * k[v] / perimeter * weighted_b_rate(i, j, l) / a4
    };
    let ppp_rate = |i: usize, j: usize, l: usize| -> Quantity<M4, f64> {
        let c_term = (0..3)
            .map(|s| {
                let (u, v) = others(s);
                k[u] * k[v] * c_rate[s][i][j][l]
            })
            .fold(Quantity::new(0.0), Add::add);
        -(hubble * product / perimeter
            * (fast_rate[i][j][l] * 3.0 / a6 + slow_rate[i][j][l] * 3.0 / a4 - c_term / a6))
    };

    fn real_block<A, B, C, D: Copy, const N: usize>(
        value: impl Fn(usize, usize, usize) -> Quantity<D, f64>,
    ) -> ThreePointBlock<A, B, C, D, Complex64, N> {
        ThreePointBlock::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|l| {
                    Quantity::new(Complex64::new(value(i, j, l).into_value(), 0.0))
                })
            })
        }))
    }

    Ok(ThreePointPhaseSpaceCorrelation {
        qqq: real_block(|i, j, l| normalization * (qqq_bracket(i, j, l) + qqq_rate(i, j, l))),
        qqp: real_block(|i, j, l| {
            normalization * (one_momentum_bracket(2, i, j, l) + one_momentum_rate(2, i, j, l))
        }),
        qpq: real_block(|i, j, l| {
            normalization * (one_momentum_bracket(1, i, j, l) + one_momentum_rate(1, i, j, l))
        }),
        pqq: real_block(|i, j, l| {
            normalization * (one_momentum_bracket(0, i, j, l) + one_momentum_rate(0, i, j, l))
        }),
        qpp: real_block(|i, j, l| {
            normalization * (two_momentum_bracket(0, i, j, l) + two_momentum_rate(0, i, j, l))
        }),
        pqp: real_block(|i, j, l| {
            normalization * (two_momentum_bracket(1, i, j, l) + two_momentum_rate(1, i, j, l))
        }),
        ppq: real_block(|i, j, l| {
            normalization * (two_momentum_bracket(2, i, j, l) + two_momentum_rate(2, i, j, l))
        }),
        ppp: real_block(|i, j, l| normalization * (ppp_bracket(i, j, l) + ppp_rate(i, j, l))),
    })
}

/// Three background states around an initial slice, used to evaluate e-fold
/// rates of the cubic tensors by central differences.
///
/// `before` and `after` lie `e_fold_step` e-folds before and after `center`
/// on the same trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackgroundNeighborhood<const N: usize> {
    /// State at `N - e_fold_step`.
    pub before: MultifieldState<N>,
    /// State on the initial slice.
    pub center: MultifieldState<N>,
    /// State at `N + e_fold_step`.
    pub after: MultifieldState<N>,
    /// Positive finite spacing in elapsed e-folds.
    pub e_fold_step: f64,
}

impl BackgroundNeighborhood<2> {
    /// Samples an immutable two-field trajectory around `e_folds`.
    pub fn from_trajectory<B: TwoFieldBackgroundProvider>(
        trajectory: &B,
        e_folds: f64,
        e_fold_step: f64,
    ) -> Self {
        Self {
            before: trajectory.state_at_e_folds(e_folds - e_fold_step),
            center: trajectory.state_at_e_folds(e_folds),
            after: trajectory.state_at_e_folds(e_folds + e_fold_step),
            e_fold_step,
        }
    }
}

/// Evaluates the covariant e-fold derivative `D_N` of every subhorizon
/// vertex tensor on the central slice of `neighborhood`.
///
/// The coordinate derivative is a central difference; the connection term
/// `Gamma^I_LM (dphi^L/dN) T^M..` is added for each of the three contravariant
/// indices, so the result is the rate seen in a parallel-transported frame.
pub fn canonical_subhorizon_vertex_rates<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    neighborhood: &BackgroundNeighborhood<N>,
    triangle: MomentumTriangle,
) -> Result<CanonicalSubhorizonInitialVertices<N>, ThreePointInitialConditionError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    let step = neighborhood.e_fold_step;
    if !step.is_finite() || step <= 0.0 {
        return Err(ThreePointInitialConditionError::InvalidEFoldStep(step));
    }
    let before = canonical_subhorizon_initial_vertices(system, &neighborhood.before, triangle)?;
    let center = canonical_subhorizon_initial_vertices(system, &neighborhood.center, triangle)?;
    let after = canonical_subhorizon_initial_vertices(system, &neighborhood.after, triangle)?;
    let state = &neighborhood.center;
    let hubble = system.hubble(state);
    let connection = system.geometry.christoffel(&state.phi);
    // flow[I][M] = Gamma^I_LM dphi^L/dN, dimensionless.
    let flow: [[f64; N]; N] = core::array::from_fn(|upper| {
        core::array::from_fn(|lower| {
            (0..N)
                .map(|velocity| {
                    connection.components()[upper][velocity][lower].value()
                        * (state.pi[velocity] / hubble).into_value()
                })
                .sum()
        })
    });
    fn rate<A, B, C, D: Copy, const N: usize>(
        before: &ThreePointBlock<A, B, C, D, f64, N>,
        center: &ThreePointBlock<A, B, C, D, f64, N>,
        after: &ThreePointBlock<A, B, C, D, f64, N>,
        step: f64,
        flow: &[[f64; N]; N],
    ) -> ThreePointBlock<A, B, C, D, f64, N> {
        let t = center.components();
        ThreePointBlock::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    let connection_term = (0..N)
                        .map(|m| {
                            t[m][j][k] * flow[i][m]
                                + t[i][m][k] * flow[j][m]
                                + t[i][j][m] * flow[k][m]
                        })
                        .fold(Quantity::new(0.0), Add::add);
                    (after.components()[i][j][k] - before.components()[i][j][k]) / (2.0 * step)
                        + connection_term
                })
            })
        }))
    }

    Ok(CanonicalSubhorizonInitialVertices {
        fast_a_times_scale_squared: rate(
            &before.fast_a_times_scale_squared,
            &center.fast_a_times_scale_squared,
            &after.fast_a_times_scale_squared,
            step,
            &flow,
        ),
        slow_a: rate(&before.slow_a, &center.slow_a, &after.slow_a, step, &flow),
        b_momentum_on_leg: core::array::from_fn(|leg| {
            rate(
                &before.b_momentum_on_leg[leg],
                &center.b_momentum_on_leg[leg],
                &after.b_momentum_on_leg[leg],
                step,
                &flow,
            )
        }),
        c_field_on_leg: core::array::from_fn(|leg| {
            rate(
                &before.c_field_on_leg[leg],
                &center.c_field_on_leg[leg],
                &after.c_field_on_leg[leg],
                step,
                &flow,
            )
        }),
    })
}

/// Evaluates the canonical `(Q, Pi)` three-point initial condition on one
/// common early slice for all three legs.
///
/// The slice is the center of `neighborhood`; its outer states supply the
/// e-fold rates of the cubic tensors through
/// [`canonical_subhorizon_vertex_rates`]. The slice is accepted only if every
/// canonical leg satisfies
/// `k_r/(aH) >= minimum_k_over_ah`; the shortest leg therefore decides
/// whether the whole triangle may start. The returned blocks use the
/// triangle's canonical leg order and the stripped Fourier normalization of
/// [`ThreePointPhaseSpaceCorrelation`]. They are canonical-momentum
/// correlations: convert them with [`velocity_three_point_from_canonical_tree`]
/// only at an explicit velocity-output boundary.
pub fn canonical_three_point_initial_condition<P, G, const N: usize>(
    system: &CovariantMultifieldSystem<P, G, N>,
    neighborhood: &BackgroundNeighborhood<N>,
    triangle: MomentumTriangle,
    minimum_k_over_ah: f64,
) -> Result<CanonicalThreePointPhaseSpaceCorrelation<N>, ThreePointInitialConditionError>
where
    P: ThriceDifferentiableMultifieldPotential<N>,
    G: ThirdOrderFieldSpaceGeometry<N> + CurvatureDerivativeGeometry<N>,
{
    if !minimum_k_over_ah.is_finite() || minimum_k_over_ah <= 1.0 {
        return Err(
            ThreePointInitialConditionError::InvalidSubhorizonRequirement(minimum_k_over_ah),
        );
    }
    let state = &neighborhood.center;
    let hubble = system.hubble(state);
    let hubble_value = hubble.into_value();
    if !hubble_value.is_finite() || hubble_value <= 0.0 {
        return Err(ThreePointInitialConditionError::InvalidHubble(hubble_value));
    }
    let scale_value = state.a.into_value();
    if !scale_value.is_finite() || scale_value <= 0.0 {
        return Err(ThreePointInitialConditionError::InvalidScaleFactor(
            scale_value,
        ));
    }
    let wavenumbers = triangle.canonical_wavenumbers().map(Mass::new);
    for (index, wavenumber) in wavenumbers.into_iter().enumerate() {
        let k_over_ah = (wavenumber / (state.a * hubble)).into_value();
        if k_over_ah < minimum_k_over_ah {
            return Err(
                ThreePointInitialConditionError::InsufficientSubhorizonDepth {
                    leg: FourierLeg::from_index(index),
                    k_over_ah,
                    required: minimum_k_over_ah,
                },
            );
        }
    }
    let vertices = canonical_subhorizon_initial_vertices(system, state, triangle)?;
    let rates = canonical_subhorizon_vertex_rates(system, neighborhood, triangle)?;

    canonical_three_point_from_subhorizon_vertices(&vertices, &rates, wavenumbers, state.a, hubble)
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
        let flows = linear_phase_space_flows(self.background, &background, self.wavenumbers());
        homogeneous_three_point_derivative(state, &flows)
    }
}

fn homogeneous_three_point_derivative<S: PhaseSpaceSlot>(
    state: &ThreePointPhaseSpaceCorrelation<2, S>,
    flows: &[LinearPhaseSpaceFlow; 3],
) -> ThreePointPhaseSpaceCorrelation<2, S> {
    ThreePointPhaseSpaceCorrelation {
        qqq: homogeneous_three_point_block(state, [PhaseSpaceSlotKind::Field; 3], flows),
        qqp: homogeneous_three_point_block(
            state,
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ],
            flows,
        ),
        qpq: homogeneous_three_point_block(
            state,
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ],
            flows,
        ),
        pqq: homogeneous_three_point_block(
            state,
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Field,
            ],
            flows,
        ),
        qpp: homogeneous_three_point_block(
            state,
            [
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
            ],
            flows,
        ),
        pqp: homogeneous_three_point_block(
            state,
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
                PhaseSpaceSlotKind::Momentum,
            ],
            flows,
        ),
        ppq: homogeneous_three_point_block(
            state,
            [
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Momentum,
                PhaseSpaceSlotKind::Field,
            ],
            flows,
        ),
        ppp: homogeneous_three_point_block(state, [PhaseSpaceSlotKind::Momentum; 3], flows),
    }
}

/// The complete tree-level canonical state evolved by
/// [`CanonicalThreePointTransportSystem`].
///
/// The source `u Sigma Sigma` needs the linear covariance of every leg at the
/// same time as `Alpha`, so the three covariances are integrated alongside it
/// on the same background. Entry `r` of `covariances` belongs to canonical leg
/// `r` of the triangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CanonicalThreePointTransportState {
    /// Canonical two-point covariances of legs one, two, and three.
    pub covariances: [CanonicalTwoPointPhaseSpaceCovariance<2>; 3],
    /// Connected canonical three-point function.
    pub alpha: CanonicalThreePointPhaseSpaceCorrelation<2>,
}

impl Add for CanonicalThreePointTransportState {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            covariances: [
                self.covariances[0] + rhs.covariances[0],
                self.covariances[1] + rhs.covariances[1],
                self.covariances[2] + rhs.covariances[2],
            ],
            alpha: self.alpha + rhs.alpha,
        }
    }
}

impl Mul<f64> for CanonicalThreePointTransportState {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            covariances: self.covariances.map(|covariance| covariance * rhs),
            alpha: self.alpha * rhs,
        }
    }
}

impl Mul<Dimensionless> for CanonicalThreePointTransportState {
    type Output = Self;

    fn mul(self, rhs: Dimensionless) -> Self::Output {
        self * rhs.into_value()
    }
}

/// Absolute error floors for the three-point blocks, each in its own
/// physical dimension.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThreePointTransportTolerance {
    /// Absolute tolerance for `QQQ`.
    pub qqq: Quantity<MMinus3, f64>,
    /// Absolute tolerance for each block with one momentum slot.
    pub one_momentum: Quantity<MMinus2, f64>,
    /// Absolute tolerance for each block with two momentum slots.
    pub two_momenta: Quantity<MMinus1, f64>,
    /// Absolute tolerance for `PPP`.
    pub ppp: Quantity<M0, f64>,
}

/// Absolute error floors for the combined canonical three-point system.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CanonicalThreePointTransportTolerance {
    /// Floors applied to each leg covariance.
    pub covariance: TwoPointTransportTolerance,
    /// Floors applied to `Alpha`.
    pub alpha: ThreePointTransportTolerance,
}

impl CanonicalThreePointTransportTolerance {
    /// Uses pure block-relative control everywhere.
    pub const fn pure_relative() -> Self {
        Self {
            covariance: TwoPointTransportTolerance::pure_relative(),
            alpha: ThreePointTransportTolerance {
                qqq: Quantity::new(0.0),
                one_momentum: Quantity::new(0.0),
                two_momenta: Quantity::new(0.0),
                ppp: Quantity::new(0.0),
            },
        }
    }
}

/// Tree-level canonical three-point transport on an immutable background.
///
/// The derivative is
/// `D_N Alpha = u(k_1) Alpha + u(k_2) Alpha + u(k_3) Alpha + S`, where `S` is
/// the operator-ordered `u^a_de Sigma Sigma` source built at every stage from
/// the model-derived cubic vertices on the interpolated background. The three
/// covariances evolve with the linear flow on the same background. Every
/// correlation uses [`CanonicalMomentumPerturbation`]; convert with
/// [`velocity_three_point_from_transport_state`] only at output.
///
/// # Panics
///
/// The derivative panics if the background supplies a non-positive or
/// non-finite Hubble rate, scale factor, or reduced Planck mass, because the
/// cubic vertices are undefined there.
pub struct CanonicalThreePointTransportSystem<'a, P, G, B> {
    background: &'a CovariantMultifieldSystem<P, G, 2>,
    trajectory: &'a B,
    triangle: MomentumTriangle,
}

impl<'a, P, G, B> CanonicalThreePointTransportSystem<'a, P, G, B> {
    /// Builds a system whose legs use the triangle's canonical order.
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

    /// Returns the canonically ordered comoving wavenumbers.
    pub fn wavenumbers(&self) -> [Mass; 3] {
        self.triangle.canonical_wavenumbers().map(Mass::new)
    }
}

impl<P, G, B> System for CanonicalThreePointTransportSystem<'_, P, G, B>
where
    P: ThriceDifferentiableMultifieldPotential<2>,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2>,
    B: TwoFieldBackgroundProvider,
{
    type IndependentVariable = Dimensionless;
    type State = CanonicalThreePointTransportState;
    type Delta = CanonicalThreePointTransportState;
    type Derivative = CanonicalThreePointTransportState;

    fn derivative(&self, e_folds: Dimensionless, state: &Self::State) -> Self::Derivative {
        let background = self.trajectory.state_at_e_folds(e_folds.into_value());
        let wavenumbers = self.wavenumbers();
        let flows = linear_phase_space_flows(self.background, &background, wavenumbers);
        let kernels =
            canonical_cubic_hamiltonian_vertices(self.background, &background, self.triangle)
                .expect("the background must admit cubic vertices")
                .e_fold_kernels(self.background.hubble(&background))
                .expect("the background must have a positive Hubble rate");
        let source = three_point_source_from_kernels(
            [&kernels[0], &kernels[1], &kernels[2]],
            [
                &state.covariances[0],
                &state.covariances[1],
                &state.covariances[2],
            ],
        );

        CanonicalThreePointTransportState {
            covariances: core::array::from_fn(|leg| {
                two_point_linear_derivative(
                    self.background,
                    &background,
                    wavenumbers[leg],
                    &state.covariances[leg],
                )
            }),
            alpha: homogeneous_three_point_derivative(&state.alpha, &flows) + source,
        }
    }
}

impl<P, G, B> AdaptiveSystem for CanonicalThreePointTransportSystem<'_, P, G, B>
where
    P: ThriceDifferentiableMultifieldPotential<2>,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2>,
    B: TwoFieldBackgroundProvider,
{
    type AbsoluteTolerance = CanonicalThreePointTransportTolerance;

    fn error_ratio(
        &self,
        current: &Self::State,
        next: &Self::State,
        error: &Self::Delta,
        tolerances: &ErrorTolerances<Self::AbsoluteTolerance>,
    ) -> f64 {
        let relative = tolerances.relative;
        let floors = tolerances.absolute;
        let covariance_ratio = (0..3)
            .map(|leg| {
                let (current, next, error) = (
                    &current.covariances[leg],
                    &next.covariances[leg],
                    &error.covariances[leg],
                );
                covariance_block_error_ratio(
                    &current.qq,
                    &next.qq,
                    &error.qq,
                    floors.covariance.qq,
                    relative,
                )
                .max(covariance_block_error_ratio(
                    &current.qp,
                    &next.qp,
                    &error.qp,
                    floors.covariance.qp,
                    relative,
                ))
                .max(covariance_block_error_ratio(
                    &current.pq,
                    &next.pq,
                    &error.pq,
                    floors.covariance.pq,
                    relative,
                ))
                .max(covariance_block_error_ratio(
                    &current.pp,
                    &next.pp,
                    &error.pp,
                    floors.covariance.pp,
                    relative,
                ))
            })
            .fold(0.0, f64::max);
        let (current, next, error) = (&current.alpha, &next.alpha, &error.alpha);
        let alpha = floors.alpha;

        [
            three_point_block_error_ratio(&current.qqq, &next.qqq, &error.qqq, alpha.qqq, relative),
            three_point_block_error_ratio(
                &current.qqp,
                &next.qqp,
                &error.qqp,
                alpha.one_momentum,
                relative,
            ),
            three_point_block_error_ratio(
                &current.qpq,
                &next.qpq,
                &error.qpq,
                alpha.one_momentum,
                relative,
            ),
            three_point_block_error_ratio(
                &current.pqq,
                &next.pqq,
                &error.pqq,
                alpha.one_momentum,
                relative,
            ),
            three_point_block_error_ratio(
                &current.qpp,
                &next.qpp,
                &error.qpp,
                alpha.two_momenta,
                relative,
            ),
            three_point_block_error_ratio(
                &current.pqp,
                &next.pqp,
                &error.pqp,
                alpha.two_momenta,
                relative,
            ),
            three_point_block_error_ratio(
                &current.ppq,
                &next.ppq,
                &error.ppq,
                alpha.two_momenta,
                relative,
            ),
            three_point_block_error_ratio(&current.ppp, &next.ppp, &error.ppp, alpha.ppp, relative),
        ]
        .into_iter()
        .fold(covariance_ratio, f64::max)
    }
}

fn three_point_block_error_ratio<A, B, C, D: Copy>(
    current: &ThreePointBlock<A, B, C, D, Complex64, 2>,
    next: &ThreePointBlock<A, B, C, D, Complex64, 2>,
    error: &ThreePointBlock<A, B, C, D, Complex64, 2>,
    absolute: Quantity<D, f64>,
    relative: f64,
) -> f64 {
    let components = |block: &ThreePointBlock<A, B, C, D, Complex64, 2>| {
        block
            .components
            .into_iter()
            .flatten()
            .flatten()
            .map(|component| component.value().norm())
            .collect::<Vec<_>>()
    };
    let reference = components(current)
        .into_iter()
        .chain(components(next))
        .fold(0.0, f64::max);
    let scale = absolute.into_value() + relative * reference;
    let errors = components(error);

    if scale == 0.0 {
        return if errors.iter().all(|error| *error == 0.0) {
            0.0
        } else {
            f64::INFINITY
        };
    }

    errors
        .into_iter()
        .map(|error| error / scale)
        .fold(0.0, f64::max)
}

/// Builds the complete initial transport state on one common subhorizon slice.
///
/// Each leg covariance is the Bunch--Davies mode state of the existing linear
/// solver, relabelled to the canonical basis at linear order, and `Alpha` is
/// [`canonical_three_point_initial_condition`]. Every canonical leg must
/// satisfy `k_r/(aH) >= minimum_k_over_ah`.
pub fn canonical_three_point_transport_initial_state<P, G>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    neighborhood: &BackgroundNeighborhood<2>,
    triangle: MomentumTriangle,
    minimum_k_over_ah: f64,
) -> Result<CanonicalThreePointTransportState, ThreePointInitialConditionError>
where
    P: ThriceDifferentiableMultifieldPotential<2> + Copy,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2> + Copy,
{
    let alpha =
        canonical_three_point_initial_condition(system, neighborhood, triangle, minimum_k_over_ah)?;
    let background = &neighborhood.center;
    let wavenumbers = triangle.canonical_wavenumbers().map(Mass::new);
    let covariance = |leg: usize| {
        crate::TwoFieldPerturbationSystem::new(*system, wavenumbers[leg])
            .bunch_davies_state(*background, minimum_k_over_ah)
            .map(|modes| {
                canonical_two_point_from_velocity_linear(two_point_covariance_from_mode_state(
                    system, &modes,
                ))
            })
            .map_err(ThreePointInitialConditionError::ModeInitialization)
    };

    Ok(CanonicalThreePointTransportState {
        covariances: [covariance(0)?, covariance(1)?, covariance(2)?],
        alpha,
    })
}

/// Converts an evolved canonical state to the velocity basis at tree level.
///
/// The vertices are evaluated on `background`, which must be the slice on
/// which `state` was evaluated.
pub fn velocity_three_point_from_transport_state<P, G>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    background: &MultifieldState<2>,
    triangle: MomentumTriangle,
    state: &CanonicalThreePointTransportState,
) -> Result<ThreePointPhaseSpaceCorrelation<2>, CubicCoefficientError>
where
    P: ThriceDifferentiableMultifieldPotential<2>,
    G: ThirdOrderFieldSpaceGeometry<2> + CurvatureDerivativeGeometry<2>,
{
    let vertices = canonical_cubic_hamiltonian_vertices(system, background, triangle)?;
    Ok(velocity_three_point_from_canonical_tree(
        state.alpha,
        &vertices,
        [
            &state.covariances[0],
            &state.covariances[1],
            &state.covariances[2],
        ],
    ))
}

/// Linear phase-space flow of each Fourier leg in elapsed e-folds.
fn linear_phase_space_flows<P, G>(
    system: &CovariantMultifieldSystem<P, G, 2>,
    background: &MultifieldState<2>,
    wavenumbers: [Mass; 3],
) -> [LinearPhaseSpaceFlow; 3]
where
    P: TwiceDifferentiableMultifieldPotential<2>,
    G: CurvedFieldSpaceGeometry<2>,
{
    let hubble = system.hubble(background).into_value();
    let connection = system.geometry.christoffel(&background.phi);
    let mass = system.perturbation_mass_matrix(background);
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
    wavenumbers.map(|wavenumber| {
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
    })
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

fn homogeneous_three_point_block<A, B, C, D, S: PhaseSpaceSlot>(
    state: &ThreePointPhaseSpaceCorrelation<2, S>,
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
        EnergyDensity, Field, FieldVelocity, FlatFieldSpace, HyperbolicFieldSpace, M0, Mass,
        MultifieldPotentialValue, MultifieldState, QuadraticMultifieldPotential, ScaleFactor,
        TwoFieldPerturbationSystem,
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
    fn source_keeps_left_external_operators_on_the_left_of_the_insertion() {
        let first_leg = TwoPointPhaseSpaceCovariance::<1> {
            qq: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            qp: TwoPointBlock::new([[Quantity::new(Complex64::new(2.0, 3.0))]]),
            pq: TwoPointBlock::new([[Quantity::new(Complex64::new(2.0, -3.0))]]),
            pp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
        };
        let third_leg = TwoPointPhaseSpaceCovariance::<1> {
            qq: TwoPointBlock::new([[Quantity::new(Complex64::new(5.0, 0.0))]]),
            qp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            pq: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            pp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
        };
        let zero_kernel = QuadraticPhaseSpaceKernel::<1>::zero();
        let middle_kernel = QuadraticPhaseSpaceKernel::<1> {
            q_pq: ThreePointBlock::new([[[Quantity::new(1.0)]]]),
            ..QuadraticPhaseSpaceKernel::zero()
        };

        let source = three_point_source_from_kernels(
            [&zero_kernel, &middle_kernel, &zero_kernel],
            [&first_leg, &third_leg, &third_leg],
        );

        // <Q_1 (u Pi_1 Q_3) Q_3> contracts as Sigma^{QP}(k_1) Sigma^{QQ}(k_3).
        assert_eq!(
            source.qqq.components()[0][0][0].into_value(),
            Complex64::new(10.0, 15.0)
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
            b_first_raised_input_exchanged: block::<M0>(6_000.0),
            b_third_raised: block::<M0>(3_000.0),
            c_first_raised: block::<MMinus1>(4_000.0),
            c_first_raised_input_exchanged: block::<MMinus1>(7_000.0),
            c_third_raised: block::<MMinus1>(5_000.0),
        };
        let kernel = CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
            coefficients,
            Hubble::new(2.0),
        )
        .unwrap();

        assert_eq!(*kernel.q_qq.components()[1][0][1].value(), -3_101.0 / 2.0);
        assert_eq!(*kernel.q_qp.components()[1][0][1].value(), -7_101.0 / 2.0);
        assert_eq!(*kernel.q_pq.components()[1][0][1].value(), -4_101.0 / 2.0);
        assert_eq!(*kernel.q_pp.components()[1][0][1].value(), 0.0);
        assert_eq!(
            *kernel.p_qq.components()[1][0][1].value(),
            3.0 * 1_101.0 / 2.0
        );
        assert_eq!(*kernel.p_qp.components()[1][0][1].value(), 2_101.0 / 2.0);
        assert_eq!(*kernel.p_pq.components()[1][0][1].value(), 6_101.0 / 2.0);
        assert_eq!(*kernel.p_pp.components()[1][0][1].value(), 5_101.0 / 2.0);
    }

    #[test]
    fn cubic_hamiltonian_map_rejects_nonphysical_hubble_rates() {
        let coefficients = CanonicalCubicHamiltonianCoefficients::<1> {
            a_first_raised: zero_real_three_point_block(),
            b_first_raised: zero_real_three_point_block(),
            b_first_raised_input_exchanged: zero_real_three_point_block(),
            b_third_raised: zero_real_three_point_block(),
            c_first_raised: zero_real_three_point_block(),
            c_first_raised_input_exchanged: zero_real_three_point_block(),
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
    fn linear_two_point_basis_conversion_is_an_explicit_round_trip() {
        let velocity = TwoPointPhaseSpaceCovariance::<1> {
            qq: TwoPointBlock::new([[Quantity::new(Complex64::new(1.0, 2.0))]]),
            qp: TwoPointBlock::new([[Quantity::new(Complex64::new(3.0, 4.0))]]),
            pq: TwoPointBlock::new([[Quantity::new(Complex64::new(5.0, 6.0))]]),
            pp: TwoPointBlock::new([[Quantity::new(Complex64::new(7.0, 8.0))]]),
        };

        let canonical = canonical_two_point_from_velocity_linear(velocity);
        assert_eq!(canonical.qq, velocity.qq);
        assert_eq!(canonical.qp.components, velocity.qp.components);
        assert_eq!(canonical.pq.components, velocity.pq.components);
        assert_eq!(canonical.pp.components, velocity.pp.components);
        assert_eq!(
            velocity_two_point_from_canonical_linear(canonical),
            velocity
        );
    }

    #[test]
    fn ordered_vertices_keep_the_output_leg_attached_to_its_momentum() {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0), Mass::new(2.0)]),
            FlatFieldSpace,
        );
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(3.0)],
            [FieldVelocity::new(0.5), FieldVelocity::new(-0.25)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(3.0, 2.0, 2.0).unwrap();

        let canonical = canonical_cubic_hamiltonian_coefficients(&model, &state, triangle).unwrap();
        let vertices = canonical_cubic_hamiltonian_vertices(&model, &state, triangle).unwrap();

        assert_eq!(*vertices.for_output(FourierLeg::One), canonical);
        assert_ne!(
            vertices.for_output(FourierLeg::One),
            vertices.for_output(FourierLeg::Two)
        );
    }

    #[test]
    fn canonical_to_velocity_three_point_adds_the_b_and_c_field_redefinition() {
        fn coefficient(b: f64, c: f64) -> CanonicalCubicHamiltonianCoefficients<1> {
            CanonicalCubicHamiltonianCoefficients {
                a_first_raised: zero_real_three_point_block(),
                b_first_raised: zero_real_three_point_block(),
                b_first_raised_input_exchanged: zero_real_three_point_block(),
                b_third_raised: ThreePointBlock::new([[[Quantity::new(b)]]]),
                c_first_raised: ThreePointBlock::new([[[Quantity::new(c)]]]),
                c_first_raised_input_exchanged: ThreePointBlock::new([[[Quantity::new(10.0 * c)]]]),
                c_third_raised: zero_real_three_point_block(),
            }
        }
        fn covariance(qq: f64, pq: f64) -> CanonicalTwoPointPhaseSpaceCovariance<1> {
            CanonicalTwoPointPhaseSpaceCovariance {
                qq: TwoPointBlock::new([[Quantity::new(Complex64::new(qq, 0.0))]]),
                qp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
                pq: TwoPointBlock::new([[Quantity::new(Complex64::new(pq, 0.0))]]),
                pp: TwoPointBlock::new([[Quantity::new(Complex64::new(0.0, 0.0))]]),
            }
        }

        let vertices = CanonicalCubicHamiltonianVertices {
            coefficients: [
                coefficient(2.0, 3.0),
                coefficient(5.0, 7.0),
                coefficient(11.0, 13.0),
            ],
        };
        let covariances = [
            covariance(17.0, 19.0),
            covariance(23.0, 29.0),
            covariance(31.0, 37.0),
        ];
        let mut canonical = CanonicalThreePointPhaseSpaceCorrelation::<1>::zero();
        canonical.pqq.components[0][0][0] = Quantity::new(Complex64::new(41.0, 0.0));

        let velocity = velocity_three_point_from_canonical_tree(
            canonical,
            &vertices,
            [&covariances[0], &covariances[1], &covariances[2]],
        );
        let expected_pqq = 41.0 - 2.0 * 23.0 * 31.0 - (30.0 * 23.0 * 37.0 + 3.0 * 29.0 * 31.0);

        assert_eq!(
            velocity.pqq.components[0][0][0].into_value(),
            Complex64::new(expected_pqq, 0.0)
        );
        assert_eq!(
            velocity.qqq.components[0][0][0].into_value(),
            Complex64::new(0.0, 0.0)
        );
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
        let rotated = canonical_cubic_c_coefficients_for_legs(
            &model,
            &state,
            triangle,
            [FourierLeg::Two, FourierLeg::Three, FourierLeg::One],
        )
        .unwrap();

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let actual = coefficients.third_raised.components()[i][j][k].into_value();
                    let expected = rotated.first_raised.components()[j][k][i].into_value();
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
        assert_eq!(
            canonical_cubic_a_coefficients(&invalid_planck, &state, triangle),
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
        assert_eq!(
            canonical_cubic_a_coefficients(&model, &zero_energy, triangle),
            Err(CubicCoefficientError::InvalidHubble(0.0))
        );

        let zero_scale = MultifieldState::new(
            [Field::new(1.0)],
            [FieldVelocity::new(0.0)],
            ScaleFactor::new(0.0),
        );
        assert_eq!(
            canonical_cubic_a_coefficients(&model, &zero_scale, triangle),
            Err(CubicCoefficientError::InvalidScaleFactor(0.0))
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
        let rotated = canonical_cubic_b_coefficients_for_legs(
            &model,
            &state,
            triangle,
            [FourierLeg::Two, FourierLeg::Three, FourierLeg::One],
        )
        .unwrap();

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let actual = coefficients.third_raised.components()[i][j][k].into_value();
                    let expected = rotated.first_raised.components()[j][k][i].into_value();
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

    #[test]
    fn canonical_cubic_a_matches_the_flat_single_field_limit() {
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

        let coefficients = canonical_cubic_a_coefficients(&model, &state, triangle).unwrap();
        let hubble = model.hubble(&state).into_value();
        let velocity: f64 = 0.5;
        let gradient = 2.0;
        let hessian = 1.0;
        let z = -3.0 * hubble * velocity - gradient + velocity.powi(3) / (2.0 * hubble);
        let expected = -velocity * hessian / (2.0 * hubble)
            + velocity.powi(2) * z / (4.0 * hubble.powi(2))
            + 3.0 * velocity * z.powi(2) / (32.0 * hubble.powi(3))
            + velocity.powi(3) * (6.0 * hubble.powi(2) - velocity.powi(2)) / (8.0 * hubble.powi(3))
            - velocity / (4.0 * hubble);
        let actual = coefficients.first_raised.components()[0][0][0].into_value();
        let scale = actual.abs().max(expected.abs()).max(1.0);

        assert!((actual - expected).abs() / scale < 1.0e-14);
    }

    #[test]
    fn canonical_cubic_a_contains_the_symmetrized_covariant_third_potential_derivative() {
        struct MixedCubicPotential;

        impl<V> MultifieldPotentialValue<V, 3> for MixedCubicPotential
        where
            V: Clone + Mul<V, Output = V>,
        {
            fn value(&self, phi: [Field<V>; 3]) -> EnergyDensity<V> {
                let [x, y, z] = phi.map(Field::into_value);
                EnergyDensity::new(x * y * z)
            }
        }

        let model = CovariantMultifieldSystem::in_planck_units(MixedCubicPotential, FlatFieldSpace);
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(3.0), Field::new(5.0)],
            [FieldVelocity::new(0.0); 3],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(3.0, 2.0, 2.0).unwrap();

        let coefficients = canonical_cubic_a_coefficients(&model, &state, triangle).unwrap();

        for permutation in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            assert_eq!(
                coefficients.first_raised.components()[permutation[0]][permutation[1]]
                    [permutation[2]]
                    .into_value(),
                -1.0 / 3.0
            );
        }
        assert_eq!(
            coefficients.first_raised.components()[0][0][0].into_value(),
            0.0
        );
    }

    #[test]
    fn canonical_cubic_a_contains_the_covariant_field_space_curvature_term() {
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

        impl ThirdOrderFieldSpaceGeometry<2> for WithoutCurvature {
            fn christoffel_derivative(
                &self,
                position: &[Field; 2],
            ) -> crate::ChristoffelDerivatives<2> {
                self.0.christoffel_derivative(position)
            }
        }

        impl CurvatureDerivativeGeometry<2> for WithoutCurvature {
            fn covariant_riemann_derivative(
                &self,
                _position: &[Field; 2],
            ) -> crate::CovariantRiemannDerivative<2> {
                crate::CovariantRiemannDerivative::zero()
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

        let actual = canonical_cubic_a_coefficients(&curved, &state, triangle).unwrap();
        let without = canonical_cubic_a_coefficients(&suppressed, &state, triangle).unwrap();
        let hubble = curved.hubble(&state).into_value();
        let metric = geometry.metric(&state.phi);
        let inverse = geometry.inverse_metric(&state.phi);
        let riemann = metric.lower_first(&geometry.riemann(&state.phi));
        let lowered_velocity = metric.lower(&ContravariantField::new(state.pi));
        let raw_curvature = |i: usize, j: usize, k: usize| {
            (0..2)
                .flat_map(|l| (0..2).map(move |m| (l, m)))
                .map(|(l, m)| {
                    -lowered_velocity.components()[k].into_value()
                        * state.pi[l].into_value()
                        * state.pi[m].into_value()
                        * (riemann.components()[l][i][j][m].into_value()
                            + riemann.components()[l][j][i][m].into_value())
                        / (4.0 * hubble)
                })
                .sum::<f64>()
        };

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let expected = (0..2)
                        .map(|contracted| {
                            let lower = (raw_curvature(contracted, j, k)
                                + raw_curvature(contracted, k, j)
                                + raw_curvature(j, contracted, k)
                                + raw_curvature(j, k, contracted)
                                + raw_curvature(k, contracted, j)
                                + raw_curvature(k, j, contracted))
                                / 6.0;
                            inverse.components()[i][contracted].into_value() * lower
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

    #[test]
    fn canonical_cubic_a_routes_the_rank_five_curvature_derivative_indices() {
        #[derive(Clone, Copy)]
        struct CurvatureDerivativeProbe;

        impl FieldSpaceGeometry<2> for CurvatureDerivativeProbe {
            fn metric(&self, _position: &[Field; 2]) -> crate::FieldMetric<2> {
                crate::FieldMetric::euclidean()
            }

            fn inverse_metric(&self, _position: &[Field; 2]) -> crate::InverseFieldMetric<2> {
                crate::InverseFieldMetric::euclidean()
            }

            fn christoffel(&self, _position: &[Field; 2]) -> crate::ChristoffelSymbols<2> {
                crate::ChristoffelSymbols::zero()
            }
        }

        impl CurvedFieldSpaceGeometry<2> for CurvatureDerivativeProbe {
            fn riemann(&self, _position: &[Field; 2]) -> crate::RiemannTensor<2> {
                crate::RiemannTensor::zero()
            }
        }

        impl ThirdOrderFieldSpaceGeometry<2> for CurvatureDerivativeProbe {
            fn christoffel_derivative(
                &self,
                _position: &[Field; 2],
            ) -> crate::ChristoffelDerivatives<2> {
                crate::ChristoffelDerivatives::zero()
            }
        }

        impl CurvatureDerivativeGeometry<2> for CurvatureDerivativeProbe {
            fn covariant_riemann_derivative(
                &self,
                _position: &[Field; 2],
            ) -> crate::CovariantRiemannDerivative<2> {
                crate::CovariantRiemannDerivative::new(core::array::from_fn(|upper| {
                    core::array::from_fn(|second| {
                        core::array::from_fn(|third| {
                            core::array::from_fn(|fourth| {
                                core::array::from_fn(|derivative| {
                                    Quantity::new(
                                        if [upper, second, third, fourth, derivative]
                                            == [0, 0, 1, 1, 0]
                                        {
                                            0.75
                                        } else {
                                            0.0
                                        },
                                    )
                                })
                            })
                        })
                    })
                }))
            }
        }

        let potential = QuadraticMultifieldPotential::new([Mass::new(1.0), Mass::new(0.75)]);
        let probe = CovariantMultifieldSystem::in_planck_units(potential, CurvatureDerivativeProbe);
        let baseline = CovariantMultifieldSystem::in_planck_units(potential, FlatFieldSpace);
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(0.5)],
            [FieldVelocity::new(0.4), FieldVelocity::new(-0.3)],
            ScaleFactor::new(1.0),
        );
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();
        let actual = canonical_cubic_a_coefficients(&probe, &state, triangle).unwrap();
        let without = canonical_cubic_a_coefficients(&baseline, &state, triangle).unwrap();
        let derivative = CurvatureDerivativeProbe.covariant_riemann_derivative(&state.phi);
        let raw = |i: usize, j: usize, k: usize| {
            (0..2)
                .flat_map(|l| (0..2).map(move |m| (l, m)))
                .map(|(l, m)| {
                    derivative.components()[i][l][m][j][k].into_value()
                        * state.pi[l].into_value()
                        * state.pi[m].into_value()
                        / 3.0
                })
                .sum::<f64>()
        };

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let expected = (raw(i, j, k)
                        + raw(i, k, j)
                        + raw(j, i, k)
                        + raw(j, k, i)
                        + raw(k, i, j)
                        + raw(k, j, i))
                        / 6.0;
                    let difference = actual.first_raised.components()[i][j][k].into_value()
                        - without.first_raised.components()[i][j][k].into_value();
                    let scale = difference.abs().max(expected.abs()).max(1.0);
                    assert!((difference - expected).abs() / scale < 1.0e-14);
                }
            }
        }
    }

    #[test]
    fn model_derived_cubic_coefficients_feed_the_e_fold_kernel() {
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
        let coefficients =
            canonical_cubic_hamiltonian_coefficients(&model, &state, triangle).unwrap();

        let kernel = CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
            coefficients,
            model.hubble(&state),
        )
        .unwrap();

        assert!(kernel.q_qq.components()[0][0][0].into_value().is_finite());
        assert!(kernel.p_qq.components()[0][0][0].into_value().is_finite());
        assert!(kernel.p_pp.components()[0][0][0].into_value().is_finite());
    }

    const SCALAR_FAST_A: f64 = 1.3;
    const SCALAR_SLOW_A: f64 = -0.7;
    const SCALAR_B: [f64; 3] = [0.11, -0.23, 0.37];
    const SCALAR_C: [f64; 3] = [1.9, -2.9, 0.53];

    fn scalar_block<D>(
        value: f64,
    ) -> ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, D, f64, 1> {
        ThreePointBlock::new([[[Quantity::new(value)]]])
    }

    fn scalar_initial_vertices() -> CanonicalSubhorizonInitialVertices<1> {
        CanonicalSubhorizonInitialVertices {
            fast_a_times_scale_squared: scalar_block(SCALAR_FAST_A),
            slow_a: scalar_block(SCALAR_SLOW_A),
            b_momentum_on_leg: SCALAR_B.map(scalar_block),
            c_field_on_leg: SCALAR_C.map(scalar_block),
        }
    }

    /// Relative e-fold rates of the scalar coefficients, of slow-roll size.
    const SCALAR_RATE_FAST_A: f64 = 0.05;
    const SCALAR_RATE_SLOW_A: f64 = -0.03;
    const SCALAR_RATE_B: [f64; 3] = [0.02, -0.04, 0.01];
    const SCALAR_RATE_C: [f64; 3] = [0.06, -0.02, 0.03];

    fn scalar_rates(drift: f64) -> CanonicalSubhorizonInitialVertices<1> {
        CanonicalSubhorizonInitialVertices {
            fast_a_times_scale_squared: scalar_block(drift * SCALAR_RATE_FAST_A),
            slow_a: scalar_block(drift * SCALAR_RATE_SLOW_A),
            b_momentum_on_leg: SCALAR_RATE_B.map(|rate| scalar_block(drift * rate)),
            c_field_on_leg: SCALAR_RATE_C.map(|rate| scalar_block(drift * rate)),
        }
    }

    /// Scalar coefficients drifting linearly by `drift` times their rates over
    /// `e_folds` e-folds from the reference slice.
    fn scalar_initial_vertices_at(
        e_folds: f64,
        drift: f64,
    ) -> CanonicalSubhorizonInitialVertices<1> {
        fn drifted<D>(
            value: f64,
            rate: f64,
            drift: f64,
            e_folds: f64,
        ) -> ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, D, f64, 1>
        {
            scalar_block(value + drift * rate * e_folds)
        }
        CanonicalSubhorizonInitialVertices {
            fast_a_times_scale_squared: drifted(SCALAR_FAST_A, SCALAR_RATE_FAST_A, drift, e_folds),
            slow_a: drifted(SCALAR_SLOW_A, SCALAR_RATE_SLOW_A, drift, e_folds),
            b_momentum_on_leg: [0, 1, 2]
                .map(|r| drifted(SCALAR_B[r], SCALAR_RATE_B[r], drift, e_folds)),
            c_field_on_leg: [0, 1, 2]
                .map(|r| drifted(SCALAR_C[r], SCALAR_RATE_C[r], drift, e_folds)),
        }
    }

    /// Scalar vertex for one output leg. Its inputs are the remaining legs in
    /// ascending order, and every tensor slot keeps the momentum of its leg.
    fn scalar_output_vertex(output: usize, scale: f64) -> CanonicalCubicHamiltonianCoefficients<1> {
        let [first_input, second_input] = match output {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        CanonicalCubicHamiltonianCoefficients {
            a_first_raised: scalar_block(SCALAR_FAST_A / (scale * scale) + SCALAR_SLOW_A),
            b_first_raised: scalar_block(SCALAR_B[second_input]),
            b_first_raised_input_exchanged: scalar_block(SCALAR_B[first_input]),
            b_third_raised: scalar_block(SCALAR_B[output]),
            c_first_raised: scalar_block(SCALAR_C[second_input]),
            c_first_raised_input_exchanged: scalar_block(SCALAR_C[first_input]),
            c_third_raised: scalar_block(SCALAR_C[output]),
        }
    }

    fn massless_de_sitter_covariance(
        wavenumber: f64,
        scale: f64,
        hubble: f64,
    ) -> CanonicalTwoPointPhaseSpaceCovariance<1> {
        let tau = -1.0 / (scale * hubble);
        let normalization = hubble / (2.0 * wavenumber.powi(3)).sqrt();
        let phase = Complex64::from_polar(1.0, -wavenumber * tau);
        let q = normalization * Complex64::new(1.0, wavenumber * tau) * phase;
        let p = -normalization * hubble * wavenumber.powi(2) * tau.powi(2) * phase;
        CanonicalTwoPointPhaseSpaceCovariance {
            qq: TwoPointBlock::new([[Quantity::new(q * q.conj())]]),
            qp: TwoPointBlock::new([[Quantity::new(q * p.conj())]]),
            pq: TwoPointBlock::new([[Quantity::new(p * q.conj())]]),
            pp: TwoPointBlock::new([[Quantity::new(p * p.conj())]]),
        }
    }

    const ALL_BLOCK_SLOTS: [[PhaseSpaceSlotKind; 3]; 8] = {
        use PhaseSpaceSlotKind::{Field as F, Momentum as P};
        [
            [F, F, F],
            [F, F, P],
            [F, P, F],
            [P, F, F],
            [F, P, P],
            [P, F, P],
            [P, P, F],
            [P, P, P],
        ]
    };

    /// Largest block-relative residual of the canonical transport equation
    /// for the scalar initial condition on massless de Sitter.
    ///
    /// The coefficients drift linearly in `N` by `drift` times the scalar
    /// rates; `rates_in_initial_condition` selects whether the initial
    /// condition is told about that drift.
    fn de_sitter_transport_residual(
        depth: f64,
        drift: f64,
        rates_in_initial_condition: bool,
    ) -> f64 {
        let (scale, hubble): (f64, f64) = (1.3, 0.7);
        let shape = [1.0, 0.8, 0.5];
        let wavenumbers = shape.map(|ratio| ratio * depth * scale * hubble);
        let rates = scalar_rates(if rates_in_initial_condition {
            drift
        } else {
            0.0
        });
        let initial = |a: f64| {
            canonical_three_point_from_subhorizon_vertices(
                &scalar_initial_vertices_at((a / scale).ln(), drift),
                &rates,
                wavenumbers.map(Mass::new),
                ScaleFactor::new(a),
                Hubble::new(hubble),
            )
            .unwrap()
        };
        let step: f64 = 1.0e-5;
        let alpha = initial(scale);
        let forward = initial(scale * step.exp());
        let backward = initial(scale * (-step).exp());
        let covariances =
            wavenumbers.map(|wavenumber| massless_de_sitter_covariance(wavenumber, scale, hubble));
        let kernels = [0, 1, 2].map(|output| {
            CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian(
                scalar_output_vertex(output, scale),
                Hubble::new(hubble),
            )
            .unwrap()
        });
        let source = three_point_source_from_kernels(
            [&kernels[0], &kernels[1], &kernels[2]],
            [&covariances[0], &covariances[1], &covariances[2]],
        );

        ALL_BLOCK_SLOTS
            .into_iter()
            .map(|slots| {
                let derivative = (forward.component(slots, [0; 3])
                    - backward.component(slots, [0; 3]))
                    / (2.0 * step);
                let mut magnitude = derivative.norm() + source.component(slots, [0; 3]).norm();
                let mut linear = Complex64::new(0.0, 0.0);
                for leg in 0..3 {
                    let mut moved = slots;
                    let term = match slots[leg] {
                        PhaseSpaceSlotKind::Field => {
                            moved[leg] = PhaseSpaceSlotKind::Momentum;
                            alpha.component(moved, [0; 3]) / hubble
                        }
                        PhaseSpaceSlotKind::Momentum => {
                            moved[leg] = PhaseSpaceSlotKind::Field;
                            -3.0 * alpha.component(slots, [0; 3])
                                - (wavenumbers[leg] / scale).powi(2) / hubble
                                    * alpha.component(moved, [0; 3])
                        }
                    };
                    magnitude += term.norm();
                    linear += term;
                }
                (derivative - linear - source.component(slots, [0; 3])).norm() / magnitude
            })
            .fold(0.0, f64::max)
    }

    #[test]
    fn vertex_rates_match_the_flat_single_field_background_derivative() {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0)]),
            FlatFieldSpace,
        );
        let (phi, velocity) = (2.0_f64, 0.5_f64);
        let at = |phi: f64, velocity: f64, e_folds: f64| {
            MultifieldState::new(
                [Field::new(phi)],
                [FieldVelocity::new(velocity)],
                ScaleFactor::new(e_folds.exp()),
            )
        };
        let center = at(phi, velocity, 0.0);
        let hubble = model.hubble(&center).into_value();
        // phi'' = -3 H phi' - m^2 phi and H' = -phi'^2 / 2 with m = M_pl = 1.
        let acceleration = -3.0 * hubble * velocity - phi;
        let hubble_rate = -velocity * velocity / 2.0;
        let step: f64 = 1.0e-4;
        let neighborhood = BackgroundNeighborhood {
            before: at(
                phi - step * velocity / hubble,
                velocity - step * acceleration / hubble,
                -step,
            ),
            center,
            after: at(
                phi + step * velocity / hubble,
                velocity + step * acceleration / hubble,
                step,
            ),
            e_fold_step: step,
        };
        let k = 3.0;
        let triangle = MomentumTriangle::try_new(k, k, k).unwrap();

        let rates = canonical_subhorizon_vertex_rates(&model, &neighborhood, triangle).unwrap();

        // Equilateral flat single field: a^2 A_fast = -phi' k^2 / (4 H) and
        // C = -phi'/H + 3 phi'^3 / (32 H^3); d/dN = H^-1 d/dt.
        let ratio_rate = (acceleration / hubble - velocity * hubble_rate / hubble.powi(2)) / hubble;
        let expected_fast = -k * k / 4.0 * ratio_rate;
        let expected_c = (-acceleration / hubble
            + velocity * hubble_rate / hubble.powi(2)
            + 9.0 / 32.0 * velocity.powi(2) * acceleration / hubble.powi(3)
            - 9.0 / 32.0 * velocity.powi(3) * hubble_rate / hubble.powi(4))
            / hubble;
        let actual_fast = rates.fast_a_times_scale_squared.components()[0][0][0].into_value();
        assert!(
            (actual_fast - expected_fast).abs() <= 1.0e-7 * expected_fast.abs(),
            "{actual_fast} != {expected_fast}"
        );
        for leg in 0..3 {
            let actual = rates.c_field_on_leg[leg].components()[0][0][0].into_value();
            assert!(
                (actual - expected_c).abs() <= 1.0e-7 * expected_c.abs(),
                "{actual} != {expected_c}"
            );
        }
        assert!(matches!(
            canonical_subhorizon_vertex_rates(
                &model,
                &BackgroundNeighborhood {
                    e_fold_step: 0.0,
                    ..neighborhood
                },
                triangle,
            ),
            Err(ThreePointInitialConditionError::InvalidEFoldStep(_))
        ));
    }

    #[test]
    fn vertex_rates_include_the_connection_in_curved_field_space() {
        let (model, center) = curved_initial_condition_model();
        let hubble = model.hubble(&center).into_value();
        let derivative = System::derivative(&model, crate::Time::new(0.0), &center);
        let step: f64 = 1.0e-4;
        let shifted = |sign: f64| {
            let dn = sign * step;
            MultifieldState::new(
                core::array::from_fn(|i| {
                    center.phi[i] + Field::new(dn * derivative.dphi_dt[i].into_value() / hubble)
                }),
                core::array::from_fn(|i| {
                    center.pi[i]
                        + FieldVelocity::new(dn * derivative.dpi_dt[i].into_value() / hubble)
                }),
                ScaleFactor::new(center.a.into_value() * dn.exp()),
            )
        };
        let neighborhood = BackgroundNeighborhood {
            before: shifted(-1.0),
            center,
            after: shifted(1.0),
            e_fold_step: step,
        };
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();
        let rates = canonical_subhorizon_vertex_rates(&model, &neighborhood, triangle).unwrap();

        // S = T^{IJK} v_I v_J v_K with v_I = G_IL phi'^L is a coordinate scalar,
        // so dS/dN = (D_N T)^{IJK} v v v + T^{IJK} D_N(v_I v_J v_K). For a
        // non-equilateral triangle T is not symmetric under index exchange
        // alone, so each slot of the product rule is kept separately.
        let lowered = |state: &MultifieldState<2>| -> [f64; 2] {
            let metric = model.geometry.metric(&state.phi);
            core::array::from_fn(|i| {
                (0..2)
                    .map(|l| metric.components()[i][l].into_value() * state.pi[l].into_value())
                    .sum()
            })
        };
        let scalar = |tensor: &[[[Quantity<M1, f64>; 2]; 2]; 2], v: [f64; 2]| -> f64 {
            (0..8)
                .map(|n| {
                    let [i, j, k] = [n / 4, (n / 2) % 2, n % 2];
                    tensor[i][j][k].into_value() * v[i] * v[j] * v[k]
                })
                .sum()
        };
        let slow_at = |state: &MultifieldState<2>| {
            *canonical_subhorizon_initial_vertices(&model, state, triangle)
                .unwrap()
                .slow_a
                .components()
        };
        let finite = (scalar(&slow_at(&neighborhood.after), lowered(&neighborhood.after))
            - scalar(
                &slow_at(&neighborhood.before),
                lowered(&neighborhood.before),
            ))
            / (2.0 * step);

        let connection = model.geometry.christoffel(&center.phi);
        let metric = model.geometry.metric(&center.phi);
        let covariant_acceleration: [f64; 2] = core::array::from_fn(|l| {
            derivative.dpi_dt[l].into_value()
                + (0..4)
                    .map(|n| {
                        let [m, q] = [n / 2, n % 2];
                        connection.components()[l][m][q].value()
                            * center.pi[m].into_value()
                            * center.pi[q].into_value()
                    })
                    .sum::<f64>()
        });
        let v = lowered(&center);
        let dv: [f64; 2] = core::array::from_fn(|i| {
            (0..2)
                .map(|l| metric.components()[i][l].into_value() * covariant_acceleration[l])
                .sum::<f64>()
                / hubble
        });
        let slow = slow_at(&center);
        let transported = scalar(rates.slow_a.components(), v)
            + (0..8)
                .map(|n| {
                    let [i, j, k] = [n / 4, (n / 2) % 2, n % 2];
                    slow[i][j][k].into_value()
                        * (dv[i] * v[j] * v[k] + v[i] * dv[j] * v[k] + v[i] * v[j] * dv[k])
                })
                .sum::<f64>();

        assert!(
            (finite - transported).abs() <= 1.0e-4 * finite.abs(),
            "{finite:e} != {transported:e}"
        );

        // The connection part of the rate (about 0.8% here) is far above that
        // tolerance, so a coordinate-only derivative would fail the check.
        let before_slow = slow_at(&neighborhood.before);
        let after_slow = slow_at(&neighborhood.after);
        let coordinate: [[[Quantity<M1, f64>; 2]; 2]; 2] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| {
                    (after_slow[i][j][k] - before_slow[i][j][k]) / (2.0 * step)
                })
            })
        });
        let connection_part = scalar(rates.slow_a.components(), v) - scalar(&coordinate, v);
        assert!(connection_part.abs() > 2.0e-3 * finite.abs());
    }

    #[test]
    fn subhorizon_rate_terms_solve_de_sitter_transport_with_drifting_coefficients() {
        let shallow = de_sitter_transport_residual(1.0e3, 1.0, true);
        let deep = de_sitter_transport_residual(1.0e4, 1.0, true);
        let frozen = de_sitter_transport_residual(1.0e4, 1.0, false);

        assert!(shallow < 2.0e-5, "residual {shallow:e} at k/(aH) ~ 1e3");
        assert!(deep < 2.0e-7, "residual {deep:e} at k/(aH) ~ 1e4");
        assert!(
            frozen > 1.0e-3,
            "frozen residual {frozen:e} should expose the drift"
        );
    }

    #[test]
    fn subhorizon_initial_condition_solves_de_sitter_transport_to_retained_order() {
        let shallow = de_sitter_transport_residual(1.0e3, 0.0, false);
        let deep = de_sitter_transport_residual(1.0e4, 0.0, false);

        assert!(shallow < 2.0e-5, "residual {shallow:e} at k/(aH) ~ 1e3");
        assert!(deep < 2.0e-7, "residual {deep:e} at k/(aH) ~ 1e4");
    }

    #[test]
    fn subhorizon_initial_condition_matches_the_independent_in_in_values() {
        // Exact rational values from the tree-level in-in integral with
        // massless de Sitter mode functions, evaluated at k = (3, 2, 3/2),
        // a = 2, H = 1/2, using the scalar vertex values above.
        let expected = [
            -83.0 / 374_400.0,
            743.0 / 1_622_400.0,
            5_801.0 / 9_734_400.0,
            943.0 / 1_081_600.0,
            -97.0 / 2_995_200.0,
            -1.0 / 166_400.0,
            73.0 / 1_497_600.0,
            -6_659.0 / 8_652_800.0,
        ];
        let alpha = canonical_three_point_from_subhorizon_vertices(
            &scalar_initial_vertices(),
            &scalar_rates(0.0),
            [Mass::new(3.0), Mass::new(2.0), Mass::new(1.5)],
            ScaleFactor::new(2.0),
            Hubble::new(0.5),
        )
        .unwrap();

        for (slots, expected) in ALL_BLOCK_SLOTS.into_iter().zip(expected) {
            let actual = alpha.component(slots, [0; 3]);
            assert_eq!(actual.im, 0.0);
            assert!(
                (actual.re - expected).abs() <= 1.0e-12 * expected.abs(),
                "{slots:?}: {} != {expected}",
                actual.re
            );
        }
    }

    fn generic_block<D>(
        offset: f64,
    ) -> ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, D, f64, 2> {
        ThreePointBlock::new(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                core::array::from_fn(|k| Quantity::new((offset + (7 * i + 3 * j + k) as f64).sin()))
            })
        }))
    }

    fn generic_initial_vertices() -> CanonicalSubhorizonInitialVertices<2> {
        CanonicalSubhorizonInitialVertices {
            fast_a_times_scale_squared: generic_block(0.1),
            slow_a: generic_block(1.7),
            b_momentum_on_leg: [generic_block(2.3), generic_block(3.1), generic_block(4.9)],
            c_field_on_leg: [generic_block(5.3), generic_block(6.7), generic_block(7.1)],
        }
    }

    /// Relabels vertices after exchanging two Fourier legs, moving each leg's
    /// field index and its leg-resolved `B` and `C` entry together.
    fn exchange_vertex_legs(
        vertices: &CanonicalSubhorizonInitialVertices<2>,
        exchange: [usize; 2],
    ) -> CanonicalSubhorizonInitialVertices<2> {
        fn relabel<D: Copy>(
            block: &ThreePointBlock<
                FieldPerturbation,
                FieldPerturbation,
                FieldPerturbation,
                D,
                f64,
                2,
            >,
            exchange: [usize; 2],
        ) -> ThreePointBlock<FieldPerturbation, FieldPerturbation, FieldPerturbation, D, f64, 2>
        {
            ThreePointBlock::new(core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    core::array::from_fn(|k| {
                        let mut indices = [i, j, k];
                        indices.swap(exchange[0], exchange[1]);
                        block.components()[indices[0]][indices[1]][indices[2]]
                    })
                })
            }))
        }
        let leg = |r: usize| {
            if r == exchange[0] {
                exchange[1]
            } else if r == exchange[1] {
                exchange[0]
            } else {
                r
            }
        };
        CanonicalSubhorizonInitialVertices {
            fast_a_times_scale_squared: relabel(&vertices.fast_a_times_scale_squared, exchange),
            slow_a: relabel(&vertices.slow_a, exchange),
            b_momentum_on_leg: [0, 1, 2]
                .map(|r| relabel(&vertices.b_momentum_on_leg[leg(r)], exchange)),
            c_field_on_leg: [0, 1, 2].map(|r| relabel(&vertices.c_field_on_leg[leg(r)], exchange)),
        }
    }

    fn assert_three_point_close<P: PhaseSpaceSlot>(
        actual: &ThreePointPhaseSpaceCorrelation<2, P>,
        expected: &ThreePointPhaseSpaceCorrelation<2, P>,
        tolerance: f64,
    ) {
        for slots in ALL_BLOCK_SLOTS {
            let scale = (0..8)
                .map(|n| {
                    expected
                        .component(slots, [n / 4, (n / 2) % 2, n % 2])
                        .norm()
                })
                .fold(0.0, f64::max);
            for n in 0..8 {
                let indices = [n / 4, (n / 2) % 2, n % 2];
                let difference =
                    (actual.component(slots, indices) - expected.component(slots, indices)).norm();
                assert!(
                    difference <= tolerance * scale,
                    "{slots:?} {indices:?}: difference {difference:e}, scale {scale:e}"
                );
            }
        }
    }

    #[test]
    fn subhorizon_initial_condition_is_covariant_under_fourier_leg_exchange() {
        let vertices = generic_initial_vertices();
        let rates = CanonicalSubhorizonInitialVertices {
            fast_a_times_scale_squared: generic_block(8.3),
            slow_a: generic_block(9.1),
            b_momentum_on_leg: [
                generic_block(10.7),
                generic_block(11.3),
                generic_block(12.9),
            ],
            c_field_on_leg: [
                generic_block(13.1),
                generic_block(14.3),
                generic_block(15.7),
            ],
        };
        let wavenumbers = [Mass::new(4.0e3), Mass::new(3.0e3), Mass::new(2.0e3)];
        let (scale, hubble) = (ScaleFactor::new(1.5), Hubble::new(0.8));
        let alpha = canonical_three_point_from_subhorizon_vertices(
            &vertices,
            &rates,
            wavenumbers,
            scale,
            hubble,
        )
        .unwrap();

        let first_second = canonical_three_point_from_subhorizon_vertices(
            &exchange_vertex_legs(&vertices, [0, 1]),
            &exchange_vertex_legs(&rates, [0, 1]),
            [wavenumbers[1], wavenumbers[0], wavenumbers[2]],
            scale,
            hubble,
        )
        .unwrap();
        assert_three_point_close(&first_second, &alpha.swap_first_second(), 1.0e-14);

        let second_third = canonical_three_point_from_subhorizon_vertices(
            &exchange_vertex_legs(&vertices, [1, 2]),
            &exchange_vertex_legs(&rates, [1, 2]),
            [wavenumbers[0], wavenumbers[2], wavenumbers[1]],
            scale,
            hubble,
        )
        .unwrap();
        assert_three_point_close(&second_third, &alpha.swap_second_third(), 1.0e-14);
        assert_eq!(alpha.conjugate(), alpha);
    }

    /// A neighborhood whose outer states equal the slice. Its coordinate rates
    /// vanish, which is sufficient for tests of symmetry and validation.
    fn static_neighborhood(state: MultifieldState<2>) -> BackgroundNeighborhood<2> {
        BackgroundNeighborhood {
            before: state,
            center: state,
            after: state,
            e_fold_step: 1.0e-3,
        }
    }

    fn curved_initial_condition_model() -> (
        CovariantMultifieldSystem<QuadraticMultifieldPotential<2>, HyperbolicFieldSpace, 2>,
        MultifieldState<2>,
    ) {
        let model = CovariantMultifieldSystem::in_planck_units(
            QuadraticMultifieldPotential::new([Mass::new(1.0e-5), Mass::new(0.75e-5)]),
            HyperbolicFieldSpace::new(Field::new(3.0)).unwrap(),
        );
        let state = MultifieldState::new(
            [Field::new(2.0), Field::new(0.5)],
            [FieldVelocity::new(4.0e-6), FieldVelocity::new(-3.0e-6)],
            ScaleFactor::new(1.0),
        );
        (model, state)
    }

    #[test]
    fn subhorizon_vertices_evaluate_each_leg_ordering_before_raising() {
        let (model, state) = curved_initial_condition_model();
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();
        let vertices = canonical_subhorizon_initial_vertices(&model, &state, triangle).unwrap();
        let inverse = *model.geometry.inverse_metric(&state.phi).components();
        let legs_with_momentum_last = [
            [FourierLeg::Two, FourierLeg::Three, FourierLeg::One],
            [FourierLeg::One, FourierLeg::Three, FourierLeg::Two],
            [FourierLeg::One, FourierLeg::Two, FourierLeg::Three],
        ];

        for (leg, coefficient_legs) in legs_with_momentum_last.into_iter().enumerate() {
            // `first_raised` in coefficient order (x, y | z) raises x; raise y, z.
            let b = raise_second_and_third(
                canonical_cubic_b_coefficients_for_legs(&model, &state, triangle, coefficient_legs)
                    .unwrap()
                    .first_raised
                    .components(),
                &inverse,
            );
            let c = raise_second_and_third(
                canonical_cubic_c_coefficients_for_legs(&model, &state, triangle, coefficient_legs)
                    .unwrap()
                    .first_raised
                    .components(),
                &inverse,
            );
            let position = coefficient_legs.map(FourierLeg::index);
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let leg_indices = [i, j, k];
                        let [x, y, z] = position.map(|p| leg_indices[p]);
                        let actual_b = vertices.b_momentum_on_leg[leg].components()[i][j][k];
                        let actual_c = vertices.c_field_on_leg[leg].components()[i][j][k];
                        let expected_b = b[x][y][z];
                        let expected_c = c[x][y][z];
                        assert!(
                            (actual_b - expected_b).into_value().abs()
                                <= 1.0e-13 * expected_b.into_value().abs().max(1.0e-30)
                        );
                        assert!(
                            (actual_c - expected_c).into_value().abs()
                                <= 1.0e-13 * expected_c.into_value().abs().max(1.0e-30)
                        );
                    }
                }
            }
        }

        let full_a = raise_second_and_third(
            canonical_cubic_a_coefficients(&model, &state, triangle)
                .unwrap()
                .first_raised
                .components(),
            &inverse,
        );
        for n in 0..8 {
            let [i, j, k] = [n / 4, (n / 2) % 2, n % 2];
            let expected = full_a[i][j][k];
            let reconstructed = vertices.slow_a.components()[i][j][k]
                + vertices.fast_a_times_scale_squared.components()[i][j][k];
            let difference = (reconstructed - expected).into_value().abs();
            assert!(difference <= 1.0e-13 * expected.into_value().abs().max(1.0e-30));
        }

        // At fixed shape, only the fast part carries physical wavenumbers.
        let doubled = canonical_subhorizon_initial_vertices(
            &model,
            &state,
            MomentumTriangle::try_new(8.0, 6.0, 4.0).unwrap(),
        )
        .unwrap();
        for n in 0..8 {
            let [i, j, k] = [n / 4, (n / 2) % 2, n % 2];
            let slow = vertices.slow_a.components()[i][j][k].into_value();
            let fast = vertices.fast_a_times_scale_squared.components()[i][j][k].into_value();
            let doubled_slow = doubled.slow_a.components()[i][j][k].into_value();
            let doubled_fast =
                doubled.fast_a_times_scale_squared.components()[i][j][k].into_value();
            assert!((doubled_slow - slow).abs() <= 1.0e-12 * slow.abs().max(fast.abs()));
            assert!((doubled_fast - 4.0 * fast).abs() <= 1.0e-13 * fast.abs().max(1.0e-30));
        }
    }

    #[test]
    fn input_exchanged_coefficients_move_momenta_with_indices() {
        let (model, state) = curved_initial_condition_model();
        let triangle = MomentumTriangle::try_new(4.0, 3.0, 2.0).unwrap();
        let canonical = canonical_cubic_hamiltonian_coefficients(&model, &state, triangle).unwrap();
        let exchanged_legs = [FourierLeg::One, FourierLeg::Three, FourierLeg::Two];
        let b = canonical_cubic_b_coefficients_for_legs(&model, &state, triangle, exchanged_legs)
            .unwrap();
        let c = canonical_cubic_c_coefficients_for_legs(&model, &state, triangle, exchanged_legs)
            .unwrap();

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    assert_eq!(
                        canonical.b_first_raised_input_exchanged.components()[i][j][k],
                        b.first_raised.components()[i][k][j]
                    );
                    assert_eq!(
                        canonical.c_first_raised_input_exchanged.components()[i][j][k],
                        c.first_raised.components()[i][k][j]
                    );
                }
            }
        }
        assert_ne!(
            canonical.c_first_raised_input_exchanged.components()[0][0][1],
            canonical.c_first_raised.components()[0][1][0]
        );
    }

    #[test]
    fn equilateral_model_initial_condition_is_invariant_under_leg_exchange() {
        let (model, state) = curved_initial_condition_model();
        let hubble = model.hubble(&state).into_value();
        let k = 50.0 * hubble;
        let triangle = MomentumTriangle::try_new(k, k, k).unwrap();
        let alpha = canonical_three_point_initial_condition(
            &model,
            &static_neighborhood(state),
            triangle,
            20.0,
        )
        .unwrap();

        assert_three_point_close(&alpha.swap_first_second(), &alpha, 1.0e-13);
        assert_three_point_close(&alpha.swap_second_third(), &alpha, 1.0e-13);
        assert_eq!(alpha.conjugate(), alpha);
    }

    #[test]
    fn initial_condition_requires_every_leg_inside_the_requested_depth() {
        let (model, state) = curved_initial_condition_model();
        let hubble = model.hubble(&state).into_value();
        let triangle =
            MomentumTriangle::try_new(80.0 * hubble, 60.0 * hubble, 30.0 * hubble).unwrap();

        assert!(
            canonical_three_point_initial_condition(
                &model,
                &static_neighborhood(state),
                triangle,
                25.0
            )
            .is_ok()
        );
        match canonical_three_point_initial_condition(
            &model,
            &static_neighborhood(state),
            triangle,
            40.0,
        ) {
            Err(ThreePointInitialConditionError::InsufficientSubhorizonDepth {
                leg,
                k_over_ah,
                required,
            }) => {
                assert_eq!(leg, FourierLeg::Three);
                assert!((k_over_ah - 30.0).abs() < 1.0e-9);
                assert_eq!(required, 40.0);
            }
            other => panic!("unexpected result {other:?}"),
        }
        for invalid in [1.0, 0.5, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                canonical_three_point_initial_condition(
                    &model,
                    &static_neighborhood(state),
                    triangle,
                    invalid
                ),
                Err(ThreePointInitialConditionError::InvalidSubhorizonRequirement(_))
            ));
        }
        assert!(matches!(
            canonical_three_point_from_subhorizon_vertices(
                &scalar_initial_vertices(),
                &scalar_rates(0.0),
                [Mass::new(1.0), Mass::new(0.0), Mass::new(1.0)],
                ScaleFactor::new(1.0),
                Hubble::new(1.0),
            ),
            Err(ThreePointInitialConditionError::InvalidWavenumber {
                leg: FourierLeg::Two,
                ..
            })
        ));
    }
}
