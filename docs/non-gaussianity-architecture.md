# Non-Gaussianity architecture

## Decision

The non-Gaussianity milestone will use tree-level phase-space transport.  It
will not begin with a slow-roll, superhorizon-only `delta N` approximation and
will not implement a separate in-in quadrature for each triangle.

This choice preserves the capabilities already demonstrated by the linear
solver: curved field space, massive entropy directions, temporary violations
of slow roll, and sharp turns.  The transport equations are ordinary
differential equations for equal-time correlation functions and can read the
same immutable background trajectory introduced by P2.

The first acceptance gate was not a number called $f_{\mathrm{NL}}$. It was a
phase-space two-point transport implementation that reproduces the existing
quantum mode solver. That P3c gate is now complete, so P3d can build
three-point evolution on a validated linear convention.

## Scope

The target is the connected tree-level primordial three-point function for
the scalar sector of a two-field sigma model,

$$
S=\int d^4x\,\sqrt{-g}\left[
\frac{M_{\mathrm{Pl}}^2}{2}R
-\frac12G_{IJ}(\phi)g^{\mu\nu}
  \partial_\mu\phi^I\partial_\nu\phi^J
-V(\phi)\right].
$$

Here $G_{IJ}$ is the field-space metric, $I,J\in\{1,2\}$ label field
coordinates, and $M_{\mathrm{Pl}}$ is the reduced Planck mass. Natural units
$c=\hbar=1$ are used throughout the crate.

The complete momentum dependence is retained.  Local, equilateral, and
squeezed values are derived evaluations of the same triangle-dependent
bispectrum, not interchangeable scalar labels.

Loops, tensor/scalar mixed correlators, non-Bunch-Davies initial states,
reheating transfer, and observational likelihood evaluation are outside this
milestone.

## Separate index domains

Three independent notions must never share a raw integer API:

1. field-space indices $I,J,K$ transform with $G_{IJ}$;
2. phase-space indices $a,b,c$ select field or momentum slots and carry a
   field-space index inside that slot;
3. Fourier-leg labels $1,2,3$ identify $\mathbf{k}_1$, $\mathbf{k}_2$, and
   $\mathbf{k}_3$.

The phase-space state is conceptually

$$
X^a=(Q^I,P^I),\qquad Q^I=U^I,\qquad P^I=D_tU^I.
$$

Here $U^I$ is the covariant field perturbation and $D_t$ is the field-space
covariant cosmic-time derivative. Therefore $P^I$ is a velocity perturbation,
not an $a^3$-weighted canonical momentum. A phase-space index is a tagged pair
such as $(Q,I)$ or $(P,I)$, not a third field-space variance marker.

At cubic order there is a second, subtler distinction. The usual Hamiltonian
transport formulas use the rescaled conjugate momentum $\Pi^I=a^{-3}p^I$.
Although $\Pi^I=D_tQ^I$ at linear order, cubic interactions make them differ
at quadratic order in perturbations. This is exactly the order entering
$u^a{}_{bc}$ and the three-point function. The derivation and the rule against
silently reinterpreting the existing slot marker are fixed in
[`cubic-kernel-conventions.md`](cubic-kernel-conventions.md).

The two-point and connected three-point objects are

$$
\Sigma^{ab}(k)
=\left\langle X^a(\mathbf{k})X^b(-\mathbf{k})\right\rangle,
$$

$$
\alpha^{abc}(k_1,k_2,k_3)
=\left\langle
X^a(\mathbf{k}_1)X^b(\mathbf{k}_2)X^c(\mathbf{k}_3)
\right\rangle_c,
$$

after stripping off the momentum-conserving delta function. The subscript
$c$ denotes the connected part, and
$\mathbf{k}_1+\mathbf{k}_2+\mathbf{k}_3=0$.

Their storage types carry physical dimensions for each slot.  Momentum-leg
permutations move both the Fourier label and its associated phase-space slot.
They are not implemented as field-index permutations.

## Momentum triangles

`MomentumTriangle` owns the three positive magnitudes
$k_r=|\mathbf{k}_r|$ and validates

$$
k_1\le k_2+k_3,\qquad
k_2\le k_1+k_3,\qquad
k_3\le k_1+k_2.
$$

The canonical representation orders `k_1 >= k_2 >= k_3` and stores scale and
shape separately:

$$
k_t=k_1+k_2+k_3,\qquad
x_2=\frac{k_2}{k_1},\qquad
x_3=\frac{k_3}{k_1}.
$$

This prevents duplicate integrations of permuted triangles while retaining a
mapping back to caller leg order.  Degenerate squeezed boundaries require an
explicit tolerance; invalid or non-finite triangles return a structured error.

The initial benchmark grid contains at least equilateral, folded-boundary,
and squeezed configurations around the same observational pivot used by the
power-spectrum benchmark.

## Transport equations

Writing the nonlinear phase-space flow schematically as

$$
D_NX^a=u^a{}_bX^b+\frac12u^a{}_{bc}X^bX^c+\cdots,
\qquad D_N\equiv H^{-1}D_t,
$$

where $u^a{}_b$ is the linear flow Jacobian and $u^a{}_{bc}$ is its quadratic
extension derived from the cubic Hamiltonian, the equal-time tree-level
hierarchy closes at three-point order:

$$
D_N\Sigma^{ab}
=u^a{}_c\Sigma^{cb}+u^b{}_c\Sigma^{ac},
$$

$$
\begin{aligned}
D_N\alpha^{abc}={}&
u^a{}_d\alpha^{dbc}
+u^b{}_d\alpha^{adc}
+u^c{}_d\alpha^{abd}\\
&+u^a{}_{de}\Sigma^{db}\Sigma^{ec}
+u^b{}_{de}\Sigma^{da}\Sigma^{ec}
+u^c{}_{de}\Sigma^{da}\Sigma^{eb}.
\end{aligned}
$$

Momentum arguments are suppressed only in this overview. They are written
explicitly in [`three-point-transport.md`](three-point-transport.md), where the
ordering of the two input legs of each $u^a{}_{bc}$ is fixed.

In the full subhorizon calculation the kernels retain their momentum
arguments and quantum real/imaginary structure.  Covariant derivatives along
the background act only on field-space content.  Connection and Riemann terms
must be derived from one declared convention; they are never reconstructed by
finite-differencing exported CSV rows.

The existing linear mode matrix remains the reference oracle.  The transport
two-point sector must reproduce its `QQ`, `QP`, `PQ`, and `PP` covariances for
shared wavenumbers before $u^a{}_{bc}$ or $\alpha^{abc}$ is accepted.

## Initial conditions and time range

Correlation functions start in the same Bunch--Davies regime as the linear
solver.  The two-point initial condition is derived from the existing mode
normalization.  The three-point initial condition and source kernels come
from the cubic action/quantum transport prescription; setting $\alpha=0$ at
horizon exit is not an acceptable substitute because it discards intrinsic
horizon-crossing interactions and general momentum shapes.

All three legs of a triangle use one common start time deep inside the horizon.
Writing $x_r=k_r/(aH)$, the start criterion must be stated in terms of all
three $x_r$ and tested by moving the common initial slice deeper. The shared
P2 background trajectory supplies every Runge--Kutta stage.

## Projection to curvature observables

Field/phase-space correlators are not yet the observable curvature
bispectrum.  At the requested output slice, a separately tested nonlinear
gauge transformation supplies first- and second-order projection coefficients
$N_a$ and $N_{ab}$. Schematically,

$$
\zeta=N_aX^a+\frac12N_{ab}X^aX^b+\cdots,
$$

so the curvature bispectrum receives both an intrinsic contribution from
$\alpha^{abc}$ and a gauge-transformation contribution quadratic in
$\Sigma^{ab}$:

$$
B_\zeta
=N_aN_bN_c\alpha^{abc}
+N_{ab}N_cN_d\bigl(\Sigma^{ac}\Sigma^{bd}
+\text{leg permutations}\bigr).
$$

The precise momentum arguments and normalization will be fixed in P3e rather
than inferred from this schematic equation. The output records:

- the field/phase-space $\alpha^{abc}$ before projection;
- $P_\zeta(k)$ after projection;
- dimensional $B_\zeta(k_1,k_2,k_3)$;
- the exact reduced-amplitude convention used for any reported
  $f_{\mathrm{NL}}$;
- the output slice, triangle convention, background/model provenance, and
  whether reheating transfer is absent.

No single $f_{\mathrm{NL}}$ column may be emitted without a shape and normalization
identifier.

## Required mathematical capabilities

The implementation sequence adds only what the transport kernels consume:

- third-order automatic differentiation for the potential;
- a symmetric covariant $V_{;(IJK)}$ with dimension $M$;
- rank-five field tensors and $\nabla_M R^I{}_{JKL}$ where required;
- phase-space tensors whose slot dimensions are explicit;
- complex Hermitian two-point blocks and permutation-aware complex
  three-point blocks;
- structured validation for Hermiticity, reality, permutation symmetry, and
  triangle closure.

Third derivative order and tensor rank remain separate concepts in the API.

## Ordered implementation milestones

1. **P3a — momentum and phase-space types.** Add validated momentum triangles,
   tagged phase-space indices, dimensioned two-/three-point containers, and
   permutation tests. No dynamics yet.
2. **P3b — third-order model geometry.** Add third-order AD, $V_{;(IJK)}$, and the
   curvature-derivative capabilities actually required by the selected
   kernels, with analytic flat and curved tests.
3. **P3c — two-point transport oracle.** Evolve the phase-space covariance on
   the shared background and reproduce every block of the existing mode
   solver across subhorizon, horizon-exit, turn, and final epochs.
4. **P3d — three-point transport.** Implement $u^a{}_{bc}$, quantum initial
   conditions, sourced $\alpha^{abc}$, and triangle/permutation regression
   tests.
5. **P3e — curvature projection and benchmark.** Implement the nonlinear
   projection, publish $B_\zeta$ with convention-labeled reduced amplitudes,
   and compare a canonical model with an independent transport result.

Each milestone has its own public types and tests.  A later milestone cannot
hide an unverified approximation made in an earlier one.

### P3a implementation status

P3a is complete.  The `transport` module provides:

- `MomentumTriangle`, which validates positive finite legs and relative
  closure tolerance, stores perimeter and canonical shape coordinates, and
  retains both directions of the caller/canonical leg permutation;
- compile-time `FieldPerturbation` and `MomentumPerturbation` slot markers,
  bounded `PhaseSpaceIndex<S, N>`, and an explicitly tagged erased index;
- dimension-carrying `TwoPointBlock` and `ThreePointBlock` containers;
- transpose and Fourier-leg permutation operations whose return types move
  the phase-space slot markers together with their component axes.

Unit tests cover equilateral/general, folded-boundary, squeezed, invalid, and
tolerance-assisted triangles, index bounds, and component permutations.  A
compile-fail case proves that a `QP` block cannot be assigned to a `QQ` block.
No transport dynamics or momentum normalization is implied by these types;
that enters through P3b/P3c.

### P3b implementation status

P3b is complete:

- `ThirdOrderJet<N>` carries a value, gradient, Hessian, and fully symmetric
  third derivative without changing the established `Jet<N>` representation;
- arithmetic, reciprocal, square root, `exp`, `ln`, `cos`, and `tanh` propagate
  derivatives through third order;
- `ThriceDifferentiableMultifieldPotential` evaluates a generic potential once
  to return dimensioned `V`, `V_,I`, `V_,IJ`, and coordinate `V_,IJK`;
- analytic tests cover mixed-index permutations, rational and radical chain
  rules, transcendental functions, a mixed cubic, and a vanishing quadratic
  third derivative.

- `ChristoffelDerivatives` represents $\partial_M\Gamma^I{}_{JK}$ explicitly as a
  coordinate-dependent non-tensor, and `ThirdOrderFieldSpaceGeometry` supplies
  it analytically;
- `CovariantMultifieldSystem` forms $\nabla_KV_{;IJ}$ and then explicitly
  symmetrizes it to the $V_{;(IJK)}$ combination selected by the cubic action;
- `FieldTensor5` and `CovariantRiemannDerivative` represent
  $\nabla_MR^I{}_{JKL}$ with variance and dimension fixed by the type system;
- `CurvatureDerivativeGeometry` is a separate capability, implemented
  analytically by the built-in flat and constant-curvature hyperbolic spaces.

Finite-difference tests independently check the hyperbolic connection
derivatives and every component of the symmetrized covariant potential third
derivative.  Flat-space reduction, $\nabla R=0$ for both built-in geometries,
full permutation symmetry, and rank-five variance rejection are also fixed by
tests. P3b is complete; P3c subsequently validated the linear transport
convention consumed by P3d.

### P3c implementation status

P3c is complete. The convention is fixed to $Q^I=U^I$ and
$P^I=D_tU^I$; $P$ is neither a conformal derivative nor an $a^3$-weighted
canonical momentum. `TwoPointPhaseSpaceCovariance` carries dimensioned complex
`QQ`, `QP`, `PQ`, and `PP` blocks, and
`two_point_covariance_from_mode_state` derives all four from one complete
mode/background state. Bunch--Davies closed forms and the identities
$\Sigma_{QQ}=\Sigma_{QQ}^\dagger$,
$\Sigma_{PP}=\Sigma_{PP}^\dagger$, and
$\Sigma_{PQ}=\Sigma_{QP}^\dagger$ are tested.

The detailed conversion and remaining acceptance stages are specified in
[`two-point-transport-oracle.md`](two-point-transport-oracle.md).
`TwoPointTransportSystem` now evolves the covariance directly on the same P2
background interpolation used by the mode solver, with no mode matrix in its
state. All four blocks agree at subhorizon, horizon-exit, maximum-turn, and
final epochs within the declared `2e-7` regression bound. Adaptive-tolerance
control preserves the physical dimension of every block; tightening the
relative tolerance from `1e-7` to `1e-10` changes the final covariance by less
than `2e-6`. Starts at 3, 4, and 5 e-folds before exit show successively
smaller differences, with each refinement shrinking by more than a factor of
two. P3d is the next milestone.

### P3d implementation status

P3d is complete. The first slice implemented the complete dimensioned
`ThreePointPhaseSpaceCorrelation` state with separate `QQQ` through `PPP`
blocks, complex conjugation, Runge--Kutta arithmetic, and Fourier-leg
transpositions that move their phase-space slots and component axes together.
A compile-fail test fixes the distinct `M^-3`, `M^-2`, `M^-1`, and `M^0`
dimensions of the delta-function-stripped three-point blocks.
`ThreePointHomogeneousTransportSystem` now applies
the P3c linear flow to all three canonical triangle legs on the shared P2
background, and is regression-tested for leg-permutation covariance. No
complete cubic kernel or physical bispectrum is claimed yet. The dimensioned
eight-block `QuadraticPhaseSpaceKernel` and the three ordered $u\Sigma\Sigma$ contractions
are implemented together with the complete model-derived coefficient layer. Its
$A_{IJK}$, $B_{IJK}$, and $C_{IJK}$ coefficients are evaluated from the model,
background, geometry, and triangle with typed dimensions, analytic flat
single-field tests, isolated Riemann-contraction regressions, and explicit
rank-five $\nabla R$ index routing. The
staged contract and next kernel-construction step are documented in
[`three-point-transport.md`](three-point-transport.md). Source-level
reconciliation found that the published eight-block cubic-Hamiltonian kernel
is expressed in $(Q,\Pi)$ rather than the already validated $(Q,D_tQ)$
variables. [`cubic-kernel-conventions.md`](cubic-kernel-conventions.md) now
records the full $A$, $B$, and $C$ formulas, ordered block map, and the explicit
implementation decision that must precede the physical kernel. That decision
is now implemented using the distinct `CanonicalMomentumPerturbation` marker
and canonical aliases for the covariance, three-point correlation, and
quadratic kernel. The typed source contraction rejects canonical/velocity
mixing at compile time.

The final P3d slices added the subhorizon quantum initial condition for all
eight canonical blocks, including first-order e-fold rates of the cubic
tensors, and `CanonicalThreePointTransportSystem`, which evolves the sourced
three-point function together with the three leg covariances. On the
reference benchmark it converges in start depth as $x_{\mathrm{init}}^{-1}$ and
in tolerance, and it transforms as a tensor under a nonlinear field
redefinition to $4\times10^{-10}$. The details and remaining limitations are
recorded in [`three-point-transport.md`](three-point-transport.md).

### P3e implementation status

P3e is complete. The projection to $\zeta$ uses the covariant "simple"
second-order gauge transformation on a superhorizon slice, and reports
$P_\zeta$, $B_\zeta$, and $f_{\mathrm{NL}}$ with an explicit sign and
normalization convention. The momentum arguments left open above are fixed
in [`curvature-projection.md`](curvature-projection.md), which also records
the validation against the mode solver, Maldacena's single-field result, and
the squeezed consistency relation. The axion-quartic model reproduces the
published result of Dias, Frazer, Mulryne, and Seery, `run_bispectrum` writes
convention-labeled CSV output, and an isosceles grid demonstrates shape
resolution.

## Acceptance boundary

The non-Gaussianity milestone is complete only when:

- two-point transport agrees with the existing mode solver within a declared
  numerical tolerance;
- three-point output obeys permutation and reality identities;
- single-field slow-roll and a documented two-field benchmark agree with
  independent analytic or published numerical references;
- start depth, integration tolerance, and triangle-grid convergence are shown;
- CSV metadata makes the gauge, output slice, shape, normalization, and
  no-reheating scope explicit.

Until then, the repository must describe the feature as transport
infrastructure, not as a prediction of primordial non-Gaussianity.

**Status (2026-09-30):** every condition above is met; see
[`curvature-projection.md`](curvature-projection.md) for the single-field,
axion-quartic, grid-convergence, and CSV-metadata evidence and for the scope
and limitations of the prediction.

## Primary references

- J.-O. Gong and T. Tanaka,
  [A covariant approach to general field space metric in multi-field inflation](https://arxiv.org/abs/1101.4809).
- J. Elliston, D. Seery, and R. Tavakol,
  [The inflationary bispectrum with curved field-space](https://arxiv.org/abs/1208.6011).
- D. J. Mulryne,
  [Transporting non-Gaussianity from sub to super-horizon scales](https://arxiv.org/abs/1302.3842).
- M. Dias, J. Frazer, and D. Seery,
  [Computing observables in curved multifield models of inflation](https://arxiv.org/abs/1502.03125).
- M. Dias, J. Frazer, D. Mulryne, and D. Seery,
  [Numerical evaluation of the bispectrum in multiple field inflation](https://arxiv.org/abs/1609.00379).
- S. Butchers and D. Seery,
  [Numerical evaluation of inflationary 3-point functions on curved field space](https://arxiv.org/abs/1803.10563).
