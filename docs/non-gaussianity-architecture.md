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

The first deliverable is not a number called `f_NL`.  It is a phase-space
two-point transport implementation that reproduces the existing quantum mode
solver.  Three-point evolution is enabled only after that equivalence is
verified.

## Scope

The target is the connected tree-level primordial three-point function for
the scalar sector of a two-field sigma model,

```text
S = integral sqrt(-g) [M_pl^2 R / 2
                       - G_IJ partial(phi^I) partial(phi^J) / 2
                       - V(phi)].
```

The complete momentum dependence is retained.  Local, equilateral, and
squeezed values are derived evaluations of the same triangle-dependent
bispectrum, not interchangeable scalar labels.

Loops, tensor/scalar mixed correlators, non-Bunch-Davies initial states,
reheating transfer, and observational likelihood evaluation are outside this
milestone.

## Separate index domains

Three independent notions must never share a raw integer API:

1. field-space indices `I, J, K` transform with `G_IJ`;
2. phase-space indices `a, b, c` select field or momentum slots and carry a
   field-space index inside that slot;
3. Fourier-leg labels `1, 2, 3` identify `k_1`, `k_2`, and `k_3`.

The phase-space state is conceptually

```text
X^a = (Q^I, P^I),
```

where the precise momentum normalization is fixed once and used by initial
conditions, evolution tensors, covariance, and output.  A phase-space index
is therefore a tagged pair, not a third field-space variance marker.

The two-point and connected three-point objects are

```text
Sigma^{ab}(k)                    = <X^a(k) X^b(-k)>,
Alpha^{abc}(k_1, k_2, k_3)      = <X^a(k_1) X^b(k_2) X^c(k_3)>_connected.
```

Their storage types carry physical dimensions for each slot.  Momentum-leg
permutations move both the Fourier label and its associated phase-space slot.
They are not implemented as field-index permutations.

## Momentum triangles

`MomentumTriangle` owns three positive normalized wavenumbers and validates

```text
k_1 <= k_2 + k_3,
k_2 <= k_1 + k_3,
k_3 <= k_1 + k_2.
```

The canonical representation orders `k_1 >= k_2 >= k_3` and stores scale and
shape separately:

```text
k_t = k_1 + k_2 + k_3,
x_2 = k_2 / k_1,
x_3 = k_3 / k_1.
```

This prevents duplicate integrations of permuted triangles while retaining a
mapping back to caller leg order.  Degenerate squeezed boundaries require an
explicit tolerance; invalid or non-finite triangles return a structured error.

The initial benchmark grid contains at least equilateral, folded-boundary,
and squeezed configurations around the same observational pivot used by the
power-spectrum benchmark.

## Transport equations

Writing the nonlinear phase-space flow schematically as

```text
D_N X^a = u^a_b X^b + 1/2 u^a_bc X^b X^c + ...,
```

the equal-time tree-level hierarchy closes at three-point order as

```text
D_N Sigma^{ab}
  = u^a_c Sigma^{cb} + u^b_c Sigma^{ac},

D_N Alpha^{abc}
  = u^a_d Alpha^{dbc} + u^b_d Alpha^{adc} + u^c_d Alpha^{abd}
    + u^a_de Sigma^{db} Sigma^{ec}
    + u^b_de Sigma^{da} Sigma^{ec}
    + u^c_de Sigma^{da} Sigma^{eb}.
```

In the full subhorizon calculation the kernels retain their momentum
arguments and quantum real/imaginary structure.  Covariant derivatives along
the background act only on field-space content.  Connection and Riemann terms
must be derived from one declared convention; they are never reconstructed by
finite-differencing exported CSV rows.

The existing linear mode matrix remains the reference oracle.  The transport
two-point sector must reproduce its `QQ`, `QP`, `PQ`, and `PP` covariances for
shared wavenumbers before `u^a_bc` or `Alpha` is accepted.

## Initial conditions and time range

Correlation functions start in the same Bunch-Davies regime as the linear
solver.  The two-point initial condition is derived from the existing mode
normalization.  The three-point initial condition and source kernels come
from the cubic action/quantum transport prescription; setting `Alpha = 0` at
horizon exit is not an acceptable substitute because it discards intrinsic
horizon-crossing interactions and general momentum shapes.

All three legs of a triangle use one common start time deep inside the horizon
of its largest relevant scale, with an explicit convergence study in
`x_init`.  The shared P2 background trajectory supplies every Runge--Kutta
stage.

## Projection to curvature observables

Field/phase-space correlators are not yet the observable curvature
bispectrum.  At the requested output slice, a separately tested nonlinear
gauge transformation supplies first- and second-order projection coefficients
`N_a` and `N_ab`.  The output records:

- the field/phase-space `Alpha` before projection;
- `P_zeta(k)` after projection;
- dimensional `B_zeta(k_1,k_2,k_3)`;
- the exact reduced-amplitude convention used for any reported `f_NL`;
- the output slice, triangle convention, background/model provenance, and
  whether reheating transfer is absent.

No single `f_NL` column may be emitted without a shape and normalization
identifier.

## Required mathematical capabilities

The implementation sequence adds only what the transport kernels consume:

- third-order automatic differentiation for the potential;
- a symmetric covariant `V_;IJK` with dimension `M`;
- rank-five field tensors and `nabla_M R^I_JKL` where required;
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
2. **P3b — third-order model geometry.** Add third-order AD, `V_;IJK`, and the
   curvature-derivative capabilities actually required by the selected
   kernels, with analytic flat and curved tests.
3. **P3c — two-point transport oracle.** Evolve the phase-space covariance on
   the shared background and reproduce every block of the existing mode
   solver across subhorizon, horizon-exit, turn, and final epochs.
4. **P3d — three-point transport.** Implement `u^a_bc`, quantum initial
   conditions, sourced `Alpha`, and triangle/permutation regression tests.
5. **P3e — curvature projection and benchmark.** Implement the nonlinear
   projection, publish `B_zeta` with convention-labeled reduced amplitudes,
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
