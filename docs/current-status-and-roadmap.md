# Current status and priority roadmap

This document is the single entry point for deciding what to do next.  The
detailed mathematical contracts remain in the design documents linked below;
this page records what is implemented, what is only designed, and the order in
which unfinished work should proceed.

## Status on 2026-09-30

The first scientific milestone is complete: `phys-type` can evolve a
covariant two-field inflationary background and its linear quantum
perturbations through a turn, then report adiabatic, entropy, and cross power
spectra.

The implemented path includes:

- compile-time mass dimensions and field-space index variance;
- first- and second-order automatic differentiation of potentials;
- covariant multifield backgrounds with metric, connection, and Riemann
  curvature contributions;
- a complex two-solution mode matrix with Bunch-Davies initialization and
  Wronskian diagnostics;
- Dormand-Prince 5(4) error control with dimensioned, per-component-group
  absolute tolerances;
- conformal-time evolution of `v = a U`, with `ln a` in the integrated
  background state;
- Hermitian covariance, curvature/entropy projection, and physical validity
  checks;
- a calibrated double-quadratic benchmark and reproducible PNG/PDF figures.

The implementation sequence 1--15 in
[`two-field-perturbation-design.md`](two-field-perturbation-design.md) is
complete.  Its conformal-time reference implementation is commit `041f2bf`.
At that commit, 133 unit tests and 17 compile-fail cases pass, as do formatting
and Clippy.
The conformal-time and former cosmic-time formulations agree in the final
curvature spectrum to a maximum relative difference of `3.1e-10` over the
41-point benchmark grid.

## Executive scientific snapshot

The repository currently has two distinct scientific boundaries:

| Capability | Status | Meaning |
| --- | --- | --- |
| Two-field background and linear perturbations | complete | produces validated adiabatic, entropy, cross, and curvature power spectra |
| Direct two-point phase-space transport | complete | reproduces the independent mode solver for all complex covariance blocks |
| Three-point transport (P3d) | complete | sourced canonical evolution converges on the reference benchmark and transforms as a tensor under field redefinitions |
| Curvature projection (P3e) | complete | $P_\zeta$, $B_\zeta$, and $f_{\mathrm{NL}}$ are validated in single-field limits and against the published axion-quartic result, with labeled CSV output and grid convergence |
| Primordial bispectrum prediction | available within scope | tree level, during inflation, no reheating; see `curvature-projection.md` |
| Nonminimally assisted inflation | in progress | all four loop-inflation branches reproduced; Class II remains |

Consequently, the crate can presently make power-spectrum predictions within
its documented inflation/no-reheating scope. It must not yet describe its
three-point output as a prediction for $B_\zeta$ or $f_{\mathrm{NL}}$.

### Next-session handoff

The implementation boundary is the branch `feat/nonminimal-assisted-inflation`,
based on merge commit `6ec9e7e` (PR #28). The validated state is 212 unit
tests plus 24 compile-fail cases, with ten expensive ignored tests that pass in
release mode (`cargo test --release -- --ignored`; they are not run by CI).
Clippy with warnings denied, rustfmt, and rustdoc pass locally.

The first external model validation after P3 is under way; see
[`nonminimal-assisted-inflation.md`](nonminimal-assisted-inflation.md):

1. `NonminimalAssistedFieldSpace` and `NonminimalLoopInflationPotential`
   implement the Einstein-frame model of arXiv:2302.05866 as ordinary
   providers. Finite-difference tests check every geometric quantity.
2. $\zeta$ is not conserved until the end of inflation, so observables are
   evaluated at the end. There, the mode solver, a numerical $\delta N$, and the
   reference's Eq. (3.1) agree on $n_s$. They agree across the coupling range
   once the exact $s_e$ is used. The slow-roll e-fold relation shifts the
   $n_s=0.965$ coupling from 0.0228 to 0.0222.
3. The reference's Eq. (4.2) typo and its Lyth--Rodriguez $f_{\mathrm{NL}}$ sign
   are recorded in [`literature-discrepancies.md`](literature-discrepancies.md),
   entries 6--9.

4. The $f_{\mathrm{NL}}$ difference is resolved: a phase-space $\delta N$ with the
   evolved covariance reproduces the transport to 0.05% for equilateral and
   squeezed shapes. The crossing-time $\delta N$ and Eq. (3.5) underestimate
   $|f_{\mathrm{NL}}|$ by about 20% because the assistant direction is tachyonic.

5. All four $(\alpha,s_*)$ branches of Sec. 4.1 follow the same pattern. The
   $n_s=0.965$ coupling moves up by 1.6--3.2% relative to the reference's
   slow-roll prescription. The transport $f_{\mathrm{NL}}$ is confirmed by
   phase-space $\delta N$ in every branch.

Resume with step 3 of the plan in that document: a general end-of-inflation
event and the Class II power-law and hybrid examples.

### Reading order for the active P3 work

1. [`non-gaussianity-architecture.md`](non-gaussianity-architecture.md)
   explains the physical target, index domains, transport hierarchy, scope,
   and P3a--P3e acceptance order.
2. [`two-point-transport-oracle.md`](two-point-transport-oracle.md) fixes the
   phase-space normalization and records why P3c is accepted.
3. [`three-point-transport.md`](three-point-transport.md) gives the detailed
   P3d equations, maps them to Rust types, and separates the implemented
   physical cubic kernel, basis conversion, quantum initial condition, and
   sourced evolution, and records their convergence and covariance tests.
4. [`field-space-curvature-design.md`](field-space-curvature-design.md) fixes
   the Riemann convention and the higher-rank geometry consumed by that
   kernel.
5. [`literature-discrepancies.md`](literature-discrepancies.md) records
   every difference from a published formula, with references and evidence.
6. [`cubic-kernel-conventions.md`](cubic-kernel-conventions.md) writes the
   curved-space cubic kernels explicitly and fixes the nonlinear distinction
   between velocity and rescaled canonical momentum.

## Release state

The release history is now explicit and internally consistent:

- `v0.1.0` points to commit `1caa3a8`, the original explicit release-preparation
  commit whose contents match the `0.1.0` changelog;
- the completed two-field milestone is versioned as `0.2.0` and has its own
  dated changelog section;
- GitHub Actions watches pushes to `main` and the first corrected run passed;
- formatting, Clippy, all tests, all-target release build, and
  `cargo package --locked` pass for the release contents.

The repository tags the verified source state, but does not publish the crate
to crates.io automatically.  Publishing remains a separate explicit action.

## Priorities

Work proceeds in this order.  A lower priority does not start while a higher
one has an unmet acceptance condition.

### P0: restore repository integrity

1. **Completed:** change the CI push branch from `master` to `main` and verify
   a run on the default branch.
2. **Completed:** run formatting, Clippy, tests, release build, and
   `cargo package --locked`.
3. **Completed:** restore `v0.1.0` at its documented release-preparation
   commit and version the two-field milestone as `v0.2.0`.

Acceptance met: `main` triggers green CI, the package contents are verified,
and the version/tag history is internally consistent.

### P1: make the scientific output self-explanatory

1. **Completed:** export both horizon-exit and requested-output spectra, with
   `N`, `k/(aH)`, and explicit evaluation-epoch metadata.
2. **Completed:** make `k/k_pivot` the canonical exported coordinate and label
   the normalization-dependent raw value as `k_internal_over_m_pl`.
3. **Completed:** add a typed observational/likelihood boundary naming Planck
   2018 legacy PR3, the baseline `TT,TE,EE+lowE` components, likelihood release
   R3.00, and the physical `0.05 Mpc^-1` pivot.  The current `A_s = 2.1e-9` is
   explicitly labeled as a benchmark target, not a fit.  See
   [`planck-2018-likelihood-boundary.md`](planck-2018-likelihood-boundary.md).

Acceptance met: a spectrum row states its evaluation epoch, horizon-exit to final
transfer is directly inspectable, and observational assumptions are explicit.

### P2: remove avoidable repeated work

**Completed:** store background checkpoints once and interpolate them for
independent mode evolution.  The former driver that re-integrated the
background for every wavenumber remains as a test-only numerical reference.

The implementation contract is fixed in
[`background-checkpoint-interpolation.md`](background-checkpoint-interpolation.md):
accepted survey states are interpolated with cubic Hermite segments in `N`,
while the current coupled driver remains the numerical reference in tests.

Acceptance met: the optimized and reference drivers agree within the declared
regression tolerance, including a held-out background point in the sharpest
turn and the complete 41-point benchmark.

### P3: extend the physics

- cubic interactions and non-Gaussianity.  The architecture decision is
  complete: use full tree-level phase-space transport, validate its two-point
  sector against the existing mode solver, then add triangle-dependent
  three-point evolution and nonlinear curvature projection.  The ordered
  P3a--P3e milestones are specified in
  [`non-gaussianity-architecture.md`](non-gaussianity-architecture.md).
  **P3a is completed:** validated canonical momentum triangles, bounded tagged
  phase-space indices, and dimensioned two-/three-point blocks are implemented.
  **P3b is completed:** third-order automatic differentiation, dimensioned
  coordinate $V_{,IJK}$, symmetrized covariant $V_{;(IJK)}$, explicit connection
  derivatives, and typed $\nabla_MR^I{}_{JKL}$ are implemented with flat and
  hyperbolic analytic tests. **P3c is completed:** the phase-space convention
  and all-block mode-solver oracle are implemented. Direct covariance evolution
  agrees with the mode solver at subhorizon, horizon-exit, maximum-turn, and
  final epochs, and passes adaptive-tolerance and initialization-depth
  convergence tests. **P3d is completed:** its first slice adds the complete
  dimensioned `QQQ`--`PPP` correlation state, permutation-safe Fourier-leg
  operations, and integration arithmetic. Homogeneous $u^a{}_b$ evolution now
  propagates all three canonical triangle legs on the shared P2 background and
  commutes with equal-leg permutations. The dimensioned $u^a{}_{bc}$ container and
  ordered three-leg $u\Sigma\Sigma$ contraction are implemented. The first
  model-derived cubic coefficients $A_{IJK}$, $B_{IJK}$, and $C_{IJK}$ are now
  evaluated from the background, field-space geometry, and canonical momentum
  triangle with typed $M^1/M^0/M^{-1}$ arithmetic and simultaneous index/leg
  symmetrization. The $A/B$ regressions isolate their Riemann contributions;
  $A$ additionally checks the symmetrized third potential derivative and the
  rank-five $\nabla R$ index routing. A single public entry point constructs all
  five raised tensors and feeds the eight-block e-fold kernel. A second public
  entry point constructs all three output-leg orderings required by the cyclic
  source, with input legs $(2,3)$, $(1,3)$, and $(1,2)$. The quantum initial
  condition and the sourced adaptive evolution are implemented.
  The source derivation has exposed a required semantic decision: the
  published cubic-Hamiltonian kernel evolves rescaled canonical momentum
  $\Pi^I$, while the accepted P3c type means $P^I=D_tQ^I$. They coincide only
  at linear order. A separate `CanonicalMomentumPerturbation` marker and
  canonical covariance, correlation, and kernel aliases now preserve this
  distinction at compile time; a compile-fail test rejects mixing a canonical
  kernel with velocity covariances. A typed five-tensor container now also
  distinguishes first- and third-index raising for $B$ and $C$, and maps the
  cosmic-time $A/B/C$ coefficients into all eight e-fold-time kernel blocks.
  Explicit linear two-point conversion relabels $V^I=\Pi^I$ only at linear
  order. The tree-level three-point conversion implements
  $V^I=\Pi^I-\tfrac12B_{JK}{}^IQ^JQ^K-C^I{}_{JK}\Pi^JQ^K$ and adds the
  corresponding ordered $B\Sigma\Sigma$ and $C\Sigma\Sigma$ terms. Typed tests
  cover the two-point round trip, preservation of `QQQ`, and exact scalar
  nonlinear corrections.
  A Fourier-normalization audit corrected the stripped three-point dimensions
  to $M^{-3}$, $M^{-2}$, $M^{-1}$, and $M^0$ for zero through three momentum
  slots. Accordingly $A$, $B$, and $C$ carry $M^1$, $M^0$, and $M^{-1}$.
  The former half-integer assignment had incorrectly treated a stripped
  three-point function as a simple product of three linear mode functions.
  Canonical triangle helpers now supply the dimensioned all-incoming dot
  products, cosines, and transverse factors used repeatedly by the cubic
  coefficients, with general, equilateral, and folded tests.
  The tree-level quantum initial condition for all eight canonical blocks is
  implemented on a common subhorizon slice and checked against an independent
  in-in derivation, the published DFMS result, and the de Sitter transport
  equation. That check corrected the operator order of the source
  covariances and the momentum arguments of the input-exchanged kernel
  blocks. The combined sourced adaptive system is implemented and converges in start
  depth and tolerance on the reference benchmark, after adding the e-fold
  rates of the cubic tensors to the initial condition. A field-redefinition test
  in sheared flat coordinates confirms that the initial condition and the
  evolved state transform as tensors. The first P3e slice projects the
  result to $P_\zeta$, $B_\zeta$, and $f_{\mathrm{NL}}$ and validates it
  against the mode solver, Maldacena's single-field result, and the squeezed
  consistency relation; see
  [`curvature-projection.md`](curvature-projection.md). See
  [`two-point-transport-oracle.md`](two-point-transport-oracle.md) and
  [`three-point-transport.md`](three-point-transport.md), with the detailed
  derivation in
  [`cubic-kernel-conventions.md`](cubic-kernel-conventions.md);
- after the original P3d/P3e path is complete, the first external model
  validation target is [*Nonminimally Assisted Inflation: A General
  Analysis*](https://arxiv.org/abs/2302.05866).
  Its Einstein-frame two-field metric and potential will be implemented as
  ordinary `FieldSpaceGeometry` and potential providers, not as exceptions in
  the solver. Acceptance requires the existing dimension types, index
  variance, Riemann identities, and covariant perturbation tests to remain
  active. The loop-inflation Class-I trajectory is the first reproduction
  target; Class-II power-law and hybrid cases additionally require a general
  end-of-inflation event;
- reheating and post-inflationary transfer;
- dynamically sized field spaces if a concrete model requires more than the
  current const-generic representation.

These are new scientific milestones, not unfinished parts of the linear
two-field solver.

## Scope boundaries

The current benchmark demonstrates numerical and physical consistency.  It is
not a parameter inference, a claim that double-quadratic inflation is favored
by data, or a general reheating prediction.  The P1 boundary makes that scope
machine-readable; those statements still require a Boltzmann/likelihood layer
and later physics extensions.
