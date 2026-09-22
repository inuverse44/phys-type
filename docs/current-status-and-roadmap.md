# Current status and priority roadmap

This document is the single entry point for deciding what to do next.  The
detailed mathematical contracts remain in the design documents linked below;
this page records what is implemented, what is only designed, and the order in
which unfinished work should proceed.

## Status on 2026-09-22

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

Store background checkpoints once and interpolate them for independent mode
evolution.  The current driver deliberately re-integrates the background for
every wavenumber because that was the simplest correct reference path.

Acceptance: the optimized and reference drivers agree within the declared
integration tolerance, including through the sharpest turn in the benchmark.

### P3: extend the physics

- cubic interactions and non-Gaussianity, following the extension boundary in
  [`field-space-curvature-design.md`](field-space-curvature-design.md);
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
