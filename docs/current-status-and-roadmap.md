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
complete.  The current `main` reference is commit `041f2bf`.  At that commit,
133 unit tests and 17 compile-fail cases pass, as do formatting and Clippy.
The conformal-time and former cosmic-time formulations agree in the final
curvature spectrum to a maximum relative difference of `3.1e-10` over the
41-point benchmark grid.

## Release state

The code is scientifically usable as a reference implementation, but the
repository is not yet in a coherent release state:

- `Cargo.toml` still declares version `0.1.0`;
- the completed two-field work remains under `Unreleased` in the changelog;
- neither the local repository nor GitHub has a release tag;
- GitHub Actions watches pushes to the former `master` branch, so pushes to
  `main` do not currently run CI.

The absence of a tag matters because the changelog currently presents
`0.1.0` as released even though there is no corresponding Git reference.
Release preparation must resolve that history explicitly instead of silently
inventing a tag.

## Priorities

Work proceeds in this order.  A lower priority does not start while a higher
one has an unmet acceptance condition.

### P0: restore repository integrity

1. Change the CI push branch from `master` to `main` and verify a run on the
   default branch.
2. Run the complete release checks: formatting, Clippy, tests, release build,
   and `cargo package --locked`.
3. Prepare the next release version and changelog.  Before publishing or
   pushing a tag, explicitly decide whether `v0.1.0` is a historical release
   to backfill or whether the next version is the first tagged release.

Acceptance: `main` triggers green CI, the package contents are verified, and
the version/tag history is internally consistent.

### P1: make the scientific output self-explanatory

1. Export both horizon-exit and final spectra instead of requiring transfer
   to be reconstructed from mode-evolution rows.
2. Normalize the scale factor at the end of inflation, or export normalized
   wavenumbers as the canonical values, so a pivot is not presented as an
   unreadable raw value near `1e9 M_pl`.
3. Add an observational/likelihood boundary that names the pivot convention
   and data release.  The current `A_s = 2.1e-9` remains a benchmark target,
   not a fit.

Acceptance: a spectrum row states its evaluation epoch, horizon-exit to final
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
by data, or a general reheating prediction.  Those statements require the P1
likelihood boundary and later physics extensions.

