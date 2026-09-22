# Background checkpoint interpolation

## Goal

The reference spectrum driver re-integrates the same background once to find
the pivot, once to find each mode start, and again together with every mode.
That is deliberately simple but makes the cost grow with the number of
wavenumbers.  P2 stores one accepted background trajectory and lets every
mode read it independently.

This is an optimization only.  It must not change the public spectrum
conventions, event definitions, tolerances, or CSV schema.

## Coordinate and checkpoints

The interpolation coordinate is elapsed e-folding number

```text
N = ln(a / a_initial).
```

Unlike a global conformal-time coordinate, `N` stays well resolved when the
scale factor becomes exponentially large.  The existing cosmic-time
background survey remains the source of truth.  Every accepted survey state,
not only every plotted sample, becomes a checkpoint containing

```text
(N, phi^I, pi^I, d phi^I/dN, d pi^I/dN).
```

The derivatives follow directly from the background equations:

```text
d phi^I/dN = pi^I/H,
d pi^I/dN = (D_t pi^I in coordinate components)/H.
```

Between adjacent checkpoints, each field and velocity component uses cubic
Hermite interpolation.  `ln a` is exact by construction,
`ln a = ln(a_initial) + N`; the ordinary scale factor is reconstructed only at
the boundary where a physical snapshot is needed.

Pivot, Bunch-Davies start, horizon exit, and requested output are all located
against the same stored trajectory.  Event roots are refined within the
bracketing checkpoint interval.

## Independent mode system

The mode state contains only `v = a U` and `D_eta v`.  With `N` as the
independent variable, the conformal-time equation is evaluated using the
interpolated background and `d/dN = (aH)^-1 d/deta`.  The field-space
connection terms are retained in the coordinate derivatives, so this is the
same covariant equation as the conformal-time reference driver.

The adaptive error norm contains only mode components.  Public step settings
remain local cosmic-time scales; at each attempt they are converted through
`dN = H dt`.  A step is capped exactly at horizon exit or the requested output
instead of overshooting and reconstructing an event from sampled rows.

## Ownership and parallelism

The stored trajectory is immutable after the survey.  Each wavenumber owns
its mode state and reads the same trajectory, so mode runs are independent.
This milestone removes repeated background work but keeps execution serial;
parallel scheduling can be added later without changing the numerical API.

## Validation

The former background-plus-mode driver remains available to unit tests as the
reference path.  Acceptance requires:

1. interpolated checkpoint endpoints reproduce the accepted background
   states exactly, within floating-point roundoff;
2. interpolation error is below the declared integration tolerance on held-out
   background points, including the benchmark's sharpest turn;
3. optimized and reference horizon-exit and final spectrum rows agree within
   the declared spectrum/integration tolerance on shared wavenumbers;
4. the optimized path performs one background survey regardless of `k_count`;
5. the full 41-point benchmark and all physical validity checks still pass.

The unit-test regression threshold is `1e-8` relative for positive spectrum
components and `1e-8` absolute for signed cross spectra.  It is deliberately
recorded in the test rather than inferred from plot resolution.

## Verified result

The implementation uses one immutable trajectory for pivot location, mode
initialization, horizon exit, and final output.  On the short coupled-driver
regression, relative differences in `P_R` and `P_S` are at most `2e-10`.

On the 41-point mass-ratio-seven benchmark, comparison with the last coupled
driver output gives:

- final `P_R`: maximum relative difference `4.0e-10`;
- horizon-exit `P_R` and `P_S`: maximum relative difference `2.0e-9`;
- final `r_C`: maximum absolute difference `6.4e-6`;
- final `P_S`: maximum absolute difference `6.1e-67` (its relative difference
  is uninformative after decay to approximately `1e-62`).

The background interpolation also has a held-out comparison against direct
integration in the checkpoint interval containing the sharpest benchmark
turn, in addition to exact checkpoint reproduction tests.
