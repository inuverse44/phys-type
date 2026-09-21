# Two-field observational benchmark

## Purpose

The original `two_field_spectra` example used unequal but otherwise arbitrary
masses.  It demonstrated the mode solver, but it was neither normalized to the
observed scalar amplitude nor arranged to expose the characteristic
heavy-to-light transfer in the pivot mode.  This benchmark gives the example a
reproducible physical meaning without presenting it as a parameter fit.

The target is a minimal canonical double-quadratic model,

```math
V(\phi_1,\phi_2)
=\frac12m_1^2\phi_1^2+\frac12m_2^2\phi_2^2,
```

with a flat field-space metric.  It is a numerical benchmark for adiabatic and
entropy transfer, not a claim that the model is preferred by current data.

## Shape parameters

The dimensionless shape is fixed by

```text
m_2 / m_1 = 7
phi_initial / M_pl = (14, 12)
N_pivot = N_end - 55
```

A background survey places the largest sampled `|eta_perp|` about `42.4`
e-folds before the end of inflation, roughly 12.6 e-folds after the pivot mode
exits.  This timing makes the pivot mode display the familiar double-inflation
evolution: the entropy spectrum develops a trough and partial rebound near the
turn while the curvature spectrum is sourced onto a higher plateau.

This statement concerns the time evolution of one mode.  It must not be
confused with a dip in the final spectrum as a function of `k`.  A localized
feature in the final `P_R(k)` requires the transition scale itself to lie in
the plotted `k` window, and is not a generic consequence of having two fields.
Moving the mass ratio or initial direction moves or removes both effects.

The spectrum grid covers `0.01 <= k/k_pivot <= 100` with 41 logarithmically
spaced points.  The wider grid is intentional: it resolves both sides of the
observable pivot neighborhood.  A separate figure follows the pivot mode in
time and is the appropriate comparison to the standard mass-ratio-seven
double-inflation evolution plot.

Mode evolution is exported every five accepted integration steps.  This is an
output-sampling choice, not an integrator step size: it resolves the sharp turn
in the plotted polyline without changing the adaptive error control or the
final spectrum.

## Amplitude normalization

The shape parameters and the overall energy scale have separate jobs.  The
mass ratio and initial fields determine the background trajectory in e-fold
time.  A common rescaling

```math
(m_1,m_2)\longrightarrow\lambda(m_1,m_2)
```

leaves that dimensionless trajectory unchanged and, for this quadratic
model, scales the scalar spectrum approximately as

```math
\mathcal P_\mathcal R\longrightarrow
\lambda^2\mathcal P_\mathcal R.
```

The example first performs a one-mode pilot run at `k = k_pivot`, starting
from `m_1 = 10^-6 M_pl`.  It then computes

```math
\lambda = \sqrt{\frac{A_s}{\mathcal P_\mathcal R^{\rm pilot}(k_*)}},
\qquad A_s=2.1\times10^{-9},
```

and repeats the complete spectrum run with both masses rescaled.  The pilot
value, target, scale factor, and calibrated masses are written into every CSV
metadata block.  The final pivot amplitude is reported separately; the target
is never substituted for a calculated spectrum value.

`A_s` here is a declared benchmark target.  A future likelihood layer should
make the pivot convention and observational data release explicit and should
fit the background parameters instead of treating this normalization as a
fit.

## Interpretation of the outputs

`P_R` in the final CSV is already the curvature spectrum after entropy-to-
curvature transfer through the requested output time.  It is not obtained by
adding `P_R` and `P_S`.  The useful comparison is instead

```text
P_R at horizon exit -> transfer during the turn -> final P_R
```

with `P_S`, `C_RS`, and `eta_perp` diagnosing the transfer.  A later output
extension should retain both horizon-exit and final spectra explicitly so that
the origin of a dip or enhancement can be read without reconstructing it from
the mode-evolution rows.

## Acceptance checks

The tracked reference run must satisfy all of the following:

- the pivot mode exits before the turn and shows entropy decay followed by
  turn-induced transfer into the curvature spectrum;
- the computed final `P_R(k_pivot)` is within one percent of `2.1e-9`;
- all spectrum rows pass the existing Hermiticity, reality, positivity, and
  correlation-range checks;
- the CSV records the shape parameters and amplitude calibration inputs;
- PNG and PDF figures are regenerated from the same CSV.

## Branch

This work is isolated on `feature/two-field-observational-benchmark`.  It is
merged only after formatting, lint, unit and compile-fail tests, the reference
run, and figure generation pass.

The time-domain comparison follows the canonical mass-ratio-seven example in
C. van de Bruck and M. Robinson,
[Power Spectra beyond the Slow Roll Approximation in Theories with
Non-Canonical Kinetic Terms](https://arxiv.org/abs/1404.7806), especially
figure 2.  That paper chooses a reference mode for its evolution plot rather
than identifying the mode with an observational 55-e-fold pivot; this benchmark
states the distinction explicitly.
