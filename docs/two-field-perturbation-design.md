# Two-field perturbations and power spectra

## Purpose

This document specifies the path from the existing canonical multifield
background solver to dimension-safe two-field scalar perturbations. The first
scientific deliverable is a reproducible figure showing

- the two coordinate-field fluctuations,
- the instantaneous adiabatic and entropy directions,
- the curvature, isocurvature, and cross spectra, and
- their evolution through horizon exit and any subsequent turn.

The conventions follow Gordon et al. for the adiabatic/entropy decomposition,
Wands et al. for transfer and correlation observables, and the covariant
field-space treatment of Langlois and Renaux-Petel. Lalak et al. is the main
reference for numerically evolving a two-field system beyond the simplest
horizon-crossing estimate.

This is a design contract. Each implementation stage gets its own branch and
must pass its tests before it is merged.

## Terminology

There are two coordinate fields, `phi^1` and `phi^2`, but only one
instantaneous adiabatic direction. It is tangent to the background trajectory.
The one independent direction perpendicular to it is the entropy, or
isocurvature, direction. Thus

```text
(Q^1, Q^2)  --field-space projection-->  (Q_sigma, Q_s)
```

The coordinate-basis spectra `P_11` and `P_22` are useful diagnostics and will
be retained. The primary physical outputs are `P_R`, `P_S`, and `C_RS`.

## Physical scope

The target action is

```math
S = \int d^4x\,\sqrt{-g}\left[
    \frac{M_{\rm pl}^2}{2}R
    -\frac12 G_{IJ}(\phi)g^{\mu\nu}\partial_\mu\phi^I\partial_\nu\phi^J
    -V(\phi)
\right].
```

The first end-to-end release uses the canonical flat metric
`G_IJ = delta_IJ`. The public tensor interfaces nevertheless distinguish
upper and lower field-space indices from the beginning. A later curved-metric
implementation can therefore replace the metric provider without changing
the meaning of vectors, covectors, projections, or spectra.

Initially out of scope are non-canonical kinetic terms, more than two fields
in the kinematic basis, nonlinear perturbations, reheating transfer to named
particle species, and CMB Boltzmann evolution.

## Typed field-space algebra

Physical mass dimension and tensor variance are independent type-level
properties:

```rust,ignore
FieldVector<Upper, D, V, N>
FieldVector<Lower, D, V, N>
FieldTensor2<Upper, Lower, D, V, N>
```

`D` is the existing mass dimension, `V` is the stored numerical value, and
`N` is the const-generic field count. `Upper` and `Lower` are zero-sized marker
types. User-facing aliases include `ContravariantField`, `CovariantField`,
`FieldMetric`, and `InverseFieldMetric`.

Only mathematically meaningful operations are exposed:

- `G_IJ v^J -> v_I` lowers an index.
- `G^IJ w_J -> w^I` raises an index.
- `w_I v^I -> scalar` contracts opposite variance.
- Same-variance contraction is a compile-time error.
- Addition requires equal variance and equal physical dimension.
- Raising and lowering preserve physical dimension.

Compile-fail tests are part of this API contract. In particular, `v^I w^I`,
using `G_IJ` to raise an index, and adding `v^I` to `w_I` must not compile.

The flat metric is an explicit value rather than an implicit dot product. This
keeps scientific code honest about which coordinate system it inhabits.

## Covariant background

Define

```math
\dot\sigma^2 = G_{IJ}\dot\phi^I\dot\phi^J,
\qquad
H^2 = \frac{\dot\sigma^2/2 + V}{3M_{\rm pl}^2},
\qquad
\epsilon_H = \frac{\dot\sigma^2}{2M_{\rm pl}^2H^2}.
```

The background equations are

```math
\mathcal D_t\dot\phi^I + 3H\dot\phi^I + G^{IJ}V_{;J}=0,
```

where

```math
\mathcal D_t A^I = \dot A^I + \Gamma^I{}_{JK}\dot\phi^J A^K.
```

For the initial flat implementation `Gamma^I_JK = 0`, but the derivative is
still represented by a covariant interface. Inflation ends at
`epsilon_H = 1`.

## Kinematic basis

The tangent unit vector and its dual are

```math
e_\sigma^I = \frac{\dot\phi^I}{\dot\sigma},
\qquad
e_{\sigma I}=G_{IJ}e_\sigma^J.
```

For two fields, `e_s^I` is the normalized vector satisfying

```math
G_{IJ}e_s^I e_s^J=1,
\qquad
G_{IJ}e_s^I e_\sigma^J=0.
```

Its sign is fixed continuously along the numerical trajectory: choose the
positive coordinate orientation at the initial point, then flip a newly
computed normal when its metric inner product with the previous normal is
negative. This prevents artificial sign changes in `Q_s` and `C_RS`.

The turn vector and dimensionless turn rate are

```math
\omega^I = \mathcal D_t e_\sigma^I,
\qquad
\eta_\perp = \frac{e_{sI}\omega^I}{H}.
```

The perturbations in the two bases are

```math
Q_\sigma=e_{\sigma I}Q^I,
\qquad
Q_s=e_{sI}Q^I.
```

## Linear perturbation equation

`Q^I` denotes the covariant field fluctuation in spatially flat gauge. For
each comoving wave number `k`, the canonical system obeys

```math
\mathcal D_t^2 Q^I + 3H\mathcal D_t Q^I
+ \left(\frac{k^2}{a^2}\delta^I{}_J + C^I{}_J\right)Q^J=0,
```

with

```math
C^I{}_J = G^{IK}V_{;KJ}
-\mathcal R^I{}_{LMJ}\dot\phi^L\dot\phi^M
-\frac{1}{M_{\rm pl}^2a^3}\mathcal D_t
 \left(\frac{a^3}{H}\dot\phi^I\dot\phi_J\right).
```

This equation also fixes the Riemann-sign convention used by the code. The
flat first release has a vanishing Riemann term. The gravitational term is
not to be omitted even in flat field space.

The potential needs a Hessian with physical dimension `M^2`. A second-order
jet value

```rust,ignore
Jet<N> { value, gradient: [f64; N], hessian: [[f64; N]; N] }
```

will evaluate `V`, `V_,I`, and `V_,IJ` from one potential formula. For curved
field space the coordinate Hessian is converted to
`V_;IJ = V_,IJ - Gamma^K_IJ V_,K`.

## Quantum initial state and mode matrix

A correlated two-field spectrum cannot be obtained by evolving one classical
pair `(Q_sigma, Q_s)`. There are two independent vacuum oscillators. The
solver therefore evolves a complex `2 x 2` mode matrix `U^I_A`, where `A`
labels the independent Bunch-Davies initial modes, together with its covariant
time derivative.

Complex values use `num_complex::Complex64` inside the existing dimensional
wrapper:

```rust,ignore
type ModeAmplitude = Quantity<MMinusHalf, Complex64>;
type ModeVelocity = Quantity<MHalf, Complex64>;
```

Physical dimension therefore remains a compile-time property of the outer
`Quantity`; using a complex numerical representation does not erase it.
Field-space variance likewise remains in `FieldVector` or tensor types and is
independent of complex conjugation.

Initial conditions are imposed in an orthonormal frame a fixed number of
e-folds before each mode leaves the horizon. The depth `x_init = k/(aH)` is
the exponential of that number, and the default is five e-folds, so
`x_init ~= 148`:

```math
U^a{}_A \simeq \frac{\delta^a{}_A}{a\sqrt{2k}}e^{-ik\tau},
\qquad
\mathcal D_t U^a{}_A \simeq
\left(-H-i\frac{k}{a}\right)U^a{}_A.
```

Only the relative phase matters for spectra, so the implementation may set
the common initial phase to zero. The normalization is checked by a conserved
Wronskian. Runs that cannot place every requested mode sufficiently inside
the horizon return a typed error instead of silently using a poor vacuum
approximation. The same applies at the other end: a mode that has not left
the horizon by the requested output time still holds its vacuum, and
reporting that as a spectrum would look entirely physical — positive,
Hermitian, and orders of magnitude too large. The driver checks `k/(aH) < 1`
at the output time and refuses otherwise.

The first correctness-oriented release evolves `Complex64` mode components
and the background together for each `k` through the existing integrator.
Because all equation coefficients are real, the numerical result is
equivalent to evolving real and imaginary parts separately, while the public
API can express conjugation, Hermitian products, and phase operations
directly. A later optimization may generate background checkpoints first and
use monotone cubic Hermite interpolation, avoiding repeated background work
without changing the mode API.

Step-size control is a property of the integrator rather than of any one
system. A system implements `AdaptiveSystem` to reduce an embedded error
estimate to a single dimensionless ratio, and the shared controller decides
acceptance and the next step size from that ratio alone. The single-field
background, both multifield backgrounds, and the two-field perturbation system
use the same controller.

Adaptive error control maps a dimensioned complex component to a real
normalized error using its modulus:

```math
E(z)=\frac{|\delta z|}{a_{\rm tol}+r_{\rm tol}|z|}.
```

Absolute tolerances retain the component's physical dimension. Relative
tolerances and the final combined norm are dimensionless. The integrator must
not compare only real parts or apply separate acceptance decisions to real and
imaginary components.

Complex numbers are split into real and imaginary columns only at explicit
interoperability boundaries such as CSV. Both columns are emitted, even when
a particular observable is expected to be real.

## Spectra and normalization

For the mode matrix, define

```math
\Sigma^{IJ}(k)=\sum_A U^I{}_A U^{J*}{}_A,
\qquad
\mathcal P^{IJ}_Q(k)=\frac{k^3}{2\pi^2}\Sigma^{IJ}(k).
```

`Sigma` is constructed with an explicit complex conjugate and checked for
Hermiticity before projection. A physically real projected spectrum may
discard an imaginary residual only after verifying that it is below a
documented tolerance.

Project this matrix with `e_sigma` and `e_s` at the requested output time:

```math
\mathcal P_{Q_\sigma}=e_{\sigma I}e_{\sigma J}\mathcal P_Q^{IJ},
\quad
\mathcal P_{Q_s}=e_{sI}e_{sJ}\mathcal P_Q^{IJ},
\quad
\mathcal C_{Q_\sigma Q_s}=e_{\sigma I}e_{sJ}\mathcal P_Q^{IJ}.
```

The dimensionless curvature and entropy variables use the same normalization,

```math
\mathcal R=\frac{H}{\dot\sigma}Q_\sigma,
\qquad
\mathcal S=\frac{H}{\dot\sigma}Q_s.
```

Consequently the reported spectra are

```math
\mathcal P_{\mathcal R}=\left(\frac{H}{\dot\sigma}\right)^2
\mathcal P_{Q_\sigma},
\quad
\mathcal P_{\mathcal S}=\left(\frac{H}{\dot\sigma}\right)^2
\mathcal P_{Q_s},
\quad
\mathcal C_{\mathcal RS}=\left(\frac{H}{\dot\sigma}\right)^2
\mathcal C_{Q_\sigma Q_s}.
```

The correlation coefficient is

```math
r_\mathcal C =
\frac{\mathcal C_{\mathcal RS}}
{\sqrt{\mathcal P_{\mathcal R}\mathcal P_{\mathcal S}}},
```

which is the `cos Delta` convention of Wands et al. Values are clamped to
`[-1, 1]` only after checking that any overshoot is consistent with floating
point tolerance.

The scale factor has an arbitrary normalization. Configuration therefore
defines a pivot `k_pivot` by `k_pivot = a_* H_*` at a selected number of
e-folds before the end of inflation, and reports both `k/k_pivot` and the
internal comoving `k`.

## Data products and figures

The Rust example writes machine-readable CSV with metadata comments recording
model parameters, tolerances, `x_init`, pivot definition, and the Git revision.
At minimum it emits:

```text
background.csv:
N, phi_1, phi_2, pi_1, pi_2, H, epsilon_H,
e_sigma_1, e_sigma_2, e_s_1, e_s_2, eta_perp

spectra.csv:
k_over_pivot, k, P_11, P_22, Re_P_12,
P_Qsigma, P_Qs, C_Qsigma_Qs, P_R, P_S, C_RS, r_C

mode_evolution.csv:
k_over_pivot, N, k_over_aH, P_R, P_S, C_RS, eta_perp
```

The principal two-field figure has four panels:

1. the trajectory in the `phi_1`--`phi_2` plane, with tangent and normal
   arrows at sparse points;
2. `P_11`, `P_22`, `P_Qsigma`, and `P_Qs` versus `k/k_pivot`;
3. `P_R`, `P_S`, and signed `C_RS` versus `k/k_pivot`;
4. `r_C` and the turn rate, or selected-mode evolution versus e-folds.

Cross spectra are signed and therefore use a symmetric-log axis when needed.
The plotting code must never take `abs(C_RS)` without also retaining its sign.

## Numerical and physical validation

Every implementation stage adds tests appropriate to its layer:

- Metric raise/lower round trips and orthonormal-basis identities.
- A concrete curved geometry, so that the covariant machinery is exercised
  with a non-vanishing connection and Riemann tensor rather than only through
  `FlatFieldSpace`, where every curvature contribution is identically zero.
  `HyperbolicFieldSpace` is checked against the Riemann symmetries and the
  first Bianchi identity, against `R = -2/L^2`, and for metric compatibility
  by finite differences. Its curvature is then shown to reach the perturbation
  mass matrix and the spectra by holding the metric and connection fixed and
  zeroing only the Riemann tensor: the mass matrix moves by `0.4%` and `P_R`
  by `6e-3`. Comparing against `FlatFieldSpace` at equal coordinates would not
  be a controlled comparison, because this metric approaches flat space in
  polar rather than Cartesian coordinates as `L` grows.
- Compile-fail tests for illegal contractions and variance mismatches.
- Potential Hessians checked against analytic quadratic models and finite
  differences.
- Straight, separable trajectories have zero turn rate and negligible cross
  spectrum.
- Deep inside the horizon the two independent modes satisfy the Wronskian and
  approach the Bunch-Davies form.
- Near horizon exit, two light canonical modes satisfy
  `P_Qsigma ~= P_Qs ~= (H_*/2pi)^2` within a documented slow-roll tolerance.
  The comparison needs a window rather than an instant. Before about two
  e-folds after exit the exact mode still carries its `1 + x^2` correction,
  which is the mirror of the `1 + 1/x^2` correction at initialization. After
  that window both spectra drift for physical reasons, at rates that are
  themselves checked: the entropy spectrum decays as a light field of its own
  mass, `P_S ~ a^(-2 m^2 / 3 H^2)`, and the adiabatic field spectrum grows
  with `2 epsilon` because `Q_sigma = (dot sigma / H) R` while `R` is
  conserved. Measured over six e-folds, the decay is `0.71` against a
  predicted `0.70` and the growth is `1.26` against a predicted `1.27`.
- The effectively single-field limit satisfies
  `P_R ~= H_*^2/(8 pi^2 epsilon_* M_pl^2)`.
- Results are stable under `x_init = 50, 100`, tighter RK tolerances, denser
  background checkpoints, and a denser `k` grid.
- The covariance matrix is Hermitian positive semidefinite and
  `|r_C| <= 1` up to numerical tolerance.
- On super-horizon scales, the transfer matrix has the triangular Wands form:
  entropy can source curvature, while a purely adiabatic mode does not create
  entropy. This statement is not applied to sub-horizon evolution. Two
  independent solutions determine the matrix as `T = M_late M_early^-1`, where
  `M = [[R_1, R_2], [S_1, S_2]]` comes from `project_modes`, so no special
  initial condition has to be constructed. The forbidden entry is not exactly
  zero, because the early time is a finite distance outside the horizon:
  measured two e-folds out, `T_SR / T_SS = 8e-5` against gradient corrections
  of order `(k/aH)^2 = 0.018`, while the permitted sourcing `T_RS` is five
  hundred times larger than `T_SR`. The triangle therefore has a direction,
  rather than every off-diagonal entry being small.

Reference curves and tolerances are stored with enough metadata to make a
failure scientifically interpretable rather than merely a changed snapshot.

Two limitations of the first implementation are recorded here rather than
left to be rediscovered. Neither blocks the current results; both shape what
the next change should be.

### The dynamic range of cosmic-time variables

The absolute error tolerance is a single number applied to raw Planck-unit
component values, and a run spans an enormous range of them. In the example
model, inflation lasts `N = 100.8`, so the scale factor reaches `3e43` while
the mode amplitude, carrying `1 / (a sqrt(2k))`, falls to `1e-51`: about
ninety-five orders of magnitude between the largest and smallest evolved
number.

An absolute floor anywhere near the mode amplitudes caps their accuracy no
matter how small the relative tolerance is. The symptom is subtle, because
nothing fails: the spectra simply stop converging as `x_init` grows, settling
on a plateau near `1e-3`.

The implemented remedy is to stop evolving a quantity with that range. The
mode matrix is integrated rescaled,

```math
v^I{}_A = a\,U^I{}_A,
```

which turns

```math
D_t^2 U + 3H D_t U + \left(M + \frac{k^2}{a^2}\right)U = 0
```

into

```math
D_t^2 v + H D_t v
+ \left(M + \frac{k^2}{a^2} - (2-\epsilon_H)H^2\right)v = 0,
```

where the subtracted term is `a''/a` written in cosmic time. The vacuum
initialization becomes `|v| = 1/sqrt(2k)` with `D_t v = -i(k/a)v`, so the
Hubble drag cancels exactly and the sub-horizon solution oscillates at fixed
size instead of decaying like `1/a`. The canonical Wronskian keeps one power
of `a` rather than three.

The measured effect is that the absolute tolerance stops mattering. Measured
with the single global floor that was in use at the time, a value of `1e-30`
destroyed convergence before the rescaling, while afterwards the spectra were
unchanged from `1e-45` up to about `1e-14`: roughly sixteen orders of
headroom. The two formulations agree to ten significant digits, which is a
strong cross-check given that the drag term, the `a''/a` term, the Wronskian
prefactor and the initialization all differ between them.

That global floor no longer exists; the next section replaces it with one
tolerance per group of components. Re-measured against the mode tolerance
alone, the curvature spectrum moves by `3e-9` at a mode floor of `1e-16` and
by `9e-7` at `1e-10`, so the setting has stopped being delicate in either
form.

This is conditioning, not a change to the type-level contract: `a` is
dimensionless, so `v` and `U` share the dimension `M^-1/2` and the type
system cannot separate them. The rescaling is therefore undone in exactly
one place, `mode_covariance`, which takes the scale factor as an argument so
that the conversion cannot be forgotten silently.

The rescaled amplitude still grows with `a` once a mode is outside the
horizon, so the rescaling moves the range into a benign region rather than
eliminating it. What removes the hazard outright is giving each group of
components its own absolute tolerance, in that group's own dimension:

```rust,ignore
pub struct TwoFieldPerturbationTolerance {
    pub background: MultifieldStateDelta<2>,
    pub rescaled_amplitude: ModeAmplitude<f64>,
    pub rescaled_velocity: ModeVelocity<f64>,
}
```

`AdaptiveSystem` carries this as an associated type, so a tolerance meant for
one quantity cannot bound another, and a compile-fail test holds that open.

The mode tolerances default to zero, which is the point rather than an
oversight. A mode error is normalized by a modulus, and the modulus of a
physical mode does not vanish, so no floor is needed to keep the ratio finite.
That argument is about the mode, not about each component of the mode matrix,
and the norm follows it: every component is measured against the largest
modulus in its matrix. An individual component can sit exactly at zero — an
axis-aligned trajectory puts exact zeros off the diagonal — and normalizing
such a component by itself would give it an unbounded ratio as soon as
rounding made its error non-zero. Any floor placed there becomes the accuracy ceiling of
the modes instead, which is exactly how the run used to lose its convergence
in `x_init` without anything failing. The background components are real and
can pass through zero, so they keep a floor, and one raw Planck-unit number
is honest for them because they stay within a few orders of each other.

One thing is deliberately left undone. The independent variable is still
cosmic time: conformal time would additionally make the sub-horizon frequency
constant at `k`, but `Time` and a conformal time would share the dimension
`M^-1` and so would be the same type, which needs its own marker to stay
safe.

Relatedly, normalizing `a = 1` at the initial state rather than at the end of
inflation is what makes `k_pivot` print as `5.5e14 M_pl`. The number is
correct but unreadable, and an end-of-inflation normalization would keep
reported wavenumbers close to unity.

### The Bunch-Davies approximation

Bunch-Davies initialization is leading order in `1/x_init`, so a spectrum
carries a systematic error that falls off as `1/x_init^2`. This is a genuine
approximation, not a numerical artifact, and is fixed by a convergence test
rather than by an absolute threshold.

Because the start is expressed in e-folds before horizon exit, the trade-off
is regular and easy to price. Each additional e-fold multiplies the work by
about `e`, since that many more sub-horizon oscillations must be resolved, and
divides the error by about `e^2`. Measured on the two-field quadratic model:

```text
e-folds before exit   x_init   relative error   time
        1                2.7        1.5e-1      0.002 s
        3               20          3.1e-4      0.003 s
        4               55          2.7e-4      0.005 s
        4.6            100          9.9e-5      0.007 s
        5              148          2.1e-5      0.010 s
        6              403          5.8e-6      0.024 s
        7             1097          3.8e-7      0.061 s
```

The first row is not an integration failure. The exact de Sitter mode has
`|u|^2 = (1 + 1/x^2)/(2 k a^2)`, and the initialization keeps only the leading
term, so at `x = 2.7` the discarded correction is `13.5%` against a measured
`14.8%`. One e-fold before exit is simply not yet Bunch-Davies.

The crossing that starts a mode, and the one the spectra are read from, are
refined by bisection inside the last accepted step using the same fifth-order
method the step was accepted on. Refining with a fourth-order sub-step instead
moved the reported spectra by `2e-5` to `5e-5`, which is the same size as the
systematic the table below is used to control, so the cheaper choice was not
worth the accuracy.

From four e-folds onward the measured error tracks `1/x_init^2` closely. The
default of five e-folds buys a systematic error around `2e-5`, well below the
other approximations in the pipeline, for a few milliseconds per mode. There
is no optimum to find: error falls as the square of the cost, so the number is
chosen from the accuracy wanted, not from a sweet spot.

## Branch sequence

Work proceeds in this order, one independently reviewable feature at a time:

1. `docs/two-field-perturbation-design` -- this contract and references.
2. `feature/field-space-tensors` -- variance markers, metric operations, and
   compile-fail tests.
3. `feature/covariant-background` -- metric provider, connection, kinematic
   basis, continuous normal orientation, and turn rate.
4. `feature/potential-hessian` -- second-order jets and covariant Hessian.
5. `feature/higher-rank-field-tensors` and `feature/field-space-curvature` --
   rank-three/rank-four tensors and the Riemann geometry contract described in
   [`field-space-curvature-design.md`](field-space-curvature-design.md).
6. `feature/two-field-perturbations` -- complex mode matrix, vacuum initial
   state, Wronskian, and adaptive evolution.
7. `docs/complex-mode-design` and `feature/complex-mode-values` -- use
   `Quantity<D, Complex64>` directly and define modulus-based error handling.
8. `feature/power-spectra` -- mode covariance, projections, reality and range
   validation of the projected spectra.
9. `feature/adaptive-system` -- an `AdaptiveSystem` norm and a reusable
   step-size controller shared by every system, including the complex modes.
10. `feature/spectrum-run` -- pivot convention, the `k` grid driver, CSV
    output, and physical validation against the analytic limits.
11. `feature/two-field-spectrum-plots` -- the four-panel scientific figure and
    reproducible generation command.
12. `feature/rescaled-mode-variable` -- integrate `v^I = a U^I` so the mode
    amplitudes stop decaying with the scale factor.
13. `feature/per-component-tolerances` -- a dimensioned absolute scale per
    component group, removing the remaining dependence on one global floor.
14. `feature/conformal-time-modes` -- conformal time as a distinct
    independent-variable type, and a background evolved in `ln a`.

Each branch starts from the updated default branch, contains one conceptual
change, passes `cargo fmt --check`, `cargo clippy --all-targets`, and
`cargo test`, and is pushed before being fast-forwarded into the default
branch.

## References

- C. Gordon, D. Wands, B. A. Bassett, and R. Maartens,
  [Adiabatic and entropy perturbations from inflation](https://arxiv.org/abs/astro-ph/0009131),
  Phys. Rev. D 63, 023506 (2001).
- D. Wands, N. Bartolo, S. Matarrese, and A. Riotto,
  [An observational test of two-field inflation](https://arxiv.org/abs/astro-ph/0205253),
  Phys. Rev. D 66, 043520 (2002).
- Z. Lalak, D. Langlois, S. Pokorski, and K. Turzynski,
  [Curvature and isocurvature perturbations in two-field inflation](https://arxiv.org/abs/0704.0212),
  JCAP 07 (2007) 014.
- D. Langlois and S. Renaux-Petel,
  [Perturbations in generalized multi-field inflation](https://arxiv.org/abs/0801.1085),
  JCAP 04 (2008) 017.

Kaloian Lozanov's reheating lectures are useful for the later transition from
inflationary entropy modes to post-inflationary dynamics, but they are not the
primary convention source for the inflationary mode solver specified here.
