# Curvature projection

## Purpose and current boundary

This document describes the first P3e slice: the projection of the canonical
phase-space correlations evolved by
[`three-point-transport.md`](three-point-transport.md) to the curvature
perturbation $\zeta$. The implementation lives in `crates/phys-type/src/curvature.rs`.

Implemented:

- the covariant first- and second-order gauge coefficients $N_a$ and
  $N_{ab}$ on a superhorizon slice;
- the stripped $P_\zeta(k)$ and $B_\zeta(k_1,k_2,k_3)$, including the
  operator-ordered gauge contribution;
- the dimensionless $\mathcal{P}_\zeta$ and a reduced $f_{\mathrm{NL}}$ with a
  stated convention;
- validation against the linear mode solver, superhorizon conservation,
  Maldacena's single-field result, the squeezed consistency relation, a field
  redefinition, the published axion-quartic result of Dias, Frazer, Mulryne,
  and Seery, and a numerical $\delta N$ calculation.

With `run_bispectrum` and its CSV writers (below), every P3e acceptance
condition of [`non-gaussianity-architecture.md`](non-gaussianity-architecture.md)
is met. Within the documented scope (tree level, the canonical transport of
[`three-point-transport.md`](three-point-transport.md), evaluation during
inflation, no reheating), the output is therefore a prediction of the
primordial bispectrum of $\zeta$ for the supplied model. Its limitations are
listed at the end of this document.

## Gauge transformation

Dias, Elliston, Frazer, Mulryne, and Seery,
[The curvature perturbation at second order](https://arxiv.org/abs/1410.3491),
give several equivalent second-order expressions for the uniform-density
curvature perturbation in terms of flat-gauge field perturbations. They are
exact except that terms which decay when every wavenumber is outside the
horizon are dropped. The implementation uses their "simple" representation,

$$
\begin{aligned}
\zeta={}&-\frac{\dot\phi_I Q^I}{2HM_{\mathrm{Pl}}^2\epsilon}
+\frac{1}{6H^2M_{\mathrm{Pl}}^2\epsilon}\biggl(
\frac{\dot\phi_I\dot\phi_J}{M_{\mathrm{Pl}}^2}\Bigl[-\frac32+\frac{9}{2\epsilon}
+\frac{3}{4\epsilon^2}\frac{V_{,K}\dot\phi^K}{M_{\mathrm{Pl}}^2H^3}\Bigr]Q^IQ^J\\
&+\frac{3}{H\epsilon}\frac{\dot\phi_I\dot\phi_J}{M_{\mathrm{Pl}}^2}Q^I\Pi^J
-3H\,\partial^{-2}\bigl[\partial_jQ_I\,\partial_j\Pi^I+\partial^2Q_I\,\Pi^I\bigr]
\biggr),
\end{aligned}
$$

with $\epsilon=\dot\phi_I\dot\phi^I/(2H^2M_{\mathrm{Pl}}^2)$, evaluated on
the output slice. Two properties make this form convenient here.

**No momentum at linear order.** The linear term contains only $Q^I$. The
velocity $V^I$ and the canonical momentum $\Pi^I$ differ only at second order,
so replacing $\delta\dot\phi$ by $\Pi$ in the quadratic terms changes $\zeta$
only at third order. The canonical correlations can therefore be projected
directly, without the tree-level velocity conversion.

**Covariance.** In Riemann normal coordinates at the background point, the
coordinate perturbation is $\delta\phi^I=Q^I+O(Q^3)$ and the connection
vanishes. Curvature then enters only through
$\delta\dot\phi^I=D_tQ^I-\tfrac12\partial_L\Gamma^I{}_{JK}\dot\phi^LQ^JQ^K$,
and $\delta\dot\phi$ appears only in quadratic terms, so this correction is
cubic. Replacing $\delta_{IJ}$ by $G_{IJ}$ therefore gives the covariant
expression. The field-redefinition test below checks this.

**Sign.** At linear order $\zeta=-(H/\dot\sigma)Q_\sigma=-\mathcal R$, where
$\mathcal R$ is the comoving curvature reported by the linear spectra. This is
the sign convention of Maldacena and of the reference above; it matters for the
sign of $B_\zeta$.

### Coefficients

Writing $\zeta(\mathbf k)=N_aX^a+\tfrac12N_{ab}X^aX^b$ with the momentum
convolution suppressed,

$$
N_I=-\frac{\dot\phi_I}{2HM_{\mathrm{Pl}}^2\epsilon},\qquad
N_{Q^IQ^J}=\frac{\dot\phi_I\dot\phi_J}{3H^2M_{\mathrm{Pl}}^4\epsilon}
\Bigl[-\frac32+\frac{9}{2\epsilon}
+\frac{3}{4\epsilon^2}\frac{V_{,K}\dot\phi^K}{M_{\mathrm{Pl}}^2H^3}\Bigr],
$$

$$
N_{Q^I\Pi^J}(\mathbf k;\mathbf k_a,\mathbf k_b)
=\frac{1}{6H^2M_{\mathrm{Pl}}^2\epsilon}\Bigl[
\frac{3\dot\phi_I\dot\phi_J}{H\epsilon M_{\mathrm{Pl}}^2}
-\frac{3H}{k^2}(\mathbf k_a\!\cdot\!\mathbf k_b+k_a^2)G_{IJ}\Bigr],
$$

and $N_{\Pi^IQ^J}(\mathbf k;\mathbf k_a,\mathbf k_b)$ is the same with $k_a^2$
replaced by $k_b^2$; the gradient term attaches the squared momentum to the
field input. $N_{\Pi\Pi}=0$. The dimensions are $M^{-1}$, $M^{-2}$, and
$M^{-3}$ in the stripped Fourier normalization of the transport.

The cross coefficients carry $1/6$, not $1/3$. The symmetric form
$\tfrac12N_{ab}X^aX^b$ counts the $Q\Pi$ and $\Pi Q$ orderings separately, so
each must be half of the real-space coefficient of $Q\Pi$. The matrix printed
in Eq. (7.4b) of Dias, Frazer, Mulryne, and Seery,
[arXiv:1609.00379](https://arxiv.org/abs/1609.00379), carries the prefactor
$1/(3H^2M_{\mathrm{Pl}}^2\epsilon)$ on these entries as well. With that choice,
the straight-trajectory benchmark below gives $f_{\mathrm{NL}}\approx0.73$
against Maldacena's $0.19$, and $\zeta$ is not conserved. With $1/6$, both
checks pass. The quadratic $QQ$ coefficient agrees with that reference. The
discrepancy, its evidence, and the exact references are recorded in
[`literature-discrepancies.md`](literature-discrepancies.md), entry 1.

### Correlators

With $\langle\zeta\zeta\rangle=(2\pi)^3\delta\,P_\zeta$ and
$\langle\zeta\zeta\zeta\rangle=(2\pi)^3\delta\,B_\zeta$,

$$
P_\zeta(k)=N_IN_J\,\mathrm{Re}\,\Sigma_{QQ}^{IJ}(k),
$$

$$
B_\zeta=N_IN_JN_K\alpha_{QQQ}^{IJK}
+\sum_{r}N_a(k_s)N_b(k_t)\,N_{cd}(\mathbf k_r;\mathbf k_s,\mathbf k_t)\,
\Sigma^{ac}(k_s)\Sigma^{bd}(k_t),
$$

where for output leg $r$ the inputs $(s,t)$ are the remaining legs in
ascending order. As in the transport source, a covariance on a leg to the
left of the quadratic insertion keeps the external operator on the left
($\Sigma^{ac}$), and one on the right uses $\Sigma^{ca}$.
`CurvatureCorrelators::bispectrum_imaginary_ratio` reports the imaginary part
of this ordered contraction, which measures residual commutator terms on the
slice; it is below $10^{-6}$ on the benchmarks. $B_\zeta$ is its real part.

The dimensionless power spectrum is $\mathcal P_\zeta=k^3P_\zeta/(2\pi^2)$,
and the reduced bispectrum follows the Maldacena/Planck convention

$$
B_\zeta=\frac65f_{\mathrm{NL}}\bigl[P_\zeta(k_1)P_\zeta(k_2)+P_\zeta(k_1)P_\zeta(k_3)
+P_\zeta(k_2)P_\zeta(k_3)\bigr]
$$

with the sign of $\zeta$ above. For single-field slow roll this gives
$f_{\mathrm{NL}}\to\tfrac{5}{12}(1-n_s)>0$ in the squeezed limit.

`curvature_correlators` rejects a slice on which any leg has
$k/(aH)$ above a caller-supplied maximum, because the simple representation
is valid only there.

## Validation

The straight trajectory is the equal-mass double-quadratic model,
$m_1=m_2=10^{-6}M_{\mathrm{Pl}}$. Its potential is rotationally symmetric, so
the entropy direction never sources $\zeta$ and the result is single-field.
For these runs the pivot exits 40 e-folds before the end of inflation, where
$\epsilon\approx0.0125$ and $\eta=d\ln\epsilon/dN\approx0.025$. The output
slice is 6 e-folds after exit, and the start is 5 e-folds before exit.

| Check | Result |
| --- | --- |
| $\mathcal P_\zeta$ from the transported covariance against $\mathcal P_{\mathcal R}$ of the mode solver, turning benchmark | relative $6\times10^{-9}$ (regular test at a start of 3 e-folds: $<10^{-5}$) |
| straight-trajectory $f_{\mathrm{NL}}$ one e-fold apart on superhorizon scales | 0.03054 and 0.03065, within the 1% start-depth error |
| against Maldacena's leading-order result, shapes $(1,1,1)$, $(1,1,0.5)$, $(1,1,0.2)$ | 0.0307, 0.0260, 0.0221 against 0.0296, 0.0255, 0.0217 (1.7--3.6%, next order in $\epsilon$) |
| squeezed consistency relation $\tfrac{5}{12}(1-n_s)=0.021001$, with $n_s$ from the mode solver | $f_{\mathrm{NL}}=0.021277$ at $k_3/k_1=0.1$ and $0.021071$ at $0.05$; the deviation falls as $(k_3/k_1)^2$ |
| $P_\zeta$ and $B_\zeta$ in sheared field coordinates against Cartesian | $<10^{-8}$ and $<10^{-7}$ relative |

Maldacena's result is written with $\eta=\dot\epsilon/(H\epsilon)$ as

$$
B_\zeta=\frac{(2\pi)^4\mathcal P_\zeta^2}{\prod_ik_i^3}\Bigl[
\frac{\eta}{8}\sum_ik_i^3+\frac{\epsilon}{8}\Bigl(-\sum_ik_i^3
+\sum_{i\neq j}k_ik_j^2+\frac{8}{k_t}\sum_{i<j}k_i^2k_j^2\Bigr)\Bigr],
$$

evaluated with $\epsilon$ and $\eta$ at the crossing of $k_t/3$. The consistency
relation is exact in slow roll, so it is the sharper check. The first
comparison runs in the regular test suite. The others are ignored release-mode
tests (about three minutes together with the P3d convergence tests); run them
with `cargo test --release -- --ignored`.

### Turning two-field benchmark: the axion-quartic model

Section 9.1 of Dias, Frazer, Mulryne, and Seery
([arXiv:1609.00379v2](https://arxiv.org/abs/1609.00379), JCAP 12 (2016) 033)
evolves the axion-quartic model of Elliston et al.
([arXiv:1203.6844](https://arxiv.org/abs/1203.6844)),

$$
V=\frac14g\phi^4+\Lambda^4\Bigl(1-\cos\frac{2\pi\chi}{f}\Bigr),\qquad
g=10^{-10},\quad\Lambda^4=\Bigl(\frac{25}{2\pi}\Bigr)^2gM_{\mathrm{Pl}}^4,\quad
f=M_{\mathrm{Pl}},
$$

from $\phi=23.5M_{\mathrm{Pl}}$ and $\chi=f/2-10^{-3}M_{\mathrm{Pl}}$ with
slow-roll velocities. The axion leaves its hilltop about 60 e-folds after the
start. That turn drives $f_{\mathrm{NL}}$ through a minimum near $-39$ and a
maximum near $+17$ before it settles at about $13$. Their Fig. 7 plots the
equilateral $f_{\mathrm{NL}}$ for $k_t/3$ leaving the horizon 14 e-folds
after the start. The figure is a vector graphic, so its plotted samples
are recovered exactly by
[`derivations/digitize_dfms_figure7.py`](derivations/digitize_dfms_figure7.py).

The same model and configuration here, started 5 e-folds before exit, give

| $N$ | Fig. 7 of DFMS | this repository |
| --- | --- | --- |
| 59.84 | $-1.5825$ | $-1.5795$ |
| 61.20 | $-22.164$ | $-22.134$ |
| 62.02 | $-39.2245$ | $-39.2236$ |
| 63.10 | $-13.321$ | $-13.340$ |
| 64.60 | $6.1135$ | $6.1067$ |
| 65.55 | $16.9603$ | $16.9593$ |
| 66.23 | $12.7325$ | $12.7328$ |
| 67.05 | $13.1427$ | $13.1430$ |
| 68.00 | $13.1162$ | $13.1161$ |

The differences are of the size expected from reading a plotted curve.
Independently of the transport, a numerical $\delta N$ calculation perturbs the
flat slice 6 e-folds after exit, integrates the background to the
uniform-density surface of the output slice, and forms
$\tfrac65f_{\mathrm{NL}}=N_IN_JN^{IJ}/(N_KN^K)^2$. Its step is 1% of the
axion's distance from the hilltop, and the result is stable to $10^{-4}$ under
changes of that step. It gives $13.106$, against $13.120$ from the transport;
the 0.1% difference is of the size of the intrinsic and slow-roll-velocity
corrections that $\delta N$ omits. Both comparisons are in the ignored test
`axion_quartic_bispectrum_matches_dfms_and_numerical_delta_n`.

## Bispectrum runs and CSV output

`CovariantMultifieldSystem::run_bispectrum` (in `crates/phys-type/src/bispectrum_run.rs`)
surveys the background once and fixes $k_{\mathrm{pivot}}=aH$ with the
convention of the spectrum run. Each `TriangleShape` gives $k_t/k_{\mathrm{pivot}}$
and the Fergusson--Shellard coordinates used in Eq. (9.4) of DFMS,

$$
k_1=\frac{k_t}{4}(1+\alpha+\beta),\qquad
k_2=\frac{k_t}{4}(1-\alpha+\beta),\qquad
k_3=\frac{k_t}{2}(1-\beta).
$$

For each shape the run starts the canonical transport when the shortest leg
has $k/(aH)=e^{\Delta N}$, evolves it to the output slice, and projects it.
It rejects an output slice on which any leg is not outside the horizon.

`BispectrumRun::write_csv` writes one row per shape with the canonically
ordered $k_i/k_{\mathrm{pivot}}$, the start and output e-folds, the largest
$k/(aH)$ on the output slice, $\mathcal P_\zeta(k_i)$,
$\mathcal B_\zeta=(k_1k_2k_3)^2B_\zeta$, $f_{\mathrm{NL}}$, the intrinsic fraction
$N_aN_bN_c\alpha^{abc}/B_\zeta$, and the imaginary residual of the ordered
contraction. `write_phase_space_csv` writes all 64 components of the canonical
$(k_1k_2k_3)^2\alpha^{abc}$ before projection. `write_bispectrum_metadata`
writes `#` lines naming:

- the model, supplied by the caller;
- the method and initial-condition prescription;
- the gauge representation, with a reference to the corrected cross
  coefficient;
- the sign of $\zeta$;
- the normalization of $P$, $B$, $\mathcal P$, and $\mathcal B$;
- the $f_{\mathrm{NL}}$ convention;
- the shape coordinates and leg ordering;
- the pivot and output slices;
- the absence of reheating;
- the units of the phase-space columns.

A regular test recomputes $f_{\mathrm{NL}}$ from the dimensionless columns alone
and checks that every one of these keys is present.

## Triangle-grid convergence

A reported shape is only meaningful if the grid resolves it. On the straight
single-field benchmark at $k_t=3k_{\mathrm{pivot}}$, the isosceles line $\alpha=0$
was sampled from $\beta=0.1$ to $0.7$ at spacing $0.05$. Along it,
$f_{\mathrm{NL}}$ rises from $0.0259$ to $0.0305$ near equilateral and falls to
$0.0240$. Linear interpolation from a grid of spacing $h$, tested at the
shapes that grid omits, misses by at most

| Spacing $h$ | Largest relative interpolation error |
| --- | --- |
| 0.2 | $2.8\times10^{-2}$ |
| 0.1 | $5.8\times10^{-3}$ |

The error falls by 4.9 as $h$ halves, as expected for a smooth shape, and at
$h=0.1$ it is comparable to the start-depth error of each point. A shape
reported on a grid of spacing 0.1 or finer is therefore resolved at the
percent level. The check is the ignored test
`isosceles_shape_is_resolved_by_the_reported_grid`.

## Scope and limitations

- Tree level only; loop corrections are not included.
- The projection uses a superhorizon slice during inflation. There is no
  reheating or post-inflationary transfer, which matters whenever $\zeta$ is
  not conserved at the end of inflation.
- The start-depth error falls only as $x_{\mathrm{init}}^{-1}$ (0.23% from 7
  e-folds before exit on the double-quadratic benchmark), because the
  odd-momentum blocks of the initial condition keep only their first
  surviving order.
- Each triangle re-evaluates the cubic vertices at every stage and costs
  seconds to a minute in release mode; squeezed triangles, whose long legs
  start much deeper, are the most expensive.
- Two fields only, through the `TwoFieldBackgroundProvider` boundary.
