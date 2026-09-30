# Nonminimally assisted inflation

## Purpose

This is the first external model validation after the non-Gaussianity
milestone. The target is Hyun, Kim, Kodama, Park, and Takahashi,
[*Nonminimally Assisted Inflation: A General Analysis*](https://arxiv.org/abs/2302.05866)
(arXiv:2302.05866v2, hereafter [HKKPT]), and within it the loop-inflation
example of Class I (Sec. 4.1). The contract from the roadmap is that the
Einstein-frame metric and potential enter as ordinary `FieldSpaceGeometry` and
potential providers, with no special case in the solvers. The dimension,
index-variance, Riemann-identity, and covariance checks must all stay active.

[HKKPT] predicts $n_s$, $r$, and $f^{\mathrm{local}}_{\mathrm{NL}}$ from slow-roll
$\delta N$ expressions, Eqs. (3.1), (3.2), and (3.5). This repository evolves the
same model exactly: the full background, the linear modes, and the tree-level
bispectrum, with no slow-roll approximation after the initial conditions. It
then compares with those expressions.

## Model

With $f(s)=\tfrac{M_{\mathrm{Pl}}^2}{2}[1+\xi_2(s/M_{\mathrm{Pl}})^2]$ (Eq. (2.5)),
the Einstein frame has, in the coordinates $(\phi,s)$ (Eqs. (2.6)--(2.9)),

$$
G=\operatorname{diag}(\mathcal K_1,\mathcal K_2),\qquad
\mathcal K_1=\frac{1}{1+\xi_2s^2},\qquad
\mathcal K_2=\mathcal K_1+6\kappa\xi_2^2s^2\mathcal K_1^2,\qquad
V=\frac{V_J(\phi)}{(1+\xi_2s^2)^2},
$$

in reduced Planck units, with $\kappa=1$ (metric) or $0$ (Palatini).
`NonminimalAssistedFieldSpace` implements this geometry, and
`NonminimalLoopInflationPotential` implements it with
$V_J=\Lambda^4(1+\alpha\ln\phi)$ (Eq. (4.1)). Both live in
`src/nonminimal.rs`.

Neither field is canonically normalized, so the geometry supplies every
connection and curvature term. Since the metric depends only on $s$ in two
dimensions, the curvature is set by the Gaussian curvature

$$
K=\frac{\mathcal K_1\mathcal K_1'\mathcal K_2'-2\mathcal K_1\mathcal K_1''\mathcal K_2
+\mathcal K_1'^2\mathcal K_2}{4\mathcal K_1^2\mathcal K_2^2},\qquad
R^I{}_{JKL}=K(\delta^I_KG_{JL}-\delta^I_LG_{JK}),
$$

and $\nabla_MR^I{}_{JKL}=(\partial_MK)(\delta^I_KG_{JL}-\delta^I_LG_{JK})$. The
derivatives of $\mathcal K_{1,2}$ to third order come from `ThirdOrderJet`.
The $K$ formula reproduces $K=-1$ for the hyperbolic metric. Unit tests check
each piece against an independent finite-difference construction:

- the Christoffel symbols against the Levi-Civita connection of the
  finite-differenced metric;
- $\partial\Gamma$ against differenced Christoffel symbols;
- the Riemann tensor against $\partial\Gamma+\Gamma\Gamma$ in the crate
  convention, where the same formula is first shown to reproduce the built-in
  hyperbolic geometry;
- the Riemann symmetries through `validate_riemann_symmetries`;
- $\nabla R$ against the covariant finite difference of $R$;
- the potential and its automatic-differentiation gradient against the
  closed form.

## Configuration

All results use the reference example $\alpha=0.1$, $s_*=10^{-3}$, $\kappa=1$,
and $\Lambda^4=10^{-10}$; $\Lambda$ does not affect $n_s$ or $f_{\mathrm{NL}}$. The
pivot leaves the horizon 60 e-folds before $\epsilon_H=1$. The initial $s$ is
tuned so that $s=s_*$ at the pivot. Linear modes and the three-point transport
start 5 e-folds before exit.

**Evaluate at the end of inflation.** The assistant field keeps rolling, so
$\zeta$ is not conserved on superhorizon scales before inflation ends. For
$\xi_2=0.022$:

| E-folds before $\epsilon_H=1$ | 10 | 5 | 2 | 1 | 0.5 | 0.1 | 0.01 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $n_s$ from the mode solver | 0.847 | 0.843 | 0.859 | 0.888 | 0.924 | 0.964 | 0.973 |

[HKKPT] evaluates the observables at the end of inflation. All comparisons
below therefore use a slice $10^{-3}$ e-folds before $\epsilon_H=1$. No
post-inflationary evolution is included.

## Spectral index

Three independent determinations agree at $\xi_2=0.022$:

| Method | $n_s$ |
| --- | --- |
| mode solver, $d\ln\mathcal P_{\mathcal R}/d\ln k$ | 0.97408 |
| numerical $\delta N$ from the horizon-crossing slice, massless crossing statistics | 0.97470 |
| [HKKPT] Eq. (3.1) with the numerical $\phi_*$, $s_*$, and $s_e$ | 0.97391 |

Eq. (3.1) follows the exact tilt across the coupling range, even though $n_s$
itself moves by 0.08:

| $\xi_2$ | 0.021 | 0.022 | 0.023 | 0.0235 | 0.024 | 0.025 |
| --- | --- | --- | --- | --- | --- | --- |
| numerical $n_s$ | 0.9789 | 0.9741 | 0.9626 | 0.9525 | 0.9387 | 0.9001 |
| Eq. (3.1) with numerical $s_e$ | 0.9786 | 0.9739 | 0.9627 | 0.9528 | 0.9391 | 0.9003 |
| Eq. (3.1) with $s_e$ from Eq. (3.9), $N_e=60$ | 0.977 | 0.969 | 0.950 | | | |
| numerical $s_e$ | 0.135 | 0.169 | 0.212 | 0.238 | 0.266 | 0.334 |

The last two rows explain the only difference from the published curves.
[HKKPT] obtains $s_e$ from the slow-roll e-fold relation, Eq. (3.9) or its
approximation (3.12). Inserting the numerical $s_*$ and $s_e$ into Eq. (3.9)
gives 58.3 e-folds instead of 60. The difference is not the choice of end
condition: on the numerical trajectory, $\epsilon^{(0)}=1$ ($\phi=0.093$) is
passed only about 0.1 e-fold before $\epsilon_H=1$. The 1.7 e-folds instead
accumulate over the whole trajectory. They are 3% of the e-fold count, of the
size expected from next-order slow-roll corrections to the $s(N)$ relation
when $\eta^{\sigma\sigma}\approx-4\xi_2$. Because $s_e\propto e^{4\xi_2N}$, the
difference shifts $s_e$ by 16% and moves $n_s$ noticeably. The numerical
$n_s=0.965$ occurs at $\xi_2\approx0.0228$, whereas the published prescription
gives $\xi_2\approx0.0222$. This is a shift of about 3% within
the range the reference quotes; the conclusion that a small nonminimal
coupling rescues loop inflation is unchanged.

## Non-Gaussianity

At $\xi_2=0.022$, with $f_{\mathrm{NL}}$ in the Maldacena/Planck sign convention
of this repository:

| Method | $f_{\mathrm{NL}}$ |
| --- | --- |
| tree-level transport, equilateral, start depths 4--7 | $-0.25141$ (stable to $10^{-5}$) |
| tree-level transport, $k_3/k_1=0.1$ | $-0.35686$ |
| numerical $\delta N$ at horizon crossing, covariant $N_{;IJ}$ | $-0.1989$ |
| [HKKPT] Eq. (3.5), numerical inputs, $\phi_e$ from $\epsilon^{(0)}=1$, sign converted | $-0.2162$ |

The sign conversion is needed because the printed Eq. (3.5) gives
$-\tfrac65f_{\mathrm{NL}}$ in the Lyth--Rodriguez convention. Its single-field
limit is $-\tfrac56(2\epsilon-\eta)$, the opposite of the Planck-sign $\delta N$
result; see [`literature-discrepancies.md`](literature-discrepancies.md),
entry 7. The end-of-inflation value of $\phi$ matters in Eq. (3.5) through
$\eta^{\phi\phi}_e$. Evaluating it at the numerical $\epsilon_H=1$ point instead,
$\phi_e=0.007$, gives $\eta^{\phi\phi}_e\approx-4\times10^3$ and
$f_{\mathrm{NL}}=-0.55$.

The two $\delta N$ values agree with each other to 9%. The transport differs from
them by 25% (equilateral) to 80% ($k_3/k_1=0.1$). Although the transport result
is converged in start depth, the difference is not yet understood, for three
reasons:

- $\delta N$ assumes massless, uncorrelated field perturbations at horizon
  crossing. The assistant direction is tachyonic, with
  $\eta^{\sigma\sigma}\approx-4\xi_2$. In a calibration run with the target
  half an e-fold before the end, applying the same formula 3 e-folds after
  crossing instead of at crossing changed $f_{\mathrm{NL}}$ by a factor of two.
  Its crossing-time error is therefore not small.
- At the end of inflation the intrinsic and gauge parts of $B_\zeta$ cancel
  to about 1%: the intrinsic fraction is 105. Such a cancellation is sensitive
  to the representation of the gauge transformation, although the result
  shows no sensitivity to start depth.
- The squeezed limit of a multifield model with a light, evolving isocurvature
  direction need not approach the shape-independent $\delta N$ value.

Resolving this requires an independent calculation of the full shape, for
example $\delta N$ with the numerically evolved field covariance and field
three-point function on a later slice, or the local gauge representation of
Eqs. (2.33a)--(2.33b) of arXiv:1410.3491 with its curvature terms. It is an
open item. Regardless of which value holds, $|f_{\mathrm{NL}}|\lesssim0.4$, well
inside the Planck bound quoted by the reference.

## Tests

The regular suite runs the geometry and potential tests in `src/nonminimal.rs`.
The ignored release-mode tests are:

- `assisted_loop_spectral_index_agrees_with_delta_n_and_eq_3_1`;
- `assisted_loop_eq_3_1_tracks_the_numerical_tilt_across_xi`;
- `assisted_loop_f_nl_is_converged_and_bracketed_by_delta_n`, which fixes the
  transport convergence and the $\delta N$/Eq. (3.5) agreement, and bounds the
  open transport difference at 40%.

## Next steps

1. Resolve the transport-versus-$\delta N$ difference in $f_{\mathrm{NL}}$ as
   described above.
2. Reproduce the $s_*=0.1$ and $\alpha=1$ branches of the reference.
3. Add a general end-of-inflation event for the Class II power-law and hybrid
   examples.
