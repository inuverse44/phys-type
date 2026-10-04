# Frames: the Jordan frame, the Einstein frame, and what is reported

This document records how a theory stated in the Jordan frame reaches the
engine, which evolves an Einstein-frame sigma model, and how quantities are
translated back. It accompanies `crates/phys-type/src/theory.rs` and the
architecture in [`architecture.md`](architecture.md). It is meant to be
checked against whenever a paper states which frame its numbers refer to.

## The theory class

A `GeneralizedNonminimal<N>` theory has the Jordan-frame action

$$
S=\int d^4x\sqrt{-g}\left[\frac{M_{\mathrm{pl}}^2}{2}\Omega^2(\phi)R
-\frac12k_{IJ}(\phi)g^{\mu\nu}\partial_\mu\phi^I\partial_\nu\phi^J-V_J(\phi)\right],
$$

with $\Omega^2=2f(\phi)/M_{\mathrm{pl}}^2$ for the coupling $f(\phi)R$, in the metric or
the Palatini formulation. Single-field Higgs inflation
($\Omega^2=1+\xi\phi^2/M_{\mathrm{pl}}^2$, $N=1$), the nonminimally assisted model of
arXiv:2302.05866 ($\Omega^2=1+\xi_2s^2/M_{\mathrm{pl}}^2$, $N=2$), and multi-field
nonminimal models with field-dependent kinetic terms are instances. Derivative
couplings (Horndeski, DHOST) are not; see the admissibility discussion in
[`architecture.md`](architecture.md).

## The map to the Einstein frame

With $g^E_{\mu\nu}=\Omega^2g^J_{\mu\nu}$ the action becomes canonical gravity plus a
sigma model with

$$
G_{IJ}=\frac{k_{IJ}}{\Omega^2}+\kappa\,\frac32M_{\mathrm{pl}}^2\,
\frac{\partial_I\Omega^2\,\partial_J\Omega^2}{\Omega^4},\qquad
V=\frac{V_J}{\Omega^4},
$$

$\kappa=1$ in the metric and $\kappa=0$ in the Palatini formulation, where the
connection is independent of the metric and the conformal transformation
generates no derivative of $\Omega$. The multi-field metric form is that of
[Kaiser, arXiv:1003.1159](https://arxiv.org/abs/1003.1159), who also shows that
more than two nonminimally coupled fields cannot in general be brought to
canonical kinetic form; the engine never needs canonical form, since it is
covariant in field space. The metric–Palatini difference for inflation is
discussed by [Bauer and Demir, arXiv:0803.2664](https://arxiv.org/abs/0803.2664).

`einstein_frame(theory)` returns the potential and a `JetGeometry` whose
metric is evaluated with automatic differentiation: $\Omega^2$ four times, so
that $\nabla R$ (three derivatives of $G$) is exact, and $k_{IJ}$ three times.
The map reproduces the hand-written assisted geometry (a core test fixture)
and its potentials to $10^{-12}$, metric and Palatini, through the potential's third
derivatives and the covariant Riemann derivative.

**Validity along the trajectory.** The map requires $\Omega^2>0$, and the engine
requires a positive-definite $G_{IJ}$ (no ghost). The background survey checks
both on every accepted step through `FieldSpaceGeometry::check_domain` and
stops with `SpectrumRunError::OutsideGeometryDomain` at the first violation.

## Correspondence of background quantities

From $g^E=\Omega^2g^J$ on a flat FRW background, $dt_E=\Omega\,dt_J$ and
$a_E=\Omega\,a_J$. `EinsteinFrame::correspondence` returns, at an Einstein-frame
state,

| Jordan-frame quantity | expression |
| --- | --- |
| scale factor | $a_J=a_E/\Omega$ |
| e-fold number | $N_J=\ln a_J=N_E-\ln\Omega$ |
| Hubble rate | $H_J=\Omega\left(H_E-\dot\phi^I\partial_I\ln\Omega\right)$, dot $=d/dt_E$ |
| field velocity | $d\phi^I/dt_J=\Omega\,\dot\phi^I$ |
| potential | $V_J=\Omega^4V$ |

The field coordinates are the same in both frames; only the metric and the
time coordinate change.

**E-fold counts depend on the frame.** An interval differs by
$\Delta N_J=\Delta N_E-\Delta\ln\Omega$. "The pivot leaves the horizon 60
e-folds before the end" must therefore name its frame. For the assisted
examples of arXiv:2302.05866, $\ln\Omega$ changes between the pivot and the end
by $\tfrac12\ln[(1+\xi_2s_e^2)/(1+\xi_2s_*^2)]$, below $10^{-3}$ e-folds for
$s_*=10^{-3}$ but about 0.04 e-folds for loop inflation with $\alpha=0.1$ and
$s_*=0.1$ ($s_e\simeq2.4$, $\xi_2\simeq0.014$). The validations in this
repository count Einstein-frame e-folds, as the reference's $\delta N$ formulas
do.

## Which curvature perturbation is reported

The engine reports the Einstein-frame comoving curvature perturbation $\zeta$
on the output slice, and $P_\zeta$, $n_s$, and $f_{\mathrm{NL}}$ built from it.

- For a single field, $\zeta$ is the same in both frames.
- With several fields it is not: [White, Minamitsuji, and Sasaki,
  arXiv:1205.0656](https://arxiv.org/abs/1205.0656) show that the Jordan- and
  Einstein-frame curvature perturbations of multi-field nonminimal models
  differ, the difference being sourced by the isocurvature perturbations, so
  that adiabaticity itself is frame dependent; see also
  [arXiv:2305.15326](https://arxiv.org/abs/2305.15326) on the comoving
  curvature perturbation in the two frames.
- The two agree once the isocurvature perturbation has decayed or the
  evolution is adiabatic. In the assisted examples the observables are taken
  at the end of inflation, where the isocurvature perturbation need not have
  decayed, so the reported numbers are Einstein-frame quantities. Whether a
  given paper's predictions are Jordan- or Einstein-frame $\zeta$ is to be
  checked when it is compared.

A Jordan-frame $\zeta$ (the Jordan-frame uniform-density or comoving slicing
of the same solution) is a future output of the observable layer; it needs
no new dynamics, only a different projection, and is listed as an open
question in [`architecture.md`](architecture.md).
