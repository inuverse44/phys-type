# Bispectrum amplitudes on the constant-turn orbit

Status: design. This is the next step after the squeezed-limit exponent of
[`constant-turn-exact-solutions.md`](constant-turn-exact-solutions.md#squeezed-bispectrum).

## Goal

The transport already reproduces the squeezed exponent $3/2-\nu$ of
Huenupi, Muñoz, Palma, and Sypsas,
[arXiv:2607.14529](https://arxiv.org/abs/2607.14529) (hereafter [HMPS2]). Its
amplitudes are not yet tested at strong mixing, where no expansion in the turn
rate exists. This step compares them with an independent in-in calculation:

1. derive the decoupling-limit cubic action of the circular orbit in
   `PolarFieldSpace` in the variables of [HMPS2];
2. evaluate the tree-level bispectrum of each cubic operator with a separate
   numerical in-in code built on the exact linear modes;
3. compare the sum with the transport in the equilateral and squeezed limits,
   and the operators of the [HMPS2] basis with its closed forms.

Item 3 tests two things independently. The transport is tested against the
in-in calculation of the same action. The paper is tested against the in-in
calculation of its own operators.

## The cubic action of the orbit

On $ds^2=d\rho^2+\rho^2d\theta^2$ with
$V=V_0-c\theta+\tfrac12m^2(\rho-\rho_0)^2$, the orbit is $\rho=R$,
$\dot\theta=\Omega$. In the decoupling limit (metric perturbations dropped,
error $O(\epsilon)$), write $\theta=\Omega t+\delta\theta$ and $\rho=R+\sigma$.
On flat slices $\zeta=-H\delta\theta/\Omega$. With
$\sqrt{2\epsilon}\,m_{\mathrm{Pl}}=R\Omega/H$, the canonical curvature
perturbation of [HMPS2] Eq. (9) is

$$
\phi\equiv\sqrt{2\epsilon}\,m_{\mathrm{Pl}}\zeta=-R\,\delta\theta .
$$

The matter Lagrangian
$\tfrac12(\dot\rho^2-(\nabla\rho)^2/a^2)+\tfrac12\rho^2(\dot\theta^2-(\nabla\theta)^2/a^2)-V$
with $\rho^2\dot\theta^2=(1+\sigma/R)^2(R\Omega-\dot\phi)^2$ expands as follows. The
linear terms vanish on the orbit, by the angular attractor and the radial
balance. The quadratic part is

$$
\mathcal L_2=\tfrac12\dot\phi^2-2\Omega\sigma\dot\phi+\tfrac12\Omega^2\sigma^2
+\tfrac12\dot\sigma^2-\tfrac12m^2\sigma^2-\tfrac{1}{2a^2}\big[(\nabla\phi)^2+(\nabla\sigma)^2\big],
$$

which is [HMPS2] Eq. (10) with

$$
\lambda=\frac{2\Omega}{H},\qquad \mu^2=m^2+3\Omega^2 .
$$

The sign of $\lambda$ follows from the sign of $\zeta$. The earlier document
quotes $|\lambda|=2|\eta_\perp|$; the sign is fixed here because the cubic
couplings below are odd in $\lambda$. The cubic part is

$$
\mathcal L_3=\frac{\sigma}{R}\Big(\dot\phi^2-\frac{(\nabla\phi)^2}{a^2}\Big)
-\frac{\Omega}{R}\sigma^2\dot\phi ,
$$

with no $\sigma^3$, because the radial potential is quadratic and the field
space is flat. In the basis of [HMPS2] Eq. (3), using
$\sqrt\epsilon\,m_{\mathrm{Pl}}=R\Omega/(\sqrt2H)$,

| operator | coefficient of the orbit |
| --- | --- |
| $\tilde\alpha\,\sqrt{\epsilon}m_{\mathrm{Pl}}\,\sigma\dot\zeta^2/\sqrt2$ | $\tilde\alpha=\lambda$ |
| $\tilde\beta H\sigma^2\dot\zeta/2$ | $\tilde\beta=-\lambda^2/2$ |
| $\tilde\gamma H^2\sigma^3/(3!\,m_{\mathrm{Pl}}\sqrt{2\epsilon})$ | $\tilde\gamma=0$ |
| $-\sigma(\nabla\phi)^2/(Ra^2)$ | outside the basis |

In the covariant basis, [HMPS2] Eq. (8), these are $\alpha=\lambda$,
$\beta=3\lambda^2/2$, and $\gamma=3\lambda^3/2$. The couplings scale as
$1/R\propto1/\sqrt\epsilon$, so $f_{\mathrm{NL}}$ is $O(1)$ in slow roll.
The decoupling limit and the nonlinear relation between $\zeta$ and $\delta\theta$
change it only at relative $O(\epsilon)$, which is $10^{-5}$ for the benchmark.

The gradient operator does not change the squeezed exponent, which is
universal ([`constant-turn-exact-solutions.md`](constant-turn-exact-solutions.md#squeezed-bispectrum)),
but it does change the amplitudes. The closed forms of [HMPS2] therefore
cannot be compared with the transport directly.

## Independent in-in code

The in-in code must share nothing with the transport except the model
parameters $(\lambda,\mu/H)$.

- **Modes.** Solve [HMPS2] Eqs. (24)–(25) in $z=-k\tau$ for the two
  Bunch–Davies doublets, started at $z_0=e^{6}$ or deeper, as plain complex
  ODEs. The mode solver of the core is not used. The power spectrum from
  these modes must reproduce [HMPS2] Eq. (147), which the core already
  matches (the tables of the earlier document); this cross-checks the two
  independent mode codes.
- **Vertices.** Each operator gives the tree-level contraction of the form of
  [HMPS2] Eq. (198): $-2\,\mathrm{Im}\int d\tau\,a^{n}\prod_i K(k_i;\tau_f,\tau)$,
  with the dressed bulk-to-boundary propagators [HMPS2] Eqs. (170)–(171)
  built from the numerical modes. The gradient operator replaces
  $D_\tau\phi(k_1)D_\tau\phi(k_2)$ by $-(\mathbf k_1\cdot\mathbf k_2)\phi(k_1)\phi(k_2)$.
- **Time integral.** The integrand oscillates as $e^{-iK\tau}$ inside the
  horizon. The $i0^+$ prescription is implemented as a damping factor
  $e^{\kappa K\tau}$, with the result extrapolated in $\kappa\to0$. Each
  integral must converge in $\kappa$, the start time, and the ODE tolerance
  separately, as the transport's start depth was.
- **Output.** $B_\zeta$ per operator, and the reduced $f_{\mathrm{NL}}$ in the
  normalization of the repository. [HMPS2] Eq. (197) has the same
  normalization at the equilateral point, because its $\Delta_\zeta$ is the power
  spectrum, Eq. (145).

The code belongs in `phys-type-models` as a reference calculation, like
`hkkpt`, and not in the core.

## Validation of the in-in code

1. **Massless limit.** For $\lambda=0$ the $\sigma$ modes are the massive
   Hankel modes; an operator with one external $\sigma$ then gives no
   $\zeta\zeta\zeta$, which tests the propagator bookkeeping.
2. **Weak mixing.** For $\lambda\ll1$ the $\tilde\alpha$ bispectrum is the
   one-transfer quasi-single-field result of
   [Chen and Wang, arXiv:0909.0496](https://arxiv.org/abs/0909.0496), first
   order in $\lambda$, which [HMPS2] states it recovers.
3. **The paper's closed forms.** The squeezed coefficients $C_\pm(\lambda,\nu)$
   of [HMPS2] Sec. X and the equilateral $f^{(\alpha)}_{\mathrm{NL}}/\alpha$ of its
   Fig. 1 are evaluated and compared with the $\alpha$, $\beta$, and $\gamma$
   operators. Any difference is recorded in
   [`literature-discrepancies.md`](literature-discrepancies.md).

## Comparison with the transport

For the orbit cases of the earlier document,
$(\lambda,\mu^2/H^2)\in\{(1,1.25),(2,1.25),(1,2),(2,2),(2,0.81)\}$ and
$\epsilon=10^{-5}$:

- equilateral $f_{\mathrm{NL}}$ from the transport against the in-in sum of
  the four operators;
- squeezed $f_{\mathrm{NL}}(r)$ for $1\ge r\ge0.01$ against the same sum,
  including the amplitude of $r^{3/2-\nu}$;
- the share of the gradient operator, reported separately.

The expected agreement is the larger of the transport's start-depth error
(about 0.2% from seven e-folds,
[`three-point-transport.md`](three-point-transport.md)) and the in-in
convergence. A difference at $O(\epsilon)$ or $O(1/z_0)$ is accepted; a larger
one is a defect of one of the two codes and is resolved before the step
closes.

## Migration

| Step | Content |
| --- | --- |
| A0 | this design |
| A1 | numerical in-in reference in `phys-type-models` with the mode and power-spectrum checks and validations 1–2 |
| A2 | the comparison with the transport, as an ignored test and a figure |
| A3 | the comparison with the closed forms of [HMPS2] (validation 3) |

## Open questions

- How to integrate the vertices at $r=0.01$, where the long mode exits about
  five e-folds after the short ones and the integrand spans both scales; the
  modes and the integral may need to be advanced together with an adaptive
  step.
- How far into the heavy regime, $\mu>3H/2$, the comparison should go. The
  squeezed oscillation needs $r\lesssim10^{-3}$, which is expensive for
  both codes.
