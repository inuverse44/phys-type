# Bispectrum amplitudes on the constant-turn orbit

Status: A1 (the in-in reference) and A2 (the comparison with the transport)
done; A3 designed. This is
the next step after the squeezed-limit exponent of
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

`phys_type_models::in_in` shares nothing with the transport except the model
parameters $(\lambda,\mu/H)$. Units are $H=1$, and a cubic coefficient is
given times $c=\sqrt{2\epsilon}\,m_{\mathrm{Pl}}$, so that the reduced
$f_{\mathrm{NL}}$ does not depend on $c$. For the orbit the coefficients of
$\sigma\dot\phi^2$, $\sigma(\nabla\phi)^2/a^2$, and $\sigma^2\dot\phi$ are $\lambda/2$,
$-\lambda/2$, and $-\lambda^2/4$; for [HMPS2] Eq. (3) they are $\tilde\alpha/2$,
$\tilde\beta/2$, and $\tilde\gamma/6$.

- **Modes.** [HMPS2] Eqs. (24)–(25) are integrated in $z$ for unit wavenumber
  by fourth-order Runge–Kutta on a uniform grid in $s=z-1$ above $z=1$ and
  $s=\ln z$ below it (step 0.01), and interpolated by cubic Hermite
  polynomials. Scale invariance gives every other wavenumber, so one table
  serves all three legs, and the long leg of a squeezed triangle starts as
  deep as the short ones.
- **Vacuum.** The doublets start from the two positive-frequency solutions
  of the coupled equations,
  $\phi,\sigma=z^{1+a}e^{iz}\sum_np_n,q_nz^{-n}$ with $a=\pm i\lambda/2$, whose
  coefficients follow from a two-by-two recursion (eight terms). The
  uncoupled modes $u_0$ and $u_\mu$ are not enough: with
  $\phi'=0$ in the entropy doublet, the curvature momentum
  $\Pi_\phi\propto D_z\phi=\lambda\sigma/z$ does not vanish, which is an $O(1)$
  error in $\dot\phi=\lambda\sigma+\dots$. It moved the power spectrum by only
  $10^{-3}$ but the orbit's $f_{\mathrm{NL}}$ by a factor of 3000. With the
  series, the power spectrum no longer depends on the start, $z_0=e^6$ or
  $e^{7.5}$.
- **Vertices.** With $\langle Q(\tau_f)\rangle=2\int dt\,\mathrm{Im}\langle Q(\tau_f)H_I(t)\rangle$
  and $H_I=-\int d^3x\,a^3\mathcal L_3$,
  $B=-2\sum g\,\mathrm{Im}\int d\tau\,a^4\sum_{\mathrm{perm}}\prod_iG_i$, where
  $G_i=\sum_b\phi_b(k_i,\tau_f)X_b^*(k_i,\tau)$ is the Wightman function of
  [HMPS2] Eq. (152) for the vertex field $X$ on leg $i$. The gradient operator
  contributes $-(\mathbf k_a\cdot\mathbf k_b)\tau^2$ for the legs of its two
  $\phi$ fields.
- **Time integral.** Eight-point Gauss–Legendre panels no wider than $1/K$
  and $0.05|\tau|$, from where the shortest leg has $z_0$ to the final slice.
  The $i0^+$ prescription is a smooth window that switches the integrand off
  between $K|\tau|=z_0K/(2k_{\max})$ and $z_0K/k_{\max}$, with every derivative
  continuous, so its error falls faster than any power of the window length.
  It leaves the horizon-crossing region untouched. A damping
  factor $e^{\kappa K\tau}$, as first designed, biases it at $O(\kappa)$, and
  Richardson extrapolation leaves $O(\kappa^3)$, about 1% for practical
  $\kappa$.
- **Superhorizon commutators.** Outside the horizon $\mathrm{Im}\,G$ is a
  commutator, smaller than $\mathrm{Re}\,G$ by about $z^3$; formed from products of
  modes it cancels to roundoff by $z\sim e^{-12}$, and for a light $\sigma$ the
  result then grew without bound with the final slice. Below $z=0.1$ it is
  therefore computed as $r(z)S(z)$: $S$ is the equal-time commutator matrix of
  $(\phi,\phi',\sigma,\sigma')$, exact from $[\phi,\Pi_\phi]=[\sigma,\Pi_\sigma]=i$, and
  $r(z)$ is the response of $\phi(z_f)$ to the state at $z$, the adjoint
  equation $dr/dz=-rA(z)$ integrated from $r(z_f)=e_\phi$. Inside the horizon the
  modes are kept: the squeezed limit relies on an exact cancellation between
  times of the short legs' $\mathrm{Im}\,G$, which holds only when its real and
  imaginary parts share their integration errors. Using the response there
  put the squeezed $\dot\phi^3$ 7% off.
- **Output.** The reduced $f_{\mathrm{NL}}=\tfrac56B/(P_1P_2+P_1P_3+P_2P_3)$,
  the convention of `CurvatureCorrelators::reduced_f_nl`. [HMPS2] Eq. (197)
  has the same normalization at the equilateral point, because its
  $\Delta_\zeta$ is the power spectrum, Eq. (145).

## Validation of the in-in code

`crates/phys-type-models/tests/in_in.rs`, with unit tests of the vacuum
series, the commutator matrix, and the response in the module:

| Check | Result |
| --- | --- |
| power spectrum against [HMPS2] Eq. (147), $(\lambda,\mu^2)$ from the earlier tables | within $2\times10^{-6}$ (heavy: $10^{-9}$); a light $\sigma$ needs the final slice at $z=e^{-30}$ |
| $g\dot\phi^3$ at $\lambda=0$ against $B=3g/(k_1k_2k_3K^3)$, equilateral, folded, $r=0.1$, $r=0.01$ | within $7\times10^{-7}$ |
| operators with an external $\sigma$ at $\lambda=0$ | exactly zero |
| [HMPS2] Eq. (3) against its covariant form, Eq. (8), with the reshuffled couplings | $10^{-9}$ |
| the orbit against $(\tilde\alpha,\tilde\beta,\tilde\gamma)=(\lambda,-\lambda^2/2,0)$ plus the gradient operator | $10^{-9}$; the gradient operator is more than 10% of the total |
| weak mixing: $f^{(\tilde\alpha)}_{\mathrm{NL}}/\lambda$ at $\lambda=10^{-3},2\times10^{-3},4\times10^{-3}$ | constant to $5\times10^{-5}$, with an $O(\lambda^2)$ correction |
| orbit $f_{\mathrm{NL}}$ for $(1,1.25)$, $(2,2)$, $(2,0.81)$, equilateral and $r=0.01$, under $z_0\to e^{7.5}$, half the step, half the panels, an earlier window, and the response from $z=0.03$ | within $10^{-6}$ (observed $5\times10^{-8}$) |

The weak-mixing check tests the order in $\lambda$ but not the value. Its value
is part of A3, because [HMPS2] states that its closed forms reduce to the
quasi-single-field result of
[Chen and Wang, arXiv:0909.0496](https://arxiv.org/abs/0909.0496).

The final slice is physical, not numerical: for a light $\sigma$ the
bispectrum keeps evolving outside the horizon until $\sigma$ has decayed. For
$\mu^2=0.81$ ($\nu=1.2$) the equilateral $f_{\mathrm{NL}}$ changes by 1.5% between
$z=e^{-16}$ and $e^{-20}$. The comparison with the transport therefore uses the
same slice in both codes.

## Comparison with the transport

The transport runs the full theory on the orbit of
`phys_type_models::constant_turn` with the pivot read `after_exit` e-folds after
its exit (20, 20, 12, 12, and 30 e-folds for the cases below). The in-in
reference uses the $\lambda$ and $\mu^2$ measured at the pivot exit and the same
slice, $z=e^{-\text{after\_exit}}$ for the short legs. Triangles are isosceles,
$k_1=k_2=k_S$, $k_3=rk_S$; $\epsilon=10^{-5}$, transport from five e-folds before
exit:

| $(\lambda,\mu^2/H^2)$ | $r$ | transport | in-in | relative difference |
| --- | --- | --- | --- | --- |
| $(1,1.25)$ | 1 | $-0.229501$ | $-0.229529$ | $-1.3\times10^{-4}$ |
| | 0.1 | $-0.066440$ | $-0.066481$ | $-6.2\times10^{-4}$ |
| | 0.01 | $-0.021281$ | $-0.021321$ | $-1.9\times10^{-3}$ |
| $(2,1.25)$ | 1 | $0.104489$ | $0.104398$ | $8.7\times10^{-4}$ |
| | 0.1 | $0.103045$ | $0.103021$ | $2.4\times10^{-4}$ |
| | 0.01 | $0.031413$ | $0.031425$ | $-3.9\times10^{-4}$ |
| $(1,2)$ | 1 | $-0.343492$ | $-0.343527$ | $-1.0\times10^{-4}$ |
| | 0.1 | $-0.063482$ | $-0.063513$ | $-4.8\times10^{-4}$ |
| | 0.01 | $-0.0077760$ | $-0.0078055$ | $-3.8\times10^{-3}$ |
| $(2,2)$ | 1 | $-0.229110$ | $-0.229131$ | $-9.2\times10^{-5}$ |
| | 0.1 | $-0.043307$ | $-0.043323$ | $-3.8\times10^{-4}$ |
| | 0.01 | $-0.0094773$ | $-0.0094940$ | $-1.8\times10^{-3}$ |
| $(2,0.81)$ | 1 | $0.332179$ | $0.331877$ | $9.1\times10^{-4}$ |
| | 0.01 | $0.135855$ | $0.135849$ | $4.1\times10^{-5}$ |

The two codes agree in sign and amplitude in every case, at strong mixing
and in the squeezed limit. The differences are roughly constant in absolute
terms, $1$–$9\times10^{-5}$, except for $3\times10^{-4}$ at $(2,0.81)$, $r=1$. That
is why the relative difference grows where $|f_{\mathrm{NL}}|$ is small. Starting
the transport seven e-folds before exit instead of five changes it by at most
$1.2\times10^{-5}$ (checked for the first three cases). With $\epsilon=10^{-6}$ the
differences fall to $0.7$–$1.9\times10^{-5}$, by a factor of 1.5 to 4 rather than
ten. Part of the difference therefore scales with $\epsilon$, as the decoupling
limit implies. The rest is at the level of the transport's start-depth error,
not a physical difference.

**The gradient coupling matters.** The in-in reference splits the orbit's
$f_{\mathrm{NL}}$ into the operators of the [HMPS2] basis and the gradient operator
$-\sigma(\nabla\phi)^2/(Ra^2)$ outside it:

| $(\lambda,\mu^2/H^2)$ | $r$ | total | gradient | [HMPS2] basis |
| --- | --- | --- | --- | --- |
| $(1,1.25)$ | 1 | $-0.229$ | $-0.162$ | $-0.067$ |
| $(2,1.25)$ | 1 | $0.104$ | $-0.136$ | $0.240$ |
| $(1,2)$ | 1 | $-0.343$ | $-0.264$ | $-0.080$ |
| $(2,2)$ | 1 | $-0.229$ | $-0.259$ | $0.030$ |
| $(2,0.81)$ | 1 | $0.331$ | $-0.078$ | $0.409$ |
| $(2,2)$ | 0.01 | $-0.0095$ | $-0.0052$ | $-0.0043$ |

In the equilateral configuration the gradient operator is often the larger
part, and for $(2,2)$ the [HMPS2] basis alone would give the wrong sign. In
the squeezed limit it remains comparable. It does not change the squeezed
exponent, which is why the earlier slope test could not see it.

`constant_turn_transport_matches_the_in_in_amplitudes` (ignored, about three
minutes in release) asserts the agreement to $5\times10^{-5}$ for $(1,1.25)$ and
$(2,2)$, equilateral and $r=0.01$, with $\epsilon=10^{-6}$. It also asserts that
dropping the gradient operator misses the transport by more than 0.1. The
figure
[`figures/constant_turn_squeezed.png`](figures/constant_turn_squeezed.png)
shows the in-in values as crosses over the transport points, with the
relative difference in its bottom panel.

## Migration

| Step | Content |
| --- | --- |
| A0 | this design |
| A1 | numerical in-in reference in `phys-type-models` and its validation (done) |
| A2 | the comparison with the transport, as an ignored test and a figure (done) |
| A3 | the comparison with the closed forms of [HMPS2]: the squeezed $C_\pm(\lambda,\nu)$ of its Sec. X, the equilateral $f^{(\alpha)}_{\mathrm{NL}}/\alpha$ of its Fig. 1, and the quasi-single-field limit; differences go to [`literature-discrepancies.md`](literature-discrepancies.md) |

## Open questions

- How far into the heavy regime, $\mu>3H/2$, the comparison should go. The
  squeezed oscillation needs $r\lesssim10^{-3}$, which is expensive for
  both codes.
