# Exact constant-turn solutions as a mode-solver benchmark

## Purpose

Huenupi, Muñoz, Palma, and Sypsas,
[*Pushing the Primordial Frontier: Exact Linear Solutions in Multifield Inflation*](https://arxiv.org/abs/2606.18248)
(arXiv:2606.18248v2, hereafter [HMPS]), solve the linear two-field system of
constant mixing $\lambda=-2\Omega/H$ and constant entropy mass $\mu$ exactly in
quasi-de Sitter space. Their closed-form $\Delta_\zeta/\Delta_0$ covers the weak and
strong mixing regimes and the massless, light, and heavy entropy fields. The
companion article [arXiv:2607.14529](https://arxiv.org/abs/2607.14529) derives
the tree-level bispectra from the same solutions.

These closed forms test the covariant mode solver in the strongly mixed
regime, where no perturbative expansion in the turn rate is available.
Nothing in the solver is specialized to this background.

## Background

[`PolarFieldSpace`](../src/field_space.rs) is flat field space in polar
coordinates, $ds^2=d\rho^2+\rho^2d\theta^2$. On it the test potential

$$
V=V_0-c\,\theta+\tfrac12 m^2(\rho-\rho_0)^2
$$

has a circular orbit $\rho=R$ with $\dot\theta=\Omega$ when
$c=3HR^2\Omega$ (the angular slow-roll attractor) and $m^2(R-\rho_0)=R\Omega^2$
(radial balance). The orbit's turn rate is $\Omega$, its entropy direction is
$\rho$, and its first slow-roll parameter is $\epsilon=R^2\Omega^2/2H^2=R^2\lambda^2/8$.
The action (1) of [HMPS] writes the $\sigma$ mass as $\mu^2=m_{\mathrm{eff}}^2+\lambda^2H^2$
([arXiv:2607.14529](https://arxiv.org/abs/2607.14529), Eq. (6)), with
$m_{\mathrm{eff}}^2=V_{NN}-\Omega^2+\epsilon H^2\mathbb{R}$. With $V_{NN}=m^2$ and
$\mathbb{R}=0$,

$$
\mu^2=m^2+3\Omega^2 .
$$

The inflaton $\theta$ never stops, so the survey ends on the field-value event
`EndOfInflation::FieldValue` for $\theta$. The mode starts $6$ e-folds before
exit ($z_0=e^6$), and $\Delta_\zeta$ is read $12$ e-folds after exit ($20$ for
$\mu^2=1.25H^2$, whose $\sigma$ decays only as $z^{1/2}$). $\Delta_0=H^2/8\pi^2\epsilon$
uses $H$ and $\epsilon$ at exit. So do the $\lambda=2|\eta_\perp|$ and
$\mu^2=m^2/H^2+3\eta_\perp^2$ inserted in the closed forms; they differ from the
design values by up to $3\times10^{-4}$ relative.

## Massive entropy field

With $\epsilon=10^{-5}$, the mode solver against Eqs. (43)–(45):

| $\lambda$ | $\mu^2/H^2$ | regime | mode solver | closed form | relative difference |
| --- | --- | --- | --- | --- | --- |
| 1 | 1.25 | light, Eq. (43) | 8.85853 | 8.85945 | $-1.0\times10^{-4}$ |
| 2 | 1.25 | light, Eq. (43) | 132.967 | 132.920 | $+3.5\times10^{-4}$ |
| 0.5 | 2 | conformal, Eq. (44) | 1.46503 | 1.46501 | $+1.1\times10^{-5}$ |
| 1 | 2 | conformal, Eq. (44) | 3.67564 | 3.67569 | $-1.5\times10^{-5}$ |
| 2 | 2 | conformal, Eq. (44) | 42.6270 | 42.6247 | $+5.6\times10^{-5}$ |
| 1 | 4 | heavy, Eq. (45) | 1.53052 | 1.53049 | $+2.2\times10^{-5}$ |
| 2 | 4 | heavy, Eq. (45) | 7.51527 | 7.51509 | $+2.4\times10^{-5}$ |
| 1 | 9 | heavy, Eq. (45) | 1.09992 | 1.09989 | $+2.9\times10^{-5}$ |

The remaining difference is the slow-roll correction that [HMPS] neglects. It
is proportional to $\epsilon$ and grows with $\lambda$. For example, at $\lambda=2$
and $\mu^2=1.25H^2$ it is $1.0\%$ for $\epsilon=2\times10^{-4}$ and $0.22\%$ for
$\epsilon=5\times10^{-5}$.

## Massless entropy field

For $\mu=0$, $\sigma$ does not decay and $\zeta$ grows on superhorizon scales. As
printed, Eq. (39) of [HMPS] contains $\ln(2/z_{\mathrm{end}})$. The mode solver, and
an independent integration of Eqs. (8)–(9) of [HMPS]
([`derivations/constant_turn_exact_solutions.py`](derivations/constant_turn_exact_solutions.py)),
agree instead with $\ln(1/2z_{\mathrm{end}})$:

$$
\frac{\Delta_\zeta}{\Delta_0}=\frac{\lambda\sinh\pi\lambda}{\pi(1+\lambda^2)}
\left[2+\ln\frac{1}{2z_{\mathrm{end}}}-2\gamma_E-\operatorname{Re}\psi\!\left(\frac{i\lambda}{2}\right)\right]^2
+\frac{\pi\lambda}{2(1+\lambda^2)}\coth\frac{\pi\lambda}{2}.
$$

This is what Eq. (38) of [HMPS] implies, since its $\mathcal H_w$ contains
$-2\ln(2z)$. At $z_{\mathrm{end}}=4.54\times10^{-5}$:

| $\lambda$ | $\epsilon$ | mode solver | corrected (39) | relative difference | printed (39) |
| --- | --- | --- | --- | --- | --- |
| 0.5 | $3\times10^{-6}$ | 34.2449 | 34.2439 | $2.9\times10^{-5}$ | $+27\%$ |
| 1 | $3\times10^{-6}$ | 202.833 | 202.778 | $2.7\times10^{-4}$ | $+28\%$ |
| 2 | $3\times10^{-7}$ | 3450.47 | 3449.17 | $3.8\times10^{-4}$ | $+29\%$ |

The massless comparison needs the exit values of $\lambda$ and $\mu$ to match
their targets, so the test corrects the orbit parameters three times. Because
the slow-roll correction accumulates over the superhorizon growth, it is
larger here: at $\lambda=2$ it is $3.9\times10^{-3}$ for $\epsilon=3\times10^{-6}$ and
$1.3\times10^{-3}$ for $\epsilon=10^{-6}$. The small-$\lambda$ limit of the corrected form
is $1+\lambda^2(N+2-\ln2-\gamma_E)^2$, with $N=\ln(1/z_{\mathrm{end}})$. This entry is also
recorded in [`literature-discrepancies.md`](literature-discrepancies.md),
entry 10.

## Tests

The regular suite runs:

- the `PolarFieldSpace` checks in `src/field_space.rs`: inverse metric,
  Levi-Civita connection and its derivative against finite differences, and
  vanishing curvature built from the connection;
- `closed_form_gamma_and_digamma_match_known_values`, for the complex
  $\ln|\Gamma|$ and $\operatorname{Re}\psi$ helpers used by the closed forms;
- `constant_turn_spectrum_matches_the_exact_massive_solutions`, the table of
  the massive entropy field, to $10^{-3}$;
- `constant_turn_spectrum_fixes_the_logarithm_of_the_massless_solution`, the
  table of the massless entropy field, to $10^{-3}$, which also asserts that
  the printed Eq. (39) is more than 20% high.

## Next steps

The companion article gives squeezed-limit bispectra for the cubic
interactions $\sigma(D_t\zeta)^2$, $\sigma^2D_t\zeta$, and $\sigma^3$. Comparing them
with the three-point transport first requires deriving the coefficients
$\alpha$, $\beta$, and $\gamma$ of that effective action from the cubic action of a
concrete sigma model, such as the orbit above.
