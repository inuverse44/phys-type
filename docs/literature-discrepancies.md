# Literature discrepancies and convention notes

This file records every place where a formula used by this repository differs
from its published source, together with notes on conventions that caused, or
nearly caused, a mistake. Each entry states the reference and equation
number, what is printed there, what the repository uses, and the evidence for
the choice. An entry does not claim the authors' codes are wrong; published
codes can be correct even when a printed formula is not.

Equation numbers refer to the arXiv versions listed below. The PDFs were
checked directly, because equation numbers derived from the TeX labels can
differ from the typeset numbering.

## References

- **[DFMS]** M. Dias, J. Frazer, D. J. Mulryne, and D. Seery, *Numerical
  evaluation of the bispectrum in multiple field inflation — the transport
  approach with code*, JCAP 12 (2016) 033,
  [arXiv:1609.00379v2](https://arxiv.org/abs/1609.00379).
- **[DEFMS]** M. Dias, J. Elliston, J. Frazer, D. Mulryne, and D. Seery, *The
  curvature perturbation at second order*, JCAP 1502 (2015) 040,
  [arXiv:1410.3491v3](https://arxiv.org/abs/1410.3491).
- **[BS]** S. Butchers and D. Seery, *Numerical evaluation of inflationary
  3-point functions on curved field space*,
  [arXiv:1803.10563](https://arxiv.org/abs/1803.10563).
- **[HKKPT]** S. C. Hyun, J. Kim, T. Kodama, S. C. Park, and T. Takahashi,
  *Nonminimally Assisted Inflation: A General Analysis*,
  [arXiv:2302.05866v2](https://arxiv.org/abs/2302.05866).
- **[HMPS]** J. Huenupi, C. Muñoz, G. A. Palma, and S. Sypsas, *Pushing the
  Primordial Frontier: Exact Linear Solutions in Multifield Inflation*,
  [arXiv:2606.18248v2](https://arxiv.org/abs/2606.18248); companion article
  *Pushing the Primordial Frontier: Cosmological Collider Signatures at
  Strong Mixing*, [arXiv:2607.14529](https://arxiv.org/abs/2607.14529).
- **[M]** J. Maldacena, *Non-Gaussian features of primordial fluctuations in
  single field inflationary models*, JHEP 05 (2003) 013,
  [arXiv:astro-ph/0210603](https://arxiv.org/abs/astro-ph/0210603).

## Discrepancies

### 1. Cross coefficients of the gauge transformation, [DFMS] Eq. (7.4b)

**Printed.** [DFMS] Eq. (7.2) writes
$\zeta(\mathbf k)=N_a\delta X^a+\tfrac12N_{ab}\delta X^a\delta X^b$, and
Eq. (7.4b) gives

$$
N_{ab}=\frac{1}{3H^2M_{\mathrm{Pl}}^2\epsilon}
\begin{pmatrix}
\dfrac{\pi_\alpha\pi_\beta}{M_{\mathrm{Pl}}^2}\Bigl[\cdots\Bigr] &
\dfrac{3}{H\epsilon}\dfrac{\pi_\alpha\pi_\beta}{M_{\mathrm{Pl}}^2}
-\dfrac{3H}{k^2}\bigl[\mathbf k_a\!\cdot\!\mathbf k_b+k_a^2\bigr]\delta_{\alpha\beta}\\[1ex]
\dfrac{3}{H\epsilon}\dfrac{\pi_\alpha\pi_\beta}{M_{\mathrm{Pl}}^2}
-\dfrac{3H}{k^2}\bigl[\mathbf k_a\!\cdot\!\mathbf k_b+k_b^2\bigr]\delta_{\alpha\beta} & 0
\end{pmatrix}.
$$

**Used here.** The same diagonal entry, but the off-diagonal entries with the
prefactor $1/(6H^2M_{\mathrm{Pl}}^2\epsilon)$.

**Reason.** The real-space source of Eq. (7.4b) is [DFMS] Eq. (7.1), identical
to $\zeta^{\mathrm{simple}}_1+\zeta^{\mathrm{simple}}_2$ of [DEFMS]
Eqs. (2.34a)–(2.34b). There, the $\delta\phi\,\delta\dot\phi$ and gradient
terms carry the overall factor $1/(6H^2M_{\mathrm{Pl}}^2\epsilon)$ and appear
once. In $\tfrac12N_{ab}X^aX^b$, the field--momentum pair is counted twice,
once as $N_{Q\Pi}$ and once as $N_{\Pi Q}$. Each off-diagonal entry must
therefore be half of the real-space coefficient, giving
$1/(6H^2M_{\mathrm{Pl}}^2\epsilon)$. The $\delta\phi\,\delta\phi$ entry already
has the factor $1/2$ built in and agrees with Eq. (7.4b).

**Evidence.** On the straight-trajectory equal-mass double-quadratic
benchmark of [`curvature-projection.md`](curvature-projection.md), where the
result must reduce to single-field inflation:

| Cross prefactor | equilateral $f_{\mathrm{NL}}$ at 3 and 4 e-folds after exit | Maldacena [M] |
| --- | --- | --- |
| $1/(3H^2M_{\mathrm{Pl}}^2\epsilon)$, as printed | 0.569, 0.729 | 0.190 |
| $1/(6H^2M_{\mathrm{Pl}}^2\epsilon)$, used here | 0.210, 0.219 | 0.190 |

This calibration run used a pivot close to the end of inflation
($\epsilon\approx0.08$) and a shallow start, so the $1/6$ row is not yet
accurate. With the pivot moved to $\epsilon\approx0.0125$, the $1/6$
coefficients conserve $f_{\mathrm{NL}}$ on superhorizon scales to 0.4%, agree
with [M] to next order in $\epsilon$, and satisfy the squeezed consistency
relation to 0.33% at $k_3/k_1=0.05$. The printed coefficient is ruled out by
both the conservation and the single-field limit.

**Where.** `src/curvature.rs`, `curvature_gauge_coefficients`; tests
`cross_coefficients_attach_the_gradient_term_to_the_field_input`,
`single_field_f_nl_is_conserved_and_matches_maldacena`, and
`squeezed_single_field_bispectrum_satisfies_the_consistency_relation`.

### 2. $B$ contribution to the two-momentum initial condition, [DFMS] Eq. (6.7c)

**Printed.** The initial condition for "$a,b$ momenta, $c$ field" ends with

$$
\cdots+k_1^2k_3^2\Bigl[\frac{a^2H}{2}B^{\alpha\beta\gamma}k_3+\text{5 perms}\Bigr].
$$

**Used here.** The $B$ contribution obtained from the tree-level in-in
integral,

$$
\alpha^{\Pi\Pi Q}\supset-\frac{H}{4\kappa k_ta^4}\Bigl(k_1^2\mathcal B_1+k_2^2\mathcal B_2
+\bigl[k_3(k_1+k_2)-k_1k_2\bigr]\mathcal B_3\Bigr),
$$

in the notation of [`three-point-transport.md`](three-point-transport.md),
after the momentum normalization of entry 3 is accounted for.

**Reason.** The block with momenta on legs 1 and 2 must be symmetric under
exchanging those two legs together with their indices. The prefactor
$k_1^2k_3^2$ outside the permutation sum is not symmetric under
$1\leftrightarrow2$. Every other block and every other contribution
(fast $A$, slow $A$, $C$) of [DFMS] Eqs. (6.7a)–(6.7d) agrees exactly with the
independent derivation. The residual difference in this term is proportional
to $(k_2-k_3)(k_2+k_3)$.

**Evidence.** `docs/derivations/three_point_initial_condition.py` computes the
in-in integral symbolically and prints the comparison. The result used here
also solves the canonical transport equation on de Sitter space to the
retained order (Rust test
`subhorizon_initial_condition_solves_de_sitter_transport_to_retained_order`).

### 3. Momentum normalization of the three-point initial conditions, [DFMS] Eqs. (6.7b)–(6.7d)

**Printed.** [DFMS] Eqs. (4.11)–(4.13) define $\delta\pi^\alpha$ as the canonical
momentum rescaled by $a^{-3}$, conjugate to cosmic time, so that
$\delta\pi=\delta\dot\phi$ at linear order. Eq. (6.2a) uses that normalization:
its $\pi\pi$ entry is $(k^2/a)/(2a^3k)$.

**Observation.** With this definition, the stripped blocks with $n$ momentum
slots have dimension $M^{-3+n}$. Eqs. (6.7b)–(6.7d) instead all have dimension
$M^{-3}$. They agree term by term with the independent derivation only after
each momentum slot is divided by $H$, that is, if $\delta\pi$ is read as
$d\delta\phi/dN$. Eqs. (6.7c) and (6.7d) carry an explicit $1/H(t_{\mathrm{init}})^2$
consistent with that reading.

**Used here.** The cosmic-time canonical momentum of Eq. (4.13). The
comparison in the derivation script multiplies the published blocks by
$H^{n}$ before comparing.

### 4. Momentum arguments of permuted cubic tensors, [DFMS] Eq. (5.17) — convention note

This is not an error in the reference. [DFMS] Eq. (5.17) lists the eight
blocks of $\tilde u^a{}_{bc}$ with index permutations such as
$C_{\gamma\beta}{}^\alpha$. The text immediately after it states that the
momenta are attached to the indices: a tensor written with indices
$\alpha,\beta,\gamma$ is evaluated at $\mathbf k_\alpha,\mathbf k_\beta,\mathbf
k_\gamma$ through the assignment $\mathbf k_a\to\mathbf k_\alpha$, and so on. A
permuted index therefore also permutes its momentum.

An earlier revision of this repository permuted only the array indices of the
`q_qp`, `p_pq`, and `p_pp` blocks and left the momentum arguments in place.
That agrees with Eq. (5.17) only for equal input wavenumbers, so equilateral
tests could not detect it. The de Sitter consistency test found it. The
corrected table with explicit momentum arguments is in
[`cubic-kernel-conventions.md`](cubic-kernel-conventions.md). [BS] use the same
cubic tensors and should be read with the same attachment rule.

### 5. Operator order of the source covariances, [DFMS] Eq. (5.14) — convention note

This is not an error in the reference. [DFMS] Eq. (5.14) contracts the source
with $\Sigma^{db}(k_b)\Sigma^{ec}(k_c)$ for the first output leg,
$\Sigma^{ad}(k_a)\Sigma^{ec}(k_c)$ for the second, and
$\Sigma^{ad}(k_a)\Sigma^{be}(k_b)$ for the third. The text notes that the
ordering matters because $\Sigma$ is not symmetric when it has an imaginary
part. (The last factor is printed with the argument $k_c$; with index $b$ it
can only mean $k_b$.)

An earlier revision of this repository transposed the covariance on legs to
the left of the insertion. The corrected source, and the same ordering in the
curvature projection, follow Eq. (5.14).

### 6. Slow-roll parameter of loop inflation, [HKKPT] Eq. (4.2)

**Printed.** For $V_J=\Lambda^4(1+\alpha\log\phi)$,
$\epsilon^{(0)}=\alpha^2/[2\phi(1+\alpha\log\phi)^2]$.

**Used here.** $\epsilon^{(0)}=\tfrac12(V_J'/V_J)^2=\alpha^2/[2\phi^2(1+\alpha\log\phi)^2]$.

**Reason.** $V_J'=\alpha\Lambda^4/\phi$. The end-of-inflation value printed
in the same section, $\phi_e=[\sqrt2\,W_0(e^{1/\alpha}/\sqrt2)]^{-1}$, solves
$\epsilon^{(0)}=1$ only with $\phi^2$. The accompanying
$\eta^{(0)}=-\alpha/[\phi^2(1+\alpha\log\phi)]$ is correct. With $\phi^2$,
Eq. (3.1) reproduces the numerical $n_s$ to $5\times10^{-4}$; see
[`nonminimal-assisted-inflation.md`](nonminimal-assisted-inflation.md).

### 7. Sign convention of $f_{\mathrm{NL}}$, [HKKPT] Eq. (3.5) — convention note

The printed left-hand side is $-\tfrac65f^{(\mathrm{local})}_{\mathrm{NL}}$, and the
single-field limit of the right-hand side is $2\epsilon-\eta$. The reference
therefore uses the Lyth--Rodriguez sign,
$\zeta=\zeta_g-\tfrac35f_{\mathrm{NL}}\zeta_g^2$. That is the opposite of the
Maldacena/Planck sign used by this repository
($f_{\mathrm{NL}}\to+\tfrac{5}{12}(1-n_s)$ in the squeezed single-field limit).
The reference compares its values with the Planck bound
$-11.1<f_{\mathrm{NL}}<9.3$, which is stated in the Planck sign. For the small
values found there the bound is not affected. Comparisons here multiply the
printed $f_{\mathrm{NL}}$ by $-1$.

### 8. The slow-roll e-fold relation, [HKKPT] Eqs. (3.9) and (3.12) — approximation note

This is not an error. [HKKPT] fixes $s_e$ from $s_*$ through the slow-roll
relation (3.9), or its small-coupling form (3.12). On the exact loop-inflation
background with $\alpha=0.1$, $s_*=10^{-3}$, and $\xi_2=0.022$, the numerical
$s_*$ and $s_e$ satisfy Eq. (3.9) with 58.3 rather than 60 e-folds. Because
$s_e\propto e^{4\xi_2N}$, this 3% difference changes $s_e$ by 16%. It moves the
coupling that gives $n_s=0.965$ from $\xi_2\approx0.0222$ to $0.0228$. Given the
numerical $s_e$, Eq. (3.1) itself agrees with the exact tilt to
$5\times10^{-4}$ for $0.021\le\xi_2\le0.025$. The other three branches of
Sec. 4.1 behave the same way: the $n_s=0.965$ coupling moves up by
1.6--3.2%, and it stays inside the quoted ranges except at the upper edge
for $(\alpha,s_*)=(1,10^{-3})$. For the Class II power-law and hybrid examples
of Secs. 4.2 and 4.3 the shift is at most 1.6% in either direction, and every
numerical coupling lies inside the quoted ranges, because $\epsilon^{(\phi)}$ is
small and nearly constant there.

### 9. Crossing-time statistics in [HKKPT] Eq. (3.5) — approximation note

This is not an error. The slow-roll $\delta N$ expression (3.5) assumes massless,
uncorrelated field perturbations at horizon crossing. For loop inflation with
$\xi_2=0.022$ and $s_*=10^{-3}$, the assistant direction is tachyonic
($\eta^{\sigma\sigma}\approx-4\xi_2$). Eq. (3.5) and a numerical crossing-time
$\delta N$ then agree with each other ($-0.216$ and $-0.199$ in the Planck sign)
but underestimate $|f_{\mathrm{NL}}|$ by about 20%. A phase-space $\delta N$ from a
later slice, which uses the evolved covariance, reproduces the tree-level
transport value $-0.2514$ to 0.05%. The conclusion of the reference that
$|f_{\mathrm{NL}}|\ll1$ is unaffected. The Class II examples with $s_*=10^{-3}$
show the same effect with a deficit of about 8%: for the power-law model with
$\lambda=0.01$ and $\xi_2=0.011$, the crossing-time $\delta N$ gives $-1.51$ and the
transport $-1.635$; for the hybrid model with $\mu=10$ and $\xi_2=0.0108$, $-12.3$
and $-13.36$.

### 10. Logarithm of the massless constant-turn spectrum, [HMPS] Eq. (39)

[HMPS] Eq. (39) gives the superhorizon $\Delta_\zeta/\Delta_0$ for a massless entropy
field with the growing factor
$\left[2+\ln(2/z_{\mathrm{end}})-2\gamma_E-\operatorname{Re}\psi(i\lambda/2)\right]^2$. The
repository uses $\ln(1/2z_{\mathrm{end}})$ in place of $\ln(2/z_{\mathrm{end}})$. The
evidence:

- the $\mathcal H_w$ of [HMPS] Eq. (38), from which (39) follows, contains
  $-2\ln(2z)$, that is $+2\ln(1/2z)$;
- an independent integration of [HMPS] Eqs. (8)–(9)
  ([`derivations/constant_turn_exact_solutions.py`](derivations/constant_turn_exact_solutions.py))
  agrees with the corrected form to $5\times10^{-4}$ for $0.02\le\lambda\le3$, the
  residual being the $\mathcal O(z_0^{-1})$ of its initial condition; the printed
  form is 1–46% high;
- the covariant mode solver on a constant-turn orbit agrees with the
  corrected form to $4\times10^{-4}$ for $\lambda=0.5$, 1, and 2
  ([`constant-turn-exact-solutions.md`](constant-turn-exact-solutions.md)).

The corrected form has the small-$\lambda$ limit
$1+\lambda^2(N+2-\ln2-\gamma_E)^2$ with $N=\ln(1/z_{\mathrm{end}})$. The massive results,
Eqs. (43)–(45), are unaffected. The equation was read from the typeset PDF.

## Agreements checked

Recording these makes clear which published results have been independently
reproduced here:

- [DFMS] Eqs. (6.7a), (6.7b), and (6.7d), every contribution, and (6.7c)
  except entry 2, after the normalization of entry 3;
- [DFMS] Eq. (6.2a)–(6.2b) subhorizon two-point initial values, reproduced by
  the Bunch--Davies mode state to the orders shown there;
- [DFMS] Eq. (7.4a) and the $QQ$ entry of Eq. (7.4b);
- [DFMS] Eq. (7.6b), the structure of the bispectrum projection, with the
  ordering note of entry 5;
- the squeezed-limit consistency relation of [M] and its leading-order
  bispectrum, as recorded in
  [`curvature-projection.md`](curvature-projection.md);
- [DFMS] Fig. 7, the equilateral $f_{\mathrm{NL}}$ of the axion-quartic model of
  Sec. 9.1, reproduced at nine plotted samples to the precision of the figure,
  with the samples recovered from its vector path by
  `docs/derivations/digitize_dfms_figure7.py`. The agreement also confirms
  that the printed $1/3$ of entry 1 is not what their code uses.

- [HKKPT] Eq. (3.1), given the numerical $s_e$ (entry 8), against the mode
  solver and a numerical $\delta N$ for loop inflation; and Eq. (3.5), with the
  sign of entry 7, against a numerical $\delta N$ to 9%.
- [HMPS] Eqs. (43)–(45), the massive constant-turn spectra, against the mode
  solver on a circular orbit to $3.5\times10^{-4}$ for $0.5\le\lambda\le2$ and
  $1.25\le\mu^2/H^2\le9$.
- [HKKPT] Eq. (3.1) for the Class II power-law and hybrid examples, given the
  numerical end surface and $s_e$, against the mode solver to
  $1.5\times10^{-3}$.

## Beyond the printed formulas

[DFMS] freeze the background coefficients on the initial slice. Here the
initial condition also keeps their first-order e-fold rates in the odd-momentum
blocks. Without them, the reference benchmark showed a start-depth floor of
about 2% in the final three-point function. This is an extension of the
prescription rather than a discrepancy; it is derived and tested in
[`three-point-transport.md`](three-point-transport.md).
