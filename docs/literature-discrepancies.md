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

## Beyond the printed formulas

[DFMS] freeze the background coefficients on the initial slice. Here the
initial condition also keeps their first-order e-fold rates in the odd-momentum
blocks. Without them, the reference benchmark showed a start-depth floor of
about 2% in the final three-point function. This is an extension of the
prescription rather than a discrepancy; it is derived and tested in
[`three-point-transport.md`](three-point-transport.md).
