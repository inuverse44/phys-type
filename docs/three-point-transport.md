# Three-point transport

## Purpose and current boundary

This document describes P3d: evolution of the connected three-point function
of scalar perturbations. It begins where the validated P3c two-point transport
ends. The immediate input is therefore a background trajectory and three
complex phase-space covariances, not a curvature bispectrum or an
$f_{\mathrm{NL}}$ value.

P3d is in progress. The following pieces are implemented:

- the complete dimensioned three-point state;
- Fourier-leg permutation and complex-conjugation operations;
- homogeneous evolution by the linear flow $u^a{}_b$;
- a dimensioned container for the quadratic flow $u^a{}_{bc}$;
- a dimensioned canonical $A/B/C$ container and its ordered conversion to all
  eight e-fold-time $u^a{}_{bc}$ blocks;
- model-derived evaluation of the fully symmetrized $A_{IJK}$ and the
  symmetrized $B_{IJK}$ and $C_{IJK}$ coefficients, including every raised
  form required by the block map;
- one model-to-kernel entry point that evaluates $A/B/C$ and feeds the eight
  canonical e-fold-time blocks;
- an output-leg-aware vertex set with input orders $(2,3)$, $(1,3)$, and
  $(1,2)$, plus contraction of the three operator-ordered $u\Sigma\Sigma$
  source terms;
- explicit linear two-point and tree-level nonlinear three-point conversion
  between the canonical $(Q,\Pi)$ and velocity $(Q,V)$ bases;
- the tree-level Bunch--Davies initial value of all eight canonical
  `QQQ`--`PPP` blocks on one common subhorizon slice.

The following pieces are not implemented:

- one adaptive system that evolves the sourced state through a physical
  benchmark;
- general-triangle convergence tests of that system in tolerance and start
  depth;
- nonlinear projection to the curvature bispectrum.

The canonical-momentum type boundary required by the first missing item is
implemented: canonical and velocity covariances, correlations, and kernels
cannot be mixed in the typed source contraction. Initial linear covariances
can be relabelled only through an explicitly named linear-order conversion;
the tree-level output conversion adds the required local $B\Sigma\Sigma$ and
$C\Sigma\Sigma$ terms before returning a velocity-marked correlation.

Therefore the current API is transport infrastructure, not yet a primordial
bispectrum prediction.

## Notation and index domains

There are three different kinds of indices:

| Symbol | Meaning |
| --- | --- |
| $I,J,K,\ldots$ | field-space coordinates, transformed with $G_{IJ}$ |
| $a,b,c,\ldots$ | phase-space indices; each means a slot $Q$ or $P$ together with a field index |
| $r,s\in\{1,2,3\}$ | Fourier-leg labels associated with $\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3$ |

Bold $\mathbf{k}_r$ denotes a wavevector and $k_r=|\mathbf{k}_r|$ its positive
magnitude. Translation invariance imposes

$$
\mathbf{k}_1+\mathbf{k}_2+\mathbf{k}_3=0,
$$

which is why the three magnitudes form a triangle. A Fourier-leg permutation
moves both $\mathbf{k}_r$ and the phase-space slot attached to that leg. It is
not a permutation of field-space indices.

The P3c phase-space convention is

$$
X^a=(Q^I,P^I),\qquad Q^I=U^I,\qquad P^I=D_tU^I,
$$

where $U^I$ is the covariant field perturbation and $D_t$ is the field-space
covariant cosmic-time derivative. The e-fold derivative is

$$
D_N\equiv H^{-1}D_t.
$$

The linear mode functions used by P3c have dimensions $M^{-1/2}$ for the
field amplitude and $M^{1/2}$ for its cosmic-time derivative. These are not
the dimensions of the Fourier-space field operators themselves. With

$$
Q(\mathbf{k})=\int d^3x\,e^{-i\mathbf{k}\cdot\mathbf{x}}Q(\mathbf{x}),
$$

the operator dimensions are $[Q(\mathbf{k})]=M^{-2}$ and
$[P(\mathbf{k})]=M^{-1}$. The annihilation operator carries the compensating
Fourier-volume dimension that makes a linear mode function proportional to
$1/\sqrt{2k}$. This distinction matters beginning at the stripped three-point
function.

Here the symbol $P$ and the Rust marker `MomentumPerturbation` mean a
covariant velocity. They must not be identified beyond linear order with the
rescaled canonical momentum used in the standard cubic-Hamiltonian
derivation. The exact relation and the eight-block canonical kernel are given
in [`cubic-kernel-conventions.md`](cubic-kernel-conventions.md).

## Correlation functions

After removing the momentum-conserving delta function, the equal-time
two-point covariance and connected three-point function are

$$
\Sigma^{ab}(k)
  =\left\langle X^a(\mathbf{k})X^b(-\mathbf{k})\right\rangle,
$$

$$
\alpha^{abc}(k_1,k_2,k_3)
  =\left\langle
      X^a(\mathbf{k}_1)X^b(\mathbf{k}_2)X^c(\mathbf{k}_3)
    \right\rangle_{\!c}.
$$

The subscript $c$ means the connected part. The Rust type
`ThreePointPhaseSpaceCorrelation<N>` represents $\alpha^{abc}$ using eight
blocks:

| Block | Meaning | Dimension |
| --- | --- | --- |
| `QQQ` | $\langle Q^I Q^J Q^K\rangle_c$ | $M^{-3}$ |
| `QQP`, `QPQ`, `PQQ` | one momentum slot | $M^{-2}$ |
| `QPP`, `PQP`, `PPQ` | two momentum slots | $M^{-1}$ |
| `PPP` | $\langle P^I P^J P^K\rangle_c$ | $M^0$ |

The slot markers and physical dimension are both part of each Rust type. For
example, a numerical array stored as `QQP` cannot be assigned to `QQQ`.

These dimensions follow after removing one
$(2\pi)^3\delta^3(\mathbf{k}_1+\mathbf{k}_2+\mathbf{k}_3)$, whose dimension is
$M^{-3}$. A stripped $n$-point function does not in general have the dimension
of a naive product of $n$ linear mode functions. For $n=2$ the two happen to
agree, which is why the error is invisible in the P3c covariance.

## Fourier-leg operations

The API supplies the two transpositions $(12)$ and $(23)$. Each operation
moves all of the following together:

1. the appropriate component-array axis;
2. the compile-time $Q/P$ slot marker;
3. the containing block among `QQQ` through `PPP`.

Each transposition squares to the identity, and their compositions generate
all six permutations in $S_3$. Complex conjugation is separate because it
changes quantum phase information without moving a Fourier leg. These
operations make it possible to define permutation and reality residuals later;
they do not force an arbitrary input state to obey those identities.

## Nonlinear phase-space flow

Expanding the perturbation equations to quadratic order gives

$$
D_NX^a(\mathbf{k})
=u^a{}_b(k)X^b(\mathbf{k})
+\frac12\int\!\frac{d^3q}{(2\pi)^3}\,
  u^a{}_{bc}(\mathbf{k};\mathbf{q},\mathbf{k}-\mathbf{q})
  X^b(\mathbf{q})X^c(\mathbf{k}-\mathbf{q})+\cdots.
$$

Here $u^a{}_b$ is the linear phase-space flow and $u^a{}_{bc}$ is the
quadratic flow in the declared phase-space variables. The commonly quoted
cubic-Hamiltonian coefficients are initially written for
$X_\Pi=(Q,\Pi)$, whereas this document uses $X=(Q,P)$ with $P=D_tQ$.
They agree linearly but require a nonlinear variable transformation at this
order. The factor $1/2$ belongs to the quadratic expansion. At tree level, the
two Wick contractions of the quadratic term cancel this factor in the
three-point source equation below.

### Homogeneous part

`ThreePointHomogeneousTransportSystem` implements

$$
\begin{aligned}
D_N\alpha^{abc}
={}&u^a{}_d(k_1)\alpha^{dbc}
   +u^b{}_d(k_2)\alpha^{adc}
   +u^c{}_d(k_3)\alpha^{abd}.
\end{aligned}
$$

The three terms propagate legs 1, 2, and 3 respectively. Every Runge--Kutta
stage reads the shared P2 background trajectory. The linear flow contains the
connection, Hubble friction, covariant perturbation mass matrix, and the
leg-dependent $k_r^2/a^2$ term already validated in P3c.

The word `Homogeneous` is intentional: if $\alpha=0$, this equation alone
keeps $\alpha=0$.

### Quadratic source

The complete tree-level source is

$$
\begin{aligned}
S^{abc}={}&
u^a{}_{de}(\mathbf{k}_1;\mathbf{k}_2,\mathbf{k}_3)
  \Sigma^{db}(k_2)\Sigma^{ec}(k_3)\\
&+u^b{}_{de}(\mathbf{k}_2;\mathbf{k}_1,\mathbf{k}_3)
  \Sigma^{ad}(k_1)\Sigma^{ec}(k_3)\\
&+u^c{}_{de}(\mathbf{k}_3;\mathbf{k}_1,\mathbf{k}_2)
  \Sigma^{ad}(k_1)\Sigma^{be}(k_2).
\end{aligned}
$$

The covariance order follows the operator order in
$\langle X^aX^bX^c\rangle$. When the quadratic insertion replaces
$X^b$, the external $X^a$ stands to its left, so the Wick contraction is
$\langle X^aX^d\rangle=\Sigma^{ad}$; an external operator to the right
contracts as $\Sigma^{db}$. Because
$\Sigma^{ad}=(\Sigma^{da})^*$ and the imaginary part is the equal-time
commutator, transposing a covariance on a left leg changes the source. This
matches the three-point transport equation in Sec. 5 of Dias, Frazer,
Mulryne, and Seery, [arXiv:1609.00379](https://arxiv.org/abs/1609.00379). An earlier revision of
this repository used $\Sigma^{da}$ on left legs; the de Sitter consistency
test below failed with an order-one residual until the ordering was corrected.

The sourced transport equation will be

$$
D_N\alpha^{abc}
=\left.D_N\alpha^{abc}\right|_{\mathrm{hom}}+S^{abc}.
$$

`QuadraticPhaseSpaceKernel<N>` stores the eight ordered $Q/P$ blocks of
$u^a{}_{bc}$. The convolution measure contributes $M^3$, so dimensional
consistency requires

$$
[u^a{}_{bc}]
=[X^a(\mathbf{k})]-[X^b(\mathbf{q})]-[X^c(\mathbf{k}-\mathbf{q})]-M^3.
$$

Therefore:

| Output | Inputs `QQ` | Inputs `QP` | Inputs `PQ` | Inputs `PP` |
| --- | ---: | ---: | ---: | ---: |
| $Q$ | $M^{-1}$ | $M^{-2}$ | $M^{-2}$ | $M^{-3}$ |
| $P$ | $M^0$ | $M^{-1}$ | $M^{-1}$ | $M^{-2}$ |

`three_point_source_from_kernels` evaluates the displayed expression for
$S^{abc}$. Its three kernel arguments have output legs 1, 2, and 3. The two
input legs are stored in ascending order among the remaining legs:

| Kernel output leg | Ordered input legs |
| --- | --- |
| 1 | $(2,3)$ |
| 2 | $(1,3)$ |
| 3 | $(1,2)$ |

Making the middle row explicit prevents confusing ascending order $(1,3)$
with a cyclic convention $(3,1)$.

## Physical kernel boundary and remaining validation

The contraction accepts correctly dimensioned numerical kernel blocks. The
library manufactures all $A_{IJK}$, $B_{IJK}$, and $C_{IJK}$
contributions from a model, background state, geometry, and momentum triangle.
The canonical route uses a distinct typed momentum $\Pi^I$; it is not silently
identified with $P^I=D_tQ^I$.

`canonical_cubic_hamiltonian_vertices` evaluates the coefficient set once for
each output leg, moving every field index with its attached Fourier momentum.
This supplies exactly the three orderings required by the source table above.
`velocity_three_point_from_canonical_tree` uses those same ordered vertices to
apply the nonlinear basis change after canonical evolution.

The $A$ evaluator includes the symmetrized covariant potential derivative
$V_{;(IJK)}$, field-space curvature, and the required covariant curvature
derivative. These are the reasons P3b implemented the third-order geometry
first.

This boundary is important: compile-time dimensions can reject a kernel in the
wrong physical slot, but they cannot prove that a coefficient, sign, momentum
argument, or curvature contraction was derived correctly. Those properties
require analytic limits, permutation tests, and comparison with an independent
calculation.

The full formulas for $A_{IJK}$, $B_{IJK}$, $C_{IJK}$, their simultaneous
index/momentum symmetrizations, the canonical eight-block $u^a{}_{bc}$ table,
and the resolved $P$--$\Pi$ type boundary are recorded in
[`cubic-kernel-conventions.md`](cubic-kernel-conventions.md). The remaining
scientific gates are the quantum initial condition, combined sourced
integration, and an independent numerical benchmark.

## Initial conditions

### Normalization and slice

The initial condition populates the stripped, canonical-basis correlation

$$
\left\langle X_\Pi^a(\mathbf{k}_1)X_\Pi^b(\mathbf{k}_2)X_\Pi^c(\mathbf{k}_3)
\right\rangle
=(2\pi)^3\delta^3(\mathbf{k}_1+\mathbf{k}_2+\mathbf{k}_3)\,
\alpha_\Pi^{abc}(k_1,k_2,k_3),
$$

with $X_\Pi=(Q^I,\Pi^I)$, comoving wavenumbers, and the block dimensions
listed above. All three legs start on one common slice. With
$x_r\equiv k_r/(aH)$, the slice is accepted only if

$$
\min_r x_r\ge x_{\mathrm{init}}>1,
$$

so the shortest leg decides whether the whole triangle may start.
`canonical_three_point_initial_condition` rejects a slice that fails this
condition and returns the blocks in the triangle's canonical leg order.

The initial value is obtained from the tree-level in-in formula,

$$
\alpha_\Pi^{abc}(\tau)
=i\int_{-\infty(1\mp i\epsilon)}^{\tau}dt'\,
\left\langle\left[H_{\mathrm{int}}(t'),
X_\Pi^a X_\Pi^b X_\Pi^c(\tau)\right]\right\rangle',
\qquad
H_{\mathrm{int}}=-\frac{a^3}{2}\left(A\,QQQ+B\,QQ\Pi+C\,\Pi\Pi Q\right),
$$

where the prime removes $(2\pi)^3\delta^3$. In the interaction picture
$\Pi^I=D_tQ^I$, while the Heisenberg operator on the left is the canonical
momentum. The result is therefore a canonical-basis correlation; it must not
be relabelled as a velocity correlation.

### Approximation

Deep inside the horizon the $k^2/a^2$ term dominates the mass matrix, so the
propagator is that of massless fields in de Sitter space,

$$
Q_k(\tau)=\frac{H}{\sqrt{2k^3}}(1+ik\tau)e^{-ik\tau},\qquad
\Pi_k(\tau)=-\frac{H^2k^2\tau^2}{\sqrt{2k^3}}e^{-ik\tau},\qquad
\tau=-\frac{1}{aH},
$$

in a parallel-transported orthonormal frame, so the field indices are raised
with $G^{IJ}$ on the slice. Background coefficients are frozen at the slice,
except for the explicit $a^{-2}$ of the gradient term

$$
A^{\mathrm{fast}}_{IJK}
=\operatorname{Sym}\left[\frac{\dot\phi_IG_{JK}\,
\mathbf{k}_J\!\cdot\!\mathbf{k}_K}{2HM_{\mathrm{Pl}}^2a^2}\right].
$$

The remaining terms of $A$ are called slow,
$A^{\mathrm{slow}}=A-A^{\mathrm{fast}}$.

For every block, and separately for each of the fast-$A$, slow-$A$, $B$, and
$C$ contributions, the implementation keeps the highest power of $x_r$ that
survives. For blocks with an even number of momentum slots, this is the naive
leading order. For blocks with one or three momentum slots, the naive leading
term is purely imaginary inside the $\operatorname{Im}$ of the in-in formula
and cancels. Their first surviving order, one power of $aH/k$ lower, is
retained instead, including the $-H$ part of $\Pi_k$ and the $1/k_t^2$ term of
the time integral. Setting these blocks to zero would violate the transport
equation below at leading order.

Massless propagators and frozen background coefficients neglect corrections
proportional to $m^2/H^2$ and to background-variation rates $\delta$, such
as $\epsilon$. A Taylor expansion of a coefficient around the slice brings
in an extra factor $i\,\delta\,aH/k_t$ relative to the retained integral.
In the even blocks this factor is imaginary inside the $\operatorname{Im}$
and first contributes at order $(aH/k)^2$. In the odd blocks the retained
term is itself one power of $aH/k$ below its natural size, so the neglected
correction can be a fraction of order $\delta$ of the retained value. These
statements were not derived for the mass correction separately.

A simple estimate indicates why this matters. An error in the initial state
excites homogeneous solutions, which fall like $a^{-3}$ inside the horizon,
whereas the sourced solution falls like $a^{-4}$. An initial relative error
$\varepsilon$ can therefore grow to about $\varepsilon\,x_{\mathrm{init}}$
by horizon exit. The retained truncation error, of relative order
$x_{\mathrm{init}}^{-2}$, then leaves a final error of order
$x_{\mathrm{init}}^{-1}$ that vanishes with increasing start depth, while an
odd-block error of order $\delta/x_{\mathrm{init}}$ could leave a final error
of order $\delta$ that does not. This is an estimate, not a result: the
start-depth convergence test of the sourced evolution must measure it, as the
DFMS numerical examples do.

### Result

Define $k_t=k_1+k_2+k_3$, $K^2=k_1k_2+k_1k_3+k_2k_3$,
$\kappa=k_1k_2k_3$, and $\mathcal{N}=1/(4\kappa k_t)$. The fully raised
tensors entering the result are

| Symbol | Meaning | Dimension |
| --- | --- | --- |
| $\mathcal{A}_f^{IJK}$ | $a^2A_{\mathrm{fast}}^{IJK}$ | $M^1$ |
| $\mathcal{A}_s^{IJK}$ | $A_{\mathrm{slow}}^{IJK}$ | $M^1$ |
| $\mathcal{B}_r^{IJK}$ | $B$ with its $\Pi$ slot on leg $r$ | $M^0$ |
| $\mathcal{C}_r^{IJK}$ | $C$ with its $Q$ slot on leg $r$ | $M^{-1}$ |

For example, $\mathcal{B}_3^{IJK}=B^{IJK}(k_1,k_2,k_3)$ and
$\mathcal{B}_1^{IJK}=B^{JKI}(k_2,k_3,k_1)$. They are three evaluations with
different Fourier-leg assignments, not index permutations of one array.
`canonical_subhorizon_initial_vertices` constructs them from the output-leg
vertices and stores every tensor as `[I][J][K]` with $I,J,K$ on legs
$1,2,3$. All components below carry these indices and the result is real.

For `QQQ`,

$$
\alpha^{QQQ}=\mathcal{N}\left[
\frac{3\mathcal{A}_f}{a^4}+\frac{3\mathcal{A}_s}{a^2}
-\frac{H}{a^2}\sum_r\mathcal{B}_r
-\frac{1}{a^4}\sum_r\frac{\kappa}{k_r}\mathcal{C}_r\right].
$$

With one momentum slot on leg $p$, and the field legs $u,v$,

$$
\begin{aligned}
\alpha^{(\Pi\text{ on }p)}=\mathcal{N}\Bigl[&
-\frac{3H(k_t+k_p)}{k_ta^4}\mathcal{A}_f
-\frac{3H(k_t-k_p)}{k_ta^2}\mathcal{A}_s
+\frac{k_p}{a^4}\sum_sk_s\mathcal{B}_s\\
&-\frac{H}{k_ta^4}\Bigl(
\bigl[k_p^2(k_u+k_v)+k_p(k_u^2+k_v^2)-k_uk_v(k_u+k_v)\bigr]\mathcal{C}_p\\
&\qquad\qquad+k_p^2(k_p+k_u)\mathcal{C}_u+k_p^2(k_p+k_v)\mathcal{C}_v
\Bigr)\Bigr].
\end{aligned}
$$

With two momentum slots on legs $u,v$ and the field leg $q$,

$$
\begin{aligned}
\alpha^{(Q\text{ on }q)}=\mathcal{N}\Bigl[&
-\frac{3k_uk_v}{a^6}\mathcal{A}_f-\frac{3k_uk_v}{a^4}\mathcal{A}_s
-\frac{H}{a^4}\Bigl(k_u^2\mathcal{B}_u+k_v^2\mathcal{B}_v
+\bigl[k_q(k_u+k_v)-k_uk_v\bigr]\mathcal{B}_q\Bigr)\\
&+\frac{1}{a^6}\Bigl(\kappa k_v\mathcal{C}_u+\kappa k_u\mathcal{C}_v
+k_u^2k_v^2\mathcal{C}_q\Bigr)\Bigr].
\end{aligned}
$$

For `PPP`,

$$
\begin{aligned}
\alpha^{PPP}=\mathcal{N}\Bigl[&
\frac{3H}{a^6}\Bigl(K^2+\frac{\kappa}{k_t}\Bigr)\mathcal{A}_f
+\frac{3H}{a^4}\Bigl(K^2-\frac{\kappa}{k_t}\Bigr)\mathcal{A}_s
-\frac{\kappa}{a^6}\sum_rk_r\mathcal{B}_r\\
&-\frac{H}{k_ta^6}\sum_r k_u^2k_v^2(k_t+k_r)\mathcal{C}_r\Bigr],
\end{aligned}
$$

where in the last sum $u,v$ are the two legs other than $r$.
`canonical_three_point_from_subhorizon_vertices` implements these formulas
with typed intermediate dimensions; the compiler checks that every bracket has
the dimension required by its block.

### Independent checks

The script
[`derivations/three_point_initial_condition.py`](derivations/three_point_initial_condition.py)
performs the in-in integral symbolically and verifies three statements:

1. the closed forms above reproduce the in-in result for every block and
   every contribution;
2. after converting their momenta, which are e-fold derivatives $\Pi/H$,
   the results agree with the four $B_{\Re}^{abc}(t_{\mathrm{init}})$
   results in Sec. 6 of Dias, Frazer, Mulryne, and Seery for all blocks and
   contributions, except the $B$ part of their two-momentum result. Its second prefactor $k_1^2k_3^2$ breaks the exchange of the two
   momentum legs that `PPQ` must obey; the result above obeys that symmetry;
3. with the operator-ordered source, the initial condition solves the
   canonical transport equation on massless de Sitter space. The residual,
   relative to the magnitudes of the terms in each block equation, is about
   $3\times10^{-6}$ at $x\sim10^3$ and falls as $x^{-2}$.

The Rust test
`subhorizon_initial_condition_solves_de_sitter_transport_to_retained_order`
repeats the third check with the library's own kernel map, source
contraction, and finite differences in $N$. It exposed two defects in the
already accepted source layer: the covariance ordering described above, and
input exchanges in the cubic kernel map that moved field indices without their
momenta. Both are fixed; the second is described in
[`cubic-kernel-conventions.md`](cubic-kernel-conventions.md).

## Verification already present

- all eight $\alpha$ blocks have compile-time slot and dimension checks;
- canonical triangle dot products carry dimension $M^2$, are symmetric, and
  reproduce equilateral and folded transverse-factor limits;
- the model-derived $C_{IJK}$ evaluator preserves physical dimensions, moves
  field indices with their Fourier legs during symmetrization, rejects invalid
  $H$ and $M_{\mathrm{Pl}}$, and reproduces the flat single-field limit;
- addition and multiplication by dimensionless Runge--Kutta coefficients
  preserve every block type;
- $(12)^2=(23)^2=1$ and complex conjugation squares to the identity;
- the `QQQ` homogeneous derivative contains the expected
  $(\mathrm{PQQ}+\mathrm{QPQ}+\mathrm{QQP})/H$ terms;
- homogeneous evolution commutes with an equal-leg exchange for an
  equilateral triangle in every block and field component;
- a scalar complex regression independently reproduces all three ordered
  $u\Sigma\Sigma$ terms;
- a nonlinear scalar regression reproduces the $B\Sigma\Sigma$ and
  $C\Sigma\Sigma$ canonical-to-velocity corrections, while `QQQ` is unchanged;
- the two-point linear basis conversion round-trips all four blocks, and a
  non-equilateral model verifies that output-leg vertices are distinct;
- compile-fail tests reject both a mismatched $\alpha$ block and an incorrectly
  dimensioned $u^P{}_{QQ}$ block;
- the subhorizon initial condition reproduces exact rational in-in values for
  all eight blocks, is real, is covariant under both leg transpositions for
  arbitrary tensors, is invariant under leg exchange for an equilateral
  curved-space model, and solves the de Sitter transport equation to the
  retained order;
- the leg-resolved $B$ and $C$ tensors agree with direct evaluations in each
  Fourier-leg order on a curved two-field geometry, and doubling all
  wavenumbers at fixed shape leaves slow-$A$ unchanged while multiplying
  fast-$A$ by four;
- the kernel map reads separately evaluated input-exchanged $B$ and $C$
  tensors, and a source regression distinguishes $\Sigma^{ad}$ from
  $\Sigma^{da}$ on a left leg.

## Remaining P3d acceptance work

1. **Completed:** select the separately typed canonical-momentum phase-space
   route while preserving the velocity types used by P3c;
2. **Completed:** map dimensioned, raised canonical $A/B/C$ coefficients to
   the eight ordered e-fold-time $u^a{}_{bc}$ blocks;
3. **Completed:** evaluate $A$, $B$, and $C$ from the model, background,
   geometry, and triangle, and assemble the canonical e-fold kernel;
4. **Completed:** evaluate the three output-leg orderings and implement the
   typed linear/two-point and nonlinear/tree-level $V\leftrightarrow\Pi$
   conversion boundary;
5. **Completed:** implement the quantum initial condition for all eight
   blocks;
6. combine homogeneous and sourced terms into one adaptive transport system;
7. test the remaining general field-coordinate covariance limit; flat-space
   and single-field limits are covered for all coefficients, and the $A/B$
   Riemann plus $A$ rank-five $\nabla R$ terms already have independent
   componentwise checks;
8. **Completed for the initial condition:** verify general-triangle
   permutation and reality identities; they must be rechecked on evolved
   sourced states;
9. demonstrate convergence in integration tolerance and start depth.

Only after these items pass does P3e project $\alpha^{abc}$ to the observable
curvature bispectrum $B_\zeta$.

## Exact restart point

The next implementation session begins at item 6. The initial condition,
cubic coefficients, operator-ordered source, and $V\leftrightarrow\Pi$
boundary should not be modified unless an independent derivation exposes a
discrepancy.

The next slice should combine the homogeneous flow and the three ordered
$u\Sigma\Sigma$ terms into one adaptive canonical-basis system. The linear
covariances in the source must be evolved alongside $\alpha$ on the same
background, starting from the Bunch--Davies mode state on the common slice
used by `canonical_three_point_initial_condition`. The homogeneous system is
currently specialized to the velocity marker and $N=2$; the sourced system
must use `CanonicalMomentumPerturbation` throughout and convert to velocity
correlations only at output.

The first acceptance tests should be: agreement with the de Sitter consistency
of the initial condition over the first e-folds of evolution, stability of the
final state when the start depth increases, and convergence under tighter
tolerances. Only then should a physical benchmark be attempted.
