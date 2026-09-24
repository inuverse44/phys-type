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
- contraction of three ordered $u\Sigma\Sigma$ source terms.

The following pieces are not implemented:

- the model-derived curved-field-space coefficients of $u^a{}_{bc}$;
- the in-in/Bunch--Davies initial value for the three-point function;
- one adaptive system that evolves the sourced state through a physical
  benchmark;
- general-triangle permutation, reality, and convergence tests;
- nonlinear projection to the curvature bispectrum.

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

In the Fourier normalization used by this crate,
$[Q]=M^{-1/2}$ and $[P]=M^{1/2}$.

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
| `QQQ` | $\langle Q^I Q^J Q^K\rangle_c$ | $M^{-3/2}$ |
| `QQP`, `QPQ`, `PQQ` | one momentum slot | $M^{-1/2}$ |
| `QPP`, `PQP`, `PPQ` | two momentum slots | $M^{1/2}$ |
| `PPP` | $\langle P^I P^J P^K\rangle_c$ | $M^{3/2}$ |

The slot markers and physical dimension are both part of each Rust type. For
example, a numerical array stored as `QQP` cannot be assigned to `QQQ`.

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
  \Sigma^{da}(k_1)\Sigma^{ec}(k_3)\\
&+u^c{}_{de}(\mathbf{k}_3;\mathbf{k}_1,\mathbf{k}_2)
  \Sigma^{da}(k_1)\Sigma^{eb}(k_2).
\end{aligned}
$$

The sourced transport equation will be

$$
D_N\alpha^{abc}
=\left.D_N\alpha^{abc}\right|_{\mathrm{hom}}+S^{abc}.
$$

`QuadraticPhaseSpaceKernel<N>` stores the eight ordered $Q/P$ blocks of
$u^a{}_{bc}$. Dimensional consistency requires

$$
[u^a{}_{bc}]=[X^a]-[X^b]-[X^c].
$$

Therefore:

| Output | Inputs `QQ` | Inputs `QP` | Inputs `PQ` | Inputs `PP` |
| --- | ---: | ---: | ---: | ---: |
| $Q$ | $M^{1/2}$ | $M^{-1/2}$ | $M^{-1/2}$ | $M^{-3/2}$ |
| $P$ | $M^{3/2}$ | $M^{1/2}$ | $M^{1/2}$ | $M^{-1/2}$ |

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

## Why the source algebra is not yet the physical kernel

The current contraction accepts correctly dimensioned numerical kernel
blocks, but does not manufacture their values. The next implementation must
derive those values from the curved-field-space cubic Hamiltonian and convert
them consistently to $P^I=D_tQ^I$, or introduce a distinct typed canonical
momentum $\Pi^I$ and keep the conversion explicit.

That conversion must include the terms involving the symmetrized covariant
potential derivative $V_{;(IJK)}$, field-space curvature, and the required
covariant curvature derivative. These are the reasons P3b implemented the
third-order geometry first.

This boundary is important: compile-time dimensions can reject a kernel in the
wrong physical slot, but they cannot prove that a coefficient, sign, momentum
argument, or curvature contraction was derived correctly. Those properties
require analytic limits, permutation tests, and comparison with an independent
calculation.

The full formulas for $A_{IJK}$, $B_{IJK}$, $C_{IJK}$, their simultaneous
index/momentum symmetrizations, the canonical eight-block $u^a{}_{bc}$ table,
and the unresolved $P$--$\Pi$ decision are recorded in
[`cubic-kernel-conventions.md`](cubic-kernel-conventions.md). That document is
the implementation gate for the next P3d slice.

## Initial conditions

All three legs begin at one common time when the triangle is sufficiently deep
inside the horizon. If

$$
x_r\equiv\frac{k_r}{aH},
$$

the start criterion must control all three $x_r$, with the smallest $k_r$
usually setting the most restrictive leg. The initial $\Sigma^{ab}$ values
come from the validated Bunch--Davies mode state.

The initial $\alpha^{abc}$ is not generally zero. It must be obtained from the
early-time in-in integral or an analytically equivalent quantum-transport
initial condition. Setting $\alpha=0$ at horizon exit would discard intrinsic
horizon-crossing interactions and distort nonlocal triangle shapes.

## Verification already present

- all eight $\alpha$ blocks have compile-time slot and dimension checks;
- addition and multiplication by dimensionless Runge--Kutta coefficients
  preserve every block type;
- $(12)^2=(23)^2=1$ and complex conjugation squares to the identity;
- the `QQQ` homogeneous derivative contains the expected
  $(\mathrm{PQQ}+\mathrm{QPQ}+\mathrm{QQP})/H$ terms;
- homogeneous evolution commutes with an equal-leg exchange for an
  equilateral triangle in every block and field component;
- a scalar complex regression independently reproduces all three ordered
  $u\Sigma\Sigma$ terms;
- compile-fail tests reject both a mismatched $\alpha$ block and an incorrectly
  dimensioned $u^P{}_{QQ}$ block.

## Remaining P3d acceptance work

1. select and type the velocity or canonical-momentum phase-space route;
2. derive and implement the corresponding curved-field-space
   $u^a{}_{bc}$ coefficients;
3. test flat-space, single-field, and field-coordinate covariance limits;
4. implement the quantum initial condition for all eight blocks;
5. combine homogeneous and sourced terms into one adaptive transport system;
6. verify general-triangle permutation and reality identities;
7. demonstrate convergence in integration tolerance and start depth.

Only after these items pass does P3e project $\alpha^{abc}$ to the observable
curvature bispectrum $B_\zeta$.
