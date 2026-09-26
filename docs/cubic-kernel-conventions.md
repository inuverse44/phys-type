# Cubic kernel conventions and the momentum boundary

## Purpose

This document fixes the analytic boundary that must be crossed before the
physical P3d kernel is implemented. It records the cubic action, its momentum
arguments, the eight quadratic phase-space blocks, and the distinction between
the velocity perturbation used by this crate and the rescaled canonical
momentum used in the standard transport derivation.

The distinction is invisible in linear theory but not in the three-point
system. Therefore the formulas below are a derivation contract, not a claim
that the physical kernel is already implemented.

The primary convention is Butchers and Seery,
[Numerical evaluation of inflationary 3-point functions on curved field
space](https://arxiv.org/abs/1803.10563), especially their Eqs. (2.9)--(2.15),
(2.24), and Appendix A.4. The covariant perturbation and curvature conventions
are cross-checked against Elliston, Seery, and Tavakol,
[The inflationary bispectrum with curved field-space](https://arxiv.org/abs/1208.6011).

## Fixed notation

Field-space indices are $I,J,K,\ldots$. Fourier legs are $1,2,3$, with

$$
\mathbf{k}_1+\mathbf{k}_2+\mathbf{k}_3=0,
\qquad k_r=|\mathbf{k}_r|>0.
$$

The association between a field index and a momentum label is

$$
I\leftrightarrow\mathbf{k}_1,
\qquad J\leftrightarrow\mathbf{k}_2,
\qquad K\leftrightarrow\mathbf{k}_3.
$$

Consequently, symmetrizing field indices also permutes their attached Fourier
momenta. The needed scalar products are determined by the three magnitudes:

$$
\mathbf{k}_1\!\cdot\!\mathbf{k}_2
=\frac{k_3^2-k_1^2-k_2^2}{2},
\quad
\mathbf{k}_1\!\cdot\!\mathbf{k}_3
=\frac{k_2^2-k_1^2-k_3^2}{2},
\quad
\mathbf{k}_2\!\cdot\!\mathbf{k}_3
=\frac{k_1^2-k_2^2-k_3^2}{2}.
$$

`MomentumTriangle::wavevector_dot` implements these relations in canonical
leg order and returns a quantity of dimension $M^2$.
`wavevector_cosine` and `transverse_factor` provide

$$
\cos\theta_{ij}
=\frac{\mathbf{k}_i\cdot\mathbf{k}_j}{k_i k_j},
\qquad
T_{ij}=1-\cos^2\theta_{ij}.
$$

The cosine is clamped only for closure-tolerance roundoff at folded
boundaries. Tests fix the general $(4,3,2)$ dot products, the equilateral
$T_{ij}=3/4$, and the folded $T_{23}=0$ limit.

The crate evolves in e-fold time, $D_N=H^{-1}D_t$. The formulas quoted from
the paper are initially cosmic-time formulas. Every block of a cosmic-time
quadratic flow must therefore be divided by $H$ before it can populate an
e-fold-time kernel.

The repository uses reduced Planck units through
`CovariantMultifieldSystem::in_planck_units`. The $M_{\mathrm{Pl}}$ factors are
kept below so that dimensions and comparison with the source remain explicit;
only the implementation boundary may set $M_{\mathrm{Pl}}=1$.

## Cubic action

After the momentum-conserving delta function is suppressed, the covariant
cubic action can be written

$$
\begin{aligned}
S_\phi=\frac12\int dt\,a^3\bigl\{&
G_{IJ}D_tQ^I D_tQ^J+M_{IJ}Q^I Q^J\\
&+A_{IJK}Q^I Q^J Q^K
+B_{IJK}Q^I Q^J D_tQ^K
+C_{IJK}D_tQ^I D_tQ^J Q^K
\bigr\}.
\end{aligned}
$$

Define

$$
Z_I=D_t\dot\phi_I
+\frac{\dot\phi_I\dot\phi_J\dot\phi^J}
       {2H M_{\mathrm{Pl}}^2}.
$$

Before the required index/momentum symmetrization, the three species kernels
are

$$
\begin{aligned}
\widetilde A_{IJK}={}&
-\frac13V_{;IJK}
-\frac{\dot\phi_I V_{;JK}}{2H M_{\mathrm{Pl}}^2}
+\frac{\dot\phi_I\dot\phi_J Z_K}{4H^2M_{\mathrm{Pl}}^4}\\
&+\frac{\dot\phi_I Z_JZ_K}{8H^3M_{\mathrm{Pl}}^4}
 \left[1-\frac{(\mathbf{k}_2\cdot\mathbf{k}_3)^2}{k_2^2k_3^2}\right]\\
&+\frac{\dot\phi_I\dot\phi_J\dot\phi_K}
        {8H^3M_{\mathrm{Pl}}^6}
 \left(6H^2M_{\mathrm{Pl}}^2-\dot\phi^2\right)\\
&-\frac{\dot\phi_K\dot\phi^L\dot\phi^M}
        {2HM_{\mathrm{Pl}}^2}R_{L(IJ)M}
+\frac13R_{(I|LM|J;K)}\dot\phi^L\dot\phi^M\\
&+\frac{\dot\phi_I G_{JK}}{2HM_{\mathrm{Pl}}^2}
  \frac{\mathbf{k}_2\cdot\mathbf{k}_3}{a^2},
\end{aligned}
$$

$$
\begin{aligned}
\widetilde B_{IJK}={}&
\frac43R_{K(IJ)L}\dot\phi^L
-\frac{\dot\phi_I Z_J\dot\phi_K}{4H^3M_{\mathrm{Pl}}^4}
 \left[1-\frac{(\mathbf{k}_2\cdot\mathbf{k}_3)^2}{k_2^2k_3^2}\right]\\
&+\frac{\dot\phi_I\dot\phi_J\dot\phi_K}
        {4H^2M_{\mathrm{Pl}}^4}
-\frac{Z_I G_{JK}}{HM_{\mathrm{Pl}}^2}
  \frac{\mathbf{k}_1\cdot\mathbf{k}_2}{k_1^2},
\end{aligned}
$$

$$
\begin{aligned}
\widetilde C_{IJK}={}&
-\frac{G_{IJ}\dot\phi_K}{2HM_{\mathrm{Pl}}^2}
+\frac{\dot\phi_I\dot\phi_J\dot\phi_K}{8H^3M_{\mathrm{Pl}}^4}
 \left[1-\frac{(\mathbf{k}_1\cdot\mathbf{k}_2)^2}{k_1^2k_2^2}\right]\\
&+\frac{\dot\phi_I G_{JK}}{HM_{\mathrm{Pl}}^2}
  \frac{\mathbf{k}_1\cdot\mathbf{k}_3}{k_1^2}.
\end{aligned}
$$

Here parentheses denote symmetrization with unit weight and indices between
vertical bars are excluded. The kernels that enter the Hamiltonian are

$$
A_{IJK}=\operatorname{Sym}_{IJK}\widetilde A_{IJK},
\qquad
B_{IJK}=\operatorname{Sym}_{IJ}\widetilde B_{IJK},
\qquad
C_{IJK}=\operatorname{Sym}_{IJ}\widetilde C_{IJK},
$$

where every permutation moves the corresponding $\mathbf{k}_r$ at the same
time. In particular, `symmetrize tensor components and leave the momentum
arguments fixed` is not the required operation.

## The two meanings of momentum

The phase-space convention already fixed and tested in P3c is

$$
X_V^a=(Q^I,V^I),
\qquad V^I\equiv D_tQ^I.
$$

The public Rust marker is currently named `MomentumPerturbation`, but its
documented physical meaning is the covariant velocity $V^I$. It is not the
canonical momentum.

The cubic Hamiltonian derivation instead introduces the conjugate momentum
$p_I=\delta S/\delta(D_tQ^I)$ and then removes its factor of $a^3$. Denote the
rescaled variable by

$$
\Pi_I\equiv a^{-3}p_I.
$$

Suppressing convolution integrals and momentum-reversal labels, its nonlinear
relation to the velocity is, through quadratic order,

$$
\Pi_I
=D_tQ_I
+\frac12B_{JKI}Q^JQ^K
+C_{IJK}\Pi^JQ^K
+O(Q^3).
$$

Equivalently,

$$
V^I
=\Pi^I
-\frac12B_{JK}{}^I Q^JQ^K
-C^I{}_{JK}\Pi^JQ^K
+O(Q^3).
$$

Therefore $V^I=\Pi^I$ at linear order, which is why P3c two-point validation
does not expose the distinction. At quadratic order they differ, and that
difference contributes at exactly the order needed for the three-point
function.

## Canonical quadratic flow and the eight-block table

For the canonical phase-space variable

$$
X_\Pi^a=(Q^I,\Pi^I),
$$

the Hamilton equations are

$$
\begin{aligned}
D_tQ^I={}&\Pi^I
-\frac12B_{JK}{}^I Q^JQ^K
-C^I{}_{JK}\Pi^JQ^K,\\
D_t\Pi^I={}&-3H\Pi^I+M^I{}_JQ^J
+\frac32A^I{}_{JK}Q^JQ^K
+B^I{}_{JK}Q^J\Pi^K
+\frac12C_{JK}{}^I\Pi^J\Pi^K.
\end{aligned}
$$

The vertical position of $I$ is not cosmetic. The required raised tensors are

$$
\begin{aligned}
A^I{}_{JK}&=G^{IL}A_{LJK},\\
B_{JK}{}^I&=G^{IL}B_{JKL},
&B^I{}_{JK}&=G^{IL}B_{LJK},\\
C_{JK}{}^I&=G^{IL}C_{JKL},
&C^I{}_{JK}&=G^{IL}C_{LJK}.
\end{aligned}
$$

Because $B$ and $C$ are not fully symmetric, raising their first or third
field index gives different contractions.

Comparing these equations with

$$
D_tX^a=u^a{}_bX^b+\frac12u^a{}_{bc}X^bX^c
$$

gives the following ordered cosmic-time blocks. The field indices $I,J,K$
belong to the output, first input, and second input respectively.

| Rust-style block | Canonical coefficient in cosmic time |
| --- | --- |
| `q_qq` | $u^{Q^I}{}_{Q^JQ^K}=-B_{JK}{}^I$ |
| `q_qp` | $u^{Q^I}{}_{Q^J\Pi^K}=-C^I{}_{KJ}$ |
| `q_pq` | $u^{Q^I}{}_{\Pi^JQ^K}=-C^I{}_{JK}$ |
| `q_pp` | $u^{Q^I}{}_{\Pi^J\Pi^K}=0$ |
| `p_qq` | $u^{\Pi^I}{}_{Q^JQ^K}=3A^I{}_{JK}$ |
| `p_qp` | $u^{\Pi^I}{}_{Q^J\Pi^K}=B^I{}_{JK}$ |
| `p_pq` | $u^{\Pi^I}{}_{\Pi^JQ^K}=B^I{}_{KJ}$ |
| `p_pp` | $u^{\Pi^I}{}_{\Pi^J\Pi^K}=C_{KJ}{}^I$ |

This table is the unambiguous reading of the two matrix slices printed in the
reference. The slices are selected by the third phase-space index; reading
them as two output-index slices exchanges several blocks.

For the crate's e-fold equation, every entry in the right column becomes
$H^{-1}$ times the displayed coefficient.

## Why these coefficients cannot populate the current type directly

`QuadraticPhaseSpaceKernel<N>` is presently parameterized by the P3c slot
marker whose meaning is $V^I=D_tQ^I$. The table above is derived for
$\Pi^I$. Both variables have the same physical dimension, so compile-time
dimension checks cannot distinguish them.

Copying the table directly would therefore create a semantic error that the
type system is not currently able to reject. In particular:

- `q_qq` and the mixed `q_qp`/`q_pq` blocks encode the nonlinear difference
  between $V$ and $\Pi$;
- transforming the $\Pi$ flow to a $V$ flow differentiates the nonlinear
  field redefinition and introduces the corresponding time-dependent terms;
- three-point correlators in the two variables also differ by local products
  of two-point covariances.

The implementation must choose one of two explicit routes:

1. introduce a distinct canonical-momentum slot marker and evolve
   $(Q,\Pi)$, including typed conversion of initial data and output
   correlators; or
2. analytically transform the full quadratic flow and initial condition to
   $(Q,V)$ before constructing the existing kernel type.

The repository selects the first route. `CanonicalMomentumPerturbation` is a
distinct public marker, and the two-point covariance, three-point correlation,
and quadratic kernel types are generic over their second phase-space slot.
Their existing one-argument forms retain `MomentumPerturbation` as the default,
so P3c remains source-compatible and continues to mean $V^I=D_tQ^I$. The
canonical aliases are:

- `CanonicalTwoPointPhaseSpaceCovariance<N>`;
- `CanonicalThreePointPhaseSpaceCorrelation<N>`;
- `CanonicalQuadraticPhaseSpaceKernel<N>`.

The source contraction requires the same marker on the kernel and both
correlator types. A compile-fail test proves that a canonical kernel cannot be
contracted with velocity covariances. Erasing a `PhaseSpaceIndex` retains only
the structural field/momentum position; semantic safety is carried by the
typed containers and must be preserved until after contraction.

The next mapping layer is also implemented.
`CanonicalCubicHamiltonianCoefficients<N>` stores the raised tensors required
by the eight-block table as five distinct fields:

| Field | Components | Dimension |
| --- | --- | --- |
| `a_first_raised` | $A^I{}_{JK}$ | $M^1$ |
| `b_first_raised` | $B^I{}_{JK}$ | $M^0$ |
| `b_third_raised` | $B_{JK}{}^I$ | $M^0$ |
| `c_first_raised` | $C^I{}_{JK}$ | $M^{-1}$ |
| `c_third_raised` | $C_{JK}{}^I$ | $M^{-1}$ |

Every array is stored in output/first-input/second-input order $[I][J][K]$.
`CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian` applies the
required $J\leftrightarrow K$ exchanges and divides every cosmic-time
coefficient by a positive finite $H$ exactly once. Tests use deliberately
asymmetric component values so that either an index exchange or a first/third
raising mistake changes the result. A compile-fail test fixes all five mass
dimensions.

These integer dimensions are also obtained term by term from the displayed
cubic action. The earlier half-integer assignment incorrectly cubed the
dimension of a linear mode function after stripping only one three-point
delta function. The Fourier-operator and convolution-measure accounting is
given explicitly in [`three-point-transport.md`](three-point-transport.md).

This layer maps already evaluated $A$, $B$, and $C$ tensors. All three physical
coefficient evaluators are now implemented:

- `canonical_cubic_c_coefficients` computes the displayed
  $\widetilde C_{IJK}$ from a `CovariantMultifieldSystem`, `MultifieldState`,
  and `MomentumTriangle`;
- it performs $\operatorname{Sym}_{IJ}$ by exchanging the field indices and
  their attached Fourier legs together;
- it returns both $C^I{}_{JK}$ and $C_{JK}{}^I$ in the storage order required
  by the eight-block map;
- dimensionful arithmetic remains in `Quantity<MMinus1, f64>`, while only
  metric components and momentum ratios are raw dimensionless scalars;
- tests cover the analytic flat single-field equilateral result, the two
  raised-index layouts in flat two-field space, and invalid $H$ or
  $M_{\mathrm{Pl}}$ boundaries.

`canonical_cubic_b_coefficients` follows the same typed symmetrization and
raised-index storage rules for $B$. It constructs

$$
Z_I=-3H\dot\phi_I-V_{,I}
+\frac{\dot\phi_I\dot\phi^2}{2HM_{\mathrm{Pl}}^2}
$$

from the covariant background equation, lowers the first Riemann index before
forming $R_{K(IJ)L}\dot\phi^L$, and keeps the completed coefficient in
`Quantity<M0, f64>`. In addition to the flat single-field and index-layout
tests, a curved-space regression compares the result with an otherwise
identical geometry whose Riemann tensor is suppressed. Their difference is
checked component by component against the typed curvature contraction.

`canonical_cubic_a_coefficients` evaluates the fully symmetrized $A$ tensor in
`Quantity<M1, f64>`. Its trait boundary requires both
`ThirdOrderFieldSpaceGeometry` and `CurvatureDerivativeGeometry`; therefore a
geometry cannot enter this evaluator without explicitly supplying the
connection derivative, Riemann tensor, and $\nabla_MR^I{}_{JKL}$. Tests isolate
the following contributions:

- the flat single-field analytic limit;
- $V_{;(IJK)}$ using a mixed cubic potential at zero background velocity;
- the Riemann term by comparison with a curvature-suppressed geometry;
- the rank-five curvature-derivative index routing using a sparse probe
  tensor.

Because $A$ contains the physical-wavenumber factor $a^{-2}\mathbf{k}_2\cdot
\mathbf{k}_3$, its evaluator also rejects non-positive or non-finite scale
factors before division. All coefficient evaluators reject invalid $H$ and
$M_{\mathrm{Pl}}$ values.

`canonical_cubic_hamiltonian_coefficients` is the complete model-facing entry
point. It evaluates $A/B/C$ and produces the five raised tensors consumed by
`CanonicalQuadraticPhaseSpaceKernel::from_cubic_hamiltonian`.

The cyclic source needs three differently ordered vertices, rather than three
references to this one canonical ordering. The
`canonical_cubic_hamiltonian_vertices` entry point returns a
`CanonicalCubicHamiltonianVertices<N>` whose output-leg order is

| Output leg | Input legs |
| --- | --- |
| $1$ | $(2,3)$ |
| $2$ | $(1,3)$ |
| $3$ | $(1,2)$ |

The field indices and the Fourier legs are permuted together before $A$, $B$,
and $C$ are evaluated. `e_fold_kernels` then converts all three vertices to
the ordered kernels expected by `three_point_source_from_kernels`.
For the nonsymmetric $B$ and $C$ tensors this rule also applies while changing
which covariant index becomes the output: the third-raised forms evaluate
$B_{JKL}$ and $C_{JKL}$ with leg order `(input 1, input 2, output)` before
raising $L$. They are not obtained by reindexing values evaluated in
`(output, input 1, input 2)` order.

## Typed canonical--velocity conversion

The public conversion boundary now implements the perturbative relations
displayed above. At linear order,
`canonical_two_point_from_velocity_linear` and
`velocity_two_point_from_canonical_linear` preserve all four numerical
covariance blocks while changing the phase-space marker. Their names make the
linear-order restriction explicit.

At tree level, `velocity_three_point_from_canonical_tree` applies

$$
\begin{aligned}
N^{V^I}{}_{Q^JQ^K}&=-B_{JK}{}^I,\\
N^{V^I}{}_{Q^J\Pi^K}&=-C^I{}_{KJ},\\
N^{V^I}{}_{\Pi^JQ^K}&=-C^I{}_{JK},\\
N^{Q^I}{}_{bc}&=N^{V^I}{}_{\Pi^J\Pi^K}=0,
\end{aligned}
$$

to each output leg. For a quadratic field redefinition
$X_V^a=X_\Pi^a+\tfrac12N^a{}_{bc}X_\Pi^bX_\Pi^c$, the connected tree-level
correlator is therefore

$$
\alpha_V^{abc}=\alpha_\Pi^{abc}
+N^a{}_{de}\Sigma_\Pi^{db}\Sigma_\Pi^{ec}
+N^b{}_{de}\Sigma_\Pi^{da}\Sigma_\Pi^{ec}
+N^c{}_{de}\Sigma_\Pi^{da}\Sigma_\Pi^{eb}.
$$

The conversion accepts only canonical correlations and covariances and returns
the velocity-marked correlation. Tests check the exact scalar $B$ and $C$
contributions, preservation of `QQQ`, the linear two-point round trip, and the
fact that a non-equilateral model produces distinct output-leg vertices.

## Acceptance checklist for implementation

The physical kernel is not accepted until all of the following are explicit
in code and tests:

1. **Implemented:** the phase-space variable is named and typed as either $V$
   or $\Pi$;
2. **Implemented:** cosmic-time coefficients are converted once, and only
   once, to e-fold coefficients by $H^{-1}$;
3. **Implemented:** every field-index symmetrization moves the attached
   Fourier legs simultaneously;
4. **Implemented:** triangle dot products use the declared all-incoming
   convention;
5. **Implemented:** $M_{\mathrm{Pl}}$ is read from the system and
   Planck units enter only through the system constructor;
6. **Implemented:** flat-field-space and single-field reductions are checked
   analytically, and the $A/B$ curvature terms plus the $A$ rank-five term have
   componentwise regressions;
7. **Implemented:** exchanging the two input legs maps each ordered block to
   its partner and all coefficients are symmetrized with their momentum labels;
8. **Implemented:** the nonlinear $V\leftrightarrow\Pi$ correlator conversion
   is tested for the selected canonical route;
9. dimensions, permutation identities, and an independent numerical
   benchmark all pass.

The model-derived physical kernel and its variable boundary are now present.
The remaining acceptance gap is item 9's independent numerical benchmark,
along with the quantum initial condition and sourced integration described in
[`three-point-transport.md`](three-point-transport.md).
