# Two-point transport oracle

## Purpose

P3c must demonstrate that direct phase-space covariance transport reproduces
the established linear mode solver before any nonlinear source is accepted.
This document fixes the phase-space convention and comparison boundary.  The
mode solver remains the independent numerical oracle; transport code must not
reuse its mode amplitudes after initialization.

## Phase-space convention

For each Fourier wavevector $\mathbf{k}$, with magnitude
$k=|\mathbf{k}|$, the transport coordinate is

$$
X^a=(Q^I,P^I),\qquad Q^I=U^I,\qquad P^I=D_tU^I.
$$

Field-space indices $I,J$ transform with $G_{IJ}$. A phase-space index $a$
means both a $Q/P$ slot and a field-space index. $D_t$ is the field-space
covariant cosmic-time derivative.

The existing solver integrates the numerically better conditioned variables
$v^I_A=aU^I_A$ and $w^I_A=D_tv^I_A$, where $A$ labels the independent
vacuum oscillator. At any oracle snapshot,

$$
Q^I_A=\frac{v^I_A}{a},\qquad
P^I_A=\frac{w^I_A-Hv^I_A}{a}.
$$

$P$ is therefore a cosmic-time covariant velocity, not a canonical momentum
with an implicit factor of $a^3$, and not a conformal-time derivative. A later
change of convention requires a new semantic type rather than reinterpretation
of `MomentumPerturbation`.

## Equal-time blocks

The unsymmetrized quantum covariance sums over the independent vacuum
oscillators $A$:

$$
\begin{aligned}
\Sigma_{QQ}^{IJ}(k)&=\sum_A Q^I_A(k)Q^{J*}_A(k),\\
\Sigma_{QP}^{IJ}(k)&=\sum_A Q^I_A(k)P^{J*}_A(k),\\
\Sigma_{PQ}^{IJ}(k)&=\sum_A P^I_A(k)Q^{J*}_A(k),\\
\Sigma_{PP}^{IJ}(k)&=\sum_A P^I_A(k)P^{J*}_A(k).
\end{aligned}
$$

The star denotes complex conjugation. These are unsymmetrized quantum
covariances; retaining their imaginary parts preserves the equal-time
commutator information needed by the three-point source.

Because $[Q]=M^{-1/2}$ and $[P]=M^{1/2}$, the block
dimensions are:

| Block | Slot type | Dimension |
| --- | --- | --- |
| `QQ` | field, field | `M^-1` |
| `QP` | field, momentum | `M^0` |
| `PQ` | momentum, field | `M^0` |
| `PP` | momentum, momentum | `M^1` |

The oracle must satisfy

$$
\Sigma_{QQ}=\Sigma_{QQ}^\dagger,\qquad
\Sigma_{PP}=\Sigma_{PP}^\dagger,\qquad
\Sigma_{PQ}=\Sigma_{QP}^\dagger.
$$

The dagger means complex conjugation together with exchange of the two field
indices.

## Ordered implementation and acceptance

1. Convert one complete mode-solver state into all four typed blocks and test
   the Bunch--Davies closed form. **Complete.**
2. Implement direct covariance evolution on the shared interpolated P2
   background, without evolving or consulting a mode matrix. **Complete.**
3. Compare every complex component of `QQ`, `QP`, `PQ`, and `PP` against mode
   snapshots in the subhorizon regime, at horizon exit, through the sharpest
   turn, and on the final slice. **Complete for the reference benchmark:** a
   shared fixed-step RK4 comparison bounds the maximum block-normalized error
   by `2e-7` at all four epochs.
4. Require convergence under tighter integration tolerance and deeper
   Bunch--Davies initialization before marking P3c complete. **Complete:**
   tightening the relative tolerance from `1e-7` to `1e-10` changes the final
   all-block result by less than `2e-6`, and the change from successively
   deeper starts at 3, 4, and 5 e-folds before exit shrinks by more than a
   factor of two.

All P3c gates pass. P3d may now build the sourced three-point hierarchy on
this two-point oracle.

## Direct evolution equation

`TwoPointTransportSystem` evolves the covariance itself in elapsed e-folds.
Let $C^I{}_J$ be the covariant perturbation mass matrix and define

$$
\mathcal{A}^I{}_J(k)
=C^I{}_J+\frac{k^2}{a^2}\delta^I{}_J.
$$

With $D_N\equiv H^{-1}D_t$, the mode-level equations are

$$
D_NQ^I=\frac{P^I}{H},\qquad
D_NP^I=-3P^I-\frac{\mathcal{A}^I{}_JQ^J}{H}.
$$

The implementation expands $D_N$ into coordinate derivatives and applies the
connection flow $\Gamma^I{}_{JK}\pi^J/H$ to both field-space indices of every
covariance block. It receives the background only through
`TwoFieldBackgroundProvider`; the P2 cubic-Hermite trajectory implements this
interface and is queried at every Runge--Kutta stage.

`TwoPointTransportTolerance` carries a correctly dimensioned absolute floor
for each block. The adaptive norm uses the largest complex modulus across a
whole block as its relative scale, so an exactly zero off-diagonal component
does not force an otherwise physical state to take vanishing steps.
