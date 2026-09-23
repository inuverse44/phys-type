# Two-point transport oracle

## Purpose

P3c must demonstrate that direct phase-space covariance transport reproduces
the established linear mode solver before any nonlinear source is accepted.
This document fixes the phase-space convention and comparison boundary.  The
mode solver remains the independent numerical oracle; transport code must not
reuse its mode amplitudes after initialization.

## Phase-space convention

For each Fourier wavenumber, the transport coordinate is

```text
X^a = (Q^I, P^I),
Q^I = U^I,
P^I = D_t U^I.
```

The existing solver integrates the numerically better conditioned variables
`v^I_A = a U^I_A` and `w^I_A = D_t v^I_A`.  At any oracle snapshot,

```text
Q^I_A = v^I_A / a,
P^I_A = (w^I_A - H v^I_A) / a.
```

`P` is therefore a cosmic-time covariant velocity, not a canonical momentum
with an implicit factor of `a^3`, and not a conformal-time derivative.  A later
change of convention requires a new semantic type rather than reinterpretation
of `MomentumPerturbation`.

## Equal-time blocks

The unsymmetrized quantum covariance sums over the independent vacuum
oscillators `A`:

```text
Sigma_QQ^{IJ} = sum_A Q^I_A Q^{J*}_A,
Sigma_QP^{IJ} = sum_A Q^I_A P^{J*}_A,
Sigma_PQ^{IJ} = sum_A P^I_A Q^{J*}_A,
Sigma_PP^{IJ} = sum_A P^I_A P^{J*}_A.
```

Because `Q` has dimension `M^-1/2` and `P` has dimension `M^1/2`, the block
dimensions are:

| Block | Slot type | Dimension |
| --- | --- | --- |
| `QQ` | field, field | `M^-1` |
| `QP` | field, momentum | `M^0` |
| `PQ` | momentum, field | `M^0` |
| `PP` | momentum, momentum | `M^1` |

The oracle must satisfy `QQ = QQ^dagger`, `PP = PP^dagger`, and
`PQ = QP^dagger`.  Complex values are retained; taking real parts would erase
the quantum commutator information carried by the mixed blocks.

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
For `A^I_J = C^I_J + (k/a)^2 delta^I_J`, the mode-level equations are

```text
D_N Q^I = P^I / H,
D_N P^I = -3 P^I - A^I_J Q^J / H.
```

The implementation expands `D_N` into coordinate derivatives and applies the
connection flow `Gamma^I_JK pi^J/H` to both field-space indices of every
covariance block. It receives the background only through
`TwoFieldBackgroundProvider`; the P2 cubic-Hermite trajectory implements this
interface and is queried at every Runge--Kutta stage.

`TwoPointTransportTolerance` carries a correctly dimensioned absolute floor
for each block. The adaptive norm uses the largest complex modulus across a
whole block as its relative scale, so an exactly zero off-diagonal component
does not force an otherwise physical state to take vanishing steps.
