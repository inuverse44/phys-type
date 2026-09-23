# Three-point transport

## Status

P3d is in progress. Its first implementation slice fixes the state and
permutation contract before adding cubic dynamics.

`ThreePointPhaseSpaceCorrelation<N>` stores the connected equal-time
correlator

```text
Alpha^{abc}(k_1,k_2,k_3)
  = <X^a(k_1) X^b(k_2) X^c(k_3)>_connected,
X^a = (Q^I, P^I),
```

as eight separate complex blocks:

| Block | Dimension |
| --- | --- |
| `QQQ` | `M^(-3/2)` |
| `QQP`, `QPQ`, `PQQ` | `M^(-1/2)` |
| `QPP`, `PQP`, `PPQ` | `M^(1/2)` |
| `PPP` | `M^(3/2)` |

These dimensions follow the P3c convention `Q^I = U^I` and
`P^I = D_t U^I`. A block's three slot markers and its dimension are both part
of its Rust type. The public API therefore rejects, for example, treating a
`QQP` block as `QQQ`, even when both contain the same numerical array shape.

## Fourier-leg operations

The two generating transpositions exchange legs `(1 2)` and `(2 3)`. They
move together:

1. the corresponding array axis;
2. the phase-space slot marker;
3. the containing block within the eight-block state.

Applying either exchange twice is the identity. Their compositions generate
all six permutations of a momentum triangle. Complex conjugation does not
move a leg and is represented separately. This is the foundation for later
permutation and reality residuals; it does not assert that an arbitrary stored
state already satisfies those identities.

## Integration boundary

The state implements addition and multiplication by dimensionless Runge--Kutta
coefficients, so it can be used directly by the existing integrators. This
first slice deliberately contains no cubic kernel and makes no physical
bispectrum claim.

The next P3d slice will add the homogeneous three-leg action of `u^a_b` using
the same P2 background trajectory and per-leg wavenumbers as P3c. After that
is verified, the sourced term

```text
u^a_de Sigma^{db}(k_2) Sigma^{ec}(k_3) + two leg permutations
```

will be connected to an explicitly typed `u^a_bc`. Quantum initial conditions
and subhorizon start-depth convergence remain acceptance requirements; zeroing
`Alpha` at horizon exit is not an allowed substitute.

## Verification in this slice

- all eight blocks have their declared compile-time dimensions;
- zero, addition, and dimensionless scaling preserve each block type;
- both generating leg exchanges move block and component axes together and
  square to the identity;
- complex conjugation squares to the identity;
- a compile-fail test rejects assigning the `M^(-1/2)` `QQP` block the
  `M^(-3/2)` dimension expected for `QQQ`.
