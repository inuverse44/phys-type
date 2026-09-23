# Three-point transport

## Status

P3d is in progress. Its state/permutation contract, homogeneous linear
transport, and typed quadratic-source contraction are implemented. The
model-derived cubic kernel and quantum initial conditions are not yet
implemented.

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

## Homogeneous transport

The state implements addition and multiplication by dimensionless Runge--Kutta
coefficients, so it can be used directly by the existing integrators.
`ThreePointHomogeneousTransportSystem` applies the P3c linear phase-space flow
to each Fourier leg,

```text
D_N Alpha^{abc}
  = u^a_d(k_1) Alpha^{dbc}
  + u^b_d(k_2) Alpha^{adc}
  + u^c_d(k_3) Alpha^{abd}.
```

It reads every Runge--Kutta stage from the shared P2 background provider. The
state legs and per-leg wavenumbers use `MomentumTriangle` canonical order. The
linear flow retains the declared P3c convention, including connection terms,
Hubble friction, the covariant perturbation mass matrix, and each leg's
`k_i^2/a^2` contribution.

The public name says `Homogeneous` deliberately: this system contains no cubic
kernel and makes no physical bispectrum claim. In particular, a zero state
remains zero.

## Quadratic source contraction

`QuadraticPhaseSpaceKernel<N>` stores the eight ordered blocks of `u^a_bc`.
Their dimensions are derived from

```text
dim(u^a_bc) = dim(X^a) - dim(X^b) - dim(X^c).
```

For a field output these are `M^(1/2)`, `M^(-1/2)`, `M^(-1/2)`, and
`M^(-3/2)` for `QQ`, `QP`, `PQ`, and `PP` inputs. For a momentum output they
are `M^(3/2)`, `M^(1/2)`, `M^(1/2)`, and `M^(-1/2)`.

`three_point_source_from_kernels` evaluates

```text
u^a_de Sigma^{db}(k_2) Sigma^{ec}(k_3) + two leg permutations
```

from three ordered kernels and three complex phase-space covariances. Kernel
`a` has output leg `a`; its input legs use ascending order among the remaining
legs: `[1,2]`, `[0,2]`, and `[0,1]`. This explicit convention prevents a
cyclic/ascending-order ambiguity for the middle output leg.

This layer performs only the hierarchy contraction. The next slice must derive
the actual curved-field-space kernel from the cubic Hamiltonian and the fixed
`P = D_t Q` convention. Quantum initial conditions and subhorizon start-depth
convergence remain acceptance requirements; zeroing `Alpha` at horizon exit is
not an allowed substitute.

## Verification in this slice

- all eight blocks have their declared compile-time dimensions;
- zero, addition, and dimensionless scaling preserve each block type;
- both generating leg exchanges move block and component axes together and
  square to the identity;
- complex conjugation squares to the identity;
- a compile-fail test rejects assigning the `M^(-1/2)` `QQP` block the
  `M^(-3/2)` dimension expected for `QQQ`.
- the `QQQ` derivative receives the expected `PQQ + QPQ + QQP` linear terms;
- for an equilateral triangle, homogeneous evolution commutes with a
  Fourier-leg exchange in all eight blocks and every field component.
- a scalar complex regression sums the three ordered `u Sigma Sigma` terms
  independently and leaves unsourced phase-space blocks zero;
- a compile-fail test rejects giving `u^P_QQ` the `M^(1/2)` dimension required
  by `u^Q_QQ` instead of its own `M^(3/2)` dimension.
