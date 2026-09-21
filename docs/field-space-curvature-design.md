# Field-space curvature and higher-rank tensor design

## Audience and purpose

This document defines the public geometry vocabulary used by `phys-type`.
It is written for users who know multifield dynamics but should not need to
know the crate's internal representation. The immediate target is the
field-space Riemann term in linear perturbations. The design also preserves a
clear route to cubic interactions and primordial non-Gaussianity.

The governing rule is that valid tensor expressions should be easy to write,
while an invalid index contraction or physical dimension should be rejected
by Rust before a numerical run begins.

## Rank is not derivative order

Two independent ideas must not be conflated:

- **Tensor rank** is the number of field-space indices. Riemann has rank four.
- **Automatic-differentiation order** is the highest derivative carried by a
  numerical value. `Jet<N>` currently carries derivatives through second
  order.

The linear perturbation mass matrix needs a rank-four Riemann tensor but only
a second derivative of the scalar potential. A future cubic calculation can
need both third derivatives of the potential and a rank-five covariant
derivative of Riemann.

## Conventions

Scalar coordinates have physical dimension `M^1`. Therefore a field-space
coordinate derivative has dimension `M^-1`. The initial convention table is:

| Object | Index form | Physical dimension | Tensor? |
| --- | --- | --- | --- |
| field velocity | `dot(phi)^I` | `M^2` | yes |
| metric | `G_IJ` | `M^0` | yes |
| inverse metric | `G^IJ` | `M^0` | yes |
| metric derivative | `partial_K G_IJ` | `M^-1` | no |
| connection | `Gamma^I_JK` | `M^-1` | no |
| Riemann curvature | `R^I_JKL` | `M^-2` | yes |
| covariant Riemann derivative | `nabla_M R^I_JKL` | `M^-3` | yes |
| potential Hessian | `V_;IJ` | `M^2` | yes |
| third potential derivative | `V_;IJK` | `M^1` | yes |

The Riemann sign and index order are fixed by

```math
[\nabla_K,\nabla_L]A^I = R^I{}_{JKL}A^J
```

and

```math
R^I{}_{JKL}
=\partial_K\Gamma^I{}_{LJ}-\partial_L\Gamma^I{}_{KJ}
+\Gamma^I{}_{KM}\Gamma^M{}_{LJ}
-\Gamma^I{}_{LM}\Gamma^M{}_{KJ}.
```

This convention is used directly in the perturbation mass matrix. A geometry
implementation using the opposite convention must convert at its boundary.

For the all-lower form, the Levi-Civita symmetries are

```math
R_{IJKL}=-R_{JIKL}=-R_{IJLK}=R_{KLIJ},
\qquad
R_{I[JKL]}=0.
```

## Public type model

Rank-one and rank-two values already use

```rust,ignore
FieldVector<I, D, V, N>
FieldTensor2<I, J, D, V, N>
```

The next implementation adds explicit fixed-rank types:

```rust,ignore
FieldTensor3<I, J, K, D, V, N>
FieldTensor4<I, J, K, L, D, V, N>
```

`I` through `L` are `Upper` or `Lower`, `D` is physical dimension, `V` is
the stored number, and `N` is the field count. Explicit rank types are chosen
over a recursive type-level index list because they produce substantially
shorter compiler diagnostics on stable Rust.

Rank-three and rank-four components are initially stored densely. Multifield
inflation normally has small `N`, and dense storage makes index order,
validation, and interoperability auditable. Symmetry compression can be added
later behind the same public semantic wrappers if profiling justifies it.

Semantic aliases and wrappers include

```rust,ignore
type RiemannTensor<N> =
    FieldTensor4<Upper, Lower, Lower, Lower, MMinus2, f64, N>;

type LoweredRiemannTensor<N> =
    FieldTensor4<Lower, Lower, Lower, Lower, MMinus2, f64, N>;
```

`ChristoffelSymbols` remains a dedicated non-tensor type. A future
`MetricDerivative` will also be a dedicated non-tensor type. Neither may be
passed to an API that requires a true `FieldTensor3`.

## Operations and index safety

Operations are named by the slot they affect. There is no ambiguous generic
`raise()` on a rank-four tensor:

```rust,ignore
let lowered: LoweredRiemannTensor<N> = metric.lower_first(&riemann);
let restored: RiemannTensor<N> = inverse_metric.raise_first(&lowered);
```

Likewise, contractions name their intended slots or expose a
domain-specific operation. The linear perturbation code should call a method
equivalent to

```rust,ignore
riemann.contract_velocities(&velocity)
```

to form `R^I_LMJ dot(phi)^L dot(phi)^M`, returning a rank-two mixed tensor of
dimension `M^2`. It should not reproduce four nested loops at each call site.

The following must fail to compile:

- contracting two upper indices;
- lowering an already lower slot;
- assigning `R_IJKL` to `R^I_JKL` without a metric;
- using a dimensionless rank-four tensor as Riemann;
- passing `ChristoffelSymbols` where a rank-three tensor is required.

The first higher-rank implementation only exposes operations needed by the
curvature and perturbation layers. A small correct API is preferable to a
large set of unverified permutations.

## Curvature provider boundary

`FieldSpaceGeometry<N>` currently supplies metric, inverse metric, and
connection. Curvature is added as a separate capability:

```rust,ignore
trait CurvedFieldSpaceGeometry<const N: usize>: FieldSpaceGeometry<N> {
    fn riemann(&self, position: &[Field; N]) -> RiemannTensor<N>;
}
```

This separation preserves existing background-only geometries. A flat
geometry implements the new trait by returning an exactly zero tensor.

The first version requires the geometry author to supply Riemann analytically.
Numerically differentiating Christoffel symbols is not an acceptable hidden
default because it can amplify noise and obscure the sign convention. A later
metric-autodiff geometry may derive inverse metric, connection, and curvature
from one metric formula, but it must be an explicit provider with independent
analytic tests.

## Validation and errors

Every public geometry implementation is testable through reusable checks:

- `G_IK G^KJ` is the identity;
- the connection is symmetric in its two lower slots for a torsion-free
  Levi-Civita geometry;
- metric compatibility `nabla_K G_IJ = 0` holds;
- lowered Riemann satisfies both pair antisymmetries, pair exchange, and the
  first Bianchi identity;
- flat Cartesian geometry has exactly zero connection and curvature;
- a known constant-curvature two-dimensional metric reproduces its analytic
  scalar curvature;
- raising and lowering Riemann round-trips within tolerance;
- the velocity-contracted curvature has dimension `M^2` and the correct mixed
  index variance.

Checked validation returns a structured `GeometryValidationError` containing
the failed identity, component indices, actual residual, and tolerance. It
does not return a bare Boolean or panic. Production evolution need not repeat
expensive global validation at every right-hand-side call.

## Documentation and usability contract

For a geometry to be considered supported, it must provide:

1. Rustdoc stating its coordinate chart, coordinate dimensions, domain, and
   Riemann convention.
2. A minimal construction example.
3. At least one analytic identity test independent of the implementation.
4. A statement of singular coordinate locations.
5. A reproducible example showing background and perturbation use.

Public tensor methods document their input and output index forms in both
mathematical notation and Rust types. CSV metadata records the geometry name,
coordinate chart, Riemann convention identifier, and crate revision. This is
part of the scientific result, not optional logging.

## Path to non-Gaussianity

Gong and Tanaka's covariant expansion shows that field-space curvature enters
systematically beyond quadratic order. Elliston, Seery, and Tavakol exhibit
curvature contributions to the inflationary three-point function, while the
transport formulation of Dias, Frazer, and Seery provides a practical route
to evolving observables in curved multifield models.

No bispectrum algorithm is selected by the present branch. The tensor layer
must nevertheless avoid blocking either an in-in/cubic-action implementation
or a transport implementation. The anticipated additional objects are:

- `V_;IJK`, a symmetric rank-three covariant tensor of dimension `M^1`;
- `nabla_M R^I_JKL`, a rank-five tensor of dimension `M^-3`;
- a rank-three field covariance/bispectrum object with momentum-triangle
  labels separate from field-space indices;
- third-order automatic differentiation, separate from tensor rank.

`FieldTensor3` and `FieldTensor4` therefore share component-access and
slot-naming conventions that a later `FieldTensor5` can follow without
changing existing APIs. Momentum indices, Fourier labels, and field-space
indices will remain distinct types.

## Branch sequence

The revised implementation order is:

1. `docs/field-space-curvature-design` -- conventions, API contract, and
   non-Gaussianity extension boundary.
2. `feature/higher-rank-field-tensors` -- dense rank-three/rank-four types,
   slot-specific raise/lower operations, and compile-fail tests.
3. `feature/field-space-curvature` -- Riemann semantic types, geometry
   capability, curvature contraction, symmetry validation, and analytic
   geometries.
4. `feature/two-field-perturbations` -- the complete linear mode equation and
   vacuum mode matrix.
5. `feature/complex-mode-values` -- dimensioned `Complex64` mode values and
   conjugation-aware diagnostics.
6. `feature/power-spectra` and `feature/two-field-spectrum-plots` -- physical
   spectra and figures already specified in the two-field design.
7. `docs/non-gaussianity-architecture` -- choose and specify the cubic/in-in,
   transport, or deliberately hybrid observable pipeline before implementing
   it.

## References

- J.-O. Gong and T. Tanaka,
  [A covariant approach to general field space metric in multi-field inflation](https://arxiv.org/abs/1101.4809),
  JCAP 03 (2011) 015; Erratum JCAP 02 (2012) E01.
- J. Elliston, D. Seery, and R. Tavakol,
  [The inflationary bispectrum with curved field-space](https://arxiv.org/abs/1208.6011),
  JCAP 11 (2012) 060.
- M. Dias, J. Frazer, and D. Seery,
  [Computing observables in curved multifield models of inflation](https://arxiv.org/abs/1502.03125),
  JCAP 12 (2015) 030.
- S. Garcia-Saenz, L. Pinol, and S. Renaux-Petel,
  [Revisiting non-Gaussianity in multifield inflation with curved field space](https://arxiv.org/abs/1907.10403),
  JHEP 01 (2020) 073.
