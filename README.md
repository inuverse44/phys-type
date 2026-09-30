# phys-type

Dimension-aware physical quantities for numerical physics, in Rust.

Physical dimensions live in the type system, so a dimensionally invalid
expression fails to compile instead of producing a wrong number. The first
target is inflationary background evolution; see
[`docs/mvp-design.md`](docs/mvp-design.md). The next scientific target is a
covariant two-field adiabatic/isocurvature spectrum, specified in
[`docs/two-field-perturbation-design.md`](docs/two-field-perturbation-design.md).
The runnable double-quadratic reference model and its scalar-amplitude
normalization are specified in
[`docs/two-field-observational-benchmark.md`](docs/two-field-observational-benchmark.md).
Its Planck 2018 data-release, likelihood-component, and physical-pivot
conventions—and the boundary between a benchmark target and an actual fit—are
fixed in
[`docs/planck-2018-likelihood-boundary.md`](docs/planck-2018-likelihood-boundary.md).
Higher-rank tensors, field-space Riemann curvature, and the extension boundary
toward non-Gaussianity are specified in
[`docs/field-space-curvature-design.md`](docs/field-space-curvature-design.md).
The selected tree-level phase-space transport architecture and its ordered
two-/three-point validation milestones are specified in
[`docs/non-gaussianity-architecture.md`](docs/non-gaussianity-architecture.md).
The active P3d state and dynamics boundary is tracked in
[`docs/three-point-transport.md`](docs/three-point-transport.md).
The cubic-action formulas and the nonlinear distinction between
$P^I=D_tQ^I$ and rescaled canonical momentum are fixed in
[`docs/cubic-kernel-conventions.md`](docs/cubic-kernel-conventions.md).
The corresponding public phase-space containers carry distinct velocity and
canonical-momentum marker types, preventing mixed source contractions. Typed
linear two-point and nonlinear tree-level three-point conversions make the
boundary explicit when entering or leaving canonical transport.
The implemented milestone, known limitations, release state, and ordered next
steps are summarized in
[`docs/current-status-and-roadmap.md`](docs/current-status-and-roadmap.md).
The shared background trajectory and independent-mode optimization are
specified and numerically compared with the reference driver in
[`docs/background-checkpoint-interpolation.md`](docs/background-checkpoint-interpolation.md).

## Running it

This crate is a library, so `cargo run` on its own has nothing to run. Use:

```sh
cargo test                      # run every test, including the compile-fail ones
cargo check                     # type-check quickly
cargo fmt                       # format
cargo clippy --all-targets      # lint
```

To integrate a quadratic inflationary background and write CSV:

```sh
cargo run --release --example quadratic_inflation > inflation.csv
```

The CSV goes to standard output and a summary to standard error. For the full
set of diagnostic CSV files, including the inflationary attractor and ten
post-inflation oscillations, run:

```sh
cargo run --release --example diagnostics
```

For a complete two-field adiabatic and entropy spectrum run, writing the
background, horizon-exit spectra, requested-output spectra, and mode-evolution
CSV files with reproducibility metadata:

```sh
cargo run --release --example two_field_spectra
```

The matching four-panel figure, showing the trajectory and its frame, the
field spectra, the curvature/entropy spectra with the signed cross spectrum,
and selected modes through horizon exit:

```sh
uv run plots/two_field_spectra.py diagnostics
```

The focused final spectrum figure puts `P_R(k)` and `P_S(k)` in two aligned
log-log panels at the end of inflation, so a strongly decayed entropy spectrum
does not flatten the visible curvature-spectrum shape:

```sh
uv run plots/final_spectra.py diagnostics/two_field_spectra.csv
```

The time-domain figure corresponding to the familiar mass-ratio-seven
double-inflation plot follows the pivot mode through entropy decay and the
turn-induced growth of curvature power:

```sh
uv run plots/pivot_mode_evolution.py diagnostics
```

A reference result generated at the end of inflation is tracked as
[`docs/figures/two_field_final_spectra.png`](docs/figures/two_field_final_spectra.png)
with a vector
[`PDF`](docs/figures/two_field_final_spectra.pdf).
The corresponding pivot-mode history is tracked as
[`docs/figures/two_field_pivot_evolution.png`](docs/figures/two_field_pivot_evolution.png)
with a vector
[`PDF`](docs/figures/two_field_pivot_evolution.pdf).

Generate every diagnostic CSV and publication-quality PNG/PDF figure with one
command:

```sh
uv run plots/generate_all.py
```

Individual plots can also be regenerated from existing CSV files:

```sh
uv run plots/background.py diagnostics/trajectory.csv
uv run plots/attractor.py diagnostics/attractor.csv
uv run plots/convergence.py diagnostics/convergence.csv
uv run plots/oscillation.py diagnostics/oscillation.csv
uv run plots/two_field_spectra.py diagnostics
uv run plots/final_spectra.py diagnostics/two_field_spectra.csv
uv run plots/pivot_mode_evolution.py diagnostics
```

The plotting scripts declare their own Python dependencies, so `uv` creates an
isolated environment automatically. The visualization design and physical
checks are documented in
[`docs/visualization-design.md`](docs/visualization-design.md).

## What the library provides

| Module | Contents |
| --- | --- |
| `autodiff` | Forward autodiff through third order with `Dual`, `Jet<N>`, and `ThirdOrderJet<N>` |
| `dimension` | `MassPower<E>` and the type-level `DimMul`, `DimDiv`, `DimSqrt` |
| `field_space` | Typed rank-one through rank-four tensors, geometry, and bases |
| `quantity` | `Quantity<D, V>` and its arithmetic |
| `time` | Semantically distinct conformal time for coordinate-safe integration |
| `state` | `State`, `StateDerivative`, `StateDelta` |
| `potential` | Automatically differentiated potentials and `QuadraticPotential` |
| `perturbation` | Cosmic/conformal two-field mode systems and Bunch-Davies evolution |
| `spectra` | Mode covariance, Hermiticity checks, and adiabatic/entropy projections |
| `spectrum_run` | The background survey, pivot convention, `k` grid driver, and CSV rows |
| `transport` | Validated momentum triangles and typed phase-space correlator blocks |
| `integrate` | The `System` and `AdaptiveSystem` traits, fixed RK4, embedded Dormand-Prince 5(4), and the shared step-size controller |
| `scalar_field` | `ScalarFieldSystem`, including the Friedmann constraint |
| `multifield` | Const-generic canonical multifield states and dynamics |
| `background` | The integration driver, termination conditions, and CSV export |

Units are natural (`c = hbar = 1`) with `M_pl = 1`, so every dimension is a
power of mass. Exponents are stored in half-power steps, which keeps
`sqrt(M^1) = M^(1/2)` representable.

```rust
use phys_type::{AdaptiveRunSettings, Field, Mass, QuadraticPotential, ScalarFieldSystem, ScaleFactor};

let system = ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(1.0e-6)));
let initial = system.slow_roll_initial_state(Field::new(16.0), ScaleFactor::new(1.0));
let run = system.run_adaptive(initial, AdaptiveRunSettings::default());
```

Each adaptive attempt records its step size, normalized local error, and
acceptance decision in `run.adaptive_steps`. Fixed-step RK4 remains available
through `system.run(initial, RunSettings::default())`.

Step-size control lives in `integrate` rather than in any one system. A system
implements `AdaptiveSystem` to reduce an embedded error estimate to a single
dimensionless ratio; `adaptive_step` and `next_step_size` then handle
acceptance and step growth for every system alike, from the single-field
background to the complex two-field modes.

Canonical multifield systems use fixed-size arrays and the same generic
integrators. `QuadraticMultifieldPotential<N>` supplies one mass per field and
obtains every gradient component through automatic differentiation.

The multifield layer includes canonical and covariant backgrounds, complex
two-field perturbations, and projected curvature/entropy spectra. Spectrum
runs initialize each mode in the Bunch-Davies vacuum and evolve `v = a U` in
conformal time. `ConformalTime` is a distinct semantic type from cosmic
`Time`, while the conformal background evolves `ln a` to keep the expanding
scale factor out of the numerical dynamic range.

Field-space variance is already available independently of physical
dimension. `FieldMetric<N>::lower` maps an upper-index vector to a lower-index
vector, `InverseFieldMetric<N>::raise` performs the inverse operation, and
only opposite-variance vectors expose contraction. Invalid index operations
are covered by compile-fail tests.

`FieldTensor3`, `FieldTensor4`, and `FieldTensor5` extend the same typed-index
model to higher rank. Metrics implement slot-specific `LowerFirst` and inverse metrics
implement `RaiseFirst`; already-lower indices are rejected at compile time.
`ChristoffelSymbols` and `ChristoffelDerivatives` deliberately remain distinct
non-tensor types.

`RiemannTensor<N>` fixes the convention $R^I{}_{JKL}$ and dimension $M^{-2}$.
`CurvedFieldSpaceGeometry<N>` supplies it analytically, while
`contract_velocities` produces the mixed `M^2` curvature contribution needed
by the perturbation mass matrix. All-lower curvature can be checked for pair
symmetries and the first Bianchi identity with structured validation errors.
Third-order geometry supplies the symmetrized $V_{;(IJK)}$ and the typed
rank-five $\nabla_MR^I{}_{JKL}$ capability needed by cubic transport.

The transport layer fixes $X^a=(Q^I,P^I)$ with $P^I=D_tQ^I$ and can
convert a complete linear mode-solver state into dimensioned complex `QQ`,
`QP`, `PQ`, and `PP` oracle blocks. `TwoPointTransportSystem` evolves those
blocks directly on the shared interpolated background and is regression-tested
against the mode oracle at four physical epochs. Dimensioned adaptive
tolerances and initialization-depth convergence complete the P3c validation.
P3d now provides the complete dimensioned `QQQ`--`PPP` state, homogeneous
three-leg evolution, the typed quadratic kernel $u^a{}_{bc}$, and the ordered
$u\Sigma\Sigma$ contraction. Model-derived $A_{IJK}$, $B_{IJK}$, and $C_{IJK}$
evaluation is also implemented with typed dimensions, simultaneous
field-index/Fourier-leg symmetrization, field-space curvature regressions, and
rank-five $\nabla R$ index checks. The cubic builder produces the three ordered
output-leg vertices required by the cyclic source. The nonlinear
$V^I\leftrightarrow\Pi^I$ field redefinition is included at tree level through
its local products of two-point covariances. The tree-level quantum initial
condition fills all eight canonical three-point blocks on a common subhorizon
slice and is checked against an independent in-in derivation and the de Sitter
transport equation. `CanonicalThreePointTransportSystem` evolves the sourced
three-point function with the three leg covariances, converges in start
depth and tolerance on the reference benchmark, and transforms as a tensor
under a nonlinear field redefinition. This completes P3d. P3e projects the
result to $P_\zeta$, $B_\zeta$, and $f_{\mathrm{NL}}$, reproduces the single-field
consistency relation and the published axion-quartic bispectrum, and
`run_bispectrum` writes convention-labeled CSV output. Within the documented
scope (tree level, during inflation, no reheating) the repository therefore
provides a primordial bispectrum prediction; see
[`docs/curvature-projection.md`](docs/curvature-projection.md).

`CovariantMultifieldSystem<P, G, N>` accepts a `FieldSpaceGeometry<N>` that
supplies $G_{IJ}$, $G^{IJ}$, and $\Gamma^I{}_{JK}$. `FlatFieldSpace` exactly reproduces
the canonical solver, and `HyperbolicFieldSpace` supplies the constant
negative curvature of alpha-attractor and hyperbolic-inflation models, which
is the smallest geometry that exercises a non-vanishing connection and Riemann
tensor at once. For two fields, `kinematics` returns the normalized
adiabatic tangent, continuously oriented entropy normal, `dot_sigma`, and the
signed dimensionless turn rate.

Potentials that can be evaluated with `Jet<N>` automatically implement
`TwiceDifferentiableMultifieldPotential<N>`. A single formula evaluation
returns $V$, $V_{,I}$, and $V_{,IJ}$; the covariant system then constructs
$V_{;IJ}=V_{,IJ}-\Gamma^K{}_{IJ}V_{,K}$ with lower tensor indices and dimension
$M^2$. The shared `Transcendental` trait provides `exp`, `ln`, `cos`, and
`tanh` for `f64`, `Dual`, and `Jet<N>`, so plateau, natural-inflation, and
alpha-attractor formulas use the same automatic-differentiation path as
polynomials.

`TwoFieldPerturbationSystem` evolves the background together with two
independent complex vacuum modes for one comoving `k`. Its mass matrix includes
the covariant potential Hessian, field-space Riemann term, and gravitational
correction. Bunch-Davies initialization checks `k/(aH)` explicitly, and the
canonical Wronskian is available as a numerical invariant.

Mode values are specified as `Quantity<D, Complex64>` so complex conjugation
and Hermitian products remain explicit without sacrificing compile-time
physical dimensions. Real and imaginary parts are separated only for output.
`normalized_mode_error` turns a complex step error into the dimensionless
ratio `|dz| / (a_tol + r_tol |z|)`, where the absolute tolerance carries the
dimension of the component it bounds and the modulus keeps one acceptance
decision for both parts.

`mode_covariance` forms the Hermitian covariance
$\Sigma^{IJ}=\sum_AU^I{}_AU^{J*}{}_A$, and `field_power_spectrum` applies the
$k^3/(2\pi^2)$ normalization, so the $M^{-1}$ covariance and the $M^2$ spectrum cannot be confused and the
normalization cannot be applied twice. `project` contracts the two upper
indices with lowered adiabatic and entropy covectors and returns `P_R`, `P_S`,
the signed `C_RS`, and the correlation coefficient. Hermiticity, the imaginary
residual of each projected spectrum, and any overshoot of `|r_C| <= 1` are
checked against explicit tolerances and reported as typed errors instead of
being silently discarded.

`run_two_field_spectra` puts those pieces together. It surveys the background
to the end of inflation, fixes `k_pivot = a_* H_*` a chosen number of e-folds
earlier, carries every wavenumber to a stated number of e-folds before its own
horizon exit before initializing the vacuum there, and reports spectra at one
output time. A mode that cannot be placed deep enough inside the horizon is a
typed error, not a silently poor approximation, and so is one that has not
left the horizon by the output time, whose vacuum would otherwise be reported
as a spectrum.

A run spans tens of orders of magnitude, so absolute error tolerances are per
group of components and carry each group's dimension: a tolerance meant for
one quantity cannot bound another. The modes are integrated rescaled as
`v = a U` and default to no absolute floor at all, because their error is
normalized by a modulus that does not vanish, and a floor there would become
their accuracy ceiling. `docs/two-field-perturbation-design.md` records the
measurements behind both choices.

## Compile-fail tests

Rejecting invalid expressions is part of the contract, so it is tested. The
cases live in `tests/compile_fail/`. After an intentional change to an error
message, refresh the expectations with:

```sh
TRYBUILD=overwrite cargo test --test compile_fail
```

## Contributing and license

Development and verification steps are documented in
[`CONTRIBUTING.md`](CONTRIBUTING.md). Changes are recorded in
[`CHANGELOG.md`](CHANGELOG.md).

This project is dual-licensed under the
[MIT license](LICENSE-MIT) or [Apache License 2.0](LICENSE-APACHE), at your
option.
