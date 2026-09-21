# phys-type

Dimension-aware physical quantities for numerical physics, in Rust.

Physical dimensions live in the type system, so a dimensionally invalid
expression fails to compile instead of producing a wrong number. The first
target is single-field inflationary background evolution; see
[`docs/mvp-design.md`](docs/mvp-design.md) for the design.

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
```

The plotting scripts declare their own Python dependencies, so `uv` creates an
isolated environment automatically. The visualization design and physical
checks are documented in
[`docs/visualization-design.md`](docs/visualization-design.md).

## What the library provides

| Module | Contents |
| --- | --- |
| `autodiff` | First-order forward-mode automatic differentiation with `Dual` |
| `dimension` | `MassPower<E>` and the type-level `DimMul`, `DimDiv`, `DimSqrt` |
| `quantity` | `Quantity<D, V>` and its arithmetic |
| `state` | `State`, `StateDerivative`, `StateDelta` |
| `potential` | Automatically differentiated potentials and `QuadraticPotential` |
| `integrate` | The `System` trait, fixed RK4, and embedded Dormand-Prince 5(4) |
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

Canonical multifield systems use fixed-size arrays and the same generic
integrators. `QuadraticMultifieldPotential<N>` supplies one mass per field and
obtains every gradient component through automatic differentiation.

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
