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

The CSV goes to standard output and a summary to standard error. Plot it with
any external tool, for example:

```sh
gnuplot -e "set datafile separator ','; \
            set key autotitle columnhead; \
            plot 'inflation.csv' using 3:4 with lines; pause -1"
```

## What the library provides

| Module | Contents |
| --- | --- |
| `dimension` | `MassPower<E>` and the type-level `DimMul`, `DimDiv`, `DimSqrt` |
| `quantity` | `Quantity<D, V>` and its arithmetic |
| `state` | `State`, `StateDerivative`, `StateDelta` |
| `potential` | The `Potential` trait and `QuadraticPotential` |
| `integrate` | The `System` trait and fixed-step `rk4_step` |
| `scalar_field` | `ScalarFieldSystem`, including the Friedmann constraint |
| `background` | The integration driver, termination conditions, and CSV export |

Units are natural (`c = hbar = 1`) with `M_pl = 1`, so every dimension is a
power of mass. Exponents are stored in half-power steps, which keeps
`sqrt(M^1) = M^(1/2)` representable.

```rust
use phys_type::{Field, Mass, QuadraticPotential, ScalarFieldSystem, ScaleFactor, RunSettings};

let system = ScalarFieldSystem::in_planck_units(QuadraticPotential::new(Mass::new(1.0e-6)));
let initial = system.slow_roll_initial_state(Field::new(16.0), ScaleFactor::new(1.0));
let run = system.run(initial, RunSettings::default());
```

## Compile-fail tests

Rejecting invalid expressions is part of the contract, so it is tested. The
cases live in `tests/compile_fail/`. After an intentional change to an error
message, refresh the expectations with:

```sh
TRYBUILD=overwrite cargo test --test compile_fail
```
