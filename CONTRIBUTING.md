# Contributing

## Branches

Create one short-lived branch per independently reviewable change:

- `feature/<name>` for user-visible capabilities
- `fix/<name>` for corrections
- `docs/<name>` for documentation-only work
- `chore/<name>` for maintenance and release preparation

Branch from the latest `main`. Keep commits focused, push the feature branch,
and merge only after the checks below pass. Prefer fast-forward integration for
local linear work; use a reviewed pull request when collaborating.

## Toolchain

Install the `llvm-tools` component together with `clippy` and `rustfmt`:

```sh
rustup component add clippy rustfmt llvm-tools
```

Since Rust 1.98, rustc strips release binaries on macOS with its bundled
`rust-objcopy`, which loads `libLLVM.dylib` from
`lib/rustlib/<host>/lib/`, the location shipped by `llvm-tools`. Without the
component it relies on the `DYLD_FALLBACK_LIBRARY_PATH` set by the rustup
proxies, so invoking a toolchain's `cargo` directly (for example with a
Homebrew-installed rustup whose proxies are not on `PATH`) prints
`stripping debug info with rust-objcopy failed: signal: 6 (SIGABRT)` on every
release build. The binaries still work, but the warning hides real ones.

## Required checks

```sh
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo build --release --all-targets
cargo package --locked
```

## Golden data

`golden/` freezes the data rows of the example CSV files. A refactor must
reproduce them. CI runs

```sh
python3 scripts/golden.py check --portable
```

on the fast examples, comparing every value to a relative tolerance of
`1e-9`; round-off diagnostics such as `Im_P_12` are compared on the scale of
their reference column. `--portable` skips files written once per accepted
integration step (background and mode histories): a last-bit difference in
the platform math library shifts the adaptive step sequence, so those rows
are comparable only on one platform. Before merging a change that touches the
engine, run the full and the identical-text checks locally, on the platform
where the golden data were frozen (currently macOS, aarch64):

```sh
python3 scripts/golden.py check --all
python3 scripts/golden.py check --exact
```

When a change is meant to alter the numbers, refreeze with
`python3 scripts/golden.py update --all` and explain the difference in the
pull request.

To verify the complete diagnostics and plotting workflow, install `uv` and run:

```sh
uv run plots/generate_all.py
```

Generated files under `diagnostics/` are intentionally not committed.
