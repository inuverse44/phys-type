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

To verify the complete diagnostics and plotting workflow, install `uv` and run:

```sh
uv run plots/generate_all.py
```

Generated files under `diagnostics/` are intentionally not committed.
