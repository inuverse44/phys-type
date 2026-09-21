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
