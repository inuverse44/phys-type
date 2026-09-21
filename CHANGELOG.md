# Changelog

All notable changes to this project will be documented in this file. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Fixed-size field-space vectors and rank-two tensors with compile-time upper
  and lower indices.
- Explicit metric raise/lower operations and opposite-variance contraction.
- Compile-fail coverage for invalid field-space index operations.
- Covariant multifield background equations with pluggable metric, inverse
  metric, and Christoffel symbols.
- An oriented two-field adiabatic/entropy basis and signed turn rate.
- Second-order `Jet<N>` automatic differentiation for a potential's value,
  complete gradient, and Hessian in one evaluation.
- Covariant potential Hessians including the Christoffel correction.
- Dense rank-three and rank-four field-space tensors with typed variance and
  physical dimension.
- Slot-specific first-index raising/lowering for higher-rank tensors, with
  compile-fail coverage for variance errors and connection/tensor confusion.

## [0.1.0] - 2026-09-21

### Added

- Compile-time mass dimensions with half-integer exponents.
- Dimension-aware quantities, states, derivatives, and state changes.
- Forward-mode automatic differentiation for potential gradients.
- Fixed-step RK4 and adaptive Dormand–Prince 5(4) integration.
- Single-field and canonical const-generic multifield backgrounds.
- Compile-fail tests for dimensionally invalid expressions.
- CSV diagnostics and reproducible Matplotlib visualization workflow.

[Unreleased]: https://github.com/inuverse44/phys-type/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/inuverse44/phys-type/releases/tag/v0.1.0
