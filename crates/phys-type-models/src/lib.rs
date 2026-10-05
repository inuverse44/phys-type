//! Models built on the `phys-type` core.
//!
//! The core knows no model: it states theories through its theory layer and
//! evolves them exactly. This crate holds particular theories, the formulas
//! of particular papers, and benchmark backgrounds, and its examples produce
//! the data behind the tracked figures. See `docs/architecture.md`.
//!
//! - [`assisted`]: nonminimally assisted inflation (arXiv:2302.05866) as a
//!   generalized nonminimal theory, with loop, power-law, and hybrid
//!   inflatons;
//! - [`hkkpt`]: that paper's slow-roll delta N prescription;
//! - [`higgs`]: single-field Higgs inflation, metric and Palatini;
//! - [`constant_turn`]: the circular orbit behind the exact constant-turn
//!   solutions of arXiv:2606.18248 and 2607.14529.

pub mod assisted;
pub mod constant_turn;
pub mod higgs;
pub mod hkkpt;
