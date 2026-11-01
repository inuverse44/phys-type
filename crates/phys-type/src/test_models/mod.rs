//! Hand-written reference models used by the engine's own tests.
//!
//! These are test fixtures, not part of the public API: the nonminimally
//! assisted geometry and loop potential of arXiv:2302.05866 written out by
//! hand. Their geometry has a nonzero covariant Riemann derivative in closed
//! form, so the generic `JetGeometry` and the theory layer are checked
//! against it. The models themselves, and their validation against the
//! literature, live in the `phys-type-models` crate.

pub mod nonminimal;

pub use nonminimal::{NonminimalAssistedFieldSpace, NonminimalLoopInflationPotential};
