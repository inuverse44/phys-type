//! Hand-written reference models used by the engine's own tests.
//!
//! These are test fixtures, not part of the public API: the nonminimally
//! assisted geometry and potentials of arXiv:2302.05866 written out by hand,
//! and that paper's slow-roll prescription. The generic theory layer is
//! checked against the hand-written geometry, and the engine's validation
//! suites use them until those suites move to the `phys-type-models` crate
//! (step 4b of `docs/architecture.md`). The public versions of these models
//! live in `phys-type-models`.

pub mod assisted_slow_roll;
pub mod nonminimal;

pub use nonminimal::{
    NonminimalAssistedFieldSpace, NonminimalHybridPotential, NonminimalLoopInflationPotential,
    NonminimalPowerLawPotential,
};
