//! Dimension-aware physical quantities for numerical physics.
//!
//! The crate encodes physical dimensions in the type system so that
//! dimensionally invalid expressions fail to compile. See `docs/mvp-design.md`
//! for the design that this implements.

pub mod autodiff;
pub mod background;
pub mod dimension;
pub mod field_space;
pub mod integrate;
pub mod multifield;
pub mod perturbation;
pub mod potential;
pub mod quantity;
pub mod scalar_field;
pub mod spectra;
pub mod spectrum_run;
pub mod state;

pub use autodiff::{Dual, Jet};
pub use background::{
    AdaptiveRunSettings, AdaptiveStep, BackgroundRun, Record, RunSettings, Termination,
};
pub use dimension::{
    DimDiv, DimMul, DimSqrt, M0, M1, M2, M3, M3Over2, M4, M5Over2, MHalf, MMinus1, MMinus2,
    MMinusHalf, MassPower,
};
pub use field_space::{
    ChristoffelSymbols, ContravariantField, CovariantField, CovariantPotentialHessian,
    CurvedFieldSpaceGeometry, FieldMetric, FieldSpaceGeometry, FieldTensor2, FieldTensor3,
    FieldTensor4, FieldVector, FlatFieldSpace, GeometryValidationError, HyperbolicFieldSpace,
    InverseFieldMetric, KinematicBasis2, KinematicBasisError, Lower, LowerFirst,
    LoweredRiemannTensor, RaiseFirst, RiemannIdentity, RiemannTensor, Tensor3Components,
    Tensor4Components, Upper, validate_riemann_symmetries,
};
pub use integrate::{
    AdaptiveAttempt, AdaptiveSystem, EmbeddedStep, ErrorTolerances, StepSizeLimits, System,
    adaptive_step, component_error_ratio, next_step_size, rk4_step, rk45_step,
};
pub use multifield::{
    CovariantMultifieldSystem, MultifieldState, MultifieldStateDelta, MultifieldStateDerivative,
    MultifieldSystem, TwoFieldKinematics,
};
pub use num_complex::Complex64;
pub use perturbation::{
    ComplexMagnitude, ComplexMatrix2, ModeAcceleration, ModeAmplitude, ModeInitializationError,
    ModeVelocity, TwoFieldModes, TwoFieldPerturbationDelta, TwoFieldPerturbationDerivative,
    TwoFieldPerturbationState, TwoFieldPerturbationSystem, TwoFieldPerturbationTolerance,
    normalized_mode_error,
};
pub use potential::{
    MultifieldPotential, MultifieldPotentialValue, Potential, PotentialDerivatives, PotentialValue,
    QuadraticMultifieldPotential, QuadraticPotential, TwiceDifferentiableMultifieldPotential,
};
pub use quantity::{Quantity, Sqrt};
pub use scalar_field::ScalarFieldSystem;
pub use spectra::{
    DimensionlessSpectrum, FieldPowerSpectrum, ModeCovariance, ModeProjections, ProjectedSpectra,
    SpectrumError, SpectrumTolerances, field_power_spectrum, hermiticity_residual, mode_covariance,
    project,
};
pub use spectrum_run::{
    BackgroundRecord, ModeEvolutionRecord, Pivot, SpectrumRecord, SpectrumRun, SpectrumRunError,
    SpectrumRunSettings, SpectrumRunStage, write_metadata,
};
pub use state::{State, StateDelta, StateDerivative};

/// A dimensionless quantity, such as the scale factor.
pub type Dimensionless<V = f64> = Quantity<M0, V>;
/// The scale factor `a`.
pub type ScaleFactor<V = f64> = Quantity<M0, V>;
/// A scalar field value `phi`.
pub type Field<V = f64> = Quantity<M1, V>;
/// A scalar field velocity `dphi/dt`.
pub type FieldVelocity<V = f64> = Quantity<M2, V>;
/// The Hubble rate `H`.
pub type Hubble<V = f64> = Quantity<M1, V>;
/// An energy density, such as `V(phi)`.
pub type EnergyDensity<V = f64> = Quantity<M4, V>;
/// The gradient of a potential, `dV/dphi`.
pub type PotentialGradient<V = f64> = Quantity<M3, V>;
/// An entry in the coordinate or covariant Hessian of a potential.
pub type PotentialHessian<V = f64> = Quantity<M2, V>;
/// A mass, such as the reduced Planck mass.
pub type Mass<V = f64> = Quantity<M1, V>;
/// Cosmic time `t`.
pub type Time<V = f64> = Quantity<MMinus1, V>;
