//! Field-space tensors in polar coordinates (rho, theta), G = diag(1, rho^2).
//!
//! At rho = 2 with velocity (0.1, 0.3) and gradient V_,I = (1e-3, 4e-3), the
//! typed code computes the kinetic term G_IJ phi_dot^I phi_dot^J and the force
//! G^IJ V_,J. The hand-written shortcuts that are right in Cartesian
//! coordinates, phi_dot^I phi_dot^I and V_,I, are printed alongside: they are
//! wrong here, and the broken examples show that they do not compile.
//!
//!     cargo run --example tensor
//!     cargo build --example forget_metric      # phi_dot^I phi_dot^I
//!     cargo build --example gradient_as_force  # V_,I added to an upper index
//!     cargo build --example raise_with_metric  # G_IJ used to raise an index

use phys_type::{
    ContravariantField, CovariantField, Field, FieldSpaceGeometry, M2, M3, M4, PolarFieldSpace,
    Quantity,
};

fn main() {
    let position = [Field::new(2.0), Field::new(0.5)]; // (rho, theta)
    let velocity = ContravariantField::<M2, f64, 2>::new([Quantity::new(0.1), Quantity::new(0.3)]);
    let gradient = CovariantField::<M3, f64, 2>::new([Quantity::new(1.0e-3), Quantity::new(4.0e-3)]);

    let metric = PolarFieldSpace.metric(&position); // G_IJ
    let inverse = PolarFieldSpace.inverse_metric(&position); // G^IJ

    // sigma_dot^2 = G_IJ phi_dot^I phi_dot^J: lower one index, then contract.
    let kinetic: Quantity<M4, f64> = metric.lower(&velocity).contract(&velocity);
    // The force term of the equation of motion has an upper index: G^IJ V_,J.
    let force: ContravariantField<M3, f64, 2> = inverse.raise(&gradient);

    let [v1, v2] = velocity.into_components().map(Quantity::into_value);
    let [g1, g2] = gradient.into_components().map(Quantity::into_value);
    let [f1, f2] = force.into_components().map(Quantity::into_value);

    println!("sigma_dot^2 = G_IJ v^I v^J   = {:.4}", kinetic.into_value());
    println!("hand-written  v^I v^I        = {:.4}", v1 * v1 + v2 * v2);
    println!("force       = G^IJ V_,J      = ({f1:.1e}, {f2:.1e})");
    println!("hand-written  V_,I           = ({g1:.1e}, {g2:.1e})");
}
