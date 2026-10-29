//! epsilon_H = phi_dot^2 / (2 H^2 M_pl^2), in reduced Planck units and in GeV.
//!
//! With M_pl = 1 the hand-written phi_dot^2 / (2 H^2) gives the same number,
//! so no test in Planck units notices the missing M_pl^2. In GeV it is off by
//! M_pl^2 = 5.9e36. The typed formula is dimensionless in both; the
//! hand-written one has dimension M^2 and does not compile
//! (--example forgot_planck_mass).

use phys_type::{Dimensionless, FieldVelocity, Hubble, Mass};

fn epsilon(phi_dot: FieldVelocity, hubble: Hubble, m_pl: Mass) -> Dimensionless {
    phi_dot * phi_dot / (hubble * hubble * m_pl * m_pl * 2.0)
}

fn main() {
    let reduced_planck_mass_gev = 2.435e18;
    for (units, scale) in [("M_pl = 1", 1.0), ("GeV     ", reduced_planck_mass_gev)] {
        let m_pl = Mass::new(scale);
        let hubble = Hubble::new(1.0e-5 * scale);
        let phi_dot = FieldVelocity::new(1.0e-6 * scale * scale);

        let typed = epsilon(phi_dot, hubble, m_pl).into_value();
        let (v, h) = (phi_dot.into_value(), hubble.into_value());
        let hand_written = v * v / (2.0 * h * h);
        println!("{units}: epsilon_H = {typed:.3e}, hand-written phi_dot^2/(2H^2) = {hand_written:.3e}");
    }
}
