//! The Planck pivot 0.05 Mpc^-1 passed where the code expects its internal
//! wavenumber k / M_pl, whose value depends on the normalization of the scale
//! factor. The observational wavenumber is a distinct type, so this does not
//! compile.

use phys_type::{ComovingWavenumberMpcInverse, Mass};

/// Evolves the mode of internal comoving wavenumber k (in units of M_pl).
fn evolve_mode(k: Mass) -> f64 {
    k.into_value()
}

fn main() {
    let pivot = ComovingWavenumberMpcInverse::try_new(0.05).unwrap();

    let result = evolve_mode(pivot);

    println!("{result}");
}
