use phys_type::{ComovingWavenumberMpcInverse, Mass};

fn takes_internal_wavenumber(_: Mass) {}

fn main() {
    takes_internal_wavenumber(ComovingWavenumberMpcInverse::try_new(0.05).unwrap());
}
