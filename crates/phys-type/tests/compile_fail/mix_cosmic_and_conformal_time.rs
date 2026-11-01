use phys_type::{ConformalTime, Time};

fn takes_cosmic_time(_: Time) {}

fn main() {
    takes_cosmic_time(ConformalTime::new(1.0));
}
