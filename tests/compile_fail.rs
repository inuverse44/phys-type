//! Compile-fail tests.
//!
//! Rejecting dimensionally invalid expressions is part of this crate's
//! contract, so the rejections are tested like any other behaviour. Each file
//! under `tests/compile_fail/` must fail to compile with the error recorded in
//! the matching `.stderr` file.
//!
//! After an intentional change to an error message, refresh the expectations
//! with `TRYBUILD=overwrite cargo test --test compile_fail`.

#[test]
fn dimensionally_invalid_programs_are_rejected() {
    let cases = trybuild::TestCases::new();
    cases.compile_fail("tests/compile_fail/*.rs");
}
