//! Provenance shared by the examples.

use std::process::Command;

/// The short Git revision, marked `-dirty` for uncommitted tracked changes.
pub fn revision() -> Option<String> {
    let output = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()?;
    let revision = output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())?;
    let status = Command::new("git")
        .args(["status", "--porcelain", "--untracked-files=no"])
        .output()
        .ok()?;
    let dirty = status.status.success() && !status.stdout.is_empty();
    Some(if dirty {
        format!("{revision}-dirty")
    } else {
        revision
    })
}
