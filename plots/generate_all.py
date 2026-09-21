# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Generate every diagnostic CSV and its PNG/PDF visualization."""

import argparse
from pathlib import Path
import subprocess
import sys


PLOTS = (
    ("background.py", "trajectory.csv", "background"),
    ("attractor.py", "attractor.csv", "attractor"),
    ("convergence.py", "convergence.csv", "convergence"),
    ("oscillation.py", "oscillation.csv", "oscillation"),
)

# The two-field figure reads three CSV files, so it takes the directory
# rather than one input file.
TWO_FIELD_PLOT = ("two_field_spectra.py", "two_field_spectra")
FINAL_SPECTRUM_PLOT = ("final_spectra.py", "two_field_final_spectra")
PIVOT_EVOLUTION_PLOT = ("pivot_mode_evolution.py", "two_field_pivot_evolution")


def run(command: list[str], cwd: Path) -> None:
    try:
        subprocess.run(command, cwd=cwd, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"command failed: {' '.join(command)}\n{error}") from error


def verify(output: Path, stem: str) -> None:
    for suffix in ("png", "pdf"):
        artifact = output / f"{stem}.{suffix}"
        if not artifact.is_file() or artifact.stat().st_size == 0:
            raise SystemExit(f"expected plot was not generated: {artifact}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=Path("diagnostics"),
        help="output directory (default: diagnostics)",
    )
    args = parser.parse_args()

    repository = Path(__file__).resolve().parents[1]
    scripts = Path(__file__).resolve().parent
    output = args.output.resolve()

    for example in ("diagnostics", "two_field_spectra"):
        run(
            ["cargo", "run", "--release", "--example", example, "--", str(output)],
            repository,
        )

    for script, csv, stem in PLOTS:
        run([sys.executable, str(scripts / script), str(output / csv)], repository)
        verify(output, stem)

    script, stem = TWO_FIELD_PLOT
    run([sys.executable, str(scripts / script), str(output)], repository)
    verify(output, stem)

    script, stem = FINAL_SPECTRUM_PLOT
    run(
        [sys.executable, str(scripts / script), str(output / "two_field_spectra.csv")],
        repository,
    )
    verify(output, stem)

    script, stem = PIVOT_EVOLUTION_PLOT
    run(
        [sys.executable, str(scripts / script), str(output)],
        repository,
    )
    verify(output, stem)

    print(f"verified {(len(PLOTS) + 3) * 2} plot files in {output}")


if __name__ == "__main__":
    main()
