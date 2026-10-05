# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
#   "scipy>=1.14,<2",
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
CONSTANT_TURN_PLOT = ("constant_turn.py", "constant_turn_spectra")
ASSISTED_PLOT = ("assisted_inflation.py", "assisted_inflation")
SQUEEZED_PLOT = ("constant_turn_squeezed.py", "constant_turn_squeezed")
HIGGS_PLOT = ("higgs_inflation.py", "higgs_inflation")


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

    # The assisted-inflation grid runs the three-point transport for 72
    # couplings and takes about five minutes; the squeezed constant-turn
    # triangles take about nine.
    for example in (
        "diagnostics",
        "two_field_spectra",
        "constant_turn_spectra",
        "assisted_inflation",
        "constant_turn_squeezed",
        "higgs_inflation",
    ):
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

    for script, stem in (CONSTANT_TURN_PLOT, ASSISTED_PLOT, SQUEEZED_PLOT, HIGGS_PLOT):
        run([sys.executable, str(scripts / script), str(output)], repository)
        verify(output, stem)

    print(f"verified {(len(PLOTS) + 7) * 2} plot files in {output}")


if __name__ == "__main__":
    main()
