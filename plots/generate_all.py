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


def run(command: list[str], cwd: Path) -> None:
    try:
        subprocess.run(command, cwd=cwd, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"command failed: {' '.join(command)}\n{error}") from error


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

    run(
        ["cargo", "run", "--release", "--example", "diagnostics", "--", str(output)],
        repository,
    )

    for script, csv, stem in PLOTS:
        run([sys.executable, str(scripts / script), str(output / csv)], repository)
        for suffix in ("png", "pdf"):
            artifact = output / f"{stem}.{suffix}"
            if not artifact.is_file() or artifact.stat().st_size == 0:
                raise SystemExit(f"expected plot was not generated: {artifact}")

    print(f"verified {len(PLOTS) * 2} plot files in {output}")


if __name__ == "__main__":
    main()
