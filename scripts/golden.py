#!/usr/bin/env python3
"""Golden-data regression check for the example outputs.

Refactors of the core must leave every number unchanged. This script runs the
examples, writes their CSV files into a scratch directory, and compares the
data rows with the files frozen under ``golden/``. Comment lines (``#``), which
carry the Git revision and other provenance, are ignored.

    python3 scripts/golden.py check            # fast examples
    python3 scripts/golden.py check --portable # as in CI, across platforms
    python3 scripts/golden.py check --all      # also the slow examples
    python3 scripts/golden.py check --exact    # require identical text
    python3 scripts/golden.py update [--all]   # refreeze after an intended change

Numbers are compared value by value: |a - b| <= rtol * max(|a|, |b|) with
rtol = 1e-9 by default, so that last-digit differences between platform math
libraries pass while any physical change fails; an absolute floor of
1e-12 times the column's largest value absorbs quantities that vanish up to
round-off. Columns that record
accumulated round-off, such as the imaginary part of an equal-time
correlator, are instead compared on the scale of a reference column. Large files are frozen with a
row stride, recorded in their first line, and the row count is checked too.
Files written once per accepted integration step are only comparable on one
platform: a last-bit difference in a math library shifts the adaptive step
sequence, and with it the e-fold of every row. ``--portable`` skips them and
compares only quantities evaluated at fixed wavenumbers or couplings.

Only the Python standard library is used.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "golden"
MAX_ROWS = 2000

# Written at accepted integration steps; skipped by --portable.
STEP_SAMPLED = {
    "trajectory.csv",
    "attractor.csv",
    "oscillation.csv",
    "two_field_background.csv",
    "two_field_modes.csv",
    "assisted_trajectories.csv",
}

# Round-off diagnostics -> the column whose scale bounds them.
NOISE_COLUMNS = {"Im_P_12": "P_11"}

# (file, column) -> absolute scale, for columns tuned to vanish: the massless
# orbit drives the measured mu^2 / H^2 to zero, leaving only round-off.
ABSOLUTE_SCALES = {("constant_turn_massless.csv", "mu2"): 1.0}

# example name -> CSV files it writes; the slow ones run only with --all.
FAST = {
    "diagnostics": [
        "trajectory.csv",
        "attractor.csv",
        "convergence.csv",
        "oscillation.csv",
    ],
    "two_field_spectra": [
        "two_field_background.csv",
        "two_field_spectra.csv",
        "two_field_horizon_exit_spectra.csv",
        "two_field_modes.csv",
    ],
    "constant_turn_spectra": [
        "constant_turn_massive.csv",
        "constant_turn_massless.csv",
    ],
}
SLOW = {
    "assisted_inflation": ["assisted_predictions.csv", "assisted_trajectories.csv"],
    "constant_turn_squeezed": ["constant_turn_squeezed.csv"],
}


def data_rows(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if not line.startswith("#")]


def frozen(rows: list[str]) -> tuple[int, list[str]]:
    """Header plus every stride-th data row, keeping at most MAX_ROWS."""
    header, body = rows[0], rows[1:]
    stride = max(1, -(-len(body) // MAX_ROWS))
    return stride, [header, *body[::stride]]


def run_examples(examples: dict[str, list[str]], output: Path) -> None:
    for example in examples:
        command = ["cargo", "run", "--quiet", "--release", "--example", example, "--", str(output)]
        print("running", " ".join(command[2:]), flush=True)
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stderr.write(result.stdout + result.stderr)
            raise SystemExit(f"example {example} failed")


def parse(field: str) -> float | str:
    try:
        return float(field)
    except ValueError:
        return field


def compare(name: str, golden: Path, fresh: Path, rtol: float, exact: bool) -> list[str]:
    with golden.open(encoding="utf-8") as handle:
        first = handle.readline()
    stride = int(first.split(":")[1]) if first.startswith("# golden-stride:") else 1
    handle_count = None
    expected = data_rows(golden)
    rows = data_rows(fresh)
    errors = []
    for line in golden.read_text(encoding="utf-8").splitlines():
        if line.startswith("# golden-rows:"):
            handle_count = int(line.split(":")[1])
    if handle_count is not None and handle_count != len(rows) - 1:
        errors.append(f"{name}: {len(rows) - 1} data rows, golden has {handle_count}")
        return errors
    actual = [rows[0], *rows[1:][::stride]]
    if len(actual) != len(expected):
        return [f"{name}: {len(actual)} compared rows, golden has {len(expected)}"]
    if actual[0] != expected[0]:
        return [f"{name}: header changed: {actual[0]!r} != {expected[0]!r}"]
    if exact:
        for index, (a, b) in enumerate(zip(actual, expected)):
            if a != b:
                errors.append(f"{name}: row {index} differs: {a!r} != {b!r}")
                break
        return errors
    table_a = [[parse(f) for f in row.split(",")] for row in actual[1:]]
    table_b = [[parse(f) for f in row.split(",")] for row in expected[1:]]
    columns = expected[0].split(",")
    for column, label in enumerate(columns):
        reference = columns.index(NOISE_COLUMNS[label]) if label in NOISE_COLUMNS else None
        floor = 1.0e-12 * max(
            (abs(row[column]) for row in table_b if isinstance(row[column], float)),
            default=0.0,
        )
        for index, (row_a, row_b) in enumerate(zip(table_a, table_b)):
            a, b = row_a[column], row_b[column]
            if isinstance(a, float) and isinstance(b, float):
                if (name, label) in ABSOLUTE_SCALES:
                    scale = ABSOLUTE_SCALES[(name, label)]
                elif reference is not None:
                    scale = abs(row_b[reference])
                else:
                    scale = max(abs(a), abs(b))
                if abs(a - b) > rtol * scale + floor:
                    errors.append(
                        f"{name}: column {label}, row {index}: {a!r} vs golden {b!r}"
                    )
                    break
            elif a != b:
                errors.append(f"{name}: column {label}, row {index}: {a!r} vs golden {b!r}")
                break
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["check", "update"])
    parser.add_argument("--all", action="store_true", help="include the slow examples")
    parser.add_argument("--exact", action="store_true", help="require identical text")
    parser.add_argument(
        "--portable", action="store_true", help="skip files sampled at integration steps"
    )
    parser.add_argument("--rtol", type=float, default=1.0e-9)
    args = parser.parse_args()

    examples = {**FAST, **(SLOW if args.all else {})}
    with tempfile.TemporaryDirectory() as scratch:
        output = Path(scratch)
        run_examples(examples, output)
        files = [name for names in examples.values() for name in names]
        if args.command == "update":
            GOLDEN.mkdir(exist_ok=True)
            for name in files:
                rows = data_rows(output / name)
                stride, kept = frozen(rows)
                text = f"# golden-stride: {stride}\n# golden-rows: {len(rows) - 1}\n"
                (GOLDEN / name).write_text(text + "\n".join(kept) + "\n", encoding="utf-8")
                print(f"froze {name}: {len(rows) - 1} rows, stride {stride}")
            return 0
        errors = []
        if args.portable:
            files = [name for name in files if name not in STEP_SAMPLED]
        for name in files:
            errors += compare(name, GOLDEN / name, output / name, args.rtol, args.exact)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        print(f"golden data reproduced for {len(files)} files")
        return 0


if __name__ == "__main__":
    sys.exit(main())
