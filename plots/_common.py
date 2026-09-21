"""Shared CSV validation and figure output for the diagnostic plots."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_csv(path: Path, required: set[str]) -> np.ndarray:
    if not path.is_file():
        raise SystemExit(f"input CSV does not exist: {path}")

    try:
        data = np.genfromtxt(
            path,
            delimiter=",",
            names=True,
            dtype=None,
            encoding="utf-8",
        )
    except (OSError, ValueError) as error:
        raise SystemExit(f"could not read {path}: {error}") from error

    if data.dtype.names is None or data.size == 0:
        raise SystemExit(f"input CSV has no data rows: {path}")

    missing = required.difference(data.dtype.names)
    if missing:
        names = ", ".join(sorted(missing))
        raise SystemExit(f"input CSV is missing required columns: {names}")

    return np.atleast_1d(data)


def require_finite(data: np.ndarray, columns: set[str]) -> None:
    for column in columns:
        if not np.all(np.isfinite(data[column])):
            raise SystemExit(f"column contains NaN or infinity: {column}")


def save_figure(figure: plt.Figure, input_path: Path, stem: str) -> None:
    output_dir = input_path.parent
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"

    try:
        figure.savefig(png, dpi=200)
        figure.savefig(pdf)
    except OSError as error:
        raise SystemExit(f"could not write plot output: {error}") from error

    print(f"wrote {png} and {pdf}")
