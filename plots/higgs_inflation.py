# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Higgs inflation, metric and Palatini, evolved exactly (arXiv:0710.3755).

Reads ``higgs_inflation.csv`` written by
``cargo run --release --example higgs_inflation``. Markers are the exact
one-field evolution with the pivot 60 Einstein-frame e-folds before the end;
thin lines are first-order slow roll from the same Einstein-frame potential;
dashed lines are the large-xi limits n_s = 1 - 2/N, r = 12/N^2 (metric) and
r = 2/(xi N^2) (Palatini). The last panel shows how many more e-folds the same
interval spans in the Jordan frame.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

COLUMNS = {"pivot_e_folds", "formulation", "xi", "n_s", "r", "n_s_slow_roll", "r_slow_roll", "jordan_e_folds"}
STYLES = {"Metric": ("tab:blue", "o"), "Palatini": ("tab:orange", "s")}
E_FOLDS = 60.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding the CSV")
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()

    path = args.directory / "higgs_inflation.csv"
    data = load_csv(path, COLUMNS)
    require_finite(data, COLUMNS - {"formulation"})
    data = data[np.isclose(data["pivot_e_folds"], E_FOLDS)]
    revision = load_metadata(path).get("revision")

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, (tilt, ratio, frames) = plt.subplots(
        1, 3, figsize=(15.0, 4.8), constrained_layout=True
    )
    grid = np.geomspace(data["xi"].min(), data["xi"].max(), 200)
    for formulation, (color, marker) in STYLES.items():
        rows = data[data["formulation"] == formulation]
        rows = rows[np.argsort(rows["xi"])]
        tilt.semilogx(rows["xi"], rows["n_s_slow_roll"], "-", color=color, linewidth=1.0)
        tilt.semilogx(
            rows["xi"], rows["n_s"], marker, color=color, markerfacecolor="white",
            markersize=6, label=f"{formulation}, exact",
        )
        ratio.loglog(rows["xi"], rows["r_slow_roll"], "-", color=color, linewidth=1.0)
        ratio.loglog(
            rows["xi"], rows["r"], marker, color=color, markerfacecolor="white",
            markersize=6, label=f"{formulation}, exact",
        )
        frames.semilogx(
            rows["xi"], rows["jordan_e_folds"] - E_FOLDS, marker + "-", color=color,
            markerfacecolor="white", markersize=6, label=formulation,
        )
    tilt.axhline(1.0 - 2.0 / E_FOLDS, color="0.4", linestyle="--", linewidth=1.0,
                 label=r"$1-2/N$")
    tilt.set_xlabel(r"$\xi$")
    tilt.set_ylabel(r"$n_s$")
    tilt.set_title(r"$n_s$ (lines: first-order slow roll)")
    tilt.legend(fontsize=8)
    ratio.loglog(grid, np.full_like(grid, 12.0 / E_FOLDS**2), "--", color="tab:blue",
                 linewidth=1.0, label=r"$12/N^2$")
    large = grid[grid >= 0.3]
    ratio.loglog(large, 2.0 / (large * E_FOLDS**2), "--", color="tab:orange", linewidth=1.0,
                 label=r"$2/(\xi N^2)$")
    ratio.set_xlabel(r"$\xi$")
    ratio.set_ylabel(r"$r$")
    ratio.set_title(r"$r=\mathcal{P}_t/\mathcal{P}_{\mathcal{R}}$")
    ratio.legend(fontsize=8)
    frames.set_xlabel(r"$\xi$")
    frames.set_ylabel(r"$N_J-N_E$ from the pivot to the end")
    frames.set_title(r"Jordan- minus Einstein-frame e-folds ($N_E=60$)")
    frames.legend(fontsize=8)

    title = "Higgs inflation evolved exactly (one field; Einstein-frame pivot 60 e-folds before the end)"
    if revision is not None:
        title = f"{title} | revision {revision}"
    figure.suptitle(title)
    save_figure(figure, path, "higgs_inflation", output_dir=args.output)


if __name__ == "__main__":
    main()
