# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Higgs inflation in the n_s-r plane against Planck 2018 and BICEP/Keck.

Reads ``higgs_ns_r.csv`` (``cargo run --release --example higgs_ns_r``) and the BICEP/Keck XIII Figure 5 contour data
(``python3 scripts/fetch_bk18.py`` puts them in ``data/bk18_rns``). The
contours are the published 68% and 95% regions for Planck 2018
TT,TE,EE+lowE+lensing and for its combination with BK18 and BAO
(arXiv:2110.00483); they are overlaid for comparison, not refitted. The
predictions are exact one-field evolutions for log10 xi from -5 to 4 (colour)
and pivots 50 to 60 Einstein-frame e-folds before the end; thin lines join
N_* = 50 and N_* = 60 along xi. The layout follows Fig. 4 of Kodama and
Takahashi, arXiv:2112.05283: metric (left) and Palatini (right), with
logarithmic (top, down to r = 1e-5) and linear (bottom) r. The contours
are lines, solid for 68% and dashed for 95%. The published axis is
r_0.002; the predictions are r at the pivot, and the tensor tilt of these
models changes r by less than a percent between the two scales. On the
logarithmic panel the contours below the first positive grid row of the
published density (r = 1.2e-3) use its density at r = 0, since the
constraint is an upper bound.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

COLUMNS = {"pivot_e_folds", "formulation", "log10_xi", "n_s", "r"}
CONTOURS = (
    ("P18", "Planck TT,TE,EE+lowE+lensing", "g"),
    ("BK18", "+BK18+BAO", "#006FED"),
)
FORMULATIONS = ("Metric", "Palatini")
LINES = {50.0: ":", 60.0: "-"}
LOG_XI_LIMITS = (-5.0, 4.0)


def contours(axes, data: Path, positive_only: bool) -> list:
    handles = []
    for prefix, label, color in CONTOURS:
        xy = np.loadtxt(data / f"{prefix}density_xy.txt")
        density = np.loadtxt(data / f"{prefix}density_P.txt")
        levels = np.loadtxt(data / f"{prefix}density_levels.txt")
        x, y = xy[:, 0], xy[:, 1]
        if positive_only:
            # The published grid is linear in r. On a logarithmic axis the
            # rows below its first positive r would be empty, although the
            # constraint is an upper bound that extends to r = 0; there the
            # published density at r = 0 is used.
            keep = y > 0.0
            extension = np.geomspace(1.0e-9, y[keep][0], 40, endpoint=False)
            y = np.concatenate([extension, y[keep]])
            density = np.vstack([np.repeat(density[:1, :], extension.size, axis=0), density[keep, :]])
        # The levels file lists the 95% and then the 68% density level.
        axes.contour(
            x, y, density, levels[:2], colors=color, linestyles=["--", "-"], linewidths=1.2, zorder=1,
        )
        handles.append(plt.Line2D([], [], color=color, label=label))
    return handles


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding higgs_ns_r.csv")
    parser.add_argument(
        "--bk18", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "bk18_rns",
        help="BICEP/Keck contour data (scripts/fetch_bk18.py)",
    )
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()
    if not (args.bk18 / "BK18density_P.txt").is_file():
        raise SystemExit(f"missing BICEP/Keck data in {args.bk18}; run scripts/fetch_bk18.py")

    path = args.directory / "higgs_ns_r.csv"
    data = load_csv(path, COLUMNS)
    require_finite(data, COLUMNS - {"formulation"})
    revision = load_metadata(path).get("revision")

    plt.style.use("seaborn-v0_8-white")
    figure, grid = plt.subplots(
        2, 2, figsize=(12.0, 9.5), sharex=True, constrained_layout=True,
    )
    for column, formulation in enumerate(FORMULATIONS):
        rows = data[data["formulation"] == formulation]
        logarithmic, linear = grid[0, column], grid[1, column]
        handles = contours(linear, args.bk18, positive_only=False)
        contours(logarithmic, args.bk18, positive_only=True)
        for axes in (logarithmic, linear):
            for e_folds, line in LINES.items():
                branch = rows[np.isclose(rows["pivot_e_folds"], e_folds)]
                branch = branch[np.argsort(branch["log10_xi"])]
                axes.plot(branch["n_s"], branch["r"], line, color="0.25", linewidth=0.8, zorder=2)
            points = axes.scatter(
                rows["n_s"], rows["r"], c=rows["log10_xi"], cmap="turbo",
                vmin=LOG_XI_LIMITS[0], vmax=LOG_XI_LIMITS[1], s=9, linewidths=0, zorder=3,
            )
        logarithmic.set_yscale("log")
        logarithmic.set_ylim(1.0e-5, 1.0)
        logarithmic.set_ylabel(r"$r$")
        logarithmic.set_title(f"{formulation} formulation")
        linear.set_ylim(0.0, 0.4)
        linear.set_ylabel(r"$r$ ($r_{0.002}$ for the contours)")
        linear.set_xlabel(r"$n_s$")
        linear.set_xlim(0.935, 1.0)
    model_handles = [
        plt.Line2D([], [], color="k", linestyle="-", linewidth=1.2, label="68%"),
        plt.Line2D([], [], color="k", linestyle="--", linewidth=1.2, label="95%"),
        plt.Line2D([], [], color="0.25", linestyle="-", linewidth=0.8, label=r"$N_*=60$"),
        plt.Line2D([], [], color="0.25", linestyle=":", linewidth=0.8, label=r"$N_*=50$"),
    ]
    grid[1, 0].legend(handles=handles + model_handles, fontsize=8, loc="upper right", frameon=False)
    grid[0, 0].text(
        0.02, 0.02,
        r"below $r=1.2\times10^{-3}$ the contours use the published density at $r=0$",
        transform=grid[0, 0].transAxes, fontsize=7, color="0.3",
    )
    bar = figure.colorbar(points, ax=grid, shrink=0.6, pad=0.01)
    bar.set_label(r"$\log_{10}\xi$")
    title = (
        r"Higgs inflation (exact, one field), $N_*=50$–$60$, against Planck 2018 and BICEP/Keck 2018"
        "\ncontours: arXiv:2110.00483 data release, overlaid not refitted; layout of arXiv:2112.05283 Fig. 4"
    )
    if revision is not None:
        title = f"{title} | revision {revision}"
    figure.suptitle(title)
    save_figure(figure, path, "higgs_ns_r", output_dir=args.output)


if __name__ == "__main__":
    main()
