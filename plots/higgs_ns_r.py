# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Higgs inflation in the n_s-r plane against Planck 2018 and BICEP/Keck.

Reads ``higgs_inflation.csv`` (``cargo run --release --example
higgs_inflation``) and the BICEP/Keck XIII Figure 5 contour data
(``python3 scripts/fetch_bk18.py`` puts them in ``data/bk18_rns``). The
contours are the published 68% and 95% regions for Planck 2018
TT,TE,EE+lowE+lensing and for its combination with BK18 and BAO
(arXiv:2110.00483); they are overlaid for comparison, not refitted. The
predictions are exact one-field evolutions with the pivot 50 (dashed) or 60
(solid) Einstein-frame e-folds before the end, for the metric and Palatini
formulations, along the nonminimal coupling xi. The published axis is
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

COLUMNS = {"pivot_e_folds", "formulation", "xi", "n_s", "r"}
CONTOURS = (
    ("P18", "Planck TT,TE,EE+lowE+lensing", ("#9ACB9A", "g"), 1.0),
    ("BK18", "+BK18+BAO", ("#8CD3F5", "#006FED"), 0.8),
)
STYLES = {"Metric": "tab:red", "Palatini": "tab:purple"}
LINES = {50.0: "--", 60.0: "-"}


def contours(axes, data: Path, positive_only: bool) -> list:
    handles = []
    for prefix, label, colors, alpha in CONTOURS:
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
        axes.contourf(x, y, density, levels, colors=colors, alpha=alpha)
        axes.contour(x, y, density, levels[:1], colors=colors[1], alpha=alpha)
        handles.append(plt.Rectangle((0, 0), 1, 1, color=colors[0], label=label))
    return handles


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding higgs_inflation.csv")
    parser.add_argument(
        "--bk18", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "bk18_rns",
        help="BICEP/Keck contour data (scripts/fetch_bk18.py)",
    )
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()
    if not (args.bk18 / "BK18density_P.txt").is_file():
        raise SystemExit(f"missing BICEP/Keck data in {args.bk18}; run scripts/fetch_bk18.py")

    path = args.directory / "higgs_inflation.csv"
    data = load_csv(path, COLUMNS)
    require_finite(data, COLUMNS - {"formulation"})
    revision = load_metadata(path).get("revision")

    plt.style.use("seaborn-v0_8-white")
    figure, (linear, logarithmic) = plt.subplots(1, 2, figsize=(12.5, 5.6), constrained_layout=True)
    handles = contours(linear, args.bk18, positive_only=False)
    contours(logarithmic, args.bk18, positive_only=True)
    for formulation, color in STYLES.items():
        for e_folds, line in LINES.items():
            rows = data[(data["formulation"] == formulation) & np.isclose(data["pivot_e_folds"], e_folds)]
            rows = rows[np.argsort(rows["xi"])]
            for axes in (linear, logarithmic):
                axes.plot(rows["n_s"], rows["r"], line, color=color, linewidth=1.4)
                axes.plot(rows["n_s"], rows["r"], "o", color=color, markersize=4)
            if e_folds == 60.0:
                # The two formulations coincide at small xi, so the metric
                # branch is labelled only at its large-xi end.
                labelled = (1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e4) if formulation == "Palatini" else (1.0e4,)
                for row in rows:
                    if row["xi"] in labelled:
                        logarithmic.annotate(
                            rf"$\xi={row['xi']:g}$", (row["n_s"], row["r"]),
                            textcoords="offset points", xytext=(5, 3), fontsize=7, color=color,
                        )
    for axes in (linear, logarithmic):
        axes.set_xlim(0.945, 1.0)
        axes.set_xlabel(r"$n_s$")
    linear.set_ylim(0.0, 0.26)
    linear.set_ylabel(r"$r_{0.002}$ (contours), $r$ at the pivot (Higgs)")
    linear.set_title("as in BICEP/Keck XIII, Fig. 5")
    logarithmic.set_yscale("log")
    logarithmic.set_ylim(1.0e-8, 0.3)
    logarithmic.set_ylabel(r"$r$")
    logarithmic.set_title(
        r"logarithmic $r$; contours below $r=1.2\times10^{-3}$ use the density at $r=0$",
        fontsize=10,
    )
    model_handles = [
        plt.Line2D([], [], color=STYLES["Metric"], label="Higgs, metric"),
        plt.Line2D([], [], color=STYLES["Palatini"], label="Higgs, Palatini"),
        plt.Line2D([], [], color="0.3", linestyle="-", label=r"$N_*=60$"),
        plt.Line2D([], [], color="0.3", linestyle="--", label=r"$N_*=50$"),
    ]
    linear.legend(handles=handles + model_handles, fontsize=8, loc="upper left", frameon=False)
    title = (
        "Higgs inflation (exact, one field) against Planck 2018 and BICEP/Keck 2018\n"
        "contours: arXiv:2110.00483 data release, overlaid not refitted"
    )
    if revision is not None:
        title = f"{title} | revision {revision}"
    figure.suptitle(title)
    save_figure(figure, path, "higgs_ns_r", output_dir=args.output)


if __name__ == "__main__":
    main()
