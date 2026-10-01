# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Reproduce the assisted-inflation figures of arXiv:2302.05866.

Reads ``assisted_predictions.csv`` and ``assisted_trajectories.csv`` written by
``cargo run --release --example assisted_inflation``. Rows are the loop
(Sec. 4.1), power-law (Sec. 4.2), and hybrid (Sec. 4.3) models. Columns are
the field trajectory from the pivot to the end at n_s = 0.965 (Figs. 9, 13,
17), n_s against xi_2, and
f_NL against n_s for s_* = 0.001 and 0.1 (Figs. 10, 14, 18). f_NL is drawn in
the sign of the reference, which is minus the Maldacena/Planck f_NL stored in
the CSV. Solid lines with markers are the exact evolution; dotted lines are
Eqs. (3.1) and (3.5) with the numerical phi_*, s_e, and end surface; dashed
lines are the reference's own slow-roll prescription.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from _common import load_csv, load_metadata, require_finite, save_figure

PREDICTION_COLUMNS = {
    "model", "s_star", "xi", "n_s", "n_s_eq_3_1", "n_s_paper",
    "f_nl", "f_nl_eq_3_5", "f_nl_paper",
}
TRAJECTORY_COLUMNS = {"model", "s_star", "xi", "n_s", "N", "phi", "s"}
FAMILIES = (
    ("loop", "Loop inflation (Class I)"),
    ("power-law", "Power-law inflation (Class II)"),
    ("hybrid", "Hybrid inflation (Class II)"),
)
COLORS = ("tab:red", "tab:blue")
PLANCK_BAND = (0.955, 0.976)


def label(model: str) -> str:
    return (
        model.replace("loop ", "")
        .replace("power-law ", "")
        .replace("hybrid ", "")
        .replace("alpha", r"$\alpha$")
        .replace("lambda", r"$\lambda$")
        .replace("mu", r"$\mu$")
        .replace("phi_e", r"$\phi_e$")
        .replace(" ", ", ")
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding both CSV files")
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()

    predictions_path = args.directory / "assisted_predictions.csv"
    predictions = load_csv(predictions_path, PREDICTION_COLUMNS)
    trajectories = load_csv(args.directory / "assisted_trajectories.csv", TRAJECTORY_COLUMNS)
    require_finite(predictions, PREDICTION_COLUMNS - {"model"})
    revision = load_metadata(predictions_path).get("revision")

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(3, 4, figsize=(17.0, 11.5), constrained_layout=True)
    for row, (family, title) in zip(axes, FAMILIES):
        trajectory, tilt, large, small = row
        models = sorted(
            {str(m) for m in predictions["model"] if str(m).startswith(family)}
        )
        for color, model in zip(COLORS, models):
            for s_star, style in ((1.0e-3, "-"), (0.1, "--")):
                path = trajectories[
                    (trajectories["model"] == model)
                    & np.isclose(trajectories["s_star"], s_star)
                ]
                if path.size:
                    # From the pivot, 60 e-folds before the end, as in the reference.
                    path = path[path["N"] >= path["N"].max() - 60.0]
                    trajectory.plot(
                        path["phi"], path["s"], style, color=color, linewidth=1.6,
                        label=rf"{label(model)}, $s_*={s_star:g}$, "
                        rf"$\xi_2={path['xi'][0]:.4f}$",
                    )
                points = predictions[
                    (predictions["model"] == model)
                    & np.isclose(predictions["s_star"], s_star)
                ]
                if not points.size:
                    continue
                points = points[np.argsort(points["xi"])]
                marker = "o" if s_star < 0.01 else "s"
                tilt.plot(
                    points["xi"], points["n_s"], "-" + marker, color=color,
                    markerfacecolor="white", markersize=5, linewidth=1.6,
                )
                tilt.plot(points["xi"], points["n_s_eq_3_1"], ":", color=color, linewidth=1.6)
                tilt.plot(points["xi"], points["n_s_paper"], "--", color=color, linewidth=1.0)
                tilt.annotate(
                    rf"$s_*={s_star:g}$", (points["xi"][-1], points["n_s"][-1]),
                    textcoords="offset points", xytext=(5, 0), fontsize=8,
                    color=color, ha="left", va="center",
                )
                target = large if s_star < 0.01 else small
                target.plot(
                    points["n_s"], -points["f_nl"], "-" + marker, color=color,
                    markerfacecolor="white", markersize=5, linewidth=1.6,
                    label=label(model),
                )
                target.plot(points["n_s_eq_3_1"], -points["f_nl_eq_3_5"], ":", color=color, linewidth=1.6)
                target.plot(points["n_s_paper"], -points["f_nl_paper"], "--", color=color, linewidth=1.0)

        if family == "hybrid":
            trajectory.set_yscale("log")
        trajectory.set_title(f"{title}: pivot to end at $n_s\\simeq0.965$", fontsize=10)
        trajectory.set_xlabel(r"$\phi/M_{\mathrm{pl}}$")
        trajectory.set_ylabel(r"$s/M_{\mathrm{pl}}$")
        trajectory.legend(fontsize=7, loc="best")
        tilt.axhline(0.965, color="0.4", linewidth=1.0, linestyle=":")
        tilt.axhspan(*PLANCK_BAND, color="0.9", zorder=0)
        tilt.set_title(r"$n_s$ at the end of inflation", fontsize=10)
        tilt.set_xlabel(r"$\xi_2$")
        tilt.set_ylabel(r"$n_s$")
        for target, s_star in ((large, 0.001), (small, 0.1)):
            target.axvspan(*PLANCK_BAND, color="0.9", zorder=0)
            target.set_title(rf"$f_{{\mathrm{{NL}}}}$, $s_*={s_star:g}$", fontsize=10)
            target.set_xlabel(r"$n_s$")
            target.set_ylabel(r"$f_{\mathrm{NL}}$, sign of the reference")
            target.legend(fontsize=7, loc="best")

    styles = [
        Line2D([], [], color="0.2", marker="o", markerfacecolor="white", label="exact (phys-type)"),
        Line2D([], [], color="0.2", linestyle=":", label="Eqs. (3.1)/(3.5), numerical inputs"),
        Line2D([], [], color="0.2", linestyle="--", linewidth=1.0, label="slow-roll prescription of the reference"),
        Line2D([], [], color="0.2", linestyle="-", label=r"trajectory: $s_*=0.001$"),
        Line2D([], [], color="0.2", linestyle="--", linewidth=1.6, label=r"trajectory: $s_*=0.1$"),
    ]
    figure.legend(handles=styles, loc="outside lower center", ncol=5)
    title = (
        "Nonminimally assisted inflation (arXiv:2302.05866): exact evolution, "
        "evaluated $10^{-3}$ e-folds before the end"
    )
    if revision is not None:
        title = f"{title}\nrevision {revision}"
    figure.suptitle(title)
    save_figure(figure, predictions_path, "assisted_inflation", output_dir=args.output)


if __name__ == "__main__":
    main()
