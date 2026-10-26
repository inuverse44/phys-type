# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Two nonminimally coupled fields off their ridge (arXiv:1210.7487, KMS).

Reads ``kms_trajectories.csv`` and ``kms_f_nl.csv`` written by
``cargo run --release --example kms_trajectories``. The panels are the
counterparts of KMS Figs. 2-3 (chi leaving the ridge), Fig. 7 (the turn
rate), and Fig. 11 (|f_NL| from the crossing-time delta N against the initial
chi). Grey points in the last panel are those where epsilon_H reaches one
during the turn, so that the survey ends at that first crossing.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

TRAJECTORY_COLUMNS = {"trajectory", "chi0", "e_folds_before_end", "chi", "epsilon_h", "eta_perp"}
F_NL_COLUMNS = {"chi0", "e_folds_total", "f_nl"}
STYLES = {1: ("tab:orange", ":"), 2: ("tab:red", "-"), 3: ("black", "--")}
# KMS Fig. 7: the turn rate peaks 34.5 and 4.9 e-folds before the end.
KMS_PEAKS = {1: 34.5, 2: 4.9}
# KMS Sec. V: |f_NL| = 43.3 for trajectory 2.
KMS_F_NL = 43.3
# Surveys that end at the first epsilon_H = 1 during the turn fall short of
# the 71.4-71.5 e-folds of an uninterrupted run.
INTERRUPTED_BELOW = 71.4


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding the CSV files")
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()

    path = args.directory / "kms_trajectories.csv"
    trajectories = load_csv(path, TRAJECTORY_COLUMNS)
    require_finite(trajectories, TRAJECTORY_COLUMNS - {"eta_perp"})
    scan = load_csv(args.directory / "kms_f_nl.csv", F_NL_COLUMNS)
    require_finite(scan, F_NL_COLUMNS)
    revision = load_metadata(path).get("revision")

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.6), constrained_layout=True)
    (ridge, turn), (slow_roll, f_nl) = axes

    for index, (color, style) in STYLES.items():
        rows = trajectories[trajectories["trajectory"] == index]
        before = rows["e_folds_before_end"]
        label = rf"trajectory {index}, $\chi_0={rows['chi0'][0]:.1e}$"
        ridge.semilogy(before, np.abs(rows["chi"]), style, color=color, label=label)
        valid = np.isfinite(rows["eta_perp"]) & (rows["eta_perp"] > 0.0)
        turn.semilogy(before[valid], rows["eta_perp"][valid], style, color=color,
                      linewidth=1.0, label=f"trajectory {index}")
        slow_roll.plot(before, rows["epsilon_h"], style, color=color, linewidth=1.0,
                       label=f"trajectory {index}")
    for index, before in KMS_PEAKS.items():
        turn.axvline(before, color=STYLES[index][0], linewidth=0.8, alpha=0.6)
    for axis in (ridge, turn, slow_roll):
        axis.set_xlim(72.0, 0.0)
        axis.set_xlabel("e-folds before the end")

    ridge.set_ylabel(r"$|\chi|/M_{\mathrm{pl}}$")
    ridge.set_title(r"leaving the ridge $\chi=0$ (KMS Figs. 2-3)")
    ridge.legend(fontsize=8)
    turn.set_ylabel(r"$|\eta_\perp|=|D_t\hat\sigma|/H$")
    turn.set_title("turn rate (KMS Fig. 7; vertical lines: KMS peaks)")
    turn.legend(fontsize=8)
    slow_roll.axhline(1.0, color="0.4", linestyle="--", linewidth=1.0)
    slow_roll.set_ylim(0.0, 1.05)
    slow_roll.set_ylabel(r"$\epsilon_H$")
    slow_roll.set_title(r"$\epsilon_H$: trajectory 2 reaches 0.99 at the turn")
    slow_roll.legend(fontsize=8)

    scan = scan[np.argsort(scan["chi0"])]
    interrupted = scan["e_folds_total"] < INTERRUPTED_BELOW
    # The line breaks where the surveys are interrupted.
    f_nl.loglog(scan["chi0"], np.where(interrupted, np.nan, np.abs(scan["f_nl"])), ".-",
                color="tab:blue", markersize=3, linewidth=0.8, label="crossing-time $\\delta N$")
    f_nl.loglog(scan["chi0"][interrupted], np.abs(scan["f_nl"][interrupted]), "o",
                color="0.6", markersize=3, label=r"$\epsilon_H=1$ during the turn")
    for index, (color, _) in STYLES.items():
        chi0 = trajectories["chi0"][trajectories["trajectory"] == index][0]
        row = scan[np.isclose(scan["chi0"], chi0)]
        f_nl.loglog(row["chi0"], np.abs(row["f_nl"]), "s", color=color, markersize=7,
                    label=f"trajectory {index}")
    chi0_2 = trajectories["chi0"][trajectories["trajectory"] == 2][0]
    f_nl.loglog([chi0_2], [KMS_F_NL], "*", color="tab:red", markersize=13,
                markerfacecolor="none", label=r"KMS: $|f_{\mathrm{NL}}|=43.3$")
    f_nl.set_xlabel(r"$\chi(\tau_0)/M_{\mathrm{pl}}$")
    f_nl.set_ylabel(r"$|f_{\mathrm{NL}}|$ (equilateral, $N_*=60$)")
    f_nl.set_title(r"$|f_{\mathrm{NL}}|$ against the initial $\chi$ (KMS Fig. 11)")
    f_nl.legend(fontsize=8)

    title = (r"KMS two-field model, $\xi_\phi=10$, $\xi_\chi=10.02$, "
             r"from rest at $\phi=3.1\,M_{\mathrm{pl}}$")
    if revision is not None:
        title = f"{title} | revision {revision}"
    figure.suptitle(title)
    save_figure(figure, path, "kms_trajectories", output_dir=args.output)


if __name__ == "__main__":
    main()
