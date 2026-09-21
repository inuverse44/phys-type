# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot the single-field inflationary background from trajectory.csv."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from _common import load_csv, require_finite, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="trajectory.csv from the diagnostics example")
    args = parser.parse_args()

    columns = {"e_folds", "phi", "dphi", "H", "epsilon_h"}
    data = load_csv(args.csv, columns)
    require_finite(data, columns)
    e_folds = data["e_folds"]

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True, constrained_layout=True)

    axes[0, 0].plot(e_folds, data["phi"], color="tab:blue")
    axes[0, 0].set_ylabel(r"$\phi/M_{\mathrm{pl}}$")
    axes[0, 0].set_title("Scalar field")

    axes[0, 1].plot(e_folds, data["dphi"], color="tab:orange")
    axes[0, 1].set_ylabel(r"$\dot{\phi}/M_{\mathrm{pl}}^2$")
    axes[0, 1].set_title("Field velocity")

    axes[1, 0].plot(e_folds, data["H"], color="tab:green")
    axes[1, 0].set_ylabel(r"$H/M_{\mathrm{pl}}$")
    axes[1, 0].set_title("Hubble rate")

    axes[1, 1].plot(e_folds, data["epsilon_h"], color="tab:red")
    axes[1, 1].axhline(1.0, color="black", linestyle="--", linewidth=1.0, label=r"$\epsilon_H=1$")
    axes[1, 1].set_ylabel(r"$\epsilon_H$")
    axes[1, 1].set_title("End-of-inflation diagnostic")
    axes[1, 1].legend(frameon=False)

    for axis in axes[1, :]:
        axis.set_xlabel(r"e-folding number $N$")

    figure.suptitle("Quadratic inflation background")
    save_figure(figure, args.csv, "background")


if __name__ == "__main__":
    main()
