# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot convergence onto the slow-roll attractor in field phase space."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, require_finite, save_figure


def ordered_labels(labels: np.ndarray) -> list[str]:
    return list(dict.fromkeys(str(label) for label in labels))


def add_direction_arrow(axis: plt.Axes, x: np.ndarray, y: np.ndarray, color: str) -> None:
    if x.size < 3:
        return
    index = min(max(x.size // 3, 1), x.size - 2)
    axis.annotate(
        "",
        xy=(x[index + 1], y[index + 1]),
        xytext=(x[index], y[index]),
        arrowprops={"arrowstyle": "->", "color": color, "linewidth": 1.2},
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="attractor.csv from the diagnostics example")
    args = parser.parse_args()

    numeric = {"phi", "dphi", "slow_roll_dphi"}
    data = load_csv(args.csv, numeric | {"label"})
    require_finite(data, numeric)

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axis = plt.subplots(figsize=(8.2, 6.2), constrained_layout=True)

    for label in ordered_labels(data["label"]):
        selected = data[data["label"] == label]
        (line,) = axis.plot(selected["phi"], selected["dphi"], linewidth=1.6, label=label)
        add_direction_arrow(axis, selected["phi"], selected["dphi"], line.get_color())

    slow_roll = data[data["label"] == "slow roll"]
    if slow_roll.size:
        order = np.argsort(slow_roll["phi"])
        axis.plot(
            slow_roll["phi"][order],
            slow_roll["slow_roll_dphi"][order],
            color="black",
            linestyle="--",
            linewidth=1.3,
            label="slow-roll estimate",
        )

    axis.set_xlabel(r"$\phi/M_{\mathrm{pl}}$")
    axis.set_ylabel(r"$\dot{\phi}/M_{\mathrm{pl}}^2$")
    axis.set_title("Convergence onto the slow-roll attractor")
    axis.legend(frameon=False, fontsize="small", ncol=2)
    save_figure(figure, args.csv, "attractor")


if __name__ == "__main__":
    main()
