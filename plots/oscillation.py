# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot coherent post-inflation oscillations and their averaged expansion."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, require_finite, save_figure


def moving_average(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if x.size < 3:
        return x, y

    spacing = float(np.median(np.diff(x)))
    window = max(3, min(x.size, round((2.0 * np.pi) / spacing)))
    kernel = np.ones(window) / window
    averaged = np.convolve(y, kernel, mode="valid")
    offset = (window - 1) // 2
    return x[offset : offset + averaged.size], averaged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="oscillation.csv from the diagnostics example")
    args = parser.parse_args()

    columns = {
        "m_t_since_end",
        "phi",
        "dphi",
        "a_over_a_end",
        "w",
        "matter_a_over_a_end",
    }
    data = load_csv(args.csv, columns)
    require_finite(data, columns)
    mt = data["m_t_since_end"]
    average_mt, average_w = moving_average(mt, data["w"])

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)

    axes[0, 0].plot(data["phi"], data["dphi"], color="tab:purple", linewidth=1.0)
    axes[0, 0].scatter(data["phi"][0], data["dphi"][0], color="black", s=24, zorder=3, label="inflation end")
    axes[0, 0].set_xlabel(r"$\phi/M_{\mathrm{pl}}$")
    axes[0, 0].set_ylabel(r"$\dot{\phi}/M_{\mathrm{pl}}^2$")
    axes[0, 0].set_title("Post-inflation phase space")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].plot(mt, data["phi"], color="tab:blue", linewidth=1.0)
    axes[0, 1].set_xlabel(r"$m(t-t_{\mathrm{end}})$")
    axes[0, 1].set_ylabel(r"$\phi/M_{\mathrm{pl}}$")
    axes[0, 1].set_title("Coherent field oscillations")

    axes[1, 0].plot(mt, data["w"], color="tab:gray", alpha=0.45, linewidth=0.8, label=r"$w$")
    axes[1, 0].plot(average_mt, average_w, color="tab:red", linewidth=1.8, label="one-period average")
    axes[1, 0].axhline(0.0, color="black", linestyle="--", linewidth=1.0)
    axes[1, 0].set_xlabel(r"$m(t-t_{\mathrm{end}})$")
    axes[1, 0].set_ylabel(r"$w$")
    axes[1, 0].set_ylim(-1.1, 1.1)
    axes[1, 0].set_title("Matter-like averaged equation of state")
    axes[1, 0].legend(frameon=False)

    axes[1, 1].plot(mt, data["a_over_a_end"], color="tab:green", linewidth=1.8, label="scalar-field solution")
    axes[1, 1].plot(
        mt,
        data["matter_a_over_a_end"],
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="matter prediction",
    )
    axes[1, 1].set_xlabel(r"$m(t-t_{\mathrm{end}})$")
    axes[1, 1].set_ylabel(r"$a/a_{\mathrm{end}}$")
    axes[1, 1].set_title("Scale-factor evolution")
    axes[1, 1].legend(frameon=False)

    figure.suptitle("Quadratic scalar after inflation")
    save_figure(figure, args.csv, "oscillation")


if __name__ == "__main__":
    main()
