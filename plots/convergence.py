# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot fixed-step RK4 error and its empirical convergence order."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, save_figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="convergence.csv from the diagnostics example")
    args = parser.parse_args()

    columns = {"h_dt", "abs_error", "order"}
    data = load_csv(args.csv, columns)
    finite = np.isfinite(data["h_dt"]) & np.isfinite(data["abs_error"])
    positive = finite & (data["h_dt"] > 0.0) & (data["abs_error"] > 0.0)
    valid_count = int(np.count_nonzero(positive))
    if valid_count < 3:
        raise SystemExit("convergence CSV needs at least three positive error rows")
    excluded = data.size - valid_count
    if excluded:
        print(f"excluded {excluded} non-finite or non-positive convergence row(s)")
    selected = data[positive]
    order = np.argsort(selected["h_dt"])
    h_dt = selected["h_dt"][order]
    error = selected["abs_error"][order]

    # Anchor in the coarse, convergent regime. Finer rows eventually reach the
    # floating-point/reference-error floor and should visibly depart from the
    # fourth-order guide.
    anchor = max(len(h_dt) - 2, 0)
    fourth_order = error[anchor] * (h_dt / h_dt[anchor]) ** 4

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)

    axes[0].loglog(h_dt, error, "o-", color="tab:blue", label=r"$|\phi-\phi_{\rm ref}|$")
    axes[0].loglog(h_dt, fourth_order, "--", color="black", label=r"reference $\propto (H_0\,\Delta t)^4$")
    axes[0].set_xlabel(r"dimensionless step $H_0\,\Delta t$")
    axes[0].set_ylabel("absolute field error")
    axes[0].set_title("Global error")
    axes[0].legend(frameon=False)

    finite_order = np.isfinite(data["order"])
    axes[1].semilogx(data["h_dt"][finite_order], data["order"][finite_order], "o-", color="tab:orange")
    axes[1].axhline(4.0, color="black", linestyle="--", linewidth=1.0, label="classical RK4")
    axes[1].set_xlabel(r"dimensionless step $H_0\,\Delta t$")
    axes[1].set_ylabel("empirical order")
    axes[1].set_title("Order under step halving")
    axes[1].legend(frameon=False)

    figure.suptitle("Fixed-step RK4 convergence")
    save_figure(figure, args.csv, "convergence")


if __name__ == "__main__":
    main()
