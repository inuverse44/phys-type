# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Compare squeezed constant-turn bispectra with the scaling of arXiv:2607.14529.

Reads ``constant_turn_squeezed.csv`` written by
``cargo run --release --example constant_turn_squeezed``. For k1 = k2 and
k3 = r k1, the leading nonanalytic squeezed shape is S ~ r^(1/2 - nu), so the
reduced f_NL ~ r^(3/2 - nu). The paper's claim is that nu is set by the entropy
mass mu, nu = sqrt(9/4 - mu^2/H^2), not by the effective mass
m_eff^2 = mu^2 - lambda^2 H^2. Both power laws are drawn through the
squeezed-most point; the lower panel shows the local logarithmic slope.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

COLUMNS = {"lambda_target", "mu2_target", "lambda", "mu2", "r", "f_nl"}
COLORS = ("tab:blue", "tab:orange", "tab:green")


def exponent(mu2: float) -> float:
    """3/2 - nu for nu = sqrt(9/4 - mu^2), continued as 3/2 for mu^2 > 9/4."""
    return 1.5 - np.sqrt(max(2.25 - mu2, 0.0))


def power_of_ten(value: float) -> str:
    """``1e-05`` as ``10^{-5}`` for mathtext."""
    mantissa, power = f"{value:.0e}".split("e")
    power = f"10^{{{int(power)}}}"
    return power if mantissa == "1" else rf"{mantissa}\times{power}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding the CSV")
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()

    path = args.directory / "constant_turn_squeezed.csv"
    data = load_csv(path, COLUMNS)
    require_finite(data, COLUMNS)
    metadata = load_metadata(path)

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, (scaling, slope_axes) = plt.subplots(
        2, 1, figsize=(7.6, 8.2), sharex=True,
        gridspec_kw={"height_ratios": [3, 2]}, constrained_layout=True,
    )
    cases = sorted({(float(l), float(m)) for l, m in zip(data["lambda_target"], data["mu2_target"])})
    for color, (lam, mu2) in zip(COLORS, cases):
        rows = data[np.isclose(data["lambda_target"], lam) & np.isclose(data["mu2_target"], mu2)]
        rows = rows[np.argsort(rows["r"])]
        r, f = rows["r"], np.abs(rows["f_nl"])
        predicted = exponent(float(rows["mu2"][0]))
        effective = exponent(float(rows["mu2"][0]) - float(rows["lambda"][0]) ** 2)
        nu = np.sqrt(2.25 - mu2)
        scaling.loglog(
            r, f, "o", color=color, markerfacecolor="white", markersize=6,
            label=rf"$\lambda={lam:g}$, $\mu^2/H^2={mu2:g}$ ($\nu={nu:.2g}$)",
        )
        grid = np.geomspace(r[0] * 0.8, 1.2, 50)
        scaling.loglog(grid, f[0] * (grid / r[0]) ** predicted, "-", color=color, linewidth=1.5)
        scaling.loglog(grid, f[0] * (grid / r[0]) ** effective, ":", color=color, linewidth=1.5)
        middle = np.sqrt(r[1:] * r[:-1])
        local = np.diff(np.log(f)) / np.diff(np.log(r))
        slope_axes.semilogx(middle, local, "o-", color=color, markerfacecolor="white", markersize=5)
        slope_axes.axhline(predicted, color=color, linewidth=1.2)
        slope_axes.axhline(effective, color=color, linewidth=1.2, linestyle=":")

    scaling.set_ylabel(r"$|f_{\mathrm{NL}}(k_S,k_S,rk_S)|$")
    handles, labels = scaling.get_legend_handles_labels()
    handles += [
        plt.Line2D([], [], color="0.3", linestyle="-", label=r"$r^{3/2-\nu}$, $\nu$ from $\mu$"),
        plt.Line2D([], [], color="0.3", linestyle=":", label=r"$r^{3/2-\nu}$, $\nu$ from $m_{\mathrm{eff}}$"),
    ]
    scaling.legend(handles=handles, fontsize=8, loc="lower left")
    scaling.set_ylim(5.0e-4, 2.0)
    slope_axes.set_xlabel(r"$r=k_L/k_S$")
    slope_axes.set_ylabel(r"$d\ln|f_{\mathrm{NL}}|/d\ln r$")
    slope_axes.set_ylim(-1.0, 1.2)

    epsilon = power_of_ten(float(metadata.get("epsilon", "nan")))
    title = (
        "Squeezed bispectrum on a constant-turn orbit (arXiv:2607.14529)\n"
        rf"tree-level transport, $\epsilon={epsilon}$"
    )
    revision = metadata.get("revision")
    if revision is not None:
        title = f"{title} | revision {revision}"
    figure.suptitle(title)
    save_figure(figure, path, "constant_turn_squeezed", output_dir=args.output)


if __name__ == "__main__":
    main()
