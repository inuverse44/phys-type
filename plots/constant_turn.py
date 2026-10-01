# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
#   "scipy>=1.14,<2",
# ]
# ///
"""Compare constant-turn spectra with the exact solutions of arXiv:2606.18248.

Reads ``constant_turn_massive.csv`` and ``constant_turn_massless.csv`` written
by ``cargo run --release --example constant_turn_spectra``. The closed forms
are evaluated at the lambda and mu^2 measured at horizon exit: Eqs. (43) and
(45) for a massive entropy field, and Eq. (39) for a massless one, both as
printed, with ln(2/z_end), and with ln(1/(2 z_end)).
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import loggamma, psi

from _common import load_csv, load_metadata, require_finite, save_figure

MASSIVE_COLUMNS = {"lambda_target", "mu2_target", "lambda", "mu2", "ratio"}
MASSLESS_COLUMNS = {"lambda_target", "lambda", "epsilon", "z_end", "ratio"}
COLORS = {0.5: "tab:blue", 1.0: "tab:orange", 2.0: "tab:green"}


def massive_ratio(lam, mu2):
    """Eq. (43) for mu^2 <= 9/4 and Eq. (45) above it."""
    lam = np.asarray(lam, dtype=float)
    mu2 = np.asarray(mu2, dtype=float)
    nu = np.sqrt(np.clip(2.25 - mu2, 0.0, None))
    rho = np.sqrt(np.clip(mu2 - 2.25, 0.0, None))
    g = lambda re, im: loggamma(re + 1j * im).real
    a, b = 0.75 - 0.5 * nu, 0.75 + 0.5 * nu
    light = 2.0 * (g(a, 0.0) + g(b, 0.0) - g(a, 0.5 * lam) - g(b, 0.5 * lam))
    heavy = (
        4.0 * g(0.75, 0.5 * rho)
        - 2.0 * g(0.75, 0.5 * (lam + rho))
        - 2.0 * g(0.75, 0.5 * (lam - rho))
    )
    return np.exp(np.where(mu2 <= 2.25, light, heavy))


def massless_ratio(lam, z_end, argument):
    """Eq. (39) with the logarithm ln(argument / z_end)."""
    lam = np.asarray(lam, dtype=float)
    growth = 2.0 + np.log(argument / z_end) - 2.0 * np.euler_gamma - psi(0.5j * lam).real
    return lam * np.sinh(np.pi * lam) / (np.pi * (1.0 + lam**2)) * growth**2 + (
        np.pi * lam / (2.0 * (1.0 + lam**2)) / np.tanh(0.5 * np.pi * lam)
    )


def power_of_ten(value: float) -> str:
    """``3e-07`` as ``3\\times10^{-7}`` for mathtext."""
    mantissa, exponent = f"{value:.0e}".split("e")
    power = f"10^{{{int(exponent)}}}"
    return power if mantissa == "1" else rf"{mantissa}\times{power}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory holding both CSV files")
    parser.add_argument("--output", type=Path, help="output directory")
    args = parser.parse_args()

    massive_path = args.directory / "constant_turn_massive.csv"
    massless_path = args.directory / "constant_turn_massless.csv"
    massive = load_csv(massive_path, MASSIVE_COLUMNS)
    massless = load_csv(massless_path, MASSLESS_COLUMNS)
    require_finite(massive, MASSIVE_COLUMNS)
    require_finite(massless, MASSLESS_COLUMNS)
    epsilon = load_metadata(massive_path).get("epsilon", "?")
    revision = load_metadata(massive_path).get("revision")

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(11.0, 6.8),
        sharex="col",
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )
    (spectrum, massless_axes), (residual, massless_residual) = axes

    mu = np.linspace(0.9, 4.2, 400)
    for lam, color in COLORS.items():
        rows = massive[np.isclose(massive["lambda_target"], lam)]
        spectrum.semilogy(
            mu, massive_ratio(lam, mu**2), color=color, linewidth=1.8,
            label=rf"exact, $\lambda={lam:g}$",
        )
        spectrum.semilogy(
            np.sqrt(rows["mu2"]), rows["ratio"], "o", color=color,
            markerfacecolor="white", markersize=6,
            label="phys-type" if lam == 0.5 else None,
        )
        residual.plot(
            np.sqrt(rows["mu2"]),
            rows["ratio"] / massive_ratio(rows["lambda"], rows["mu2"]) - 1.0,
            "o-", color=color, markerfacecolor="white", markersize=5,
        )
    spectrum.axvline(1.5, color="0.6", linestyle=":", linewidth=1.0)
    spectrum.text(1.53, 300.0, r"$\mu=3H/2$", color="0.4", fontsize=9)
    spectrum.set_ylabel(r"$\Delta_\zeta/\Delta_0$")
    spectrum.set_title(
        "massive entropy field: Eqs. (43)–(45)\n"
        rf"$\epsilon={power_of_ten(float(epsilon))}$, $z_{{\mathrm{{end}}}}\to0$"
    )
    spectrum.legend(loc="upper right")
    residual.axhline(0.0, color="0.5", linewidth=1.0)
    residual.set_xlabel(r"$\mu/H$")
    residual.set_ylabel("phys-type / exact − 1")
    residual.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

    lam = np.linspace(0.05, 2.5, 300)
    z_end = float(np.median(massless["z_end"]))
    massless_axes.semilogy(
        lam, massless_ratio(lam, z_end, 0.5), color="tab:purple", linewidth=1.8,
        label=r"Eq. (39) with $\ln(1/2z_{\mathrm{end}})$",
    )
    massless_axes.semilogy(
        lam, massless_ratio(lam, z_end, 2.0), "--", color="tab:red", linewidth=1.6,
        label=r"Eq. (39) as printed, $\ln(2/z_{\mathrm{end}})$",
    )
    massless_axes.semilogy(
        massless["lambda"], massless["ratio"], "o", color="tab:purple",
        markerfacecolor="white", markersize=6, label="phys-type",
    )
    massless_axes.set_ylabel(r"$\Delta_\zeta/\Delta_0$")
    massless_axes.set_title(
        "massless entropy field: Eq. (39)\n"
        rf"$z_{{\mathrm{{end}}}}=e^{{-10}}$, $\epsilon={power_of_ten(float(massless['epsilon'][0]))}$"
    )
    massless_axes.legend(loc="upper left")
    corrected = massless_ratio(massless["lambda"], massless["z_end"], 0.5)
    printed = massless_ratio(massless["lambda"], massless["z_end"], 2.0)
    massless_residual.plot(
        massless["lambda"], massless["ratio"] / corrected - 1.0, "o-",
        color="tab:purple", markerfacecolor="white", markersize=5,
        label=r"vs $\ln(1/2z_{\mathrm{end}})$",
    )
    massless_residual.axhline(0.0, color="0.5", linewidth=1.0)
    massless_residual.set_xlabel(r"$\lambda=-2\Omega/H$")
    massless_residual.set_ylabel("phys-type / exact − 1")
    massless_residual.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    massless_residual.text(
        0.02, 0.8,
        rf"printed Eq. (39) is {100 * np.min(printed / massless['ratio'] - 1):.0f}–"
        rf"{100 * np.max(printed / massless['ratio'] - 1):.0f}% high",
        transform=massless_residual.transAxes, color="tab:red", fontsize=9,
    )

    title = "Constant-turn power spectra against the exact solutions of arXiv:2606.18248"
    if revision is not None:
        title = f"{title}\nrevision {revision}"
    figure.suptitle(title)
    save_figure(figure, massive_path, "constant_turn_spectra", output_dir=args.output)


if __name__ == "__main__":
    main()
