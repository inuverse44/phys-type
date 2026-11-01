# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot final two-field curvature and isocurvature power spectra."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

REQUIRED_COLUMNS = {"k_over_pivot", "P_R", "P_S", "r_C"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="two_field_spectra.csv")
    parser.add_argument(
        "--output",
        type=Path,
        help="output directory (default: beside the input CSV)",
    )
    args = parser.parse_args()

    spectra = load_csv(args.csv, REQUIRED_COLUMNS)
    require_finite(spectra, REQUIRED_COLUMNS)
    if np.any(spectra["k_over_pivot"] <= 0.0):
        raise SystemExit("k_over_pivot must be positive")
    if np.any(spectra["P_R"] <= 0.0) or np.any(spectra["P_S"] <= 0.0):
        raise SystemExit("P_R and P_S must be positive")

    order = np.argsort(spectra["k_over_pivot"])
    spectra = spectra[order]
    k = spectra["k_over_pivot"]

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, (curvature, entropy) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.4),
        sharex=True,
        constrained_layout=True,
    )

    curvature.loglog(
        k, spectra["P_R"], linewidth=2.5, label=r"$\mathcal{P}_{\mathcal{R}}$"
    )
    metadata = load_metadata(args.csv)
    target = float(metadata.get("target_scalar_amplitude", "nan"))
    if np.isfinite(target) and target > 0.0:
        curvature.axhline(
            target,
            color="tab:gray",
            linestyle=":",
            linewidth=1.5,
            label=r"$A_s$ target",
        )
    entropy.loglog(
        k,
        spectra["P_S"],
        "--",
        color="tab:orange",
        linewidth=2.0,
        label=r"$\mathcal{P}_{\mathcal{S}}$",
    )
    for axes in (curvature, entropy):
        axes.axvline(1.0, color="0.75", linewidth=1.0, zorder=0)
        axes.set_ylabel("dimensionless power")
        axes.legend(loc="best")
    entropy.set_xlabel(r"$k/k_{\mathrm{pivot}}$")
    pivot_index = np.argmin(np.abs(k - 1.0))
    curvature.text(
        0.03,
        0.08,
        rf"$r_{{\mathcal{{C}}}}(k_*)={spectra['r_C'][pivot_index]:+.3f}$",
        transform=curvature.transAxes,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )

    title = "Final curvature and isocurvature power spectra"
    output_before_end = float(metadata.get("output_e_folds_before_end", "nan"))
    if np.isfinite(output_before_end):
        when = (
            "at the end of inflation"
            if output_before_end == 0.0
            else rf"{output_before_end:g} e-folds before inflation ends"
        )
        revision = metadata.get("revision")
        caption = when if revision is None else f"{when} | revision {revision}"
        title = f"{title}\n{caption}"
    figure.suptitle(title)

    save_figure(
        figure,
        args.csv,
        "two_field_final_spectra",
        output_dir=args.output,
    )


if __name__ == "__main__":
    main()
