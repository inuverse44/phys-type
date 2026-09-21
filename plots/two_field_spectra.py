# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot the two-field adiabatic and entropy spectra as one four-panel figure.

Reads the three CSV files written by the `two_field_spectra` example from the
directory given as the only argument.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

BACKGROUND = "two_field_background.csv"
SPECTRA = "two_field_spectra.csv"
MODES = "two_field_modes.csv"

BACKGROUND_COLUMNS = {
    "N",
    "phi_1",
    "phi_2",
    "e_sigma_1",
    "e_sigma_2",
    "e_s_1",
    "e_s_2",
    "eta_perp",
}
SPECTRA_COLUMNS = {
    "k_over_pivot",
    "P_11",
    "P_22",
    "P_Qsigma",
    "P_Qs",
    "P_R",
    "P_S",
    "C_RS",
    "r_C",
}
MODE_COLUMNS = {"k_over_pivot", "N", "k_over_aH", "P_R", "P_S", "eta_perp"}

ARROWS = 12


def trajectory_panel(axes: plt.Axes, background: np.ndarray) -> None:
    """The field-space path with its adiabatic and entropy frame."""
    phi_1 = background["phi_1"]
    phi_2 = background["phi_2"]
    axes.plot(phi_1, phi_2, color="tab:blue", zorder=2)
    axes.plot(phi_1[0], phi_2[0], "o", color="tab:blue", zorder=3, label="start")

    # Arrows at sparse points, so the frame stays readable.
    step = max(len(phi_1) // ARROWS, 1)
    sparse = slice(None, None, step)
    for components, color, label in (
        (("e_sigma_1", "e_sigma_2"), "tab:red", r"$e_\sigma$"),
        (("e_s_1", "e_s_2"), "tab:green", r"$e_s$"),
    ):
        axes.quiver(
            phi_1[sparse],
            phi_2[sparse],
            background[components[0]][sparse],
            background[components[1]][sparse],
            color=color,
            label=label,
            angles="xy",
            width=0.005,
            zorder=4,
        )

    axes.set_xlabel(r"$\phi^1/M_{\mathrm{pl}}$")
    axes.set_ylabel(r"$\phi^2/M_{\mathrm{pl}}$")
    axes.set_title("Background trajectory and frame")
    axes.legend(loc="best", fontsize="small")


def field_spectra_panel(axes: plt.Axes, spectra: np.ndarray) -> None:
    """The field-space and projected field spectra."""
    k = spectra["k_over_pivot"]
    for column, style, label in (
        ("P_11", "-", r"$\mathcal{P}^{11}$"),
        ("P_22", "-", r"$\mathcal{P}^{22}$"),
        ("P_Qsigma", "--", r"$\mathcal{P}_{Q_\sigma}$"),
        ("P_Qs", "--", r"$\mathcal{P}_{Q_s}$"),
    ):
        axes.plot(k, spectra[column], style, label=label)

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlabel(r"$k/k_{\mathrm{pivot}}$")
    axes.set_ylabel(r"$\mathcal{P}/M_{\mathrm{pl}}^2$")
    axes.set_title("Field spectra")
    axes.legend(loc="best", fontsize="small", ncols=2)


def curvature_panel(axes: plt.Axes, spectra: np.ndarray) -> None:
    """Curvature, entropy, and the signed cross spectrum."""
    k = spectra["k_over_pivot"]
    cross = spectra["C_RS"]
    axes.plot(k, spectra["P_R"], color="tab:blue", label=r"$\mathcal{P}_\mathcal{R}$")
    axes.plot(k, spectra["P_S"], color="tab:orange", label=r"$\mathcal{P}_\mathcal{S}$")
    # The cross spectrum is signed, so it keeps its sign on a symmetric-log
    # axis instead of being folded through an absolute value.
    axes.plot(k, cross, color="tab:purple", label=r"$\mathcal{C}_{\mathcal{RS}}$")

    # Put the linear region below every series, so each one is displayed
    # logarithmically and the smallest does not flatten against the threshold.
    magnitudes = np.abs(np.concatenate([spectra["P_R"], spectra["P_S"], cross]))
    magnitudes = magnitudes[magnitudes > 0.0]
    threshold = magnitudes.min() / 10.0 if magnitudes.size else 1.0
    axes.set_xscale("log")
    axes.set_yscale("symlog", linthresh=threshold)
    axes.set_xlabel(r"$k/k_{\mathrm{pivot}}$")
    axes.set_ylabel("dimensionless spectrum")
    axes.set_title("Curvature and entropy")

    correlation = axes.twinx()
    correlation.plot(k, spectra["r_C"], ":", color="tab:gray", label=r"$r_\mathcal{C}$")
    correlation.set_ylim(-1.05, 1.05)
    correlation.set_ylabel(r"$r_\mathcal{C}$")
    correlation.grid(False)

    handles, labels = axes.get_legend_handles_labels()
    extra_handles, extra_labels = correlation.get_legend_handles_labels()
    axes.legend(
        handles + extra_handles,
        labels + extra_labels,
        loc="best",
        fontsize="small",
        ncols=2,
    )


def evolution_panel(axes: plt.Axes, modes: np.ndarray) -> None:
    """Selected modes evolving through horizon exit, against the turn rate."""
    wavenumbers = np.unique(modes["k_over_pivot"])
    selected = wavenumbers[:: max(len(wavenumbers) // 3, 1)]
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(selected)))

    for wavenumber, color in zip(selected, colors):
        rows = modes[modes["k_over_pivot"] == wavenumber]
        order = np.argsort(rows["N"])
        rows = rows[order]
        axes.plot(
            rows["N"],
            rows["P_R"],
            color=color,
            label=rf"$k/k_*={wavenumber:.2g}$",
        )
        axes.plot(rows["N"], rows["P_S"], "--", color=color)

        # Mark where each mode leaves the horizon.
        exit_index = np.argmin(np.abs(rows["k_over_aH"] - 1.0))
        axes.plot(rows["N"][exit_index], rows["P_R"][exit_index], "o", color=color, ms=4)

    axes.set_yscale("log")
    axes.set_xlabel(r"$N$")
    axes.set_ylabel(r"$\mathcal{P}_\mathcal{R}$ (solid), $\mathcal{P}_\mathcal{S}$ (dashed)")
    axes.set_title("Mode evolution and turn rate")

    turn = axes.twinx()
    first = modes[modes["k_over_pivot"] == wavenumbers[0]]
    turn.plot(
        first["N"][np.argsort(first["N"])],
        first["eta_perp"][np.argsort(first["N"])],
        ":",
        color="tab:red",
        label=r"$\eta_\perp$",
    )
    turn.set_ylabel(r"$\eta_\perp$")
    turn.grid(False)

    handles, labels = axes.get_legend_handles_labels()
    extra_handles, extra_labels = turn.get_legend_handles_labels()
    axes.legend(
        handles + extra_handles,
        labels + extra_labels,
        loc="best",
        fontsize="small",
    )


def provenance(metadata: dict[str, str]) -> str:
    """A one-line caption identifying the run that produced the rows."""
    parts = []
    if "revision" in metadata:
        parts.append(f"revision {metadata['revision']}")
    if "e_folds_before_exit" in metadata:
        start = float(metadata["e_folds_before_exit"])
        parts.append(rf"start ${start:.1f}$ e-folds before exit")
    elif "x_init" in metadata:
        parts.append(rf"$x_{{\mathrm{{init}}}}={metadata['x_init']}$")
    if "e_folds_at_pivot" in metadata and "e_folds_at_end" in metadata:
        pivot = float(metadata["e_folds_at_pivot"])
        end = float(metadata["e_folds_at_end"])
        parts.append(rf"pivot at $N={pivot:.1f}$ of {end:.1f}")

    return " | ".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory",
        type=Path,
        help="directory holding the CSV files from the two_field_spectra example",
    )
    args = parser.parse_args()

    background = load_csv(args.directory / BACKGROUND, BACKGROUND_COLUMNS)
    spectra = load_csv(args.directory / SPECTRA, SPECTRA_COLUMNS)
    modes = load_csv(args.directory / MODES, MODE_COLUMNS)
    require_finite(background, BACKGROUND_COLUMNS)
    require_finite(spectra, SPECTRA_COLUMNS)
    require_finite(modes, MODE_COLUMNS)

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 9.0), constrained_layout=True)

    trajectory_panel(axes[0, 0], background)
    field_spectra_panel(axes[0, 1], spectra)
    curvature_panel(axes[1, 0], spectra)
    evolution_panel(axes[1, 1], modes)

    caption = provenance(load_metadata(args.directory / SPECTRA))
    if caption:
        figure.suptitle(caption, fontsize="small", y=0.005, va="bottom")

    save_figure(figure, args.directory / SPECTRA, "two_field_spectra")


if __name__ == "__main__":
    main()
