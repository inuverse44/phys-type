# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Plot the calibrated pivot mode through the double-inflation turn."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from _common import load_csv, load_metadata, require_finite, save_figure

MODE_COLUMNS = {"k_over_pivot", "N", "k_over_aH", "P_R", "P_S", "eta_perp"}
SPECTRUM_COLUMNS = {"k_over_pivot", "P_R"}


def required_float(metadata: dict[str, str], key: str) -> float:
    try:
        value = float(metadata[key])
    except (KeyError, ValueError) as error:
        raise SystemExit(f"metadata must contain a numeric {key}") from error
    if not np.isfinite(value):
        raise SystemExit(f"metadata value is not finite: {key}")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory",
        type=Path,
        help="directory holding two_field_modes.csv and two_field_spectra.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output directory (default: the input directory)",
    )
    args = parser.parse_args()

    mode_path = args.directory / "two_field_modes.csv"
    spectrum_path = args.directory / "two_field_spectra.csv"
    modes = load_csv(mode_path, MODE_COLUMNS)
    require_finite(modes, MODE_COLUMNS)
    spectra = load_csv(spectrum_path, SPECTRUM_COLUMNS)
    require_finite(spectra, SPECTRUM_COLUMNS)
    metadata = load_metadata(spectrum_path)

    wavenumbers = np.unique(modes["k_over_pivot"])
    pivot_wavenumber = wavenumbers[np.argmin(np.abs(np.log(wavenumbers)))]
    if not np.isclose(pivot_wavenumber, 1.0, rtol=1.0e-8):
        raise SystemExit("mode CSV does not contain the pivot mode k/k_pivot = 1")
    pivot = modes[modes["k_over_pivot"] == pivot_wavenumber]
    pivot = pivot[np.argsort(pivot["N"])]

    pivot_n = required_float(metadata, "e_folds_at_pivot")
    hubble = required_float(metadata, "hubble_at_pivot_over_m_pl")
    epsilon_h = required_float(metadata, "epsilon_h_at_pivot")
    horizon_power = hubble**2 / (8.0 * np.pi**2 * epsilon_h)
    if horizon_power <= 0.0:
        raise SystemExit("single-field horizon normalization must be positive")

    relative_n = pivot["N"] - pivot_n
    turn_index = np.argmax(np.abs(pivot["eta_perp"]))
    turn_n = relative_n[turn_index]
    window = (relative_n >= -5.0) & (relative_n <= 25.0)
    final_pivot_index = np.argmin(np.abs(np.log(spectra["k_over_pivot"])))
    final_pivot_power = spectra["P_R"][final_pivot_index]

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, power = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    power.semilogy(
        relative_n[window],
        pivot["P_R"][window] / horizon_power,
        linewidth=2.4,
        label=r"$\mathcal{P}_{\mathcal{R}}/\mathcal{P}_*$",
    )
    power.semilogy(
        relative_n[window],
        pivot["P_S"][window] / horizon_power,
        "--",
        linewidth=2.2,
        label=r"$\mathcal{P}_{\mathcal{S}}/\mathcal{P}_*$",
    )
    power.axvline(0.0, color="0.5", linestyle=":", linewidth=1.2)
    power.axvline(turn_n, color="tab:red", linestyle=":", linewidth=1.2)
    power.set_xlabel(r"$N-N_*$")
    power.set_ylabel(r"power normalized by $\mathcal{P}_*=H_*^2/(8\pi^2\epsilon_*)$")
    power.set_title("Pivot-mode evolution through the double-inflation turn")

    turn = power.twinx()
    turn.plot(
        relative_n[window],
        pivot["eta_perp"][window],
        color="tab:red",
        alpha=0.55,
        label=r"$\eta_\perp$",
    )
    turn.set_ylabel(r"$\eta_\perp$")
    turn.grid(False)

    handles, labels = power.get_legend_handles_labels()
    turn_handles, turn_labels = turn.get_legend_handles_labels()
    power.legend(handles + turn_handles, labels + turn_labels, loc="best")
    power.text(
        0.03,
        0.05,
        rf"turn at $N-N_*={turn_n:.1f}$"
        + "\n"
        + rf"final $\mathcal{{P}}_\mathcal{{R}}(k_*)={final_pivot_power:.2e}$",
        transform=power.transAxes,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )

    save_figure(
        figure,
        mode_path,
        "two_field_pivot_evolution",
        output_dir=args.output,
    )


if __name__ == "__main__":
    main()
