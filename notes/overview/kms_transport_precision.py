# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""Where the transport f_NL of KMS trajectory 2 loses precision.

Run from the repository root:

    uv run notes/overview/kms_transport_precision.py

It reads ``kms_transport_precision.csv`` next to this script, written by

    cargo run --release --example kms_transport_probe -- transport 1 \
        1e-10,1e-11,1e-12,1e-13 <slices>

(the slices are 20, 16, 14, then 12 to 0.1 in steps of 0.1, then 0.05 and
0.01 e-folds before the end), and writes ``kms_transport_precision.{pdf,png}``.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
# The first maximum of the turn rate, 5.15 e-folds before the end
# (docs/general-nonminimal-inflation.md, Results).
TURN = 5.15
# Phase-space delta N from the slice 50 e-folds before the end, with the chi
# step extrapolated to zero (issue #67).
DELTA_N = 38.30
COLORS = {1e-10: "tab:blue", 1e-11: "tab:green", 1e-12: "tab:orange", 1e-13: "tab:red"}


def main():
    data = np.genfromtxt(HERE / "kms_transport_precision.csv", delimiter=",", names=True)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)

    finest = data[data["tol"] == 1e-13]
    final = finest["p_zeta"][np.argmin(finest["before_end"])]
    top.semilogy(finest["before_end"], finest["p_zeta"] / final, color="black")
    top.set_ylabel(r"$P_\zeta(N)\,/\,P_\zeta(\mathrm{end})$")

    for tol, color in COLORS.items():
        rows = data[data["tol"] == tol]
        exponent = int(round(np.log10(tol)))
        bottom.semilogy(
            rows["before_end"],
            np.abs(rows["f_nl"]),
            color=color,
            label=rf"transport, $\varepsilon=10^{{{exponent}}}$",
        )
    bottom.axhline(DELTA_N, color="gray", ls="--", label=r"phase-space $\delta N$ (38.30)")
    bottom.set_ylabel(r"$|f_{\mathrm{NL}}|$ on the slice")
    bottom.set_xlabel("e-folds before the end of inflation")
    bottom.legend(fontsize=8, loc="upper left")

    for axis in (top, bottom):
        axis.axvline(TURN, color="gray", lw=0.8)
        axis.axvspan(3.0, 0.0, color="tab:red", alpha=0.08)
        axis.set_xlim(20.0, 0.0)
    top.text(TURN, 3e0, " turn", color="gray", fontsize=8)
    top.text(1.5, 1e4, "valley", color="tab:red", fontsize=8, ha="center")

    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(HERE / f"kms_transport_precision.{suffix}", dpi=200)


if __name__ == "__main__":
    main()
