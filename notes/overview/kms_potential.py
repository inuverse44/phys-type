# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib>=3.10,<4",
#   "numpy>=2,<3",
# ]
# ///
"""The Einstein-frame potential of the KMS model with its three trajectories.

Run from the repository root:

    uv run notes/overview/kms_potential.py

It reads ``golden/kms_trajectories.csv`` and writes ``kms_potential.{pdf,png}``
next to this script.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

XI_PHI, XI_CHI = 10.0, 10.02
LAMBDA_PHI = 1.0
LAMBDA_CHI, G = 0.5 * LAMBDA_PHI, 1.0 * LAMBDA_PHI
STYLES = {1: ("tab:orange", ":"), 2: ("tab:red", "-"), 3: ("black", "--")}
# Growth rate of chi on the ridge from the covariant tachyonic mass
# m^2 / H^2 = -0.24 of the potential at phi = 3.1 (curvature and gravity
# terms omitted): -3/2 + sqrt(9/4 + 0.24).
RIDGE_RATE = 0.078
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def potential(phi, chi):
    """V = V_J / Omega^4 in units of lambda_phi M_pl^4."""
    omega2 = 1.0 + XI_PHI * phi**2 + XI_CHI * chi**2
    jordan = LAMBDA_PHI * phi**4 / 4 + G * phi**2 * chi**2 / 2 + LAMBDA_CHI * chi**4 / 4
    return jordan / omega2**2


def canonical_radius(r):
    """rho(r) = int_0^r sqrt(G_rr) dr with G_rr = 1/Omega^2 + 6 xi^2 r^2/Omega^4, xi = xi_phi."""
    grid = np.linspace(0.0, max(float(np.max(r)), 1e-9), 4001)
    omega2 = 1.0 + XI_PHI * grid**2
    integrand = np.sqrt(1.0 / omega2 + 6.0 * XI_PHI**2 * grid**2 / omega2**2)
    rho = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(grid))])
    return np.interp(r, grid, rho)


def load_trajectories():
    # Two "# golden-..." lines precede the header row.
    data = np.genfromtxt(ROOT / "golden" / "kms_trajectories.csv", delimiter=",",
                         skip_header=2, names=True)
    return {i: data[data["trajectory"] == i] for i in STYLES}


def main() -> None:
    trajectories = load_trajectories()

    figure = plt.figure(figsize=(12.5, 9.6), constrained_layout=True)
    grid = figure.add_gridspec(2, 2)
    surface = figure.add_subplot(grid[0, 0], projection="3d")
    top = figure.add_subplot(grid[0, 1])
    ball = figure.add_subplot(grid[1, 0])
    growth = figure.add_subplot(grid[1, 1])

    # (a) The landscape in the field coordinates of the action.
    phi = np.linspace(-1.0, 3.3, 241)
    chi = np.linspace(-0.6, 2.4, 241)
    p, c = np.meshgrid(phi, chi)
    v = potential(p, c)
    surface.plot_surface(p, c, v, cmap="viridis", linewidth=0, antialiased=True,
                         alpha=0.7, rcount=120, ccount=120)
    for i, (color, style) in STYLES.items():
        rows = trajectories[i]
        surface.plot(rows["phi"], rows["chi"], potential(rows["phi"], rows["chi"]) + 1e-4,
                     style, color=color, linewidth=1.6, label=f"trajectory {i}")
    surface.set_xlabel(r"$\phi/M_{\rm pl}$")
    surface.set_ylabel(r"$\chi/M_{\rm pl}$")
    surface.set_zlabel(r"$V/\lambda_\phi M_{\rm pl}^4$", labelpad=8)
    surface.view_init(elev=30, azim=-125)
    surface.set_title(r"(a) Einstein-frame potential $V=\tilde V/\Omega^4$")
    surface.legend(fontsize=8, loc="upper left")

    # (b) Top view: a plateau at large field, a valley along phi = 0.
    levels = np.linspace(0.0, v.max(), 31)
    filled = top.contourf(p, c, v, levels=levels, cmap="viridis")
    figure.colorbar(filled, ax=top, label=r"$V/\lambda_\phi M_{\rm pl}^4$")
    for i, (color, style) in STYLES.items():
        rows = trajectories[i]
        top.plot(rows["phi"], rows["chi"], style, color=color, linewidth=1.4,
                 label=f"trajectory {i}")
    top.plot([3.1], [0.0], "w*", markersize=12, label=r"start at rest, $\phi=3.1$")
    top.annotate("plateau", xy=(2.6, 1.2), color="0.15", fontsize=11)
    top.annotate("valley $\\phi=0$", xy=(0.02, 0.6), xytext=(0.45, 0.25), color="white",
                 arrowprops={"arrowstyle": "->", "color": "white"})
    top.set_xlabel(r"$\phi/M_{\rm pl}$")
    top.set_ylabel(r"$\chi/M_{\rm pl}$")
    top.set_title("(b) top view: the fields roll from the star to the origin")
    top.legend(fontsize=8, loc="upper right")

    # (c) Coordinates in which field space is nearly flat: the angle and the
    # canonical radial distance. The ridge theta = 0 and the valleys
    # theta = +-90 deg are now visible as such.
    theta = np.radians(np.linspace(-100.0, 100.0, 401))
    radius = np.linspace(1e-3, 3.3, 400)
    t, r = np.meshgrid(theta, radius)
    vb = potential(r * np.cos(t), r * np.sin(t))
    rho = canonical_radius(r)
    filled = ball.contourf(np.degrees(t), rho, vb, levels=31, cmap="viridis")
    for i, (color, style) in STYLES.items():
        rows = trajectories[i]
        rr = np.hypot(rows["phi"], rows["chi"])
        tt = np.degrees(np.arctan2(rows["chi"], rows["phi"]))
        ball.plot(tt, canonical_radius(rr), style, color=color, linewidth=1.4,
                  label=f"trajectory {i}")
    ball.axvline(0.0, color="white", linewidth=0.8, linestyle="--")
    ball.text(3.0, 4.45, "ridge", color="white")
    ball.text(70.0, 4.45, "valley", color="white")
    ball.set_xlim(-100.0, 100.0)
    ball.set_xlabel(r"angle $\theta=\arctan(\chi/\phi)$ [deg]")
    ball.set_ylabel(r"canonical radial distance $\rho/M_{\rm pl}$")
    ball.set_title(r"(c) the ball's view: a flat ridge at $\theta=0$")
    ball.legend(fontsize=8, loc="lower left")

    # (d) chi grows at a constant rate per e-fold on the ridge: a tachyonic
    # mass across it.
    for i, (color, style) in STYLES.items():
        rows = trajectories[i]
        growth.semilogy(rows["e_folds"], np.abs(rows["chi"]), style, color=color,
                        linewidth=1.2, label=rf"trajectory {i}, $\chi_0={rows['chi0'][0]:.1e}$")
        n = np.linspace(0.0, 60.0, 2)
        growth.semilogy(n, rows["chi0"][0] * np.exp(RIDGE_RATE * n), color=color,
                        linewidth=0.6, alpha=0.6)
    growth.set_ylim(5e-5, 5.0)
    growth.set_xlabel("e-folds from the start, $N$")
    growth.set_ylabel(r"$|\chi|/M_{\rm pl}$")
    growth.set_title(r"(d) off the ridge: $\chi\propto e^{0.078N}$ (thin lines) until the fall")
    growth.legend(fontsize=8, loc="lower right")

    figure.suptitle(r"KMS: $\Omega^2=1+10\,\phi^2+10.02\,\chi^2$, "
                    r"$\tilde V=\lambda_\phi(\phi^4+2\phi^2\chi^2+\chi^4/2)/4$, metric formulation")
    figure.savefig(HERE / "kms_potential.pdf")
    figure.savefig(HERE / "kms_potential.png", dpi=110)


if __name__ == "__main__":
    main()
