"""Independent check of the constant-turn power spectra of arXiv:2606.18248.

Integrates the mode equations (8)-(9) of Huenupi, Munoz, Palma, and Sypsas in
z = k/(aH), from Bunch--Davies initial conditions at z0 to z_end, and compares
sum_b |phi_b|^2 / |u_0(0)|^2 with the closed forms (39), (43), and (45).
Eq. (39) is evaluated both as printed, with ln(2/z_end), and with
ln(1/(2 z_end)), which is what Eq. (38) implies.

Usage:
    python3 -m venv venv && venv/bin/pip install numpy scipy mpmath
    venv/bin/python docs/derivations/constant_turn_exact_solutions.py
"""

import numpy as np
from mpmath import coth, digamma, euler, gamma, log, mpf, pi, re, sinh, sqrt
from scipy.integrate import solve_ivp
from scipy.special import hankel1


def numerical_ratio(lam, mu2, z0=3000.0, z_end=1.0e-5):
    """Delta_zeta / Delta_0 from Eqs. (8)-(9), for 0 <= mu^2 <= 9/4."""
    nu = np.sqrt(9.0 / 4.0 - mu2)
    theta = 0.5 * np.pi * (nu - 1.5)

    def u0(z):
        return 1j * (1 - 1j * z) * np.exp(1j * z)

    def du0(z):
        return 1j * z * np.exp(1j * z)

    def um(z):
        return -np.sqrt(np.pi / 2) * np.exp(1j * theta) * z**1.5 * hankel1(nu, z)

    def dum(z, h=1.0e-6):
        return (um(z + h) - um(z - h)) / (2 * h)

    # y = (phi, D_z phi, sigma, d sigma / dz), integrated in x = ln z.
    def rhs(x, y):
        z = np.exp(x)
        phi, d, s, ds = y[0] + 1j * y[1], y[2] + 1j * y[3], y[4] + 1j * y[5], y[6] + 1j * y[7]
        derivative = np.array(
            [
                d - lam * s / z,
                2 * d / z - phi,
                ds,
                2 * ds / z - (1 + mu2 / z**2) * s + lam * d / z,
            ]
        ) * z
        return np.concatenate([[v.real, v.imag] for v in derivative])

    total = 0.0
    for initial in [(u0(z0), du0(z0), 0.0, 0.0), (0.0, 0.0, um(z0), dum(z0))]:
        phi, dphi, s, ds = (complex(v) for v in initial)
        start = [phi, dphi + lam * s / z0, s, ds]
        y0 = np.concatenate([[v.real, v.imag] for v in start])
        solution = solve_ivp(
            rhs, [np.log(z0), np.log(z_end)], y0, method="DOP853", rtol=1e-11, atol=1e-13
        )
        y = solution.y[:, -1]
        total += abs(y[0] + 1j * y[1]) ** 2
    return total


def massive_ratio(lam, mu2):
    """Eq. (43) for mu^2 <= 9/4 and Eq. (45) above it."""
    lam, mu2 = mpf(lam), mpf(mu2)
    if mu2 <= mpf(9) / 4:
        nu = sqrt(mpf(9) / 4 - mu2)
        a, b = mpf(3) / 4 - nu / 2, mpf(3) / 4 + nu / 2
        return abs(gamma(a) * gamma(b) / (gamma(a + 0.5j * lam) * gamma(b + 0.5j * lam))) ** 2
    rho = sqrt(mu2 - mpf(9) / 4)
    return abs(gamma(mpf(3) / 4 + 0.5j * rho)) ** 4 / (
        abs(gamma(mpf(3) / 4 + 0.5j * (lam + rho))) ** 2
        * abs(gamma(mpf(3) / 4 + 0.5j * (lam - rho))) ** 2
    )


def massless_ratio(lam, z_end, argument):
    """Eq. (39) with the logarithm ln(argument / z_end)."""
    lam = mpf(lam)
    growth = 2 + log(argument / mpf(z_end)) - 2 * euler - re(digamma(0.5j * lam))
    return lam * sinh(pi * lam) / (pi * (1 + lam**2)) * growth**2 + pi * lam / (
        2 * (1 + lam**2)
    ) * coth(pi * lam / 2)


def main():
    print("massive: lambda, mu^2/H^2, numerical, closed form, relative difference")
    for lam, mu2 in [(1.0, 2.0), (2.0, 2.0), (1.0, 1.25)]:
        z_end = 1.0e-9 if mu2 < 2 else 1.0e-5
        numerical = numerical_ratio(lam, mu2, z_end=z_end)
        exact = float(massive_ratio(lam, mu2))
        print(f"  {lam} {mu2} {numerical:.6f} {exact:.6f} {numerical / exact - 1:.1e}")

    print("massless: lambda, z_end, numerical, (39) printed, (39) with ln(1/(2 z_end))")
    for lam in [0.02, 0.5, 1.0, 2.0, 3.0]:
        for z_end in [1.0e-3, 1.0e-5]:
            numerical = numerical_ratio(lam, 0.0, z_end=z_end)
            printed = float(massless_ratio(lam, z_end, 2))
            corrected = float(massless_ratio(lam, z_end, mpf(1) / 2))
            print(
                f"  {lam} {z_end:.0e} {numerical:.6g} {printed / numerical - 1:+.3f}"
                f" {corrected / numerical - 1:+.1e}"
            )


if __name__ == "__main__":
    main()
