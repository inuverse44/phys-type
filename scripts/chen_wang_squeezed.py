# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy>=2,<3",
#   "scipy>=1.14,<2",
# ]
# ///
"""The squeezed coefficient s(nu) of quasi-single-field inflation.

Chen and Wang, arXiv:0911.3380, Eqs. (5.10)--(5.11): for p3 << p1 = p2,

    <. . .> -> s(nu) c2^3 c3 / (H R^6) p1^(nu - 7/2) p2^-1 p3^(-3/2 - nu),

with the transfer vertex c2 = 2 R thetadot and the cubic coupling
c3 = V''' / 6, at first order in each transfer. With x = -s, Eq. (5.11) is

    s(nu) = 3 pi^2 Gamma(nu) / 2^(3 - nu) * (T1 + T2) * Y,
    T1 = int_{0 < s1 < s2 < s3} s1^-1/2 s2^(1/2-nu) s3^-1/2 sin s1
             Im(H1(s1) H2(s2)) Im(H1(s2) H2(s3) e^{-i s3}),
    T2 = int_{0 < s1 < s2 < s3} s1^-1/2 s2^-1/2 s3^(1/2-nu) sin s1 sin s2
             Im(H2(s1) H2(s2) H1(s3)^2),
    Y  = int_0^inf ds s^-1/2 Re(H1(s) e^{i s}),

with H1, H2 the Hankel functions of order nu. The nested integrals are
evaluated as tails from infinity on a grid that is logarithmic below s = 1
and uniform above, with a smooth window
that switches the integrands off at large s, which converges for
1/2 < nu < 3/2 as the paper notes (and, with the window, below 1/2).

The grid starts at s = e^LOWER (environment variable CW_LOWER, default -40).
For nu close to 3/2 the integrands converge at s -> 0 only as s^(3/2 - nu)
and their imaginary parts cancel to roundoff below about e^-20; for nu = 1.2
the value is extrapolated geometrically in the cutoff from LOWER = -14 to -20,
which gives 84.4 to 84.5.
"""

import argparse
import os

import numpy as np
from scipy.special import gamma, hankel1, hankel2


LOWER = float(os.environ.get("CW_LOWER", "-40"))


def window(s: np.ndarray, end: float) -> np.ndarray:
    """1 below end/2, 0 above end, with every derivative continuous."""
    t = np.clip((s - 0.5 * end) / (0.5 * end), 0.0, 1.0)
    f = lambda x: np.where(x > 0, np.exp(-1.0 / np.maximum(x, 1e-300)), 0.0)
    return f(1 - t) / (f(1 - t) + f(t))


def tail(values: np.ndarray, s: np.ndarray) -> np.ndarray:
    """int_s^end of a sampled function, by the trapezoid rule from the end."""
    pieces = 0.5 * (values[1:] + values[:-1]) * np.diff(s)
    return np.concatenate((np.cumsum(pieces[::-1])[::-1], [0.0]))


def integral(values: np.ndarray, s: np.ndarray) -> complex:
    return np.sum(0.5 * (values[1:] + values[:-1]) * np.diff(s))


def s_coefficient(nu: float, end: float = 4000.0, step: float = 0.005) -> float:
    # Logarithmic below s = 1, for the integrable powers at s -> 0; uniform
    # above, for the oscillations.
    s = np.concatenate((np.exp(np.arange(LOWER, 0.0, 0.01)), np.arange(1.0, end + step / 2, step)))
    w = window(s, end)
    h1, h2 = hankel1(nu, s), hankel2(nu, s)

    # T1: the s3 tail G(s2) = int s3^-1/2 H2(s3) e^{-i s3}, then the s2 tail
    # of s2^(1/2-nu) Im(H1(s1) H2(s2)) Im(H1(s2) G(s2)).
    g = tail(s**-0.5 * h2 * np.exp(-1j * s) * w, s)
    q = (h1 * g).imag
    weight2 = s ** (0.5 - nu) * q * w
    re_part = tail(weight2 * h2.imag, s)
    im_part = tail(weight2 * h2.real, s)
    # Im(H1(s1) H2(s2)) = Re H1(s1) Im H2(s2) + Im H1(s1) Re H2(s2).
    middle = h1.real * re_part + h1.imag * im_part
    t1 = integral(s**-0.5 * np.sin(s) * middle * w, s)

    # T2: G2(s2) = int s3^(1/2-nu) H1(s3)^2, then W = H2(s2) G2(s2) and
    # Im(H2(s1) W(s2)) under int s2^-1/2 sin s2.
    g2 = tail(s ** (0.5 - nu) * h1**2 * w, s)
    big_w = h2 * g2
    weight2 = s**-0.5 * np.sin(s) * w
    re_w = tail(weight2 * big_w.real, s)
    im_w = tail(weight2 * big_w.imag, s)
    # Im(H2(s1) W) = Re H2(s1) Im W + Im H2(s1) Re W.
    middle2 = h2.real * im_w + h2.imag * re_w
    t2 = integral(s**-0.5 * np.sin(s) * middle2 * w, s)

    y = integral(s**-0.5 * (h1 * np.exp(1j * s)).real * w, s)
    return 3 * np.pi**2 * gamma(nu) / 2 ** (3 - nu) * (t1 + t2) * y


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nu", type=float, nargs="+", required=True)
    parser.add_argument("--end", type=float, default=4000.0)
    parser.add_argument("--step", type=float, default=0.005)
    args = parser.parse_args()
    for nu in args.nu:
        print(f"nu {nu}: s = {s_coefficient(nu, args.end, args.step):.8e}")


if __name__ == "__main__":
    main()
