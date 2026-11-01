# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "mpmath>=1.3,<2",
#   "numpy>=2,<3",
# ]
# ///
"""Squeezed-limit coefficients of arXiv:2607.14529, Sec. X.

Evaluates the closed forms of the leading nonanalytic squeezed bispectra of
the alpha, beta, and gamma interactions, Eqs. (236)--(242), (245)--(250), and
(252)--(257), for given lambda and nu. A curvature leg of the double
integrals J diverges at t -> 0 as t^(-1 + i w lambda / 2) and is defined by
analytic continuation in a = i w lambda / 2; with t = e^u,

    int_0^inf dt t^(-1 + a) g(t)
      = int du e^(a u) [g(e^u) - g(0) / (1 + e^u)] + g(0) pi / sin(pi a),

applied in both variables. The subtracted integrand is smooth and decays at
both ends, so the trapezoid rule in u converges exponentially. An entropy leg, (t / (1 + t))^a = t^(-1 + (1 + a)) (1 + t)^(-a), converges at
t = 0 and is integrated without subtraction. Beyond
nu < 1 (beta) and nu < 1/2 (gamma) the integrals also diverge at large t and
would need continuation in nu, which is not implemented.

Prints c_+- and the prediction for r^2 B(1, 1, r) / g, g the coupling, in the
units of phys_type_models::in_in (a cubic coefficient carries
c = sqrt(2 epsilon) m_Pl, H = 1). Reading Eqs. (240), (248), and (255) as
f_NL S, Eq. (197) and B_zeta = B / c^4 give

    r^2 B / g = s (28.8 / (2^n sqrt 2)) [c_+ r^(1/2 - nu) + c_- r^(1/2 + nu)],

with (s, n) = (-1, 10), (-1, 8), and (+1, 8). The printed prediction also
carries the factor 10/9 by which those equations fall short of their own
Eqs. (197) and (199); see docs/literature-discrepancies.md.
"""

import argparse
import itertools

import mpmath as mp
import numpy as np


def hyp(nu: float, lam: float, w: int, t: np.ndarray) -> np.ndarray:
    """D_nu^(w)(-t) = 2F1(1/2 - nu, 1/2 + nu; 1 + i w lambda; -t)."""
    c = 1 + 1j * w * lam
    return np.array([complex(mp.hyp2f1(0.5 - nu, 0.5 + nu, c, -x)) for x in t])


def c_coefficient(nu: float, lam: float, w: int, sign: int) -> complex:
    """C_+-^(w), Eq. (119)."""
    s = sign
    return complex(
        mp.gamma(1 + 1j * w * lam)
        * mp.gamma(2 * s * nu)
        / (mp.gamma(0.5 + s * nu) * mp.gamma(0.5 + s * nu + 1j * w * lam))
    )


def f0(nu: float, lam: float, w: int) -> complex:
    """F_0^(w), Eq. (121)."""
    a = 1j * w * lam / 2
    return complex(
        1j * mp.exp(mp.pi * w * lam / 4) / mp.sqrt(2 * mp.pi)
        * mp.gamma(0.5 + a) * mp.gamma(0.75 - nu / 2) * mp.gamma(0.75 + nu / 2)
        / (mp.gamma(0.75 - nu / 2 + a) * mp.gamma(0.75 + nu / 2 + a))
    )


class Grid:
    """Uniform nodes u_j, t_j = e^(u_j); index 0 of a value array holds t = 0."""

    def __init__(self, lower: float, upper: float, step: float):
        self.u = np.arange(lower, upper + step / 2, step)
        self.step = step
        self.t = np.concatenate(([0.0], np.exp(self.u)))

    def regularized(self, a: complex, values: np.ndarray, axis: int) -> np.ndarray:
        """The continued int dt t^(-1+a) g(t) along `axis`, with g(0) at index 0."""
        values = np.moveaxis(values, axis, 0)
        at_zero, g = values[0], values[1:]
        shape = (-1,) + (1,) * (g.ndim - 1)
        kernel = (np.exp(a * self.u) * self.step).reshape(shape)
        damping = (1.0 / (1.0 + np.exp(self.u))).reshape(shape)
        return (kernel * (g - at_zero * damping)).sum(axis=0) + at_zero * complex(
            mp.pi / mp.sin(mp.pi * a)
        )

    def plain(self, b: complex, values: np.ndarray, axis: int) -> np.ndarray:
        """int dt t^(-1+b) g(t) for Re b = 1, which converges at t = 0: no
        subtraction, whose 1/(1+t) tail would not decay for Re b = 1."""
        values = np.moveaxis(values, axis, 0)[1:]
        shape = (-1,) + (1,) * (values.ndim - 1)
        kernel = (np.exp(b * self.u) * self.step).reshape(shape)
        return (kernel * values).sum(axis=0)


# Per interaction: which legs are curvature legs, the sign s and power n of
# the prefactor, and N^(w) of Eqs. (242), (250), (257) without the Gammas.
UPPER = 60.0

INTERACTIONS = {
    "alpha": ((True, True, False), -1, 10, lambda w1, w2, w3: -1j * w3),
    "beta": ((True, False, False), -1, 8, lambda w1, w2, w3: w2 * w3),
    "gamma": ((False, False, False), 1, 8, lambda w1, w2, w3: -1j * w1 * w2 * w3),
}


def j_integral(nu, lam, w1, w2, sign, grid, d, curvature):
    """J_{+-}^(w1, w2), Eqs. (238), (247), (254)."""
    exponents, factors = [], []
    for w, is_curvature, axis in ((w1, curvature[0], 0), (w2, curvature[1], 1)):
        a = 1j * w * lam / 2
        t = grid.t if axis == 0 else grid.t
        if is_curvature:
            # t^-2 (t / (1 + t))^(1 + a) = t^(-1 + a) (1 + t)^(-1 - a).
            exponents.append(a)
            factors.append((1 + t) ** (-1 - a) * d[w])
        else:
            # (t / (1 + t))^a = t^(-1 + (1 + a)) (1 + t)^(-a).
            exponents.append(1 + a)
            factors.append((1 + t) ** (-a) * d[w])
    t1, t2 = grid.t[:, None], grid.t[None, :]
    f = factors[0][:, None] * factors[1][None, :] * (1 + t1 + t2) ** (-2.5 + sign * nu)
    integrate = [
        grid.regularized if is_curvature else grid.plain for is_curvature in curvature[:2]
    ]
    inner = integrate[0](exponents[0], f, axis=0)
    return integrate[1](exponents[1], inner, axis=0)


def coefficients(
    lam: float, nu: float, interaction: str = "alpha", step: float = 0.04
) -> tuple[float, float]:
    curvature, _, _, phase = INTERACTIONS[interaction]
    grid = Grid(-40.0, UPPER, step)
    d = {w: hyp(nu, lam, w, grid.t) for w in (1, -1)}
    result = []
    for sign in (1, -1):
        total = 0j
        for w1, w2, w3 in itertools.product((1, -1), repeat=3):
            a = [1j * w * lam / 2 for w in (w1, w2, w3)]
            gammas = mp.mpf(1)
            for a_i, is_curvature in zip(a, curvature):
                gammas *= mp.gamma(a_i) if is_curvature else mp.gamma(1 + a_i)
            normalization = complex(phase(w1, w2, w3) * mp.exp(mp.pi * (w1 + w2 + w3) * lam / 4) / gammas)
            long_leg = complex(mp.gamma(2.5 - sign * nu) * mp.gamma(0.5 + sign * nu))
            reduced = c_coefficient(nu, lam, w3, sign) * long_leg * j_integral(
                nu, lam, w1, w2, sign, grid, d, curvature
            )
            amplitudes = np.conj(f0(nu, lam, w1)) * np.conj(f0(nu, lam, w2)) * np.conj(f0(nu, lam, w3))
            total += amplitudes * normalization * reduced
        result.append(total.imag)
    return result[0], result[1]


def prediction(lam: float, nu: float, interaction: str, step: float = 0.04) -> tuple[float, float]:
    """r^2 B / g = A_+ r^(1/2 - nu) + A_- r^(1/2 + nu), with the factor 10/9."""
    _, s, n, _ = INTERACTIONS[interaction]
    plus, minus = coefficients(lam, nu, interaction, step)
    prefactor = s * 28.8 / (2**n * np.sqrt(2.0)) * 10.0 / 9.0
    return prefactor * plus, prefactor * minus


def exact_bispectrum(interaction: str, lam: float, nu: float, k: tuple[float, float, float]) -> float:
    """B(k) of the exact representations, Eqs. (199), (207), (212) with
    (203), (209), (214), in the units of phys_type_models::in_in with unit
    coupling and the late-time limit tau_f -> 0.

    Each operator gives, with X_D = -(i / sqrt(2 k^3)) z^2 e^{iz} g and
    X_S = (2 w / sqrt(2 k^3)) z^2 e^{iz} h,
    B = 2 kappa sum_perm Im int dtau a^4 sum_w prod F^(w)*(tau_f) X^(w)(tau).
    The time integral is int dtau tau^2 e^{-i Omega tau} = -2i / Omega^3 with
    Omega = K + 2 sum k_i t_i, and 2 / Omega^3 = int ds s^2 e^{-s Omega}
    factorizes the three t integrals.
    """
    u_grid = Grid(-40.0, 60.0, 0.05)
    s = np.exp(np.arange(-55.0, 8.0, 0.05))
    ds = 0.05 * s
    k_total = sum(k)
    d = {w: hyp(nu, lam, w, u_grid.t) for w in (1, -1)}
    f0s = {w: f0(nu, lam, w) for w in (1, -1)}

    def leg(kind: str, w: int, kk: float) -> np.ndarray:
        """F(s) = int dt weight(t) e^{-2 s k t} for a curvature (g) or entropy (h) leg."""
        a = 1j * w * lam / 2
        t = u_grid.t[:, None]
        damping = np.exp(-2.0 * s[None, :] * kk * t)
        if kind == "D":
            values = (1 + t) ** (-1 - a) * d[w][:, None] * damping
            integral = u_grid.regularized(a, values, axis=0)
            prefactor = complex(mp.exp(mp.pi * w * lam / 4) / (mp.sqrt(2) * mp.gamma(a)))
        else:
            values = (1 + t) ** (-a) * d[w][:, None] * damping
            integral = u_grid.plain(1 + a, values, axis=0)
            prefactor = complex(mp.exp(mp.pi * w * lam / 4) / (mp.sqrt(2) * mp.gamma(1 + a)))
        return prefactor * integral

    slots = {"alpha": ("S", "D", "D"), "beta": ("D", "S", "S"), "gamma": ("S", "S", "S")}[interaction]
    kappa = {"alpha": 0.5, "beta": 0.5, "gamma": 1.0 / 6.0}[interaction]
    cache: dict[tuple[str, int, int], np.ndarray] = {}
    total = 0.0
    for permutation in itertools.permutations(range(3)):
        kinds = [None, None, None]
        for slot, leg_index in zip(slots, permutation):
            kinds[leg_index] = slot
        for ws in itertools.product((1, -1), repeat=3):
            product = np.ones_like(s, dtype=complex)
            factor = 1.0 + 0j
            for i, (kind, w) in enumerate(zip(kinds, ws)):
                key = (kind, w, i)
                if key not in cache:
                    cache[key] = leg(kind, w, k[i])
                product = product * cache[key]
                # F*(tau_f) X(tau) without the z^2 e^{iz}, which the time
                # integral absorbs: prod_i z_i^2 = (k1 k2 k3)^2 tau^6, a^4 = tau^-4.
                vertex = -1j if kind == "D" else 2 * w
                factor *= np.conj(f0s[w]) * vertex / (2 * k[i] ** 3)
            time_integral = -1j * np.sum(ds * s**2 * np.exp(-s * k_total) * product)
            total += (2 * kappa * factor * (k[0] * k[1] * k[2]) ** 2 * time_integral).imag
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lambda", dest="lam", type=float, required=True)
    parser.add_argument("--nu", type=float, required=True)
    parser.add_argument("--interaction", choices=INTERACTIONS, default="alpha")
    parser.add_argument("--step", type=float, default=0.04)
    parser.add_argument("--exact", type=float, nargs="*", help="r values for the exact representation")
    args = parser.parse_args()
    if args.exact:
        for r in args.exact:
            value = exact_bispectrum(args.interaction, args.lam, args.nu, (1.0, 1.0, r))
            print(f"exact {args.interaction} lambda {args.lam} nu {args.nu} r {r}: B {value:.10e}")
        return
    plus, minus = coefficients(args.lam, args.nu, args.interaction, args.step)
    a_plus, a_minus = prediction(args.lam, args.nu, args.interaction, args.step)
    print(f"{args.interaction} lambda {args.lam} nu {args.nu}: c_+ {plus:.10e} c_- {minus:.10e}")
    print(f"r^2 B / g = {a_plus:.10e} r^(1/2-nu) + {a_minus:.10e} r^(1/2+nu)")


if __name__ == "__main__":
    main()
