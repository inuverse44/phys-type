"""Reproducible derivation of the subhorizon canonical three-point initial condition.

Run with ``python3 docs/derivations/three_point_initial_condition.py`` (requires
SymPy). The script performs four independent checks that are summarized in
``docs/three-point-transport.md``:

1. evaluates the tree-level in-in integral with massless de Sitter mode
   functions and frozen background coefficients, and keeps the highest
   surviving power of k/(aH) separately for every block and every
   contribution (fast A, slow A, three leg-resolved B, three leg-resolved C);
2. compares those results with the closed-form weights implemented in
   ``canonical_three_point_from_subhorizon_vertices``;
3. compares them with the four B_Re^{abc}(t_init) results in Sec. 6 of
   Dias, Frazer, Mulryne & Seery, arXiv:1609.00379, whose momenta are e-fold
   derivatives (our P / H);
4. checks numerically that the initial condition solves the canonical
   transport equation on massless de Sitter with the operator-ordered source,
   with a relative residual that decreases as (k/aH)^-2, also when a
   coefficient drifts in N (rate terms included) and when the linear flow
   carries a small mass.

Conventions: H_int = -(a^3/2)[A QQQ + B QQ Pi + C Pi Pi Q]; alpha^{abc} is the
stripped <X^a(k1) X^b(k2) X^c(k3)>; bPi_r is B with its momentum slot on leg r,
cQ_r is C with its field slot on leg r, a_f = a^2 A_fast, a_s = A - A_fast;
d_c is the e-fold rate dc/dN of coefficient c.
"""
import cmath
import itertools
from collections import defaultdict
from fractions import Fraction

import sympy as sp

STATIC = ["a_f", "a_s", "bPi1", "bPi2", "bPi3", "cQ1", "cQ2", "cQ3"]
# "d_c" is the e-fold rate dc/dN; the coefficient on the slice is expanded as
# c(tau') = c + (dc/dN)(N' - N) with N' - N = 1 - tau'/tau.
RATES = ["d_" + c for c in STATIC]
COEFFS = STATIC + RATES
BLOCKS = ["".join(b) for b in itertools.product("QP", repeat=3)]
TRUNC = 6

# ---------------------------------------------------------------------------
# 1. In-in integral as Laurent polynomials.
#    Monomial exponents: (j, k1, k2, k3, Ki, tau, H), j^2 = -1, Ki = 1/k_t.
# ---------------------------------------------------------------------------
NAMES = ["j", "k1", "k2", "k3", "Ki", "tau", "H"]
KN = ["k1", "k2", "k3"]


def mono(**e):
    v = [0] * len(NAMES)
    for n, x in e.items():
        v[NAMES.index(n)] = x
    return tuple(v)


def poly(d):
    return {k: Fraction(v) for k, v in d.items() if v != 0}


def add(p, q):
    r = defaultdict(Fraction, p)
    for k, v in q.items():
        r[k] += v
    return {k: v for k, v in r.items() if v != 0}


def mul(p, q):
    r = defaultdict(Fraction)
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            r[tuple(a + b for a, b in zip(k1, k2))] += v1 * v2
    return {k: v for k, v in r.items() if v != 0}


ONE = {mono(): Fraction(1)}


def external(slot, r):
    # y(tau) e^{+ik tau} / (H / sqrt(2k^3)): Q -> 1 + i k tau, Pi -> -H k^2 tau^2
    if slot == "Q":
        return poly({mono(): 1, mono(j=1, **{KN[r]: 1}, tau=1): 1})
    return poly({mono(H=1, **{KN[r]: 2}, tau=2): -1})


def internal(slot, r):
    # y*(tau') e^{-ik tau'} / (H / sqrt(2k^3)) as [(coefficient, tau' power)]
    if slot == "Q":
        return [(ONE, 0), (poly({mono(j=1, **{KN[r]: 1}): -1}), 1)]
    return [(poly({mono(H=1, **{KN[r]: 2}): -1}), 2)]


def static_vertex(sigma):
    n_momenta = sigma.count("P")
    if n_momenta == 0:
        # A_fast = a_f H^2 tau'^2 = a_f / a'^2
        return {"a_s": [(ONE, 0)], "a_f": [(poly({mono(H=2): 1}), 2)]}
    if n_momenta == 1:
        return {f"bPi{sigma.index('P') + 1}": [({mono(): Fraction(1, 3)}, 0)]}
    if n_momenta == 2:
        return {f"cQ{sigma.index('Q') + 1}": [({mono(): Fraction(1, 3)}, 0)]}
    return {}


def vertex(sigma):
    out = {}
    for c, terms in static_vertex(sigma).items():
        out[c] = terms
        rate = []
        for p, n in terms:
            rate.append((p, n))
            rate.append((mul(p, poly({mono(tau=-1): -1})), n + 1))
        out["d_" + c] = rate
    return out


def laurent_mul(left, right):
    out = defaultdict(dict)
    for pa, na in left:
        for pb, nb in right:
            out[na + nb] = add(out[na + nb], mul(pa, pb))
    return [(p, n) for n, p in out.items()]


def derive():
    integrand = {c: defaultdict(dict) for c in COEFFS}
    for sigma in itertools.product("QP", repeat=3):
        for c, terms in vertex(sigma).items():
            acc = terms
            for r in range(3):
                acc = laurent_mul(acc, internal(sigma[r], r))
            for p, n in acc:
                # a'^4 = 1 / (H^4 tau'^4)
                integrand[c][n - 4] = add(integrand[c][n - 4], mul(p, poly({mono(H=-4): 1})))

    def integral(g):
        # int_{-inf}^{tau} g e^{i k_t tau'} = e^{i k_t tau} sum_m (-1)^m g^(m)(tau) (-i/k_t)^(m+1)
        res, cur = {}, dict(g)
        for m in range(TRUNC):
            for n, p in cur.items():
                term = mul(p, poly({mono(tau=n, Ki=m + 1, j=m + 1): (-1) ** m * (-1) ** (m + 1)}))
                res = add(res, term)
            nxt = defaultdict(dict)
            for n, p in cur.items():
                if n != 0:
                    nxt[n - 1] = add(nxt[n - 1], {k: v * n for k, v in p.items()})
            cur = nxt
        return res

    k1, k2, k3, H, a = sp.symbols("k1 k2 k3 H a", positive=True)

    def to_sym(p):
        e = 0
        for kk, v in p.items():
            _, e1, e2, e3, eK, et, eH = kk
            e += (sp.Rational(v.numerator, v.denominator) * k1**e1 * k2**e2 * k3**e3
                  * (k1 + k2 + k3) ** (-eK) * (-1 / (a * H)) ** et * H**eH)
        return e

    integrals = {c: integral(integrand[c]) for c in COEFFS}
    results = {}
    for block in BLOCKS:
        ext = ONE
        for r in range(3):
            ext = mul(ext, external(block[r], r))
        results[block] = {}
        for c in COEFFS:
            imaginary = {}
            for kk, v in mul(ext, integrals[c]).items():
                phase = kk[0] % 4
                if phase in (1, 3):
                    imaginary = add(imaginary, {(0,) + kk[1:]: v if phase == 1 else -v})
            by_degree = defaultdict(dict)
            for kk, v in imaginary.items():
                by_degree[kk[1] + kk[2] + kk[3] - kk[4]][kk] = v

            def at(degree):
                # alpha = -6 Im[...] times prod_r H^2/(2 k_r^3)
                return sp.factor(sp.together(
                    -6 * to_sym(by_degree.get(degree, {})) * H**6 / (8 * k1**3 * k2**3 * k3**3)))

            if c in STATIC:
                for degree in sorted(by_degree, reverse=True):
                    expr = at(degree)
                    if expr != 0:
                        results[block][c] = expr
                        results[block]["degree_" + c] = degree
                        break
            else:
                # Keep the rate term at the power of k/(aH) of the retained
                # static term of the same contribution: it is a relative
                # correction of the order of the background-variation rate.
                results[block][c] = at(results[block]["degree_" + c[2:]])
    return results, (k1, k2, k3, H, a)


# ---------------------------------------------------------------------------
# 2. Closed-form weights implemented in Rust.
# ---------------------------------------------------------------------------
def implemented(block, c, k, H, a):
    kt = sum(k)
    product = k[0] * k[1] * k[2]
    pair_sum = k[0] * k[1] + k[0] * k[2] + k[1] * k[2]
    norm = 1 / (4 * product * kt)
    others = {0: (1, 2), 1: (0, 2), 2: (0, 1)}
    moments = [i for i, s in enumerate(block) if s == "P"]
    kind, leg = c[:-1], (int(c[-1]) - 1 if c[-1].isdigit() else None)
    if len(moments) == 0:
        w = {"a_f": 3 / a**4, "a_s": 3 / a**2}.get(c)
        if kind == "bPi":
            w = -H / a**2
        if kind == "cQ":
            u, v = others[leg]
            w = -k[u] * k[v] / a**4
    elif len(moments) == 1:
        p = moments[0]
        u, v = others[p]
        w = {"a_f": -3 * H * (kt + k[p]) / kt / a**4,
             "a_s": -3 * H * (kt - k[p]) / kt / a**2}.get(c)
        if kind == "bPi":
            w = k[p] * k[leg] / a**4
        if kind == "cQ":
            if leg == p:
                w = -H * (k[p]**2 * (k[u] + k[v]) + k[p] * (k[u]**2 + k[v]**2)
                          - k[u] * k[v] * (k[u] + k[v])) / kt / a**4
            else:
                w = -H * k[p]**2 * (k[p] + k[leg]) / kt / a**4
    elif len(moments) == 2:
        q = block.index("Q")
        u, v = others[q]
        w = {"a_f": -3 * k[u] * k[v] / a**6, "a_s": -3 * k[u] * k[v] / a**4}.get(c)
        if kind == "bPi":
            w = (-H * (k[q] * (k[u] + k[v]) - k[u] * k[v]) if leg == q
                 else -H * k[leg]**2) / a**4
        if kind == "cQ":
            if leg == q:
                w = k[u]**2 * k[v]**2 / a**6
            else:
                other = v if leg == u else u
                w = product * k[other] / a**6
    else:
        w = {"a_f": 3 * H * (pair_sum + product / kt) / a**6,
             "a_s": 3 * H * (pair_sum - product / kt) / a**4}.get(c)
        if kind == "bPi":
            w = -product * k[leg] / a**6
        if kind == "cQ":
            u, v = others[leg]
            w = -H * k[u]**2 * k[v]**2 * (kt + k[leg]) / kt / a**6
    return norm * w


def implemented_rate(block, c, k, H, a):
    """Closed-form rate weights used in the Rust assembly (static name c)."""
    kt = sum(k)
    product = k[0] * k[1] * k[2]
    norm = 1 / (4 * product * kt)
    others = {0: (1, 2), 1: (0, 2), 2: (0, 1)}
    moments = [i for i, s in enumerate(block) if s == "P"]
    kind, leg = c[:-1], (int(c[-1]) - 1 if c[-1].isdigit() else None)
    w = 0
    if len(moments) == 0:
        if kind == "bPi":
            w = -H * k[leg] / (kt * a**2)
    elif len(moments) == 1:
        p = moments[0]
        w = {"a_f": 3 / a**4, "a_s": 3 / a**2}.get(c, 0)
        if kind == "cQ":
            w = -(product / k[leg]) / a**4
        w *= H * k[p] / kt
    elif len(moments) == 2:
        q = block.index("Q")
        u, v = others[q]
        if kind == "bPi":
            w = H * k[u] * k[v] * k[leg] / (kt * a**4)
    else:
        w = {"a_f": 3 / a**6, "a_s": 3 / a**4}.get(c, 0)
        if kind == "cQ":
            u, v = others[leg]
            w = -k[u] * k[v] / a**6
        w *= -H * product / kt
    return norm * w


# ---------------------------------------------------------------------------
# 3. Dias, Frazer, Mulryne & Seery (2016), Sec. 6 initial conditions.
# ---------------------------------------------------------------------------
def dfms(symbols):
    k1, k2, k3, H, a = symbols
    ks = [k1, k2, k3]
    af, as_ = sp.symbols("a_f a_s")
    bP = sp.symbols("bPi1:4")
    cQ = sp.symbols("cQ1:4")
    kt = k1 + k2 + k3
    K2 = k1 * k2 + k1 * k3 + k2 * k3

    def perms(f):
        return sum(f([ks[i] for i in p], cQ[p[2]], bP[p[2]], p)
                   for p in itertools.permutations(range(3)))

    fast = af / 2          # six permutations of fast A give 3 a^2 A_fast
    slow = a**2 * as_ / 2
    fff = 1 / (4 * a**4 * k1 * k2 * k3 * kt) * perms(
        lambda q, C, B, p: fast - C * q[0] * q[1] / 2 + slow
        + a**2 * H / 2 * B * ((q[0] + q[1]) * q[2] - K2) / (q[0] * q[1]))
    ffp = 1 / (4 * a**4 * (k1 * k2 * k3)**2 * kt) * (
        k1**2 * (k2 + k3) * perms(lambda q, C, B, p: fast + slow - C * q[0] * q[1] / 2)
        + k1 * perms(lambda q, C, B, p: q[0]**2 * q[1]**2 / 2 * (1 + q[2] / kt) * C
                     + B * q[0] * q[1] * q[2]**2 / (2 * H)
                     - fast * (K2 + k1 * k2 * k3 / kt) - slow * (K2 - k1 * k2 * k3 / kt)))
    fpp = 1 / (4 * a**6 * (k1 * k2 * k3)**2 * kt * H**2) * (
        k1**2 * k2**2 * k3 * perms(lambda q, C, B, p: C * q[0] * q[1] / 2 - fast - slow
                                   - a**2 * H / 2 * B * (q[0] + q[1]) * q[2] / (q[0] * q[1]))
        + k1**2 * k3**2 * perms(lambda q, C, B, p: a**2 * H / 2 * B * q[2]))
    ppp = 1 / (4 * a**6 * k1 * k2 * k3 * kt * H**2) * perms(
        lambda q, C, B, p: fast * (K2 + k1 * k2 * k3 / kt) - bP[p[0]] * q[0]**2 * q[1] * q[2] / (2 * H)
        - q[0]**2 * q[1]**2 / 2 * (1 + q[2] / kt) * C + slow * (K2 - k1 * k2 * k3 / kt))
    return {"QQQ": (fff, 0), "PQQ": (ffp, 1), "PPQ": (fpp, 2), "PPP": (ppp, 3)}


# ---------------------------------------------------------------------------
# 4. Canonical transport residual on massless de Sitter.
# ---------------------------------------------------------------------------
def transport_residual(results, symbols, coefficient, depth, drift=0.0, rates=True,
                       mass_squared_over_h2=0.0):
    """Residual of the canonical transport equation on de Sitter space.

    The coefficient drifts as c(N) = 1 + drift (N - N_0); ``rates`` selects
    whether the initial condition includes the rate terms. The linear flow
    may carry a mass m^2 = mass_squared_over_h2 H^2; the initial condition
    and the covariances remain massless.
    """
    k1, k2, k3, H, a = symbols
    rate_name = "d_" + coefficient

    def alpha_expr(b):
        e = results[b][coefficient]
        if rates:
            e += drift * results[b][rate_name]
        return e

    alpha_f = {b: sp.lambdify((k1, k2, k3, a, H), alpha_expr(b)) for b in BLOCKS}
    # Cosmic-time derivative d/dt = H a d/da at fixed coefficients, plus the
    # explicit drift dc/dt = H dc/dN of the static coefficient.
    dalpha_f = {b: sp.lambdify((k1, k2, k3, a, H),
                               H * a * sp.diff(alpha_expr(b), a)
                               + H * drift * results[b][coefficient])
                for b in BLOCKS}
    scale, hubble = 1.3, 0.7
    k = [x * depth * scale * hubble for x in (1.0, 0.8, 0.5)]
    values = {c: float(c == coefficient) for c in COEFFS}
    A = values["a_f"] / scale**2 + values["a_s"]
    bP = [values[f"bPi{r}"] for r in (1, 2, 3)]
    cQ = [values[f"cQ{r}"] for r in (1, 2, 3)]

    def modes(kk):
        tau = -1 / (scale * hubble)
        f = hubble / (2 * kk**3) ** 0.5 * (1 + 1j * kk * tau) * cmath.exp(-1j * kk * tau)
        p = -hubble**2 * kk**2 * tau**2 / (2 * kk**3) ** 0.5 * cmath.exp(-1j * kk * tau)
        return {"Q": f, "P": p}

    m = [modes(kk) for kk in k]

    def sigma(r, left, right):
        return m[r][left] * m[r][right].conjugate()

    def u2(out, s_out, s1, s2, in1, in2):
        return {("Q", "Q", "Q"): -bP[out], ("Q", "Q", "P"): -cQ[in1],
                ("Q", "P", "Q"): -cQ[in2], ("Q", "P", "P"): 0.0,
                ("P", "Q", "Q"): 3 * A, ("P", "Q", "P"): bP[in2],
                ("P", "P", "Q"): bP[in1], ("P", "P", "P"): cQ[out]}[(s_out, s1, s2)]

    alpha = {b: alpha_f[b](*k, scale, hubble) for b in BLOCKS}
    worst = 0.0
    for b in BLOCKS:
        derivative = dalpha_f[b](*k, scale, hubble)
        magnitude = abs(derivative)
        rhs = 0.0
        for r in range(3):
            if b[r] == "Q":
                term = alpha[b[:r] + "P" + b[r + 1:]]
            else:
                restoring = k[r]**2 / scale**2 + mass_squared_over_h2 * hubble**2
                term = -3 * hubble * alpha[b] - restoring * alpha[b[:r] + "Q" + b[r + 1:]]
            rhs += term
            magnitude += abs(term)
        for out, (i1, i2) in enumerate([(1, 2), (0, 2), (0, 1)]):
            for s1 in "QP":
                for s2 in "QP":
                    # An external leg to the left of the insertion contracts as Sigma^{a d}.
                    left = sigma(i1, b[i1], s1) if i1 < out else sigma(i1, s1, b[i1])
                    right = sigma(i2, b[i2], s2) if i2 < out else sigma(i2, s2, b[i2])
                    term = u2(out, b[out], s1, s2, i1, i2) * left * right
                    rhs += term
                    magnitude += abs(term)
        worst = max(worst, abs(derivative - rhs) / magnitude)
    return worst


def main():
    results, symbols = derive()
    k1, k2, k3, H, a = symbols

    for block in BLOCKS:
        for c in STATIC:
            difference = sp.simplify(results[block][c]
                                     - implemented(block, c, [k1, k2, k3], H, a))
            assert difference == 0, (block, c, difference)
            difference = sp.simplify(results[block]["d_" + c]
                                     - implemented_rate(block, c, [k1, k2, k3], H, a))
            assert difference == 0, ("rate", block, c, difference)
    print("closed-form Rust weights reproduce the in-in integral: ok")

    static_only = {b: {c: results[b][c] for c in STATIC} for b in BLOCKS}
    for block, (paper, n_momenta) in dfms(symbols).items():
        ours = sum(static_only[block][c] * sp.Symbol(c) for c in STATIC) / H**n_momenta
        difference = sp.expand(sp.together(ours - paper))
        for c in STATIC:
            term = sp.factor(sp.together(difference.coeff(sp.Symbol(c))))
            status = "agrees" if term == 0 else f"differs by {term}"
            print(f"DFMS {block} {c}: {status}")

    for c in STATIC:
        shallow = transport_residual(results, symbols, c, 1.0e3)
        deep = transport_residual(results, symbols, c, 1.0e4)
        drift = transport_residual(results, symbols, c, 1.0e4, drift=0.05)
        frozen = transport_residual(results, symbols, c, 1.0e4, drift=0.05, rates=False)
        massive = [transport_residual(results, symbols, c, x, mass_squared_over_h2=0.05)
                   for x in (1.0e3, 1.0e4)]
        print(f"transport residual {c}: {shallow:.2e} (x~1e3), {deep:.2e} (x~1e4); "
              f"drift with rates {drift:.2e}, frozen {frozen:.2e}; "
              f"m^2=0.05H^2: {massive[0]:.2e}, {massive[1]:.2e}")
        assert shallow < 2.0e-5 and deep < 2.0e-7 and drift < 2.0e-7


if __name__ == "__main__":
    main()
