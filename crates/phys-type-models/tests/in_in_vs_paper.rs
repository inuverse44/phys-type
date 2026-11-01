//! The in-in reference against the closed forms of arXiv:2607.14529. See
//! `docs/constant-turn-bispectrum-amplitudes.md`.

use phys_type_models::in_in::{
    MixingTheory, ModeSettings, ModeTable, VertexSettings, bispectrum, covariant_paper_terms,
};

/// `r^2 B(1, 1, r)` of the covariant interactions with `couplings` on the slice
/// `z_final` of the short legs.
fn squeezed(lambda: f64, nu: f64, couplings: [f64; 3], ratios: &[f64], z_final: f64) -> Vec<f64> {
    let table = ModeTable::new(
        MixingTheory {
            lambda,
            mu_squared: 2.25 - nu * nu,
        },
        ModeSettings {
            z_end: z_final * ratios.iter().copied().fold(1.0, f64::min) * 0.5,
            ..ModeSettings::default()
        },
    );
    ratios
        .iter()
        .map(|&r| {
            r * r
                * bispectrum(
                    &table,
                    &covariant_paper_terms(couplings[0], couplings[1], couplings[2]),
                    [1.0, 1.0, r],
                    z_final,
                    &VertexSettings::default(),
                )
                .bispectrum
        })
        .collect()
}

/// Least squares for `y = sum_j c_j x^{p_j}` by normal equations.
fn fit(x: &[f64], y: &[f64], powers: &[f64]) -> Vec<f64> {
    let n = powers.len();
    let mut a = vec![vec![0.0; n + 1]; n];
    for (&xi, &yi) in x.iter().zip(y) {
        let basis: Vec<f64> = powers.iter().map(|&p| xi.powf(p)).collect();
        // Weight each point by its leading term so that small r count.
        let weight = 1.0 / (yi * yi);
        for i in 0..n {
            for j in 0..n {
                a[i][j] += weight * basis[i] * basis[j];
            }
            a[i][n] += weight * basis[i] * yi;
        }
    }
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&p, &q| a[p][column].abs().total_cmp(&a[q][column].abs()))
            .unwrap();
        a.swap(column, pivot);
        let pivot_row = a[column].clone();
        for (row, entries) in a.iter_mut().enumerate() {
            if row != column {
                let factor = entries[column] / pivot_row[column];
                for (entry, pivot_entry) in entries.iter_mut().zip(&pivot_row).skip(column) {
                    *entry -= factor * pivot_entry;
                }
            }
        }
    }
    (0..n).map(|i| a[i][n] / a[i][i]).collect()
}

/// The leading squeezed coefficient `A_+` of `r^2 B(1, 1, r) = A_+ r^{1/2 - nu} + ...`
/// from the closed forms of arXiv:2607.14529, Eqs. (240), (248), and (255)
/// with the factor 10/9 that they miss (`docs/literature-discrepancies.md`
/// entry 12), evaluated by
/// `uv run scripts/hmps2_squeezed.py --interaction NAME --lambda L --nu NU`.
const PAPER_LEADING: [(&str, [f64; 3], f64, f64, f64); 15] = [
    ("alpha", [1.0, 0.0, 0.0], 0.05, 0.8, 2.3197251407e-3),
    ("alpha", [1.0, 0.0, 0.0], 0.5, 0.6, 6.6676967947e-2),
    ("alpha", [1.0, 0.0, 0.0], 1.0, 0.6, 1.9481076737e0),
    ("alpha", [1.0, 0.0, 0.0], 1.0, 0.8, 2.7956950811e0),
    ("alpha", [1.0, 0.0, 0.0], 2.0, 0.8, 1.2178470358e3),
    ("alpha", [1.0, 0.0, 0.0], 1.0, 1.2, 2.0106616226e1),
    ("alpha", [1.0, 0.0, 0.0], 2.0, 1.2, 7.9351103115e3),
    ("alpha", [1.0, 0.0, 0.0], 3.0, 0.8, 5.3468978345e5),
    ("alpha", [1.0, 0.0, 0.0], 1.0, 0.4, 1.8618749945e0),
    ("beta", [0.0, 1.0, 0.0], 1.0, 0.4, -5.2068045175e0),
    ("beta", [0.0, 1.0, 0.0], 1.0, 0.8, -7.1846384711e0),
    ("beta", [0.0, 1.0, 0.0], 2.0, 0.8, -1.6800854349e3),
    ("gamma", [0.0, 0.0, 1.0], 1.0, 0.4, 2.1549435360e0),
    ("gamma", [0.0, 0.0, 1.0], 1.0, 0.3, 2.6276526945e0),
    ("gamma", [0.0, 0.0, 1.0], 2.0, 0.3, 7.2919664385e2),
];

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn squeezed_in_in_matches_the_corrected_closed_forms() {
    // The in-in r^2 B over 1e-6 <= r <= 1e-2 is fitted with the leading and
    // subleading nonanalytic powers and the analytic r and r^2; the leading
    // coefficient must equal the paper's, times 10/9.
    let ratios: Vec<f64> = (0..=32)
        .map(|i| 10f64.powf(-2.0 - 0.125 * f64::from(i)))
        .collect();
    for (name, couplings, lambda, nu, paper) in PAPER_LEADING {
        let y = squeezed(lambda, nu, couplings, &ratios, (-30.0_f64).exp());
        let leading = fit(&ratios, &y, &[0.5 - nu, 0.5 + nu, 1.0, 2.0])[0];
        assert!(
            (leading / paper - 1.0).abs() < 5.0e-4,
            "{name} lambda {lambda} nu {nu}: in-in {leading} vs paper {paper}"
        );
        // Without the factor 10/9 the paper is 10% low everywhere.
        assert!((leading / (0.9 * paper) - 1.0).abs() > 0.1);
    }
}

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn in_in_matches_the_exact_integral_representations() {
    // Eqs. (199)/(203), (207)/(209), and (212)/(214) of arXiv:2607.14529 at
    // finite r, before any squeezed approximation, evaluated by
    // `uv run scripts/hmps2_squeezed.py --interaction NAME --lambda L --nu NU --exact 1 0.1`.
    // (interaction, couplings, lambda, nu, B at r = 1, B at r = 0.1)
    type Case = (&'static str, [f64; 3], f64, f64, f64, f64);
    let exact: [Case; 5] = [
        (
            "gamma",
            [0.0, 0.0, 1.0],
            1.0,
            0.3,
            7.3183777132e-1,
            9.4991246145e1,
        ),
        (
            "gamma",
            [0.0, 0.0, 1.0],
            2.0,
            0.3,
            9.0142812949e1,
            1.8316158496e4,
        ),
        (
            "beta",
            [0.0, 1.0, 0.0],
            1.0,
            0.4,
            -2.6637892249e0,
            -3.1362557454e2,
        ),
        (
            "beta",
            [0.0, 1.0, 0.0],
            1.0,
            0.8,
            -7.9221972295e0,
            -1.3856239080e3,
        ),
        (
            "alpha",
            [1.0, 0.0, 0.0],
            1.0,
            0.8,
            5.1571483604e0,
            6.1307268500e2,
        ),
    ];
    for (name, couplings, lambda, nu, at_one, at_tenth) in exact {
        for (r, paper) in [(1.0, at_one), (0.1, at_tenth)] {
            let b = squeezed(lambda, nu, couplings, &[r], (-60.0_f64).exp())[0] / (r * r);
            assert!(
                (b / paper - 1.0).abs() < 1.0e-6,
                "{name} lambda {lambda} nu {nu} r {r}: in-in {b} vs paper {paper}"
            );
        }
    }
}

/// The squeezed coefficient `s(nu)` of Chen and Wang, arXiv:0911.3380,
/// Eq. (5.11), evaluated by `uv run scripts/chen_wang_squeezed.py --nu NU`,
/// with the relative tolerance of each value. For `nu = 1.2` the integrand
/// near `s = 0` converges only as `s^{3/2 - nu}` and roundoff limits the
/// cutoff, so the value is extrapolated in the cutoff and uncertain at 1e-3.
const CHEN_WANG: [(f64, f64, f64); 5] = [
    (0.4, 3.107_12, 5.0e-5),
    (0.6, 3.328_87, 5.0e-5),
    (0.8, 5.438_05, 5.0e-5),
    (1.0, 14.395_7, 5.0e-5),
    (1.2, 84.48, 3.0e-3),
];

#[test]
#[ignore = "expensive; run with `cargo test --release -- --ignored`"]
fn weak_mixing_sigma_cubed_matches_chen_wang() {
    // Quasi-single-field inflation (arXiv:0911.3380) is the constant-turn
    // theory at first order in each transfer. For p3 << p1 = p2 its Eq. (5.10)
    // gives <dtheta^3> = s(nu) c2^3 c3 / (H R^6) p1^(nu-7/2) p2^-1 p3^(-3/2-nu),
    // with c2 = 2 R thetadot = lambda R and L_3 = -c3 sigma^3. With
    // phi = -R dtheta and c = R thetadot, a unit coefficient of sigma^3 then
    // has r^2 B = s(nu) lambda^3 r^(1/2 - nu) + ..., with B the in-in
    // bispectrum of phi times c. (The paper labels the left side of Eq. (5.10)
    // zeta; only <dtheta^3>, as in its letter version, arXiv:0909.0496,
    // Eq. (13), scales as lambda^3.)
    let ratios: Vec<f64> = (0..=32)
        .map(|i| 10f64.powf(-2.0 - 0.125 * f64::from(i)))
        .collect();
    for (nu, s, tolerance) in CHEN_WANG {
        let leading = |lambda: f64| {
            // A unit coefficient of sigma^3 is gamma = 6.
            let y = squeezed(lambda, nu, [0.0, 0.0, 6.0], &ratios, (-40.0_f64).exp());
            fit(&ratios, &y, &[0.5 - nu, 0.5 + nu, 1.0, 2.0])[0] / lambda.powi(3)
        };
        // Remove the O(lambda^2) correction.
        let (small, larger) = (leading(0.01), leading(0.02));
        let limit = (4.0 * small - larger) / 3.0;
        assert!(
            (limit / s - 1.0).abs() < tolerance,
            "nu {nu}: in-in {limit} vs Chen-Wang {s}"
        );
    }
}
