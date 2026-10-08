//! A numerical in-in reference for the constant-mixing two-field theory.
//!
//! The theory is the quadratic action of arXiv:2607.14529, Eq. (10), with
//! constant mixing `lambda` and entropy mass `mu` in exact de Sitter space,
//! plus cubic interactions. Everything here is independent of the engine:
//! the modes solve the paper's Eqs. (24)--(25) directly as complex ODEs in
//! `z = -k tau`, and the tree-level bispectrum is the in-in time integral of
//! the dressed Wightman functions, Eqs. (152) and (198). It is the reference
//! against which the transport amplitudes on the constant-turn orbit are
//! checked; see `docs/constant-turn-bispectrum-amplitudes.md`.
//!
//! Units: `H = 1`. The fields are the canonical curvature perturbation
//! `phi = c zeta`, `c = sqrt(2 epsilon) m_Pl`, and the entropy perturbation
//! `sigma`. A cubic term is given by its coefficient in the Lagrangian
//! `a^3 L_3` times `c`; the reduced `f_NL` is then independent of `c`.

use phys_type::Complex64;

/// The quadratic theory: mixing `lambda` and entropy mass `mu^2 / H^2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixingTheory {
    /// The mixing `lambda` of `D_t phi = phidot - lambda H sigma`.
    pub lambda: f64,
    /// `mu^2 / H^2`.
    pub mu_squared: f64,
}

/// The `z` range and step of a [`ModeTable`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModeSettings {
    /// The Bunch--Davies start, deep inside the horizon.
    pub z_start: f64,
    /// The latest `z` that will be evaluated.
    pub z_end: f64,
    /// The grid step: in `z` above `z = 1` and in `ln z` below it.
    pub step: f64,
}

impl Default for ModeSettings {
    fn default() -> Self {
        Self {
            z_start: 7.5_f64.exp(),
            z_end: (-20.0_f64).exp(),
            step: 1.0e-2,
        }
    }
}

/// `[phi, d phi / dz, sigma, d sigma / dz]` of one mode doublet.
pub type Doublet = [Complex64; 4];

/// The two Bunch--Davies mode doublets of unit wavenumber on a grid in `z`.
///
/// The doublet `b` of wavenumber `k` is `k^{-3/2}` times the tabulated one,
/// evaluated at `z = -k tau`. The grid coordinate is `s = z - 1` above
/// `z = 1` and `s = ln z` below it, so that both the subhorizon oscillation
/// and the superhorizon decay are resolved with a uniform step.
#[derive(Clone, Debug)]
pub struct ModeTable {
    theory: MixingTheory,
    step: f64,
    /// `s` of the first node.
    s_min: f64,
    /// Values and `d/ds` at the nodes, in increasing `s`.
    values: Vec<[Doublet; 2]>,
    rates: Vec<[Doublet; 2]>,
}

fn z_of(s: f64) -> f64 {
    if s >= 0.0 { 1.0 + s } else { s.exp() }
}

fn s_of(z: f64) -> f64 {
    if z >= 1.0 { z - 1.0 } else { z.ln() }
}

/// `d/dz` of a doublet, Eqs. (24)--(25) of arXiv:2607.14529.
fn mode_rate(theory: &MixingTheory, z: f64, y: &Doublet) -> Doublet {
    let [phi, phi_z, sigma, sigma_z] = *y;
    let lambda = theory.lambda;
    let d_phi = phi_z + sigma * (lambda / z);
    let phi_zz = d_phi * (2.0 / z) - phi - sigma_z * (lambda / z) + sigma * (lambda / (z * z));
    let sigma_zz =
        sigma_z * (2.0 / z) - sigma * (1.0 + theory.mu_squared / (z * z)) + d_phi * (lambda / z);
    [phi_z, phi_zz, sigma_z, sigma_zz]
}

fn s_rate(theory: &MixingTheory, s: f64, y: &[Doublet; 2]) -> [Doublet; 2] {
    let z = z_of(s);
    let dz_ds = if s >= 0.0 { 1.0 } else { z };
    y.map(|doublet| mode_rate(theory, z, &doublet).map(|value| value * dz_ds))
}

/// The two positive-frequency solutions at large `z`, the Bunch--Davies
/// doublets. Substituting `phi = z e^{iz} z^a sum_n p_n z^{-n}` and the same
/// for `sigma` with `q_n` into Eqs. (24)--(25), the leading order fixes
/// `a = i s lambda / 2` and `q_0 = -i s p_0` for `s = +-1`, and order `n >= 1`
/// gives
///
/// ```text
/// 2i(a - n) p_n + i lambda q_n = -[(a^2 - (2n-1) a + c_n) p_{n-1} + lambda (a - n - 1) q_{n-1}]
/// 2i(a - n) q_n - i lambda p_n = -[(a^2 - (2n-1) a + c_n - lambda^2 + mu^2) q_{n-1}
///                                  + lambda (n - 2 - a) p_{n-1}]
/// ```
///
/// with `c_n = n (n - 1) - 2`. With `p_0 = 1/2`, each doublet carries the
/// Bunch--Davies norm `|phi|^2 + |sigma|^2 -> z^2 / 2`, and the two are
/// orthogonal. Truncating the series leaves an error `O(z^{-terms})`.
fn bunch_davies(theory: &MixingTheory, z: f64) -> [Doublet; 2] {
    const TERMS: usize = 8;
    let i = Complex64::new(0.0, 1.0);
    let (lambda, mu_squared) = (theory.lambda, theory.mu_squared);
    [1.0, -1.0].map(|sign: f64| {
        let a = i * (0.5 * sign * lambda);
        let mut p = [Complex64::new(0.0, 0.0); TERMS];
        let mut q = p;
        p[0] = Complex64::new(0.5, 0.0);
        q[0] = -i * sign * 0.5;
        for n in 1..TERMS {
            let nf = n as f64;
            let c = nf * (nf - 1.0) - 2.0;
            let common = a * a - a * (2.0 * nf - 1.0) + c;
            let rp = -(common * p[n - 1] + (a - nf - 1.0) * lambda * q[n - 1]);
            let rq = -((common - lambda * lambda + mu_squared) * q[n - 1]
                + (nf - 2.0 - a) * lambda * p[n - 1]);
            // [[d, i lambda], [-i lambda, d]] (p_n, q_n) = (rp, rq).
            let d = (a - nf) * (2.0 * i);
            let determinant = d * d - lambda * lambda;
            p[n] = (d * rp - i * lambda * rq) / determinant;
            q[n] = (d * rq + i * lambda * rp) / determinant;
        }
        // phi = z^{1+a} e^{iz} S_p with S_p = sum p_n z^{-n}, and
        // d phi / dz = z^a e^{iz} [(1 + a + iz) S_p + z S_p'].
        let envelope = (i * z + a * z.ln()).exp();
        let sum = |c: &[Complex64; TERMS]| -> (Complex64, Complex64) {
            let mut value = Complex64::new(0.0, 0.0);
            let mut slope = Complex64::new(0.0, 0.0);
            for (n, coefficient) in c.iter().enumerate() {
                let power = z.powi(-(n as i32));
                value += coefficient * power;
                slope -= coefficient * (n as f64 * power / z);
            }
            (value, slope)
        };
        let ((sp, sp_z), (sq, sq_z)) = (sum(&p), sum(&q));
        [
            envelope * z * sp,
            envelope * ((1.0 + a + i * z) * sp + z * sp_z),
            envelope * z * sq,
            envelope * ((1.0 + a + i * z) * sq + z * sq_z),
        ]
    })
}

impl ModeTable {
    /// Integrates both doublets from `z_start` to `z_end` with the classical
    /// fourth-order Runge--Kutta method on the uniform `s` grid.
    pub fn new(theory: MixingTheory, settings: ModeSettings) -> Self {
        let step = settings.step;
        let first = (s_of(settings.z_end) / step).floor() as i64;
        let last = (s_of(settings.z_start) / step).ceil() as i64;
        let count = (last - first + 1) as usize;
        let mut values = vec![[[Complex64::new(0.0, 0.0); 4]; 2]; count];
        let mut rates = values.clone();
        let mut y = bunch_davies(&theory, z_of(last as f64 * step));
        let add = |y: &[Doublet; 2], k: &[Doublet; 2], scale: f64| -> [Doublet; 2] {
            core::array::from_fn(|b| core::array::from_fn(|c| y[b][c] + k[b][c] * scale))
        };
        for index in (0..count).rev() {
            let s = (first + index as i64) as f64 * step;
            values[index] = y;
            rates[index] = s_rate(&theory, s, &y);
            if index == 0 {
                break;
            }
            // One step of -step in s.
            let h = -step;
            let k1 = rates[index];
            let k2 = s_rate(&theory, s + 0.5 * h, &add(&y, &k1, 0.5 * h));
            let k3 = s_rate(&theory, s + 0.5 * h, &add(&y, &k2, 0.5 * h));
            let k4 = s_rate(&theory, s + h, &add(&y, &k3, h));
            y = core::array::from_fn(|b| {
                core::array::from_fn(|c| {
                    y[b][c] + (k1[b][c] + (k2[b][c] + k3[b][c]) * 2.0 + k4[b][c]) * (h / 6.0)
                })
            });
        }
        Self {
            theory,
            step,
            s_min: first as f64 * step,
            values,
            rates,
        }
    }

    /// The theory the table solves.
    pub fn theory(&self) -> MixingTheory {
        self.theory
    }

    /// The doublets of unit wavenumber at `z`, by cubic Hermite
    /// interpolation in `s`.
    pub fn modes(&self, z: f64) -> [Doublet; 2] {
        let position = (s_of(z) - self.s_min) / self.step;
        assert!(
            position >= 0.0 && position <= (self.values.len() - 1) as f64,
            "z = {z} lies outside the mode table"
        );
        let index = (position.floor() as usize).min(self.values.len() - 2);
        let t = position - index as f64;
        let (h00, h10) = ((1.0 + 2.0 * t) * (1.0 - t).powi(2), t * (1.0 - t).powi(2));
        let (h01, h11) = (t * t * (3.0 - 2.0 * t), t * t * (t - 1.0));
        let (left, right) = (&self.values[index], &self.values[index + 1]);
        let (left_rate, right_rate) = (&self.rates[index], &self.rates[index + 1]);
        core::array::from_fn(|b| {
            core::array::from_fn(|c| {
                left[b][c] * h00
                    + left_rate[b][c] * (self.step * h10)
                    + right[b][c] * h01
                    + right_rate[b][c] * (self.step * h11)
            })
        })
    }

    /// `Delta_zeta / Delta_0` at `z`, with `Delta_0 = H^2 / (8 pi^2 epsilon)`:
    /// `2 k^3 sum_b |phi_b|^2` from Eq. (145).
    pub fn power_ratio(&self, z: f64) -> f64 {
        2.0 * self.modes(z).iter().map(|d| d[0].norm_sqr()).sum::<f64>()
    }
}

/// The real matrix `A(z)` of `dy/dz = A y` for `y = (phi, phi', sigma, sigma')`.
fn rate_matrix(theory: &MixingTheory, z: f64) -> [[f64; 4]; 4] {
    let columns: [Doublet; 4] = core::array::from_fn(|column| {
        let unit: Doublet =
            core::array::from_fn(|c| Complex64::new(if c == column { 1.0 } else { 0.0 }, 0.0));
        mode_rate(theory, z, &unit)
    });
    core::array::from_fn(|row| core::array::from_fn(|column| columns[column][row].re))
}

/// The response of `phi(z_f)` to the state at `z`: the row
/// `r(z) = e_phi T(z_f, z)` of the transfer matrix, `y(z_f) = T(z_f, z) y(z)`.
///
/// The imaginary part of a Wightman function is a commutator, which is
/// exact and small on superhorizon scales; forming it from products of the
/// complex modes cancels to roundoff there. Instead
/// `Im sum_b phi_b(z_f) y_b(z)^* = r(z) S(z)` with the equal-time commutator
/// matrix `S(z)` of [`commutator_matrix`], and `r` solves the adjoint
/// equation `dr/dz = -r A(z)` from `r(z_f) = e_phi`.
#[derive(Clone, Debug)]
struct Response {
    step: f64,
    /// `s` of the first stored node, `z_f` itself.
    s_first: f64,
    values: Vec<[f64; 4]>,
    rates: Vec<[f64; 4]>,
}

impl Response {
    fn new(theory: &MixingTheory, z_final: f64, z_start: f64, step: f64) -> Self {
        let rate = |s: f64, r: &[f64; 4]| -> [f64; 4] {
            let z = z_of(s);
            let dz_ds = if s >= 0.0 { 1.0 } else { z };
            let a = rate_matrix(theory, z);
            core::array::from_fn(|column| {
                -dz_ds * (0..4).map(|row| r[row] * a[row][column]).sum::<f64>()
            })
        };
        let rk4 = |s: f64, r: &[f64; 4], h: f64| -> [f64; 4] {
            let shift = |r: &[f64; 4], k: &[f64; 4], c: f64| -> [f64; 4] {
                core::array::from_fn(|i| r[i] + c * k[i])
            };
            let k1 = rate(s, r);
            let k2 = rate(s + 0.5 * h, &shift(r, &k1, 0.5 * h));
            let k3 = rate(s + 0.5 * h, &shift(r, &k2, 0.5 * h));
            let k4 = rate(s + h, &shift(r, &k3, h));
            core::array::from_fn(|i| r[i] + h / 6.0 * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]))
        };
        // The first node is z_f; later nodes lie on the grid of `ModeTable`,
        // so that no step crosses s = 0.
        let s_final = s_of(z_final);
        let first_node = (s_final / step).floor() as i64 + 1;
        let last_node = (s_of(z_start) / step).ceil() as i64;
        let mut values = vec![[1.0, 0.0, 0.0, 0.0]];
        let mut rates = vec![rate(s_final, &values[0])];
        let mut s = s_final;
        let mut r = values[0];
        for node in first_node..=last_node {
            let target = node as f64 * step;
            r = rk4(s, &r, target - s);
            s = target;
            values.push(r);
            rates.push(rate(s, &r));
        }
        Self {
            step,
            s_first: s_final,
            values,
            rates,
        }
    }

    fn at(&self, z: f64) -> [f64; 4] {
        let s = s_of(z);
        // Node 0 is s_first; node j >= 1 is (floor(s_first / step) + j) step.
        let base = (self.s_first / self.step).floor();
        let (index, left, width) = if s < (base + 1.0) * self.step {
            (0, self.s_first, (base + 1.0) * self.step - self.s_first)
        } else {
            let j = ((s / self.step).floor() - base) as usize;
            let j = j.min(self.values.len() - 2);
            (j, (base + j as f64) * self.step, self.step)
        };
        assert!(
            index + 1 < self.values.len(),
            "z = {z} lies outside the response"
        );
        let t = ((s - left) / width).clamp(0.0, 1.0);
        let (h00, h10) = ((1.0 + 2.0 * t) * (1.0 - t).powi(2), t * (1.0 - t).powi(2));
        let (h01, h11) = (t * t * (3.0 - 2.0 * t), t * t * (t - 1.0));
        core::array::from_fn(|c| {
            self.values[index][c] * h00
                + self.rates[index][c] * (width * h10)
                + self.values[index + 1][c] * h01
                + self.rates[index + 1][c] * (width * h11)
        })
    }
}

/// `S_ac(z) = Im sum_b y_a,b y_c,b^*` for unit wavenumber, from the canonical
/// commutators. With `H = 1` and `a = 1/z`, `Pi_phi = a^3 D_t phi =
/// -(phi' + lambda sigma / z) / z^2` and `Pi_sigma = -sigma' / z^2`, and
/// `Im sum_b phi_b Pi_phi,b^* = Im sum_b sigma_b Pi_sigma,b^* = 1/2`.
fn commutator_matrix(theory: &MixingTheory, z: f64) -> [[f64; 4]; 4] {
    // y = J Z for Z = (phi, Pi_phi, sigma, Pi_sigma).
    let j = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, -z * z, -theory.lambda / z, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, -z * z],
    ];
    let c = [
        [0.0, 0.5, 0.0, 0.0],
        [-0.5, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.5],
        [0.0, 0.0, -0.5, 0.0],
    ];
    core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            (0..16)
                .map(|m| {
                    let [p, q] = [m / 4, m % 4];
                    j[a][p] * c[p][q] * j[b][q]
                })
                .sum()
        })
    })
}

/// A field in a cubic vertex.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Slot {
    Phi,
    PhiDot,
    /// `D_t phi = phidot - lambda sigma`.
    CovariantPhiDot,
    Sigma,
}

/// A cubic operator of `L_3`, with `H = 1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CubicOperator {
    /// `sigma phidot^2`.
    SigmaPhiDotSquared,
    /// `sigma (nabla phi)^2 / a^2`.
    SigmaGradientPhiSquared,
    /// `sigma^2 phidot`.
    SigmaSquaredPhiDot,
    /// `sigma^3`.
    SigmaCubed,
    /// `sigma (D_t phi)^2`.
    SigmaCovariantPhiDotSquared,
    /// `sigma^2 D_t phi`.
    SigmaSquaredCovariantPhiDot,
    /// `phidot^3`, which has a closed form at `lambda = 0`.
    PhiDotCubed,
}

impl CubicOperator {
    fn slots(self) -> [Slot; 3] {
        use Slot::*;
        match self {
            Self::SigmaPhiDotSquared => [Sigma, PhiDot, PhiDot],
            Self::SigmaGradientPhiSquared => [Sigma, Phi, Phi],
            Self::SigmaSquaredPhiDot => [Sigma, Sigma, PhiDot],
            Self::SigmaCubed => [Sigma, Sigma, Sigma],
            Self::SigmaCovariantPhiDotSquared => [Sigma, CovariantPhiDot, CovariantPhiDot],
            Self::SigmaSquaredCovariantPhiDot => [Sigma, Sigma, CovariantPhiDot],
            Self::PhiDotCubed => [PhiDot, PhiDot, PhiDot],
        }
    }
}

/// One term `coefficient * operator` of `L_3`, the coefficient multiplied
/// by `c = sqrt(2 epsilon) m_Pl`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicTerm {
    /// The operator.
    pub operator: CubicOperator,
    /// Its coefficient times `c`.
    pub coefficient: f64,
}

/// The cubic terms of the circular orbit of `PolarFieldSpace` in the
/// decoupling limit: `(sigma / R)(phidot^2 - (nabla phi)^2 / a^2) -
/// (Omega / R) sigma^2 phidot` with `c / R = Omega = lambda / 2`.
pub fn orbit_terms(lambda: f64) -> [CubicTerm; 3] {
    [
        CubicTerm {
            operator: CubicOperator::SigmaPhiDotSquared,
            coefficient: 0.5 * lambda,
        },
        CubicTerm {
            operator: CubicOperator::SigmaGradientPhiSquared,
            coefficient: -0.5 * lambda,
        },
        CubicTerm {
            operator: CubicOperator::SigmaSquaredPhiDot,
            coefficient: -0.25 * lambda * lambda,
        },
    ]
}

/// The interactions of arXiv:2607.14529, Eq. (3), with couplings
/// `(alpha~, beta~, gamma~)`: `c` times their coefficients is
/// `(alpha~ / 2, beta~ / 2, gamma~ / 6)`.
pub fn paper_terms(alpha: f64, beta: f64, gamma: f64) -> [CubicTerm; 3] {
    [
        CubicTerm {
            operator: CubicOperator::SigmaPhiDotSquared,
            coefficient: 0.5 * alpha,
        },
        CubicTerm {
            operator: CubicOperator::SigmaSquaredPhiDot,
            coefficient: 0.5 * beta,
        },
        CubicTerm {
            operator: CubicOperator::SigmaCubed,
            coefficient: gamma / 6.0,
        },
    ]
}

/// The same interactions in the covariant basis of arXiv:2607.14529,
/// Eq. (8), `sigma (D_t zeta)^2`, `sigma^2 D_t zeta`, and `sigma^3`, with
/// couplings `(alpha, beta, gamma)` and the normalizations of Eq. (3).
pub fn covariant_paper_terms(alpha: f64, beta: f64, gamma: f64) -> [CubicTerm; 3] {
    [
        CubicTerm {
            operator: CubicOperator::SigmaCovariantPhiDotSquared,
            coefficient: 0.5 * alpha,
        },
        CubicTerm {
            operator: CubicOperator::SigmaSquaredCovariantPhiDot,
            coefficient: 0.5 * beta,
        },
        CubicTerm {
            operator: CubicOperator::SigmaCubed,
            coefficient: gamma / 6.0,
        },
    ]
}

/// The tree-level bispectrum of `phi` (in units where the coefficients
/// carry `c`) and the power spectra of its legs on the final slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InInBispectrum {
    /// The wavenumbers.
    pub k: [f64; 3],
    /// `c B_phi`.
    pub bispectrum: f64,
    /// `P_phi(k_i) = sum_b |phi_b|^2`.
    pub power: [f64; 3],
}

impl InInBispectrum {
    /// `(5/6) B / (P_1 P_2 + P_1 P_3 + P_2 P_3)`, the convention of
    /// `phys_type::CurvatureCorrelators::reduced_f_nl`.
    pub fn reduced_f_nl(&self) -> f64 {
        let [p1, p2, p3] = self.power;
        5.0 / 6.0 * self.bispectrum / (p1 * p2 + p1 * p3 + p2 * p3)
    }
}

/// Integration settings of [`bispectrum`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VertexSettings {
    /// Panel width at most this times `1/K`.
    pub oscillation_width: f64,
    /// Panel width at most this times `|tau|`.
    pub logarithmic_width: f64,
    /// The window switches off between this fraction of the earliest
    /// `K |tau|` and the earliest one.
    pub window_onset: f64,
    /// Below this `z` of a leg, the imaginary part of its Wightman function
    /// comes from the exact commutator instead of the modes.
    pub response_below: f64,
}

impl Default for VertexSettings {
    fn default() -> Self {
        Self {
            oscillation_width: 1.0,
            logarithmic_width: 0.05,
            window_onset: 0.5,
            response_below: 0.1,
        }
    }
}

/// Eight-point Gauss--Legendre nodes and weights on `[-1, 1]`.
const GAUSS_LEGENDRE: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_2, 0.101_228_536_290_376_26),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_47),
    (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_47),
    (0.960_289_856_497_536_2, 0.101_228_536_290_376_26),
];

/// A smooth step from 1 at `t <= 0` to 0 at `t >= 1`, with every derivative
/// continuous, so that switching off the early-time integrand converges
/// faster than any power of the window length.
fn window(t: f64) -> f64 {
    if t <= 0.0 {
        return 1.0;
    }
    if t >= 1.0 {
        return 0.0;
    }
    let f = |x: f64| (-1.0 / x).exp();
    f(1.0 - t) / (f(1.0 - t) + f(t))
}

const PERMUTATIONS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// The tree-level bispectrum `<phi_k1 phi_k2 phi_k3>` of `terms` on the
/// slice where the shortest leg has `z = z_final`; the table must reach the
/// longest leg there, `z_final k_min / k_max`.
///
/// `<Q(tau_f)> = 2 Im int dt <Q(tau_f) H_I(t)>` with `H_I = -int d^3x a^3 L_3`
/// gives `B = -2 sum_terms g Im int dtau a^4 sum_perm prod_i G_i(slot)`, where
/// `G_i(X) = sum_b phi_b(k_i, tau_f) X_b^*(k_i, tau)` is the Wightman function
/// of Eq. (152). The integral starts where the shortest leg has the table's
/// `z_start`; a smooth window replaces the `i0^+` prescription there.
pub fn bispectrum(
    table: &ModeTable,
    terms: &[CubicTerm],
    k: [f64; 3],
    z_final: f64,
    settings: &VertexSettings,
) -> InInBispectrum {
    let lambda = table.theory.lambda;
    let k_max = k.iter().copied().fold(0.0, f64::max);
    let total = k[0] + k[1] + k[2];
    let z_start = z_of((table.values.len() - 1) as f64 * table.step + table.s_min);
    let tau_start = -z_start / k_max;
    let tau_final = -z_final / k_max;
    let x_end = total * tau_start.abs();
    let x_onset = settings.window_onset * x_end;

    let leg_values = |leg: usize, tau: f64| -> [Doublet; 2] { table.modes(-k[leg] * tau) };
    let external: [[Complex64; 2]; 3] = core::array::from_fn(|leg| {
        let modes = leg_values(leg, tau_final);
        let scale = k[leg].powf(-1.5);
        [modes[0][0] * scale, modes[1][0] * scale]
    });
    let power: [f64; 3] =
        core::array::from_fn(|leg| external[leg].iter().map(Complex64::norm_sqr).sum::<f64>());
    let responses: [Response; 3] = core::array::from_fn(|leg| {
        Response::new(&table.theory, -k[leg] * tau_final, z_start, table.step)
    });

    // The integrand a^4 w sum_terms g sum_perm prod G at one time.
    let integrand = |tau: f64| -> f64 {
        let propagators: [[Complex64; 4]; 3] = core::array::from_fn(|leg| {
            let z = -k[leg] * tau;
            let modes = leg_values(leg, tau);
            let scale = k[leg].powf(-1.5);
            let slot_value = |doublet: &Doublet, slot: usize| -> Complex64 {
                let [phi, phi_z, sigma, _] = *doublet;
                match slot {
                    0 => phi,
                    1 => -phi_z * z,
                    2 => -(phi_z + sigma * (lambda / z)) * z,
                    _ => sigma,
                }
            };
            // Inside the horizon G comes from the modes. Outside, Im G is
            // far below Re G and the products of modes cancel to roundoff, so
            // it comes from the exact commutator, r(z) S(z) contracted with
            // the slot. The modes stay in use inside the horizon because the
            // squeezed limit relies on a cancellation between times of the
            // short legs' Im G that is exact only when Re and Im share their
            // integration errors.
            let response = responses[leg].at(z);
            let commutator = commutator_matrix(&table.theory, z);
            let imaginary: [f64; 4] =
                core::array::from_fn(|c| (0..4).map(|a| response[a] * commutator[a][c]).sum());
            let slot_row = |slot: usize| -> [f64; 4] {
                match slot {
                    0 => [1.0, 0.0, 0.0, 0.0],
                    1 => [0.0, -z, 0.0, 0.0],
                    2 => [0.0, -z, -lambda, 0.0],
                    _ => [0.0, 0.0, 1.0, 0.0],
                }
            };
            core::array::from_fn(|slot| {
                let real: f64 = (0..2)
                    .map(|b| (external[leg][b] * (slot_value(&modes[b], slot) * scale).conj()).re)
                    .sum();
                let row = slot_row(slot);
                let im: f64 = if z < settings.response_below {
                    (0..4).map(|c| row[c] * imaginary[c]).sum::<f64>() / k[leg].powi(3)
                } else {
                    (0..2)
                        .map(|b| {
                            (external[leg][b] * (slot_value(&modes[b], slot) * scale).conj()).im
                        })
                        .sum()
                };
                Complex64::new(real, im)
            })
        });
        let index = |slot: Slot| match slot {
            Slot::Phi => 0,
            Slot::PhiDot => 1,
            Slot::CovariantPhiDot => 2,
            Slot::Sigma => 3,
        };
        let dot = |a: usize, b: usize| {
            let c = 3 - a - b;
            0.5 * (k[c] * k[c] - k[a] * k[a] - k[b] * k[b])
        };
        let mut sum = Complex64::new(0.0, 0.0);
        for term in terms {
            let slots = term.operator.slots();
            let mut contraction = Complex64::new(0.0, 0.0);
            for permutation in PERMUTATIONS {
                let mut product = Complex64::new(1.0, 0.0);
                for (slot, &leg) in slots.iter().zip(permutation.iter()) {
                    product *= propagators[leg][index(*slot)];
                }
                if term.operator == CubicOperator::SigmaGradientPhiSquared {
                    // (nabla phi)^2 / a^2 -> -(k_a . k_b) tau^2 for the legs
                    // of the two phi slots.
                    product *= -dot(permutation[1], permutation[2]) * tau * tau;
                }
                contraction += product;
            }
            sum += contraction * term.coefficient;
        }
        let x = total * tau.abs();
        let weight = window((x - x_onset) / (x_end - x_onset));
        sum.im * weight / tau.powi(4)
    };

    let mut integral = 0.0;
    let mut left = tau_start;
    while left < tau_final {
        let width = (settings.oscillation_width / total)
            .min(settings.logarithmic_width * left.abs())
            .min(tau_final - left);
        let (middle, half) = (left + 0.5 * width, 0.5 * width);
        integral += GAUSS_LEGENDRE
            .iter()
            .map(|(node, weight)| weight * integrand(middle + half * node))
            .sum::<f64>()
            * half;
        left += width;
    }

    InInBispectrum {
        k,
        bispectrum: -2.0 * integral,
        power,
    }
}

/// `B = 3 g / (k_1 k_2 k_3 K^3)` for `g phidot^3` at `lambda = 0`, from
/// `G(phidot) = -z^2 e^{-iz} / (2 k^3)` and `int tau^2 e^{iK tau} = 2i / K^3`.
pub fn phi_dot_cubed_bispectrum(coefficient: f64, k: [f64; 3]) -> f64 {
    let total = k[0] + k[1] + k[2];
    3.0 * coefficient / (k[0] * k[1] * k[2] * total.powi(3))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residual(theory: &MixingTheory, z: f64) -> f64 {
        // Eqs. (24)--(25) on the asymptotic solutions, by central differences.
        let h = 1.0e-4;
        let at = |z: f64| bunch_davies(theory, z);
        let (minus, centre, plus) = (at(z - h), at(z), at(z + h));
        let mut worst: f64 = 0.0;
        for b in 0..2 {
            let [phi, phi_z, sigma, sigma_z] = centre[b];
            let first_phi = (plus[b][0] - minus[b][0]) / (2.0 * h);
            let first_sigma = (plus[b][2] - minus[b][2]) / (2.0 * h);
            let rate = mode_rate(theory, z, &centre[b]);
            let second_phi = (plus[b][1] - minus[b][1]) / (2.0 * h);
            let second_sigma = (plus[b][3] - minus[b][3]) / (2.0 * h);
            let scale = phi.norm() + sigma.norm();
            for error in [
                first_phi - phi_z,
                first_sigma - sigma_z,
                second_phi - rate[1],
                second_sigma - rate[3],
            ] {
                worst = worst.max(error.norm() / scale);
            }
        }
        worst
    }

    #[test]
    fn the_asymptotic_doublets_solve_the_mode_equations() {
        for (lambda, mu_squared) in [(0.0, 0.0), (1.0, 1.25), (2.0, 4.0)] {
            let theory = MixingTheory { lambda, mu_squared };
            assert!(residual(&theory, 300.0) < 1.0e-7, "{lambda} {mu_squared}");
        }
    }

    #[test]
    fn the_asymptotic_doublets_are_bunch_davies_normalized() {
        let theory = MixingTheory {
            lambda: 1.5,
            mu_squared: 2.0,
        };
        let z = 1.0e4;
        let doublets = bunch_davies(&theory, z);
        for d in &doublets {
            let norm = d[0].norm_sqr() + d[2].norm_sqr();
            assert!((norm / (0.5 * z * z) - 1.0).abs() < 1.0e-3);
        }
        let overlap =
            doublets[0][0].conj() * doublets[1][0] + doublets[0][2].conj() * doublets[1][2];
        assert!(overlap.norm() < 1.0e-3 * z * z);
    }

    #[test]
    fn massless_modes_without_mixing_freeze_at_the_bunch_davies_amplitude() {
        let theory = MixingTheory {
            lambda: 0.0,
            mu_squared: 0.0,
        };
        let table = ModeTable::new(
            theory,
            ModeSettings {
                z_start: 400.0,
                z_end: 1.0e-4,
                step: 1.0e-2,
            },
        );
        // Each doublet carries half of u_0, whose frozen |u_0|^2 = 1/2.
        assert!((table.power_ratio(1.0e-4) - 1.0).abs() < 1.0e-7);
        let early = bunch_davies(&theory, 150.0);
        let tabulated = table.modes(150.0);
        for b in 0..2 {
            for c in 0..4 {
                assert!((tabulated[b][c] - early[b][c]).norm() < 1.0e-7 * 150.0);
            }
        }
    }

    #[test]
    fn the_window_is_a_smooth_step() {
        assert_eq!(window(-0.1), 1.0);
        assert_eq!(window(1.2), 0.0);
        assert!((window(0.5) - 0.5).abs() < 1.0e-15);
        assert!((window(0.3) + window(0.7) - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn the_commutator_matrix_matches_the_modes() {
        let theory = MixingTheory {
            lambda: 1.7,
            mu_squared: 1.1,
        };
        let table = ModeTable::new(theory, ModeSettings::default());
        for z in [40.0, 3.0, 0.5] {
            let modes = table.modes(z);
            let exact = commutator_matrix(&theory, z);
            for a in 0..4 {
                for c in 0..4 {
                    let direct: f64 = (0..2).map(|b| (modes[b][a] * modes[b][c].conj()).im).sum();
                    assert!(
                        (direct - exact[a][c]).abs() < 1.0e-7 * (1.0 + z * z),
                        "{z} {a} {c}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_response_reproduces_the_wightman_imaginary_part() {
        let theory = MixingTheory {
            lambda: 1.7,
            mu_squared: 1.1,
        };
        let settings = ModeSettings::default();
        let table = ModeTable::new(theory, settings);
        let z_final = 0.05;
        let response = Response::new(&theory, z_final, settings.z_start, settings.step);
        let external = table.modes(z_final);
        for z in [30.0, 2.0, 0.3, 0.06] {
            let modes = table.modes(z);
            let r = response.at(z);
            let s = commutator_matrix(&theory, z);
            for c in 0..4 {
                let direct: f64 = (0..2)
                    .map(|b| (external[b][0] * modes[b][c].conj()).im)
                    .sum();
                let via: f64 = (0..4).map(|a| r[a] * s[a][c]).sum();
                assert!(
                    (direct - via).abs() < 1.0e-6 * (1.0 + z * z),
                    "{z} {c}: {direct} vs {via}"
                );
            }
        }
    }
}
