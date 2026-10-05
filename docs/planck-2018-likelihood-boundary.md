# Planck 2018 observational boundary

## Purpose

`phys-type` computes primordial spectra.  It does not evolve the spectra
through recombination, calculate CMB angular spectra, evaluate a released
likelihood, or infer posterior distributions.  This document fixes the
observational reference attached to the two-field benchmark without crossing
that boundary.

The benchmark's `A_s = 2.1e-9` is a calibration target.  It is not a value
fitted by this repository and must not be reported as a likelihood result.

## Named baseline

The reference is the Planck 2018 legacy PR3 baseline combination
`TT,TE,EE+lowE`:

- high-multipole Plik `TT`, `TE`, and `EE` likelihoods for `ell >= 30`;
- low-multipole Commander temperature likelihood;
- low-multipole SimAll `EE` likelihood (`lowE`).

The associated released likelihood archive is
`COM_Likelihood_Data-baseline_R3.00`.  The CMB lensing-reconstruction
likelihood is not part of this named baseline.  A future inference layer must
add lensing explicitly if it uses it.

The scalar-spectrum pivot is

```text
k_0 = 0.05 Mpc^-1.
```

This physical observational pivot is not numerically identified with the
solver's `k_pivot_internal_over_m_pl`.  The latter depends on the arbitrary
normalization of the simulated scale factor.  The public
`ComovingWavenumberMpcInverse` newtype keeps the observational coordinate
separate from the internal Planck-unit `Mass` type.

Primary references:

- [Planck 2018 results V: CMB power spectra and likelihoods](https://arxiv.org/abs/1907.12875)
- [Planck 2018 results VI: Cosmological parameters](https://arxiv.org/abs/1807.06209)
- [Planck 2018 results X: Constraints on inflation](https://arxiv.org/abs/1807.06211)

## Machine-readable contract

Every CSV from `two_field_spectra` records:

```text
# observational_reference: planck_2018_tt_te_ee_low_e
# observational_data_release: planck_2018_legacy_pr3
# likelihood_release: COM_Likelihood_Data-baseline_R3.00
# likelihood_components: plik_high_l_tt_te_ee,commander_low_l_tt,simall_low_l_ee
# scalar_pivot_mpc_inverse: 5e-2
# amplitude_role: benchmark_target_not_likelihood_fit
```

The typed source of this block is `PLANCK_2018_TT_TE_EE_LOW_E`.  Lensing is a
separate `LikelihoodComponent`, so it cannot appear in the baseline by an
undocumented convention.

## Published contours as an overlay

`plots/higgs_ns_r.py` draws the 68% and 95% $n_s$–$r_{0.002}$ contours from
the BICEP/Keck 2018 data release
([arXiv:2110.00483](https://arxiv.org/abs/2110.00483), `rns_code.tgz`), for the
combinations "Planck TT,TE,EE+lowE+lensing" and "+BK18+BAO".
`scripts/fetch_bk18.py` fetches the release and pins its SHA-256; the grids
are not redistributed. The overlay places predicted points against published
posteriors that were derived for $\Lambda$CDM+$r$. It does not evaluate any
likelihood, so the boundary above still applies: falling inside a contour is
a benchmark comparison, not a fit.

## Future inference boundary

A genuine fit requires at least these additional, separately validated
layers:

1. convert the primordial scalar, entropy, cross, and (when implemented)
   tensor spectra into the convention required by a Boltzmann solver;
2. propagate them through the thermal and recombination history;
3. calculate lensed or unlensed CMB angular spectra as required;
4. evaluate the named released likelihood, including nuisance parameters and
   priors;
5. sample or optimize model parameters and preserve the complete inference
   configuration.

Until those steps exist, the observational metadata is provenance for a
benchmark comparison, not evidence of agreement with or preference by the
Planck likelihood.
