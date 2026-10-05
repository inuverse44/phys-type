# Step 5: perturbations for any number of fields

Status: implemented in steps 5a–5f, except the sum-separable $N=3$ $\delta N$ comparison. It is step 5 of
[`architecture.md`](architecture.md).

## Goal

The theory layer states theories with any number of fields $N$, and the
background is already integrated for any $N$. The perturbation pipeline, from
the mode solver to the observables, is fixed to $N=2$. Step 5 removes that
restriction so that the same engine predicts

- single-field theories, $N=1$, beginning with Higgs inflation in the metric
  and Palatini formulations;
- the existing two-field results, unchanged;
- three or more fields, such as N-flation-like models and several
  nonminimally coupled fields.

## What is already general, and what is not

| Layer | Status | Items |
| --- | --- | --- |
| geometry, potentials, theory layer | any $N$ | `FieldSpaceGeometry<N>` and its extensions, `JetGeometry`, `MultifieldPotential<N>`, `GeneralizedNonminimal<N>` |
| background | any $N$ | `CovariantMultifieldSystem<P, G, N>`, `MultifieldState<N>` |
| cubic action | any $N$ | `canonical_cubic_{a,b,c}_coefficients`, the Hamiltonian vertices, the subhorizon initial vertices and their rates, `three_point_source_from_kernels`, the phase-space covariance types |
| gauge coefficients | any $N$ | `curvature_gauge_coefficients` |
| kinematics | $N=2$ | `TwoFieldKinematics`, `KinematicBasis2` (tangent $e_\sigma$, normal $e_s$, signed $\eta_\perp$) |
| mode solver | $N=2$ | `perturbation.rs`: `TwoFieldModes`, `ComplexMatrix2`, the cosmic- and conformal-time systems, Bunch–Davies initialization |
| projections | $N=2$ | `spectra.rs`: `mode_covariance`, `project` onto $(\mathcal R,\mathcal S)$ |
| transport systems | $N=2$ | `TwoPointTransportSystem`, the three-point systems and `CanonicalThreePointTransportState`, `BackgroundNeighborhood<2>`, `TwoFieldBackgroundProvider` |
| observables and runs | $N=2$ | `curvature_correlators`, `survey`, `run_two_field_spectra`, `mode_spectrum`, `run_bispectrum`, `three_point_state`, the record types and CSV writers |

The general layers were written with `const N: usize` from the start, so the
work is concentrated in the time evolution and the projections.

## Design decisions

**1. Const generics, not heap matrices.** Every pipeline type takes
`const N: usize`, like the layers that are already general. Field-space
indices keep their type-level variance and dimension, no allocation enters
the inner loops, and $N=1,2,3$ are separate monomorphizations. Stable Rust
cannot form `[T; N - 1]`, so nothing may be sized by the number of entropy
directions (see decision 3).

**2. Names lose "TwoField", in one break.** `TwoFieldModes` becomes
`Modes<N>`, `ConformalTwoFieldPerturbationSystem` becomes
`ConformalPerturbationSystem<P, G, N>`, `run_two_field_spectra` becomes
`run_spectra`, and so on. The crate is pre-1.0 with no outside users, so the
renaming is a clean break released as 0.3.0, without deprecated aliases.

**3. The adiabatic direction is a vector; the entropy sector is a
projector.** (The alternative is the Frenet basis: $e_1=e_\sigma$, $e_2$ along
the turn $D_te_1$, $e_3$ along the part of $D_te_2$ orthogonal to $e_1$ and $e_2$,
and so on, each with its own turn rate, as in the Frenet–Serret frame of a
space curve. It names $N-1$ entropy directions but is undefined wherever a
turn rate vanishes, for example on straight segments.) For any $N$ the tangent $e_\sigma^I=\dot\phi^I/\dot\sigma$ is unique.
The entropy sector is the orthogonal complement, of dimension $N-1$, and has
no preferred basis for $N\ge3$ (a Frenet basis degenerates on straight
segments). The kinematics therefore return

- the tangent $e_\sigma$ and the entropy projector
  $\perp^{IJ}=G^{IJ}-e_\sigma^Ie_\sigma^J$, both $N$-index objects;
- the turn vector $\omega^I=D_te_\sigma^I/H$ and its norm $|\eta_\perp|$;
- for $N=2$ only, the oriented normal $e_s$ and the signed $\eta_\perp$, as now.

**4. Observables are defined without an entropy basis.** With the field
power spectrum $\mathcal P^{IJ}(k)$, the tangent $e_\sigma$, the projector
$\perp^{IJ}=G^{IJ}-e_\sigma^Ie_\sigma^J$, and the normalization $H/\dot\sigma$,

$$
\mathcal P_{\mathcal R}=\Big(\tfrac{H}{\dot\sigma}\Big)^2e_{\sigma I}e_{\sigma J}\mathcal P^{IJ},\qquad
S^{IJ}=\Big(\tfrac{H}{\dot\sigma}\Big)^2\perp^I{}_K\perp^J{}_L\mathcal P^{KL},\qquad
C^I=\Big(\tfrac{H}{\dot\sigma}\Big)^2\perp^I{}_Ke_{\sigma L}\,\mathrm{Re}\,\mathcal P^{KL}.
$$

$S^{IJ}$ is the entropy covariance, an $N\times N$ tensor of rank at most $N-1$
(it vanishes along $e_\sigma$). The reported quantities are invariant under any
change of field coordinates and need no entropy basis:

| quantity | definition | meaning |
| --- | --- | --- |
| $\mathcal P_{\mathcal R}$ | above | curvature power |
| $\mathcal P_{\mathcal S}$ | $S^I{}_I$ | total entropy power |
| $\mathcal P_{\mathcal S,1}\ge\dots\ge\mathcal P_{\mathcal S,N-1}$ | eigenvalues of $S^I{}_J$, dropping the zero along $e_\sigma$ | principal entropy spectra, along axes set by the perturbations themselves |
| $\lvert\mathcal C_{\mathcal{RS}}\rvert$ | $\sqrt{G_{IJ}C^IC^J}$ | strength of the curvature–entropy correlation |
| $r_C$ | $\lvert\mathcal C_{\mathcal{RS}}\rvert/\sqrt{\mathcal P_{\mathcal R}\mathcal P_{\mathcal S}}$ | correlation coefficient, in $[0,1]$ |

The eigenvalues are computed from the $N\times N$ matrix, so nothing is sized
$N-1$; unlike a Frenet basis, the principal axes exist on straight segments
too. The CSV columns are

| $N$ | columns |
| --- | --- |
| 1 | `P_R` |
| 2 | `P_R, P_S, C_RS, r_C`, as now, with the signed `C_RS` |
| $\ge3$ | `P_R, P_S, P_S_1, …, P_S_{N-1}, C_RS_norm, r_C` |

For $N=2$ the single eigenvalue equals the trace and $\lvert\mathcal C_{\mathcal{RS}}\rvert$
is the magnitude of today's signed $\mathcal C_{\mathcal{RS}}$, so the general
definitions reproduce the current outputs, which is itself a test. Within the
inflation-only scope these entropy quantities are diagnostics of how much
non-adiabatic power remains, that is, whether $\zeta$ can still evolve;
observable isocurvature of a matter species needs reheating.

**5. Tensor modes join the observables.** Single-field validation needs the
tensor-to-scalar ratio, which the engine does not compute today. A tensor
mode obeys $v_t''+(k^2-a''/a)v_t=0$ in the Einstein frame, independent of the
fields, and $h_{ij}$ is unchanged by the conformal map to the Jordan frame.
It is evolved with the same conformal-time integrator, from Bunch–Davies to
the output slice, giving $\mathcal P_t$ and $r=\mathcal P_t/\mathcal P_{\mathcal R}$.

**6. Bit-for-bit two-field results are not promised; $10^{-12}$ is.** Loops
over $N$ reorder floating-point sums, so a general implementation can change
the last digits of today's two-field numbers. Each step must reproduce the
golden data to a relative $10^{-12}$; the golden data are refrozen once, at
the end, with the measured differences recorded in the pull request.

**7. The single-field background stays.** `ScalarFieldSystem` remains the
light single-field diagnostic of the MVP; single-field predictions go through
`CovariantMultifieldSystem<P, G, 1>` like every other $N$.

## Validation

| Target | Check |
| --- | --- |
| $N=2$ | golden data to $10^{-12}$; every existing test, including the ignored validations in both crates |
| $N=1$, analytic | quadratic and quartic single-field inflation: $\mathcal P_{\mathcal R}$ against the slow-roll amplitude to $O(\epsilon)$, $n_s$ and $r$ against their slow-roll values, $f_{\mathrm{NL}}$ against Maldacena's result and the squeezed consistency relation (already validated as a two-field limit) |
| $N=1$, Higgs inflation | $\Omega^2=1+\xi\phi^2$, $V_J=\lambda\phi^4/4$ through the theory layer, metric formulation: at large $\xi$, $n_s\simeq1-2/N_*$ and $r\simeq12/N_*^2$, with $N_*$ in the frame the reference uses ([Bezrukov and Shaposhnikov, arXiv:0710.3755](https://arxiv.org/abs/0710.3755)); Palatini: $r\simeq2/(\xi N_*^2)$ ([Bauer and Demir, arXiv:0803.2664](https://arxiv.org/abs/0803.2664)); the exact numbers differ from these at $O(1/N_*)$ and the differences are recorded, not tuned away; with a reproduction figure of $n_s$ and $r$ against $\xi$ |
| $N=1$ embedded | a single field run as $N=2$ with a decoupled heavy spectator reproduces the $N=1$ run |
| $N=3$, separable | (a case with a known answer) three uncoupled quadratic fields: $\mathcal P_{\mathcal R}$ and $f_{\mathrm{NL}}$ against the sum-separable $\delta N$ results of [Battefeld and Easther, arXiv:astro-ph/0610296](https://arxiv.org/abs/astro-ph/0610296) and a numerical $\delta N$ |
| $N=3$, covariance | (the answer must not depend on coordinates) a constant rotation of the three fields and a nonlinear field redefinition (generic three-field metric from `JetGeometry`) leave the projected spectra and $f_{\mathrm{NL}}$ unchanged, and the phase-space correlators transform as tensors, extending the existing two-field covariance tests |
| $N=3$ embedding | (an extra decoupled field must change nothing) the assisted model plus a decoupled heavy third field reproduces the two-field assisted results |

## Migration

One pull request per step; the two-field golden data are checked to
$10^{-12}$ at every step.

| Step | Content |
| --- | --- |
| 5a | kinematics for any $N$ (decision 3): `Kinematics<N>` with tangent, entropy projector, and turn vector; the oriented two-field basis becomes `oriented_kinematics` |
| 5b | the mode solver and projections for any $N$ (decisions 1, 2, 4): `PerturbationSystem<P, G, N>`, `ConformalPerturbationSystem<P, G, N>`, `Modes<N>`, a Bunch–Davies frame from Gram–Schmidt (the two-field normal for $N=2$, so the two-field numbers are unchanged bit for bit), and `invariant_spectra`; checked against the oriented two-field projection, rotation invariance for $N=3$, and a one-field run against two fields with a decoupled spectator. The analytic single-field checks need the run API and move to 5e |
| 5c | tensor modes and $r$ (decision 5): `CovariantMultifieldSystem::tensor_spectrum`, evolved on the surveyed trajectory with the scalar integrator; single-field quadratic inflation reproduces $\mathcal P_t=(2H_*^2/\pi^2)[1-2(C+1)\epsilon_1]$ and $r=16\epsilon_1(1+C\epsilon_2)$ (Stewart–Lyth, $C=-2+\ln2+\gamma_E$) to $2\epsilon_1^2$; the vacuum must start well inside the horizon (two e-folds before exit leaves a 1.3% error, six leave $4\times10^{-4}$). $r$ enters the records in 5e |
| 5d | two- and three-point transport for any $N$: `BackgroundProvider<N>`, `BackgroundNeighborhood<N>`, the transport systems and `CanonicalThreePointTransportState<N>`, `curvature_correlators`; two-field numbers unchanged bit for bit; the three-field two-point transport matches the evolution of three-field modes on flat and curved spaces. The three-point checks for $N\ne2$ need a surveyed trajectory and move to 5e |
| 5e | the run APIs for any $N$: `survey` (with `BackgroundSample<N>`), `BackgroundTrajectory<N>` (defaulting to two fields), `mode_spectrum` (returning `ModeSpectrum<N>` with `InvariantSpectra<N>`), `tensor_spectrum`, `three_point_state`, and `run_bispectrum`. `run_two_field_spectra` stays the two-field run with its oriented records and CSV files (the golden data). A one-field run equals a two-field run with a decoupled spectator ($\mathcal P_{\mathcal R}$, $\mathcal P_t$ to $10^{-9}$, $f_{\mathrm{NL}}$ to $10^{-6}$), and a three-field run equals the two-field run it embeds ($10^{-8}$, $10^{-6}$) |
| 5f | Higgs inflation in `phys-type-models` (metric and Palatini) with its figure, [`higgs-inflation.md`](higgs-inflation.md). The $N=3$ checks (rotation invariance, embedding of a two-field run) landed with 5b and 5e; the sum-separable $\delta N$ comparison remains |

## Open questions for review

None remain; see below.

## Resolved in review

- **Entropy reporting for $N\ge3$.** The trace, the principal entropy spectra
  (eigenvalues of the entropy covariance), the norm of the cross vector, and
  $r_C$ (decision 4); a Frenet basis is not used.
- **Deprecation.** A clean break at 0.3.0 (decision 2).
- **Frames and $n_s$ for Higgs inflation.** The physics does not depend on
  the frame. For a single field, $\zeta$ is frame invariant, and $\Omega\to1$
  after inflation, so $n_s(k)$ at a physical wavenumber $k$ is the same in both
  frames. What depends on the frame is the label "the pivot left the horizon
  $N_*$ e-folds before the end": the e-fold count between the pivot and the
  end differs by the change of $\ln\Omega$, with $\xi\phi_*^2\simeq4N_*/3\simeq80$
  and $\xi\phi_e^2\simeq1.15$ about $\tfrac12\ln(81/2.15)\simeq1.8$. "$N_*=60$" in the
  Einstein frame and in the Jordan frame therefore name different physical
  scales, whose $n_s$ differ by about $2\Delta N/N_*^2\simeq10^{-3}$. Predictions
  are reported at a physical $k$, with $N_*$ given in both frames, and a
  comparison with a paper first establishes which frame its $N_*$ uses.
