# Architecture: core, theories, and models

This document fixes the boundary between the general machinery of
`phys-type` and the particular models it is validated on. It was written
after the nonminimally assisted inflation work showed that model-specific
code had grown inside the core. Each migration step below is a separate pull
request.

## Philosophy

The repository has two cores.

1. **Physics in the types.** Mass dimensions, index variance, field-space
   tensors, and their derivatives are types, so that dimensionally or
   covariantly invalid expressions do not compile.
2. **From an action to predictions without approximation.** Given the
   action of a theory, such as a generalized nonminimal scalar-tensor theory,
   the background, the linear perturbations, and the tree-level three-point
   function are evolved exactly, and observables such as $P_\zeta$, $n_s$, and
   $f_{\mathrm{NL}}$ are returned with their definitions stated.

Models are users of these cores: single-field Higgs inflation, multi-field
nonminimally assisted inflation, constant-turn benchmarks, and later others.
No core module may know about a particular model or paper.

## Audit of the current boundary

At merge commit `6d6a7dc`:

| Item | Where | Problem |
| --- | --- | --- |
| No action layer | — | Users hand-write the Einstein-frame metric $G_{IJ}$ and potential $V$; nothing derives them from an action. |
| `NonminimalAssistedFieldSpace` | `src/nonminimal.rs` | One model, $f(s)=1+\xi_2s^2$ with a diagonal $(\phi,s)$ metric and a Gaussian-curvature formula valid only for that form, sits beside the generic geometries. |
| `NonminimalLoopInflationPotential`, `…PowerLaw…`, `…Hybrid…` | `src/nonminimal.rs`, root exports | Example potentials of one paper are core exports. |
| `assisted_slow_roll` | `src/assisted_slow_roll.rs`, a root module | One paper's slow-roll prescription is a top-level public module. |
| Model validation tests | `src/spectrum_run.rs` test module | About a thousand lines of assisted and constant-turn validation live in a core file. |
| Two-field perturbations | `perturbation`, `transport`, `spectrum_run`, `bispectrum_run` | The background is $N$-field, but the perturbation, transport, and observable pipeline is fixed to $N=2$, so single-field Higgs inflation and three or more fields cannot use it. |

The type layer (`quantity`, `dimension`, `autodiff`, `field_space`,
`potential`) is already model-independent.

## Target layers

```
phys-type (core crate) — knows no model
├─ type layer:        quantity, dimension, autodiff, field_space, potential
├─ theory layer:      theory classes and their frame maps          (new)
│    GeneralizedNonminimal<N>:  Omega^2(phi) R, k_IJ(phi), V_J(phi),
│                               metric or Palatini formulation
│    (later) Horndeski, DHOST, ... each with an admissibility statement
├─ engine layer:      exact evolution in a computational representation
│    EinsteinSigmaModel<N>: background, linear modes, two- and three-point
│                           transport for G_IJ(phi), V(phi), canonical gravity
└─ observable layer:  zeta projection, spectra, bispectra, f_NL, with the
                      frame and slice of every quantity recorded

phys-type-models (models crate) — depends only on the public core API
├─ higgs, assisted (loop, power-law, hybrid), constant_turn, ...
├─ references/: the formulas of particular papers (e.g. HKKPT Eqs. (3.1)-(3.9))
├─ tests/:      validation against the literature
└─ examples/:   the data behind the tracked figures
```

The dependency runs one way, from `phys-type-models` to `phys-type`. A
reference from the core to a model is a compile error, which is the same kind
of guarantee the type layer gives for dimensions.

## Frames are a representation, not a principle

The core does not require every theory to be computed in the Einstein frame.
It separates three things.

- **The theory is stated in the user's language**, for example the Jordan
  frame of a nonminimal coupling. A theory class says what it is and which
  computational representations it admits.
- **The computation happens in a representation the engine supports.** At
  present the only validated engine is the Einstein-frame sigma model with
  canonical gravity, because the cubic transport was derived and tested for
  it. For a generalized nonminimal theory the map is exact:
  $G_{IJ}=k_{IJ}/\Omega^2+\kappa\,\tfrac32M_{\mathrm{pl}}^2\,\partial_I\Omega^2\,\partial_J\Omega^2/\Omega^4$
  and $V=V_J/\Omega^4$, with $\kappa=1$ for the metric and $\kappa=0$ for the
  Palatini formulation.
- **Observables carry their frame and slice.** The theory layer returns the
  correspondence along the trajectory: the conformal factor, the e-fold
  numbers $N_J=N_E+\ln\Omega$, the Jordan-frame Hubble rate, and the field
  map. A pivot "60 e-folds before the end" or an end-of-inflation condition
  must name its frame. $\zeta$ is invariant under a field-dependent conformal
  transformation once the evolution is adiabatic; before that, and whenever
  matter couples to the Jordan metric, the definition used is stated
  explicitly. These correspondences are written up in `docs/frames.md` with
  references when the theory layer lands.

Horndeski and DHOST theories may later need their own engine. They plug in
at the same boundary: a theory class, an engine that accepts it, and
observables with stated definitions.

## Can admissibility of the Einstein frame be tested?

Partly, at three levels.

1. **By type, at compile time.** Only theory classes that implement the map
   to `EinsteinSigmaModel` can be handed to that engine. For
   `GeneralizedNonminimal<N>` the map exists by construction. Whether a class
   admits the map is a property of the class, decided once with a derivation
   and references, not per run.
2. **Along the trajectory, at run time.** The map has validity conditions
   that depend on the fields: $\Omega^2>0$ everywhere visited, and a positive
   definite Einstein-frame metric $G_{IJ}$ (no ghost). These are checked on
   every step of the background, and the run fails with the violated
   condition and its location.
3. **For Horndeski and DHOST, a supplied transformation is verifiable, its
   existence is a classification question.** Given a candidate conformal and
   disformal transformation $\tilde g_{\mu\nu}=C(\phi,X)g_{\mu\nu}+D(\phi,X)\partial_\mu\phi\partial_\nu\phi$,
   checking it is mechanical: transform the action and test that the result
   is in the target class, with the transformation invertible on the
   trajectory. Whether such a transformation exists is answered by the
   classification of degenerate theories: disformal transformations map
   Horndeski theories outside the Horndeski class
   ([Zumalacárregui and García-Bellido, arXiv:1308.4685](https://arxiv.org/abs/1308.4685)),
   the degenerate theories are classified by
   [Langlois and Noui, arXiv:1510.06930](https://arxiv.org/abs/1510.06930),
   and the classes are related by conformal and disformal transformations in
   [Ben Achour, Langlois, and Noui, arXiv:1602.08398](https://arxiv.org/abs/1602.08398).
   Before a class is admitted, its admissibility conditions will be derived
   and recorded with these references; a theory outside every admitted class
   needs an engine of its own.

## Public API the models need

The assisted validation currently calls crate-private functions of
`spectrum_run`, such as the background survey and the single-mode evolution.
Moving the models out forces a decision for each: either it becomes a
documented public API of the engine, or the test that needs it stays in the
core as an engine test. The intent is a small public surface: a background
survey returning a queryable trajectory, single-mode and spectrum runs on it,
and bispectrum runs, all with explicit settings.

## Migration plan

Each step is one pull request. A refactoring step must leave every number
unchanged, which the golden data of step 1 checks bit for bit.

| Step | Pull request | Content | Safety net |
| --- | --- | --- | --- |
| 0 | design | this document | review |
| 1 | golden data | freeze the CSV output of the four examples and test that it is reproduced | — |
| 2 | workspace | move the crate to `crates/phys-type`, no code change | CI, golden data |
| 3 | theory layer | `GeneralizedNonminimal<N>` and its exact map, the frame correspondence, run-time admissibility checks, `docs/frames.md` | the generic map reproduces `NonminimalAssistedFieldSpace` to round-off |
| 4 | models crate | assisted as an instance of step 3; papers' formulas, example potentials, validation tests, and examples move to `crates/phys-type-models` | golden data |
| 5+ | $N$ fields | the perturbation, transport, and observable pipeline for any $N$, starting with $N=1$ and Higgs inflation | the two-field results unchanged |

Expensive validations remain ignored tests run in release mode; CI checks the
golden data of the fast examples with `scripts/golden.py` (see
`CONTRIBUTING.md`), and the slow examples are checked locally with `--all`.

## Open questions

- Which frame is the default for reported e-fold numbers when a theory is
  given in the Jordan frame.
- How the Palatini formulation enters the theory class beyond the $\kappa=0$
  metric, for example with torsion or non-metricity in more general actions.
- Whether the three-point engine should also accept non-canonical kinetic
  terms $P(X^{IJ},\phi)$ before Horndeski, as an intermediate class.
