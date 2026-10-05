# Higgs inflation

The first single-field target of the $N$-field pipeline
([`n-field-perturbations.md`](n-field-perturbations.md)), stated through the
theory layer: $\Omega^2=1+\xi\phi^2$, a canonical Jordan-frame kinetic term, and
$V_J=\lambda\phi^4/4$ ([Bezrukov and Shaposhnikov,
arXiv:0710.3755](https://arxiv.org/abs/0710.3755)), in the metric and the
Palatini formulation ([Bauer and Demir,
arXiv:0803.2664](https://arxiv.org/abs/0803.2664)). The model is
`phys_type_models::higgs::HiggsTheory`; `phys_type::einstein_frame` gives
the one-field Einstein-frame system, and nothing in the engine is specific to
it.

## What is computed

The pivot leaves the horizon 60 Einstein-frame e-folds before $\epsilon_H=1$;
the spectra are read one e-fold before the end. $n_s$ is the logarithmic
derivative of the exact $\mathcal P_{\mathcal R}$ over $k e^{\pm0.02}$, and
$r=\mathcal P_t/\mathcal P_{\mathcal R}$ with the tensor mode of the same
wavenumber. $\lambda$ only sets the amplitude. The comparison values are
first-order slow roll from the same Einstein-frame potential and metric,
$\epsilon_V=\tfrac12(\partial_\chi\ln V)^2$ and $\eta_V=V_{,\chi\chi}/V$ with
$d\chi=\sqrt G\,d\phi$, and the large-$\xi$ limits.

![Higgs inflation](figures/higgs_inflation.png)

[`figures/higgs_inflation.png`](figures/higgs_inflation.png)
([PDF](figures/higgs_inflation.pdf)) is generated, in a few seconds, by

```sh
cargo run --release --example higgs_inflation
uv run plots/higgs_inflation.py diagnostics
```

## Results

At $\xi=10^4$:

| | exact | first-order slow roll | large-$\xi$ limit |
| --- | --- | --- | --- |
| metric $n_s$ | 0.96775 | 0.96783 | $1-2/N=0.96667$ |
| metric $r$ | $3.007\times10^{-3}$ | $2.964\times10^{-3}$ | $12/N^2=3.33\times10^{-3}$ |
| Palatini $n_s$ | 0.96643 | 0.96667 | $1-2/N=0.96667$ |
| Palatini $r$ | $5.70\times10^{-8}$ | $5.55\times10^{-8}$ | $2/(\xi N^2)=5.56\times10^{-8}$ |

- The exact $n_s$ follows first-order slow roll to $2.5\times10^{-4}$ for
  $\xi\ge10^{-2}$ and to $9\times10^{-4}$ near quartic inflation ($\xi=10^{-3}$,
  where $\epsilon$ is larger); $r$ follows it to 3%. These are the sizes of
  second-order slow-roll terms.
- The metric $r$ sits about 10% below $12/N^2$ at $N=60$, an $O(1/N)$
  correction of the large-$\xi$ limit (effectively a shift of $N$ by a few);
  $n_s$ is $1.1\times10^{-3}$ above $1-2/N$ for the same reason.
- In the Palatini formulation $\xi r$ is constant to 1% between $\xi=10^3$ and
  $10^4$, as $r\propto1/\xi$ requires.
- **E-folds in the two frames.** The same interval, 60 e-folds in the Einstein
  frame, spans $N_J-N_E\simeq2$ more Jordan-frame e-folds in the metric
  formulation at large $\xi$ and up to about 7 more in the Palatini
  formulation at $\xi=10^4$, where $\Omega$ changes far more between the pivot
  and the end. $n_s$ at a physical scale is frame invariant
  ([`frames.md`](frames.md)), but "$N_*=60$" in a paper must be read in the
  frame it was counted in; for Palatini Higgs inflation the choice moves $n_s$
  by about $2\Delta N/N^2\simeq4\times10^{-3}$.

## Tests

`crates/phys-type-models/tests/higgs.rs` (ignored, release mode): the exact
$n_s$ and $r$ against first-order slow roll for $\xi=10^{-2}$, 1, and $10^3$ in
both formulations, and the large-$\xi$ limits. The slow-roll functions are
checked against quartic inflation at $\xi=0$ and the large-$\xi$ limits in the
crate's unit tests.
