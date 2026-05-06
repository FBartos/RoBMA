README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/FBartos/RoBMA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FBartos/RoBMA/actions/workflows/R-CMD-check.yaml)
[![R-unit-tests](https://github.com/FBartos/RoBMA/actions/workflows/R-CMD-tests.yaml/badge.svg)](https://github.com/FBartos/RoBMA/actions/workflows/R-CMD-tests.yaml)
[![test-coverage](https://github.com/FBartos/RoBMA/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/FBartos/RoBMA/actions/workflows/test-coverage.yaml)
[![Codecov](https://codecov.io/gh/FBartos/RoBMA/graph/badge.svg)](https://app.codecov.io/gh/FBartos/RoBMA)
[![CRAN
status](https://www.r-pkg.org/badges/version/RoBMA)](https://CRAN.R-project.org/package=RoBMA)
<!-- badges: end -->

# Robust Bayesian Meta-Analysis (RoBMA)

The `RoBMA` R package fits a comprehensive collection of Bayesian
meta-analytic models (Bartoš & Maier, 2020). It provides single-model
fits, meta-regression, location-scale models, multilevel models, fully
Bayesian model-averaged ensembles, and built-in publication-bias
adjustment using selection models, PET-PEESE, and Robust Bayesian
Meta-Analysis (RoBMA). Posterior sampling is performed by JAGS through a
custom C++ module that implements the weighted distributions used in
selection models. Prior handling, plotting, and Bayesian diagnostics are
delegated to the companion `BayesTools` package.

This introduction maps the package at a high level and points to the
more detailed vignettes that follow.

## System Requirements

The `RoBMA` R package requires JAGS 4.3.1 or newer (Plummer, 2003).
Install JAGS from <https://mcmc-jags.sourceforge.io/> before installing
the R package.

``` r
install.packages("RoBMA")
```

**Note on backwards compatibility.** The 4.0 release is not backwards
compatible with earlier versions of the package. The *Manuscript
Companions* section accompanies the methodological papers in which the
`RoBMA` R package was originally developed and illustrated; the
vignettes there have been updated to the 4.0 version of the package, and
the numerical results differ from the published values because the
underlying algorithms have changed. To reproduce the published analyses
exactly, install `RoBMA` version 3.6.1 together with `BayesTools`
version 0.2.23:

``` r
remotes::install_version("BayesTools", version = "0.2.23")
remotes::install_version("RoBMA",      version = "3.6.1")
```

## Main Functions

The package exposes a small set of high-level fitting functions.

| Function | What it fits |
|:---|:---|
| [`brma()`](man/brma.Rd), [`brma.norm()`](man/brma.Rd) | Bayesian random-effects meta-analysis |
| [`brma.glmm()`](man/brma.glmm.Rd) | Bayesian GLMM meta-analysis (binomial / log OR or Poisson / log IRR) |
| [`bselmodel()`](man/bselmodel.Rd) | Bayesian weight-function selection model |
| [`bPET()`](man/bPET.Rd), [`bPEESE()`](man/bPEESE.Rd) | Bayesian PET / PEESE publication-bias adjustment |
| [`BMA()`](man/BMA.Rd), [`BMA.norm()`](man/BMA.Rd) | Bayesian model-averaging across presence / absence of effect and heterogeneity (no bias adjustment) |
| [`BMA.glmm()`](man/BMA.glmm.Rd) | Bayesian model-averaging for GLMM meta-analysis |
| [`RoBMA()`](man/RoBMA.Rd) | Robust Bayesian model-averaging including publication-bias models |

`brma()` and `RoBMA()` are the main workhorses. `brma()` fits a single
model in the same spirit as `metafor::rma()`; `RoBMA()` averages across
an ensemble that includes publication-bias adjustment. `BMA()` is
`RoBMA()` without the bias-adjustment models.

All fitting functions

- share an effect-size interface (`yi`, `sei` / `vi`, `measure`
  arguments; apart from glmm),
- can be extended into a meta-regression (`mods` argument),
- can be extended into a location-scale model (`scale` argument),
- can be extended into a multilevel model (`cluster` argument),
- and all produce objects that work with the same set of inference
  helpers (`summary()`, `plot()`, `predict()`, `loo()`, `funnel()`,
  `regplot()`, residuals and influence diagnostics, …).

## Vignette Map

The vignettes are organized into four sections.

### Foundations

| Vignette | Topic |
|:---|:---|
| [*Introduction to RoBMA*](vignettes/v00-introduction.Rmd) | Package overview, fitting functions, and vignette map |
| [*Prior Distributions*](vignettes/v01-prior-distributions.Rmd) | Default, informed, and custom prior distributions; rescaling guidance |
| [*Bayesian Meta-Analysis*](vignettes/v02-bayesian-meta-analysis.Rmd) | Random-effects meta-analysis with `brma()`, compared with `metafor::rma()` |
| [*Feature Coverage*](vignettes/v03-feature-coverage.Rmd) | Overview of available functionality across model families |

The foundations section gives the package overview, prior specification,
a baseline `brma()` workflow, and current feature coverage. The
[*Bayesian Meta-Analysis*](vignettes/v02-bayesian-meta-analysis.Rmd)
vignette walks through the BCG-vaccine example from the `metafor`
package in `brma()`, covering summaries, meta-regression, marginal
means, residuals, influence, LOO, and the standard meta-analytic plots.

### Correspondence with `metafor`

| Vignette | Topic |
|:---|:---|
| [*Multilevel Meta-Analysis*](vignettes/v10-metafor-parity-multilevel.Rmd) | `brma()` with `cluster` and `metafor::rma.mv()` for 3-level data |
| [*Publication-Bias Adjustment*](vignettes/v11-metafor-parity-publication-bias.Rmd) | `bselmodel()` / `bPET()` / `bPEESE()` and `metafor::selmodel()` |
| [*Location-Scale Models*](vignettes/v12-metafor-parity-location-scale.Rmd) | `brma()` with `scale` and `metafor::rma.ls()` |
| [*Generalized Linear Mixed-Effects Meta-Analysis*](vignettes/v13-metafor-parity-glmm.Rmd) | `brma.glmm()` and `metafor::rma.glmm()` for binomial and Poisson outcomes |

Each vignette starts from a `metafor` package analysis and shows the
matching `RoBMA` R package syntax, output, and diagnostics.

### Bayesian Model Averaging

| Vignette | Topic |
|:---|:---|
| [*Bayesian Model Averaging*](vignettes/v20-bayesian-model-averaging.Rmd) | Accounting for model uncertainty across presence and absence of effect and heterogeneity, with posterior model probabilities and inclusion Bayes factors |
| [*Robust Bayesian Meta-Analysis*](vignettes/v21-robust-bayesian-meta-analysis.Rmd) | Extending the ensemble to publication-bias models used in `RoBMA()`, with `PSMA`, `PP`, and bespoke ensemble specifications |

The first vignette introduces model averaging on an ensemble of models
with and without effect and heterogeneity. The second extends the
ensemble to publication-bias models used in `RoBMA()`.

### Manuscript Companions

| Vignette | Companion paper / dataset |
|:---|:---|
| [*Adjusting for Publication Bias Tutorial*](vignettes/v30-tutorial.Rmd) | Bartoš et al. (2022), JASP / R tutorial on `Lui2015` |
| [*Robust Bayesian Meta-Regression*](vignettes/v31-robma-metaregression.Rmd) | Bartoš et al. (2025), RoBMA-reg on `Andrews2021` |
| [*Multilevel Robust Bayesian Meta-Analysis*](vignettes/v32-robma-multilevel.Rmd) | Bartoš et al. (2026), Multilevel RoBMA on `Johnides2025` |
| [*Multilevel Robust Bayesian Meta-Regression*](vignettes/v33-robma-multilevel-metaregression.Rmd) | Multilevel RoBMA-reg on `Kroupova2021` |
| [*Informed Bayesian Meta-Analysis in Medicine*](vignettes/v34-bma-norm-medicine.Rmd) | Bartoš et al. (2021), informed prior distributions for medical meta-analysis (continuous outcomes) |
| [*Informed Bayesian Meta-Analysis with Binary Outcomes*](vignettes/v35-bma-glmm-medicine.Rmd) | Bartoš et al. (2023), informed prior distributions for medical meta-analysis (binary and time-to-event outcomes) |
| [*Zplot Publication-Bias Diagnostics*](vignettes/v36-zplot.Rmd) | Bartoš & Schimmack (2025), zplot diagnostics on `Hoppen2025` |

These vignettes reproduce or update analyses from published papers and
serve as references when citing the corresponding methodological work.

## Where to Start

- For readers familiar with the `metafor` package: read [*Bayesian
  Meta-Analysis*](vignettes/v02-bayesian-meta-analysis.Rmd), then the
  *Correspondence with `metafor`* vignettes for the model family they
  use most.
- For readers new to Bayesian meta-analysis: read [*Prior
  Distributions*](vignettes/v01-prior-distributions.Rmd) and then
  [*Bayesian Meta-Analysis*](vignettes/v02-bayesian-meta-analysis.Rmd).
- For publication-bias adjustment: read [*Robust Bayesian
  Meta-Analysis*](vignettes/v21-robust-bayesian-meta-analysis.Rmd), with
  [*Adjusting for Publication Bias
  Tutorial*](vignettes/v30-tutorial.Rmd) as a worked example.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-bartos2021bayesian" class="csl-entry">

Bartoš, F., Gronau, Q. F., Timmers, B., Otte, W. M., Ly, A., &
Wagenmakers, E.-J. (2021). Bayesian model-averaged meta-analysis in
medicine. *Statistics in Medicine*, *40*(30), 6743–6761.
<https://doi.org/10.1002/sim.9170>

</div>

<div id="ref-RoBMA" class="csl-entry">

Bartoš, F., & Maier, M. (2020). *RoBMA: An R package for robust Bayesian
meta-analyses*. <https://CRAN.R-project.org/package=RoBMA>

</div>

<div id="ref-bartos2020adjusting" class="csl-entry">

Bartoš, F., Maier, M., Quintana, D. S., & Wagenmakers, E.-J. (2022).
Adjusting for publication bias in JASP and R — Selection models,
PET-PEESE, and robust Bayesian meta-analysis. *Advances in Methods and
Practices in Psychological Science*, *5*(3), 1–19.
<https://doi.org/10.1177/25152459221109259>

</div>

<div id="ref-bartos2023robust" class="csl-entry">

Bartoš, F., Maier, M., Stanley, T., & Wagenmakers, E.-J. (2025). Robust
Bayesian meta-regression: Model-averaged moderation analysis in the
presence of publication bias. *Psychological Methods*.
<https://doi.org/10.1037/met0000737>

</div>

<div id="ref-bartos2025robust" class="csl-entry">

Bartoš, F., Maier, M., & Wagenmakers, E.-J. (2026). Robust Bayesian
multilevel meta-analysis: Adjusting for publication bias in the presence
of dependent effect sizes. *Behavior Research Methods*.
<https://doi.org/10.31234/osf.io/9tgp2_v1>

</div>

<div id="ref-bartos2023empirical" class="csl-entry">

Bartoš, F., Otte, W. M., Gronau, Q. F., Timmers, B., Ly, A., &
Wagenmakers, E.-J. (2023). *Empirical prior distributions for Bayesian
meta-analyses of binary and time-to-event outcomes*.
<https://doi.org/10.48550/arXiv.2306.11468>

</div>

<div id="ref-bartos2025zcurve" class="csl-entry">

Bartoš, F., & Schimmack, U. (2025). Zplot: A visual diagnostic for
publication bias in meta-analysis. In *arXiv*.
<https://doi.org/10.48550/arXiv.2509.07171>

</div>

<div id="ref-plummer2003jags" class="csl-entry">

Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical
models using Gibbs sampling. *Proceedings of the 3rd International
Workshop on Distributed Statistical Computing (DSC 2003)*, 1–10.

</div>

</div>
