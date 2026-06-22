The `RoBMA` R package fits a comprehensive collection of Bayesian meta-analytic models [@RoBMA].
It provides single-model fits, meta-regression, location-scale models, multilevel models, fully Bayesian model-averaged ensembles, and built-in publication-bias adjustment using selection models, PET-PEESE, and Robust Bayesian Meta-Analysis (RoBMA).
Posterior sampling is performed by JAGS through a custom C++ module that implements the weighted distributions used in selection models.
Prior handling, plotting, and Bayesian diagnostics are delegated to the companion `BayesTools` package.

This introduction maps the package at a high level and points to the more detailed vignettes that follow.

## System Requirements

The `RoBMA` R package requires JAGS 4.3.1 or newer [@plummer2003jags].
Install JAGS from <https://mcmc-jags.sourceforge.io/> before installing the R package.

```r
install.packages("RoBMA")
```

**Note on backwards compatibility.**
The 4.0 release is not backwards compatible with earlier versions of the package.
The *Manuscript Companions* section accompanies the methodological papers in which the `RoBMA` R package was originally developed and illustrated; the vignettes there have been updated to the 4.0 version of the package, and the numerical results differ from the published values because the underlying algorithms have changed.
To reproduce the published analyses exactly, install `RoBMA` version 3.6.1 together with `BayesTools` version 0.2.23:

```r
remotes::install_version("BayesTools", version = "0.2.23")
remotes::install_version("RoBMA",      version = "3.6.1")
```

## Main Functions

The package exposes a small set of high-level fitting functions.

| Function              | What it fits                                                               |
|:----------------------|:---------------------------------------------------------------------------|
| [`brma()`](`r reference_url("brma")`), [`brma.norm()`](`r reference_url("brma")`) | Bayesian random-effects meta-analysis                                      |
| [`brma.glmm()`](`r reference_url("brma.glmm")`) | Bayesian GLMM meta-analysis (binomial / log OR or Poisson / log IRR)       |
| [`bselmodel()`](`r reference_url("bselmodel")`) | Bayesian weight-function selection model                                   |
| [`bPET()`](`r reference_url("bPET")`), [`bPEESE()`](`r reference_url("bPEESE")`) | Bayesian PET / PEESE publication-bias adjustment                           |
| [`BMA()`](`r reference_url("BMA")`), [`BMA.norm()`](`r reference_url("BMA")`) | Bayesian model-averaging across presence / absence of effect and heterogeneity (no bias adjustment) |
| [`BMA.glmm()`](`r reference_url("BMA.glmm")`) | Bayesian model-averaging for GLMM meta-analysis                            |
| [`RoBMA()`](`r reference_url("RoBMA")`) | Robust Bayesian model-averaging including publication-bias models          |

`brma()` and `RoBMA()` are the main workhorses.
`brma()` fits a single model in the same spirit as `metafor::rma()`; `RoBMA()` averages across an ensemble that includes publication-bias adjustment.
`BMA()` is `RoBMA()` without the bias-adjustment models.

All fitting functions

- share an effect-size interface (`yi`, `sei` / `vi`, `measure` arguments; apart from glmm),
- can be extended into a meta-regression (`mods` argument),
- can be extended into a location-scale model (`scale` argument),
- can be extended into a multilevel model (`cluster` argument),
- and all produce objects that work with the same set of inference helpers (`summary()`, `plot()`, `predict()`, `loo()`, `funnel()`, `regplot()`, residuals and influence diagnostics, ...).

## Vignette Map

The vignettes are organized into four sections.

### Foundations

| Vignette                  | Topic                                                                 |
|:--------------------------|:----------------------------------------------------------------------|
| [*Introduction to RoBMA*](`r article_url("v00-introduction")`) | Package overview, fitting functions, and vignette map |
| [*Prior Distributions*](`r article_url("v01-prior-distributions")`) | Default, informed, and custom prior distributions; rescaling guidance |
| [*Bayesian Meta-Analysis*](`r article_url("v02-bayesian-meta-analysis")`) | Random-effects meta-analysis with `brma()`, compared with `metafor::rma()` |
| [*Feature Coverage*](`r article_url("v03-feature-coverage")`) | Overview of available functionality across model families             |

The foundations section gives the package overview, prior specification, a baseline `brma()` workflow, and current feature coverage.
The [*Bayesian Meta-Analysis*](`r article_url("v02-bayesian-meta-analysis")`) vignette walks through the BCG-vaccine example from the `metafor` package in `brma()`, covering summaries, meta-regression, marginal means, residuals, influence, LOO, and the standard meta-analytic plots.

### Correspondence with `metafor`

| Vignette                                          | Topic                                                |
|:--------------------------------------------------|:-----------------------------------------------------|
| [*Multilevel Meta-Analysis*](`r article_url("v10-metafor-parity-multilevel")`) | `brma()` with `cluster` and `metafor::rma.mv()` for 3-level data |
| [*Publication-Bias Adjustment*](`r article_url("v11-metafor-parity-publication-bias")`) | `bselmodel()` / `bPET()` / `bPEESE()` and `metafor::selmodel()` |
| [*Location-Scale Models*](`r article_url("v12-metafor-parity-location-scale")`) | `brma()` with `scale` and `metafor::rma.ls()`          |
| [*Generalized Linear Mixed-Effects Meta-Analysis*](`r article_url("v13-metafor-parity-glmm")`) | `brma.glmm()` and `metafor::rma.glmm()` for binomial and Poisson outcomes |
| [*Multivariate and Multilevel Meta-Analysis*](`r article_url("v14-metafor-parity-multivariate")`) | `brma.mv()` with known sampling covariance matrices and structured random effects |

Each vignette starts from a `metafor` package analysis and shows the matching `RoBMA` R package syntax, output, and diagnostics.

### Bayesian Model Averaging

| Vignette                       | Topic                                                                  |
|:-------------------------------|:-----------------------------------------------------------------------|
| [*Bayesian Model Averaging*](`r article_url("v20-bayesian-model-averaging")`) | Accounting for model uncertainty across presence and absence of effect and heterogeneity, with posterior model probabilities and inclusion Bayes factors |
| [*Robust Bayesian Meta-Analysis*](`r article_url("v21-robust-bayesian-meta-analysis")`) | Extending the ensemble to publication-bias models used in `RoBMA()`, with `PSMA`, `PP`, and bespoke ensemble specifications |

The first vignette introduces model averaging on an ensemble of models with and without effect and heterogeneity.
The second extends the ensemble to publication-bias models used in `RoBMA()`.

### Manuscript Companions

| Vignette                                  | Companion paper / dataset                                          |
|:------------------------------------------|:-------------------------------------------------------------------|
| [*Adjusting for Publication Bias Tutorial*](`r article_url("v30-tutorial")`) | @bartos2020adjusting, JASP / R tutorial on `Lui2015`              |
| [*Robust Bayesian Meta-Regression*](`r article_url("v31-robma-metaregression")`) | @bartos2023robust, RoBMA-reg on `Andrews2021`                     |
| [*Multilevel Robust Bayesian Meta-Analysis*](`r article_url("v32-robma-multilevel")`) | @bartos2025robust, Multilevel RoBMA on `Johnides2025`             |
| [*Multilevel Robust Bayesian Meta-Regression*](`r article_url("v33-robma-multilevel-metaregression")`) | Multilevel RoBMA-reg on `Kroupova2021`                            |
| [*Informed Bayesian Meta-Analysis in Medicine*](`r article_url("v34-bma-norm-medicine")`) | @bartos2021bayesian, informed prior distributions for medical meta-analysis (continuous outcomes) |
| [*Informed Bayesian Meta-Analysis with Binary Outcomes*](`r article_url("v35-bma-glmm-medicine")`) | @bartos2023empirical, informed prior distributions for medical meta-analysis (binary and time-to-event outcomes) |
| [*Zplot Publication-Bias Diagnostics*](`r article_url("v36-zplot")`) | @bartos2025zcurve, zplot diagnostics on `Hoppen2025`              |

These vignettes reproduce or update analyses from published papers and serve as references when citing the corresponding methodological work.

## Where to Start

- For readers familiar with the `metafor` package: read [*Bayesian Meta-Analysis*](`r article_url("v02-bayesian-meta-analysis")`), then the *Correspondence with `metafor`* vignettes for the model family they use most.
- For readers new to Bayesian meta-analysis: read [*Prior Distributions*](`r article_url("v01-prior-distributions")`) and then [*Bayesian Meta-Analysis*](`r article_url("v02-bayesian-meta-analysis")`).
- For publication-bias adjustment: read [*Robust Bayesian Meta-Analysis*](`r article_url("v21-robust-bayesian-meta-analysis")`), with [*Adjusting for Publication Bias Tutorial*](`r article_url("v30-tutorial")`) as a worked example.

## References
