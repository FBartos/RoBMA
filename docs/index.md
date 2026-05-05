# RoBMA

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
| [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md), [`brma.norm()`](https://fbartos.github.io/RoBMA/reference/brma.md) | Bayesian random-effects meta-analysis |
| [`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md) | Bayesian GLMM meta-analysis (binomial / log OR or Poisson / log IRR) |
| [`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md) | Bayesian weight-function selection model |
| [`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md), [`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md) | Bayesian PET / PEESE publication-bias adjustment |
| [`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md), [`BMA.norm()`](https://fbartos.github.io/RoBMA/reference/BMA.md) | Bayesian model-averaging across presence / absence of effect and heterogeneity (no bias adjustment) |
| [`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md) | Bayesian model-averaging for GLMM meta-analysis |
| [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) | Robust Bayesian model-averaging including publication-bias models |

[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) and
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) are the
main workhorses.
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) fits a
single model in the same spirit as
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html);
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) averages
across an ensemble that includes publication-bias adjustment.
[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md) is
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) without
the bias-adjustment models.

All fitting functions

- share an effect-size interface (`yi`, `sei` / `vi`, `measure`
  arguments; apart from glmm),
- can be extended into a meta-regression (`mods` argument),
- can be extended into a location-scale model (`scale` argument),
- can be extended into a multilevel model (`cluster` argument),
- and all produce objects that work with the same set of inference
  helpers ([`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`loo()`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
  [`funnel()`](https://fbartos.github.io/RoBMA/reference/funnel.md),
  [`regplot()`](https://fbartos.github.io/RoBMA/reference/regplot.md),
  residuals and influence diagnostics, …).

## Vignette Map

The vignettes are organized into four sections.

### Foundations

| Vignette | Topic |
|:---|:---|
| [*Introduction to RoBMA*](https://fbartos.github.io/RoBMA/articles/00-introduction.md) | Package overview, fitting functions, and vignette map |
| [*Prior Distributions*](https://fbartos.github.io/RoBMA/articles/01-prior-distributions.md) | Default, informed, and custom prior distributions; rescaling guidance |
| [*Bayesian Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/02-bayesian-meta-analysis.md) | Random-effects meta-analysis with [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md), compared with [`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html) |
| [*Feature Coverage*](https://fbartos.github.io/RoBMA/articles/03-feature-coverage.md) | Overview of available functionality across model families |

The foundations section gives the package overview, prior specification,
a baseline [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md)
workflow, and current feature coverage. The [*Bayesian
Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/02-bayesian-meta-analysis.md)
vignette walks through the BCG-vaccine example from the `metafor`
package in
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md), covering
summaries, meta-regression, marginal means, residuals, influence, LOO,
and the standard meta-analytic plots.

### Correspondence with `metafor`

| Vignette | Topic |
|:---|:---|
| [*Multilevel Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/10-metafor-parity-multilevel.md) | [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) with `cluster` and [`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html) for 3-level data |
| [*Publication-Bias Adjustment*](https://fbartos.github.io/RoBMA/articles/11-metafor-parity-publication-bias.md) | [`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md) / [`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md) / [`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md) and [`metafor::selmodel()`](https://wviechtb.github.io/metafor/reference/selmodel.html) |
| [*Location-Scale Models*](https://fbartos.github.io/RoBMA/articles/12-metafor-parity-location-scale.md) | [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) with `scale` and `metafor::rma.ls()` |
| [*Generalized Linear Mixed-Effects Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/13-metafor-parity-glmm.md) | [`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md) and [`metafor::rma.glmm()`](https://wviechtb.github.io/metafor/reference/rma.glmm.html) for binomial and Poisson outcomes |

Each vignette starts from a `metafor` package analysis and shows the
matching `RoBMA` R package syntax, output, and diagnostics.

### Bayesian Model Averaging

| Vignette | Topic |
|:---|:---|
| [*Bayesian Model Averaging*](https://fbartos.github.io/RoBMA/articles/20-bayesian-model-averaging.md) | Accounting for model uncertainty across presence and absence of effect and heterogeneity, with posterior model probabilities and inclusion Bayes factors |
| [*Robust Bayesian Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/21-robust-bayesian-meta-analysis.md) | Extending the ensemble to publication-bias models used in [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md), with `PSMA`, `PP`, and bespoke ensemble specifications |

The first vignette introduces model averaging on an ensemble of models
with and without effect and heterogeneity. The second extends the
ensemble to publication-bias models used in
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md).

### Manuscript Companions

| Vignette | Companion paper / dataset |
|:---|:---|
| [*Adjusting for Publication Bias Tutorial*](https://fbartos.github.io/RoBMA/articles/30-tutorial.md) | Bartoš et al. (2022), JASP / R tutorial on `Lui2015` |
| [*Robust Bayesian Meta-Regression*](https://fbartos.github.io/RoBMA/articles/31-robma-metaregression.md) | Bartoš et al. (2025), RoBMA-reg on `Andrews2021` |
| [*Multilevel Robust Bayesian Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/32-robma-multilevel.md) | Bartoš et al. (2026), Multilevel RoBMA on `Johnides2025` |
| [*Multilevel Robust Bayesian Meta-Regression*](https://fbartos.github.io/RoBMA/articles/33-robma-multilevel-metaregression.md) | Multilevel RoBMA-reg on `Kroupova2021` |
| [*Informed Bayesian Meta-Analysis in Medicine*](https://fbartos.github.io/RoBMA/articles/34-bma-norm-medicine.md) | Bartoš et al. (2021), informed prior distributions for medical meta-analysis (continuous outcomes) |
| [*Informed Bayesian Meta-Analysis with Binary Outcomes*](https://fbartos.github.io/RoBMA/articles/35-bma-glmm-medicine.md) | Bartoš et al. (2023), informed prior distributions for medical meta-analysis (binary and time-to-event outcomes) |
| [*Zplot Publication-Bias Diagnostics*](https://fbartos.github.io/RoBMA/articles/36-zplot.md) | Bartoš & Schimmack (2025), zplot diagnostics on `Hoppen2025` |

These vignettes reproduce or update analyses from published papers and
serve as references when citing the corresponding methodological work.

## Where to Start

- For readers familiar with the `metafor` package: read [*Bayesian
  Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/02-bayesian-meta-analysis.md),
  then the *Correspondence with `metafor`* vignettes for the model
  family they use most.
- For readers new to Bayesian meta-analysis: read [*Prior
  Distributions*](https://fbartos.github.io/RoBMA/articles/01-prior-distributions.md)
  and then [*Bayesian
  Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/02-bayesian-meta-analysis.md).
- For publication-bias adjustment: read [*Robust Bayesian
  Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/21-robust-bayesian-meta-analysis.md),
  with [*Adjusting for Publication Bias
  Tutorial*](https://fbartos.github.io/RoBMA/articles/30-tutorial.md) as
  a worked example.

## References

Bartoš, F., Gronau, Q. F., Timmers, B., Otte, W. M., Ly, A., &
Wagenmakers, E.-J. (2021). Bayesian model-averaged meta-analysis in
medicine. *Statistics in Medicine*, *40*(30), 6743–6761.
<https://doi.org/10.1002/sim.9170>

Bartoš, F., & Maier, M. (2020). *RoBMA: An R package for robust Bayesian
meta-analyses*. <https://CRAN.R-project.org/package=RoBMA>

Bartoš, F., Maier, M., Quintana, D. S., & Wagenmakers, E.-J. (2022).
Adjusting for publication bias in JASP and R — Selection models,
PET-PEESE, and robust Bayesian meta-analysis. *Advances in Methods and
Practices in Psychological Science*, *5*(3), 1–19.
<https://doi.org/10.1177/25152459221109259>

Bartoš, F., Maier, M., Stanley, T., & Wagenmakers, E.-J. (2025). Robust
Bayesian meta-regression: Model-averaged moderation analysis in the
presence of publication bias. *Psychological Methods*.
<https://doi.org/10.1037/met0000737>

Bartoš, F., Maier, M., & Wagenmakers, E.-J. (2026). Robust Bayesian
multilevel meta-analysis: Adjusting for publication bias in the presence
of dependent effect sizes. *Behavior Research Methods*.
<https://doi.org/10.31234/osf.io/9tgp2_v1>

Bartoš, F., Otte, W. M., Gronau, Q. F., Timmers, B., Ly, A., &
Wagenmakers, E.-J. (2023). *Empirical prior distributions for Bayesian
meta-analyses of binary and time-to-event outcomes*.
<https://doi.org/10.48550/arXiv.2306.11468>

Bartoš, F., & Schimmack, U. (2025). Zplot: A visual diagnostic for
publication bias in meta-analysis. In *arXiv*.
<https://doi.org/10.48550/arXiv.2509.07171>

Plummer, M. (2003). JAGS: A program for analysis of Bayesian graphical
models using Gibbs sampling. *Proceedings of the 3rd International
Workshop on Distributed Statistical Computing (DSC 2003)*, 1–10.
