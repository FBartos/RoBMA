# Changelog

## version 4.0.0

### Breaking changes

- rewrites the package around the unified `brma` class hierarchy.
  Single-model fits now use
  [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md),
  [`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md),
  [`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md),
  [`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md), and
  [`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md);
  model-averaged fits use
  [`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md),
  [`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md),
  and [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md).
- removes the legacy `RoBMA.reg()`, `NoBMA()`, `NoBMA.reg()`, `BiBMA()`,
  and `BiBMA.reg()` constructors. Use `mods`, `scale`, and `cluster` in
  the new constructors,
  [`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md) for
  no-bias normal-likelihood model averaging, and
  [`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md)
  for GLMM model averaging.
- replaces old input aliases such as `d`, `r`, `logOR`, `OR`, `z`, `y`,
  `se`, `v`, `n`, `study_names`, `study_ids`, `weight`, and
  `transformation` with `yi`, `vi`/`sei`, `ni`, `slab`, `cluster`,
  `weights`, `measure`, `output_measure`, and `transform`.
- removes legacy helper APIs including `combine_data()`,
  `check_setup()`, `extract_posterior()`, `marginal_summary()`,
  `marginal_plot()`, `plot_models()`, `adjusted_effect()`,
  `as_zcurve()`, and the old z-curve plotting methods.
- normal-likelihood fitting functions now require an explicit `measure`
  for fitted models. Use `measure = "GEN"` for generic effect sizes
  without a known unit-information scale.
- [`update()`](https://rdrr.io/r/stats/update.html) for `brma` objects
  now focuses on extending MCMC samples, updating labels, and refreshing
  cached quantities, not changing model structure.
- [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  no longer accepts the old `remove_failed` and `balance_probability`
  arguments.

### Features

- adds [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) /
  [`brma.norm()`](https://fbartos.github.io/RoBMA/reference/brma.md) for
  single normal-likelihood Bayesian meta-analysis, including
  random-effects, meta-regression, multilevel, and location-scale
  models.
- adds
  [`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
  for binomial-normal and Poisson-normal GLMM meta-analysis from raw
  two-arm counts (`measure = "OR"` and `"IRR"`).
- adds single-model publication-bias constructors
  [`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md),
  [`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md), and
  [`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md).
- adds [`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md) /
  [`BMA.norm()`](https://fbartos.github.io/RoBMA/reference/BMA.md) for
  Bayesian model averaging without publication-bias adjustment.
- adds
  [`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md)
  for Bayesian model averaging of GLMM meta-analyses without
  publication-bias adjustment.
- rewrites
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) as a
  product-space model-averaged ensemble over effect, heterogeneity,
  moderator, scale, and publication-bias components.
- adds formula/data-frame input handling for effect sizes, moderators,
  scale predictors, clusters, labels, subsets, likelihood weights, and
  raw GLMM counts.
- adds default prior construction from standardized effect-size
  measures, estimated or manually supplied unit-information standard
  deviations, and informed empirical priors.
- adds
  [`prior_weightfunction()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md),
  [`wf_cumulative()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md),
  [`wf_fixed()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md),
  and
  [`wf_independent()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md)
  for BayesTools-backed selection-weightfunction priors.
- adds
  [`prior_PET()`](https://fbartos.github.io/RoBMA/reference/prior_PET.md),
  [`prior_PEESE()`](https://fbartos.github.io/RoBMA/reference/prior_PEESE.md),
  [`prior_none()`](https://fbartos.github.io/RoBMA/reference/prior_none.md),
  [`prior_factor()`](https://fbartos.github.io/RoBMA/reference/prior_factor.md),
  [`prior_informed()`](https://fbartos.github.io/RoBMA/reference/prior_informed.md),
  and BayesTools contrast helpers as package-level prior utilities.
- adds `posterior` package interfaces via
  [`as_draws()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md),
  [`as_draws_array()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md),
  [`as_draws_df()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md),
  [`as_draws_list()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md),
  [`as_draws_matrix()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md),
  and
  [`as_draws_rvars()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  for fitted models and `brma_samples`.
- adds the `brma_samples` posterior-sample class with print, summary,
  matrix, and `posterior` conversion methods.
- adds
  [`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
  for posterior predictions of fixed terms, cluster effects, latent true
  effects, observed responses, and scale terms, with `newdata`,
  `conditional`, `bias_adjusted`, `output_measure`, and `transform`
  support.
- adds convenience wrappers
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md),
  [`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md),
  [`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md),
  [`true_effects()`](https://fbartos.github.io/RoBMA/reference/true_effects.md),
  and [`ranef()`](https://fbartos.github.io/RoBMA/reference/ranef.md)
  for `brma` objects.
- adds model-comparison helpers
  [`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
  [`loo()`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
  [`loo_compare()`](https://fbartos.github.io/RoBMA/reference/loo_compare.brma.md),
  [`loo_weights()`](https://fbartos.github.io/RoBMA/reference/loo_weights.brma.md),
  [`check_loo()`](https://fbartos.github.io/RoBMA/reference/check_loo.brma.md),
  [`add_waic()`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md),
  [`waic()`](https://fbartos.github.io/RoBMA/reference/waic.brma.md),
  and [`logLik()`](https://rdrr.io/r/stats/logLik.html) using the `loo`
  package.
- adds bridge-sampling marginal likelihood support for single-model
  `brma` fits via
  [`add_marglik()`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md),
  [`bridge_sampler()`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html),
  [`logml()`](https://rdrr.io/pkg/bridgesampling/man/logml.html),
  [`bf()`](https://rdrr.io/pkg/bridgesampling/man/bf.html),
  [`bayes_factor()`](https://rdrr.io/pkg/bridgesampling/man/bf.html),
  and
  [`post_prob()`](https://rdrr.io/pkg/bridgesampling/man/post_prob.html).
- adds residual and influence diagnostics:
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`rstandard()`](https://rdrr.io/r/stats/influence.measures.html),
  [`rstudent()`](https://rdrr.io/r/stats/influence.measures.html) /
  `LOO-PIT`,
  [`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html),
  [`influence()`](https://rdrr.io/r/stats/lm.influence.html),
  [`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html),
  [`dffits()`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md),
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html),
  [`covratio()`](https://fbartos.github.io/RoBMA/reference/covratio.brma.md),
  and [`vif()`](https://fbartos.github.io/RoBMA/reference/vif.md).
- adds plotting methods for `brma` objects: posterior/prior plots,
  [`funnel()`](https://fbartos.github.io/RoBMA/reference/funnel.md),
  [`regplot()`](https://fbartos.github.io/RoBMA/reference/regplot.md),
  [`qqnorm()`](https://rdrr.io/r/stats/qqnorm.html),
  [`radial()`](https://fbartos.github.io/RoBMA/reference/radial.md) /
  [`galbraith()`](https://fbartos.github.io/RoBMA/reference/radial.md),
  MCMC diagnostic plots, weightfunction plots, and PET-PEESE plots.
- adds
  [`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md)
  with summary and plotting methods for moderator models.
- adds
  [`summary_models()`](https://fbartos.github.io/RoBMA/reference/summary_models.md)
  for marginal and individual model-weight summaries of product-space
  `RoBMA`, `BMA`, and `BMA.glmm` objects.
- adds
  [`interpret()`](https://fbartos.github.io/RoBMA/reference/interpret.md)
  for concise textual interpretation of fitted `brma` and model-averaged
  objects.
- renames the zplot diagnostic API to
  [`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
  and adds the direct plotting wrapper
  [`zplot()`](https://fbartos.github.io/RoBMA/reference/zplot.brma.md),
  with [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`hist()`](https://rdrr.io/r/graphics/hist.html),
  [`lines()`](https://rdrr.io/r/graphics/lines.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and print methods
  for zplot objects.
- adds
  [`RoBMA.options()`](https://fbartos.github.io/RoBMA/reference/RoBMA_options.md)
  and
  [`RoBMA.get_option()`](https://fbartos.github.io/RoBMA/reference/RoBMA_options.md)
  package options for defaults such as core count, automatic
  LOO/WAIC/marginal-likelihood computation, prior scaling defaults, and
  selection-bias defaults.

### Changes

- renames the multilevel clustering argument to `cluster`.
- renames study labels to `slab`, matching `metafor` naming.
- renames likelihood weights to `weights` and applies them consistently
  to posterior fitting, log-likelihoods, LOO, WAIC, and diagnostics.
- uses `measure`, `output_measure`, and `transform` for effect-size
  scale handling. Supported conversions include `SMD`, `COR`, `ZCOR`,
  and `OR`; `transform = "EXP"` exponentiates log ratio measures for
  display.
- standardizes continuous predictors by default and transforms reported
  coefficients back to the original scale unless standardized
  coefficients are requested.
- uses treatment contrasts by default for single-model constructors and
  mean-difference contrasts by default for model-averaged constructors.
- changes
  [`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
  default to `type = "terms"`. GLMM `type = "response"` predictions
  return continuity-corrected effect-size estimators by default via
  `as_measure = TRUE`.
- separates output `unit` from `conditioning_depth` for residuals,
  fitted values, LOO, WAIC, and related diagnostics.
- supports estimate-level and, for multilevel models, cluster-level
  LOO/WAIC targets with target metadata to prevent invalid comparisons.
- keeps bridge-sampling marginal likelihoods for single-model `brma`
  objects; product-space `RoBMA`, `BMA`, and `BMA.glmm` objects relly on
  product-space only.
- routes selection-weightfunction priors through the BayesTools
  selection backend and selected-normal kernel, removing legacy
  weighted-normal mapping paths.
- uses `bias_indicator` and branch-aware selected-normal contexts for
  RoBMA publication-bias mixtures instead of inferring selection
  branches from `omega`.
- increases zplot default posterior thinning controls to `10000` samples
  and accepts `Inf` where full posterior evaluation is requested.
- adds `max_samples` controls to expensive funnel, regplot, and zplot
  summaries.
- updates the package startup message to point users to
  [`vignette("v00-introduction", package = "RoBMA")`](https://fbartos.github.io/RoBMA/articles/v00-introduction.md).
- requires BayesTools 0.3.0 for forward API and selection-backend
  support.
- adds `bridgesampling`, `loo`, `MASS`, and `parallel` as imports and
  `posterior` as a suggested package.

### Fixes

- fixes loading and runtime checks for the RoBMA JAGS module and native
  R routines.

### Performance and internals

- moves fitting to JAGS product-space models with mixture-prior
  indicators for model averaging.
- replaces legacy weighted-normal and multivariate-normal native code
  with selected-normal kernels shared by JAGS and R-native calls.
- adds native selected-normal routines for log likelihoods, normalizers,
  CDFs, moments, RNG, weighted summaries, funnel contours, regplot
  intervals, and zplot densities/threshold summaries.
- adds native GLMM marginal and cluster log-likelihood helpers for
  binomial and Poisson models.
- caches selected-normal normalizers and uses telescoping selection
  probabilities with log-space fallbacks for better numerical stability.
- relocates selected-normal C++ code to `src/selnorm/` and updates
  `Makevars*`, native registration, cleanup rules, and JAGS distribution
  registration.
- removes unused native matrix/LAPACK helper sources and older
  source-level transformation helpers.

### Documentation and tests

- reorganizes vignettes into numbered workflows covering introduction,
  prior distributions, baseline Bayesian meta-analysis, feature
  coverage, metafor parity, model averaging, RoBMA, multilevel models,
  medicine examples, and zplot diagnostics.
- regenerates roxygen documentation for the new constructors, priors,
  predictions, summaries, diagnostics, plots, model-comparison methods,
  and datasets.
- refreshes the README and pkgdown site for the 4.0.0 API.
- adds cached model fits under numbered vignette/model directories.
- refactors tests into ordered input, fitting, prediction, plotting,
  diagnostics, model-comparison, selected-normal kernel, and
  vignette-cache coverage.
- adds regression tests for selected-normal telescope probabilities,
  native/R fallback parity, posterior-row alignment, GLMM response
  conversion, LOO/WAIC targets, bridge sampling, and visual outputs.

## version 3.6.1

CRAN release: 2025-12-17

### Features

- `Explanation` vignette that helps navigate users through the vignettes
- two vignettes demonstrating robust Bayesian meta-analysis and
  meta-regressions
- [`summary()`](https://rdrr.io/r/base/summary.html) function now
  provides publication bias model type summary (`type = "models"`) for
  models fitted using `algorithm = "ss"`
- improves control over zplot diagnostics (i.e., specifying col, border,
  etc for the individual elements)

## version 3.6

### Features

- [`funnel()`](https://fbartos.github.io/RoBMA/reference/funnel.md) plot
  to visualize residuals vs the expected sampling distribution for
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  `RoBMA.reg()` models when using the `algorithm = "ss"`
- [`residuals()`](https://rdrr.io/r/stats/residuals.html) method for
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  `RoBMA.reg()` models when using the `algorithm = "ss"`
- [`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
  function to transform meta-analytic models into a zplot object, only
  available for
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  `RoBMA.reg()` fitted using the `algorithm = "ss"`
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`print()`](https://rdrr.io/r/base/print.html) functions for the
  `as_zplot` objects

## version 3.5.1

CRAN release: 2025-07-28

### Features

- [`summary()`](https://rdrr.io/r/base/summary.html) function now
  supports a `standardized_coefficients` argument to report either
  standardized (default) or raw meta-regression coefficients
- `extract()` function to extract the posterior samples of the model
  parameters
- [`true_effects()`](https://fbartos.github.io/RoBMA/reference/true_effects.md)
  function to summarize the true effect size estimates of
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  `RoBMA.reg()` models when using the `algorithm = "ss"`
- [`predict()`](https://rdrr.io/r/stats/predict.html) method for
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  `RoBMA.reg()` models when using the `algorithm = "ss"`

### Fixes

- fitting a meta-regression using predictors with missing values result
  in a clear error message

### Changes

- improving the speed of unit tests

## version 3.5

### Features

- approximate and computationally feasibly 3lvl selection models via the
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  `RoBMA.reg()` functions with the `cluster` argument when using
  `algorithm = "ss"`
- 3lvl binomial-normal models for binary data via the `BiBMA` and
  `BiBMA.reg` functions with the `cluster` argument when using
  `algorithm = "ss"`
- [`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)
  function to compute the pooled effect size from the `RoBMA.reg`,
  `NoBMA.reg`, and `BiBMA.reg` models
- `adjusted_effect()` function to compute the adjusted effect size from
  the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- enables
  [`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.md)
  for BiBMA models

### Fixes

- passing and checks of the `cluster` and `study_labels` arguments
- PEESE prior distribution now scale as 1/scale instead of 1/scale^2
  with the `rescale_priors` argument\
- the conditional prediction interval based on
  [`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.md)
  is now conditional on the presence of the effect
- additional minor prior handling fixes (i.e., missing marginal
  estimates when only alternative prior distributions were specified
  etc)
- diagnostics with mixture baseline priors when using `algorithm = "ss"`
- [`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.md)
  with only a single study does not produce relative heterogeneity
  instead of crashing

## version 3.4

### Features

- adding binomial-normal meta-regression models for binary data via the
  `BiBMA.reg` function
- the spike and slab algorithm for faster model estimation via the
  `algorithm = "ss"` argument for BiBMA models
- default prior distributions for all parameters of BiBMA models are now
  set via the `set_default_binomial_priors()` function

## version 3.3

### Features

- the spike and slab algorithm for faster model estimation via the
  `algorithm = "ss"` argument (see a new vignette for more details)
- refactoring of the JAGS C++ code of weighted distributions and
  exporting of the lpdfs into JAGS (maintenance)
- weights_mix JAGS prior distribution to sample a mixture of weight
  functions directly

### Fixes

- incorrectly omitting models with more than one predictor when
  computing conditional marginal summary

## version 3.2.1

### Features

- default prior distributions for all parameters are now set via the
  `set_default_priors()` function
- `rescale_priors` argument allows to conveniently re-scale the prior
  distributions for the effect, heterogeneity, and bias simultaneously

## version 3.2

### Features

- [`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.md)
  function to summarize the heterogeneity of the RoBMA models
  (prediction interval, tau, tau^2, I^2, and H^2)
- `check_RoBMA_convergence()` function to check the convergence of the
  RoBMA models
- adds informed prior distributions for binary and time-to-event
  outcomes via BayesTools 0.2.17

### Fixes

- checking and fixing the number of available cores upon loading the
  package (hopefully fixes some parallelization issues)
- [`update()`](https://rdrr.io/r/stats/update.html) function
  re-evaluates convergence checks of individual models
  (<https://github.com/FBartos/RoBMA/issues/34>)
- typos and minor issues in the vignettes

## version 3.1

### Features

- binomial-normal models for binary data via the `BiBMA` function
- `NoBMA` and `NoBMA.reg()` functions as wrappers around `RoBMA`
  `RoBMA.reg()` functions for simpler specification of publication bias
  unadjusted Bayesian model-averaged meta-analysis
- adding odds ratios output transformation\`
- extending (instead of a complete refitting) of models via the
  `update.RoBMA()` function (only non-converged models by default or all
  by setting `extend_all = TRUE`)

### Fixes

- handling of non-converged models

## version 3.0.1

CRAN release: 2023-06-02

### Fixes (thanks to Don & Rens)

- compilation issues with Clang
  (<https://github.com/FBartos/RoBMA/issues/28>)
- lapack path specifications
  (<https://github.com/FBartos/RoBMA/issues/24>)

## version 3.0

### Features

- meta-regression with `RoBMA.reg()` function
- posterior marginal summary and plots for the `RoBMA.reg` models with
  `summary_marginal()` and `plot_marginal()` functions
- new vignette on hierarchical Bayesian model-averaged meta-analysis
- new vignette on robust Bayesian model-averaged meta-regression
- adding vignette from AMPPS tutorial
- faster implementation of JAGS multivariate normal distribution (based
  on the BUGS JAGS module)
- incorporating `weight` argument in the `RoBMA` and `combine_data`
  functions in order to pass `custom` likelihood weights
- ability to use inverse square weights in the weighted meta-analysis by
  setting a `weighted_type = "inverse_sqrt"` argument

### Changes

- reworked interface for the hierarchical models. Prior distributions
  are now specified via the `priors_hierarchical` and
  `priors_hierarchical_null` arguments instead of `priors_rho` and
  `priors_rho_null`. The model summary now shows `Hierarchical`
  component summary.

## version 2.3.2

CRAN release: 2023-03-13

### Fixes

- suppressing start-up message
- cleaning up imports

## version 2.3.1

CRAN release: 2022-07-16

### Fixes

- fixing weighted meta-analysis parameterization

## version 2.3

### Features

- weighted meta-analysis by specifying `cluster` argument in
  [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) and
  setting `weighted = TRUE`. The likelihood contribution of estimates
  from each study is down-weighted proportionally to the number of
  estimates in that study. Note that this experimental feature is
  supposed to provide a conservative alternative for estimating RoBMA in
  cases with multiple estimates from a study where the multivariate
  option is not computationally feasible.

## version 2.2.3

### Fixes

- updating the Makevars to install with R 4.2 and JAGS 4.3.1

## version 2.2.2

CRAN release: 2022-04-20

### Fixes

- updating the C++ to compile on M1 Mac

## version 2.2.1

CRAN release: 2022-04-06

### Changes

- message about the effect size scale of parameter estimates is always
  shown
- compatibility with BayesTools 0.2.0+

## version 2.2

### Features

- three-level meta-analysis by specifying `cluster` argument in `RoBMA`.
  However, note that this is (1) an experimental feature and (2) the
  computational expense of fitting selection models with clustering is
  extreme. As of now, it is almost impossible to have more than 2-3
  estimates clustered within a single study).

## version 2.1.2

CRAN release: 2022-01-12

### Fixes

- adding Windows ucrt patch (thanks to Tomas Kalibera)
- adding BayesTools version check

## version 2.1.1

CRAN release: 2021-11-03

### Fixes

- incorrectly formatted citations in vignettes and capitalization

### Features

- adding `informed_prior()` function (from the BayesTools package) that
  allows specification of various informed prior distributions from the
  field of medicine and psychology
- adding a vignette reproducing the example of dentine sensitivity with
  the informed Bayesian model-averaged meta-analysis from Bartoš et al.,
  2021
  ([open-access](https://onlinelibrary.wiley.com/doi/10.1002/sim.9170)),
- further reductions of fitted object size when setting `save = "min"`

## version 2.1

### Fixes

- more informative error message when the JAGS module fails to load
- correcting wrong PEESE transformation for the individual models
  summaries (issue [\#12](https://github.com/FBartos/RoBMA/issues/12))
- fixing error message for missing conditional PET-PEESE
- fixing incorrect lower bound check for log(OR)

### Features

- adding
  [`interpret()`](https://fbartos.github.io/RoBMA/reference/interpret.md)
  function (issue [\#11](https://github.com/FBartos/RoBMA/issues/11))
- adding effect size transformation via `output_scale` argument to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  `plot_models()` functions
- better handling of effect size transformations and scaling -
  BayesTools style back-end functions with Jacobian transformations

## version 2.0

Please notice that this is a major release that breaks backwards
compatibility.

### Changes

- naming of the arguments specifying prior distributions for the
  different parameters/components of the models changed (`priors_mu` -\>
  `priors_effect`, `priors_tau` -\> `priors_heterogeneity`, and
  `priors_omega` -\> `priors_bias`),
- prior distributions for specifying weight functions now use a
  dedicated function
  (`prior(distribution = "two.sided", parameters = ...)` -\>
  `prior_weightfunction(distribution = "two.sided", parameters = ...)`),
- new dedicated function for specifying no publication bias adjustment
  component / no heterogeneity component
  ([`prior_none()`](https://fbartos.github.io/RoBMA/reference/prior_none.md)),
- new dedicated functions for specifying models with the PET and PEESE
  publication bias adjustments
  (`prior_PET(distribution = "Cauchy", parameters = ...)` and
  `prior_PEESE(distribution = "Cauchy", parameters = ...)`),
- new default prior distribution specification for the publication bias
  adjustment part of the models (corresponding to the RoBMA-PSMA model
  from Bartoš et al., 2021
  [manuscript](https://doi.org/10.1002/jrsm.1594)),
- new `model_type` argument allowing to specify different “pre-canned”
  models (`"PSMA"` = RoBMA-PSMA, `"PP"` = RoBMA-PP, `"2w"` =
  corresponding to Maier et al., in press ,
  [manuscript](https://doi.org/10.1037/met0000405)),
- `combine_data` function allows combination of different effect sizes /
  variability measures into a common effect size measure (also used from
  within the `RoBMA` function),
- better and improved automatic fitting procedure now enabled by default
  (can be turned of with `autofit = FALSE`)
- prior distributions can be specified on the different scale than the
  supplied effect sizes (the package fits the model on Fisher’s z scale
  and back transforms the results back to the scale that was used for
  prior distributions specification, Cohen’s d by default, but both of
  them can be overwritten with the `prior_scale` and `transformation`
  arguments),
- new prior distributions, e.g., beta or fixed weight functions,
- estimates from individual models are now plotted with the
  `plot_models()` function and the forest plot can be obtained with the
  [`forest()`](https://wviechtb.github.io/metafor/reference/forest.html)
  function,
- the posterior distribution plots for the individual weights are no
  able supported, however, the weightfunction and the PET-PEESE
  publication bias adjustments can be visualized with the `plot.RoBMA()`
  function and `parameter = "weightfunction"` and
  `parameter = "PET-PEESE"`.

## version 1.2.1

CRAN release: 2021-02-16

### Fixes

- check_setup function not working at all

## version 1.2.0

CRAN release: 2021-01-21

### Changes

- the studies’s true effects are now marginalized out of the random
  effects models and are no longer estimated (see Appendix A of our
  [manuscript](https://doi.org/10.1037/met0000405) for more details). As
  a results, arguments referring to the true effects are now disabled.
- all models are now being estimated using the likelihood of effect
  sizes (instead of test-statistics as usually defined). We reproduced
  the simulation study that we used to evaluate the method performance
  and it achieved identical results (up to MCMC error, before
  marginalizing out the true effects). A big advantage of using the
  normal likelihood for effect sizes is a considerable speed up of the
  whole estimation process.
- as a results of these two changes, the results of the models will
  differ to those of pre 1.2.0 version

### Fixes

- autofit being turn on if any control argument was specified

## version 1.1.2

CRAN release: 2020-12-10

### Fixes

- vdiffr not being used conditionally in unit tests

## version 1.1.1

CRAN release: 2020-11-10

### Fixes

- inability to fit a model without specifying a seed
- inability to produce individual model plots due to incompatibility
  with the newer versions of ggplot2

## version 1.1.0

CRAN release: 2020-10-30

### Features

- parallel within and between model fitting using the parallel package
  with ‘parallel = TRUE’ argument

## version 1.0.5

CRAN release: 2020-10-13

### Fixes:

- models being fitted automatically until reaching R-hat lower than 1.05
  without setting max_rhat and autofit control parameters
- bug preventing to draw a bivariate plot of mu and tau
- range for parameter estimates from individual models no containing 0
  (or 1 in case of OR measured effect sizes)
- inability to fit a model with only null mu distributions if
  correlation or OR measured effect sizes were specified
- ordering of the estimated and observed effects when both of them are
  requested simultaneously
- formatting of this file (NEWS.md)

### Improvements:

- priors plot: parameter specification, default plotting range, clearer
  x-axis labels in cases when the parameter is defined on transformed
  scale
- parameters plots: probability scale always ends at the same spot as is
  the last tick on the density scale
- adding warnings if any of the specified models has Rhat higher than
  1.05 or the specified value
- grouping the same warnings messages together

## version 1.0.4

CRAN release: 2020-08-07

### Fixes:

- inability to run models without the silent = TRUE control

## version 1.0.3

CRAN release: 2020-08-06

### Features:

- x-axis rescaling for the weight function plot (by setting ‘rescale_x =
  TRUE’ in the ‘plot.RoBMA’ function)
- setting expected direction of the effect in for RoBMA function

### Fixes:

- marginal likelihood calculation for models with spike prior
  distribution on mean parameter which location was not set to 0
- some additional error messages

### CRAM requested changes:

- changing information messages from ‘cat’ to ‘message’ from plot
  related functions
- saving and returning the ‘par’ settings to the user defined one in the
  base plot functions

## version 1.0.2

### Fixes:

- the summary and plot function now shows quantile based confidence
  intervals for individual models instead of the HPD provided before
  (this affects only ‘summary’/‘plot’ with ‘type = “individual”’, all
  other confidence intervals were quantile based before)

## version 1.0.1

### Fixes:

- summary function returning median instead of mean

## version 1.0.0 (vs the osf version)

### Fixes:

- incorrectly weighted theta estimates
- models with non-zero point prior distribution incorrectly plotted
  using when “models” option in case that the mu parameter was
  transformed

### Additional features:

- analyzing OR
- distributions implemented using boost library (helps with convergence
  issues)
- ability to mute the non-suppressible “precision not achieved” warning
  messages by using “silent” = TRUE inside of the control argument
- vignettes

### Notable changes:

- the way how the seed is set before model fitting (the simulation study
  will not be reproducible with the new version of the package)
