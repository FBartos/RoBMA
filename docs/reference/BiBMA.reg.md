# Estimate a Robust Bayesian Meta-Analysis Meta-Regression

`RoBMA` is used to estimate a robust Bayesian meta-regression. The
interface allows a complete customization of the ensemble with different
prior (or list of prior) distributions for each component.

## Usage

``` r
BiBMA.reg(
  formula,
  data,
  test_predictors = TRUE,
  study_names = NULL,
  cluster = NULL,
  standardize_predictors = TRUE,
  priors = NULL,
  rescale_priors = 1,
  priors_effect = set_default_binomial_priors("effect", rescale = rescale_priors),
  priors_heterogeneity = set_default_binomial_priors("heterogeneity", rescale =
    rescale_priors),
  priors_effect_null = set_default_binomial_priors("effect", null = TRUE),
  priors_heterogeneity_null = set_default_binomial_priors("heterogeneity", null = TRUE),
  prior_covariates = set_default_binomial_priors("covariates", rescale = rescale_priors),
  prior_covariates_null = set_default_binomial_priors("covariates", null = TRUE),
  prior_factors = set_default_binomial_priors("factors", rescale = rescale_priors),
  prior_factors_null = set_default_binomial_priors("factors", null = TRUE),
  priors_hierarchical = set_default_binomial_priors("hierarchical"),
  priors_hierarchical_null = set_default_binomial_priors("hierarchical", null = TRUE),
  priors_baseline = set_default_binomial_priors("baseline"),
  priors_baseline_null = set_default_binomial_priors("baseline", null = TRUE),
  algorithm = "bridge",
  chains = 3,
  sample = 5000,
  burnin = 2000,
  adapt = 500,
  thin = 1,
  parallel = FALSE,
  autofit = TRUE,
  autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),
  save = "all",
  seed = NULL,
  silent = TRUE,
  ...
)
```

## Arguments

- formula:

  a formula for the meta-regression model

- data:

  a data.frame containing the data for the meta-regression. Note that
  the column names have to correspond to the effect sizes (`d`, `logOR`,
  `OR`, `r`, `z`), a measure of sampling variability (`se`, `v`, `n`,
  `lCI`, `uCI`, `t`), and the predictors. See
  [`combine_data()`](https://fbartos.github.io/RoBMA/reference/combine_data.md)
  for a complete list of reserved names and additional information about
  specifying input data.

- test_predictors:

  vector of predictor names to test for the presence of moderation
  (i.e., assigned both the null and alternative prior distributions).
  Defaults to `TRUE`, all predictors are tested using the default prior
  distributions (i.e., `prior_covariates`, `prior_covariates_null`,
  `prior_factors`, and `prior_factors_null`). To only estimate and
  adjust for the effect of predictors use `FALSE`. If `priors` is
  specified, any settings in `test_predictors` is overridden.

- study_names:

  an optional argument with the names of the studies

- cluster:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

- standardize_predictors:

  whether continuous predictors should be standardized prior to
  estimating the model. Defaults to `TRUE`. Continuous predictor
  standardization is important for applying the default prior
  distributions for continuous predictors. Note that the resulting
  output corresponds to standardized meta-regression coefficients.

- priors:

  named list of prior distributions for each predictor (with names
  corresponding to the predictors). It allows users to specify both the
  null and alternative hypothesis prior distributions for each predictor
  by assigning the corresponding element of the named list with another
  named list (with `"null"` and `"alt"`). If only one prior is specified
  for a given parameter, it is assumed to correspond to the alternative
  hypotheses and the default null hypothesis is specified (i.e.,
  `prior_covariates_null` or `prior_factors_null`). If a named list with
  only one named prior distribution is provided (either `"null"` or
  `"alt"`), only this prior distribution is used and no default
  distribution is filled in. Parameters without specified prior
  distributions are assumed to be only adjusted for using the default
  alternative hypothesis prior distributions (i.e., `prior_covariates`
  or `prior_factors`). If `priors` is specified, `test_predictors` is
  ignored.

- rescale_priors:

  a re-scaling factor for the prior distributions. The re-scaling factor
  allows to adjust the width of all default priors simultaneously.
  Defaults to `1`.

- priors_effect:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to
  `prior(distribution = "student", parameters = list(location = 0, scale = 0.58, df = 4))`,
  based on logOR meta-analytic estimates from the Cochrane Database of
  Systematic Reviews (Bartoš et al. 2023) .

- priors_heterogeneity:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to
  `prior(distribution = "invgamma", parameters = list(shape = 1.77, scale = 0.55))`
  that is based on heterogeneities of logOR estimates from the Cochrane
  Database of Systematic Reviews (Bartoš et al. 2023) .

- priors_effect_null:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the null hypothesis. Defaults to a
  point null hypotheses at zero,
  `prior(distribution = "point", parameters = list(location = 0))`.

- priors_heterogeneity_null:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the null hypothesis. Defaults to a
  point null hypotheses at zero (a fixed effect meta-analytic models),
  `prior(distribution = "point", parameters = list(location = 0))`.

- prior_covariates:

  a prior distributions for the regression parameter of continuous
  covariates on the effect size under the alternative hypothesis (unless
  set explicitly in `priors`). Defaults to a relatively wide normal
  distribution
  `prior(distribution = "normal", parameters = list(mean = 0, sd = 0.25))`.

- prior_covariates_null:

  a prior distributions for the regression parameter of continuous
  covariates on the effect size under the null hypothesis (unless set
  explicitly in `priors`). Defaults to a no effect
  `prior("spike", parameters = list(location = 0))`.

- prior_factors:

  a prior distributions for the regression parameter of categorical
  covariates on the effect size under the alternative hypothesis (unless
  set explicitly in `priors`). Defaults to a relatively wide
  multivariate normal distribution specifying differences from the mean
  contrasts
  `prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25), contrast = "meandif")`.

- prior_factors_null:

  a prior distributions for the regression parameter of categorical
  covariates on the effect size under the null hypothesis (unless set
  explicitly in `priors`). Defaults to a no effect
  `prior("spike", parameters = list(location = 0))`.

- priors_hierarchical:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the alternative
  hypothesis. This setting allows users to fit a hierarchical
  (three-level) meta-analysis when `cluster` are supplied. Note that
  this is an experimental feature and see News for more details.
  Defaults to a beta distribution
  `prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))`.

- priors_hierarchical_null:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the null
  hypothesis. Defaults to `NULL`.

- priors_baseline:

  prior distributions for the alternative hypothesis about intercepts
  (`pi`) of each study. Defaults to `NULL`.

- priors_baseline_null:

  prior distributions for the null hypothesis about intercepts (`pi`)
  for each study. Defaults to an independent uniform prior distribution
  for each intercept
  `prior("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent")`.

- algorithm:

  a string specifying the algorithm used for the model averaging.
  Defaults to `"bridge"` which results in estimating individual models
  using JAGS and computing the marginal likelihood using bridge
  sampling. An alternative is `"ss"` which uses spike and slab like
  parameterization to approximate the Bayesian model averaging with a
  single model. Note that significantly more `sample`, `burnin`, and
  `adapt` iterations are needed for the `"ss"` algorithm.

- chains:

  a number of chains of the MCMC algorithm.

- sample:

  a number of sampling iterations of the MCMC algorithm. Defaults to
  `5000`.

- burnin:

  a number of burnin iterations of the MCMC algorithm. Defaults to
  `2000`.

- adapt:

  a number of adaptation iterations of the MCMC algorithm. Defaults to
  `500`.

- thin:

  a thinning of the chains of the MCMC algorithm. Defaults to `1`.

- parallel:

  whether the individual models should be fitted in parallel. Defaults
  to `FALSE`. The implementation is not completely stable and might
  cause a connection error.

- autofit:

  whether the model should be fitted until the convergence criteria
  (specified in `autofit_control`) are satisfied. Defaults to `TRUE`.

- autofit_control:

  allows to pass autofit control settings with the
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  function. See
  [`?set_autofit_control`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  for options and default settings.

- convergence_checks:

  automatic convergence checks to assess the fitted models, passed with
  [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  function. See
  [`?set_convergence_checks`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  for options and default settings.

- save:

  whether all models posterior distributions should be kept after
  obtaining a model-averaged result. Defaults to `"all"` which does not
  remove anything. Set to `"min"` to significantly reduce the size of
  final object, however, some model diagnostics and further manipulation
  with the object will not be possible.

- seed:

  a seed to be set before model fitting, marginal likelihood
  computation, and posterior mixing for reproducibility of results.
  Defaults to `NULL` - no seed is set.

- silent:

  whether all print messages regarding the fitting process should be
  suppressed. Defaults to `TRUE`. Note that `parallel = TRUE` also
  suppresses all messages.

- ...:

  additional arguments.

## Value

`RoBMA.reg` returns an object of class 'RoBMA.reg'.

## Details

`BiBMA.reg()` function estimates the Bayesian model-averaged binomial
meta-regression. See
[`vignette("/MetaRegression", package = "RoBMA")`](https://fbartos.github.io/RoBMA/doc/MetaRegression.md)
vignette describes how to use the similar
[`RoBMA.reg()`](https://fbartos.github.io/RoBMA/reference/RoBMA.reg.md)
function to fit Bayesian meta-regression ensembles. See Bartoš et al.
(2025) for more details about the methodology and
[`BiBMA()`](https://fbartos.github.io/RoBMA/reference/BiBMA.md) for more
details about the function options. By default, the function
standardizes continuous predictors. As such, the output should be
interpreted as standardized meta-regression coefficients.

Generic
[`summary.RoBMA()`](https://fbartos.github.io/RoBMA/reference/summary.RoBMA.md),
[`print.RoBMA()`](https://fbartos.github.io/RoBMA/reference/print.RoBMA.md),
and
[`plot.RoBMA()`](https://fbartos.github.io/RoBMA/reference/plot.RoBMA.md)
functions are provided to facilitate manipulation with the ensemble. A
visual check of the individual model diagnostics can be obtained using
the
[`diagnostics()`](https://fbartos.github.io/RoBMA/reference/diagnostics.md)
function. The fitted model can be further updated or modified by
[`update.RoBMA()`](https://fbartos.github.io/RoBMA/reference/update.RoBMA.md)
function. Estimated marginal means can be computed by
[`marginal_summary()`](https://fbartos.github.io/RoBMA/reference/marginal_summary.md)
function and visualized by the
[`marginal_plot()`](https://fbartos.github.io/RoBMA/reference/marginal_plot.md)
function.

## References

Bartoš F, Maier M, Stanley TD, Wagenmakers E (2025). “Robust Bayesian
meta-regression: Model-averaged moderation analysis in the presence of
publication bias.” *Psychological Methods*.
[doi:10.1037/met0000737](https://doi.org/10.1037/met0000737) .\
\
Bartoš F, Otte WM, Gronau QF, Timmers B, Ly A, Wagenmakers E (2023).
“Empirical prior distributions for Bayesian meta-analyses of binary and
time-to-event outcomes.”
[doi:10.48550/arXiv.2306.11468](https://doi.org/10.48550/arXiv.2306.11468)
, Preprint available at https://doi.org/10.48550/arXiv.2306.11468.

## See also

[`BiBMA()`](https://fbartos.github.io/RoBMA/reference/BiBMA.md)
[`summary.RoBMA()`](https://fbartos.github.io/RoBMA/reference/summary.RoBMA.md),
[`update.BiBMA()`](https://fbartos.github.io/RoBMA/reference/update.BiBMA.md),
[`check_setup.reg()`](https://fbartos.github.io/RoBMA/reference/check_setup.reg.md)
