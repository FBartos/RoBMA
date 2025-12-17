# Estimate a Bayesian Model-Averaged Meta-Regression

`NoBMA.reg` is a wrapper around [`RoBMA.reg()`](reference/RoBMA.reg.md)
that can be used to estimate a publication bias unadjusted Bayesian
model-averaged meta-regression. The interface allows a complete
customization of the ensemble with different prior (or list of prior)
distributions for each component.

## Usage

``` r
NoBMA.reg(
  formula,
  data,
  test_predictors = TRUE,
  study_names = NULL,
  study_ids = NULL,
  transformation = if (any(colnames(data) != "y")) "fishers_z" else "none",
  prior_scale = if (any(colnames(data) != "y")) "cohens_d" else "none",
  standardize_predictors = TRUE,
  priors = NULL,
  model_type = NULL,
  rescale_priors = 1,
  priors_effect = set_default_priors("effect", rescale = rescale_priors),
  priors_heterogeneity = set_default_priors("heterogeneity", rescale = rescale_priors),
  priors_effect_null = set_default_priors("effect", null = TRUE),
  priors_heterogeneity_null = set_default_priors("heterogeneity", null = TRUE),
  priors_hierarchical = set_default_priors("hierarchical"),
  priors_hierarchical_null = set_default_priors("hierarchical", null = TRUE),
  prior_covariates = set_default_priors("covariates", rescale = rescale_priors),
  prior_covariates_null = set_default_priors("covariates", null = TRUE),
  prior_factors = set_default_priors("factors", rescale = rescale_priors),
  prior_factors_null = set_default_priors("factors", null = TRUE),
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
  [`combine_data()`](reference/combine_data.md) for a complete list of
  reserved names and additional information about specifying input data.

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

- study_ids:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

- transformation:

  transformation to be applied to the supplied effect sizes before
  fitting the individual models. Defaults to `"fishers_z"`. We highly
  recommend using `"fishers_z"` transformation since it is the only
  variance stabilizing measure and does not bias PET and PEESE style
  models. The other options are `"cohens_d"`, correlation coefficient
  `"r"` and `"logOR"`. Supplying `"none"` will treat the effect sizes as
  unstandardized and refrain from any transformations.

- prior_scale:

  an effect size scale used to define priors. Defaults to `"cohens_d"`.
  Other options are `"fishers_z"`, correlation coefficient `"r"`, and
  `"logOR"`. The prior scale does not need to match the effect sizes
  measure - the samples from prior distributions are internally
  transformed to match the `transformation` of the data. The
  `prior_scale` corresponds to the effect size scale of default output,
  but can be changed within the summary function.

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

- model_type:

  string specifying the RoBMA ensemble. Defaults to `NULL`. The other
  options are `"PSMA"`, `"PP"`, and `"2w"` which override settings
  passed to the `priors_effect`, `priors_heterogeneity`,
  `priors_effect`, `priors_effect_null`, `priors_heterogeneity_null`,
  `priors_bias_null`, and `priors_effect`. See details for more
  information about the different model types.

- rescale_priors:

  a re-scaling factor for the prior distributions. The re-scaling factor
  allows to adjust the width of all default priors simultaneously.
  Defaults to `1`.

- priors_effect:

  list of prior distributions for the effect size (`mu`) parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to a standard normal distribution
  `prior(distribution = "normal", parameters = list(mean = 0, sd = 1))`.

- priors_heterogeneity:

  list of prior distributions for the heterogeneity `tau` parameter that
  will be treated as belonging to the alternative hypothesis. Defaults
  to
  `prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))`
  that is based on heterogeneities estimates from psychology (van Erp et
  al. 2017) .

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

- priors_hierarchical:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the alternative
  hypothesis. This setting allows users to fit a hierarchical
  (three-level) meta-analysis when `study_ids` are supplied. Note that
  this is an experimental feature and see News for more details.
  Defaults to a beta distribution
  `prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))`.

- priors_hierarchical_null:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the null
  hypothesis. Defaults to `NULL`.

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
  [`set_autofit_control()`](reference/RoBMA_control.md) function. See
  [`?set_autofit_control`](reference/RoBMA_control.md) for options and
  default settings.

- convergence_checks:

  automatic convergence checks to assess the fitted models, passed with
  [`set_convergence_checks()`](reference/RoBMA_control.md) function. See
  [`?set_convergence_checks`](reference/RoBMA_control.md) for options
  and default settings.

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

`NoBMA.reg` returns an object of class 'RoBMA'.

## Details

See [`RoBMA.reg()`](reference/RoBMA.reg.md) for more details.

Note that these default prior distributions are relatively wide and more
informed prior distributions for testing for the presence of moderation
should be considered.

## See also

[`RoBMA()`](reference/RoBMA.md),
[`RoBMA.reg()`](reference/RoBMA.reg.md),
[`summary.RoBMA()`](reference/summary.RoBMA.md),
[`update.RoBMA()`](reference/update.RoBMA.md),
[`check_setup()`](reference/check_setup.md)
