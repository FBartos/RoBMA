# Estimate a Bayesian Model-Averaged Meta-Analysis

`NoBMA` is a wrapper around
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) that can
be used to estimate a publication bias unadjusted Bayesian
model-averaged meta-analysis. The interface allows a complete
customization of the ensemble with different prior (or list of prior)
distributions for each component.

## Usage

``` r
NoBMA(
  d = NULL,
  r = NULL,
  logOR = NULL,
  OR = NULL,
  z = NULL,
  y = NULL,
  se = NULL,
  v = NULL,
  n = NULL,
  lCI = NULL,
  uCI = NULL,
  t = NULL,
  study_names = NULL,
  cluster = NULL,
  data = NULL,
  weight = NULL,
  transformation = if (is.null(y)) "fishers_z" else "none",
  prior_scale = if (is.null(y)) "cohens_d" else "none",
  model_type = NULL,
  rescale_priors = 1,
  priors_effect = set_default_priors("effect", rescale = rescale_priors),
  priors_heterogeneity = set_default_priors("heterogeneity", rescale = rescale_priors),
  priors_effect_null = set_default_priors("effect", null = TRUE),
  priors_heterogeneity_null = set_default_priors("heterogeneity", null = TRUE),
  priors_hierarchical = set_default_priors("hierarchical"),
  priors_hierarchical_null = set_default_priors("hierarchical", null = TRUE),
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

- d:

  a vector of effect sizes measured as Cohen's d / Hedges' g
  (standardized mean differences)

- r:

  a vector of effect sizes measured as correlations

- logOR:

  a vector of effect sizes measured as log odds ratios

- OR:

  a vector of effect sizes measured as odds ratios

- z:

  a vector of effect sizes measured as Fisher's z

- y:

  a vector of unspecified effect sizes (note that effect size
  transformations are unavailable with this type of input)

- se:

  a vector of standard errors of the effect sizes

- v:

  a vector of variances of the effect sizes

- n:

  a vector of overall sample sizes

- lCI:

  a vector of lower bounds of confidence intervals

- uCI:

  a vector of upper bounds of confidence intervals

- t:

  a vector of t/z-statistics

- study_names:

  an optional argument with the names of the studies

- cluster:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

- data:

  a data object created by the `combine_data` function. This is an
  alternative input entry to specifying the `d`, `r`, `y`, etc...
  directly. I.e., RoBMA function does not allow passing a data.frame and
  referencing to the columns.

- weight:

  specifies likelihood weights of the individual estimates. Notes that
  this is an untested experimental feature.

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
  (three-level) meta-analysis when `cluster` are supplied. Note that
  this is an experimental feature and see News for more details.
  Defaults to a beta distribution
  `prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))`.

- priors_hierarchical_null:

  list of prior distributions for the correlation of random effects
  (`rho`) parameter that will be treated as belonging to the null
  hypothesis. Defaults to `NULL`.

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

`NoBMA` returns an object of class 'RoBMA'.

## Details

See [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) for
more details.

Note that these default prior distributions are relatively wide and more
informed prior distributions for testing for the presence of moderation
should be considered.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`summary.RoBMA()`](https://fbartos.github.io/RoBMA/reference/summary.RoBMA.md),
[`update.RoBMA()`](https://fbartos.github.io/RoBMA/reference/update.RoBMA.md),
[`check_setup()`](https://fbartos.github.io/RoBMA/reference/check_setup.md)
