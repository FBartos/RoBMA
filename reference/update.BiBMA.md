# Updates a fitted BiBMA object

`update.BiBMA` can be used to

1.  add an additional model to an existing `"BiBMA"` object by
    specifying either a null or alternative prior for each parameter and
    the prior odds of the model (`prior_weights`), see the
    [`vignette("CustomEnsembles")`](articles/CustomEnsembles.md)
    vignette,

2.  change the prior odds of fitted models by specifying a vector
    `prior_weights` of the same length as the fitted models,

3.  refitting models that failed to converge with updated settings of
    control parameters,

4.  or changing the convergence criteria and recalculating the ensemble
    results by specifying new `control` argument and setting
    `refit_failed == FALSE`.

## Usage

``` r
# S3 method for class 'BiBMA'
update(
  object,
  refit_failed = TRUE,
  extend_all = FALSE,
  prior_effect = NULL,
  prior_heterogeneity = NULL,
  prior_baseline = NULL,
  prior_weights = NULL,
  prior_effect_null = NULL,
  prior_heterogeneity_null = NULL,
  prior_baseline_null = NULL,
  study_names = NULL,
  chains = NULL,
  adapt = NULL,
  burnin = NULL,
  sample = NULL,
  thin = NULL,
  autofit = NULL,
  parallel = NULL,
  autofit_control = NULL,
  convergence_checks = NULL,
  save = "all",
  seed = NULL,
  silent = TRUE,
  ...
)
```

## Arguments

- object:

  a fitted BiBMA object

- refit_failed:

  whether failed models should be refitted. Relevant only if new priors
  or `prior_weights` are not supplied. Defaults to `TRUE`.

- extend_all:

  extend sampling in all fitted models based on `"sample_extend"`
  argument in [`set_autofit_control()`](reference/RoBMA_control.md)
  function. Defaults to `FALSE`.

- prior_effect:

  prior distribution for the effect size (`mu`) parameter that will be
  treated as belonging to the alternative hypothesis. Defaults to
  `NULL`.

- prior_heterogeneity:

  prior distribution for the heterogeneity `tau` parameter that will be
  treated as belonging to the alternative hypothesis. Defaults to
  `NULL`.

- prior_baseline:

  prior distribution for the intercepts (`pi`) of each study that will
  be treated as belonging to the alternative hypothesis. Defaults to
  `NULL`.

- prior_weights:

  either a single value specifying prior model weight of a newly
  specified model using priors argument, or a vector of the same length
  as already fitted models to update their prior weights.

- prior_effect_null:

  prior distribution for the effect size (`mu`) parameter that will be
  treated as belonging to the null hypothesis. Defaults to `NULL`.

- prior_heterogeneity_null:

  prior distribution for the heterogeneity `tau` parameter that will be
  treated as belonging to the null hypothesis. Defaults to `NULL`.

- prior_baseline_null:

  prior distribution for the intercepts (`pi`) of each study that will
  be treated as belonging to the null hypothesis. Defaults to `NULL`.

- study_names:

  an optional argument with the names of the studies

- chains:

  a number of chains of the MCMC algorithm.

- adapt:

  a number of adaptation iterations of the MCMC algorithm. Defaults to
  `500`.

- burnin:

  a number of burnin iterations of the MCMC algorithm. Defaults to
  `2000`.

- sample:

  a number of sampling iterations of the MCMC algorithm. Defaults to
  `5000`.

- thin:

  a thinning of the chains of the MCMC algorithm. Defaults to `1`.

- autofit:

  whether the model should be fitted until the convergence criteria
  (specified in `autofit_control`) are satisfied. Defaults to `TRUE`.

- parallel:

  whether the individual models should be fitted in parallel. Defaults
  to `FALSE`. The implementation is not completely stable and might
  cause a connection error.

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

`BiBMA` returns an object of class 'BiBMA'.

## Details

See [`BiBMA()`](reference/BiBMA.md) for more details.

## See also

[`BiBMA()`](reference/BiBMA.md),
[`summary.RoBMA()`](reference/summary.RoBMA.md),
[`prior()`](reference/prior.md),
[`check_setup()`](reference/check_setup.md)
