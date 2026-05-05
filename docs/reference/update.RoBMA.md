# Updates a fitted RoBMA object

`update.RoBMA` can be used to

1.  add an additional model to an existing `"RoBMA"` object by
    specifying either a null or alternative prior for each parameter and
    the prior odds of the model (`prior_weights`), see the
    [`vignette("CustomEnsembles")`](https://fbartos.github.io/RoBMA/articles/CustomEnsembles.md)
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
# S3 method for class 'RoBMA'
update(
  object,
  refit_failed = TRUE,
  extend_all = FALSE,
  prior_effect = NULL,
  prior_heterogeneity = NULL,
  prior_bias = NULL,
  prior_hierarchical = NULL,
  prior_weights = NULL,
  prior_effect_null = NULL,
  prior_heterogeneity_null = NULL,
  prior_bias_null = NULL,
  prior_hierarchical_null = NULL,
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

  a fitted RoBMA object

- refit_failed:

  whether failed models should be refitted. Relevant only if new priors
  or `prior_weights` are not supplied. Defaults to `TRUE`.

- extend_all:

  extend sampling in all fitted models based on `"sample_extend"`
  argument in
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  function. Defaults to `FALSE`.

- prior_effect:

  prior distribution for the effect size (`mu`) parameter that will be
  treated as belonging to the alternative hypothesis. Defaults to
  `NULL`.

- prior_heterogeneity:

  prior distribution for the heterogeneity `tau` parameter that will be
  treated as belonging to the alternative hypothesis. Defaults to
  `NULL`.

- prior_bias:

  prior distribution for the publication bias adjustment component that
  will be treated as belonging to the alternative hypothesis. Defaults
  to `NULL`.

- prior_hierarchical:

  prior distribution for the correlation of random effects (`rho`)
  parameter that will be treated as belonging to the alternative
  hypothesis. This setting allows users to fit a hierarchical
  (three-level) meta-analysis when `cluster` are supplied. Note that
  this is an experimental feature and see News for more details.
  Defaults to a beta distribution
  `prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))`.

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

- prior_bias_null:

  prior distribution for the publication bias adjustment component that
  will be treated as belonging to the null hypothesis. Defaults to
  `NULL`.

- prior_hierarchical_null:

  prior distribution for the correlation of random effects (`rho`)
  parameter that will be treated as belonging to the null hypothesis.
  Defaults to `NULL`.

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

`RoBMA` returns an object of class 'RoBMA'.

## Details

See [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) for
more details.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`summary.RoBMA()`](https://fbartos.github.io/RoBMA/reference/summary.RoBMA.md),
[`prior()`](https://fbartos.github.io/RoBMA/reference/prior.md),
[`check_setup()`](https://fbartos.github.io/RoBMA/reference/check_setup.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Bem 2011 and fitting the default (RoBMA-PSMA) model
fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study)

# the update function allows us to change the prior model weights of each model
fit1 <- update(fit, prior_weights = c(0, rep(1, 35)))

# add an additional model with different priors specification
# (see '?prior' for more information)
fit2 <- update(fit,
               priors_effect_null = prior("point", parameters = list(location = 0)),
               priors_heterogeneity = prior("normal",
                                  parameters = list(mean = 0, sd = 1),
                                  truncation = list(lower = 0, upper = Inf)),
               priors_bias = prior_weightfunction("one-sided",
                                    parameters = list(cuts = c(.05, .10, .20),
                                                      alpha = c(1, 1, 1, 1))))

# update the models with an increased number of sample iterations
fit3 <- update(fit, autofit_control = set_autofit_control(sample_extend = 1000), extend_all = TRUE)
} # }

```
