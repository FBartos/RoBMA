# Estimate a Bayesian Model-Averaged Meta-Analysis of Binomial Data

`BiBMA` estimate a binomial-normal Bayesian model-averaged
meta-analysis. The interface allows a complete customization of the
ensemble with different prior (or list of prior) distributions for each
component.

## Usage

``` r
BiBMA(
  x1,
  x2,
  n1,
  n2,
  study_names = NULL,
  cluster = NULL,
  rescale_priors = 1,
  priors_effect = set_default_binomial_priors("effect", rescale = rescale_priors),
  priors_heterogeneity = set_default_binomial_priors("heterogeneity", rescale =
    rescale_priors),
  priors_effect_null = set_default_binomial_priors("effect", null = TRUE),
  priors_heterogeneity_null = set_default_binomial_priors("heterogeneity", null = TRUE),
  priors_hierarchical = set_default_binomial_priors("hierarchical"),
  priors_hierarchical_null = set_default_binomial_priors("hierarchical", null = TRUE),
  priors_baseline = set_default_binomial_priors("baseline"),
  priors_baseline_null = set_default_binomial_priors("baseline", null = TRUE),
  chains = 3,
  sample = 5000,
  burnin = 2000,
  adapt = 500,
  thin = 1,
  parallel = FALSE,
  autofit = TRUE,
  autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),
  algorithm = "bridge",
  save = "all",
  seed = NULL,
  silent = TRUE,
  ...
)
```

## Arguments

- x1:

  a vector with the number of successes in the first group

- x2:

  a vector with the number of successes in the second group

- n1:

  a vector with the number of observations in the first group

- n2:

  a vector with the number of observations in the second group

- study_names:

  an optional argument with the names of the studies

- cluster:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

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

- algorithm:

  a string specifying the algorithm used for the model averaging.
  Defaults to `"bridge"` which results in estimating individual models
  using JAGS and computing the marginal likelihood using bridge
  sampling. An alternative is `"ss"` which uses spike and slab like
  parameterization to approximate the Bayesian model averaging with a
  single model. Note that significantly more `sample`, `burnin`, and
  `adapt` iterations are needed for the `"ss"` algorithm.

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

The `BiBMA()` function estimates the binomial-normal Bayesian
model-averaged meta-analysis described in Bartoš et al. (2023) . See
[[`vignette("MedicineBiBMA", package = "RoBMA")`](https://fbartos.github.io/RoBMA/articles/MedicineBiBMA.md)](https://fbartos.github.io/RoBMA/doc/MedicineBiBMA.md)
vignette for a reproduction of the Oduwole et al. (2018) example. Also
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) for
additional details.

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
function.

## References

Bartoš F, Otte WM, Gronau QF, Timmers B, Ly A, Wagenmakers E (2023).
“Empirical prior distributions for Bayesian meta-analyses of binary and
time-to-event outcomes.”
[doi:10.48550/arXiv.2306.11468](https://doi.org/10.48550/arXiv.2306.11468)
, Preprint available at https://doi.org/10.48550/arXiv.2306.11468.\
\
Oduwole O, Udoh EE, Oyo-Ita A, Meremikwu MM (2018). “Honey for acute
cough in children.” *Cochrane Database of Systematic Reviews*.
[doi:10.1002/14651858.CD007094.pub5](https://doi.org/10.1002/14651858.CD007094.pub5)
.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`summary.RoBMA()`](https://fbartos.github.io/RoBMA/reference/summary.RoBMA.md),
[`update.RoBMA()`](https://fbartos.github.io/RoBMA/reference/update.RoBMA.md),
[`check_setup()`](https://fbartos.github.io/RoBMA/reference/check_setup.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Oduwole (2018) and reproducing the example from
# Bartos et al. (2023) with domain specific informed prior distributions

fit <- BiBMA(
  x1          = c(5, 2),
  x2          = c(0, 0),
  n1          = c(35, 40),
  n2          = c(39, 40),
  priors_effect        = prior_informed(
      "Acute Respiratory Infections",
      type = "logOR", parameter = "effect"),
  priors_heterogeneity = prior_informed(
      "Acute Respiratory Infections",
      type = "logOR", parameter = "heterogeneity")
 )

 summary(fit)

 # produce summary on OR scale
 summary(fit, output_scale = "OR")

} # }
```
