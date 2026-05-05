# Prints summary of `"RoBMA"` ensemble implied by the specified priors

`check_setup` prints summary of `"RoBMA"` ensemble implied by the
specified prior distributions. It is useful for checking the ensemble
configuration prior to fitting all of the models.

## Usage

``` r
check_setup(
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  models = FALSE,
  silent = FALSE
)

check_setup.RoBMA(
  model_type = NULL,
  priors_effect = prior(distribution = "normal", parameters = list(mean = 0, sd = 1)),
  priors_heterogeneity = prior(distribution = "invgamma", parameters = list(shape = 1,
    scale = 0.15)),
  priors_bias = list(prior_weightfunction(distribution = "two.sided", parameters =
    list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.1)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights =
    1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha =
    c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), 
    
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1,
    1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution =
    "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)),
    prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0, 1),
    truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy",
    parameters = list(0, 5), truncation = list(0, Inf), prior_weights = 1/4)),
  priors_effect_null = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null = prior(distribution = "point", parameters = list(location =
    0)),
  priors_bias_null = prior_none(),
  priors_hierarchical = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null = NULL,
  models = FALSE,
  silent = FALSE
)
```

## Arguments

- model_type:

  string specifying the RoBMA ensemble. Defaults to `NULL`. The other
  options are `"PSMA"`, `"PP"`, and `"2w"` which override settings
  passed to the `priors_effect`, `priors_heterogeneity`,
  `priors_effect`, `priors_effect_null`, `priors_heterogeneity_null`,
  `priors_bias_null`, and `priors_effect`. See details for more
  information about the different model types.

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

- priors_bias:

  list of prior distributions for the publication bias adjustment
  component that will be treated as belonging to the alternative
  hypothesis. Defaults to
  `list( prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.10)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.5)), prior_weights = 1/12), prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12), prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf), prior_weights = 1/4), prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf), prior_weights = 1/4) )`,
  corresponding to the RoBMA-PSMA model introduce by Bartoš et
  al. (2023) .

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

- priors_bias_null:

  list of prior weight functions for the `omega` parameter that will be
  treated as belonging to the null hypothesis. Defaults no publication
  bias adjustment,
  [`prior_none()`](https://fbartos.github.io/RoBMA/reference/prior_none.md).

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

- models:

  should the models' details be printed.

- silent:

  do not print the results.

## Value

`check_setup` invisibly returns list of summary tables.

## See also

[`check_setup.reg()`](https://fbartos.github.io/RoBMA/reference/check_setup.reg.md)
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)
