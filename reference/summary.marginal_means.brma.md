# Summarize Estimated Marginal Means

Summarizes estimated marginal means stored in a `marginal_means.brma`
object using
[`BayesTools::marginal_estimates_table()`](https://fbartos.github.io/BayesTools/reference/BayesTools_ensemble_tables.html).

## Usage

``` r
# S3 method for class 'marginal_means.brma'
summary(
  object,
  type = NULL,
  probs = c(0.025, 0.5, 0.975),
  logBF = FALSE,
  BF01 = FALSE,
  bf = NULL,
  output_measure = NULL,
  transform = NULL,
  ...
)
```

## Arguments

- object:

  a `marginal_means.brma` object.

- type:

  for RoBMA product-space objects, whether to summarize model-averaged
  (`"averaged"`) or conditional (`"conditional"`) marginal means.
  Defaults to `"averaged"` and is available only for RoBMA marginal
  means.

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .50, .975)`.

- logBF:

  whether to show inclusion Bayes factors on the log scale. Defaults to
  `FALSE`.

- BF01:

  whether to show inverse inclusion Bayes factors. Defaults to `FALSE`.

- bf:

  whether to show inclusion Bayes factors. Defaults to the setting
  stored by
  [`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md).

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- ...:

  additional arguments (currently ignored).

## Value

A `BayesTools_table` of class `summary.marginal_means.brma`.
