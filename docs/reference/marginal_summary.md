# Summarize marginal estimates of a fitted RoBMA regression object

`marginal_summary` creates summary tables for marginal estimates of a
RoBMA regression model.

## Usage

``` r
marginal_summary(
  object,
  conditional = FALSE,
  output_scale = NULL,
  probs = c(0.025, 0.975),
  logBF = FALSE,
  BF01 = FALSE
)
```

## Arguments

- object:

  a fitted RoBMA regression object

- conditional:

  show the conditional estimates (assuming that the alternative is
  true).

- output_scale:

  transform the meta-analytic estimates to a different scale. Defaults
  to `NULL` which returns the same scale as the model was estimated on.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- logBF:

  show log of Bayes factors. Defaults to `FALSE`.

- BF01:

  show Bayes factors in support of the null hypotheses. Defaults to
  `FALSE`.

## Value

`marginal_summary` returns a list of tables of class 'BayesTools_table'.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`summary.RoBMA()`](https://fbartos.github.io/RoBMA/reference/summary.RoBMA.md),
[`diagnostics()`](https://fbartos.github.io/RoBMA/reference/diagnostics.md),
[`check_RoBMA()`](https://fbartos.github.io/RoBMA/reference/check_RoBMA.md)
