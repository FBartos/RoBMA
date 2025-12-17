# Summarize fitted zcurve_RoBMA object

`summary.zcurve_RoBMA` creates summary tables for a zcurve_RoBMA object.

## Usage

``` r
# S3 method for class 'zcurve_RoBMA'
summary(object, conditional = FALSE, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  a fitted RoBMA object

- conditional:

  show the conditional estimates (assuming that the alternative is
  true). Defaults to `FALSE`. Only available for `type == "ensemble"`.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- ...:

  additional arguments

## Value

`summary.RoBMA` returns a list of tables of class 'BayesTools_table'.

## See also

[`as_zcurve()`](reference/as_zcurve.md)
