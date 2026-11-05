# Summarize Zplot Results

Creates summary tables for zplot estimates including EDR, Soric FDR, and
the conventional missing-count diagnostic when applicable.

## Usage

``` r
# S3 method for class 'zplot_brma'
summary(object, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  a zplot_brma object.

- probs:

  quantiles of the posterior distribution to display. Defaults to
  `c(.025, .975)`.

- ...:

  additional arguments (currently unused).

## Value

An object of class `"summary.zplot_brma"` containing the estimates
table.

## Details

The summary includes:

- EDR:

  Expected Discovery Rate - the posterior average probability of
  significance for a new outcome from the bias-adjusted reference

- Soric FDR:

  Expected proportion of false discoveries among significant results,
  computed from EDR following Soric (1989)

- Missing N:

  The retained reference-bin-normalized inverse-weight diagnostic for
  independent univariate and conditional-cluster product models. See
  [`as_zplot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/as_zplot.brma.md)
  for the additional assumptions required by a literal count
  interpretation. Unavailable for multivariate models, retained
  estimate-level or sampling-error sources, applicable integrated
  contexts, or best-rule selection models because relative weights do
  not identify a missing count

The footer reports the Observed Discovery Rate (ODR) with 95% CI for
comparison with the model-estimated EDR.

## See also

[`as_zplot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/as_zplot.brma.md),
[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.zplot_brma.md)
