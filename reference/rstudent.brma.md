# Externally Standardized (Studentized) Residuals for brma Objects

Computes externally standardized residuals (also called studentized
residuals or standardized deleted residuals) from a fitted brma object
using LOO-PIT (Leave-One-Out Probability Integral Transform). Returns a
data frame with raw residuals, standard errors, and standardized
residuals (z-values).

## Usage

``` r
# S3 method for class 'brma'
rstudent(model, unit = "estimate", conditioning_depth = "marginal", ...)
```

## Arguments

- model:

  a fitted brma object.

- unit:

  output unit. Only `"estimate"` is available for LOO-PIT residuals.

- conditioning_depth:

  unused for LOO-PIT residuals. LOO-PIT residuals always use the
  estimate-unit LOO target.

- ...:

  additional arguments (currently ignored)

## Value

A data frame with columns:

- `resid`: Raw residuals

- `se`: Standard errors of the residuals

- `z`: Externally standardized residuals (LOO-PIT)

## Details

This function returns a data frame with three columns matching the
output of
[`metafor::rstudent`](https://wviechtb.github.io/metafor/reference/residuals.rma.html):

- `resid`: LOO predictive residuals (observed - fitted values)

- `se`: LOO predictive standard errors when available

- `z`: Externally standardized residuals (LOO-PIT transformed)

LOO-PIT residuals are the Bayesian equivalent of studentized deleted
residuals. They are computed via leave-one-out probability integral
transformation using Pareto smoothed importance sampling. For each
observation, the LOO-weighted CDF is computed and transformed to a
standard normal quantile.

Under a correctly specified model, LOO-PIT residuals should follow a
standard normal distribution. Large absolute values may indicate
outliers or model misspecification.

The `z` column is the primary standardized diagnostic. The `resid` and
`se` columns are raw-scale companions computed from LOO predictive
moments using the normalized PSIS weights. For selection models, these
moments are computed from the fitted selected-normal predictive
distribution. For GLMMs, they are computed on the approximate
effect-size scale used by the LOO-PIT diagnostic; they are not exact PIT
diagnostics for the raw count likelihood.

Unlike
[`rstandard.brma`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md)
(which uses the hat matrix), LOO-PIT residuals properly account for
estimation uncertainty and leverage without requiring explicit hat
matrix computation. This makes `rstudent.brma` suitable for all model
types including selection models and GLMMs.

## See also

[`rstandard.brma()`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md),
[`residuals.brma()`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md),
[`loo.brma()`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`blup.brma()`](https://fbartos.github.io/RoBMA/reference/blup.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  # externally standardized residuals
  rstudent(fit)

  # check Pareto k values
  plot(loo(fit))
}
} # }
```
