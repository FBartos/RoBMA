# Internally Standardized Residuals for brma Objects

Computes internally standardized residuals from a fitted brma object
using the hat matrix. Returns a data frame with raw residuals, standard
errors, and standardized residuals (z-values). Available for normal
outcome models only.

## Usage

``` r
# S3 method for class 'brma'
rstandard(model, unit = "estimate", conditioning_depth = "marginal", ...)
```

## Arguments

- model:

  a fitted brma object.

- unit:

  output unit. Only `"estimate"` is implemented currently.

- conditioning_depth:

  conditioning depth. Options are:

  - `"marginal"` (default): Residuals from fixed effects predictions
    (\\observed - X\beta\\).

  - `"cluster"`: Residuals from cluster-level predictions (observed -
    (\\X\beta + gamma\\)). Only available for multilevel (3-level)
    models.

  - `"estimate"`: Residuals from BLUPs, i.e., deviations of the observed
    effect sizes from the best linear unbiased predictions of the
    estimate-specific true effects (observed - theta).

- ...:

  additional arguments (currently ignored)

## Value

A data frame with columns:

- `resid`: Raw residuals

- `se`: Standard errors of the residuals

- `z`: Standardized residuals

## Details

This function returns a data frame with three columns matching the
output of
[`metafor::rstandard`](https://wviechtb.github.io/metafor/reference/residuals.rma.html):

- `resid`: Raw residuals (observed - fitted values)

- `se`: Standard errors of the residuals

- `z`: Standardized residuals (resid / se)

Internally standardized residuals divide the observed residuals by their
corresponding standard errors computed using the hat matrix. For a
correctly specified model, these residuals should approximately follow a
standard normal distribution.

This function is only available for normal outcome models without
selection (weightfunction) bias adjustment. For other model types, use
[`rstudent.brma`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md)
which uses LOO-PIT.

## See also

[`rstudent.brma()`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md),
[`residuals.brma()`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md),
[`blup.brma()`](https://fbartos.github.io/RoBMA/reference/blup.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  # marginal internally standardized residuals (default)
  rstandard(fit)

  # estimate-level (BLUP-based) residuals
  rstandard(fit, conditioning_depth = "estimate")
}
} # }
```
