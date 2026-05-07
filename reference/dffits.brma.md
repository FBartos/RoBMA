# DFFITS for brma Objects

Computes DFFITS (Difference in FITS, standardized) for a fitted brma
object. DFFITS measures how much the fitted value for observation \\i\\
changes if observation \\i\\ is removed, standardized by the estimated
standard error of the fit.

## Usage

``` r
# S3 method for class 'brma'
dffits(model, ...)
```

## Arguments

- model:

  a fitted normal-outcome `brma` object without a weightfunction
  component.

- ...:

  additional arguments (currently ignored).

## Value

A named numeric vector of DFFITS values, one for each observation.

## Details

DFFITS values are computed as a PSIS leave-one-out deletion diagnostic.
For each observation \\i\\, the leave-one-out posterior mean fitted
value at that observation is estimated with normalized PSIS weights and
compared to the full-posterior fitted value: \$\$DFFITS_i =
\frac{\hat{\mu}\_i - \hat{\mu}\_{i(-i)}}{SD\_{(-i)}(\mu_i)}\$\$

This targets deletion influence on fitted values directly. It does not
use LOO-PIT residuals, which are predictive outlier diagnostics rather
than fitted-value deletion diagnostics.

Estimate-unit LOO must first be computed with
`model <- add_loo(model, unit = "estimate")`. If the leave-one-out
posterior SD of a fitted value is near zero, the corresponding DFFITS
value is returned as `NA`.

## See also

[`influence.brma`](https://fbartos.github.io/RoBMA/reference/influence.brma.md),
[`cooks.distance.brma`](https://fbartos.github.io/RoBMA/reference/cooks.distance.brma.md),
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  dffits(fit)
}
} # }
```
