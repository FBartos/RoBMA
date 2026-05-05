# Cook's Distance for brma Objects

Computes Cook's distance for a fitted brma object. Cook's distance
measures the aggregate influence of each observation on the model
coefficients.

## Usage

``` r
# S3 method for class 'brma'
cooks.distance(model, ...)
```

## Arguments

- model:

  a fitted brma object.

- ...:

  additional arguments (currently ignored).

## Value

A numeric vector of Cook's distance values, one for each observation.

## Details

Cook's distance is computed as a PSIS leave-one-out deletion diagnostic.
For each observation \\i\\, normalized PSIS weights estimate the fitted
values under the leave-one-out posterior. The distance is the posterior
Mahalanobis distance between the full-data and leave-one-out
fitted-value vectors: \$\$D_i = \frac{\Delta_i' V\_\mu^+
\Delta_i}{P}\$\$

where \\\Delta_i = \hat{\mu} - \hat{\mu}\_{(-i)}\\, \\V\_\mu^+\\ is the
generalized inverse of the full-posterior fitted-value covariance, and
\\P\\ is the rank of the fixed-effect model matrix.

## See also

[`influence.brma`](https://fbartos.github.io/RoBMA/reference/influence.brma.md),
[`dffits.brma`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md),
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  cooks.distance(fit)
}
} # }
```
