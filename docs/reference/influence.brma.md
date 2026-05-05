# Measure Influence for brma Objects

Computes DFFITS, Cook's distance, COVRATIO, and other influence
diagnostics for a fitted brma object.

## Usage

``` r
# S3 method for class 'brma'
influence(model, ...)
```

## Arguments

- model:

  a fitted brma object.

- ...:

  additional arguments (currently ignored).

## Value

An object of class `"infl.brma"`, which corresponds to the structure of
[`metafor::influence`](https://wviechtb.github.io/metafor/reference/influence.rma.uni.html)
objects. It is a list containing:

- inf:

  A data frame with columns: `rstudent`, `dffits`, `cook.d`, `cov.r`,
  `tau.del`, and `hat`, where available. The `tau.del` column is the
  PSIS leave-one-out posterior mean of the aggregate heterogeneity. For
  scale models, aggregation is over the remaining scale-model rows after
  each deletion.

- dfbs:

  A data frame with DFBETAS values for the model coefficients.

Undefined determinant- or variance-standardized diagnostics are reported
as `NaN` and printed with an explanatory note.

## See also

[`dffits.brma`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md),
[`cooks.distance.brma`](https://fbartos.github.io/RoBMA/reference/cooks.distance.brma.md),
[`covratio.brma`](https://fbartos.github.io/RoBMA/reference/covratio.brma.md),
[`dfbetas.brma`](https://fbartos.github.io/RoBMA/reference/dfbetas.brma.md),
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )
  fit <- add_loo(fit)

  inf <- influence(fit)
  print(inf)
}
} # }
```
