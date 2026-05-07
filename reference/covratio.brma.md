# COVRATIO for brma Objects

Computes COVRATIO for a fitted brma object. COVRATIO measures the change
in the determinant of the covariance matrix of the estimates when
observation \\i\\ is removed.

## Usage

``` r
# S3 method for class 'brma'
covratio(model, type = "mods", ...)
```

## Arguments

- model:

  a fitted brma object.

- type:

  type of parameters to be summarized. Defaults to `"mods"` (for the
  effect size and meta-regression coefficients). Use `"scale"` for
  heterogeneity and scale-regression coefficients.

- ...:

  additional arguments. The internal `.weights` argument can supply
  precomputed PSIS weights for callers that already extracted them.

## Value

A named numeric vector of COVRATIO values, one for each observation.

## Details

COVRATIO is computed using importance sampling weights to approximate
the leave-one-out covariance matrices without refitting the model.
Estimate-unit LOO must first be computed with
`model <- add_loo(model, unit = "estimate")`, unless internal weights
are supplied.

\$\$COVRATIO_i = \frac{\det(Cov(\beta)\_{-i})}{\det(Cov(\beta))}\$\$

Values \> 1 indicate that the observation improves precision (decreases
variance), while values \< 1 indicate that the observation decreases
precision (increases variance). If any included parameter has zero
posterior variance, or if a full or LOO covariance determinant is zero
or non-finite, COVRATIO is undefined. In that case, values are reported
as `NaN` with a printed note when available.

## See also

[`influence.brma`](https://fbartos.github.io/RoBMA/reference/influence.brma.md),
[`dffits.brma`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md),
[`cooks.distance.brma`](https://fbartos.github.io/RoBMA/reference/cooks.distance.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
  fit <- add_loo(fit)

  covratio(fit)
}
} # }
```
