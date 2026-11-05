# COVRATIO for brma Objects

Computes COVRATIO for a fitted brma object. COVRATIO measures the change
in the determinant of the covariance matrix of the estimates when
observation \\i\\ is removed.

## Usage

``` r
# S3 method for class 'brma'
covratio(model, type = "mods", component = NULL, parameter = NULL, ...)
```

## Arguments

- model:

  a fitted brma object.

- type:

  type of parameters to be summarized. Defaults to `"mods"` (for the
  effect size and meta-regression coefficients). Use `"scale"` for
  heterogeneity and scale-regression coefficients.

- component:

  optional parameter namespace. Use `"random"` with one explicitly
  selected semantic random-effect quantity.

- parameter:

  semantic random-effect quantity used when `component = "random"`.

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
precision (increases variance). For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
objects, COVRATIO uses estimate-unit PSIS weights. With correlated
known-`V` data, this is influence under conditional estimate deletion,
not independent-study deletion. A requested random-effect correlation
must be defined in every posterior draw. If either SD is zero in some
draws, the complete influence request is unavailable; explicitly select
an everywhere-defined quantity with `component = "random"` and
`parameter` instead. Parameters with zero posterior variance are
excluded from the covariance determinant and reported in a printed note
when available. Ordered terms retain their fitted totals and omit one
share per metadata-declared Dirichlet simplex. Aliases of a shared
allocation source contribute its independent coordinates only once. If
no estimable parameters remain, or if a full or LOO covariance
determinant is zero or non-finite after exclusion, COVRATIO is undefined
and values are reported as `NaN`.

## See also

[`influence.brma`](https://fbartos.github.io/RoBMA/dev/reference/influence.brma.md),
[`dffits.brma`](https://fbartos.github.io/RoBMA/dev/reference/dffits.brma.md),
[`cooks.distance.brma`](https://fbartos.github.io/RoBMA/dev/reference/cooks.distance.brma.md)

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
