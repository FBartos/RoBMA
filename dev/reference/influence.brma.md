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
  PSIS leave-one-out posterior mean of the root mean square of total
  row-level heterogeneity over the remaining scale-model rows after each
  deletion, the root mean square over the observed design that
  [`summary_heterogeneity`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md)
  reports. It differs from the arithmetic mean of the
  observation-specific tau that 4.0.0 used, which is never larger. It is
  not a pooled log-scale design prediction.

- dfbs:

  A data frame with DFBETAS values for the model coefficients.

Undefined COVRATIO and coefficient DFBETAS values are reported as `NaN`
with an explanatory note. DFFITS returns `NA` when its leave-one-out
posterior fitted-value SD is zero, without a separate note. For binomial
and Poisson GLMMs the `rstudent` component is the descriptive normal
approximation of
[`rstudent.brma`](https://fbartos.github.io/RoBMA/dev/reference/rstudent.brma.md)
(and the DFFITS, Cook's distance, and hat columns are not available).
For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
objects, influence returns only components with implemented
estimate-unit targets. DFFITS and Cook's distance report influence on
fixed-location fitted values \\\mu = X\beta\\; `tau.del` is currently
deferred. With correlated known-`V` data, the returned PSIS influence
components use conditional estimate deletion rather than
independent-study deletion.

## Details

The `hat` column uses the design of each draw's member for ensembles
whose members differ in their PET and PEESE terms, as described in
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/dev/reference/hatvalues.brma.md).

Selection
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
models are currently unavailable through this aggregate diagnostic.
Standalone DFFITS and Cook's distance are available for normal outcomes
without a weightfunction component.

For known-`V` models, the hat values require applying the inverse of the
marginal covariance. DFFITS and Cook's distance use PSIS fixed-location
moments and do not apply this inverse. A covariance that is too
ill-conditioned for that (a reciprocal condition number of about `2e-10`
or less in any posterior draw used) is unavailable: the call stops with
an error of class `RoBMA_known_v_precision_unavailable`, whose parent
class is `RoBMA_numerical_unavailable`, instead of returning a rounded
number. The condition carries the draw positions (`rows`) and their
`reciprocal_condition`.

## See also

[`dffits.brma`](https://fbartos.github.io/RoBMA/dev/reference/dffits.brma.md),
[`cooks.distance.brma`](https://fbartos.github.io/RoBMA/dev/reference/cooks.distance.brma.md),
[`covratio.brma`](https://fbartos.github.io/RoBMA/dev/reference/covratio.brma.md),
[`dfbetas.brma`](https://fbartos.github.io/RoBMA/dev/reference/dfbetas.brma.md),
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/dev/reference/hatvalues.brma.md)

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
