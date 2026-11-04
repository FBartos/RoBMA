# Internally Standardized Residuals for brma Objects

Computes internally standardized residuals from a fitted brma object
using the hat matrix or the corresponding known-`V` GLS/BLUP projection.
Returns a data frame with raw residuals, standard errors, and
standardized residuals (z-values). Available for normal outcome models
only.

## Usage

``` r
# S3 method for class 'brma'
rstandard(
  model,
  unit = "estimate",
  conditioning_depth = "marginal",
  max_samples = Inf,
  ...
)
```

## Arguments

- model:

  a fitted brma object.

- unit:

  output unit. Only `"estimate"` is implemented currently.

- conditioning_depth:

  conditioning depth. Options are:

  - `"marginal"` (default): GLS projection residuals \\(I - H)y\\.

  - `"cluster"`: Projection residuals after the fitted cluster
    contribution is removed. Only available for specialized multilevel
    models.

  - `"estimate"`: Sampling residual projection after the fitted
    random-effect contribution is removed (a BLUP-style target).

- max_samples:

  maximum posterior draws used for known-`V`
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  marginal covariance computations. Defaults to `Inf`. Finite values
  deterministically thin draws across posterior row order.

- ...:

  additional arguments (currently ignored)

## Value

A data frame with columns:

- `resid`: Posterior mean projection residuals at the requested depth

- `se`: Standard errors of the residuals

- `z`: Standardized residuals

## Details

This function returns a data frame with three columns matching the
output of
[`metafor::rstandard`](https://wviechtb.github.io/metafor/reference/residuals.rma.html):

- `resid`: Projection residuals: marginal depth uses \\(I - H)y\\;
  cluster and estimate depths retain their corresponding random-effect
  target. Known-`V` estimate depth uses the sampling residual
  projection. Residuals and SDs come from the same projection. These are
  not the posterior-mean residuals returned by
  `residuals(type = "outcome")`, which subtract the posterior mean of
  the fitted values and therefore carry the prior's shrinkage; the two
  agree as the effect prior becomes diffuse.

- `se`: Standard errors of the residuals

- `z`: Posterior mean of draw-wise standardized residuals

Internally standardized residuals are computed draw by draw by dividing
the observed residual by its corresponding residual standard error from
the hat matrix. The returned `z` column is the posterior mean of those
draw-wise ratios, so it is not generally equal to the displayed
`resid / se` summaries. For correlated known-`V` models, marginal
residuals use the GLS residual covariance and estimate-depth residuals
use the BLUP-style sampling residual projection. For a correctly
specified model, these residuals should approximately follow a standard
normal distribution. A design-exactly fitted row is set to zero only
when its covariance multiplier has no other partners; fully saturated
designs remain exactly zero under correlated covariance as well. Such a
row has residual and standard error zero and an undefined standardized
residual, so its `z` is `NA` (also when it is exactly fitted in only
some posterior draws); [`qqnorm()`](https://rdrr.io/r/stats/qqnorm.html)
omits these rows with a message.

For ensembles whose members differ in their PET and PEESE terms, each
draw uses the design of its member, as described in
[`hatvalues.brma`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md).

This function is only available for normal outcome models without
selection (weightfunction) bias adjustment. Selection models can use
[`rstudent.brma`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md),
which uses LOO-PIT. Internally standardized residuals are unavailable
for binomial and Poisson GLMMs because no count-likelihood
standardization has been implemented. Use
[`rstudent()`](https://rdrr.io/r/stats/influence.measures.html) (the
LOO-PIT residuals) for GLMMs.

For known-`V` models, the standardized residuals require applying the
inverse of the marginal covariance. A covariance that is too
ill-conditioned for that (a reciprocal condition number of about `2e-10`
or less in any posterior draw used) is unavailable: the call stops with
an error of class `RoBMA_known_v_precision_unavailable`, whose parent
class is `RoBMA_numerical_unavailable`, instead of returning a rounded
number. The condition carries the draw positions (`rows`) and their
`reciprocal_condition`. These positions refer to the current request
after deterministic thinning, including estimate-depth requests.

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
