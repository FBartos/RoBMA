# Pooled Heterogeneity for brma Objects

Computes the pooled (aggregated) heterogeneity estimate (tau) from a
fitted brma object by averaging across the scale model matrix.

## Usage

``` r
# S3 method for class 'brma'
pooled_heterogeneity(object, probs = c(0.025, 0.975), conditional = FALSE, ...)
```

## Arguments

- object:

  a fitted brma object

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- conditional:

  whether to return the pooled heterogeneity conditional on the
  heterogeneity component for RoBMA product-space objects. Defaults to
  `FALSE`.

- ...:

  additional arguments passed to
  [`predict.brma`](https://fbartos.github.io/RoBMA/reference/predict.brma.md);
  wrapper arguments such as `newdata`, `type`, and `quiet` are fixed.

## Value

A `brma_samples` object containing posterior samples. When printed,
displays a summary table. Use
[`summary()`](https://rdrr.io/r/base/summary.html) to obtain the summary
table directly. The samples can be converted to posterior draws formats
using
[`as_draws()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md).

## Details

This function is a convenience wrapper around
`predict.brma(..., type = "terms.scale", newdata = TRUE)`.

For location-scale models (with scale regression), the pooled
heterogeneity averages tau across the scale model matrix proportionately
to the levels observed in the data.

For models without scale regression, this returns the single tau
parameter.

For multilevel (3-level) models, the returned tau is the total
heterogeneity: `tau = sqrt(tau_within^2 + tau_between^2)`.

## See also

[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md),
[`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md)

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

  pooled_heterogeneity(fit)
}
} # }
```
