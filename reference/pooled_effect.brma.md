# Pooled Effect Size for brma Objects

Computes the pooled (aggregated) effect size estimate from a fitted brma
object by averaging across the moderator model matrix.

## Usage

``` r
# S3 method for class 'brma'
pooled_effect(
  object,
  bias_adjusted = TRUE,
  output_measure = NULL,
  transform = NULL,
  probs = c(0.025, 0.975),
  conditional = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted brma object

- bias_adjusted:

  whether to adjust for publication bias. Defaults to `TRUE`. For
  PET/PEESE models this removes the regression bias term from the pooled
  location effect. Selection-model weighting affects response
  predictions, not this `type = "terms"` wrapper.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- probs:

  quantiles of the posterior distribution to be displayed. Defaults to
  `c(.025, .975)` for 95% credible intervals.

- conditional:

  whether to return the pooled effect conditional on the effect
  component for RoBMA product-space objects. Defaults to `FALSE`.

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
`predict.brma(..., type = "terms", newdata = TRUE, bias_adjusted = TRUE, quiet = TRUE)`.

For meta-regression models, the pooled effect averages the effect size
estimate across moderator levels proportionately to the levels observed
in the data. This provides an estimate representative of the sample of
studies.

For models without moderators, this returns the single mu parameter.

## See also

[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md),
[`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md),
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

  pooled_effect(fit)
}
} # }
```
