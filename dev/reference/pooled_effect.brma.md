# Pooled Effect Size for brma Objects

Computes the pooled effect and prediction interval at the average
expanded moderator design of a fitted brma object.

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
  location effect. This reference declares no new-publication standard
  error or selection event; explicit unadjusted
  [`predict()`](https://rdrr.io/r/stats/predict.html) rows instead
  target a selected new publication. Selection-model weighting affects
  response predictions, not this fixed-location posterior summary.

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
  [`predict.brma`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md)
  (`as_measure`, `V_new`, and `conditioning_depth`); wrapper arguments
  such as `newdata`, `type`, and `quiet` are fixed (supplying one is R's
  error for an argument matched twice, as for
  [`predict()`](https://rdrr.io/r/stats/predict.html)), and other
  unknown arguments are unused and trigger a warning.

## Value

A `brma_samples` object containing posterior samples. When printed,
displays a one-row summary table whose `PI` columns contain posterior
prediction quantiles. They are computed exactly from the posterior draws
of the pooled location and heterogeneity, as quantiles of the normal
mixture over the draws, so they involve no random numbers. Use
[`summary()`](https://rdrr.io/r/base/summary.html), or
`as.data.frame(summary(x))`, to obtain the table; `as.data.frame(x)`
returns the posterior draws. The samples can be converted to posterior
draws formats using
[`as_draws()`](https://fbartos.github.io/RoBMA/dev/reference/as_draws.brma.md).

## Details

This function is a convenience wrapper around
`predict.brma(..., type = "terms", newdata = TRUE, bias_adjusted = TRUE, quiet = TRUE)`
that adds the prediction interval columns.

The fixed-location model is evaluated at the average expanded moderator
design. Because the location predictor is linear in its design, this
equals averaging its fitted-design posterior draws directly.

For meta-regression models, the pooled effect averages the effect size
estimate across moderator levels proportionately to the levels observed
in the data. This provides an estimate representative of the sample of
studies.

For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
models and models with factor-specific or random-slope terms, the
average expanded design can be synthetic and need not correspond to a
single realizable outcome or study. The pooled effect and its prediction
interval then describe that linear composite and should be interpreted
only when the composite has a meaningful scientific interpretation. Use
`predict(..., newdata = ...)` for a specified realizable design.

For models without moderators, this returns the single mu parameter. The
prediction interval is marginal over one new true effect at the same
average location and random-effect design, with the scale given by the
mean of the observation-specific tau, as in
[`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_heterogeneity.md).
It therefore uses `pooled_heterogeneity(component = "total")`, rather
than the observed-design RMS heterogeneity reported by
[`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md).
For known-\\R\\ random effects this is a unit-kernel reference interval;
use [`predict()`](https://rdrr.io/r/stats/predict.html) for an interval
tied to a fitted group and its known covariance multiplier.

## See also

[`predict.brma()`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md),
[`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_heterogeneity.md),
[`blup()`](https://fbartos.github.io/RoBMA/dev/reference/blup.md)

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
