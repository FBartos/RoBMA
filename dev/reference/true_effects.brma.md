# True Effects for brma Objects

Computes the estimated true effects (theta) for a fitted brma object.
This is an alias for
[`blup.brma`](https://fbartos.github.io/RoBMA/dev/reference/blup.brma.md).

## Usage

``` r
# S3 method for class 'brma'
true_effects(
  object,
  bias_adjusted = FALSE,
  output_measure = NULL,
  transform = NULL,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- object:

  a fitted brma object

- bias_adjusted:

  whether to adjust for publication bias. Defaults to `FALSE`, which
  returns estimates including publication bias effects (i.e., what we
  expect the true effects to be given the biased observations). Set to
  `TRUE` to obtain bias-corrected estimates.

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

- ...:

  additional arguments passed to
  [`predict.brma`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md);
  wrapper arguments such as `newdata`, `type`, `quiet`,
  `output_measure`, and `transform` are fixed by this method (supplying
  one is R's error for an argument matched twice, as for
  [`predict()`](https://rdrr.io/r/stats/predict.html)), and other
  unknown arguments are unused and trigger a warning.

## Value

A `brma_samples` object containing posterior draws of BLUP or
empirical-Bayes true-effect summaries with one column per estimate. For
existing normal data, these are conditional BLUP means, not simulated
latent-effect draws. When printed, displays a summary table. Use
[`summary()`](https://rdrr.io/r/base/summary.html), or
`as.data.frame(summary(x))`, to obtain the summary table;
`as.data.frame(x)` returns the posterior draws. The samples can be
converted to posterior draws formats using
[`as_draws()`](https://fbartos.github.io/RoBMA/dev/reference/as_draws.brma.md).

## Details

This function is identical to
[`blup.brma`](https://fbartos.github.io/RoBMA/dev/reference/blup.brma.md).
See that function for full details on how true effects are computed.

## See also

[`blup.brma()`](https://fbartos.github.io/RoBMA/dev/reference/blup.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_effect.md),
[`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_heterogeneity.md)

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

  true_effects(fit)
}
} # }
```
