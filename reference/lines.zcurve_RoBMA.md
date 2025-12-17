# Add Lines With Posterior Predictive Distribution of Z-Statistics

Adds lines to a plot of a fitted zcurve_RoBMA object. This function is
typically used to overlay additional information or model fits on an
existing plot.

## Usage

``` r
# S3 method for class 'zcurve_RoBMA'
lines(
  x,
  conditional = FALSE,
  plot_type = "base",
  probs = c(0.025, 0.975),
  max_samples = 500,
  extrapolate = FALSE,
  plot_CI = TRUE,
  from = -6,
  to = 6,
  by = 0.05,
  length.out = NULL,
  col = if (extrapolate) "blue" else "black",
  ...
)
```

## Arguments

- x:

  A RoBMA object

- conditional:

  whether conditional estimates should be plotted. Defaults to `FALSE`
  which plots the model-averaged estimates. Note that both
  `"weightfunction"` and `"PET-PEESE"` are always ignoring the other
  type of publication bias adjustment.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- probs:

  quantiles of the posterior samples to be displayed. Defaults to
  `c(.025, .975)`

- max_samples:

  Maximum number of samples from the posterior distribution that will be
  used for estimating z-curve estimates.

- extrapolate:

  Logical indicating whether to extrapolate values beyond the observed
  data range.

- plot_CI:

  Should credible intervals be included in the plot? Defaults to TRUE.

- from:

  Lower bound of the z-value range for plotting. Defaults to -6.

- to:

  Upper bound of the z-value range for plotting. Defaults to 6.

- by:

  Numeric value specifying the increment for the sequence.

- length.out:

  Optional integer specifying the desired length of the output sequence.

- col:

  Color of the plotted line.

- ...:

  Additional graphical parameters for the line and CI ribbon. For base R
  lines: `lwd`, `col`, `lty`. For base R ribbon: `alpha` (alpha
  transparency for the ribbon). For ggplot2 lines: `linewidth`, `color`,
  `linetype`. For ggplot2 ribbon: `alpha`.

## See also

[`as_zcurve()`](reference/as_zcurve.md),
[`plot.zcurve_RoBMA()`](reference/plot.zcurve_RoBMA.md),
[`hist.zcurve_RoBMA()`](reference/hist.zcurve_RoBMA.md)
