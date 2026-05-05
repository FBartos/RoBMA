# Plot Zplot Results

Plots a zplot visualization showing the histogram of observed
z-statistics overlaid with model-implied densities.

## Usage

``` r
# S3 method for class 'zplot_brma'
plot(
  x,
  plot_type = "base",
  probs = c(0.025, 0.975),
  max_samples = 500,
  plot_fit = TRUE,
  plot_extrapolation = TRUE,
  plot_ci = TRUE,
  plot_thresholds = TRUE,
  from = -6,
  to = 6,
  by.hist = 0.5,
  length.out.hist = NULL,
  by.lines = 0.05,
  length.out.lines = NULL,
  dots_hist = NULL,
  dots_fit = NULL,
  dots_extrapolation = NULL,
  dots_thresholds = NULL,
  ...
)
```

## Arguments

- x:

  a zplot_brma object.

- plot_type:

  graphics system: `"base"` or `"ggplot"`. Defaults to `"base"`.

- probs:

  quantiles for credible intervals. Defaults to `c(.025, .975)`.

- max_samples:

  maximum posterior samples for density estimation. Defaults to 500.

- plot_fit:

  whether to show fitted density (with bias adjustments). Defaults to
  `TRUE`.

- plot_extrapolation:

  whether to show extrapolated density (bias removed). Defaults to
  `TRUE`.

- plot_ci:

  whether to show credible interval bands. Defaults to `TRUE`.

- plot_thresholds:

  whether to show significance threshold lines. Defaults to `TRUE`.

- from, to:

  z-value range for plotting. Defaults to `-6` and `6`.

- by.hist:

  bin width for histogram. Defaults to 0.5.

- length.out.hist:

  number of histogram bins (alternative to `by.hist`).

- by.lines:

  step size for density lines. Defaults to 0.05.

- length.out.lines:

  number of density points (alternative to `by.lines`).

- dots_hist:

  graphical parameters for histogram (list).

- dots_fit:

  graphical parameters for fit lines (list).

- dots_extrapolation:

  graphical parameters for extrapolation lines (list).

- dots_thresholds:

  graphical parameters for threshold lines (list).

- ...:

  additional graphical parameters passed to components.

## Value

`NULL` invisibly for base graphics, or a ggplot2 object.

## Details

The plot displays two density curves:

- Fit (black):

  Model-implied density including publication bias adjustments. This
  represents the expected distribution of z-statistics given the
  estimated selection process.

- Extrapolation (blue):

  Bias-corrected curve representing the hypothetical distribution
  without selective reporting, scaled by the expected suppressed studies
  under selection models. This curve is not normalized to integrate to
  one when selection implies missing studies.

## See also

[`hist.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/hist.zplot_brma.md),
[`lines.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/lines.zplot_brma.md),
[`summary.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/summary.zplot_brma.md)
