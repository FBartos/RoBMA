# Create Z-Curve Meta-Analytic Plot

Plots a fitted `zcurve_RoBMA` object, visualizing the z-curve, model
fit, extrapolation, and credible intervals.

## Usage

``` r
# S3 method for class 'zcurve_RoBMA'
plot(
  x,
  conditional = FALSE,
  plot_type = "base",
  probs = c(0.025, 0.975),
  max_samples = 500,
  plot_fit = TRUE,
  plot_extrapolation = TRUE,
  plot_CI = TRUE,
  plot_thresholds = TRUE,
  from = -6,
  to = 6,
  by.hist = 0.5,
  length.out.hist = NULL,
  by.lines = 0.05,
  length.out.lines = NULL,
  dots_hist = NULL,
  dots_extrapolation = NULL,
  dots_thresholds = NULL,
  ...
)
```

## Arguments

- x:

  A zcurve_RoBMA object to be plotted.

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

  Maximum number of posterior samples to use for plotting credible
  intervals. Defaults to 500.

- plot_fit:

  Should the model fit be included in the plot? Defaults to TRUE.

- plot_extrapolation:

  Should model extrapolation be included in the plot? Defaults to TRUE.

- plot_CI:

  Should credible intervals be included in the plot? Defaults to TRUE.

- plot_thresholds:

  Should significance thresholds be displayed in the plot? Defaults to
  TRUE.

- from:

  Lower bound of the z-value range for plotting. Defaults to -6.

- to:

  Upper bound of the z-value range for plotting. Defaults to 6.

- by.hist:

  Bin width for the histogram of observed z-values. Defaults to 0.5.

- length.out.hist:

  Number of bins for the histogram. If NULL, determined by by.hist.
  Defaults to NULL.

- by.lines:

  Step size for plotting model fit and extrapolation lines. Defaults to
  0.05.

- length.out.lines:

  Number of points for plotting lines. If NULL, determined by by.lines.
  Defaults to NULL.

- dots_hist:

  List of additional graphical parameters for the histogram. For base R:
  `border`, `col`, `density`, `angle`. For ggplot2: `color`, `fill`,
  `alpha`.

- dots_extrapolation:

  List of additional graphical parameters for the extrapolation
  line/ribbon. For base R lines: `lwd`, `col`, `lty`. For base R ribbon:
  `col` (will be alpha-blended), `border`. For ggplot2 lines:
  `linewidth`, `color`, `linetype`. For ggplot2 ribbon: `fill`, `alpha`.

- dots_thresholds:

  List of additional graphical parameters for the threshold lines. For
  base R: `col`, `lty`, `lwd`. For ggplot2: `color`, `linetype`,
  `linewidth`.

- ...:

  Additional graphical arguments for the main fit line and ribbon
  (passed to lines.zcurve_RoBMA), as well as basic plotting arguments.
  For base R lines: `lwd`, `col`, `lty`. For base R ribbon: `col` (will
  be alpha-blended), `border`. For ggplot2 lines: `linewidth`, `color`,
  `linetype`. For ggplot2 ribbon: `fill`, `alpha`. Basic plotting
  arguments (both base R and ggplot2): `xlab` (x-axis label, default:
  "Z-Statistic"), `ylab` (y-axis label, default: "Density"), `main`
  (plot title, default: ""), `ylim` (y-axis limits). For base R only:
  `xaxt` (x-axis type, default: "s"), `yaxt` (y-axis type, default:
  "s").

## Value

Returns `NULL` if `plot_type = "base"`, or a `ggplot2` object if
`plot_type = "ggplot2"`.

## See also

[`as_zcurve()`](reference/as_zcurve.md),
[`lines.zcurve_RoBMA()`](reference/lines.zcurve_RoBMA.md),
[`hist.zcurve_RoBMA()`](reference/hist.zcurve_RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n,
             study_names = Anderson2010$labels, algorithm = "ss")

zcurve_fit <- as_zcurve(fit)
plot(zcurve_fit)
} # }
```
