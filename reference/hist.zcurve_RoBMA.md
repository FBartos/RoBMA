# Create Histogram of Z-Statistics

Plots a histogram of observed z-values for a fitted `zcurve_RoBMA`
object, with options to customize the plotting range, bin width, and
display of significance thresholds.

## Usage

``` r
# S3 method for class 'zcurve_RoBMA'
hist(
  x,
  plot_type = "base",
  from = -6,
  to = 6,
  by = 0.5,
  length.out = NULL,
  add = FALSE,
  plot_thresholds = TRUE,
  dots_thresholds = NULL,
  ...
)
```

## Arguments

- x:

  A `zcurve_RoBMA` object containing the fitted model.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- from:

  Lower bound of the z-value range for plotting. Defaults to -6.

- to:

  Upper bound of the z-value range for plotting. Defaults to 6.

- by:

  Numeric value specifying the bin width for the histogram. Defaults to
  0.5.

- length.out:

  Optional integer specifying the number of bins. If NULL, determined by
  `by`. Defaults to NULL.

- add:

  Logical; if TRUE, adds histogram bars to an existing plot without
  creating a new canvas. Only applies to base R graphics. Defaults to
  FALSE.

- plot_thresholds:

  Logical; should significance thresholds be displayed on the plot?
  Defaults to TRUE.

- dots_thresholds:

  List of additional graphical parameters for the threshold lines. For
  base R: `col`, `lty`, `lwd`. For ggplot2: `color`, `linetype`,
  `linewidth`.

- ...:

  Additional graphical parameters for the histogram and basic plotting
  arguments. For base R histogram: `border`, `col`, `density`, `angle`.
  For ggplot2 histogram: `color`, `fill`, `alpha`. Basic plotting
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
[`plot.zcurve_RoBMA()`](reference/plot.zcurve_RoBMA.md),
`hist.zcurve_RoBMA()`
