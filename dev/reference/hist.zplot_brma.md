# Histogram of Z-Statistics

Plots a histogram of observed z-values from the meta-analysis.

## Usage

``` r
# S3 method for class 'zplot_brma'
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
  dots_hist = NULL,
  dots_all = NULL,
  ...
)
```

## Arguments

- x:

  a zplot_brma object.

- plot_type:

  graphics system: `"base"` or `"ggplot"`. Defaults to `"base"`.

- from, to:

  z-value range for plotting. Defaults to `-6` and `6`.

- by:

  bin width. Defaults to 0.5.

- length.out:

  number of bins (alternative to `by`).

- add:

  whether to add to an existing base plot. Defaults to `FALSE`. For
  ggplot, combine the histogram and density layers with `+`;
  `add = TRUE` is unsupported.

- plot_thresholds:

  whether to show significance threshold lines. Defaults to `TRUE`.

- dots_thresholds:

  graphical parameters for threshold lines (list).

- dots_hist:

  graphical parameters for histogram (list).

- dots_all:

  graphical parameters passed to all components (list).

- ...:

  additional graphical parameters, such as `xlim` and `ylim`. Use
  `as_data = TRUE` to return histogram data without drawing.

## Value

`NULL` invisibly for base graphics, a ggplot2 object, or a data frame
with columns `x`, `density`, and `breaks` (bin widths) when
`as_data = TRUE`.

## Details

Z-statistics are computed as effect size divided by standard error
(`yi / sei`). Histogram bins are adjusted to align with significance
thresholds when selection model priors are present. The bins cover
`from` to `to`; when `by` does not divide that range, the last bin ends
at `to`. `xlim` and `ylim` set the displayed axis ranges without
changing the bins. Standalone ggplot histograms scale automatically to
added layers unless `ylim` is supplied.

## See also

[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.zplot_brma.md),
[`lines.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/lines.zplot_brma.md)
