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

  whether to add to existing plot. Defaults to `FALSE`.

- plot_thresholds:

  whether to show significance threshold lines. Defaults to `TRUE`.

- dots_thresholds:

  graphical parameters for threshold lines (list).

- dots_hist:

  graphical parameters for histogram (list).

- dots_all:

  graphical parameters passed to all components (list).

- ...:

  additional graphical parameters.

## Value

`NULL` invisibly for base graphics, or a ggplot2 object.

## Details

Z-statistics are computed as effect size divided by standard error
(`yi / sei`). Histogram bins are adjusted to align with significance
thresholds when selection model priors are present.

## See also

[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md),
[`lines.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/lines.zplot_brma.md)
