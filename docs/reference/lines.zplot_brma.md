# Add Zplot Density Lines

Adds model-implied density lines to an existing zplot.

## Usage

``` r
# S3 method for class 'zplot_brma'
lines(
  x,
  plot_type = "base",
  probs = c(0.025, 0.975),
  max_samples = 500,
  plot_ci = TRUE,
  extrapolate = FALSE,
  from = -6,
  to = 6,
  by = 0.05,
  length.out = NULL,
  col = "black",
  as_data = FALSE,
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

  maximum posterior samples for density. Defaults to 500.

- plot_ci:

  whether to show credible interval bands. Defaults to `TRUE`.

- extrapolate:

  whether to remove bias adjustments. Defaults to `FALSE`.

- from, to:

  z-value range for density. Defaults to `-6` and `6`.

- by:

  step size for density points. Defaults to 0.05.

- length.out:

  number of density points (alternative to `by`).

- col:

  line color. Defaults to `"black"`.

- as_data:

  whether to return data instead of plotting. Defaults to `FALSE`.

- ...:

  additional graphical parameters.

## Value

`NULL` invisibly for base graphics, ggplot2 layers for ggplot, or a data
frame with columns `x`, `y`, `y_lCI`, `y_uCI` if `as_data = TRUE`.

## Details

When `extrapolate = FALSE`, the density includes all bias adjustments
(PET/PEESE regression, selection weights) representing the fitted model.
When `extrapolate = TRUE`, bias adjustments are removed to show the
hypothetical distribution without publication bias. Under selection
models, the curve is scaled by inverse selection probability and need
not integrate to one.

## See also

[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md),
[`hist.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/hist.zplot_brma.md)
