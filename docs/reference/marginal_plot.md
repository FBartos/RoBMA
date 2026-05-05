# Plots marginal estimates of a fitted RoBMA regression object

`marginal_plot` allows to visualize prior and posterior distributions of
marginal estimates of a RoBMA regression model.

## Usage

``` r
marginal_plot(
  x,
  parameter,
  conditional = FALSE,
  plot_type = "base",
  prior = FALSE,
  output_scale = NULL,
  dots_prior = NULL,
  ...
)
```

## Arguments

- x:

  a fitted RoBMA regression object

- parameter:

  regression parameter to be plotted

- conditional:

  whether conditional marginal estimates should be plotted. Defaults to
  `FALSE` which plots the model-averaged estimates.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- prior:

  whether prior distribution should be added to figure. Defaults to
  `FALSE`.

- output_scale:

  transform the effect sizes and the meta-analytic effect size estimate
  to a different scale. Defaults to `NULL` which returns the same scale
  as the model was estimated on.

- dots_prior:

  list of additional graphical arguments to be passed to the plotting
  function of the prior distribution. Supported arguments are `lwd`,
  `lty`, `col`, and `col.fill`, to adjust the line thickness, line type,
  line color, and fill color of the prior distribution respectively.

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `main`, `xlim`, `ylim` to adjust the line thickness,
  line type, line color, fill color, x-label, y-label, title, x-axis
  range, and y-axis range respectively.

## Value

`plot.RoBMA` returns either `NULL` if `plot_type = "base"` or an object
object of class 'ggplot2' if `plot_type = "ggplot2"`.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)
