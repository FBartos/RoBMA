# Plots brma Object

`plot.brma` visualizes posterior (and prior) distribution a brma object.

## Usage

``` r
# S3 method for class 'brma'
plot(
  x,
  parameter,
  parameter_mods,
  parameter_scale,
  prior = FALSE,
  standardized_coefficients = FALSE,
  conditional = FALSE,
  output_measure = NULL,
  transform = NULL,
  plot_type = "base",
  dots_prior = NULL,
  ...
)
```

## Arguments

- x:

  a fitted `brma`, `BMA`, or `RoBMA` object.

- parameter:

  a parameter to be plotted. Defaults to `"mu"` for the effect size, or
  to the meta-regression intercept when moderators are present.
  Additional options are `"tau"`, `"rho"` for multilevel models,
  `"PET"`, `"PEESE"`, and `"omega"` or `"weightfunction"` for selection
  models. Use
  [`plot_pet_peese()`](https://fbartos.github.io/RoBMA/reference/plot_pet_peese.md)
  for PET/PEESE regression plots.

- parameter_mods:

  character. Moderator term to plot. Use `"intercept"` for the adjusted
  effect in meta-regression models.

- parameter_scale:

  character. Scale-regression term to plot. Use `"intercept"` for the
  heterogeneity intercept in location-scale models.

- prior:

  whether prior distribution should be added to figure. Defaults to
  `FALSE`.

- standardized_coefficients:

  whether to plot moderator and scale-regression coefficients on the
  standardized predictor scale. Defaults to `FALSE`.

- conditional:

  whether to plot the conditional posterior distribution for RoBMA
  product-space objects. Defaults to `FALSE`.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

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

`plot.brma` returns either `NULL` if `plot_type = "base"` or a `ggplot2`
object if `plot_type = "ggplot"`.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  plot(fit, parameter = "mu")
  plot(fit, parameter = "tau", prior = TRUE)
  plot(fit, parameter = "PET")
}
} # }

```
