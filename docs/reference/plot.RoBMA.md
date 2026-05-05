# Plots a fitted RoBMA object

`plot.RoBMA` allows to visualize different `"RoBMA"` object parameters
in various ways. See `type` for the different model types.

## Usage

``` r
# S3 method for class 'RoBMA'
plot(
  x,
  parameter = "mu",
  conditional = FALSE,
  plot_type = "base",
  prior = FALSE,
  output_scale = NULL,
  rescale_x = FALSE,
  show_data = TRUE,
  dots_prior = NULL,
  ...
)
```

## Arguments

- x:

  a fitted RoBMA object

- parameter:

  a parameter to be plotted. Defaults to `"mu"` (for the effect size).
  The additional options are `"tau"` (for the heterogeneity),
  `"weightfunction"` (for the estimated weightfunction), or
  `"PET-PEESE"` (for the PET-PEESE regression).

- conditional:

  whether conditional estimates should be plotted. Defaults to `FALSE`
  which plots the model-averaged estimates. Note that both
  `"weightfunction"` and `"PET-PEESE"` are always ignoring the other
  type of publication bias adjustment.

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

- rescale_x:

  whether the x-axis of the `"weightfunction"` should be re-scaled to
  make the x-ticks equally spaced. Defaults to `FALSE`.

- show_data:

  whether the study estimates and standard errors should be show in the
  `"PET-PEESE"` plot. Defaults to `TRUE`.

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

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)

### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
# the 'plot' function allows to visualize the results of a fitted RoBMA object, for example;
# the model-averaged effect size estimate
plot(fit, parameter = "mu")

# and show both the prior and posterior distribution
plot(fit, parameter = "mu", prior = TRUE)

# conditional plots can by obtained by specifying
plot(fit, parameter = "mu", conditional = TRUE)

# plotting function also allows to visualize the weight function
plot(fit, parameter = "weightfunction")

# re-scale the x-axis
plot(fit, parameter = "weightfunction", rescale_x = TRUE)

# or visualize the PET-PEESE regression line
plot(fit, parameter = "PET-PEESE")
} # }

```
