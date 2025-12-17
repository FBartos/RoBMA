# Models plot for a RoBMA object

`plot_models` plots individual models' estimates for a `"RoBMA"` object.

## Usage

``` r
plot_models(
  x,
  parameter = "mu",
  conditional = FALSE,
  output_scale = NULL,
  plot_type = "base",
  order = "decreasing",
  order_by = "model",
  ...
)
```

## Arguments

- x:

  a fitted RoBMA object

- parameter:

  a parameter to be plotted. Defaults to `"mu"` (for the effect size).
  The additional option is `"tau"` (for the heterogeneity).

- conditional:

  whether conditional estimates should be plotted. Defaults to `FALSE`
  which plots the model-averaged estimates. Note that both
  `"weightfunction"` and `"PET-PEESE"` are always ignoring the other
  type of publication bias adjustment.

- output_scale:

  transform the effect sizes and the meta-analytic effect size estimate
  to a different scale. Defaults to `NULL` which returns the same scale
  as the model was estimated on.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- order:

  how the models should be ordered. Defaults to `"decreasing"` which
  orders them in decreasing order in accordance to `order_by` argument.
  The alternative is `"increasing"`.

- order_by:

  what feature should be use to order the models. Defaults to `"model"`
  which orders the models according to their number. The alternatives
  are `"estimate"` (for the effect size estimates), `"probability"` (for
  the posterior model probability), and `"BF"` (for the inclusion Bayes
  factor).

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `main`, `xlim`, `ylim` to adjust the line thickness,
  line type, line color, fill color, x-label, y-label, title, x-axis
  range, and y-axis range respectively.

## Value

`plot_models` returns either `NULL` if `plot_type = "base"` or an object
object of class 'ggplot2' if `plot_type = "ggplot2"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)

### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
# the plot_models function creates a plot for of the individual models' estimates, for example,
# the effect size estimates from the individual models can be obtained with
plot_models(fit)

# and effect size estimates from only the conditional models
plot_models(fit, conditional = TRUE)
} # }

```
