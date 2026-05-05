# Forest plot for a RoBMA object

`forest` creates a forest plot for a `"RoBMA"` object.

## Usage

``` r
forest(
  x,
  conditional = FALSE,
  plot_type = "base",
  output_scale = NULL,
  order = NULL,
  ...
)
```

## Arguments

- x:

  a fitted RoBMA object

- conditional:

  whether conditional estimates should be plotted. Defaults to `FALSE`
  which plots the model-averaged estimates. Note that both
  `"weightfunction"` and `"PET-PEESE"` are always ignoring the other
  type of publication bias adjustment.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- output_scale:

  transform the effect sizes and the meta-analytic effect size estimate
  to a different scale. Defaults to `NULL` which returns the same scale
  as the model was estimated on.

- order:

  order of the studies. Defaults to `NULL` - ordering as supplied to the
  fitting function. Studies can be ordered either `"increasing"` or
  `"decreasing"` by effect size, or by labels `"alphabetical"`.

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `main`, `xlim`, `ylim` to adjust the line thickness,
  line type, line color, fill color, x-label, y-label, title, x-axis
  range, and y-axis range respectively.

## Value

`forest` returns either `NULL` if `plot_type = "base"` or an object
object of class 'ggplot2' if `plot_type = "ggplot2"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)

### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
# the forest function creates a forest plot for a fitted RoBMA object, for example,
# the forest plot for the individual studies and the model-averaged effect size estimate
forest(fit)

# the conditional effect size estimate
forest(fit, conditional = TRUE)

# or transforming the effect size estimates to Fisher's z
forest(fit, output_scale = "fishers_z")
} # }
```
