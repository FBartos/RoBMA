# Funnel plot for a RoBMA object

`funnel` creates a funnel plot for a `"RoBMA"` object. Only available
for normal-normal models estimated using the spike-and-slab algorithm
(i.e., `algorithm = "ss"`). This function uses several simplifications
to visualize the sampling distribution, see Details for more
information.

## Usage

``` r
funnel(
  x,
  conditional = FALSE,
  plot_type = "base",
  output_scale = NULL,
  incorporate_heterogeneity = TRUE,
  incorporate_publication_bias = TRUE,
  max_samples = 500,
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

- incorporate_heterogeneity:

  Whether heterogeneity should be incorporated into the sampling
  distribution. Defaults to `TRUE`.

- incorporate_publication_bias:

  Whether publication bias should be incorporated into the sampling
  distribution. Defaults to `TRUE`.

- max_samples:

  Maximum number of samples from the posterior distribution that will be
  used for estimating the funnel plot under publication bias. Defaults
  to `500`.

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `main`, `xlim`, `ylim` to adjust the line thickness,
  line type, line color, fill color, x-label, y-label, title, x-axis
  range, and y-axis range respectively.

## Value

`funnel` returns either `NULL` if `plot_type = "base"` or an object
object of class 'ggplot2' if `plot_type = "ggplot2"`.

## Details

The `funnel` function differs from the corresponding
[funnel](https://wviechtb.github.io/metafor/reference/funnel.html)
function in two regards; 1) The heterogeneity is by default incorporated
into the model and 2) the sampling distribution under publication bias
is approximate by sampling from the estimated weighted normal
distribution under the pooled effect. This approximation may distort the
true sampling distribution (but that would be impossible) to visualize
in the usual funnel plot style anyway?.

The sampling distribution is drawn under the mean effect size and
heterogeneity estimates (the uncertainty about those values is not
incorporated into the figure).

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n,
             study_names = Anderson2010$labels, algorithm = "ss")

funnel(fit)
} # }
```
