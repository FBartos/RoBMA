# Funnel Plot for brma Object

`funnel.brma` creates a funnel plot for a fitted brma object. For
intercept-only models without scale regression, the default outcome mode
displays observed effect sizes against the fitted sampling distribution.
For models with location or scale moderators, the default residual mode
displays residuals against a standard-error funnel.

## Usage

``` r
# S3 method for class 'brma'
funnel(
  x,
  residual,
  type = "LOO-PIT",
  unit = "estimate",
  conditioning_depth = "marginal",
  sampling_heterogeneity = TRUE,
  sampling_bias = TRUE,
  max_samples = 10000,
  plot_type = "base",
  ...
)
```

## Arguments

- x:

  a fitted brma object

- residual:

  whether to use residual mode. Defaults to not specified, which means
  the function automatically determines the mode:

  - not specified: For intercept-only models without scale regression,
    displays observed effect sizes against the fitted sampling
    distribution funnel; for models with moderators or scale regression,
    automatically uses residual mode.

  - `FALSE`: explicitly requests outcome mode.

  - `TRUE`: explicitly requests residual mode, displaying residuals on
    the x-axis and using `type` to determine how those residuals are
    computed.

- type:

  the type of residuals to use when in residual mode. Options are:

  - `"LOO-PIT"` (alias: `"rstudent"`; default): Leave-one-out
    probability integral transform residuals returned by
    [`rstudent.brma`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md).

  - `"rstandard"`: Internally standardized residuals using
    [`rstandard.brma`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md).
    Only available for normal outcome models.

  - `"outcome"`: Raw outcome residuals from
    [`residuals.brma`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md)
    with `type = "outcome"`.

  Only used when funnel is in residual mode.

- unit:

  output unit for residual mode. Only `"estimate"` is implemented in
  this pass.

- conditioning_depth:

  residual conditioning depth for residual mode. Options are
  `"marginal"`, `"cluster"`, and `"estimate"`. The default LOO-PIT
  residual path is available only with marginal conditioning depth.

- sampling_heterogeneity:

  whether heterogeneity should be incorporated into the sampling
  distribution funnel. Defaults to `TRUE`. Only used in outcome mode and
  ignored in residual mode.

- sampling_bias:

  whether publication bias should be incorporated into the sampling
  distribution funnel. Defaults to `TRUE`. Only used when
  `residual = FALSE` or when automatic mode selects outcome mode.
  Ignored in residual mode. When `TRUE` and the model includes selection
  models (weightfunction), uses selected-normal quantiles. When `TRUE`
  and the model includes PET/PEESE, incorporates the expected skew from
  these regression adjustments.

- max_samples:

  maximum number of posterior samples used for model-averaged
  publication-bias funnel contours. Defaults to `10000`. Use `Inf` to
  use all posterior samples.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- ...:

  additional graphical arguments to customize the plot appearance:

  xlim, ylim

  :   numeric vectors of length 2 specifying axis limits

  xlab, ylab

  :   character strings for axis labels

  main

  :   character string for plot title (default: no title)

  pch

  :   point symbol (default: 21, filled circle). Use standard R pch
      values.

  col

  :   point border color (default: "black")

  bg

  :   point fill/background color (default: "#A6A6A6")

  cex

  :   point size multiplier for base graphics (default: 1)

  size

  :   point size for ggplot2 (default: 2)

  las

  :   axis-label style for base graphics (default: 1)

  back

  :   background region color (default: "grey"). Set to `NA` to
      suppress.

  shade

  :   funnel region fill color (default: "white"). Set to `NA` to
      suppress.

  lty

  :   line type for funnel edges and center line (default: "dotted")

  col.line

  :   color for funnel edge lines (default: "black")

  refline

  :   numeric override for the reference line. By default, residual mode
      uses 0, while outcome mode uses the center of the fitted sampling
      distribution, which may be curved when PET/PEESE bias adjustment
      is incorporated.

  col.refline

  :   color of vertical reference line (default: "black")

  as_data

  :   if `TRUE`, returns plot data instead of creating a plot

## Value

If `as_data = TRUE`, `funnel.brma` returns a list with the data used for
plotting, including the plotted points, funnel polygons, plotting
limits, labels, and reference line. Otherwise, it returns `NULL`
invisibly if `plot_type = "base"` or a ggplot object if
`plot_type = "ggplot"`.

## Details

The funnel plot has two modes. If `residual` is not specified, the mode
is chosen automatically from the fitted model: intercept-only models
without scale regression use outcome mode, whereas models with location
or scale moderators use residual mode.

**Outcome mode** (intercept-only models without scale regression):
Displays observed effect sizes on the x-axis and standard errors on the
y-axis. The reference line follows the center of the fitted sampling
distribution. When `sampling_bias = FALSE`, this center is the pooled
effect; when PET/PEESE bias adjustment is incorporated, the center line
can vary with the standard error. The funnel region represents the
central 95\\ region of the sampling distribution, optionally
incorporating heterogeneity and publication bias.

**Residual mode** (models with moderators or scale regression): Displays
residuals on the x-axis and standard errors on the y-axis. The funnel
region represents the central 95\\ \\N(0, \mathrm{SE}^2)\\. With
`type = "LOO-PIT"`, the plotted residuals and standard errors are the
raw-scale LOO predictive companions returned by
[`rstudent.brma`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md);
the PIT-normalized `z` values are used by
[`qqnorm.brma`](https://fbartos.github.io/RoBMA/reference/qqnorm.brma.md)
and influence diagnostics. With `type = "rstandard"`, the plotted values
are internally standardized residual companions from
[`rstandard.brma`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md).
With `type = "outcome"`, these are raw outcome residuals. Under a
correctly specified model, most points should fall within this region.

The `type` argument controls how residuals are computed in residual
mode. See
[`residuals.brma`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md)
for details on each type. The `sampling_heterogeneity` and
`sampling_bias` arguments are ignored in residual mode.

For GLMM models, observed effect sizes are computed from the raw
frequency data using formulas equivalent to
[`metafor::escalc`](https://wviechtb.github.io/metafor/reference/escalc.html).
Residual-mode GLMM funnels use approximate effect-size-scale
residual/PIT companions, not exact PIT diagnostics for the raw count
likelihood.

## See also

[`residuals.brma()`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md),
[`rstandard.brma()`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md),
[`rstudent.brma()`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE) &&
    requireNamespace("metafor", quietly = TRUE)) {
  data(dat.bcg, package = "metadat")
  dat <- metafor::escalc(
    measure = "RR",
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg
  )

  fit <- brma(yi = yi, vi = vi, data = dat, measure = "RR")
  funnel(fit)
  funnel(fit, pch = 19, col = "blue", bg = "lightblue")

  fit_reg <- brma(
    yi      = yi,
    vi      = vi,
    mods    = ~ ablat + year,
    data    = dat,
    measure = "RR"
  )
  fit_reg <- add_loo(fit_reg)
  funnel(fit_reg)
  funnel(fit_reg, type = "outcome")

  funnel_data <- funnel(fit, as_data = TRUE)
  funnel(fit, plot_type = "ggplot")
}
} # }
```
