# Regression Plot (Bubble Plot) for brma Object

`regplot.brma` creates a regression plot (also known as bubble plot) for
a fitted brma object with moderators. The plot displays observed effect
sizes against a moderator variable, with point sizes proportional to
study precision.

## Usage

``` r
# S3 method for class 'brma'
regplot(
  x,
  mod = NULL,
  pred = TRUE,
  ci = TRUE,
  pi = FALSE,
  si = FALSE,
  level = 95,
  at = NULL,
  digits = 2,
  transf = NULL,
  atransf = NULL,
  targs = NULL,
  refline = NULL,
  psize = NULL,
  plim = c(0.5, 3),
  by = NULL,
  legend = TRUE,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  sampling_bias = TRUE,
  sei = NULL,
  plot_type = "base",
  as_data = FALSE,
  ...
)
```

## Arguments

- x:

  a fitted brma object with moderators

- mod:

  index or name of the moderator variable to plot on the x-axis. If not
  specified and only one moderator is present, that moderator is used.
  If multiple moderators are present, this argument is required.

- pred:

  logical; whether to show the prediction line. Defaults to `TRUE`.

- ci:

  logical; whether to show credible interval bands. Defaults to `TRUE`.

- pi:

  logical; whether to show prediction interval bands. Defaults to
  `FALSE`.

- si:

  logical; whether to show sampling interval bands. Defaults to `FALSE`.
  The sampling interval shows the expected range of observed effect
  sizes, incorporating both heterogeneity (tau) and a representative
  level of sampling error (median SE across studies by default). When
  the model includes publication bias adjustments and
  `sampling_bias = TRUE`, the sampling interval incorporates the
  expected distortion from the selection process.

- level:

  numeric; credible/prediction interval level in percent. Defaults to
  `95`.

- at:

  numeric vector; for continuous moderators, values at which to evaluate
  the prediction. If not specified, uses a sequence across the observed
  range.

- digits:

  integer; number of decimal places for labels. Defaults to `2`.

- transf:

  function; transformation to apply to the y-axis (effect sizes).
  Defaults to `NULL` (no transformation).

- atransf:

  reserved for axis-label transformations. Currently not implemented and
  must be `NULL`.

- targs:

  reserved for additional transformation arguments. Currently not
  implemented and must be `NULL`.

- refline:

  numeric; position of horizontal reference line. Defaults to `NULL` (no
  reference line).

- psize:

  numeric vector or `NULL`; point sizes for each study. If scalar, it is
  recycled to all studies. If `NULL` (default), sizes are computed based
  on inverse sampling variance.

- plim:

  numeric vector of length 2; range for point sizes. Defaults to
  `c(0.5, 3)`.

- by:

  character; name of a moderator variable to use for separate
  lines/colors. Defaults to `NULL`. If omitted and the selected
  moderator enters exactly one two-way interaction, the other variable
  in the interaction is used automatically. Continuous `by` moderators
  are shown at mean - SD, mean, and mean + SD.

- legend:

  logical; whether to show legend when `by` is specified. Defaults to
  `TRUE`.

- xlab:

  character; x-axis label. Defaults to the moderator name.

- ylab:

  character; y-axis label. Defaults to "Observed Effect Size".

- xlim:

  numeric vector of length 2; x-axis limits. Defaults to data range.

- ylim:

  numeric vector of length 2; y-axis limits. Defaults to data range.

- sampling_bias:

  whether publication bias should be incorporated into plotted
  predictions and sampling intervals. Defaults to `TRUE`. For PET/PEESE
  models, this includes the expected regression bias in predictions. For
  selection models, sampling-bias adjustment applies to sampling
  intervals when `si = TRUE`; the mean prediction is unchanged.

- sei:

  single positive numeric value used as the reference standard error for
  sampling-bias and sampling-interval calculations. Defaults to the
  median observed standard error.

- plot_type:

  character; whether to use base R graphics (`"base"`) or ggplot2
  (`"ggplot"`). Defaults to `"base"`.

- as_data:

  logical; if `TRUE`, returns plot data instead of creating plot.
  Defaults to `FALSE`.

- ...:

  additional graphical arguments:

  main

  :   character string for plot title

  pch

  :   point symbol (default: 21)

  col

  :   point border color (default: "black")

  bg

  :   point fill color (default: "#A6A6A6")

  lcol

  :   line color (default: "black")

  lwd

  :   line width (default: 2)

  shade

  :   CI/PI/SI band shading (default: TRUE)

  col.ci

  :   CI band color (default: "gray70")

  col.pi

  :   PI band color (default: "gray85")

  col.si

  :   SI band color (default: "gray92")

  alpha.ci

  :   CI band transparency (default: 0.4)

  alpha.pi

  :   PI band transparency (default: 0.2)

  alpha.si

  :   SI band transparency (default: 0.15)

  jitter

  :   jitter amount for categorical moderators (default: 0.2)

  box.width

  :   box width for categorical interval summaries (default: 0.5)

## Value

`regplot.brma` returns `NULL` invisibly if `plot_type = "base"` or a
ggplot object if `plot_type = "ggplot"`. If `as_data = TRUE`, returns a
list with plot data components.

## Details

The regression plot (bubble plot) is a standard visualization for
meta-regression results. It displays:

- Observed effect sizes (y-axis) against moderator values (x-axis)

- Point sizes proportional to study precision (inverse variance)

- Prediction line showing the estimated regression relationship

- Confidence bands showing uncertainty in the mean prediction

- Optional prediction bands showing expected range of true effects

- Optional sampling interval bands showing expected range of observed
  outcomes

For continuous moderators, predictions are computed across the observed
range of the moderator. For categorical moderators (factors),
predictions are computed at each factor level with optional jittering of
points.

The `by` argument allows displaying separate regression lines for
different levels of a second moderator, useful for visualizing
interactions.

## See also

[`funnel.brma()`](https://fbartos.github.io/RoBMA/reference/funnel.md),
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

  fit <- brma(
    yi      = yi,
    vi      = vi,
    mods    = ~ ablat + year,
    data    = dat,
    measure = "RR"
  )
  regplot(fit, mod = "ablat")
  regplot(fit, mod = "year", pi = TRUE, si = TRUE)
  regplot(fit, mod = "ablat", plot_type = "ggplot")

  fit_cat <- brma(yi = yi, vi = vi, mods = ~ alloc, data = dat, measure = "RR")
  regplot(fit_cat)
}
} # }
```
