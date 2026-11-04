# Regression Plot (Bubble Plot) for brma Object

`regplot.brma` creates a regression plot (also known as bubble plot) for
a fitted brma object with location or scale moderators. The plot
displays observed effect sizes against a moderator variable, with point
sizes proportional to study precision.

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
  max_samples = 10000,
  plot_type = "base",
  as_data = FALSE,
  ...
)
```

## Arguments

- x:

  a fitted brma object with location or scale moderators

- mod:

  index or name of the moderator variable to plot on the x-axis. If not
  specified and only one moderator is present, that moderator is used.
  If multiple moderators are present, this argument is required.

- pred:

  logical; whether to show the prediction line. Defaults to `TRUE`.

- ci:

  logical; whether to show credible interval bands. Defaults to `TRUE`.
  See Details for the quantile estimator of the bands.

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

  numeric; credible/prediction interval level in percent, from `50` to
  `99.9`. Defaults to `95`.

- at:

  numeric vector; for continuous moderators, values at which to evaluate
  the prediction. If not specified, uses a sequence across the observed
  range. Continuous evaluation values are sorted with all associated
  predictions and intervals; duplicate evaluations are retained.

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
  recycled to all studies. If `NULL` (default), sizes are computed from
  the inverse sampling variance, mapped to `plim`. Equal precisions use
  the lower endpoint; zero SEs receive the established maximum priority
  when positive SEs are also present.

- plim:

  numeric vector of length 2; range for point sizes. Defaults to
  `c(0.5, 3)`.

- by:

  character; name of a moderator variable to use for separate
  lines/colors. Defaults to `NULL`. If omitted and the selected
  moderator enters exactly one two-way interaction, the other variable
  in the interaction is used automatically. Continuous `by` moderators
  are shown at mean - SD, mean, and mean + SD. Grouping by `sei` or `vi`
  supplies the pointwise precision when neither precision measure is
  plotted; an explicit `sei` reference is then not allowed. The
  companion precision measure (`sei` versus `vi`) is not selected
  automatically; an explicitly inconsistent companion grouping grid is
  unavailable.

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
  intervals when `si = TRUE`; the mean prediction is unchanged. Selected
  sampling intervals require a supported scalar selection law; joint
  configurations needing retained-source or multivariate-event
  integration are unavailable. Use `sampling_bias = FALSE` for
  bias-adjusted intervals. Conditioning every applicable selection
  source cancels positive weights and uses the ordinary Gaussian
  interval law. This Gaussian case permits zero pointwise SE; active
  selection kernels still require strictly positive pointwise SE.

- sei:

  single positive numeric value used as the reference standard error for
  sampling-bias and sampling-interval calculations. Defaults to the
  median observed standard error. When `mod = "sei"` or `mod = "vi"`,
  the plotted moderator supplies the pointwise standard error and an
  explicit `sei` override is not allowed. A grouping moderator
  `by = "sei"` or `by = "vi"` also supplies pointwise precision and
  disallows an explicit `sei` override. The companion precision column
  is synchronized with the plotted or grouped values. When precision is
  neither plotted nor grouped and no override is supplied, its existing
  representative predictor value supplies the precision.

- max_samples:

  maximum number of posterior samples used for prediction summaries and
  interval bands. Defaults to `10000`. Use `Inf` to use all posterior
  samples.

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

- Optional marginal prediction bands showing the expected range of new
  true effects

- Optional sampling interval bands showing expected range of observed
  outcomes

For continuous moderators, predictions are computed across the observed
range of the moderator. For categorical moderators (factors),
predictions are computed at each factor level with optional jittering of
points. The other predictors are held at a reference design: continuous
predictors of the location and scale formulas at their means and factors
at their reference levels, predictors of random slopes in the same way
(a random slope whose predictor is held at zero, such as a centered
predictor at its mean, adds no variance to the prediction and sampling
bands), and grouping variables at new row-unique levels. The line and
bands describe the regression at this reference design, not an average
over the observed predictor values, and the prediction and sampling
bands are pointwise rather than a joint band for one random regression
curve. For PET/PEESE models with `sampling_bias = TRUE`, the line and
its credible band include the expected bias at the reference standard
error (`sei`); `sampling_bias = FALSE` shows the bias-adjusted line. A
moderator that appears only in the scale formula produces a constant
mean prediction while prediction and sampling intervals reflect its
effect on heterogeneity.

The credible interval band of a continuous moderator (`ci = TRUE`) is
formed by pointwise Harrell-Davis quantiles (at the probabilities
implied by `level`) of the posterior draws of the mean prediction at
each grid value; they are computed by
[`BayesTools::harrell_davis_quantile()`](https://fbartos.github.io/BayesTools/reference/harrell_davis_quantile.html)
on the model scale, before `transf` is applied. The Harrell-Davis
quantiles estimate the same quantiles as the empirical quantiles of the
draws but smooth the kinks that the empirical quantiles show along the
moderator. The bands remain pointwise, and, because every draw
contributes, they can differ visibly from the empirical quantiles when
few draws are used (`max_samples` below a few hundred). Near the edge of
a point mass in the draws (e.g., draws of models that fix the plotted
quantity in a model-averaged posterior) the band blends the point mass
and the remaining draws, and extreme probabilities (e.g.,
`level = 99.9`) are not estimated more accurately than by the empirical
quantiles; see
[`BayesTools::harrell_davis_quantile()`](https://fbartos.github.io/BayesTools/reference/harrell_davis_quantile.html)
for these caveats. The credible intervals of categorical moderators use
the empirical quantiles of the draws of each level. Prediction and
sampling intervals are quantiles of a different distribution, the
posterior predictive mixture over the draws, and are computed by a
different estimator (mixture quantiles, not Harrell-Davis quantiles of
the draws).

Prediction intervals marginalize newly realized random effects. Sampling
intervals use the corresponding marginal response distribution by also
adding representative sampling error; evaluating a fitted moderator
value does not condition these bands on its fitted latent effect. For
marginalized known-`R` random-effect blocks, `pi` and `si` are
unavailable because the group-covariance multipliers belong to the
fitted layout and `R_new` is unsupported. Use `pi = FALSE, si = FALSE`
for the mean and credible-interval plot. For binomial and Poisson GLMMs,
prediction intervals retain this latent-effect interpretation, but
sampling intervals are only descriptive normal approximations on the
derived effect-size scale. They are not intervals from the fitted
discrete likelihood and do not use the prior-nuisance marginal GLMM
response target from
[`predict()`](https://rdrr.io/r/stats/predict.html).

The `by` argument allows displaying separate regression lines for
different levels of a second moderator, useful for visualizing
interactions. For
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
known-`V` models, prediction and sampling bands are pointwise
descriptive bands. They summarize row-wise marginal random-effect
variation and, for sampling intervals, a representative sampling
standard error. They are not full joint predictive regions for
correlated known-`V` outcomes.

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
