# Frequentist-Style and Bayesian Funnel Plots for brma Objects

`funnel.brma` creates fast plug-in funnel contours. Location and bias
parameters are fixed at their within-model posterior means;
heterogeneity uses the root mean square SD, \\\sqrt{E(\tau^2)}\\. For
model-averaged fits, complete joint models retain their posterior model
probabilities and their plug-in sampling CDFs are averaged before the
contours are inverted. This makes the ordinary normal-model construction
directly comparable to a conventional plug-in funnel such as
`metafor::funnel(..., addtau2 = TRUE)`.

`bfunnel.brma` instead averages the sampling CDF over posterior draws
before inversion. It therefore incorporates continuous parameter
uncertainty as well as model uncertainty.

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

bfunnel(x, ...)

# Default S3 method
bfunnel(x, ...)

# S3 method for class 'brma'
bfunnel(
  x,
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

  - not specified: Displays outcomes only for intercept-only models with
    one common marginal heterogeneity distribution. Otherwise it
    automatically uses LOO-PIT residual mode.

  - `FALSE`: explicitly requests outcome mode. This is available only
    for intercept-only models with common marginal heterogeneity.

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
    with `type = "outcome"`. Its funnel bands are descriptive
    sampling-error-only reference bands.

  Only used when funnel is in residual mode.

- unit:

  output unit for residual mode. Only `"estimate"` is currently
  implemented.

- conditioning_depth:

  conditioning depth for non-LOO residual mode. Options are
  `"marginal"`, `"cluster"`, and `"estimate"`. Omit this argument for
  `type = "LOO-PIT"` (the default) or its `"rstudent"` alias: the
  estimate-deletion context determines the LOO predictive target.

- sampling_heterogeneity:

  whether heterogeneity should be incorporated into the sampling
  distribution funnel. Defaults to `TRUE`. Only used in outcome mode and
  ignored in residual mode.

- sampling_bias:

  whether publication bias should be incorporated into the sampling
  distribution funnel. Defaults to `TRUE`. Only used when
  `residual = FALSE` or when automatic mode selects outcome mode.
  Ignored in residual mode. When `TRUE` and the model includes selection
  models (weightfunction), uses selected-normal quantiles where the
  fitted source and publication structure admit a scalar calculation.
  Unsupported joint-selection configurations fail explicitly; use
  `sampling_bias = FALSE` for bias-adjusted contours. Exact-zero
  integrated selected cells require an existing positive-weight
  cancellation certificate; otherwise the selected contour is
  unavailable. `sampling_bias = FALSE` remains an explicitly unselected
  reference. Ordinary and constant-positive step-selection rows are
  evaluated as exact Gaussian components before conditioned-source
  quadrature. Remaining selected rows are unavailable when the required
  odd quadrature order exceeds the existing 401-node cap or is
  nonfinite; the call then stops with an error of class
  `RoBMA_funnel_quadrature_rejected`, whose parent class is
  `RoBMA_numerical_unavailable`, carrying the standard deviation
  `ratio`, the `required_order` and the `maximum` order. `max_samples`
  controls posterior rows and does not increase this quadrature budget.
  When `TRUE` and the model includes PET/PEESE, incorporates the
  expected skew from these regression adjustments.

- max_samples:

  maximum number of posterior draws used by `bfunnel()`. Defaults to
  `10000`; use `Inf` for all draws. `funnel()` uses all draws to
  estimate conditional plug-in values and posterior model probabilities,
  then evaluates only one CDF per joint model.

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

If `as_data = TRUE`, `funnel.brma` and `bfunnel.brma` return a list with
the data used for plotting, including the plotted points, funnel
polygons, plotting limits, labels, reference line, and outcome-funnel
estimand. GLMM outcome-mode data of `funnel.brma` also contain an
`approximation` record identifying the descriptive normal effect-size
approximation and stating that it is not fitted discrete- likelihood
coverage (`bfunnel.brma` is available only for normal outcome models).
Otherwise, it returns `NULL` invisibly if `plot_type = "base"` or a
ggplot object if `plot_type = "ggplot"`.

## Details

The ordinary `funnel()` has outcome and residual modes. If `residual` is
not specified, outcomes are used only when the fitted model is
intercept-only and its marginal heterogeneity is identical across
studies. The latter check uses fitted row-marginal random-effect
variances, including complete \\ZGZ'\\ contributions for
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
random- formula models. All other models use residual mode.

**Outcome mode** (intercept-only models with common heterogeneity):
Displays observed effect sizes on the x-axis and standard errors on the
y-axis. The reference line follows the center of the fitted sampling
distribution. When `sampling_bias = FALSE`, this center is the plug-in
mixture median (the pooled effect in a single model); when PET/PEESE
bias adjustment is incorporated, the center line can vary with the
standard error. The funnel region represents the central 95\\ region of
the plug-in sampling distribution, optionally incorporating
heterogeneity and publication bias. Location and bias parameters use
posterior means, and heterogeneity uses the square root of the posterior
mean variance, separately within every complete joint model. Integrated
and conditioned variance components are averaged on the variance scale.
The resulting model-specific CDFs are averaged with posterior model
probabilities and only then inverted; neither parameters nor contour
endpoints are averaged across models. For normal models, these contours
use the marginal response target: heterogeneity represents newly
realized latent effects and fitted study effects are not retained. For
correlated known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, outcome-mode funnels use the diagonal of `V` for Gaussian
marginal contours. The contour describes a new, independent estimate at
each hypothetical standard error, so the fitted rows' correlations do
not enter it: under `weight_rule = "product"` such an estimate forms its
own dependency block and its selected law is the scalar selected normal
at that standard error, mixed over the population law of any conditioned
source. Nonconstant joint selection under a best weight rule instead
weighs an estimate against the others of its publication; it requires a
full publication design at each hypothetical standard error and remains
unavailable in the scalar contour calculation. Conditioning every
applicable selection source cancels positive weights and therefore
retains the Gaussian contours. The observed estimates' own selected
marginals inside a correlated block differ from this new-study contour;
use [`zplot()`](https://fbartos.github.io/RoBMA/reference/zplot.brma.md)
to display them at the fitted publication designs. Residual mode follows
the fitted estimate-unit residual target.

`bfunnel()` is available only for intercept-only normal models with
common marginal heterogeneity. It has no residual mode: explicitly
supplied `residual`, `type`, `unit`, or `conditioning_depth` arguments
are refused. It averages the same sampling CDF over posterior draws and
then inverts the averaged CDF. Thus it incorporates posterior
uncertainty in the effect, heterogeneity, PET/PEESE coefficients, and
selection parameters in addition to model uncertainty. For models
without a common outcome-scale target, use the default LOO-PIT residual
funnel from `funnel()`.

**Residual mode** (models without a common outcome-scale target):
Displays residuals on the x-axis and standard errors on the y-axis. The
funnel region is bounded by the 0.025 and 0.975 normal quantiles times
the displayed standard error. It is a descriptive reference band, not a
posterior-predictive coverage interval. With `type = "LOO-PIT"`, the
plotted x-coordinate is the PIT-normalized `z` value returned by
[`rstudent.brma`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md)
multiplied by its LOO predictive standard deviation. Consequently, each
plotted `x / se` is exactly its LOO-PIT `z`. This reduces to the
ordinary raw deleted residual when the LOO predictive distribution is
normal and extends the same funnel geometry to non-normal continuous
predictive mixtures. With `type = "rstandard"`, the plotted values are
internally standardized residual companions from
[`rstandard.brma`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md).
With `type = "outcome"`, these are raw outcome residuals and the band
uses sampling error only; it does not incorporate fitted heterogeneity.

The `type` argument controls how residuals are computed in residual
mode. See
[`residuals.brma`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md)
for details on each type. The `sampling_heterogeneity` and
`sampling_bias` arguments are ignored in residual mode.

For GLMM models, outcome-mode observed effect sizes are computed from
the raw frequency data using formulas equivalent to
[`metafor::escalc`](https://wviechtb.github.io/metafor/reference/escalc.html).
Their contours use the corresponding continuity-corrected effect-size
estimates and approximate standard errors in a normal reference
calculation. The default x-axis label and a warning identify this as a
descriptive approximation; it is not a coverage interval from the fitted
binomial or Poisson likelihood and does not use the prior-nuisance
marginal GLMM response target from
[`predict()`](https://rdrr.io/r/stats/predict.html). Exact
discrete-likelihood coverage cannot be indexed by a scalar standard
error alone because it also depends on arm totals or exposures and
nuisance rates. LOO-PIT residual mode of GLMMs uses the descriptive
normal approximation on the effect-size scale of
[`residuals.brma`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md)
(the fitted predictive distribution is discrete); the default x-axis
label identifies it. The sampling-error-only bands of the residual
funnel are not exact diagnostics for the fitted count likelihood.

Outcome-mode plot data carry the canonical `RoBMA_target` attribute. Its
`reference_law` identifies the independently published new-estimate
product-rule law; `estimand` separately records plug-in versus
posterior-predictive averaging. Declared source roles retain their
provenance.

## See also

[`residuals.brma()`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md),
[`rstandard.brma()`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md),
[`rstudent.brma()`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md),
[`zplot.brma()`](https://fbartos.github.io/RoBMA/reference/zplot.brma.md)

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
  bfunnel(fit)
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
