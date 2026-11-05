# Normal QQ Plot for brma Object

`qqnorm.brma` creates a normal QQ plot for a fitted brma object. The
plot displays sorted standardized residuals on the vertical axis against
the theoretical quantiles of a standard normal distribution on the
horizontal axis. Under a correctly specified model, the points should
fall approximately along the diagonal reference line.

## Usage

``` r
# S3 method for class 'brma'
qqnorm(
  y,
  type = "rstudent",
  unit = "estimate",
  conditioning_depth = "marginal",
  envelope = TRUE,
  conf_level = 95,
  bonferroni = FALSE,
  reps = 1000,
  smooth = TRUE,
  xlim,
  ylim,
  xlab,
  ylab,
  plot_type = "base",
  max_samples = Inf,
  ...
)
```

## Arguments

- y:

  a fitted brma object.

- type:

  the type of standardized residuals to use. Options are:

  - `"rstudent"` (alias: `"LOO-PIT"`; default): Leave-one-out
    probability integral transform residuals. Works for normal outcome
    models, including selection models, and, as a descriptive normal
    approximation, for binomial and Poisson GLMMs. This is the
    recommended type for Bayesian models as it properly accounts for
    estimation uncertainty and leverage. Note: requires that loo has
    been computed (see
    [`add_loo`](https://fbartos.github.io/RoBMA/dev/reference/add_loo.brma.md)).

  - `"rstandard"`: Internally standardized residuals using the hat
    matrix. Only available for normal outcome models without selection
    (weightfunction) bias adjustment.

- unit:

  output unit. Only `"estimate"` is implemented currently.

- conditioning_depth:

  residual conditioning depth. Options are `"marginal"`, `"cluster"`,
  and `"estimate"`. It does not apply to LOO-PIT residuals, whose
  metadata record the deletion unit and retained context instead; leave
  it omitted for those residuals.

- envelope:

  logical indicating whether to add a pseudo confidence envelope to the
  plot. The deterministic envelope uses Beta order-statistic quantiles
  transformed to the standard normal scale, providing a reference under
  a correctly specified model. Defaults to `TRUE`.

- conf_level:

  finite numeric value strictly between 0 and 100 specifying the
  confidence level for the envelope bounds. Defaults to `95`.

- bonferroni:

  logical indicating whether to apply Bonferroni correction to the
  envelope, so it can be regarded as a simultaneous confidence region
  for all \\k\\ residuals. Defaults to `FALSE`.

- reps:

  retained for compatibility with earlier simulation-based envelopes.
  The current envelope is deterministic and does not use this argument.

- smooth:

  logical indicating whether to smooth the envelope bounds using
  Friedman's SuperSmoother
  ([`supsmu`](https://rdrr.io/r/stats/supsmu.html)). Defaults to `TRUE`.

- xlim:

  x-axis limits. If not specified, limits are computed from data.

- ylim:

  y-axis limits. If not specified, limits are computed from data.

- xlab:

  title for the x-axis. If not specified, defaults to
  `"Theoretical Quantiles"`.

- ylab:

  title for the y-axis. If not specified, defaults to
  `"Sample Quantiles"`.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- max_samples:

  maximum posterior draws used to compute known-`V` internally
  standardized residuals. Defaults to `Inf`. Finite values select draws
  deterministically and are ignored for other residual types.

- ...:

  additional graphical arguments to customize the plot appearance:

  pch

  :   point symbol (default: 21, filled circle)

  col

  :   point border color (default: "black")

  bg

  :   point fill/background color (default: "#A6A6A6")

  cex

  :   point size multiplier for base graphics (default: 1)

  size

  :   point size for ggplot2 (default: 2)

  shade

  :   fill color for the envelope region (default: "grey90"). Set to
      `NA` to suppress shading.

  col.envelope

  :   color for envelope boundary lines (default: "grey50")

  col.refline

  :   color of the diagonal reference line (default: "black")

  main

  :   character string for plot title (default: no title)

  as_data

  :   if `TRUE`, returns plot data instead of creating the plot

## Value

For `plot_type = "base"`, returns an invisible list with components `x`
(theoretical quantiles) and `y` (sorted residuals). For
`plot_type = "ggplot"`, returns a ggplot object.

If `as_data = TRUE`, returns a list with all computed plot data
including: `x`, `y`, `points`, `refline`, `envelope`, `xlim`, `ylim`,
`xlab`, and `ylab`.

## Details

The normal QQ plot is a diagnostic tool for assessing whether the
standardized residuals follow a standard normal distribution, as
expected under a correctly specified model. Each point represents one
estimate, plotted at coordinates \\(\Phi^{-1}(p_i), z\_{(i)})\\ where
\\p_i\\ are the plotting positions and \\z\_{(i)}\\ are the sorted
standardized residuals. Each retained input row contributes one point,
irrespective of its likelihood weight. Weights affect fitting and the
residual calculation; they do not replicate points or multiply residual
values in this plot.

The default residual type is `"rstudent"` (LOO-PIT), which differs from
[`metafor::qqnorm`](https://wviechtb.github.io/metafor/reference/qqnorm.rma.html)
(which defaults to `"rstandard"`). LOO-PIT residuals are preferred for
normal outcome models because they properly account for estimation
uncertainty via leave-one-out cross-validation.

The confidence envelope is computed in closed form from the distribution
of normal order statistics. For the \\i\\-th sorted residual,
\\\Phi(Z\_{(i)})\\ follows a beta distribution with shape parameters
\\i\\ and \\k + 1 - i\\. Pointwise quantile bounds are transformed back
with \\\Phi^{-1}\\. When `smooth = TRUE`, bounds are smoothed using
Friedman's SuperSmoother. When `bonferroni = TRUE`, the confidence level
is adjusted for simultaneous inference across all \\k\\ comparisons.

For correlated known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
`"rstandard"` residuals, this envelope remains an independent
order-statistic reference band. It is descriptive and does not model
residual correlation induced by known sampling covariance.

For binomial and Poisson GLMMs the LOO-PIT residuals are a descriptive
normal approximation on the effect-size scale rather than a PIT of the
fitted discrete likelihood (see
[`residuals.brma`](https://fbartos.github.io/RoBMA/dev/reference/residuals.brma.md));
the default y-axis label says so.

## See also

[`rstudent.brma()`](https://fbartos.github.io/RoBMA/dev/reference/rstudent.brma.md),
[`rstandard.brma()`](https://fbartos.github.io/RoBMA/dev/reference/rstandard.brma.md),
[`residuals.brma()`](https://fbartos.github.io/RoBMA/dev/reference/residuals.brma.md),
[`funnel.brma()`](https://fbartos.github.io/RoBMA/dev/reference/funnel.md),
[`radial.brma()`](https://fbartos.github.io/RoBMA/dev/reference/radial.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )
  fit <- add_loo(fit)

  qqnorm(fit)
  qqnorm(fit, type = "rstandard")
  qqnorm(fit, envelope = FALSE)
  qqnorm(fit, plot_type = "ggplot")
  qqnorm(fit, pch = 19, col = "blue", shade = "lightblue")
}
} # }
```
