# Radial (Galbraith) Plot for brma Object

`radial.brma` creates a radial (Galbraith) plot for a fitted brma
object. The plot displays each study as a point with precision on the
x-axis and the standardized effect size on the z-axis. Without
centering, a line from the origin through any point has slope equal to
the observed effect size. With centering, slopes are relative to the
pooled estimate. An arc on the right side maps standardized values back
to the effect size scale.

## Usage

``` r
# S3 method for class 'brma'
radial(
  x,
  center = FALSE,
  xlim,
  zlim,
  xlab,
  zlab,
  atz,
  aty,
  steps = 7,
  level = 95,
  digits = 2,
  transf,
  targs,
  plot_type = "base",
  ...
)

# S3 method for class 'brma'
galbraith(x, ...)
```

## Arguments

- x:

  a fitted brma object. Must be an intercept-only model (no moderators
  or scale regression).

- center:

  logical indicating whether to center the plot at the pooled estimate.
  When `TRUE`, the pooled estimate is shifted to zero on the z-axis.
  Defaults to `FALSE`.

- xlim:

  x-axis limits. If not specified, limits are computed from data.

- zlim:

  z-axis limits. If not specified, limits are computed from data.

- xlab:

  title for the x-axis. If not specified, a default label is used.

- zlab:

  title for the z-axis. If not specified, a default label is used.

- atz:

  numeric vector of positions for tick marks on the z-axis (left
  vertical axis). If not specified, positions are chosen automatically.

- aty:

  numeric vector of effect size values to mark on the y-axis arc (right
  side). If not specified, values are chosen automatically.

- steps:

  integer specifying the number of tick marks on the y-axis arc when
  `aty` is not specified. Defaults to `7`.

- level:

  numeric value between 0 and 100 specifying the confidence level for
  the confidence band and pooled-effect CI arc. Defaults to `95`.

- digits:

  integer specifying the number of decimal places for y-axis arc labels.
  Defaults to `2`.

- transf:

  optional transformation function applied to the y-axis arc labels
  (e.g., `transf = exp` when effect sizes are on a log scale).

- targs:

  optional list passed to `transf` as its second argument,
  `transf(x, targs)`, as in metafor (e.g.,
  `transf = metafor::transf.exp.int, targs = list(tau2 = 0)`).

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

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

  las

  :   axis-label style for base graphics (default: 1)

  back

  :   background color for the confidence band (default: "grey90"). Set
      to `NA` to suppress.

  arc.res

  :   integer specifying the number of line segments used to draw arcs
      (default: 100)

  main

  :   character string for plot title (default: no title)

  as_data

  :   if `TRUE`, returns plot data instead of creating the plot

## Value

`radial.brma` returns `NULL` invisibly if `plot_type = "base"` or a
ggplot object if `plot_type = "ggplot"`.

If `as_data = TRUE`, returns a list with the computed plot data
including: `points`, `refline`, `band`, `arc`, `arc_ticks`, `ci_arc`,
`ci_ticks`, `xlim`, and `zlim`. Returned arc coordinates use plain
trigonometry; both renderers recompute arcs for their aspect ratio, so
the data do not encode final device geometry.

## Details

The radial (Galbraith) plot transforms each study's effect size and
standard error into a point in precision-standardized space:

- x-axis: precision = \\1/\sqrt{v_i + \hat{\tau}\_i^2}\\

- z-axis: standardized effect = \\y_i/\sqrt{v_i + \hat{\tau}\_i^2}\\

For ordinary and specialized multilevel models, the marginal
heterogeneity is common and \\\hat{\tau}\_i = \hat{\tau}\\.
Random-formula models use each row's marginal random-effect standard
deviation. Known sampling covariance uses its diagonal for \\v_i\\;
plotted points can therefore be dependent.

For a single normal model, the sloped parallelogram is the usual
marginal Gaussian reference band around the pooled-effect line. For
GLMMs it is a continuity-corrected effect-size approximation; for
selection, PET/PEESE, and product-space models it remains a
bias-adjusted posterior-mean reference rather than an exact envelope of
the fitted observation law. The arc on the right side allows reading
individual effect sizes by projecting from the origin through a point to
the arc; when `center = TRUE`, the plotted slope is relative to the
pooled effect.

This function requires an intercept-only model without scale regression
or likelihood weights. Likelihood-weighted fits do not have a single
radial precision target because the observational variance and
power-likelihood information scales differ.

`galbraith()` is a same-argument alias for `radial()`.

## References

Galbraith, R. F. (1988). Graphical display of estimates having differing
standard errors. *Technometrics*, 30(3), 271-281.

Galbraith, R. F. (1994). Some applications of radial plots. *Journal of
the American Statistical Association*, 89(428), 1232-1242.

## See also

[`funnel.brma()`](https://fbartos.github.io/RoBMA/dev/reference/funnel.md),
[`pooled_effect.brma()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_effect.brma.md),
[`pooled_heterogeneity.brma()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_heterogeneity.brma.md)

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
  radial(fit)
  radial(fit, center = TRUE)
  radial(fit, plot_type = "ggplot")
  galbraith(fit)
}
} # }
```
