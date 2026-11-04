# Forest Plot for brma Objects Using metafor

[`metafor::forest()`](https://wviechtb.github.io/metafor/reference/forest.html)
creates a forest plot for a fitted brma object by adapting RoBMA results
to metafor's public forest plotting API. Study rows are drawn by
[`forest.default`](https://wviechtb.github.io/metafor/reference/forest.default.html)
and the RoBMA pooled posterior summary is added with
[`addpoly.default`](https://wviechtb.github.io/metafor/reference/addpoly.default.html).

## Usage

``` r
# S3 method for class 'brma'
forest(
  x,
  addfit = TRUE,
  addpred = FALSE,
  newdata = NULL,
  predstyle = "line",
  mlab = NULL,
  row = -1,
  level = 95,
  bias_adjusted = TRUE,
  conditional = FALSE,
  as_data = FALSE,
  predlim,
  border,
  constarea,
  ...
)

# S3 method for class 'metafor_forest.brma'
forest(
  x,
  addfit = TRUE,
  predstyle = "line",
  mlab = NULL,
  row = -1,
  predlim,
  border,
  constarea,
  ...
)
```

## Arguments

- x:

  a fitted brma object or an object returned by
  [`as_metafor_forest()`](https://fbartos.github.io/RoBMA/reference/as_metafor_forest.md).

- addfit:

  logical; whether to add the RoBMA pooled-effect summary polygon.
  Defaults to `TRUE`.

- addpred:

  logical; whether to add a posterior prediction interval for one new
  true effect. Defaults to `FALSE`.

- newdata:

  optional one-row prediction design used when `addpred = TRUE`. See
  [`as_metafor_forest()`](https://fbartos.github.io/RoBMA/reference/as_metafor_forest.md).

- predstyle:

  character; prediction interval style passed to
  [`addpoly.default`](https://wviechtb.github.io/metafor/reference/addpoly.default.html).
  Defaults to `"line"`.

- mlab:

  label for the summary row. Defaults to `"Pooled Effect"` for a pooled
  summary and `"Predicted Effect"` when `newdata` requests a
  moderator-specific prediction.

- row:

  numeric row position for the pooled-effect summary polygon. Defaults
  to `-1`.

- level:

  credible/confidence level in percent, greater than 0 and less
  than 100. A value below 1 is read as a proportion (`0.95` is 95
  percent), so `1` means 1 percent. Defaults to `95`.

- bias_adjusted:

  logical; whether pooled and prediction summaries should adjust for
  publication bias. Defaults to `TRUE`.

- conditional:

  logical; whether to return conditional posterior summaries for RoBMA
  product-space objects. Defaults to `FALSE`.

- as_data:

  logical; if `TRUE`, returns the prepared data instead of creating a
  plot.

- predlim:

  argument passed to
  [`addpoly.default`](https://wviechtb.github.io/metafor/reference/addpoly.default.html)
  for prediction interval rendering.

- border, constarea:

  arguments passed to
  [`addpoly.default`](https://wviechtb.github.io/metafor/reference/addpoly.default.html)
  for summary polygon rendering.

- ...:

  additional arguments passed to
  [`forest.default`](https://wviechtb.github.io/metafor/reference/forest.default.html).
  This includes auxiliary-column arguments such as `ilab`, `ilab.lab`,
  `ilab.xpos`, and `ilab.pos`. Shared graphical arguments such as
  `annotate`, `digits`, `transf`, `atransf`, `efac`, `col`, `lty`,
  `fonts`, and `cex` are also forwarded to
  [`addpoly.default`](https://wviechtb.github.io/metafor/reference/addpoly.default.html).

## Value

Invisibly returns the object returned by
[`forest.default`](https://wviechtb.github.io/metafor/reference/forest.default.html).
If `as_data = TRUE`, returns a `metafor_forest.brma` object.

## Details

This method deliberately uses metafor's public vector-based forest
plotting API instead of constructing a fake metafor model object. The
observed study rows are shown as usual Wald intervals based on the
supplied effect sizes and standard errors. The model summary row uses
the posterior mean and credible interval from
[`pooled_effect`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md).
When `addpred = TRUE` without `newdata`, the displayed center, credible
interval, and marginal prediction interval use the average-design target
returned by
[`pooled_effect`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md).
With explicit `newdata`, all three refer to that one prediction row and
the prediction interval uses
[`predict.brma`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
with `type = "estimate"` and `conditioning_depth = "marginal"`. For
metafor's `predstyle = "shade"` and `predstyle = "dist"`, the adapter
supplies the predictive density. Without `newdata` it is the posterior
predictive density, a normal mixture over the posterior draws, smoothed
by a normal kernel for display (a kernel density estimate of predictive
draws without their simulation noise); the smoothing keeps point masses
(draws without heterogeneity) visible. With `newdata` it is a kernel
density estimate of the predictive draws from
[`predict.brma`](https://fbartos.github.io/RoBMA/reference/predict.brma.md).
The displayed prediction-interval endpoints are the quantiles reported
by
[`pooled_effect`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)
or
[`predict.brma`](https://fbartos.github.io/RoBMA/reference/predict.brma.md).

Additional study-level columns can be added with metafor's `ilab`
interface. When supplied to
[`as_metafor_forest()`](https://fbartos.github.io/RoBMA/reference/as_metafor_forest.md)
or to `metafor::forest(..., as_data = TRUE)`, those arguments are stored
in the adapter object's `forest_args`; presentation options such as
`addfit`, `row`, and `predstyle` are retained by
`forest(..., as_data = TRUE)` for later rendering. Explicit plot-time
arguments override stored defaults. When supplied to
[`metafor::forest()`](https://wviechtb.github.io/metafor/reference/forest.html)
at plot time, they are forwarded directly.

Use
[`as_metafor_forest()`](https://fbartos.github.io/RoBMA/reference/as_metafor_forest.md)
to prepare the adapter object without requiring a metafor installation.

## See also

[`as_metafor_forest()`](https://fbartos.github.io/RoBMA/reference/as_metafor_forest.md),
[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md),
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
  metafor::forest(fit)
  metafor::forest(fit, addpred = TRUE, atransf = exp)

  count_columns <- as.matrix(dat[, c("tpos", "tneg", "cpos", "cneg")])
  metafor::forest(
    fit,
    addpred   = TRUE,
    slab      = paste(dat$author, dat$year),
    ilab      = count_columns,
    ilab.lab  = c("T+", "T-", "C+", "C-"),
    ilab.xpos = c(-8.5, -7.5, -6.5, -5.5),
    xlim      = c(-10, 4),
    atransf   = exp
  )
}
} # }
```
