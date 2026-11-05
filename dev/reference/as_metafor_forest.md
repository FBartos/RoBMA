# Prepare brma Data for metafor Forest Plots

`as_metafor_forest` extracts observed effect sizes, standard errors,
study labels, and RoBMA posterior summary intervals in a structure that
can be passed to metafor's forest plotting functions.

## Usage

``` r
as_metafor_forest(x, ...)

# S3 method for class 'brma'
as_metafor_forest(
  x,
  level = 95,
  addpred = FALSE,
  newdata = NULL,
  bias_adjusted = TRUE,
  conditional = FALSE,
  ...
)
```

## Arguments

- x:

  a fitted brma object.

- ...:

  additional named arguments stored in `forest_args` and later passed to
  [`forest.default`](https://wviechtb.github.io/metafor/reference/forest.default.html).
  Shared graphical arguments are also stored in `addpoly_args`. This
  includes metafor's auxiliary-column arguments `ilab`, `ilab.lab`,
  `ilab.xpos`, and `ilab.pos`.

- level:

  credible/confidence level in percent, greater than 0 and less
  than 100. A value below 1 is read as a proportion (`0.95` is 95
  percent), so `1` means 1 percent. Defaults to `95`.

- addpred:

  logical; whether to compute a marginal posterior prediction interval
  for one new true effect. Defaults to `FALSE`.

- newdata:

  optional one-row prediction design used when `addpred = TRUE`. When
  omitted, the interval is evaluated at the average expanded design
  using
  [`pooled_effect()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_effect.md).
  When supplied, each row passed to
  [`predict.brma()`](https://fbartos.github.io/RoBMA/dev/reference/predict.brma.md)
  is one new true effect; with `bias_adjusted = FALSE`, it targets a
  selected new publication. The pooled reference declares no
  new-publication SE or selection event. The forest adapter accepts
  exactly one because it draws one summary polygon.

- bias_adjusted:

  logical; whether pooled and prediction summaries should adjust for
  publication bias. Defaults to `TRUE`.

- conditional:

  logical; whether to return conditional posterior summaries for RoBMA
  product-space objects. Defaults to `FALSE`.

## Value

A list of class `metafor_forest.brma` with study-level data, pooled
posterior summaries, and ready-to-use argument lists for
[`forest.default`](https://wviechtb.github.io/metafor/reference/forest.default.html)
and
[`addpoly.default`](https://wviechtb.github.io/metafor/reference/addpoly.default.html).
The `forest_args` include enough vertical space for the default summary
row. The returned object can also be plotted directly with
[`forest`](https://wviechtb.github.io/metafor/reference/forest.html).

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
  count_columns <- as.matrix(dat[, c("tpos", "tneg", "cpos", "cneg")])
  forest_data <- as_metafor_forest(
    fit,
    addpred   = TRUE,
    slab      = paste(dat$author, dat$year),
    ilab      = count_columns,
    ilab.lab  = c("T+", "T-", "C+", "C-"),
    ilab.xpos = c(-7, -6, -5, -4),
    xlim      = c(-10, 4),
    atransf   = exp
  )
  metafor::forest(forest_data, predstyle = "bar")

  # Or draw the metafor calls directly:
  do.call(metafor::forest.default, forest_data$forest_args)
  do.call(metafor::addpoly.default, forest_data$addpoly_args)
}
} # }
```
