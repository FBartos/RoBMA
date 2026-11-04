# Add Zplot Density Lines

Adds model-implied density lines to an existing zplot.

## Usage

``` r
# S3 method for class 'zplot_brma'
lines(
  x,
  plot_type = "base",
  probs = c(0.025, 0.975),
  max_samples = if (inherits(x, "brma.mv")) 1000 else 10000,
  plot_ci = TRUE,
  extrapolate = FALSE,
  from = -6,
  to = 6,
  by = 0.05,
  length.out = NULL,
  col = "black",
  as_data = FALSE,
  integration_control = x[["zplot"]][["data"]][["integration_control"]],
  parallel = FALSE,
  cores = min(4, RoBMA.get_option("max_cores")),
  ...
)
```

## Arguments

- x:

  a zplot_brma object.

- plot_type:

  graphics system: `"base"` or `"ggplot"`. Defaults to `"base"`.

- probs:

  quantiles for credible intervals. Defaults to `c(.025, .975)`. The
  bands are pointwise Harrell-Davis quantiles of the posterior density
  draws (see Details).

- max_samples:

  maximum posterior samples for density. Defaults to 1000 for `brma.mv`
  models and 10000 for other models. Use `Inf` to use all posterior
  samples. This plotting budget does not change the EDR and
  missing-study summaries stored by
  [`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md).

- plot_ci:

  whether to show credible interval bands. Defaults to `TRUE`.

- extrapolate:

  whether to remove bias adjustments. Defaults to `FALSE`.

- from, to:

  z-value range for density. Defaults to `-6` and `6`.

- by:

  step size for density points. Defaults to 0.05.

- length.out:

  number of density points (alternative to `by`).

- col:

  line color. Defaults to `"black"`.

- as_data:

  whether to return data instead of plotting. Defaults to `FALSE`.

- integration_control:

  numerical integration settings created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md).
  Defaults to the settings stored by
  [`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md);
  an explicit value overrides them for these lines. Legacy objects
  without stored settings require an explicit control, or recreation
  with
  [`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md),
  when their selected density needs integration.

- parallel:

  whether to distribute multivariate selection-model density
  calculations across local workers. Defaults to `FALSE`. Only posterior
  draws are divided; the full publication event and integration settings
  are retained. Scalar models, Gaussian-reference-only curves, and EDR
  summaries remain serial.

- cores:

  number of local workers when `parallel = TRUE`. Defaults to the
  smaller of 4 and `RoBMA.get_option("max_cores")`. One core or fewer
  posterior draws than workers uses serial computation. Parallel workers
  require matching installed package/native builds and share the
  configured selection-cache budget. Worker count is reported when
  parallel work starts.

- ...:

  additional graphical parameters.

## Value

`NULL` invisibly for base graphics, ggplot2 layers for ggplot, or a data
frame with columns `x`, `y`, `y_lCI`, `y_uCI` if `as_data = TRUE`.

## Details

The density uses the `conditioning_depth` stored by
[`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md).
Recreate the zplot object to change the predictive target.

The line is the posterior mean of the density at each z value, and its
credible band (`plot_ci = TRUE`, or `y_lCI` and `y_uCI` of the returned
data) is formed by pointwise Harrell-Davis quantiles of the posterior
density draws at the probabilities `probs`, computed by
[`BayesTools::harrell_davis_quantile()`](https://fbartos.github.io/BayesTools/reference/harrell_davis_quantile.html).
They estimate the same quantiles as the empirical quantiles of the draws
but smooth the kinks that the empirical quantiles show along z. The
bands remain pointwise, and, because every draw contributes, they can
differ visibly from the empirical quantiles when few draws are used
(`max_samples` below a few hundred). Near the edge of a point mass in
the draws (e.g., draws of models that fix the plotted quantity in a
model-averaged posterior) the band blends the point mass and the
remaining draws, and extreme probabilities (e.g.,
`probs = c(.0005, .9995)`) are not estimated more accurately than by the
empirical quantiles; see
[`BayesTools::harrell_davis_quantile()`](https://fbartos.github.io/BayesTools/reference/harrell_davis_quantile.html)
for these caveats.

When `extrapolate = FALSE`, the density includes all bias adjustments
(PET/PEESE regression, selection weights) representing the fitted model.
When `extrapolate = TRUE`, bias adjustments are removed to show the
hypothetical distribution without publication bias. Multivariate models,
retained estimate-level or sampling-error sources, applicable integrated
contexts, and best-rule selection models use a normalized pre-selection
Gaussian reference. Released univariate product targets retain
inverse-weight scaling and their extrapolated curve need not integrate
to one.

## See also

[`plot.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md),
[`hist.zplot_brma()`](https://fbartos.github.io/RoBMA/reference/hist.zplot_brma.md)
