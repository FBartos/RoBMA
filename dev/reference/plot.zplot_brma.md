# Plot Zplot Results

Plots a zplot visualization showing the histogram of observed
z-statistics overlaid with model-implied densities.

## Usage

``` r
# S3 method for class 'zplot_brma'
plot(
  x,
  plot_type = "base",
  probs = c(0.025, 0.975),
  max_samples = if (inherits(x, "brma.mv")) 1000 else 10000,
  plot_fit = TRUE,
  plot_extrapolation = FALSE,
  plot_ci = TRUE,
  plot_thresholds = TRUE,
  from = -6,
  to = 6,
  by.hist = 0.5,
  length.out.hist = NULL,
  by.lines = 0.05,
  length.out.lines = NULL,
  dots_hist = NULL,
  dots_fit = NULL,
  dots_extrapolation = NULL,
  dots_thresholds = NULL,
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

  maximum posterior samples for density estimation. Defaults to 1000 for
  `brma.mv` models and 10000 for other models. Use `Inf` to use all
  posterior samples. This plotting budget does not change the EDR and
  missing-study summary budget in
  [`as_zplot()`](https://fbartos.github.io/RoBMA/dev/reference/as_zplot.brma.md).

- plot_fit:

  whether to show fitted density (with bias adjustments). Defaults to
  `TRUE`.

- plot_extrapolation:

  whether to show extrapolated density (bias removed). Defaults to
  `FALSE`.

- plot_ci:

  whether to show credible interval bands. Defaults to `TRUE`.

- plot_thresholds:

  whether to show significance threshold lines. Defaults to `TRUE`.

- from, to:

  z-value range for plotting. Defaults to `-6` and `6`.

- by.hist:

  bin width for histogram. Defaults to 0.5.

- length.out.hist:

  number of histogram bins (alternative to `by.hist`).

- by.lines:

  step size for density lines. Defaults to 0.05.

- length.out.lines:

  number of density points (alternative to `by.lines`).

- dots_hist:

  graphical parameters for histogram (list).

- dots_fit:

  graphical parameters for fit lines (list).

- dots_extrapolation:

  graphical parameters for extrapolation lines (list).

- dots_thresholds:

  graphical parameters for threshold lines (list).

- integration_control:

  numerical integration settings created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/dev/reference/set_selection_likelihood_control.md).
  Defaults to the settings stored by
  [`as_zplot()`](https://fbartos.github.io/RoBMA/dev/reference/as_zplot.brma.md);
  an explicit value overrides them for this plot. Legacy objects without
  stored settings require an explicit control, or recreation with
  [`as_zplot()`](https://fbartos.github.io/RoBMA/dev/reference/as_zplot.brma.md),
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

  additional graphical parameters passed to components, such as `xlim`
  and `ylim` for the displayed axis ranges.

## Value

`NULL` invisibly for base graphics, or a ggplot2 object.

## Details

All model-implied curves use the `conditioning_depth` stored by
[`as_zplot()`](https://fbartos.github.io/RoBMA/dev/reference/as_zplot.brma.md).
Recreate the zplot object to change the predictive target.

The curves are posterior means of the density at each z value. Their
credible bands (`plot_ci = TRUE`) are pointwise Harrell-Davis quantiles
of the posterior density draws at the probabilities `probs`, computed by
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

The plot can display two density curves:

- Fit (black):

  Model-implied density including publication bias adjustments. This
  represents the expected distribution of z-statistics given the
  estimated selection process. Shown by default.

- Extrapolation (blue):

  The bias-adjusted reference without selective reporting. Multivariate
  models, retained estimate-level or sampling-error sources, applicable
  integrated contexts, and best-rule selection models use the normalized
  pre-selection Gaussian marginal. Released univariate product targets
  retain their inverse-weight scaling convention. Shown only when
  `plot_extrapolation = TRUE`, since the extrapolated curve rescales the
  display by an inferred selection convention rather than describing the
  observed literature.

## See also

[`hist.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/hist.zplot_brma.md),
[`lines.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/lines.zplot_brma.md),
[`summary.zplot_brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.zplot_brma.md)
