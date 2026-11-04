# Plot Estimated Marginal Means

Plots estimated marginal means stored in a `marginal_means.brma` object
using
[`BayesTools::plot_marginal()`](https://fbartos.github.io/BayesTools/reference/plot_marginal.html).

## Usage

``` r
# S3 method for class 'marginal_means.brma'
plot(
  x,
  parameter,
  type = NULL,
  prior = FALSE,
  plot_type = "base",
  dots_prior = NULL,
  output_measure = NULL,
  transform = NULL,
  density_method = NULL,
  density_control = NULL,
  add = FALSE,
  ...
)

# S3 method for class 'marginal_means.brma'
lines(x, parameter, prior = FALSE, ...)
```

## Arguments

- x:

  a `marginal_means.brma` object.

- parameter:

  moderator term to plot. Use the original term name, for example
  `"measure"`, `"intercept"` for the intercept when available, `"mu"` as
  an intercept alias, or the internal parameter name, for example
  `"mu_measure"`.

- type:

  for RoBMA product-space objects, whether to plot model-averaged
  (`"averaged"`) or conditional (`"conditional"`) marginal means.
  Defaults to `"averaged"` and is available only for RoBMA marginal
  means.

- prior:

  whether the marginal prior distribution should be added to the plot.
  Defaults to `FALSE`.

- plot_type:

  whether to use base R graphics (`"base"`) or ggplot2 (`"ggplot"`).
  Defaults to `"base"`.

- dots_prior:

  list of additional graphical arguments passed to the prior plotting
  function.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- density_method:

  posterior density method. The default `NULL` reuses the method stored
  on `x`; specify `"KDE"` to override a stored qCMDE/IWMDE method.
  `"qCMDE"` and `"IWMDE"` compute any missing densities for the plotted
  marginal means before plotting and never silently mix with KDE.
  Plot-time qCMDE/IWMDE computation can be slow and is used only for the
  current plot call; precompute with
  `marginal_means(..., density_method = "qCMDE")` or
  `marginal_means(..., density_method = "IWMDE")` to reuse densities
  across plots. Prefer qCMDE for likelihood-aware final results; IWMDE
  can be faster but is more sensitive to its fitted conditional weights.
  A qCMDE/IWMDE calculation for selected weightfunction models near zero
  cumulative ordered effects can take more than 10 minutes. A
  qCMDE/IWMDE request that the source model does not support stops
  before any density is estimated with an error of class
  `RoBMA_density_method_unavailable` and the class naming its cause (see
  [`marginal_means.brma()`](https://fbartos.github.io/RoBMA/reference/marginal_means.brma.md)).

- density_control:

  named list of qCMDE/IWMDE density-estimation settings. Supported
  entries are `n_points` (default `100`), `samples` (default `500` for
  qCMDE and `1000` for IWMDE density curves), `target_relative_mcse`
  (default `0.05`), `display_grid` (default `"adaptive"`),
  `normalization_points` (default `NULL`, resolved to
  `max(50, n_points)`), and `normalization_prob` (default `0.999`; see
  [`hypothesis()`](https://fbartos.github.io/RoBMA/reference/hypothesis.md)
  for the per-row truncation it permits). `integration_control` (default
  `NULL`) accepts
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md)
  settings for Gaussian selection models; see
  [`marginal_means.brma()`](https://fbartos.github.io/RoBMA/reference/marginal_means.brma.md)
  for an example. The fitted object and posterior draws are unchanged.
  `samples` controls the fixed posterior-row budget for the density
  curve. `target_relative_mcse` is a point-ordinate diagnostic target
  and does not alter this fixed-budget plot. Curve diagnostics use the
  empirical 5–95 percent bulk, report the 5 and 95 percent tail
  checkpoints, and retain a whole-curve absolute MCSE safeguard relative
  to the density peak. Supplying `density_control` forces recomputation
  for the plotted marginal means instead of reusing stored densities.
  Increase row and normalization budgets when diagnostics report
  concentrated contributions or unstable normalization.

- add:

  whether to add the densities to an existing plot. This is used by
  `lines.marginal_means.brma()` and defaults to `FALSE`.

- ...:

  additional graphical arguments passed to
  [`BayesTools::plot_marginal()`](https://fbartos.github.io/BayesTools/reference/plot_marginal.html).

## Value

`NULL` invisibly if `plot_type = "base"` or a ggplot object if
`plot_type = "ggplot"`.

## Details

`lines.marginal_means.brma()` adds posterior densities to an existing
base plot. With `plot_type = "ggplot"`, it returns ggplot2 layer(s) that
can be added to a marginal-means plot with `+`.
