# Plots brma Object

`plot.brma` visualizes posterior (and prior) distribution a brma object.

## Usage

``` r
# S3 method for class 'brma'
plot(
  x,
  parameter,
  parameter_mods,
  parameter_scale,
  prior = FALSE,
  standardized_coefficients = FALSE,
  conditional = FALSE,
  output_measure = NULL,
  transform = NULL,
  plot_type = "base",
  dots_prior = NULL,
  density_method = c("KDE", "qCMDE", "IWMDE"),
  density_control = NULL,
  component = "auto",
  ...
)

# S3 method for class 'brma'
lines(
  x,
  parameter = NULL,
  parameter_mods = NULL,
  parameter_scale = NULL,
  prior = FALSE,
  standardized_coefficients = FALSE,
  conditional = FALSE,
  output_measure = NULL,
  transform = NULL,
  plot_type = "base",
  dots_prior = NULL,
  density_method = c("KDE", "qCMDE", "IWMDE"),
  density_control = NULL,
  component = "auto",
  ...
)
```

## Arguments

- x:

  a fitted `brma`, `BMA`, or `RoBMA` object.

- parameter:

  a parameter to be plotted. Defaults to `"mu"` for the effect size, or
  to the meta-regression intercept when moderators are present.
  Additional options are `"tau"`, `"rho"` for multilevel models,
  `"PET"`, `"PEESE"`, and `"omega"` or `"weightfunction"` for selection
  models. Use
  [`plot_pet_peese()`](https://fbartos.github.io/RoBMA/reference/plot_pet_peese.md)
  for PET/PEESE regression plots. Factor terms select all coefficient
  cells; use a semantic selector such as `"group[level]"` to plot one
  cell. Structural point masses are retained. A quantity of the fitted
  model that is no model parameter (an inclusion indicator such as
  `"mu_x_indicator"`, `"bias_indicator"`, a weight-function coordinate
  such as `"omega[1]"`, or a latent cluster effect) stops with the
  classes `BayesTools_hypothesis_target` and
  `BayesTools_hypothesis_unavailable`, with which
  [`hypothesis()`](https://fbartos.github.io/RoBMA/reference/hypothesis.md)
  refuses it. Ordered level cells retain their semantic marginal draws
  and prior metadata on the requested predictor scale, including
  conditional row alignment.

- parameter_mods:

  legacy moderator selector. Prefer `parameter` with
  `component = "mods"`. Use `"intercept"` for the adjusted effect in
  meta-regression models.

- parameter_scale:

  legacy scale-regression selector. Prefer `parameter` with
  `component = "scale"`. Use `"intercept"` for the heterogeneity
  intercept in location-scale models.

- prior:

  whether prior distribution should be added to figure. Defaults to
  `FALSE`.

- standardized_coefficients:

  whether to plot moderator and scale-regression coefficients on the
  standardized predictor scale. Defaults to `FALSE`.

- conditional:

  whether to plot the conditional posterior distribution for RoBMA
  product-space objects. Defaults to `FALSE`.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional plotting transformation. `"EXP"` exponentiates effect-size
  location and individual meta-regression coefficients for fitted
  log-scale OR, RR, HR, and IRR models. For scale-regression
  coefficients, `"EXP"` displays the multiplicative change in
  heterogeneity. `"LOG"` displays a positive heterogeneity intercept or
  random-effect variance or SD multiplier on the log scale. The
  transformation is applied to KDE and precomputed qCMDE/IWMDE densities
  with the corresponding Jacobian. For automatically unscaled
  scale-formula intercepts, use `standardized_coefficients = TRUE` with
  qCMDE/IWMDE.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- dots_prior:

  list of additional graphical arguments to be passed to the plotting
  function of the prior distribution. Supported arguments are `lwd`,
  `lty`, `col`, and `col.fill`, to adjust the line thickness, line type,
  line color, and fill color of the prior distribution respectively. For
  ordered predictor targets, omitted line colors and types match each
  posterior level's color with a dashed prior curve, including selected
  cells. For selected ordered cells with `prior = TRUE`,
  `color`/`colour` aliases resolve to `col`, with priority `col`,
  `color`, then `colour`; `linetype` takes precedence over `lty` in both
  renderers. Explicit `col`, `lty`, or `linetype` values take
  precedence; the ordinary prior defaults are retained for other
  targets. Ordered metadata validation is restricted to the selected
  terms, including selected group members; unrelated ordered terms do
  not affect prior styles.

- density_method:

  posterior density method. `"KDE"` uses the standard BayesTools kernel
  density estimate. `"qCMDE"` attaches RoBMA row-normalized q-grid
  conditional densities. `"IWMDE"` attaches Chen-style moment-matched
  IWMDE densities. qCMDE is preferred when its additional normalization
  cost is acceptable; IWMDE can be faster but is more sensitive to its
  fitted conditional weights. Matching is case-insensitive. Ordered
  effects use the same normalized conditional-source estimator for both
  methods, with their respective posterior-row budgets. Selected
  weightfunction models near zero cumulative ordered effects can take
  more than 10 minutes with qCMDE or IWMDE. Exact regular ordered
  offsets that need a finite supplier integral use a positive
  conditional block; offsets classified as nonfinite or unsupported are
  omitted only at that exact value from a density curve. Point
  hypotheses retain their structural or precision refusal. For semantic
  random-effect quantities, qCMDE/IWMDE support direct scalar fitted
  sources, allocated component SDs backed by a scalar aggregate
  (`tau_total` or `tau_common`), and allocation proportions or
  multipliers backed by a fitted simplex coordinate. Other nonlinear
  derived quantities remain KDE-only. Shared inclusion gates are
  supported for aggregate and component SDs and allocation proportions.
  SDs retain their excluded zero branches; proportions condition on
  positive total heterogeneity. Allocated component SDs in a
  predictor-scaled block require `standardized_coefficients = TRUE` or
  `density_method = 'KDE'` for an original-scale request. This
  block-wide limitation also covers a centered intercept whose value
  happens to be unchanged by scaling. A positive point source multiplied
  by stochastic allocation weights has a continuous target law, but its
  source-held-weight chart cannot sweep that target; use
  `density_method = 'KDE'`. Source-zero and inclusion-gate atoms retain
  their declared locations. In independently gated allocations, realized
  `tau_total`/`tau2_total` and `tau2_prop(...)` are KDE-only because
  they combine multiple gate and allocation coordinates. Their
  structural point masses are displayed separately from the continuous
  density; `tau2_prop(...)` is conditioned on positive realized total
  heterogeneity. qCMDE/IWMDE are not available for non-known-`V`
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  random-formula models,
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  models with component-specific scale formulas (a named `scale` list),
  or selection-weightfunction coordinates requiring joint replacement.
  IWMDE is also unavailable for binomial and Poisson GLMMs; use qCMDE
  for GLMM density plots. qCMDE/IWMDE densities are evaluated on the
  fitted coefficient coordinate. If automatic predictor scaling changes
  the requested display coordinate, use
  `standardized_coefficients = TRUE`; RoBMA does not infer a coefficient
  transformation from posterior draws.

- density_control:

  named list of density-estimation settings. Supported entries are
  `n_points` (default `100`), `samples` (default `500` for qCMDE and
  `1000` for IWMDE density curves), `target_relative_mcse` (default
  `0.05`), `display_grid` (default `"adaptive"`), `normalization_points`
  (default `NULL`, resolved to `max(50, n_points)`), and
  `normalization_prob` (default `0.999`; see
  [`hypothesis()`](https://fbartos.github.io/RoBMA/reference/hypothesis.md)
  for the per-row truncation it permits). `integration_control` (default
  `NULL`) retains the fitted selection-integration settings. Supply a
  control created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md)
  to change those settings for this post-fit calculation, for example
  `list(integration_control = set_selection_likelihood_control(max_points_per_scramble = 32768))`.
  This entry is available for Gaussian selection models with a fitted
  integration plan. The fitted object and posterior draws are unchanged.
  The maximum point budget controls factor QMC fallback; analytic and
  deterministic quadrature rules remain unchanged. `samples` controls
  the fixed posterior-row budget for the density curve.
  `target_relative_mcse` is a point-ordinate diagnostic target and does
  not alter this fixed-budget density plot. The normalization entries
  are used with `density_method = "qCMDE"` and
  `density_method = "IWMDE"`. Curve diagnostics apply local
  relative-MCSE, effective-sample-size, and contribution-concentration
  gates over the empirical 5–95 percent bulk and record the 5 and 95
  percent tail checkpoints. The entire display retains an absolute MCSE
  safeguard relative to the density peak. Requested point ordinates used
  for Bayes factors retain separate strict local diagnostics. Increase
  the row and normalization budgets and compare results when density
  diagnostics report low effective sample size, concentrated
  contributions, or unstable normalization.

- component:

  parameter component. Defaults to `"auto"`, which infers the component
  when possible. Use `"mods"` (alias `"location"`), `"scale"`,
  `"random"`, or `"bias"` to disambiguate terms used in multiple model
  components. The random component selects semantic standard deviation,
  correlation, and allocation parameters from
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  random formulas. A `parameter_mods` or `parameter_scale` selection
  with another `component` stops with an error of class
  `RoBMA_component_mismatch`.

- ...:

  list of additional graphical arguments to be passed to the plotting
  function. Supported arguments are `lwd`, `lty`, `col`, `col.fill`,
  `xlab`, `ylab`, `ylab2`, `main`, `xlim`, `ylim`, and `ylim2` to adjust
  the line thickness, line type, line color, fill color, x-label,
  density-axis label, probability-axis label, title, x-axis range,
  density-axis range, and probability-axis range respectively. Set
  `ylim2` on the initial
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) call;
  subsequent base [`lines()`](https://rdrr.io/r/graphics/lines.html)
  calls reuse that probability mapping.

## Value

`plot.brma` returns either `NULL` if `plot_type = "base"` or a `ggplot2`
object (a named list of plots for a whole factor term) if
`plot_type = "ggplot"`. `lines.brma` returns `NULL` for
`plot_type = "base"` and ggplot2 layer(s) for `plot_type = "ggplot"`.

## Details

A qCMDE/IWMDE request that the fitted model does not support stops
before any density is estimated with an error of class
`RoBMA_density_method_unavailable` and one class naming the cause:
`RoBMA_density_method_random_unknown_v`
([`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
random-formula models without known `V`),
`RoBMA_density_method_scale_components` (component-specific scale
formulas), or `RoBMA_density_method_glmm` (IWMDE for binomial and
Poisson GLMMs). Requests that are unavailable for the plotted quantity
stop with the same parent class and
`RoBMA_density_method_conditional_random` (conditional random-effect
plots, which are KDE-only), `RoBMA_density_method_random_target` (a
random-effect quantity without a supported scalar random-component
coordinate), `RoBMA_density_method_bias_target` (a composite
weightfunction target; use `density_method = 'KDE'`), or
`RoBMA_density_method_original_scale` (a coefficient or factor cell
whose original-scale value is not a linear combination of its fitted
coefficients, such as the exponentiated intercept of a scale formula
with standardized predictors; use `standardized_coefficients = TRUE`). A
density whose construction fails stops with an error of class
`RoBMA_density_construction_error` (a `RoBMA_density_error`); use
`density_method = "KDE"` to estimate it from the posterior draws. A
fitted chain with no eligible row of the target adds the class
`RoBMA_density_chain_coverage`, and the message names longer chains and
their convergence as the other remedy.
[`hypothesis()`](https://fbartos.github.io/RoBMA/reference/hypothesis.md)
refuses qCMDE/IWMDE for the same causes with these classes after its
own.
[`plot_prior()`](https://fbartos.github.io/RoBMA/reference/plot_prior.md)
stops with `RoBMA_density_method_original_scale` for the prior of a term
with several fitted coordinates whose original-scale coordinates the
standardization changes (see there).

A returned qCMDE/IWMDE estimate that fails the density availability
checks raises a `RoBMA_density_plot_error` identifying the plotted
parameter. Its `density_diagnostics` field contains the returned
diagnostic records, including computed density values and selected rows,
without executable estimate plans that could retain the fitted object.

`lines.brma()` adds the posterior density to an existing base plot. Omit
`add`: this method always adds, and supplying `add` triggers an error.
With `plot_type = "ggplot"`, it returns ggplot2 layer(s) that can be
added to a `plot.brma(..., plot_type = "ggplot")` object with `+`. For
base plots containing point masses, the initial `plot.brma()` call
establishes the secondary probability axis and `lines.brma()` reuses it
across fitted objects. An overlaid point mass outside the initial
`ylim2` is clipped with a warning.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  plot(fit, parameter = "mu")
  lines(fit, parameter = "mu", col = "blue", lwd = 2)

  ggplot_fit <- plot(fit, parameter = "mu", plot_type = "ggplot")
  ggplot_fit + lines(fit, parameter = "mu", plot_type = "ggplot", col = "blue")

  plot(fit, parameter = "tau", prior = TRUE)
  plot(fit, parameter = "PET")
}
} # }

```
