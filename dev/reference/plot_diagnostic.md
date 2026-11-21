# Plot MCMC Diagnostics

`plot_diagnostic` creates visual MCMC diagnostics for a fitted brma
object. Convenience wrappers are available for trace, density, and
autocorrelation plots.

## Usage

``` r
plot_diagnostic(x, ...)

# S3 method for class 'brma'
plot_diagnostic(
  x,
  parameter = NULL,
  parameter_mods = NULL,
  parameter_scale = NULL,
  type,
  plot_type = "base",
  lags = 30,
  component = "auto",
  ...,
  statistic = "value",
  batch_size = NULL
)

plot_diagnostic_autocorrelation(x, ...)

# S3 method for class 'brma'
plot_diagnostic_autocorrelation(
  x,
  parameter = NULL,
  parameter_mods = NULL,
  parameter_scale = NULL,
  type = "autocorrelation",
  plot_type = "base",
  lags = 30,
  component = "auto",
  ...
)

plot_diagnostic_trace(x, ...)

# S3 method for class 'brma'
plot_diagnostic_trace(
  x,
  parameter = NULL,
  parameter_mods = NULL,
  parameter_scale = NULL,
  type = "trace",
  plot_type = "base",
  lags = 30,
  component = "auto",
  ...,
  statistic = "value",
  batch_size = NULL
)

plot_diagnostic_density(x, ...)

# S3 method for class 'brma'
plot_diagnostic_density(
  x,
  parameter = NULL,
  parameter_mods = NULL,
  parameter_scale = NULL,
  type = "density",
  plot_type = "base",
  lags = 30,
  component = "auto",
  ...
)
```

## Arguments

- x:

  a fitted brma object

- ...:

  additional graphical arguments passed through RoBMA's diagnostic setup
  to
  [`BayesTools::JAGS_diagnostics()`](https://fbartos.github.io/BayesTools/reference/JAGS_diagnostics.html).

- parameter:

  base parameter to plot. Defaults to `NULL`, which uses `"mu"` or the
  meta-regression intercept. Valid values include `"mu"`, `"tau"`,
  `"rho"`, `"PET"`, `"PEESE"`, and `"omega"` or `"weightfunction"` when
  present. A quantity of the fitted model that is no model parameter
  (e.g. an inclusion indicator or a weight-function coordinate such as
  `"omega[1]"`) stops with the classes `BayesTools_hypothesis_target`
  and `BayesTools_hypothesis_unavailable`, as in
  [`plot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.brma.md).

- parameter_mods:

  legacy moderator selector. Prefer `parameter` with
  `component = "mods"`.

- parameter_scale:

  legacy scale-regression selector. Prefer `parameter` with
  `component = "scale"`.

- type:

  diagnostic plot type. Convenience wrappers set a type-specific default
  but still forward this argument to `plot_diagnostic.brma()`.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- lags:

  number of lags for autocorrelation plots. Defaults to 30.

- component:

  parameter component. Defaults to `"auto"`, which infers the component
  when possible. Use `"mods"` (alias `"location"`), `"scale"`,
  `"random"`, or `"bias"` to disambiguate terms used in multiple model
  components. The random component selects semantic standard deviation,
  correlation, and allocation parameters from
  [`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
  random formulas. A `parameter_mods` or `parameter_scale` selection
  with another `component` stops with an error of class
  `RoBMA_component_mismatch`.

- statistic:

  trace statistic. `"value"` shows the sampled parameter. `"inclusion"`
  shows the fraction of each batch spent in the alternative components
  and needs every component declared null or alternative, with both
  present. `"component"` shows the fraction spent in each mixture
  component in separate plots (one plot for a binary null/alternative
  choice). `"switching"` shows the fraction of adjacent draws within
  each batch that changed component. Use `parameter = "bias"` for the
  publication-bias mixture or select a whole moderator term. A
  random-effect quantity must have a unique indicator from its prior or
  allocation gate. Otherwise the diagnostic raises
  `RoBMA_indicator_unavailable` (also a `BayesTools_plot_condition`).

- batch_size:

  optional positive integer giving the number of retained draws per
  batch for indicator traces. Switching needs at least two retained
  draws per batch and chain; a final one-draw batch is not drawn. By
  default the shortest chain supplies approximately its square root in
  batches. Each chain is batched separately and plotted at the batches'
  midpoint iterations.

## Value

`plot_diagnostic` returns the object returned by
[`BayesTools::JAGS_diagnostics()`](https://fbartos.github.io/BayesTools/reference/JAGS_diagnostics.html),
including a list when several coordinate plots are produced, invisibly
for base graphics.

## See also

[`summary.brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  plot_diagnostic_trace(fit, parameter = "mu", statistic = "inclusion")
  plot_diagnostic_trace(fit, parameter = "bias", statistic = "component")
  plot_diagnostic(fit, parameter = "bias", type = "trace", statistic = "switching")
}
} # }
```
