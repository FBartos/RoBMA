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
  ...
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
  ...
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
  present.

- parameter_mods:

  moderator term for location regression.

- parameter_scale:

  term for scale regression.

- type:

  diagnostic plot type. Convenience wrappers set a type-specific default
  but still forward this argument to `plot_diagnostic.brma()`.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- lags:

  number of lags for autocorrelation plots. Defaults to 30.

## Value

`plot_diagnostic` returns the object returned by
[`BayesTools::JAGS_diagnostics()`](https://fbartos.github.io/BayesTools/reference/JAGS_diagnostics.html),
invisibly for base graphics.

## See also

[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)
