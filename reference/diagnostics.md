# Checks a fitted RoBMA object

`diagnostics` creates visual checks of individual models convergence.
Numerical overview of individual models can be obtained by
`summary(object, type = "models", diagnostics = TRUE)`, or even more
detailed information by `summary(object, type = "individual")`.

## Usage

``` r
diagnostics(
  fit,
  parameter,
  type,
  plot_type = "base",
  show_models = NULL,
  lags = 30,
  title = is.null(show_models) | length(show_models) > 1,
  ...
)

diagnostics_autocorrelation(
  fit,
  parameter = NULL,
  plot_type = "base",
  show_models = NULL,
  lags = 30,
  title = is.null(show_models) | length(show_models) > 1,
  ...
)

diagnostics_trace(
  fit,
  parameter = NULL,
  plot_type = "base",
  show_models = NULL,
  title = is.null(show_models) | length(show_models) > 1,
  ...
)

diagnostics_density(
  fit,
  parameter = NULL,
  plot_type = "base",
  show_models = NULL,
  title = is.null(show_models) | length(show_models) > 1,
  ...
)
```

## Arguments

- fit:

  a fitted RoBMA object

- parameter:

  a parameter to be plotted. Either `"mu"`, `"tau"`, `"omega"`, `"PET"`,
  or `"PEESE"`.

- type:

  type of MCMC diagnostic to be plotted. Options are `"chains"` for the
  chains' trace plots, `"autocorrelation"` for autocorrelation of the
  chains, and `"densities"` for the overlaying densities of the
  individual chains. Can be abbreviated to first letters.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- show_models:

  MCMC diagnostics of which models should be plotted. Defaults to `NULL`
  which plots MCMC diagnostics for a specified parameter for every model
  that is part of the ensemble.

- lags:

  number of lags to be shown for `type = "autocorrelation"`. Defaults to
  `30`.

- title:

  whether the model number should be displayed in title. Defaults to
  `TRUE` when more than one model is selected.

- ...:

  additional arguments to be passed to
  [par](https://rdrr.io/r/graphics/par.html) if `plot_type = "base"`.

## Value

`diagnostics` returns either `NULL` if `plot_type = "base"` or an
object/list of objects (depending on the number of parameters to be
plotted) of class 'ggplot2' if `plot_type = "ggplot2"`.

## Details

The visualization functions are based on
[stan_plot](https://mc-stan.org/rstan/reference/stan_plot.html) function
and its color schemes.

## See also

[`RoBMA()`](reference/RoBMA.md),
[`summary.RoBMA()`](reference/summary.RoBMA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# using the example data from Anderson et al. 2010 and fitting the default model
# (note that the model can take a while to fit)
fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)

### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
# diagnostics function allows to visualize diagnostics of a fitted RoBMA object, for example,
# the trace plot for the mean parameter in each model model
diagnostics(fit, parameter = "mu", type = "chain")

# in order to show the trace plot only for the 11th model, add show_models parameter
diagnostics(fit, parameter = "mu", type = "chain", show_models = 11)

# furthermore, the autocorrelations
diagnostics(fit, parameter = "mu", type = "autocorrelation")

# and overlying densities for each plot can also be visualize
diagnostics(fit, parameter = "mu", type = "densities")
} # }

```
