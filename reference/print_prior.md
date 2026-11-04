# Print Prior Distributions

`print_prior` prints prior distributions stored in `brma`, `BMA`, and
`RoBMA` objects. This is a user-facing helper for inspecting priors
without extracting them from the internal `$priors` list.

## Usage

``` r
print_prior(x, ...)

# S3 method for class 'prior'
print_prior(x, ...)

# S3 method for class 'brma'
print_prior(
  x,
  parameter,
  parameter_mods,
  parameter_scale,
  component = "auto",
  ...
)
```

## Arguments

- x:

  a `brma`, `BMA`, or `RoBMA` object, or a prior distribution object.

- ...:

  additional arguments passed to the prior printing method. Use
  `silent = TRUE` for programmatic inspection without console output.

- parameter:

  character. Base parameter to print. If omitted, all stored prior
  distributions are printed. Common options are `"mu"`, `"tau"`,
  `"rho"`, `"PET"`, `"PEESE"`, `"omega"`, and `"bias"`, with aliases
  `"effect"` = `"mu"`, `"heterogeneity"` = `"tau"`, and
  `"weightfunction"` = `"omega"`. GLMM outcome priors `"pi"` and `"phi"`
  are available when present. Moderator and scale terms can also be
  selected by name when unambiguous. For
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  objects with random-effect formulas, `"random"` prints the structured
  random-effect prior specification, and the random-effect quantities of
  fitted objects (such as `"tau"`, `"tau_total"`, `"tau2"`,
  `"tau2_prop(study)"`, or `"study: tau"`) print their declared prior,
  or, for a quantity that follows from several declared priors, one line
  giving the derivation (for example `tau2 = tau^2` or
  `study: tau = sigma_total * sqrt(w[1])`, in the notation of the
  `"random"` specification) followed by the declared priors it refers
  to. A quantity without a prior density stops as in
  [`plot_prior()`](https://fbartos.github.io/RoBMA/reference/plot_prior.md).
  A character vector requests multiple base parameters. A quantity of
  the fitted model that is no model parameter (e.g. an inclusion
  indicator or a weight-function coordinate such as `"omega[1]"`) stops
  with the classes `BayesTools_hypothesis_target` and
  `BayesTools_hypothesis_unavailable`, as in
  [`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md).

- parameter_mods:

  legacy moderator selector. Prefer `parameter` with
  `component = "mods"`. Use `"intercept"` for the adjusted effect in
  meta-regression models.

- parameter_scale:

  legacy scale-regression selector. Prefer `parameter` with
  `component = "scale"`. Use `"intercept"` for the heterogeneity
  intercept in location-scale models.

- component:

  parameter component. Defaults to `"auto"`, which infers the component
  when possible. Use `"mods"` (alias `"location"`), `"scale"`, `"bias"`,
  or `"random"` (the random-effect quantities of fitted objects) to
  disambiguate terms used in multiple prior components. A
  `parameter_mods` or `parameter_scale` selection with another
  `component` stops with an error of class `RoBMA_component_mismatch`.

## Value

`print_prior` invisibly returns the selected prior distribution (for a
random-effect quantity that is not one declared prior, its prior
density, a BayesTools `prior_linear_density`). If multiple parameters
are requested, a named list of prior distributions is returned
invisibly.

## See also

[`plot_prior()`](https://fbartos.github.io/RoBMA/reference/plot_prior.md)
[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md)
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md)
[`prior()`](https://fbartos.github.io/RoBMA/reference/prior.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  priors <- BMA(
    yi          = yi,
    vi          = vi,
    mods        = ~ Preregistered,
    data        = dat.lehmann2018,
    measure     = "SMD",
    only_priors = TRUE
  )

  print_prior(priors)
  print_prior(priors, parameter = "mu")
  print_prior(priors, parameter = "tau")
  print_prior(priors, parameter = "Preregistered", component = "mods")
}
} # }
```
