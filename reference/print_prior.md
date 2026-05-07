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
print_prior(x, parameter, parameter_mods, parameter_scale, ...)
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
  selected by name when unambiguous. A character vector requests
  multiple base parameters.

- parameter_mods:

  character. Moderator term to print. Use `"intercept"` for the adjusted
  effect in meta-regression models.

- parameter_scale:

  character. Scale-regression term to print. Use `"intercept"` for the
  heterogeneity intercept in location-scale models.

## Value

`print_prior` invisibly returns the selected prior distribution. If
multiple parameters are requested, a named list of prior distributions
is returned invisibly.

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
  print_prior(priors, parameter_mods = "Preregistered")
}
} # }
```
