# Plot Prior Distributions

`plot_prior` visualizes prior distributions stored in `brma`, `BMA`, and
`RoBMA` objects. This is especially useful for objects created with
`only_priors = TRUE`.

## Usage

``` r
plot_prior(x, ...)

# S3 method for class 'prior'
plot_prior(x, plot_type = "base", ...)

# S3 method for class 'brma'
plot_prior(
  x,
  parameter = "mu",
  parameter_mods,
  parameter_scale,
  standardized_coefficients = TRUE,
  output_measure = NULL,
  transform = NULL,
  plot_type = "base",
  ...
)
```

## Arguments

- x:

  a `brma`, `BMA`, or `RoBMA` object, or a prior distribution object.

- ...:

  additional arguments passed to the prior plotting method.

- plot_type:

  whether to use a base plot `"base"` or ggplot2 `"ggplot"` for
  plotting. Defaults to `"base"`.

- parameter:

  character. Base parameter to plot. Defaults to `"mu"`. Common options
  are `"mu"`, `"tau"`, `"rho"`, `"PET"`, `"PEESE"`, `"omega"`, and
  `"bias"`, with aliases `"effect"` = `"mu"`, `"heterogeneity"` =
  `"tau"`, and `"weightfunction"` = `"omega"`. `"bias"` plots only
  non-mixed or homogeneous bias priors; for mixed weightfunction and
  PET/PEESE mixtures use `"omega"`, `"PET"`, or `"PEESE"`. Moderator and
  scale terms can also be selected by name when unambiguous. A character
  vector requests multiple base parameters.

- parameter_mods:

  character. Moderator term to plot. Use `"intercept"` for the adjusted
  effect in meta-regression models.

- parameter_scale:

  character. Scale-regression term to plot. Use `"intercept"` for the
  heterogeneity intercept in location-scale models.

- standardized_coefficients:

  whether to plot moderator and scale-regression priors on the
  standardized predictor scale. Defaults to `TRUE`, which shows the
  priors as specified. Set to `FALSE` to transform them to the original
  predictor scale when continuous predictors were standardized.

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

## Value

`plot_prior` returns either `NULL` invisibly if `plot_type = "base"` or
a `ggplot2` object if `plot_type = "ggplot"`. If multiple parameters are
requested, a named list is returned, invisibly for base plots.

## Details

`output_measure` and `transform` transform the prior plotting scale only
for effect-size location priors (`"mu"` or the meta-regression
intercept).

## See also

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
    yi           = yi,
    vi           = vi,
    mods         = ~ Preregistered,
    data         = dat.lehmann2018,
    measure      = "SMD",
    only_priors  = TRUE
  )

  plot_prior(priors, parameter = "mu")
  plot_prior(priors, parameter = "tau")
  plot_prior(priors, parameter_mods = "Preregistered")
}
} # }
```
