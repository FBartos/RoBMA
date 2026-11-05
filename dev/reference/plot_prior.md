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
  component = "auto",
  ...
)
```

## Arguments

- x:

  a `brma`, `BMA`, or `RoBMA` object, or a prior distribution object.

- ...:

  additional arguments passed to the prior plotting method. For mixed
  continuous and point priors, `ylim` controls the density axis, `ylim2`
  controls the probability-mass axis, and `ylab2` controls its label.

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
  vector requests multiple base parameters. A quantity of the fitted
  model that is no model parameter (e.g. an inclusion indicator or a
  weight-function coordinate such as `"omega[1]"`) stops with the
  classes `BayesTools_hypothesis_target` and
  `BayesTools_hypothesis_unavailable`, as in
  [`plot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.brma.md).

- parameter_mods:

  legacy moderator selector. Prefer `parameter` with
  `component = "mods"`. Use `"intercept"` for the adjusted effect in
  meta-regression models.

- parameter_scale:

  legacy scale-regression selector. Prefer `parameter` with
  `component = "scale"`. Use `"intercept"` for the heterogeneity
  intercept in location-scale models.

- standardized_coefficients:

  whether to plot moderator and scale-regression priors on the
  standardized predictor scale. Defaults to `TRUE`, which shows the
  priors as specified. Set to `FALSE` to transform them to the original
  predictor scale when continuous predictors were standardized (see
  Details for terms without an original-scale prior).

- output_measure:

  effect-size measure for location/effect predictions. Defaults to the
  fitted measure. Supported conversions are among `"SMD"`, `"COR"`,
  `"ZCOR"`, and `"OR"`; `"RR"`, `"HR"`, `"IRR"`, `"RD"`, and `"GEN"` can
  only be returned on their fitted measure. Use `transform = "EXP"` for
  ratio-scale output from log-scale measures.

- transform:

  optional display transformation. Currently `"EXP"` exponentiates
  log-scale measures `"OR"`, `"RR"`, `"HR"`, and `"IRR"`.

- component:

  parameter component. Defaults to `"auto"`, which infers the component
  when possible. Use `"mods"` (alias `"location"`), `"scale"`, `"bias"`,
  or `"random"` (the random-effect quantities of fitted objects) to
  disambiguate terms used in multiple prior components. A
  `parameter_mods` or `parameter_scale` selection with another
  `component` stops with an error of class `RoBMA_component_mismatch`.

## Value

`plot_prior` returns the selected prior plotting method's result,
invisibly for base plots (which may include plot-data metadata), or a
`ggplot2` object for `plot_type = "ggplot"`. If multiple parameters are
requested, a named list is returned, invisibly for base plots.

## Details

`output_measure` and `transform` transform the prior plotting scale only
for effect-size location priors (`"mu"` or the meta-regression
intercept). For base plots, the probability mapping created by
`plot_prior()` is reused by subsequent posterior
[`lines.brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.brma.md)
calls, including overlays from other fitted objects. PET and PEESE
priors, single or mixed, show bias-only effect-versus-standard-error
curves by default, with a zero location contribution, the effect
direction of the object (override it with `effect_direction`), and the
declared model weights. Use `individual = TRUE` for coefficient
densities, or supply `prior_list_mu` to request full regression curves
with explicit location priors, as supported by
[`BayesTools::plot_prior_list()`](https://fbartos.github.io/BayesTools/reference/plot_prior_list.html).

The random-effect quantities of fitted `brma.mv`, `BMA.mv`, and
`RoBMA.mv` objects are selected by their names in `parameter`, such as
`"tau"`, `"tau_total"`, `"tau2"`, `"tau2_prop(study)"`, or
`"study: tau"` (`"tau"` stops as ambiguous when several blocks have an
SD). A quantity that is one declared prior (the SD of a block, the total
SD of a variance allocation, or a correlation with its own prior) shows
that prior. Any other quantity (a variance, a share of the total
variance, the SD of a component of an allocation, and the quantities of
gated allocations) shows the prior law it follows from the declared
priors, the density that
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) overlays with
`prior = TRUE`, drawn without a posterior; point masses (for example of
inclusion gates) are drawn as arrows. For these laws, `xlim`,
`n_points`, `par_name`, and graphical arguments apply as for other
priors, while `x_seq` and `x_range_quant` stop. A quantity without a
prior density stops with the classes `BayesTools_hypothesis_target` and
`BayesTools_hypothesis_unavailable`. `standardized_coefficients` selects
the fitted (default) or the original predictor scale of these laws.
Objects created with `only_priors = TRUE` have no fitted random-effect
quantities; print their random-effect priors with
`print_prior(x, parameter = "random")`.

Inconsistent random-allocation parent metadata stop with
`RoBMA_random_prior_derivation_error`; cycles additionally carry
`RoBMA_random_allocation_cycle`, and unknown parents carry
`RoBMA_random_allocation_parent_unavailable`. A parent chain beyond the
supported depth of 16 allocations carries
`RoBMA_random_allocation_depth_unavailable`. The condition retains its
`reason` and visited `allocation_chain`. Other source errors propagate
with their original classes and causes.

A term with several fitted coordinates (a factor with three or more
levels or with independent contrasts, and the interactions of such a
factor) has an original-scale prior only when the standardization of
continuous predictors leaves its coordinates unchanged. With
`standardized_coefficients = FALSE`, the prior of such a term whose
coordinates the standardization changes (a factor or factor-by-factor
interaction next to its interaction with a standardized predictor, and
that interaction term) stops with an error of class
`RoBMA_density_method_original_scale` and its parent
`RoBMA_density_method_unavailable`; use
`standardized_coefficients = TRUE`. The other terms (including factors
and factor-by-factor interactions without an interaction with a
standardized predictor) are plotted as before.

## See also

[`BMA()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.md)
[`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
[`brma()`](https://fbartos.github.io/RoBMA/dev/reference/brma.md)
[`prior()`](https://fbartos.github.io/RoBMA/dev/reference/prior.md)

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
  plot_prior(priors, parameter = "Preregistered", component = "mods")
}
} # }
```
