# Robust Bayesian Model-Averaged Meta-Analysis

Fits a robust Bayesian model-averaged meta-analysis. The default
ensemble averages across models with and without an effect,
heterogeneity, and publication-bias adjustment.

## Usage

``` r
RoBMA(
  yi,
  vi,
  sei,
  weights,
  ni,
  mods,
  scale,
  cluster,
  data,
  slab,
  subset,
  measure,
  effect_direction = "detect",
  prior_effect,
  prior_heterogeneity,
  prior_mods,
  prior_scale,
  prior_heterogeneity_allocation,
  prior_bias,
  prior_effect_null,
  prior_heterogeneity_null,
  prior_mods_null,
  prior_scale_null,
  prior_heterogeneity_allocation_null,
  prior_bias_null,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "meandif",
  prior_unit_information_sd,
  rescale_priors = 1,
  prior_informed_field,
  prior_informed_subfield,
  model_type = "PSMA",
  sample = 5000,
  burnin = 2000,
  adapt = 500,
  chains = 3,
  thin = 1,
  parallel = FALSE,
  autofit = FALSE,
  autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),
  seed = NULL,
  silent = TRUE,
  ...
)
```

## Arguments

- yi:

  a vector of effect sizes, or a formula with the effect size on the
  left-hand side and location moderators on the right-hand side (for
  example `yi ~ x1 + x2`). If a formula is supplied, `mods` must not be
  specified.

- vi:

  a vector of sampling variances. Either `vi` or `sei` must be supplied
  for normal models.

- sei:

  a vector of standard errors. Either `vi` or `sei` must be supplied for
  normal models.

- weights:

  an optional vector of positive likelihood weights. For
  normal/effect-size models, each weight powers the estimate likelihood.
  For constructors with GLMM raw-count input, each weight powers the
  paired two-arm likelihood for one study.

- ni:

  an optional vector of sample sizes. Used for `measure = "GEN"` or when
  estimating `"UISD"`).

- mods:

  an optional matrix, data.frame, or formula specifying location
  moderators (meta-regressors). Formula input is evaluated in `data`.

- scale:

  an optional matrix, data.frame, or formula specifying scale predictors
  for location-scale models. Formula input is evaluated in `data`.

- cluster:

  an optional vector of cluster identifiers for multilevel
  meta-analysis.

- data:

  an optional data frame containing the variables.

- slab:

  an optional vector of study labels.

- subset:

  an optional logical or numeric vector specifying a subset of data to
  be used.

- measure:

  a character string specifying the effect size measure.
  Normal/effect-size constructors require an explicit value and support
  `"SMD"`, `"ZCOR"`, `"RR"`, `"OR"`, `"HR"`, `"RD"`, `"IRR"`, and
  `"GEN"`. Use `"GEN"` only for general effect sizes without a known
  unit information standard deviation. GLMM raw-count constructors
  support only `"OR"` and `"IRR"` and default to `"OR"`.

- effect_direction:

  direction used by publication-bias adjustments. `"positive"` assumes
  statistically significant positive estimates are more likely to be
  selected; `"negative"` mirrors the selection direction; `"detect"`
  infers the direction from the fitted data.

- prior_effect:

  prior distribution(s) for the alternative effect component(s).

- prior_heterogeneity:

  prior distribution(s) for the alternative heterogeneity component(s).

- prior_mods:

  prior distribution(s) for alternative moderator components. A single
  prior applies to all terms; a named list can specify term-specific
  components.

- prior_scale:

  prior distribution(s) for alternative scale-regression components. A
  single prior applies to all terms; a named list can specify
  term-specific components.

- prior_heterogeneity_allocation:

  prior distribution(s) for the alternative cluster-level heterogeneity
  allocation component(s).

- prior_bias:

  prior distribution(s) for alternative publication-bias component(s),
  such as weight functions, PET, or PEESE.

  Alternative prior arguments can be:

  - A single prior distribution object (creates a mixture with one
    alternative)

  - A list of prior distributions (creates a mixture with multiple
    alternatives)

  - `NULL` or `FALSE` (omits the alternative hypothesis component)

  See
  [`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/reference/publication_bias_prior_specification.md)
  for details on specifying publication-bias priors and
  [`prior_specification`](https://fbartos.github.io/RoBMA/reference/prior_specification.md)
  for details on specifying meta-analytic parameter priors.

- prior_effect_null:

  prior distribution(s) for the null effect component(s).

- prior_heterogeneity_null:

  prior distribution(s) for the null heterogeneity component(s).

- prior_mods_null:

  prior distribution(s) for null moderator components. A single prior
  applies to all terms; a named list can specify term-specific
  components.

- prior_scale_null:

  prior distribution(s) for null scale-regression components. A single
  prior applies to all terms; a named list can specify term-specific
  components.

- prior_heterogeneity_allocation_null:

  prior distribution(s) for the null cluster-level heterogeneity
  allocation component(s).

- prior_bias_null:

  prior distribution(s) for null publication-bias component(s), usually
  [`prior_none()`](https://fbartos.github.io/RoBMA/reference/prior_none.md).
  See
  [`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/reference/publication_bias_prior_specification.md).

  Null prior arguments can be:

  - A single prior distribution object (creates a mixture with one null)

  - A list of prior distributions (creates a mixture with multiple
    nulls)

  - `NULL` or `FALSE` (omits the null hypothesis component)

  Defaults to a point mass (spike) at zero for effect and heterogeneity
  parameters.

- standardize_continuous_predictors:

  logical. Whether to standardize continuous predictors. Defaults to
  `TRUE`.

- set_contrast_factor_predictors:

  character. How to set contrast for factor predictors. Defaults are
  constructor-specific and shown in each function usage; single-model
  constructors use `"treatment"`, while model-averaging constructors use
  `"meandif"`.

- prior_unit_information_sd:

  numeric. The unit information standard deviation (\\\sigma\_{unit}\\).
  Cannot be used together with `prior_informed_field`.

- rescale_priors:

  numeric. A scaling factor for supported prior distributions. Point and
  none priors are unchanged. For constructors with publication-bias
  prior distributions, `rescale_priors` does not rescale them except for
  the default PEESE prior's UISD adjustment. Defaults to 1.

- prior_informed_field:

  character. The field of the informed prior distributions. Omit to use
  the standard default prior specification; explicit `NULL` is invalid.

- prior_informed_subfield:

  character. The subfield of the informed prior distributions. Omit to
  use the field-specific default, such as `"Cochrane"` for
  `prior_informed_field = "medicine"`; explicit `NULL` is invalid.

- model_type:

  character string specifying predefined publication-bias model
  ensembles. One of:

  - `"PSMA"` (default): Full RoBMA-PSMA ensemble with 6 weight
    functions + PET + PEESE

  - `"6w"`: Six weight function models

  - `"2w"`: Two weight function models

  - `"PP"`: PET-PEESE models only

  Custom `prior_bias` replaces the preset alternative bias components.
  If `prior_bias` is omitted, `model_type` determines the default
  alternatives even when `prior_bias_null` is customized or omitted.

- sample:

  numeric. Number of MCMC samples to save. Defaults to `5000`.

- burnin:

  numeric. Number of burn-in iterations. Defaults to `2000`.

- adapt:

  numeric. Number of adaptation iterations. Defaults to `500`.

- chains:

  numeric. Number of MCMC chains. Defaults to `3`.

- thin:

  numeric. Thinning interval. Defaults to `1`.

- parallel:

  logical. Whether to run MCMC chains in parallel. Defaults to `FALSE`.

- autofit:

  logical. Whether to automatically extend the MCMC chains if
  convergence is not met. Defaults to `FALSE`.

- autofit_control:

  list of autofit control settings. See
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  for details.

- convergence_checks:

  list of convergence check settings. See
  [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  for details.

- seed:

  numeric. Random seed for reproducibility. Defaults to `NULL`.

- silent:

  logical. Whether to suppress output. Constructors with no explicit
  default use `RoBMA.get_option("silent")` when `silent` is omitted.
  Model-averaging wrappers default to `TRUE` unless explicitly changed.

- ...:

  additional advanced arguments. Fitting functions reject unused
  arguments; currently recognized internal arguments include
  `only_data`, `only_priors`, `is_JASP`, and `is_JASP_prefix`.

## Value

A fitted object of class `c("RoBMA", "brma")`. The object contains
checked `data`, checked `priors`, the JAGS `fit`, cached `summary`, and
cached `coefficients`. It can be passed to
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`funnel()`](https://fbartos.github.io/RoBMA/reference/funnel.md),
[`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
and related methods.

## Details

`RoBMA()` uses product-space Bayesian model averaging. Inclusion Bayes
factors and model-averaged estimates are obtained from mixture priors
for effect, heterogeneity, moderators, scale regression, and
publication-bias components.

By default, `model_type = "PSMA"` includes selection-model weight
functions together with PET and PEESE publication-bias adjustments. Use
[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md) for model
averaging without publication-bias adjustment, or
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) for
fitting a single meta-analytic model.

`RoBMA()` uses normal/effect-size input (`yi` with `vi` or `sei`).
Raw-count GLMM model averaging is provided by
[`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md).

Product-space objects support predictive comparison with
[`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md)
and
[`add_waic()`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md).
Bridge-sampling marginal likelihood via
[`add_marglik()`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md)
is not available for product-space model-averaging objects.

## See also

[publication_bias_prior_specification](https://fbartos.github.io/RoBMA/reference/publication_bias_prior_specification.md),
[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md),
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md),
[`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md),
[`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md),
[`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- RoBMA(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
  plot(fit)
}
} # }
```
