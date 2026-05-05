# Bayesian Model-Averaged Generalized Meta-Analysis

Function for fitting Bayesian model-averaged meta-analytic models
directly to binary or count data using a generalized linear mixed model
(GLMM) framework. Unlike
[`RoBMA`](https://fbartos.github.io/RoBMA/reference/RoBMA.md), this
function does not adjust for publication bias, as weight function and
regression-based bias adjustment methods are not available for GLMM
outcomes.

## Usage

``` r
BMA.glmm(
  ai,
  bi,
  ci,
  di,
  n1i,
  n2i,
  x1i,
  x2i,
  t1i,
  t2i,
  weights,
  mods,
  scale,
  cluster,
  data,
  slab,
  subset,
  measure = "OR",
  prior_effect,
  prior_heterogeneity,
  prior_mods,
  prior_scale,
  prior_heterogeneity_allocation,
  prior_baserate,
  prior_lograte,
  prior_effect_null,
  prior_heterogeneity_null,
  prior_mods_null,
  prior_scale_null,
  prior_heterogeneity_allocation_null,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd,
  rescale_priors = 1,
  prior_informed_field,
  prior_informed_subfield,
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

- ai:

  a vector of the number of events in the treatment or experimental
  group for binomial GLMM models.

- bi:

  a vector of the number of non-events in the treatment or experimental
  group for binomial GLMM models.

- ci:

  a vector of the number of events in the control group for binomial
  GLMM models.

- di:

  a vector of the number of non-events in the control group for binomial
  GLMM models.

- n1i:

  a vector of the sample size in the treatment or experimental group. If
  omitted for binomial GLMMs, it is computed as `ai + bi`.

- n2i:

  a vector of the sample size in the control group. If omitted for
  binomial GLMMs, it is computed as `ci + di`.

- x1i:

  a vector of the number of events in the treatment/experimental group
  (for Poisson data).

- x2i:

  a vector of the number of events in the control group (for Poisson
  data).

- t1i:

  a vector of the person-time in the treatment/experimental group.

- t2i:

  a vector of the person-time in the control group.

- weights:

  an optional vector of positive likelihood weights. For
  normal/effect-size models, each weight powers the estimate likelihood.
  For constructors with GLMM raw-count input, each weight powers the
  paired two-arm likelihood for one study.

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

- prior_baserate:

  prior distribution for the estimate-specific midpoint base-rate
  probability in binomial GLMM models. If omitted or `NULL`, defaults to
  independent `Beta(1, 1)` priors.

- prior_lograte:

  prior distribution for the estimate-specific midpoint log-rate in
  Poisson GLMM models. If omitted or `NULL`, a data-based
  unit-information normal prior is used independently for each estimate.

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
  none priors are unchanged; publication-bias priors are not rescaled
  except for the default PEESE prior's UISD adjustment. Defaults to 1.

- prior_informed_field:

  character. The field of the informed prior distributions. Omit to use
  the standard default prior specification; explicit `NULL` is invalid.

- prior_informed_subfield:

  character. The subfield of the informed prior distributions. Omit to
  use the field-specific default, such as `"Cochrane"` for
  `prior_informed_field = "medicine"`; explicit `NULL` is invalid.

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

A fitted object of class `c("BMA.glmm", "RoBMA", "brma.glmm", "brma")`.
The object contains checked `data`, checked mixture `priors`, the JAGS
`fit`, cached `summary`, and cached `coefficients`.

## Details

`BMA.glmm` combines the data input style of
[`brma.glmm`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
with the mixture prior specification of
[`RoBMA`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) for
Bayesian model-averaging.

Model for odds ratios (`measure = "OR"`) uses a binomial-normal model as
described in Jackson et al. (2018) .

Model for incidence rate ratios (`measure = "IRR"`) uses a
Poisson-normal model as described in Bagos and Nikolopoulos (2009) .

When `weights` are supplied, they are treated as likelihood weights on
the paired two-arm study contribution.

## See also

[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md)
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)
[`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.bcg, package = "metadat")

  fit <- BMA.glmm(
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    data    = dat.bcg,
    measure = "OR",
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
}
} # }
```
