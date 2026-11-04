# Bayesian Generalized Meta-Analysis

Function for fitting random-effects, meta-regression, multilevel, and
location-scale meta-analytic models directly to either binary or count
data.

## Usage

``` r
brma.glmm(
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
  silent,
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

  an optional vector of finite positive likelihood weights. For
  normal/effect-size models, each weight powers the estimate likelihood.
  For constructors with GLMM raw-count input, each weight powers the
  paired two-arm likelihood for one study. An integer weight \\w\\
  repeats that conditional likelihood contribution \\w\\ times, with the
  declared latent-effect and grouping structure held fixed. Positive
  fractional likelihood powers are also supported.

- mods:

  an optional matrix, data.frame, or formula specifying location
  moderators (meta-regressors). Formula input is evaluated in `data` and
  supports explicit columns and standard formula operators, including
  `^` as R's crossing operator with a literal whole exponent of at least
  2: `(x + z)^2` expands to `x + z + x:z`, and `x^2` is not a quadratic
  term (`x + x^2` is `x`). Precompute transformed predictors, such as a
  squared predictor, as columns in `data`; inline calls such as
  [`I()`](https://rdrr.io/r/base/AsIs.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html), spline constructors,
  and user functions are not supported.

- scale:

  an optional matrix, data.frame, or formula specifying scale predictors
  for location-scale models. Formula input is evaluated in `data` and
  follows the same explicit-column formula grammar as `mods`.
  Random-formula models also accept a named list targeting top-level
  random components or concrete random-effect blocks; see
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).

- cluster:

  an optional vector of cluster identifiers for multilevel
  meta-analysis.

- data:

  an optional data frame containing the variables.

- slab:

  an optional vector of study labels.

- subset:

  an optional logical or numeric vector specifying a subset of data to
  be used. Subsetting precedes missing-row removal; retained-row outcome
  and predictor checks follow these steps. Type, shape, and existing
  extraction checks can still reject malformed inputs before subsetting.
  Row selection (by `subset` or missing-row removal) drops factor levels
  left without rows from `mods`, `scale`, and random-effect grouping
  factors. The index variable of a structured random-effect term
  (`cs()`, `hcs()`, `ar1()`, `har()` with one index variable) keeps its
  declared levels, which define the covariance structure; the variables
  of a multi-variable `cs(a + b | group)` index follow the ordinary
  rule, because that index is built from their observed combinations. If
  an index variable that keeps its levels is also used in `mods` or
  `scale` and loses one of its levels to the selection, the fit stops
  with an error of class `RoBMA_structured_index_level_removed` (also
  `RoBMA_data_error`); remove the level from the data with
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html) or use a
  separate column for the moderator.

- measure:

  a character string specifying the effect size measure.
  Normal/effect-size constructors require an explicit value and support
  `"SMD"`, `"ZCOR"`, `"RR"`, `"OR"`, `"HR"`, `"RD"`, `"IRR"`, and
  `"GEN"`. Use `"GEN"` only for general effect sizes without a known
  unit information standard deviation. GLMM raw-count constructors
  support only `"OR"` and `"IRR"` and default to `"OR"`.

- prior_effect:

  prior distribution for the effect size (\\\mu\\) parameter (the
  intercept). If omitted and consumed, a default prior is constructed.
  In single-model functions, explicit `NULL` or `FALSE` sets a spike at
  zero.

- prior_heterogeneity:

  prior distribution for the heterogeneity (\\\tau\\) parameter. If
  omitted and consumed, a default prior is constructed. In single-model
  functions, explicit `NULL` or `FALSE` sets a spike at zero. For
  [`brma.mv`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  random-effect formulas, this can also be a partial or complete
  [`prior_random`](https://fbartos.github.io/RoBMA/reference/prior_random.md)
  object; see
  [`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/reference/random_effect_prior_specification.md).

- prior_mods:

  prior distribution for the moderators (\\\beta\\) parameters. A single
  prior applies to all terms; a named list can specify term-specific
  priors. If omitted or `NULL`, default priors are used. This argument
  can be specified only when `mods` is specified. A formula without an
  intercept (`~ 0 + x`) fixes the location intercept at zero, so an
  `intercept` entry of the list is accepted only as a point prior at
  zero; any other `intercept` prior stops with an error of class
  `RoBMA_prior_intercept_conflict` (also of class `RoBMA_prior_error`).
  A named `intercept` entry supplies the complete intercept prior.
  Explicit competing global effect declarations are refused when this
  would discard them. Named `NULL` or `FALSE` intercept entries are
  invalid; specify the complete prior in one place.

- prior_scale:

  prior distribution for the scale (\\\delta\\) parameters. A single
  prior applies to all terms; a named list can specify term-specific
  priors. If omitted or `NULL`, default priors are used. This argument
  can be specified only when `scale` is specified. A named `intercept`
  entry supplies the complete positive SD prior used by
  `log(intercept)`. Global heterogeneity declarations remain meaningful
  when they supply another random root or its model-inclusion weights;
  declarations fully replaced without a consumer raise
  `RoBMA_prior_intercept_conflict`.

- prior_heterogeneity_allocation:

  prior distribution for the fraction of heterogeneity allocated to the
  cluster-level component in multilevel models (\\\rho\\, a variance
  fraction rather than a random-coefficient correlation). If omitted or
  `NULL`, defaults to `Beta(1, 1)`. This argument can be specified only
  when `cluster` is specified.

- prior_baserate:

  prior distribution for the estimate-specific midpoint base-rate
  probability in binomial GLMM models (`measure = "OR"`). If omitted or
  `NULL`, defaults to independent `Beta(1, 1)` priors. Continuous priors
  are restricted to `[0, 1]` (see
  [`prior_specification`](https://fbartos.github.io/RoBMA/reference/prior_specification.md))
  and point priors must lie strictly inside `(0, 1)`. The estimate-unit
  log-likelihood of the model
  ([`log_lik()`](https://fbartos.github.io/RoBMA/reference/log_lik.md),
  LOO, WAIC) integrates the base rate over any of these priors: point
  priors and untruncated beta priors adaptively (a uniform(0, 1) prior
  is a Beta(1, 1) prior), other families on a refined prior-CDF grid.
  The grid is accepted only when its log-likelihood changes by less than
  `1e-6` between refinements; any prior without an adaptive route (other
  than an untruncated beta prior; a truncated beta prior takes the grid
  too) can fail to reach this, especially with few events or very large
  arms, and then stops with an error of class
  `RoBMA_glmm_grid_nonconvergence`, as a lognormal base-rate prior does
  on cells with few events. A beta prior avoids it. The automatic
  computation of LOO/WAIC after fitting warns instead and returns the
  fit without them.

- prior_lograte:

  prior distribution for the estimate-specific midpoint log-rate in
  Poisson GLMM models (`measure = "IRR"`). If omitted or `NULL`, a
  normal prior centered at the pooled crude log-rate is used
  independently for each estimate. Its standard deviation defaults to
  `1` and can be changed via `RoBMA.options(default_lograte.sd = ...)`.
  Rescaling all person-times shifts the prior mean by the corresponding
  log conversion while leaving its standard deviation unchanged. The
  estimate-unit log-likelihood of the model
  ([`log_lik()`](https://fbartos.github.io/RoBMA/reference/log_lik.md),
  LOO, WAIC) integrates the log-rate over the prior: point priors and
  untruncated normal priors adaptively, other families on a refined
  prior-CDF grid. The grid is accepted only when its log-likelihood
  changes by less than `1e-6` between refinements; any prior without an
  adaptive route (other than an untruncated normal prior; a truncated
  normal prior takes the grid too) can fail to reach this, especially
  with few events or very large arms, and then stops with an error of
  class `RoBMA_glmm_grid_nonconvergence`, as gamma, lognormal, and
  inverse-gamma log-rate priors do on cells with few events. A normal
  prior avoids it. The automatic computation of LOO/WAIC after fitting
  warns instead and returns the fit without them.

- standardize_continuous_predictors:

  logical. Whether to standardize continuous predictors. Defaults to
  `TRUE`.

- set_contrast_factor_predictors:

  character. How to set contrast for factor predictors. Defaults are
  constructor-specific and shown in each function usage. Single-model
  constructors use `"treatment"`; model-averaging constructors use
  `"meandif"`. Accepted values are `"treatment"`, `"meandif"`,
  `"orthonormal"`, and `"independent"`. Independent coding creates one
  coefficient per factor level. Random-coefficient blocks reuse the
  resolved fixed-factor contrast unless `random_block(contrasts = ...)`
  overrides it. A warning is issued when default or inherited
  independent contrasts make a fixed or random design overspecified.
  With `"meandif"` or `"orthonormal"` contrasts, a formula term that
  interacts a factor with a predictor whose own term is missing, such as
  `g:x` in `~ g + g:x` or `~ g / x`, stops the fit: the term codes `g`
  by level indicators, with one coefficient per level, on which
  mean-difference and orthonormal priors are not defined. Specify an
  independent prior for that term, e.g.,
  `prior_mods = list("g:x" = prior_factor("normal", list(0, 1), contrast = "independent"))`
  (or `prior_scale` for a scale formula), include the missing term
  (`~ g * x`), or use `set_contrast_factor_predictors = "treatment"`.

- prior_unit_information_sd:

  numeric. The unit information standard deviation (\\\sigma\_{unit}\\).
  Cannot be used together with `prior_informed_field`.

- rescale_priors:

  numeric. A scaling factor for supported prior distributions. Point and
  none priors are unchanged; supplied
  [`prior_random()`](https://fbartos.github.io/RoBMA/reference/prior_random.md)
  scale architectures and publication-bias priors are not rescaled. The
  default PEESE prior separately uses its documented UISD adjustment,
  independently of `rescale_priors`. Defaults to 1.

- prior_informed_field:

  character. The field of the informed prior distributions. Omit to use
  the standard default prior specification; explicit `NULL` is invalid.

- prior_informed_subfield:

  character. The subfield of the informed prior distributions. Can only
  be specified together with `prior_informed_field`. Omit to use the
  field-specific default, such as `"Cochrane"` for
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

A fitted object of class `c("brma.glmm", "brma")`. The object contains
checked `data`, checked `priors`, the JAGS `fit`, cached `summary`, and
cached `coefficients`. If the corresponding package options are enabled,
it can also contain cached LOO, WAIC, or marginal likelihood results.

## Details

Model for odds ratios (`measure = "OR"`) corresponds to Model 4
described in Jackson et al. (2018) . `logit(pi[i])` is the
study-specific midpoint of the two arm logits. `prior_baserate` defines
the estimate-specific prior distribution on `pi[i]`.

Model for incidence rate ratios (`measure = "IRR"`) corresponds to Bagos
and Nikolopoulos (2009) . `phi[i]` is the study-specific midpoint of the
two arm log incidence rates. `prior_lograte` defines the
estimate-specific prior distribution on `phi[i]`. If unspecified, a
unit-information prior is based on the data and used independently for
each estimate.

When `weights` are supplied, they are treated as likelihood weights on
the paired two-arm study contribution.

## See also

[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md),
[`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.bcg, package = "metadat")

  fit <- brma.glmm(
    ai      = tpos,
    bi      = tneg,
    ci      = cpos,
    di      = cneg,
    mods    = ~ alloc,
    data    = dat.bcg,
    measure = "OR",
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
}
} # }
```
