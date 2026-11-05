# Bayesian Model-Averaged Generalized Meta-Analysis

Function for fitting Bayesian model-averaged meta-analytic models
directly to binary or count data using a generalized linear mixed model
(GLMM) framework. Unlike
[`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md), this
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
  set_contrast_factor_predictors = "meandif",
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
  [`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md).

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

  prior distribution(s) for the alternative effect component(s).

- prior_heterogeneity:

  prior distribution(s) for the alternative heterogeneity component(s).

- prior_mods:

  prior distribution(s) for alternative moderator components. A single
  prior applies to all terms; a named list can specify term-specific
  components. This argument can be specified only when `mods` is
  specified. For a formula without an intercept (`~ 0 + x`), an
  `intercept` entry other than a point prior at zero stops with an error
  of class `RoBMA_prior_intercept_conflict` (also of class
  `RoBMA_prior_error`). A named alternative `intercept` must be a single
  simple or point prior; it is the complete intercept prior and replaces
  the global effect mixture. Specify it in one place: null-list
  intercepts and explicit global effect declarations left without a
  consumer are refused. Ordinary term entries replace only their
  alternative component.

- prior_scale:

  prior distribution(s) for alternative scale-regression components. A
  single prior applies to all terms; a named list can specify
  term-specific components. This argument can be specified only when
  `scale` is specified. A named alternative `intercept` must be a single
  positive simple or point prior; it is the complete SD prior used by
  `log(intercept)`. Global heterogeneity declarations remain meaningful
  when they supply another random root or its inclusion weights;
  declarations replaced without a consumer are refused. Null-list
  intercept declarations are unavailable.

- prior_heterogeneity_allocation:

  prior distribution(s) for the alternative cluster-level heterogeneity
  allocation component(s). This argument can be specified only when
  `cluster` is specified.

- prior_baserate:

  prior distribution for the estimate-specific midpoint base-rate
  probability in binomial GLMM models (`measure = "OR"`). If omitted or
  `NULL`, defaults to independent `Beta(1, 1)` priors. Continuous priors
  are restricted to `[0, 1]` (see
  [`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md))
  and point priors must lie strictly inside `(0, 1)`. The estimate-unit
  log-likelihood of the model
  ([`log_lik()`](https://fbartos.github.io/RoBMA/dev/reference/log_lik.md),
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
  ([`log_lik()`](https://fbartos.github.io/RoBMA/dev/reference/log_lik.md),
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

- prior_effect_null:

  prior distribution(s) for the null effect component(s).

- prior_heterogeneity_null:

  prior distribution(s) for the null heterogeneity component(s).

- prior_mods_null:

  prior distribution(s) for null moderator components. A named
  `intercept` entry, including `NULL` or `FALSE`, is unavailable: put
  the complete intercept prior in `prior_mods`, or use global prior
  arguments without a formula-list intercept. The same rule applies to
  `prior_scale_null`. A single prior applies to all terms; a named list
  can specify term-specific components. This argument can be specified
  only when `mods` is specified. The null and the alternative prior of a
  factor term need the same contrast; a conflict stops with an error of
  class `RoBMA_prior_contrast_conflict` (also of class
  `RoBMA_prior_error`) that names the term and the remedies. The same
  rule applies to `prior_scale_null`.

- prior_scale_null:

  prior distribution(s) for null scale-regression components. A single
  prior applies to all terms; a named list can specify term-specific
  components. This argument can be specified only when `scale` is
  specified.

- prior_heterogeneity_allocation_null:

  prior distribution(s) for the null cluster-level heterogeneity
  allocation component(s). This argument can be specified only when
  `cluster` is specified.

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
  [`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
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
  [`set_autofit_control()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
  for details.

- convergence_checks:

  list of convergence check settings. See
  [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_control.md)
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
[`brma.glmm`](https://fbartos.github.io/RoBMA/dev/reference/brma.glmm.md)
with the mixture prior specification of
[`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) for
Bayesian model-averaging.

Model for odds ratios (`measure = "OR"`) uses a binomial-normal model as
described in Jackson et al. (2018) .

Model for incidence rate ratios (`measure = "IRR"`) uses a
Poisson-normal model as described in Bagos and Nikolopoulos (2009) .

When `weights` are supplied, they are treated as likelihood weights on
the paired two-arm study contribution.

For factor moderators, `BMA.glmm()` defaults to
`set_contrast_factor_predictors = "meandif"`, matching the other
model-averaging constructors. This differs from
[`brma.glmm()`](https://fbartos.github.io/RoBMA/dev/reference/brma.glmm.md),
which defaults to `"treatment"` as a single-model constructor.

## See also

[`brma.glmm()`](https://fbartos.github.io/RoBMA/dev/reference/brma.glmm.md)
[`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
[`summary.brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma.md)
[`plot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.brma.md)

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
