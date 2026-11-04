# Bayesian Model-Averaged Multivariate Meta-Analysis

Fits Bayesian model-averaged multivariate and multilevel meta-analytic
models with known sampling covariance and formula random effects. Random
structures supplied through `random` are model-averaged with independent
component gates.

## Usage

``` r
BMA.mv(
  yi,
  V,
  ni,
  mods,
  scale,
  random,
  R = NULL,
  Rscale = "cor",
  data,
  slab,
  subset,
  measure,
  prior_effect,
  prior_heterogeneity,
  prior_mods,
  prior_scale,
  prior_effect_null,
  prior_heterogeneity_null,
  prior_mods_null,
  prior_scale_null,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "meandif",
  prior_unit_information_sd,
  rescale_priors = 1,
  prior_informed_field,
  prior_informed_subfield,
  known_v_parameterization = "auto",
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
  ...,
  vi = NULL,
  sei = NULL
)
```

## Arguments

- yi:

  a vector of effect sizes, or a formula with the effect size on the
  left-hand side and location moderators on the right-hand side (for
  example `yi ~ x1 + x2`). If a formula is supplied, `mods` must not be
  specified.

- V:

  a known working variance-covariance matrix, a list of block
  variance-covariance matrices, or an exact diagonal-plus-factor
  declaration created by
  [`known_v_factor()`](https://fbartos.github.io/RoBMA/reference/known_v_factor.md),
  used by
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).
  Ordinary covariance matrices from
  [`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
  can be supplied directly. Product selection requires no publication
  groups. For best-p-value selection, declare groups separately with
  `selection_model(weight_rule = "best", group = ...)`. Matrix symmetry
  uses base R's numerical tolerance. Accepted matrices are stored
  symmetrically by averaging off-diagonal pairs; the supplied object and
  its diagonal are unchanged. Dependency blocks that are block constant
  on nested partitions are represented exactly as a diagonal-plus-factor
  structure without a
  [`known_v_factor()`](https://fbartos.github.io/RoBMA/reference/known_v_factor.md)
  declaration. Dimensions and row alignment are checked on the complete
  input; when rows are dropped (`subset` or missing values), the
  numerical covariance checks apply to the retained rows only, so
  acceptance says nothing about the entries of the dropped rows.

- ni:

  an optional vector of sample sizes. Used for `measure = "GEN"` or when
  estimating the unit information standard deviation.

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

- random:

  Optional BayesTools random-effect formula or list of formulas.
  Omission or `NULL` means a fixed-effect model with no implicit
  heterogeneity term. Known sampling covariance `V` remains unchanged.
  Declare `random` before supplying heterogeneity priors or scale
  formulas.

- R:

  optional known random-effect group covariance/correlation matrix, or a
  named list of matrices. Matrix row and column names identify grouping
  levels. This is distinct from the known sampling covariance `V`.

- Rscale:

  scaling applied to each supplied `R` matrix before fitting. Supported
  values are `"cor"`, `"none"`, `"cor0"`, `"cov0"`, with
  metafor-compatible aliases `TRUE`, `FALSE`, `0`, `1`, `2`, and `3`. An
  unnamed scalar is recycled across all supplied matrices. A vector or
  list can specify one scaling value per named `R` entry, with exactly
  matching names. An unnamed one-element list is supported only for a
  single `R` matrix.

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

- known_v_parameterization:

  known-`V` backend used by
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).
  `"auto"` chooses an exact backend when feasible; `"latent"` uses a
  latent `D + BB'` decomposition; `"whitened"` uses an eigen-rotated
  normal likelihood; `"block_mvn"` uses an exact native block
  multivariate-normal likelihood.

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

- vi:

  a vector of sampling variances. Either `vi` or `sei` must be supplied
  for normal models.

- sei:

  a vector of standard errors. Either `vi` or `sei` must be supplied for
  normal models.

## Value

A fitted object of class
`c("BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")`.

## Details

Omitted or NULL `random` specifies no heterogeneity and no random-effect
inclusion mixture. Heterogeneity priors and scale formulas require an
explicit random structure; the known sampling covariance remains
unchanged.

`BMA.mv()` combines the product-space model-averaging workflow of
[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md) with the
likelihood, formula, random-effect, and known-covariance machinery of
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md). For
multiple top-level random components, a Dirichlet prior allocates the
slab total variance and an independent Bernoulli gate multiplies each
allocated contribution. The allocation is not renormalized when a gate
is off. Thus component \\j\\ contributes \$\$I_j w_j \tau^2,\$\$ where
\\I_j\\ is its inclusion indicator, \\w_j\\ is its slab allocation, and
\\\tau\\ is the slab total SD. With one component, no artificial
Dirichlet weight is created and the contribution is \\I_1 \tau^2\\.

By default, every top-level component has an independent 0.5 inclusion
probability. `prior_heterogeneity` defines the shared positive slab and
`prior_heterogeneity_null` defines exclusion at exactly zero. Their
prior weights determine the common gate probability; multiple
alternative priors retain their relative weights within the slab.
Setting `prior_heterogeneity_null = NULL` or `FALSE` fixes all component
gates on, while `prior_heterogeneity = NULL` or `FALSE` fixes them off.
Nested random terms are gated as one top-level component before their
internal Dirichlet split. Component-specific or partial scale models
instead have independent SD roots, with allocation within each
component; unscaled components keep their own base SD prior. See
[`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/reference/random_effect_prior_specification.md).

A partial
[`BayesTools::prior_random()`](https://fbartos.github.io/BayesTools/reference/prior_random.html)
supplied through `prior_heterogeneity` may override contrasts,
correlations, monitoring, and parameterization. `BMA.mv()` owns the
gated scale architecture, so a custom
[`prior_random()`](https://fbartos.github.io/RoBMA/reference/prior_random.md)
must not supply an SD, SD source, term-specific SD, covariance-owned SD,
or variance allocation.

Random-effect inclusion probabilities are reported separately from fixed
effect and scale-regression inclusion. Averaged random-component SDs
retain the excluded zero branch; conditional random summaries condition
each component on its own gate. `tau_total` and `tau2_total` are the
realized gated aggregate, including zero when every component is
excluded. `tau2_prop(j)` is the realized share \\I_j w_j / \sum_k I_k
w_k\\, conditional on positive total heterogeneity. Excluded components
therefore have zero share, and the all-off branch is omitted only from
variance-proportion summaries. The positive slab scale and raw Dirichlet
weights remain internal coordinates.

Product-space marginal likelihood and bridge-sampling methods are not
available. Predictive comparison through
[`loo.brma()`](https://fbartos.github.io/RoBMA/reference/loo.brma.md)
and
[`waic.brma()`](https://fbartos.github.io/RoBMA/reference/waic.brma.md)
remains available. If `V` is singular, every allowed product-space
branch must contain structural variance sufficient to regularize its
null space; an allowed all-components-off branch is therefore rejected.

The function does not add publication-bias models. Use
[`RoBMA.mv()`](https://fbartos.github.io/RoBMA/reference/RoBMA.mv.md) to
add model-averaged selection, PET, and PEESE publication-bias
adjustments.

## See also

[`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md),
[`RoBMA.mv()`](https://fbartos.github.io/RoBMA/reference/RoBMA.mv.md),
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`summary_models()`](https://fbartos.github.io/RoBMA/reference/summary_models.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("dat.assink2016", package = "metadat")
V <- metafor::vcalc(
  vi,
  cluster = study,
  type = deltype,
  obs = esid,
  rho = c(0.7, 0.5),
  data = dat.assink2016
)
fit <- BMA.mv(
  yi = yi,
  V = V,
  mods = ~ deltype,
  random = ~ 1 | study / esid,
  data = dat.assink2016,
  measure = "SMD",
  seed = 1,
  silent = TRUE
)
summary(fit)
summary_models(fit)
} # }
```
