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
  selection = BayesTools::selection_model(),
  selection_control = set_selection_likelihood_control(),
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

  an optional vector of finite positive likelihood weights. For
  normal/effect-size models, each weight powers the estimate likelihood.
  For constructors with GLMM raw-count input, each weight powers the
  paired two-arm likelihood for one study. An integer weight \\w\\
  repeats that conditional likelihood contribution \\w\\ times, with the
  declared latent-effect and grouping structure held fixed. Positive
  fractional likelihood powers are also supported.

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
  [`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md)
  for details on specifying publication-bias priors and
  [`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md)
  for details on specifying meta-analytic parameter priors.

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

- prior_bias_null:

  prior distribution(s) for null publication-bias component(s), usually
  [`prior_none()`](https://fbartos.github.io/RoBMA/dev/reference/prior_none.md).
  See
  [`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md).

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

- selection:

  specification created by
  [`selection_model()`](https://fbartos.github.io/RoBMA/dev/reference/selection_model.md)
  for automatically constructed weightfunction priors. The default
  integrates estimate-level random effects and the complete sampling
  error, and conditions on other random effects. Each source can instead
  be conditioned upon or integrated through the corresponding
  selection-model setting. The default `weight_rule = "product"` needs
  no publication groups and ignores `group`. Publication grouping is
  active only for `weight_rule = "best"`. Explicit `prior_bias` objects
  retain their own selection specification; this argument does not
  overwrite it. Conditioning on sampling variation retains the entire
  sampling-error realization, including in ordinary univariate models.
  This sampling setting requires omitted or unit observation `weights`.

- selection_control:

  numerical integration settings created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/dev/reference/set_selection_likelihood_control.md).
  These affect numerical evaluation, without changing the conditioning
  model or the vector weighting rule.

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
[`funnel()`](https://fbartos.github.io/RoBMA/dev/reference/funnel.md),
[`add_loo()`](https://fbartos.github.io/RoBMA/dev/reference/add_loo.brma.md),
and related methods.

## Details

`RoBMA()` uses product-space Bayesian model averaging. Inclusion Bayes
factors and model-averaged estimates are obtained from mixture priors
for effect, heterogeneity, moderators, scale regression, and
publication-bias components.

By default, `model_type = "PSMA"` includes selection-model weight
functions together with PET and PEESE publication-bias adjustments. Use
[`BMA()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.md) for
model averaging without publication-bias adjustment, or
[`brma()`](https://fbartos.github.io/RoBMA/dev/reference/brma.md) for
fitting a single meta-analytic model.

`RoBMA()` uses normal/effect-size input (`yi` with `vi` or `sei`).
Raw-count GLMM model averaging is provided by
[`BMA.glmm()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.glmm.md).

`selection` configures every generated default weightfunction. Its
default integrates estimate-level random effects and the complete
sampling error, and conditions on other random effects. Explicit priors
retain their nested
[`selection_model()`](https://fbartos.github.io/RoBMA/dev/reference/selection_model.md)
settings. Active selection branches must share the same conditioning
choices. Only `weight_rule = "best"` branches require and must share a
publication partition. Product branches ignore `group`; an ensemble
without `"best"` needs no publication groups and accepts `group = NULL`.
Branches can differ in weight-height priors, bins, and `weight_rule`.
These fixed model choices do not introduce new inclusion parameters.

Product-space objects support predictive comparison with
[`add_loo()`](https://fbartos.github.io/RoBMA/dev/reference/add_loo.brma.md)
and
[`add_waic()`](https://fbartos.github.io/RoBMA/dev/reference/add_waic.brma.md).
Bridge-sampling marginal likelihood via
[`add_marglik()`](https://fbartos.github.io/RoBMA/dev/reference/add_marglik.brma.md)
is not available for product-space model-averaging objects.

## See also

[publication_bias_prior_specification](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md),
[`RoBMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md),
[`BMA()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.md),
[`brma()`](https://fbartos.github.io/RoBMA/dev/reference/brma.md),
[`bselmodel()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.md),
[`bPET()`](https://fbartos.github.io/RoBMA/dev/reference/bPET.md),
[`bPEESE()`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/dev/reference/summary.brma.md),
[`plot.brma()`](https://fbartos.github.io/RoBMA/dev/reference/plot.brma.md)

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
