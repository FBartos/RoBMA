# Bayesian Selection Model

Function for fitting random-effects, meta-regression, multilevel, and
location-scale meta-analytic selection models.

## Usage

``` r
bselmodel(
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
  prior_effect,
  prior_heterogeneity,
  prior_mods,
  prior_scale,
  prior_heterogeneity_allocation,
  prior_bias,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd,
  rescale_priors = 1,
  prior_informed_field,
  prior_informed_subfield,
  effect_direction = "detect",
  steps,
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

- prior_bias:

  selection-model bias prior, usually created by
  [`prior_weightfunction()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md).
  If omitted or `NULL`, a default one-sided weightfunction prior is
  constructed from `steps`.

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

- effect_direction:

  direction used by publication-bias adjustments. `"positive"` assumes
  statistically significant positive estimates are more likely to be
  selected; `"negative"` mirrors the selection direction; `"detect"`
  infers the direction from the fitted data.

- steps:

  numeric vector of one-sided p-value cut points for the default
  selection model. If `prior_bias` is supplied, the prior carries its
  own side, steps, and weights. If omitted, the default is `0.025`,
  yielding intervals `[0, .025]` and `(.025, 1]`.

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

- selection:

  specification created by
  [`selection_model()`](https://fbartos.github.io/RoBMA/reference/selection_model.md)
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
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md).
  These affect numerical evaluation, without changing the conditioning
  model or the vector weighting rule.

- ...:

  additional advanced arguments. Fitting functions reject unused
  arguments; currently recognized internal arguments include
  `only_data`, `only_priors`, `is_JASP`, and `is_JASP_prefix`.

## Value

A fitted object of class `c("bselmodel", "brma")` containing a single
Bayesian selection model fit.

## Details

`bselmodel()` is a normal/effect-size selection-model constructor.
Custom `prior_bias` can be a weightfunction prior or a supported
BayesTools selection-kernel prior; p-hacking kernels are not supported
in active RoBMA. By default, the model conditions on contextual cluster
effects and uses the product of estimate weights. For independent rows
conditional on these effects, this gives the usual selected-normal
likelihood with within-cluster heterogeneity integrated. Use
`selection = selection_model(...)` to choose `weight_rule = "best"`,
which uses the weight at the smallest p-value in each publication group.
For `"best"`, `group` identifies a data column; when `group = NULL`, the
specialized `cluster` argument supplies publication groups, or an
unclustered model uses one group per estimate. Product models use
neither input to define selection groups.

All conditioning choices use the same Gaussian source model. Conditioned
sources remain latent; their population distributions stay outside
selection normalization. Integrated sources may be reweighted by
selection. With all sources conditioned and positive weights, the
weights cancel and the observed law is the ordinary Gaussian model. See
[`bselmodel.mv()`](https://fbartos.github.io/RoBMA/reference/bselmodel.mv.md)
for the complete source and covariance contract.

Product normalizers factorize only for conditionally independent rows.
Dependent events use supported covariance-factor quadrature or fixed
randomized quasi-Monte Carlo integration with explicit error
diagnostics. Non-unit observation `weights` require
`known_sampling_variance = "integrate"` and independent product factors.
Unit weights are equivalent to omitting `weights`.

## See also

[publication_bias_prior_specification](https://fbartos.github.io/RoBMA/reference/publication_bias_prior_specification.md),
[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md),
[`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`funnel.brma()`](https://fbartos.github.io/RoBMA/reference/funnel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  fit <- bselmodel(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    steps   = 0.025,
    seed    = 1,
    silent  = TRUE
  )

  summary(fit)
  funnel(fit)
}
} # }
```
