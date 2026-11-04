# Bayesian Multivariate and Multilevel Selection Model

Fits estimate-level Bayesian selection models with known sampling
covariance and BayesTools random-effect formulas.

## Usage

``` r
bselmodel.mv(
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
  prior_bias,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd,
  rescale_priors = 1,
  prior_informed_field,
  prior_informed_subfield,
  effect_direction = "detect",
  steps,
  selection = BayesTools::selection_model(),
  selection_control = set_selection_likelihood_control(),
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
  silent,
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
`c("bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma")`.

## Details

Omitted or NULL `random` specifies a fixed-effect selection model with
no implicit heterogeneity. Both random-effect selection axes are then
inapplicable; the sampling-error selection axis still applies.
Heterogeneity priors and scale formulas require an explicit random
structure.

`bselmodel.mv()` combines the selection-model interface of
[`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md)
with the known-covariance and random-formula model of
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).
Each weightfunction carries a
[`selection_model()`](https://fbartos.github.io/RoBMA/reference/selection_model.md)
specifying whether estimate-level random effects, other random effects,
and the complete sampling error are conditioned upon or integrated
before selection normalization. The default integrates estimate-level
random effects and sampling error, conditions on other random effects,
and uses the product of estimate weights. The default product rule needs
no publication groups, so adding a three-level structure requires only
`random = ~ 1 | study / esid`. `weight_rule = "best"` instead uses the
weight at the smallest actual p-value in the publication group.

With one cutoff and all three sources integrated, `"best"` corresponds
to the relaxed report-all rule of van Aert et al. (in press) : a
favorable result can select the publication vector, after which all
measured outcomes are reported. The other conditioning cells and
multiple-cutoff best rules are extensions of that specification.

Publication grouping is active only for `weight_rule = "best"`. For
grouped multivariate models, specify it with
`selection = selection_model(weight_rule = "best", group = study)`, or
in the selection model of an explicit weightfunction prior. Independent
models without a random structure can use singleton groups
automatically. Covariance matrices, including standard
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
results, and random-effect grouping factors do not establish publication
identities. With `"product"`, `group` is ignored: it is neither bound to
a data column nor checked as a publication partition. Integration
follows the dependencies in `V` and the integrated random effects,
independently of publication grouping.

Sampling error is \\e \sim N(0,V)\\. With positive sampling standard
errors it can equivalently be written as \\e = Sz\\, with \\z \sim
N(0,R)\\ and \\V = SRS\\. Here \\S\\ contains the original sampling
standard errors and \\R\\ preserves their complete correlation
structure. Conditioning retains the entire realized vector \\e\\ during
selection normalization. No residual sampling variation remains inside
that normalizer. Integration instead contributes the full \\V\\ to
candidate variation. These choices apply equally to univariate,
diagonal, and correlated sampling covariance, including `vi` and `sei`
inputs.

The covariance `V` remains authoritative. A diagonal-plus-factor input
from
[`known_v_factor()`](https://fbartos.github.io/RoBMA/reference/known_v_factor.md)
can support more efficient calculations, but equivalent representations
define the same selection model with matching source settings and, for
`"best"`, publication groups. Selection thresholds always use the
original `sqrt(diag(V))`. The same faster routes are reached from an
ordinary matrix whose dependency blocks are block constant on nested
partitions, such as compound symmetry or nested
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
`type`/`obs` structures: that structure is recovered from the supplied
entries and reproduces them to working precision. It changes how the
same integral is evaluated, not which sampling sources are conditioned
upon.

The likelihood and post-fit methods use the same resolved conditioning
model. `known_v_parameterization` and numerical integration settings
cannot change it. The `estimate_random_effects` setting controls
estimate-level true-effect variation. Shared covariance-factor and dense
numerical paths use explicit integration diagnostics.

The source settings define the reporting model. An optional comparison
with integrating retained contexts can be requested with
[`selection_sensitivity_diagnostics()`](https://fbartos.github.io/RoBMA/reference/selection_sensitivity_diagnostics.md)
or stored with
[`add_selection_sensitivity_diagnostics()`](https://fbartos.github.io/RoBMA/reference/selection_sensitivity_diagnostics.md).
Fitting does not compute this separate comparison automatically.

`estimate_random_effects`, `other_random_effects`, and
`known_sampling_variance` each accept `"condition"` or `"integrate"`.
Both `"product"` and `"best"` support the resulting choices for
supported Gaussian sources, subject to the dependency restrictions
below.

Retained contexts remain unknown and are estimated. Their population
mixing laws stay outside the selection normalizer; integrating a context
instead allows selection to reweight its mixing law. An absent context
makes its two settings coincide. BayesTools identifies an estimate-level
random term by a distinct grouping level for every retained estimate. At
most one such term is allowed. Every other declared term is an
other-level source, including terms with independent coefficient
supports inside repeated groups. These roles do not split coefficient
families or discard known group covariance; each choice applies to the
complete declared random term.

With every outcome-generating source conditioned upon, candidate
outcomes are deterministic within a retained context. If weights are
positive almost surely, they cancel under conditional normalization: the
observed law is the ordinary Gaussian model and carries no information
about the selection weights. A zero acceptance probability on a
positive-probability set of retained contexts makes this selection
process invalid. It must not be handled by silently discarding those
contexts. Ordinary population-level publication selection instead
integrates all outcome-generating sources.

Product normalization retains all supported Gaussian dependencies. Best
weights require integrated dependencies to remain within each
publication group; an unsupported source is identified before fitting.
Conditioned sources may connect groups. Integrated candidate covariance
can be singular while the observed law remains nondegenerate after
averaging over retained sources. No residual noise is added to change
this boundary. Ensembles require a common conditioning cell across
active selection branches and a common publication partition across
their `"best"` branches only. Product branches impose no grouping
requirements; an ensemble without `"best"` accepts the default
`group = NULL`. Bins, weight priors, and supported weight rules may
differ across branches. Like the other
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)-family
constructors, `bselmodel.mv()` does not accept observation `weights`.

Joint likelihoods, single-model bridge evidence, conditional-density
estimation, selected prediction, latent summaries, and z-plots use this
resolved model. Estimate-unit LOO, WAIC, and LOO-PIT integrate deleted
outcomes within the original publication event. Their conditional row
scores must not be summed and interpreted as a joint group likelihood.
Cluster or block deletion remains unavailable for arbitrary
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
random formulas. Product-space ensembles do not provide bridge marginal
likelihoods. With conditioned sampling, estimate deletion integrates the
deleted sampling error conditional on the retained sampling errors of
the remaining estimates.

Marginal selected predictions draw new retained contexts from their
original mixing laws before drawing selected outcomes and reconstructing
latent true effects. At estimate depth, latent predictions draw fitted
true effects;
[`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md) and
[`fitted.brma()`](https://fbartos.github.io/RoBMA/reference/fitted.brma.md)
instead return conditional means. Response prediction at estimate depth
describes a new reporting event around fitted truth, with new sampling
error. Prediction's `conditioning_depth` is a separate choice from the
selection-model source settings. Cluster depth is available only for the
specialized `cluster` interface. Known-covariance latent and response
predictions with `newdata` require `V_new` and, for `"best"`, the
resolved publication-group reference; cross-covariance with fitted
outcomes and non-marginal `newdata` are unavailable. See
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
for prediction targets and labels.

When the complete sampling error is conditioned upon, total fitted truth
equals the observed estimate minus its fitted sampling error.
Conditional means and estimate-depth latent draws then coincide within
each posterior draw; uncertainty remains across posterior draws.
Blockwise
[`ranef()`](https://fbartos.github.io/RoBMA/reference/ranef.md) retains
its conditional-mean target for each random-effect block.

Partial-vector calculations retain the full selection event. Densities
need a nondegenerate observed Gaussian law, except for supported
one-coordinate projections of explicitly declared rank-one sources.
Multicoordinate singular Lebesgue densities and general singular
conditional latent posteriors are unavailable; zero integrated
true-effect covariance gives deterministic latent effects. No general
public multivariate CDF is provided. LOO-PIT and z-plots use their
respective deletion and marginal projection targets.

## References

van Aert RCM, Riley R, Jackson D (in press). “Correcting for publication
bias in multivariate and multilevel meta-analysis: A multivariate step
function selection model approach.” *Psychological Methods*.
[doi:10.31234/osf.io/rk4mn_v1](https://doi.org/10.31234/osf.io/rk4mn_v1)
.

## See also

[`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md),
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md),
[`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(yi = c(0.10, 0.20))
V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
fit <- bselmodel.mv(
  yi      = yi,
  V       = V,
  data    = dat,
  measure = "GEN",
  seed    = 1,
  silent  = TRUE
)
summary(fit)
} # }
```
