# Bayesian Meta-Analysis With Known Sampling Covariance

Fits normal-likelihood Bayesian meta-analytic models with a known
working variance-covariance matrix `V`. The interface supports location
moderators, BayesTools random-effect formulas using [random-effect
formula structure
tags](https://fbartos.github.io/RoBMA/reference/random_effect_formula_tags.md),
and `brma.mv()`-specific known-`V` backends.

## Usage

``` r
brma.mv(
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
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
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
  used by `brma.mv()`. Ordinary covariance matrices from
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
  random components or concrete random-effect blocks; see `brma.mv()`.

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
  `brma.mv` random-effect formulas, this can also be a partial or
  complete
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

  known-`V` backend used by `brma.mv()`. `"auto"` chooses an exact
  backend when feasible; `"latent"` uses a latent `D + BB'`
  decomposition; `"whitened"` uses an eigen-rotated normal likelihood;
  `"block_mvn"` uses an exact native block multivariate-normal
  likelihood.

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

A fitted object of class `c("brma.mv", "brma.norm", "brma")`. The
advanced internal `only_data = TRUE` path returns checked data before
prior validation and random-effect marginalization decisions; use
`only_priors = TRUE` when inspecting the compiled random-effect
structure.

## Details

`brma.mv()` represents known sampling dependence through a latent
`D + BB'` decomposition of `V` (`known_v_parameterization = "latent"`),
through an eigen-rotated likelihood
(`known_v_parameterization = "whitened"`), or through an exact native
block multivariate-normal likelihood
(`known_v_parameterization = "block_mvn"`). The default
`known_v_parameterization = "auto"` chooses the whitened backend
whenever the extra diagonal variance is scalar and no scale regression
is present, otherwise chooses block-MVN when the largest dense
covariance block is at or below the configured feasibility threshold,
and falls back to the latent backend for large blocks. The automatic
block-MVN threshold defaults to 128 rows per covariance block and can be
changed with `options(RoBMA.known_v_block_mvn_max_block_size = value)`,
where `value` is a positive integer or `Inf`. Correlated `V` with scale
regression cannot use `"whitened"`; request `"block_mvn"` explicitly if
an automatic latent decomposition is unavailable for a nearly singular
block. Covariance inputs retain their symmetric representation.
Numerical positive-semidefinite checks use standardized correlations so
small marginal variances are not mistaken for null directions. At a
numerically singular boundary, the computational factor uses the
existing eigensolver roundoff tolerance; sampling, whitening, and
conditional calculations share its support. This numerical rule is
separate from the selection model's condition/integrate source choices.
Exact rank-one dependency blocks are represented by one latent sampling
factor, including when another backend was requested, so correlations of
one and any separate positive diagonal variance retain their
mathematical structure. Singular positive-semidefinite `V` matrices are
supported only when integrated conditional variance is structurally
positive on every retained null direction. This can come from strictly
positive heterogeneity or a supported marginalized random-effect
variance. The check is performed after formulas and priors are compiled.
Sampled random effects change the conditional mean and therefore do not
regularize the likelihood covariance. For the block-MVN backend, a fixed
point regularizer must also be large enough to make each block
numerically positive definite at double precision; ineffective fixed
values fail during preflight rather than inside JAGS. The latent
known-`V` backend is unavailable for other singular `V` matrices; use
`"auto"`, `"whitened"`, or `"block_mvn"`. Unsupported configurations
fail before fitting rather than invoking a degenerate-normal likelihood.
As a convenience for diagonal known-`V` inputs, `brma.mv()` also accepts
named `vi` or `sei` arguments when `V` is omitted; these are converted
to `diag(vi)` or `diag(sei^2)`. Do not supply `V` together with `vi` or
`sei`.

Omitted or NULL `random` declares no random-effects structure. No
implicit `tau` or heterogeneity mixture is added. Do not supply
`prior_heterogeneity` or `prior_heterogeneity_null` without `random`,
even as NULL or a zero prior; omit these irrelevant arguments instead.
Scale formulas also require an explicit random structure. This rule
applies to all .mv constructors; specialized univariate models retain
their usual implicit heterogeneity and clustered `tau`/`rho`/`I2`
behavior.

The `random` argument uses BayesTools [random-effect formula structure
tags](https://fbartos.github.io/RoBMA/reference/random_effect_formula_tags.md).
Plain `(expr | group)` syntax creates an unstructured random-coefficient
block, and `expr || group` creates a diagonal block. Therefore
`~ 1 | study`, `~ 1 + x | study`, and nested forms such as
`~ 1 | study/effect` are supported. The explicit `us()` / `un()`,
[`diag()`](https://rdrr.io/r/base/diag.html), and `id()` tags use the
same coefficient-formula semantics. The structured `cs()` / `hcs()`,
`ar1()` / [`ar()`](https://rdrr.io/r/stats/ar.html) / `har()`, and
`car()` tags instead own their index basis; see the linked tag
documentation before specifying factor slopes or index variables.
BayesTools stores fully explicit canonical random-effect names of the
form `(formula) owner: quantity(arguments)`. RoBMA uses their simplified
display aliases through RoBMA's meta-analytic naming map consistently in
summaries, plotting, density estimation, hypotheses, and posterior
draws. A sole random intercept drops its redundant argument, so a bare
block is displayed and accepted as `tau`, while an explicitly named
block is `study: tau`; the unique shorthand `tau` is also accepted. A
non-intercept coefficient remains explicit, for example `study: tau(x)`.
Bare formulas and unnamed one-entry lists omit a redundant owner;
explicitly named one-entry lists retain it, and missing names in longer
lists become `component 1`, `component 2`, and so on. Allocation names
include `tau_total`, `tau_common`, and `tau2_prop(study)`. BayesTools
retains its general `sd`/`var`/`cor` vocabulary internally; none of
these public RoBMA names expose backend coordinates. The optional `R`
argument supplies known covariance or correlation matrices across
random-effect grouping levels, following
[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
naming. `R` is separate from the known sampling covariance `V`: `V`
describes sampling-error dependence between estimates, while `R`
describes the latent group-axis covariance of random intercepts. Current
`R` support is limited to random-intercept blocks. `Rscale` controls how
each supplied `R` matrix is scaled before fitting and accepts `"cor"`,
`"none"`, `"cor0"`, `"cov0"`, plus metafor-compatible aliases `TRUE`/`1`
for `"cor"`, `FALSE`/`0` for `"none"`, `2` for `"cor0"`, and `3` for
`"cov0"`. In non-selection Gaussian models with `random`, shared sampled
random-effect blocks enter the conditional mean through BayesTools. A
supported one-to-one random-intercept block is automatically integrated
into the likelihood variance. Other terms remain sampled with their full
covariance. With a single top-level random component, `scale` models the
row-wise total random-effect SD consumed by that component. A named list
of `scale` formulas can instead target top-level random components or
concrete random-effect blocks. Block names and grouping labels are
accepted, and the terminal grouping variable is an alias for an
unambiguous nested block; for example, `scale = list(esid = ~ x)`
targets the `esid` block in `random = ~ 1 | study / esid`. Untargeted
blocks retain their ordinary SD priors. Component names can be ordinary
text; RoBMA normalizes them internally for JAGS. A single `scale`
formula is rejected for multiple top-level random components because the
target is ambiguous. In model- averaged models, block-specific scales
inside one grouped component require the blocks to be supplied as
separately named `random` entries so their inclusion gates remain
well-defined. When `scale` is a named list, `prior_scale` must either be
omitted or use the same named-list shape with one prior list per scale
component, using the same user-facing names. Post-fit scale predictions
for random-formula models keep these component scales separate rather
than collapsing them into one total \\\tau\\. In non-selection Gaussian
models, a supported random-intercept block that maps one-to-one to
estimates is automatically compiled as a marginalized variance component
instead of sampled latent coefficients. This preserves estimate-level
heterogeneity semantics while keeping shared higher-level random effects
conditional. For known `R`, this automatic marginalization is available
only when BayesTools validates a one-column random-intercept block whose
`Z K Z'` row-space contribution is diagonal and one-to-one with the
fitted rows; RoBMA then adds the BayesTools-prepared row multiplier to
the known-`V` extra variance. Other known-`R` blocks remain sampled or
are rejected by BayesTools validation. Slopes and other unsupported
variance forms remain sampled. Multiple declared terms with one grouping
level per estimate are rejected, regardless of their coefficient basis
or known group covariance.

The implementation intentionally omits the `weights` and `cluster`
arguments from
[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html).
Use `random` for multilevel structures. Estimate-unit
[`log_lik()`](https://fbartos.github.io/RoBMA/reference/log_lik.md),
[`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
[`add_waic()`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`rstudent()`](https://rdrr.io/r/stats/influence.measures.html) are
available for known-`V` models. Bridge-sampling marginal likelihoods via
[`add_marglik()`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md)
are available for fitted single-model known-`V` objects.
[`rstandard()`](https://rdrr.io/r/stats/influence.measures.html) is
available for known-`V` models at supported estimate and marginal
conditioning depths.

## Diagnostic target semantics

The main `brma.mv()` post-fit methods intentionally use different
targets:

- `log_lik(unit = "estimate")`,
  [`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md)/[`loo()`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
  [`add_waic()`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md)/[`waic()`](https://fbartos.github.io/RoBMA/reference/waic.brma.md),
  and [`rstudent()`](https://rdrr.io/r/stats/influence.measures.html)
  use estimate-unit log-score targets. When known `V` induces
  dependence, these contributions are Schur-complement conditionals
  \\p(y_i \mid y\_{-i}, \theta)\\ within dependency blocks. They are
  existing-estimate diagnostics, not independent new-estimate prediction
  targets. Every local Gaussian random effect is integrated, so the
  target is the same whether the fitting parameterization sampled or
  marginalized the effect. When known random-effect covariance `R` is
  supplied, the known `R` matrix and the fitted random-effect covariance
  enter the marginal \\ZGZ'\\ covariance of the Schur conditional
  through BayesTools' compiled random-effect metadata.

- [`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html) and
  [`covratio()`](https://fbartos.github.io/RoBMA/reference/covratio.brma.md)
  use estimate-unit PSIS weights. With correlated known `V`, their
  interpretation is parameter influence under conditional estimate
  deletion, not independent-study deletion.

- [`add_marglik()`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md)/[`bridge_sampler()`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html)
  uses the full fitted likelihood corresponding to the selected
  known-`V` backend. It does not use the estimate-wise conditional
  target from LOO/WAIC. During bridge evaluation, sampled Gaussian
  location random effects are integrated exactly as \\ZGZ'\\ while
  retaining every covariance parameter and prior; fitted diagonal
  marginalized blocks remain in the row variance. Selection likelihoods
  instead preserve the normalization boundary declared by
  [`selection_model()`](https://fbartos.github.io/RoBMA/reference/selection_model.md):
  retained contextual sources remain outside the normalizer and
  integrated sources enter the conditional covariance.

- [`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html),
  marginal
  [`rstandard()`](https://rdrr.io/r/stats/influence.measures.html), and
  [`vif()`](https://fbartos.github.io/RoBMA/reference/vif.md) use a
  marginal GLS covariance target based on `V + ZGZ'`, where formula
  random effects are marginalized through the BayesTools covariance
  metadata.

- [`dffits()`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md)
  and
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  use fixed-location fitted-value influence targets. Cluster/block LOO
  is deferred for `brma.mv()` because its deletion target is a joint
  dependency block, not an independent row.

Without `random`, the model estimates no heterogeneity. Explicit
random-formula models report component-specific standard deviation,
correlation and allocation summaries. A single homogeneous random
component uses the RoBMA names `tau` and `tau2`.
`summary_heterogeneity(component = ...)` describes component
heterogeneity over the observed design, whereas
`pooled_heterogeneity(component = ...)` evaluates absolute component
heterogeneity at the average expanded random design. Both report the
fitted heterogeneity scale without fixed known-`R` group covariance
multipliers; row-specific scale prediction and covariance consumers
retain the full known `R`. General known-`V` models do not report scalar
I2/H2.

## See also

[`random_effect_formula_tags`](https://fbartos.github.io/RoBMA/reference/random_effect_formula_tags.md),
[`random_effect_prior_specification`](https://fbartos.github.io/RoBMA/reference/random_effect_prior_specification.md),
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md),
[`BMA.mv()`](https://fbartos.github.io/RoBMA/reference/BMA.mv.md),
[`RoBMA.mv()`](https://fbartos.github.io/RoBMA/reference/RoBMA.mv.md),
[`bPET.mv()`](https://fbartos.github.io/RoBMA/reference/bPET.mv.md),
[`bPEESE.mv()`](https://fbartos.github.io/RoBMA/reference/bPEESE.mv.md),
[`bselmodel.mv()`](https://fbartos.github.io/RoBMA/reference/bselmodel.mv.md),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
fit <- brma.mv(
  yi                        = c(0.10, 0.20),
  V                         = V,
  measure                   = "GEN",
  prior_unit_information_sd = 1,
  seed                      = 1,
  silent                    = TRUE
)
summary(fit)
} # }
```
