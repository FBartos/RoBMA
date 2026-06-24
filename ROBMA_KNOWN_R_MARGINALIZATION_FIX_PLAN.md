# RoBMA Plan: Marginalized Known-R with Known-V

## Status

BayesTools now owns the random-effect design support needed for known group
covariance in marginalized one-column random-intercept blocks. RoBMA should
consume that support instead of carrying its own covariance scaling, level
matching, or formula-design logic.

Required BayesTools-side APIs:

- `random_group_covariance()`
- `random_effects_formula(..., group_covariance = ...)`
- `random_effects_compile(..., marginalized = ...)`
- `random_effects_marginal_vcov()`
- `random_effects_marginal_variance_factors()`

The RoBMA side should stay responsible for the metafor-facing interface:

- `R`
- `Rscale`
- `known_v_parameterization`
- `known_v_residual_fraction`
- `marginalize_estimate_level`

Do not expose RoBMA/metafor names in BayesTools-facing metadata beyond the
existing adapter that turns `R`/`Rscale` into BayesTools
`random_group_covariance()` objects.

## Goal

Allow `brma.mv()` known-V normal models to marginalize supported known-R
random-intercept blocks when `marginalize_estimate_level = TRUE`.

The supported MVP is deliberately narrow:

- one random-effect column, i.e. random intercept only
- group levels match the known-R kernel through BayesTools design metadata
- row-space contribution is diagonal
- each observed row maps one-to-one to a random-effect group level
- no row-indexed external SD source
- no new levels in prediction/reconstruction

For this case, the known-V likelihood can consume the marginalized random
effect as a diagonal extra variance contribution:

```text
tau^2 * row_multiplier[i]
```

where `row_multiplier` is provided by BayesTools from the prepared known-R
metadata. RoBMA must not silently convert `Rscale = "none"` to a correlation
matrix or otherwise rescale the kernel itself.

## Semantics

Known-V and known-R remain distinct:

- Known-V is observed estimate covariance.
- Known-R is group-axis covariance for a latent random-effect block.

This feature is not a general dense likelihood implementation for
`V + tau^2 * Z K Z'`. RoBMA may use the diagonal extra-variance path only when
BayesTools confirms that diagonal consumption is valid. Use
`require_diagonal = TRUE` and `require_one_to_one = TRUE` when asking
BayesTools for row-level marginal variance factors.

Sampled known-R behavior must remain unchanged:

- latent random effects are sampled when the block is not marginalized
- estimate-depth log-likelihood, LOO, WAIC, CDF, residual, and prediction
targets stay conditional on sampled known-R effects
- known-R is not treated as a pointwise likelihood term for sampled blocks

Marginalized known-R behavior should follow the existing marginalized random
effect semantics in RoBMA's known-V backend:

- the extra diagonal variance is added to the known-V estimate likelihood
- bridge-sampling likelihood evaluation uses the same variance contribution
- post-fit log-likelihood and predictive summaries use the same contribution
- GLS and marginal covariance diagnostics may continue to use the full
  BayesTools marginal covariance where that is already the diagnostic target

Do not add a nugget to the known-V matrix to make an unsupported design pass.
Unsupported dense row-space covariance should error clearly.

## Non-Goals

Do not implement any of the following in this pass:

- dense known-R likelihoods for off-diagonal `tau^2 * Z K Z'`
- known-R random slopes
- multi-column known-R random-effect blocks
- row-indexed external SD sources for known-R blocks
- new-level prediction for known-R blocks
- GLMM known-R marginalization
- selection-model known-R marginalization
- cluster/unit-depth predictive scoring changes
- RoBMA-side copies of BayesTools kernel scaling or level matching
- BayesTools API changes or BayesTools-facing `R`/`Rscale` names

## Implementation Areas

### 1. Random-Effect Compile Selection

Update `.prepare_brma_mv_random_effects_compile()` in
`R/random-effects-compile.R` so supported known-R random-intercept blocks may
be compiled as marginalized when `marginalize_estimate_level = TRUE`.

The selection rule should be conservative:

- keep `marginalize_estimate_level = FALSE` sampled behavior unchanged
- keep unsupported known-R designs sampled only if existing semantics require
  sampled behavior and no known-V diagonal consumption is requested
- error when the user explicitly requests a marginalized known-R path that
  cannot be consumed safely by RoBMA's known-V backend
- rely on BayesTools validation for random slopes, multi-column blocks,
  mismatched levels, and unsupported compile modes

Prefer a small helper that decides whether a term is eligible for RoBMA's
known-V diagonal marginalization path. Keep it separate from BayesTools
formula-design internals.

### 2. Row Multiplier Metadata

After the BayesTools formula design is available, call:

```r
BayesTools::random_effects_marginal_variance_factors(
  formula_design,
  require_diagonal = TRUE,
  require_one_to_one = TRUE
)
```

Use the returned block metadata to store row-aligned variance multipliers in
RoBMA's existing marginalized-random metadata. Do not recompute:

- kernel scaling
- group-level ordering
- row-to-group matching
- diagonal validity

The stored metadata should be sufficient for:

- JAGS data preparation
- JAGS variance expressions
- post-fit extra variance evaluation
- bridge likelihood evaluation
- same-data prediction/reconstruction where already supported

Keep naming RoBMA-native in RoBMA internals, but make it obvious that the value
is a multiplier, not a variance parameter. Suggested names include
`row_multiplier`, `marginal_variance_multiplier`, or a similar existing local
pattern.

### 3. JAGS Syntax and Data

Extend the existing marginalized random variance expression path so known-R
marginalized blocks contribute:

```text
pow(sd_parameter, 2) * row_multiplier[i]
```

or the equivalent expression used by the current known-V parameterization.

This must work for the existing normal known-V backends:

- scalar/latent known-V path
- whitened known-V path
- block-MVN known-V path

For marginalized known-R terms, do not generate:

- latent group-level `dmnorm` syntax
- sampled random-effect coefficient matrices
- latent random-effect monitors
- coefficient monitors for the marginalized block

Sampled known-R terms should continue to generate the existing latent syntax,
data, and monitor names.

### 4. Post-Fit, Log-Likelihood, and Bridge Evaluation

Update the downstream variance evaluation helpers so marginalized known-R
blocks use the stored BayesTools row multipliers:

- `.evaluate_marginalized_random_variance()`
- `.known_v_extra_variance_from_setup()`
- `.marglik_known_v_extra_variance()`
- any related helper that currently assumes a unit row multiplier

For a posterior draw, each supported marginalized known-R block should add:

```text
tau_draw^2 * row_multiplier
```

Bridge-sampling likelihood evaluation must use the same expression as fitted
object log-likelihood evaluation. Avoid two independent implementations of the
same calculation.

BayesTools bridge rebuild validation should remain active. If the known-R
metadata changes between fit and rebuild, validation must fail rather than
silently evaluating the wrong model.

### 5. Prediction and Reconstruction

Preserve existing known-R new-level restrictions. For the MVP, prediction with
new levels for a marginalized known-R block should error clearly.

Same-data prediction or reconstruction may use stored row multipliers when the
existing prediction path already supports marginalized random effects. Do not
invent a broader newdata matching system in this pass.

### 6. Diagnostics and Marginal Covariance

Keep the existing diagnostic distinction:

- likelihood extra variance uses the diagonal row multiplier path only when
  BayesTools confirms diagonal one-to-one consumption
- marginal covariance diagnostics may use
  `BayesTools::random_effects_marginal_vcov()` to obtain the full
  `tau^2 * Z K Z'` contribution where that is already the intended diagnostic
  target

Do not replace full diagnostic covariance with diagonal-only approximations.

### 7. Documentation and NEWS

Update `brma.mv()` documentation and NEWS to state the new behavior precisely:

- `R`/`Rscale` still define known group covariance for random-intercept blocks
- sampled known-R behavior remains supported
- when `marginalize_estimate_level = TRUE`, RoBMA can marginalize known-R
  random-intercept blocks in known-V normal models only when the BayesTools
  diagonal one-to-one checks pass
- unsupported known-R designs remain unsupported for marginalization

Remove or revise documentation that says known-R blocks are never
auto-marginalized.

If RoBMA needs a development dependency on the BayesTools version containing
`random_effects_marginal_variance_factors()`, update the minimum version only
after the target BayesTools release version is chosen.

## Testing Plan

Add focused tests before broad release-profile runs.

### Compile and Metadata Tests

In `tests/testthat/test-00-brma-mv-known-r.R` or a nearby existing file:

- known-R plus known-V with `marginalize_estimate_level = TRUE` compiles the
  supported random-intercept block as marginalized
- the stored RoBMA marginalized-random metadata includes the BayesTools row
  multiplier in row order
- sampled behavior remains unchanged when `marginalize_estimate_level = FALSE`
- sampled behavior remains unchanged for paths that are intentionally not
  eligible for the diagonal known-V marginalization path
- no latent known-R random-effect nodes or monitors are generated for the
  marginalized block
- latent syntax and monitor names remain unchanged for sampled known-R blocks

### Scaling Tests

Use small deterministic kernels to verify row multipliers for:

- `Rscale = "none"`
- `Rscale = "cor"`
- `Rscale = "cor0"`
- `Rscale = "cov0"`

These tests should assert RoBMA consumes BayesTools-prepared multipliers and
does not silently normalize `"none"` to a correlation matrix.

### Error Tests

Assert clear errors for:

- off-diagonal row-space covariance when diagonal known-V consumption would be
  required
- repeated group levels when one-to-one mapping is required
- random slopes
- multi-column known-R random effects
- row-indexed external SD sources
- prediction with new known-R levels
- changed known-R metadata during bridge rebuild validation

### JAGS and Likelihood Tests

Add or update tests around known-V syntax/data:

- scalar/latent known-V uses the known-R row multiplier in the marginalized
  extra variance expression
- whitened known-V uses the same multiplier semantics
- block-MVN known-V uses the same multiplier semantics
- `.evaluate_marginalized_random_variance()` equals
  `tau^2 * row_multiplier` for posterior draws
- bridge/marginal-likelihood known-V extra variance equals the fitted-object
  log-likelihood extra variance

### Diagnostics and Prediction Tests

Where affected, verify:

- known-V GLS/marginal covariance diagnostics still use the full
  BayesTools marginal covariance contribution
- same-data prediction uses the stored marginalized variance contribution
  where the existing path supports it
- new-level prediction errors for marginalized known-R blocks

## Suggested Test Commands

Start with focused tests:

```powershell
Rscript -e "Sys.setenv(AGENT='1'); devtools::test(filter='00-brma-mv-known-r|00-input-data-mv', reporter=testthat::LlmReporter$new())"
```

If likelihood, bridge, or fitted-object behavior changes, run the fitting and
downstream tests according to RoBMA's test policy:

```powershell
Rscript -e "Sys.setenv(AGENT='1'); devtools::test(filter='01-', reporter=testthat::LlmReporter$new())"
Rscript -e "Sys.setenv(AGENT='1'); devtools::test(filter='02-|03-', reporter=testthat::LlmReporter$new())"
```

When stable, run the release profile if feasible:

```powershell
Rscript tools/test-profile.R release
```

If the full profile cannot run for environmental reasons, report the exact
blocker and the narrower profiles that passed.

## Review Checklist

Before considering the RoBMA integration complete, perform several senior
developer review rounds. Check specifically for:

- accidental RoBMA-side duplication of BayesTools covariance scaling
- accidental conversion of `Rscale = "none"` to a correlation matrix
- inconsistent variance calculations between JAGS, logLik, bridge, and
  prediction
- conditional sampled known-R semantics being changed unintentionally
- dense off-diagonal row covariance being silently dropped
- over-broad newdata support
- confusing names that blur known-V and known-R concepts
- brittle dependence on BayesTools internal object structure
- missing documentation for the narrow diagonal one-to-one limitation

After each review round, clean issues found and rerun the relevant focused
tests. Do not stop at the first passing test run if the implementation can be
simplified or made more robust.
