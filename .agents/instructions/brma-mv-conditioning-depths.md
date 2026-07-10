# brma.mv Conditioning Depths and Predictive Targets

Use this note when changing `log_lik()`, `add_loo()`, `add_waic()`,
`rstudent()`, `add_marglik()`, or known-`V`/known-`R` random-effect
diagnostics. It records the current design decisions so the same target
semantics do not need to be reconstructed from scratch.

## Core Vocabulary

- `unit` is the deletion/output unit requested by the user. Current values are
  `"estimate"` and `"cluster"`.
- `conditioning_depth` records which fitted latent effects are conditioned on
  by the target. Current values are `"marginal"`, `"cluster"`, and `"estimate"`.
- LOO/WAIC expose `unit`; they derive `conditioning_depth` with
  `.loo_conditioning_depth_from_unit()`.
- `log_lik()`, LOO, WAIC, and LOO-PIT residuals use predictive log-score
  targets. They are not marginal covariance diagnostics.
- `hatvalues()`, marginal `rstandard()`, Pearson residual SEs, and `vif()` use
  marginal GLS covariance targets.
- same-data `predict.brma.mv(type = "response")` also uses marginal component
  generation, so response predictions are comparable to `brma()` predictions
  rather than conditional fitted-random-effect BLUPs.
- `add_marglik()`/`bridge_sampler()` use the full joint fitted likelihood plus
  the fitted model prior density. They are not pointwise predictive scores.

The main target metadata lives in `attr(x, "RoBMA_target")`. For `brma.mv()`
diagnostic documentation, keep `.brma_mv_diagnostic_target_table()` in
`R/brma-mv-targets.R` synchronized with any semantic change.

## Current Predictive Targets

| Model/path | User target | Current interpretation |
| --- | --- | --- |
| `brma()`/`brma.norm()`, no cluster | `unit = "estimate"` | One scalar contribution per estimate. Normal models integrate estimate-level heterogeneity through `tau^2 + se_i^2`. |
| `brma()` 3-level, estimate unit | `unit = "estimate"` | One contribution per estimate, conditional on fitted cluster effects and marginal over estimate-level heterogeneity. |
| `brma()` 3-level, cluster unit | `unit = "cluster"` | One joint held-out-cluster contribution. The held-out cluster effect is integrated by analytic covariance, quadrature, or GLMM integration. |
| `brma.mv()`, known `V`, no random formula | `unit = "estimate"` | One contribution per estimate. Correlated known `V` uses Schur conditionals `p(y_i | y_-i, theta)` inside dependency blocks. |
| `brma.mv()`, known `V` plus sampled random formula | `unit = "estimate"` | One contribution per estimate, conditional on sampled fitted random effects. Correlated known `V` still uses Schur conditionals around the conditional mean. |
| `brma.mv()`, known `V` plus marginalized random block | `unit = "estimate"` | Same known-`V` estimate target, with marginalized estimate-level random variance entering as diagonal extra variance. |
| `brma.mv()`, known `V` plus sampled known `R` | `unit = "estimate"` | Same as sampled random formula: known-`R` sampled random effects are conditioned on. `R` shapes their posterior/prior but is not added as pointwise `ZGZ'`. |
| `brma.mv()`, known `V` plus marginalized known `R` | `unit = "estimate"` | Same known-`V` estimate target, with BayesTools-prepared known-`R` row multipliers entering as diagonal extra variance. |
| `brma.mv()` random formula, no known `V` | `unit = "estimate"` | Deferred for post-fit log-likelihood/LOO/WAIC. Availability is guarded in `.check_log_lik_target_available()`. |
| `brma.mv()` known `V` | `unit = "cluster"` | Deferred. There is no implemented joint dependency-block LOO target yet. |

Important caveat: `conditioning_depth = "cluster"` in fitted values means
"condition on cluster-level effects". `unit = "cluster"` for LOO/log_lik means
"score a held-out cluster as one joint unit". For normal multilevel LOO, the
held-out cluster effect is integrated, not conditioned on its fitted value.

## Known V

Known sampling covariance `V` is about dependence between observed estimates.
It is not a random-effect prior and it is not known group covariance.

Implementation paths:

- Availability: `.check_log_lik_target_available()` in `R/unit_level.R`.
- Estimate setup: `.estimate_likelihood_setup_from_parts()` in `R/log-lik.R`.
- Backend decision: `.known_v_estimate_target_uses_backend()` in
  `R/log-lik-known-v.R`.
- Schur conditionals: `.log_lik_known_v_estimate_target_from_setup()` and
  `.known_v_component_conditional_distribution()` in `R/log-lik-known-v.R`.
- LOO/WAIC metadata: `R/loo.R` and `.add_loo_target_metadata()` in
  `R/unit_level.R`.

For correlated known `V`, estimate-unit predictive scores are not independent
new-study scores. They are existing-estimate diagnostics:

```text
p(y_i | y_-i, theta)
```

within each known-`V` dependency block. Singleton blocks reduce to the scalar
normal likelihood.

## Random Formula Effects in brma.mv Known-V Scores

For `brma.mv()` random-formula models, estimate-depth LOO/log_lik/WAIC currently
exist only for known-`V` fits. The setup intentionally does this:

1. It evaluates fixed effects.
2. It adds sampled formula random effects to the mean with
   `.evaluate.brma.random_effects()` when `unit = "estimate"`.
3. It sets ordinary `tau` contributions to zero for sampled formula random
   effects, because those effects are already conditioned in the mean.
4. It adds only explicitly marginalized random-effect variance as diagonal
   extra variance through `.known_v_extra_variance_from_setup()`.
5. It evaluates the known-`V` scalar or Schur target around that conditional
   mean.

This is why adding a sampled random formula changes the target from
"marginal over that formula random effect" to "conditional on the fitted random
effect". That is not an accidental change from known `V`; it follows the
established estimate-depth convention.

## Known R

Known `R` is group-axis latent random-effect covariance/correlation. It is
separate from known sampling covariance `V`.

Conceptually:

```text
u = tau * z
z ~ N(0, K)
```

where `K` comes from `R` after `Rscale` handling in BayesTools. For marginal GLS
covariance consumers, the random-effect contribution is:

```text
tau^2 * Z K Z'
```

Current known-`R` rules:

- Public RoBMA API is `R`/`Rscale` on `brma.mv()`.
- BayesTools remains the implementation authority via group covariance metadata.
- Current support is random intercepts only.
- Supported known-`R` blocks may be auto-marginalized by
  `marginalize_estimate_level = TRUE` only when BayesTools validates a
  one-column random-intercept block with diagonal `Z K Z'` row-space
  contribution and one-to-one fitted row/group mapping.
- RoBMA obtains row multipliers with
  `BayesTools::random_effects_marginal_variance_factors(...,
  require_diagonal = TRUE, require_one_to_one = TRUE)` and uses them as
  `tau^2 * row_multiplier[i]` in known-`V` extra variance.
- Known-`R` blocks with off-diagonal row-space covariance, repeated group rows,
  random slopes, multi-column blocks, or row-indexed external SD sources are
  sampled or rejected by BayesTools/RoBMA validation; RoBMA does not implement a
  dense known-`R` likelihood.
- New prediction levels for known-`R` blocks are rejected.
- In estimate-unit `log_lik()`, LOO, WAIC, and LOO-PIT residuals, sampled
  known-`R` random effects are conditioned on; marginalized known-`R` blocks
  enter as diagonal extra variance.
- Known `R` is not split into pointwise likelihood terms and is not added again
  as pointwise marginal `ZGZ'` covariance for estimate-depth LOO/WAIC.

This is deliberate. Adding sampled `R` as pointwise likelihood would mix prior
density with observation likelihood and would double count sampled random
effects in the estimate-depth target. The supported marginalized known-`R` path
is different: the latent random effects are not sampled, so their validated
diagonal row-space contribution belongs in the observation variance.

Known-`R` target metadata is attached by
`.known_r_log_lik_target_metadata()` in `R/log-lik.R`:

- `known_r`
- `known_r_blocks`
- `known_r_semantics`

These metadata are informational. They should not by themselves make
`loo_compare()` reject two models. Different known-`R` structures are model
differences, not different deletion targets, as long as data, unit,
conditioning depth, and likelihood target match.

## Marginal Covariance Diagnostics

The following methods use a marginal GLS covariance target rather than a
pointwise predictive log-score:

- `hatvalues()`
- `rstandard(conditioning_depth = "marginal")`
- Pearson residual SEs with marginal conditioning
- `vif()`
- `predict(..., type = "response")` generative distribution for known-`V`
  `brma.mv()` paths

These paths use known `V` plus random-effect covariance:

```text
M = V + Z G Z'
```

For known `R`, BayesTools contributes the group-axis covariance inside
`ZGZ'`. RoBMA should keep delegating that construction to
`BayesTools::random_effects_marginal_vcov()` instead of duplicating matrix
scaling, ordering, or definiteness logic.

For same-data `brma.mv()` response predictions, the posterior mean is evaluated
from fixed/location terms only. Sampling noise is generated from known `V`, and
each applicable random-effect component is generated marginally. Sampled fitted
random effects are not added as conditional BLUPs. This intentionally mirrors
the univariate `brma()` response target, where between-study heterogeneity is
represented through posterior predictive generation rather than point BLUPs.
`type = "estimate"` remains the fitted latent-effect target for `newdata = NULL`
and conditions on sampled random effects when they are part of the fitted model.

Explicit prediction rows represent new true effects. Random-formula components
are generated marginally for each row; deliberately repeated ordinary grouping
IDs share a generated effect. Known-`R` rows must use fitted levels because an
`R_new` interface is deferred. The former unreleased logical aggregate
prediction mode does not exist. `pooled_effect()` averages fitted-design fixed
location draws directly and is not a prediction target.

Key files:

- `R/known-v-gls.R`
- `R/hat_matrix.R`
- `R/residuals.R`
- `R/vif.R`
- `R/predict-brma-mv.R`
- `R/heterogeneity-mv.R`

Location-scale pooled heterogeneity uses posterior root-mean-square
heterogeneity across fitted rows. The univariate and multivariate wrappers
therefore both average variances first and take the square root afterward.

## Marginal Likelihood

`add_marglik()` and `bridge_sampler()` use the full joint fitted model target.
They do not use the estimate-wise LOO/WAIC target.

For known `V`:

- the full fitted likelihood uses the selected known-`V` backend;
- this is not Schur estimate deletion.

For known `R`:

- the likelihood is conditional on latent random effects where the fitted JAGS
  model is conditional;
- the known group covariance enters through the latent random-effect prior
  density handled by BayesTools/JAGS bridge machinery;
- supported marginalized known-`R` blocks are not latent `dmnorm` blocks and
  instead enter through the same diagonal known-`V` extra variance used by JAGS;
- it is valid that `R` is not a pointwise likelihood term.

Key files:

- `R/marglik.R`
- `R/brma-mv-targets.R`
- BayesTools bridge/JAGS formula machinery

## LOO, WAIC, and log_lik Consistency

LOO and WAIC reuse the same log-likelihood matrix as `log_lik()` for the same
unit. Therefore, if a semantic change is made to the estimate log-likelihood,
it also changes LOO, WAIC, LOO-PIT residuals, and influence diagnostics that
use PSIS weights.

`log_lik()` returns pointwise posterior draws, not a scalar maximized
log-likelihood. For correlated known `V`, summing its Schur-conditional columns
produces a composite score rather than the full joint fitted likelihood.
`logLik.brma()` is intentionally absent, and `AIC()`/`BIC()` are unsupported.

At estimate depth, deletion means holding out one estimate within the retained
grouping and known-`V` dependency structure. It does not mean holding out a
whole new group. Sampled, analytically marginalized, and mixed random-effect
representations remain comparable when their data hash, deletion unit,
conditioning depth, and likelihood target agree. `random_effects` and
`estimate_level_random` record draw provenance and are not comparison keys.

This policy does not imply equal numerical reliability. PSIS reweights draws
toward the leave-one-estimate-out posterior, but sampled local effects can
produce high Pareto-k values. WAIC can be especially sensitive to whether a
latent effect is sampled or marginalized. Check diagnostics and compare matched
sampled/marginalized representations when changing these paths.

Current consistency checks:

- `.check_loo_compare_targets()` rejects mismatched data hashes, units,
  conditioning depths, or likelihood target labels.
- The four-field compatibility key is applied separately within LOO and WAIC;
  LOO and WAIC objects are not mixed in one comparison table.
- `data_hash` describes the observed outcome target, including known `V`.
- Cached LOO/WAIC extraction recomputes this four-field target key and rejects
  stale results after outcome, known-`V`, unit, or target changes.
- Known `R` metadata is preserved for transparency but does not define a
  separate deletion target.
- WAIC target metadata is copied from the same log-likelihood target used by
  LOO.

## If Adding a Dense Fully Marginalized Known-R Predictive Score

A dense fully marginalized known-`R` predictive score for sampled or
off-diagonal known-`R` blocks would be a new target, not a silent replacement
for estimate-depth LOO. It would integrate sampled known-`R` random effects and
score with a covariance containing `ZGZ'`. This is broader than the current MVP,
which only adds BayesTools-validated diagonal `tau^2 * row_multiplier[i]`
variance for blocks that are not sampled.

To add it safely:

1. Add an explicit public or internal target selector.
2. Add new `RoBMA_target` labels so comparisons cannot mix conditional and
   marginalized known-`R` scores accidentally.
3. Define whether the target is estimate-wise Schur conditional, whole-block
   joint, or another deletion unit.
4. Revisit PSIS reliability because high-dimensional latent random-effect
   marginalization changes the importance-sampling problem.
5. Add direct tests against analytic or BayesTools covariance oracles.

Do not implement this by changing the current `unit = "estimate"` target in
place. The existing target is intentionally conditional on sampled fitted
random effects.

## Tests That Guard This Design

Current relevant tests include:

- `tests/testthat/test-00-input-data-mv.R`
  - diagnostic target registry
  - known-`V` Schur target
  - target metadata availability
  - known-`R` row multipliers in known-`V` syntax/data
  - log_lik/bridge extra variance from `tau^2 * row_multiplier`
  - backend restrictions for row-varying known-`R` variance
- `tests/testthat/test-00-brma-mv-known-r.R`
  - known-`R` metadata
  - JAGS/formula design preservation
  - `Rscale = "none"` marginal covariance semantics
  - supported known-`R` auto-marginalization and sampled fallbacks
- `tests/testthat/test-03-loo.R`
  - direct known-`R` `log_lik()`, LOO, and WAIC metadata
  - known-`V` estimate-unit target metadata
  - unified no-random, sampled, marginalized, and mixed comparison policy
  - matched sampled/marginalized reliability and stale-cache checks
- `tests/testthat/test-03-bridgesampling.R`
  - marginal likelihood target metadata and bridge availability
- `tests/testthat/test-02-vif.R`, `test-02-residuals.R`,
  `test-02-hatvalues.R`, and `test-02-predict-mv.R`
  - marginal covariance consumers and known-`R`/known-`V` behavior

When changing these semantics, run the affected `test-00-*` tests, regenerate
cached `test-01-*` fits if the fitted object shape changes, and rerun the
post-fit tests that consume LOO/WAIC, marginal covariance, prediction, and
bridge sampling.
