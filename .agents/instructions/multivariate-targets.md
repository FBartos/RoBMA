# Multivariate Predictive Targets

Use this guide when changing `brma.mv()`, `log_lik()`, LOO/WAIC,
residuals, marginal covariance diagnostics, prediction, or marginal likelihood.

## Keep Targets Distinct

- `unit` is the deletion/output unit: `"estimate"` or `"cluster"`.
- `conditioning_depth` records which fitted latent effects are retained by
  prediction and residual targets: `"marginal"`, `"cluster"`, or
  `"estimate"`. LOO target metadata instead records the deletion `unit` and
  `retained_context`; do not equate these axes.
- `log_lik()`, LOO, WAIC, and LOO-PIT use predictive log-score targets.
- `hatvalues()`, marginal/Pearson residual scaling, and `vif()` use marginal
  GLS covariance.
- `add_marglik()` and `bridge_sampler()` use the full joint fitted
  likelihood and prior, not a pointwise predictive score.

Do not reuse one target's covariance or likelihood merely because dimensions
match.

## Known Sampling Covariance V

For correlated known `V`, estimate-unit log scores are Schur conditionals
`p(y_i | y_-i, theta)` within dependency blocks. Their column sum is a
composite score, not the full joint likelihood.

For Gaussian estimate-unit scores, integrate every local Gaussian random
effect through BayesTools' compiled marginal covariance plan. The target is
`p(y_i | y_-i, theta)`: retained rows from the same dependency block remain
available after deleting estimate `i`. This is the same target whether the
fitting parameterization sampled or marginalized a local effect. Never infer
dependencies from posterior draws or reconstruct a random structure in RoBMA.

The specialized `brma(..., cluster = ...)` interface uses the same Gaussian
covariance scorer for estimate deletion. Cluster-unit scoring deletes the
whole cluster jointly and retains its distinct new-cluster target. Non-Gaussian,
selection-model, and weighted likelihoods retain their supported conditional
representations; do not silently approximate their local-effect integrals.

## Known Random-Effect Correlation R

Known `R` describes group-axis latent random-effect covariance, not sampling
covariance. BayesTools owns its scaling and metadata.

- Sampled and marginalized known-`R` effects enter Gaussian estimate scores
  through the same BayesTools metadata-defined marginal `ZGZ'` covariance.
- Delegate covariance factor plans and states to BayesTools. RoBMA may combine
  them with known sampling covariance and validated row variance, but must not
  recompile or special-case individual random structures.
- The full joint fitted likelihood used by marginal likelihood remains distinct
  from the estimate-wise Schur score.

## Diagnostics and Prediction

Marginal covariance consumers use `M = V + ZGZ'`. Same-data
`predict.brma.mv()` follows the same two-axis contract as `brma()`:

- `type = "terms"` returns fixed location, `type = "estimate"` returns latent
  true-effect draws, and `type = "response"` adds sampling error;
- `conditioning_depth = "marginal"` draws all applicable random effects anew,
  irrespective of `newdata = NULL` versus an equivalent explicit design;
- `conditioning_depth = "estimate"` draws fitted latent effects from their
  conditional posterior, including conditional uncertainty rather than only
  Gaussian BLUP means;
- cluster depth is unavailable for arbitrary random formulas because their
  hierarchy does not identify one canonical cluster level.

`newdata` selects design and identity, not conditioning. Non-marginal explicit
rows must be rejected unless fitted identities can be validated. Marginal
matching labels preserve joint new-draw dependence and never reuse fitted BLUPs.
Response prediction preserves full `V`/`V_new` and joint `ZGZ'` dependence.
Conditional means remain available through `blup()` and `fitted()`.

Keep target metadata in `attr(x, "RoBMA_target")`. LOO comparisons must reject
mismatched data, unit, retained context, or likelihood target. Known-`R`
metadata is informative and is not itself a comparison key.

## Relevant Files

- Target availability and metadata: `R/unit_level.R`,
  `R/brma-mv-targets.R`.
- Known-`V` scores: `R/log-lik.R`, `R/log-lik-known-v.R`.
- Marginal covariance: `R/known-v-gls.R`,
  `R/random-marginal-vcov.R`.
- Prediction: `R/predict-brma-mv.R`.
- Marginal likelihood: `R/marglik.R`.

Any semantic change must update the target registry, metadata, documentation,
and the focused consumers of the same target.
