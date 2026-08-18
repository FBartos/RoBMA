# Multivariate Predictive Targets

Use this guide when changing `brma.mv()`, `log_lik()`, LOO/WAIC,
residuals, marginal covariance diagnostics, prediction, or marginal likelihood.

## Keep Targets Distinct

- `unit` is the deletion/output unit: `"estimate"` or `"cluster"`.
- `conditioning_depth` records which fitted latent effects are conditioned on:
  `"marginal"`, `"cluster"`, or `"estimate"`.
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

At estimate depth, sampled fitted random effects are conditioned on in the
mean. Only explicitly marginalized, validated diagonal random-effect variance
is added to the observation variance.

Cluster-unit known-`V` scoring and estimate-unit scoring for random-formula
models without known `V` are currently unavailable. Do not silently substitute
a different target.

## Known Random-Effect Correlation R

Known `R` describes group-axis latent random-effect covariance, not sampling
covariance. BayesTools owns its scaling and metadata.

- Sampled known-`R` effects are conditioned on in estimate-unit predictive
  scores.
- Supported marginalized one-to-one random-intercept blocks contribute their
  validated diagonal variance.
- Never add sampled `R` again as pointwise `ZGZ'`; that mixes the prior with
  observation likelihood and double counts the effect.
- Delegate marginal random-effect covariance to
  `BayesTools::random_effects_marginal_vcov()` and validated marginal variance
  factors to BayesTools.

Dense or off-diagonal marginalized known-`R` predictive scoring would be a new
explicit target. It must not replace the current estimate target in place.

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
mismatched data, unit, conditioning depth, or likelihood target. Known-`R`
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
