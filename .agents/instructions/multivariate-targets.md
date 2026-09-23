# Predictive and Multivariate Targets

Use this guide when changing prediction, `brma.mv()`, `log_lik()`, LOO/WAIC,
residuals, marginal covariance diagnostics, or marginal likelihood.

## Keep Targets Distinct

- `unit` is the deletion/output unit: `"estimate"` or `"cluster"`.
- `conditioning_depth` records which fitted latent effects are retained by
  prediction and residual targets: `"marginal"`, `"cluster"`, or
  `"estimate"`. LOO target metadata instead records the deletion `unit` and
  `retained_context`; do not equate these axes.
- `log_lik()`, LOO, WAIC, and LOO-PIT use predictive log-score targets.
- `hatvalues()` and `vif()` use the fitted GLS projection; Pearson and
  standardized residual scaling use the marginal outcome covariance.
- `add_marglik()` and `bridge_sampler()` use the full joint fitted
  likelihood and prior, not a pointwise predictive score.

Do not reuse one target's covariance or likelihood merely because dimensions
match.

An integer likelihood weight `w` repeats the conditional likelihood contribution
`w` times with the declared latent-effect and grouping structure held fixed.
Fractional likelihood powers remain supported. This fitting interpretation does
not choose between per-copy and whole-row deletion or define fractional-weight
diagnostic conventions.

Likelihood weights affect fitting and PSIS deletion, not the observation law
used to scale diagnostics. Pearson and LOO-PIT/rstudent use the original outcome
variance/CDF. Internally standardized residuals use the fitted weighted
projection with the original covariance in `(I-H) M (I-H)'`; the ordinary GLS
subtraction identity applies only without extra likelihood weights. Keep
Bayesian posterior averaging and LOO-PIT distinct from metafor's refit z-score.

## Known Sampling Covariance V

Covariance classification, sampling support, and numerical representations are
distinct. Preserve the accepted sampling factor and PSD/PD status when a
positive direction cannot be represented by a Cholesky pivot or squared
singular value. Operations requiring that unavailable representation fail with
their target context; they must not treat underflow as a structural zero, add
jitter, or silently use a lower-rank inverse. Genuine PSD targets retain their
supported spectral paths, and sampling may use the accepted factor directly.

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

The specialized `brma(..., cluster = ...)` interface remains supported. Its
two-level allocation permits within- and between-cluster I2 even though its
covariance can also be represented by `brma.mv()`. Use
`y = X beta + u_cluster + u_estimate + epsilon` as its canonical notation;
random-formula models replace the named effects with their fitted Gaussian
blocks.

Marginal covariance consumers use `M = V + ZGZ'`. Same-data
`predict.brma.mv()` follows the same two-axis contract as `brma()`:

- `type = "terms"` returns fixed location, `type = "estimate"` returns latent
  true-effect draws, and `type = "response"` adds sampling error;
- `conditioning_depth = "marginal"` draws all applicable random effects anew,
  irrespective of `newdata = NULL` versus an equivalent explicit design;
- `conditioning_depth = "estimate"` draws fitted latent effects from their
  conditional posterior, including conditional uncertainty rather than only
  Gaussian BLUP means;
- `conditioning_depth = "cluster"` retains the specialized multilevel model's
  fitted cluster effect and predicts a new estimate within it;
- cluster depth is unavailable for arbitrary random formulas because their
  hierarchy does not identify one canonical cluster level.

Marginal is the canonical prediction default. Response prediction adds only
the sampling layer to the corresponding latent target. For GLMM responses,
estimate depth uses fitted posterior `pi`/`phi`; marginal and cluster depths
draw a new nuisance rate from its prior and are partly prior predictive.
Document this wherever those GLMM targets are exposed.

`newdata` selects design and identity, not conditioning. Non-marginal explicit
rows must be rejected unless fitted identities can be validated. Marginal
matching labels preserve joint new-draw dependence and never reuse fitted BLUPs.
Response prediction preserves full `V`/`V_new` and joint `ZGZ'` dependence.
Conditional means remain available through `blup()` and
`fitted(..., conditioning_depth = "estimate")`.

Normal-model funnel/sampling bands, regression/forest prediction intervals,
default z-plots, and default residual diagnostics use explicit marginal
targets. GLMM funnel and regression sampling bands are descriptive normal
effect-size approximations and must be labeled accordingly. LOO and LOO-PIT
keep deletion-conditioned predictive-score targets and must not be routed
through generic response prediction.

Keep target metadata in `attr(x, "RoBMA_target")`. LOO comparisons must reject
mismatched data, unit, retained context, or likelihood target. Known-`R`
metadata is informative and is not itself a comparison key.
Cluster-unit scores also retain the ordered row partition; matching cluster
counts or labels alone do not make pointwise scores comparable. Recompute old
cluster diagnostics that lack this partition metadata.

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
