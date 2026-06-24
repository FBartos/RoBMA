# Goal Prompt: RoBMA Marginalized Known-R with Known-V

```text
/goal Create and pursue a goal to fully implement the RoBMA marginalized known-R with known-V integration described in:

  C:\R-Packages\RoBMA\ROBMA_KNOWN_R_MARGINALIZATION_FIX_PLAN.md

Work in:

  C:\R-Packages\RoBMA

Implement the feature end to end, not just the narrow syntax path. Treat `R`/`Rscale` as RoBMA/metafor-facing adapters that are converted to BayesTools known group covariance metadata. Do not duplicate BayesTools kernel scaling, level matching, diagonal checks, or formula-design logic in RoBMA.

Before editing, read the plan, AGENTS.md, `.agents/instructions/engineer-mode.md`, `.agents/instructions/r-development.md`, `.agents/instructions/testing.md`, `.agents/instructions/brma-mv-conditioning-depths.md`, relevant `brma.mv` known-V/known-R/random-effects/JAGS/logLik/bridge/prediction/diagnostic code, and existing tests. Match RoBMA style: snake_case, `<-`, two-space indentation, BayesTools validators where used locally, concise `stop(..., call. = FALSE)` messages, roxygen for user-facing API docs, and scoped changes.

Implementation requirements:
- Allow supported known-R random-intercept blocks in `brma.mv()` known-V normal models to be marginalized when `marginalize_estimate_level = TRUE`.
- Keep the supported MVP narrow: one random-effect column, random intercept only, BayesTools-validated known-R metadata, diagonal row-space contribution, one-to-one row/group mapping, no row-indexed external SD sources, and no new-level prediction.
- Preserve sampled known-R behavior exactly when `marginalize_estimate_level = FALSE` or when a known-R block remains sampled: existing latent syntax, JAGS data, monitor names, conditional estimate-depth logLik/LOO/WAIC/CDF/residual/prediction semantics, and bridge behavior.
- Consume `BayesTools::random_effects_marginal_variance_factors(..., require_diagonal = TRUE, require_one_to_one = TRUE)` for marginalized known-R row multipliers after the BayesTools formula design is available.
- Do not recompute known-R kernel scaling, group-level ordering, row matching, diagonal validity, or one-to-one validity in RoBMA.
- Store BayesTools-prepared row multipliers in RoBMA marginalized-random metadata in a way usable by JAGS data, JAGS variance expressions, fitted-object logLik, bridge likelihoods, and same-data prediction/reconstruction where existing marginalized-random paths support them.
- Ensure marginalized known-R terms add `tau^2 * row_multiplier[i]` to the known-V extra variance path.
- Ensure scalar/latent, whitened, and block-MVN known-V parameterizations use the same multiplier semantics.
- Do not generate latent group-level `dmnorm` syntax, sampled coefficient matrices, latent random-effect monitors, or coefficient monitors for marginalized known-R terms.
- Keep sampled known-R terms on the existing sampled path.
- Update `.evaluate_marginalized_random_variance()`, known-V logLik setup, and marginal-likelihood/bridge helpers so fitted-object likelihood and bridge likelihood evaluate the same extra variance.
- Keep BayesTools bridge rebuild validation active so changed known-R metadata fails validation instead of silently rebuilding the wrong likelihood.
- Preserve the diagnostic distinction: likelihood extra variance may use the diagonal row multiplier only after BayesTools confirms it is valid, while GLS/marginal covariance diagnostics may continue using `BayesTools::random_effects_marginal_vcov()` for the full `tau^2 * Z K Z'` contribution where that is the intended diagnostic target.
- Preserve known-R new-level prediction restrictions; for this MVP, error clearly on new levels for marginalized known-R blocks.
- Update `brma.mv()` documentation and NEWS to replace any statement that known-R blocks are never auto-marginalized with the narrower supported behavior.
- If RoBMA must depend on a newer BayesTools release for `random_effects_marginal_variance_factors()`, update the minimum dependency only after checking the intended BayesTools release version.
- Do not implement dense known-R likelihoods, GLMM/selection-model known-R marginalization, cluster/unit-depth predictive scoring changes, or BayesTools-facing `R`/`Rscale` APIs.

Testing requirements:
- Add focused tests for known-R plus known-V with `marginalize_estimate_level = TRUE` compiling the supported random-intercept block as marginalized.
- Test that RoBMA stores and consumes BayesTools row multipliers in row order.
- Test that marginalized known-R emits no latent known-R random-effect nodes or monitors.
- Test that sampled known-R behavior remains unchanged when `marginalize_estimate_level = FALSE` or when sampled behavior is intended.
- Test row multipliers and downstream variance for `Rscale = "none"`, `"cor"`, `"cor0"`, and `"cov0"` without silently converting `"none"` to a correlation matrix.
- Test scalar/latent, whitened, and block-MVN known-V syntax/data or variance-expression behavior as applicable.
- Test `.evaluate_marginalized_random_variance()` and bridge/marginal-likelihood known-V extra variance equal `tau^2 * row_multiplier`.
- Test unsupported designs still error or remain sampled as appropriate: off-diagonal row-space covariance under diagonal consumption, repeated group levels under one-to-one consumption, random slopes, multi-column known-R blocks, row-indexed external SD sources, and prediction new levels.
- Test BayesTools bridge rebuild validation detects changed known-R metadata for marginalized terms.
- Test diagnostics/prediction paths touched by the implementation, including full marginal covariance diagnostics if changed.
- Run focused tests first, for example:
  `Rscript -e "Sys.setenv(AGENT='1'); devtools::test(filter='00-brma-mv-known-r|00-input-data-mv', reporter=testthat::LlmReporter$new())"`
- If likelihood, bridge, or fitted-object behavior changes, run fitting and downstream tests:
  `Rscript -e "Sys.setenv(AGENT='1'); devtools::test(filter='01-', reporter=testthat::LlmReporter$new())"`
  `Rscript -e "Sys.setenv(AGENT='1'); devtools::test(filter='02-|03-', reporter=testthat::LlmReporter$new())"`
- When stable, run `Rscript tools/test-profile.R release` if feasible.
- If any test profile cannot run for environmental reasons, report the blocker precisely.

Review requirements:
- After tests pass, perform at least three senior-developer review rounds.
- In each round, inspect correctness, API alignment, known-V versus known-R semantics, sampled versus marginalized predictive targets, duplicated BayesTools logic, overengineering, brittle assumptions, error messages, documentation, and under-tested edge cases.
- Clean up issues found in each round and rerun the relevant tests.
- Do not stop after the first passing test run if the implementation can still be simplified or made more robust.

Git permission for this task:
- You may run non-mutating inspection commands only if explicitly allowed by the current repository instructions and current user turn.
- Do not run mutating git commands, do not commit, and do not revert unrelated user changes.

Finish with a concise report listing:
- implemented RoBMA integration points
- BayesTools APIs consumed
- key files changed
- tests run and results
- senior review rounds performed
- remaining limitations intentionally left for future work
```
