# Pending Maintainer Decisions

Keep only unresolved choices here. Once decided, implement the decision and
remove the item.

## Visual Baseline Refresh

- Issue: numerical corrections and graphics-toolchain changes produced ignored
  `.new.svg` files across selected-normal, radial, forest, and other plot
  families. Some forest output also contains stochastic prediction intervals.
- Impact: affected certification snapshots cannot be treated as verified
  baselines.
- Recommendation: review intended differences manually, make stochastic
  snapshots deterministic, and accept only verified baselines. Do not hide
  differences with broad coordinate tolerances.
- Decision: I will review them

## Exact Selection Normalizer for Large Multilevel Fits

- Issue: the finite-vector exact selection likelihood checks the randomized-QMC
  normalizer at every evaluated parameter state. For the 412-estimate,
  128-study `Johnides2025` RoBMA model (largest block 18), 4,096 points times 8
  scrambles failed the 0.5% relative-MCSE requirement. Increasing to 16,384
  points times 8 scrambles still failed across three independent
  initializations. A current-binary diagnostic probe measured 3.3959% relative
  MCSE at the first rejected state. The multilevel `Weingarten2018` RoBMA fit
  in v36 also exhausted its ordinary restarts: observed failures were 0.5639%
  at 512 points and 0.5912% at 1,024 points. Minimum-length one-restart probes
  failed at 0.5085% with 2,048 points and 0.5637% with 4,096 points, showing
  that fixed-budget doubling is not uniformly or monotonically resolving the
  states visited by the sampler. A minimum-length probe using v33's actual
  priors and the default 512-point control failed at 3.3304% relative MCSE.
- Attempts: increased only the integration-point budget while preserving the
  exact likelihood and tolerance; removed duplicate JAGS transfer of identical
  QMC designs across blocks of the same size; installed current source binaries
  in an isolated library so PSOCK workers used the instrumented JAGS module;
  reproduced the failure with a one-restart, minimum-length probe. The data
  deduplication reduced preparation from about four hours to about 13 minutes,
  but does not address estimator variance.
- Impact: three current-source exact-selection cache refreshes cannot be
  certified: `v32-robma-multilevel/fit`,
  `v33-robma-multilevel-metaregression/fit_reg`, and
  `v36-zplot/fit_RoBMA_Weingarten2018`. The other 45 declared vignette caches
  were regenerated and validated with current parameter maps. The three old
  cached objects cannot render with the current package because they predate
  parameter-map metadata. Treating a failed
  integration check as a rejected MCMC proposal, loosening the tolerance, or
  switching silently to the approximate likelihood would change the numerical
  or statistical contract.
- Alternatives: (1) a deterministic adaptive QMC normalizer that adds points
  until the requested error is met and fails only at an explicit maximum; (2)
  a lower-variance selected-Gaussian probability algorithm, such as a
  mathematically justified ordering/pivoting extension of the current
  sequential estimator; or (3) explicitly mark the historical v32/v33
  manuscript reproductions as using `selection_likelihood = "approximate"`.
  Brute-force fixed budgets appear unsuitable: under ordinary Monte Carlo
  scaling, reducing 3.3959% to 0.5% would require roughly 46 times the current
  independent effort for the rejected state.
- Recommendation: obtain mathematical review of alternatives (1) and (2).
  Do not weaken the 0.5% criterion or alter proposal handling. Use alternative
  (3) only if the maintainer decides those two vignettes should reproduce the
  previously released approximate model rather than demonstrate the new exact
  likelihood.
- Decision:

## Correlation in Homogeneous Structured Heterogeneity Summaries

- Issue: the Ishak scenario shows `cor` in `summary_heterogeneity()` for
  `har(time | study)`, but `summary_heterogeneity()` for the homogeneous
  `ar(time | study)` model prints only `sd` and `var`, even though `cor` is a
  fitted public parameter and is shown by `summary()`.
- Impact: otherwise parallel HAR and AR output exposes different parameter
  families and makes the nested comparison less systematic.
- Recommendation: include `cor` in `summary_heterogeneity()` for homogeneous
  scalar-correlation structures such as AR, CS, and CAR. Keep the existing
  `sd`/`var` naming and do not add an RoBMA-side naming workaround.
- Decision:
