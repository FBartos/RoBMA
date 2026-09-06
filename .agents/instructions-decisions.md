# Pending Maintainer Decisions

Keep only unresolved choices here. Once decided, implement the decision and
remove the item.

## Existing Covariance Spectral Truncation

- Issue: `R/covariance-factorization.R` symmetrizes its input and replaces
  spectral eigenvalues whose absolute value is within its roundoff envelope
  by zero, including small positive eigenvalues. The selected-response native
  fallback in `src/r-selnorm-mv.cc.inc` has the same truncation policy after
  Cholesky fails. This policy predates the selection-speed work.
- Impact: the spectral factor can differ from the submitted covariance. The
  R helper is consumed by ordinary known-V latent decomposition and covariance
  sampling, not just selection models. Its comments distinguish preserved
  covariance storage from altered spectral factors; storage preservation alone
  does not establish that downstream sampling preserves the target.
- Recommendation: review the shared policy separately. Preserve a successful
  Cholesky factor and valid positive eigenvalues; decide how to handle
  roundoff-sign ambiguity for genuinely singular inputs before changing the
  spectral fallback. Do not silently add jitter or covariance repair.
- Decision: pending; the current performance changes leave this policy intact.

## Existing Known-V Residual-Decomposition Policy

- Issue: `.known_v_decompose_block()` limits the latent backend's residual
  fraction to `0.99 * lambda_min(cov2cor(V))` and rejects a maximum fraction
  at or below `sqrt(.Machine$double.eps)`. It warns when the requested
  fraction is reduced; these policies predate the current work.
- Impact: the decomposition preserves the Gaussian covariance (subject to
  the spectral policy above), but changes which factors the approximate
  selection likelihood conditions on. The 0.99 margin and near-singularity
  rejection also exclude some mathematically valid decompositions. The
  Assink ordinary-matrix target uses the requested 0.1 fraction without this
  reduction, so these constants do not explain the current Assink timing.
- Recommendation: review the residual-fraction contract and singular-input
  policy together before changing these constants. Do not change the
  approximate target or add covariance repair as a performance shortcut.
- Maintainer proposal: reject approximate selection when the requested
  decomposition is infeasible instead of reducing its residual fraction, and
  suggest 'selection_likelihood = "exact"' where supported. Base rejection on
  the requested decomposition, not a generic high-dependence cutoff: another
  declared factorization can be feasible. This does not resolve the shared
  spectral policy or guarantee that exact selection supports every singular V.
- Decision: pending; left unchanged in this work.

## Existing BayesTools Truncated-Normal CDF Clipping

- Issue: BayesTools `R/priors-methods-joint.R` clips its truncated-normal CDF
  and complementary CDF to `[0, 1]` with `pmin`/`pmax`. This is preexisting
  code adjacent to the log-density implementation inspected during profiling.
- Impact: overshoot from its interval-probability calculation is hidden rather
  than diagnosed. The bridge-prior optimization uses log densities, not these
  CDF clipping branches; this finding does not explain the measured bottleneck.
- Recommendation: review this policy with the other numerical-boundary
  decisions. No inference-changing error was established here; do not alter
  unrelated prior-CDF behavior as part of the performance change.
- Decision: pending; left unchanged in this work.

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
