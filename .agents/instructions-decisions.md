# Pending Maintainer Decisions

Keep only unresolved choices here. Remove an item when its decision is
implemented or superseded.

## Selected Funnel and Zplot Target Consistency

- Issue: the marginal new-effect contract is implemented by selected zplot as
  an outer average over newly realized retained contexts of conditionally
  normalized full publication-event marginals. In contrast,
  `R/funnel-quantiles.R` combines random sources into total marginal SD before
  calling the scalar selected-CDF API in `R/regplot-quantiles.R`. Source roles,
  publication groups and vector selection rules do not reach that native API.
- Impact: `funnel()` and `bfunnel()` differ correctly in plug-in versus
  posterior averaging, but their selected distributions generally disagree
  with the fitted source-normalization order. This already occurs for diagonal
  sampling covariance with a retained study effect. The documented descriptive
  scalar-SE convention for full V additionally omits the original multivariate
  selection event. Previous post-fit runtime passes did not verify this
  consistency and must not be interpreted as such.
- Recommendation: retain zplot's agreed marginal replicated-literature target
  and share its source-aware selected marginal law with contour calculations.
  For arbitrary hypothetical SE grids under full V, define the accompanying
  multivariate design/publication event explicitly; a scalar SE alone does not
  determine that law. Options are full-event contours at validated designs or
  a separately identified descriptive scalar display. Do not silently replace
  full-event selection by scalar selection or condition on fitted BLUPs.
- Current implementation: fully conditioned positive-weight contours use the
  Gaussian law. Other configurations unsupported by the scalar calculation
  fail explicitly, rather than returning a different selected distribution.
- Decision: pending discussion of contour design semantics. The proposed
  convention scales the complete original publication's SE vector proportionally
  at each hypothetical focal SE, preserving correlations, relative SEs and random
  design, then averages over observed designs. The maintainer question is pending.

## Covariance Working-Precision Boundary

- Issue: dense covariance inputs with standardized eigenvalues within the
  eigensolver roundoff envelope cannot be classified as exactly singular,
  positive definite or indefinite from ordinary floating-point output alone.
- Current implementation: the inherited numerical PSD rule now operates in SRS
  correlation coordinates. One accepted factor defines sampling, whitening and
  conditional support; a rounded positive Cholesky pivot cannot override it.
  Small independent variances are preserved. Known-group kernels requiring a
  nonsingular density fail explicitly when positive definiteness is unresolved.
  External numerical-symmetry canonicalization remains necessary for ordinary
  metafor::vcalc output and is documented. No diagonal jitter was introduced.
- Evidence: exact integer rank-two matrices had eight support failures over 30
  row orderings before the shared-support correction and none afterward. Keeping
  every computed positive eigenvalue or using zero-tolerance pivoted Cholesky
  did not solve the general support problem.
- Decision: any replacement of the inherited within-roundoff PSD convention
  with stricter unavailability for all ambiguous dense inputs remains a separate
  maintainer choice. The current review makes the existing convention consistent;
  it does not claim exact binary64 PSD certification.

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
