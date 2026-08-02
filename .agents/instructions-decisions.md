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
- Decision:

## Cached Vignette Fits

- Issue: the 16 committed cache directories under `models/` contain about
  220 MB of fits created before BayesTools' current parameter registry.
- Impact: current summaries and plots reject them, preventing a normal vignette
  rebuild.
- Recommendation: regenerate one vignette cache at a time with current RoBMA and
  BayesTools. Do not add an in-memory migration for stale fitted objects.
- Decision:

## Covariance Validation Scope

- Issue: `V` is validated after explicit subsetting and missing-row removal, so
  a wholly excluded invalid block is not inspected.
- Impact: the fitted covariance is validated, but unused supplied entries are
  not certified.
- Alternatives: define `V` as describing the submitted dataset and validate it
  before row removal, or define it as describing retained observations and
  document current behavior.
- Recommendation: choose the public contract before changing validation order.
- Decision:

## STEP Fast-Route Selection

- Issue: `.selection_telescope_probability_check()` selects equivalent numerical
  algorithms using fixed `z = 8` and `omega = 100` thresholds.
- Impact: the statistical model is unchanged, but the accuracy/performance
  contract is implicit.
- Recommendation: use one canonical log-partition route or a documented,
  operation-specific error certificate. Do not retune constants without such a
  contract.
- Decision:

## Cluster Selected-Normal Quadrature

- Issue: cluster selected-normal AGHQ uses fixed finite-difference and curvature
  floors without a complete public accuracy/nonconvergence policy.
- Impact: difficult boundary or near-singular cluster likelihoods can depend on
  heuristic numerical choices.
- Recommendation: define quadrature order, convergence checks, failure behavior,
  and boundary handling before changing the implementation.
- Decision:
