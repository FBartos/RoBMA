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
