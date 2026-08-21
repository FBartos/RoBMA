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

## Cached Vignette Fits

- Issue: the 16 committed cache directories under `models/` contain about
  220 MB of fits created before BayesTools' current coordinate-registry and
  semantic-catalog contracts.
- Impact: current summaries and plots reject them, preventing a normal vignette
  rebuild.
- Recommendation: regenerate one vignette cache at a time with current RoBMA and
  BayesTools. Do not add an in-memory migration for stale fitted objects.
- Decision: let's fully ignore vignette work till we are fully finished with all other work. I will explicitly say when we should deal with vignettes

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
