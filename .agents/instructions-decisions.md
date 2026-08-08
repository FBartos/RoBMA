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
  220 MB of fits created before BayesTools' current parameter registry.
- Impact: current summaries and plots reject them, preventing a normal vignette
  rebuild.
- Recommendation: regenerate one vignette cache at a time with current RoBMA and
  BayesTools. Do not add an in-memory migration for stale fitted objects.
- Decision: let's fully ignore vignette work till we are fully finished with all other work. I will explicitly say when we should deal with vignettes

## Cross-Level Point Equalities

- Issue: BayesTools point syntax intentionally requires a finite numeric literal
  on the right of `=`. Consequently `alloc[random] = alloc[systematic]` is not
  parsed as a point null. The explicit difference form
  `alloc[random] - alloc[systematic] = 0` reaches RoBMA, but direct fitted-model
  hypotheses currently reject all compound point expressions so that derived
  expressions cannot bypass declared product-space atoms. Marginal-means KDE
  hypotheses already support the explicit difference form because a
  single-model object has one common averaged posterior.
- Impact: accepting the shorthand expands the public hypothesis grammar.
  Supporting the derived point for fitted objects also requires an explicit
  atom-free certification rule. KDE can estimate an atom-free difference from
  draws; qCMDE/IWMDE additionally need a certified joint linear-target route and
  cannot reuse either level's marginal ordinate.
- Recommendation: extend the BayesTools AST to normalize symbolic equality to
  an explicit difference, and allow the resulting KDE point test only when
  posterior metadata certify that the derived target is atom-free. Continue to
  reject qCMDE/IWMDE for this compound target until its joint replacement and
  prior-ordinate route is certified. Keep model-averaged marginal means
  fail-closed because their levels can require different conditioning events.
- Investigation: the explicit difference briefly worked in commit `e348346b`
  because that implementation accepted compound KDE point expressions after
  checking only that the result looked scalar. Commit `114ac376` replaced it
  with the current fail-closed rule after showing that a derived expression
  could bypass declared product-space atoms. The single-model marginal-means
  route remains supported because both levels share one averaged, atom-free
  posterior. The fitted-model route was therefore removed for correctness, not
  lost accidentally during catalog work.
- Decision:
