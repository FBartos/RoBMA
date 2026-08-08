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
- Decision:

## Fixed qCMDE/IWMDE Ordinate Budget

- Issue: the review replaced contribution-dependent adaptive row stopping with
  one state-independent fixed sample of 500 rows. This preserves a valid
  simple-random-sampling design, but exposes substantially higher MCSE in the
  BCG regression. On the 15,000-row `ablat = 0` ordinate, IWMDE MCSE was 10.9%,
  8.28%, 6.13%, 5.20%, and 3.79% for 500, 1,000, 2,000, 4,000, and all rows;
  qCMDE MCSE was 9.12%, 7.67%, 5.55%, 4.94%, and 3.73%. Corresponding elapsed
  times were 5.7/2.9, 7.1/4.6, 9.0/8.9, 15.6/20.8, and 83/157 seconds for
  IWMDE/qCMDE.
- Impact: raising the public default to 4,000 multiplies representative runtime
  by roughly 3--7 and still leaves IWMDE just above the 5% target. Using all
  rows removes row-sampling error but costs 1.4--2.6 minutes per ordinate here.
  Restoring adaptive stopping would make the chosen sample size depend on the
  evaluated contributions and undo the reviewed sampling design.
- Recommendation: retain the public 500-row default and its explicit warning.
  Give the BCG scenario an explicit larger fixed budget (or the all-row census
  when it is intended as numerical certification) after choosing its acceptable
  runtime, then review the changed text baselines. Do not restore the old
  adaptive loop.
- Decision:
