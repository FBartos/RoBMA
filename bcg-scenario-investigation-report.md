# BCG Scenario TODO Investigation

- Date: 2026-08-03
- Worktree: `C:\R-Packages\RoBMA-bcg-investigation`
- Branch: `codex/bcg-scenario-investigation`
- Base: `dev-4-1-0` at `47f9740c`

## Scope

The five TODO comments in `tests/scenarios/bcg.R` reduce to three questions:

1. Why are qCMDE/IWMDE location-density overlays sometimes taller than KDE?
2. Why do the regression location, heterogeneity, and allocation plots report
   high Monte Carlo standard errors (MCSEs)?
3. Can regression Bayes factors be checked against independently fitted null
   models?

This investigation did not regenerate scenario baselines or edit the scenario.
Those changes require maintainer review of intentional visual changes.

## Findings

### 1. Taller simple-model density overlays

With a current-profile refit of the exact BCG simple model, the location curves
were:

| Method | Numerical area | Peak density |
|---|---:|---:|
| KDE | 0.9999998 | 2.290366 |
| qCMDE | 1.001593 | 2.313534 |
| IWMDE, 1,000 rows | 1.000601 | 2.289634 |

The qCMDE line is about 1% taller than KDE but has the correct total mass. It is
a conditional-density estimate rather than a kernel-smoothed estimate, so equal
samples and equal mass do not imply equal peaks. A KDE with half the default
bandwidth reached a peak of 2.3456, confirming that a difference of this size is
ordinary smoothing behavior, not a scale error.

The old 500-row IWMDE curve had area 1.0226 and peak 2.3536 for this seeded fit.
Increasing grid points did not remove that lift. Increasing estimator rows did:

| IWMDE rows | Numerical area | Peak density |
|---:|---:|---:|
| 500 | 1.02263 | 2.3536 |
| 1,000 | 1.00060 | 2.2896 |
| 2,000 | 1.00319 | 2.2938 |
| 5,000 | 1.00240 | 2.3047 |

Cause: expected KDE smoothing for qCMDE, plus visible finite-row Monte Carlo
variation for IWMDE at 500 rows. There is no normalization or parameter-scale
bug.

### 2. Regression plot MCSE warnings

The old gate took the worst relative MCSE, ESS, and importance-weight share over
every display-grid ordinate. The failures came from extreme display tails near
the 0.999 quantile, where density was only about 0.1%--0.4% of the curve peak.
At 500 rows the worst raw diagnostics were:

| Target/method | Worst relative MCSE | Minimum ESS | Largest weight share | Maximum MCSE / peak |
|---|---:|---:|---:|---:|
| location, qCMDE | 0.4517 | 5.3 | 0.3932 | 0.0249 |
| heterogeneity, qCMDE | 0.1049 | 78.9 | 0.0290 | 0.0093 |
| allocation random, qCMDE | 0.3291 | 9.2 | 0.2772 | 0.0262 |
| allocation systematic, qCMDE | 0.2194 | 19.6 | 0.1584 | 0.0200 |
| location, IWMDE | 0.4316 | 5.4 | 0.3813 | 0.0335 |
| heterogeneity, IWMDE | 0.1084 | 71.3 | 0.0435 | 0.0203 |
| allocation random, IWMDE | 0.3584 | 8.0 | 0.3016 | 0.0378 |
| allocation systematic, IWMDE | 0.2026 | 22.3 | 0.1308 | 0.0243 |

Thus the worst pointwise relative error looked large, while the largest
absolute MCSE occupied at most 3.8% of the plot's vertical scale. Raising the
budget to 2,000 rows cost roughly four to five times as much and still did not
make every negligible tail pass the old ESS gate.

The prototype uses plot-scale density diagnostics:

- MCSE: `max(mcse) / peak_density`;
- ESS: `min(ess / (density / peak_density)^2)`;
- largest-weight influence: `max(weight_share * density / peak_density)`.

Raw pointwise diagnostics remain attached for inspection. Only the gate used to
decide whether an entire curve is usable changes. Requested-ordinate Bayes
factor checks retain their stricter local diagnostics.

On a current-profile refit, all six public BCG regression plot combinations
(location, heterogeneity, allocation crossed with qCMDE/IWMDE) returned ggplot
objects with zero warnings. Their individual evaluation times were 2.38--5.80
seconds.

### 3. Regression Bayes-factor null-model comparisons

An independent bridge-marginal-likelihood comparison is feasible. The tested
results were:

| Hypothesis | Bridge BF | qCMDE BF | IWMDE BF | qCMDE log difference | IWMDE log difference |
|---|---:|---:|---:|---:|---:|
| intercept = 0 | 2.21768 | 2.24401 | 2.43420 | 0.0118 | 0.0932 |
| allocation random = 0 | 0.99195 | 0.99669 | 0.98148 | 0.0048 | -0.0106 |
| allocation systematic = 0 | 0.68841 | 0.70337 | 0.70466 | 0.0215 | 0.0233 |
| heterogeneity = 0 | 2.12035e22 | 2.02623e22 | 2.08810e22 | -0.0454 | -0.0153 |

The allocation-coefficient oracle requires replacing the factor by equivalent
numeric indicator columns and setting
`standardize_continuous_predictors = FALSE`. Leaving standardization enabled
changes the joint coefficient prior, so it is not the same null-model oracle.

The heterogeneity boundary used all 15,000 allowed rows, took about 62 seconds
for qCMDE and 39 seconds for IWMDE, and retained about 22% estimated relative
BF error. It is useful certification evidence but too expensive and noisy for
the ordinary maintainer scenario. The intercept and allocation comparisons are
appropriate routine scenario checks.

During this work, the current BayesTools hypothesis AST exposed a separate bug:
primitive expressions are symbols, but the warning matcher attempted to subset
them as calls. The prototype safely converts both symbols and calls to text.
This affects warning attachment after a BF has already been computed; it does
not change the BF calculation.

## Prototype changes

- Accept symbol-valued primitive nodes in the hypothesis warning matcher.
- Evaluate density-curve MCSE, ESS, and largest-weight gates on the plot's
  vertical scale while retaining raw diagnostic fields.
- Increase the default IWMDE density-curve budget from 500 to 1,000 rows;
  qCMDE remains 500 and requested ordinates remain adaptive up to all rows.
- Bump the IWMDE algorithm provenance version from 4 to 5.
- Update roxygen documentation and add focused regression tests for method-
  specific defaults, plot-scale metric calculation/gating, AST symbols, and
  provenance invalidation.
- Update shared density-diagnostic test fixtures to the algorithm-5 format and
  assert that the marginalized-random-SD plot reaches the renderer with a
  finite qCMDE curve.

## Verification

- Focused IWMDE, plotting, and marginal-means tests: 783 passed, 0 failed,
  0 warnings, 28 documented skips.
- Standard profile using 38 validated shared fits: 8,980 passed, 0 failed,
  0 warnings, 156 documented skips for unavailable/stale optional fixtures and
  certification or extended-gallery cases.
- Fresh exact BCG regression fit: all six qCMDE/IWMDE public plots completed
  without warnings in 2.38--5.80 seconds per call.
- Fresh exact BCG simple fit: KDE, qCMDE, and IWMDE curve areas and peaks were
  measured through the public ggplot path as reported above.
- Source-package check, excluding tests and vignettes already handled
  separately: 0 errors, 0 warnings, 0 notes. The check used `LC_ALL=C`,
  `LC_CTYPE=C`, and `LANG=C` because this Windows R installation rejects the
  shell's `C.UTF-8` locale.
- A vignette-building attempt reproduced the existing stale-cache failure:
  committed vignette fits lack BayesTools' canonical parameter registry. No
  cache migration or baseline regeneration was attempted.

## Consequences

- Valid central curves are no longer rejected solely by negligible display
  tails.
- A genuinely unstable visible curve still fails through absolute plot-scale
  MCSE, scaled ESS, or scaled weight influence.
- Strict local checks for requested BF ordinates are unchanged.
- IWMDE plotting uses twice as many posterior rows by default. In the BCG
  regression checks, total per-call times remained below six seconds, but the
  estimator portion can approach twice the previous cost in other models.
- Algorithm-version 5 prevents reuse of density attributes produced under the
  previous reliability/default policy.
- Existing human-reviewed scenario SVGs are intentionally untouched.

## Recommendation

1. Accept the AST-symbol fix independently; it is a compatibility correctness
   fix required by the current BayesTools contract.
2. Accept plot-scale reliability gates and the 1,000-row IWMDE curve default.
   This is the smallest tested change that addresses both misleading regression
   warnings and the visible 500-row IWMDE fluctuation without weakening BF
   validation.
3. After accepting the policy, update `tests/scenarios/bcg.R` to remove the
   resolved comments and add bridge comparisons for the intercept and both
   allocation coefficients. Review any resulting SVG changes manually.
4. Add the heterogeneity-zero bridge comparison only to the certification
   profile, with its runtime and Monte Carlo uncertainty recorded. Do not add it
   to the standard suite or ordinary BCG scenario.
5. Defer package-version and `NEWS.md` changes until the maintainer selects the
   prototype changes for `dev-4-1-0`; this branch is an investigation candidate,
   not an accepted feature branch.
