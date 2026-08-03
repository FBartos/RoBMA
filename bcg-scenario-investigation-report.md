# BCG Scenario TODO Investigation

- Date: 2026-08-03
- Implementation worktree: `C:\R-Packages\RoBMA-bulk-tail-diagnostics`
- Branch: `codex/bulk-tail-density-diagnostics`
- Base: updated `dev-4-1-0` at `830f06a6`
- Original investigation worktree: `C:\R-Packages\RoBMA-bcg-investigation`

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
every display-grid ordinate. The display range starts from the empirical 0.5%
and 99.5% quantiles and adds 20% of that span on each side. For an approximately
normal posterior, this padding moves the endpoints to about the 0.015% and
99.985% quantiles. The separate `normalization_prob = 0.999` setting controls
the integration grid; it is not what made the display diagnostic fail.

At 500 rows, the padded endpoints had density only about 0.1%--0.4% of the curve
peak. The original worst raw diagnostics were:

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

Stan distinguishes rank-normalized bulk ESS from tail ESS and defines its tail
quantity from the 5% and 95% quantiles. RoBMA adopts those probability boundaries,
not Stan's ESS estimator: RoBMA's local ESS measures density-estimator importance
contributions rather than MCMC autocorrelation. See the
[Stan Reference Manual](https://mc-stan.org/docs/2_37/reference-manual-2_37.pdf)
and the [`posterior::ess_tail()` documentation](https://mc-stan.org/posterior/reference/ess_tail.html).

The implemented density-curve policy is:

- apply local relative MCSE, importance-contribution ESS, and largest-weight
  gates only over the empirical 5%--95% continuous posterior bulk;
- record the local diagnostics at the exact 5% and 95% tail checkpoints;
- retain `max(mcse) / peak_density` over the complete display as an absolute
  whole-curve safeguard;
- keep worst raw full-grid diagnostics attached for information, without letting
  padded extreme endpoints reject the curve; and
- leave exact requested-ordinate Bayes-factor gates unchanged and strict.

A fresh exact BCG allocation meta-regression produced:

| Target/method | MCSE / peak | Bulk rel. MCSE | Bulk ESS | Bulk max share | Raw rel. MCSE | Raw ESS | Raw max share |
|---|---:|---:|---:|---:|---:|---:|---:|
| location, qCMDE | 0.0264 | 0.0777 | 141.4 | 0.0166 | 0.4103 | 7.3 | 0.2922 |
| location, IWMDE | 0.0229 | 0.0507 | 273.0 | 0.0135 | 0.3070 | 8.9 | 0.2005 |
| heterogeneity, qCMDE | 0.0132 | 0.0323 | 355.9 | 0.0105 | 0.0970 | 96.0 | 0.0572 |
| heterogeneity, IWMDE | 0.0096 | 0.0253 | 627.0 | 0.0105 | 0.0714 | 184.1 | 0.0187 |
| allocation random, qCMDE | 0.0183 | 0.0618 | 176.3 | 0.0132 | 0.2434 | 18.6 | 0.1209 |
| allocation systematic, qCMDE | 0.0156 | 0.0512 | 233.9 | 0.0114 | 0.1751 | 37.2 | 0.0691 |
| allocation random, IWMDE | 0.0191 | 0.0432 | 325.7 | 0.0073 | 0.2206 | 18.4 | 0.1515 |
| allocation systematic, IWMDE | 0.0133 | 0.0358 | 443.1 | 0.0095 | 0.1346 | 46.3 | 0.0737 |

All six public plot combinations returned ggplot objects with zero warnings.
Their individual evaluation times were 2.33--5.96 seconds.

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

## Implemented changes

- Accept symbol-valued primitive nodes in the hypothesis warning matcher.
- Define the continuous density bulk by empirical 5% and 95% quantiles and
  include both cutpoints exactly in the deterministic display grid.
- Protect already inserted cutpoints and requested ordinates from being
  overwritten when multiple required coordinates share the same nearest grid
  slot.
- Gate local relative MCSE, importance-contribution ESS, and largest-weight
  share over the bulk; report both tail checkpoints separately.
- Retain whole-display peak-scaled absolute MCSE and raw full-grid diagnostics.
- Increase the default IWMDE density-curve budget from 500 to 1,000 rows;
  qCMDE remains 500 and requested ordinates remain adaptive up to all rows.
- Remove the superseded algorithm-5 scaled-ESS and scaled-weight fields.
- Bump the IWMDE algorithm provenance version from 5 to 6 and the package
  version to 4.1.3, preventing reuse of attributes under the old policy.
- Update roxygen documentation and add focused regression tests for method-
  specific defaults, bulk/tail calculation and every gate, AST symbols,
  deterministic checkpoint placement, and provenance invalidation.
- Update shared density-diagnostic test fixtures to the algorithm-6 format and
  assert that the marginalized-random-SD plot reaches the renderer with a
  finite qCMDE curve.

## Verification

- Standard cache refresh: 279 passed, 0 failed, 0 warnings, 21 documented
  certification skips, and all 38 standard cached fits validated. No new fit or
  catalog entry was added.
- Focused IWMDE, plotting, provenance, and marginal-means tests: 798 passed,
  0 failed, 0 warnings, and 28 documented certification, gallery, or unavailable
  optional-fixture skips.
- Standard profile using the 38 validated shared fits: 8,993 passed, 0 failed,
  0 warnings, and 156 documented unavailable optional-fixture, certification,
  or extended-gallery skips; runtime was 382.7 seconds.
- Fresh exact BCG regression fit: all six qCMDE/IWMDE public plots completed
  without warnings in 2.33--5.96 seconds per call; bulk, tail-checkpoint,
  full-grid, and peak-scale metrics were inspected as reported above.
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
- A genuinely unstable curve still fails through whole-display absolute MCSE
  or through relative MCSE, ESS, or concentrated importance contributions in
  the 5%--95% bulk.
- The 5% and 95% tail checkpoints are deterministic grid points and remain
  directly inspectable; only the more extreme padded display points are
  informational for local relative diagnostics.
- Strict local checks for requested BF ordinates are unchanged.
- IWMDE plotting uses twice as many posterior rows by default. In the BCG
  regression checks, total per-call times remained below six seconds, but the
  estimator portion can approach twice the previous cost in other models.
- Algorithm-version 6 prevents reuse of density attributes produced under the
  previous reliability/default policy.
- Existing human-reviewed scenario SVGs are intentionally untouched.

## Recommendation

1. Merge the algorithm-6 bulk/tail diagnostic policy with the already reviewed
   AST-symbol fix and 1,000-row IWMDE curve default. The BCG evidence does not
   justify increasing qCMDE to 5,000 rows or increasing IWMDE beyond 1,000.
2. Do not add adaptive curve budgets yet. The current fixed budgets pass the
   meaningful bulk/tail and absolute-error checks within the existing runtime;
   adaptation would add complexity without an observed failure.
3. After merging the policy, update `tests/scenarios/bcg.R` to remove the
   resolved comments and add bridge comparisons for the intercept and both
   allocation coefficients. Review any resulting SVG changes manually.
4. Add the heterogeneity-zero bridge comparison only to the certification
   profile, with its runtime and Monte Carlo uncertainty recorded. Do not add it
   to the standard suite or ordinary BCG scenario.
