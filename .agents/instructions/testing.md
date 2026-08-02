# Testing

Use this guide for changes under `tests/testthat/`, `tests/scenarios/`, or the
test runners under `tools/`.

## Development Workflow

Always use the LLM reporter:

```r
devtools::test(filter = "topic", reporter = "llm")
```

Run the narrowest relevant test first. Populate affected cached fits only when
needed:

```r
devtools::test(filter = "01-", reporter = "llm")
```

Do not use repeated full-suite runs as an iteration loop.

## Test Profiles

- `Rscript tools/test-profile.R standard`: routine unit, integration,
  representative metafor, and human-reviewed visual tests using an already
  valid fit cache; maximum 15 minutes on the reference machine. It never fits
  models.
- `Rscript tools/test-profile.R refresh-standard`: create only missing or stale
  standard cached fits. Use `--clean` only for an intentional full refresh.
- `Rscript tools/test-profile.R certification --list`: list independently
  runnable expensive numerical and extended-regression cases. Run named cases
  by appending their names; omitting names runs every case in sequence.
- `Rscript tools/test-profile.R release`: refresh standard fits, run the
  standard suite, run all certification cases, then call `devtools::check()`.

Every certification case has a hard one-hour limit. Certification has no total
limit because cases execute in separate processes and can be selected or rerun
independently.

## Test Organization

- `test-00-*`: deterministic validation and focused low-level tests.
- `test-01-*`: model fits and cache generation.
- `test-02-*`: post-fit methods using cached fits.
- `test-03-*`: higher-level integration.

Only `test-01-*` files create ordinary cached model fits. Reuse existing fits
elsewhere. Check `fit_catalog()` before adding a fit, and update the catalog in
the same change when a new cached fit is necessary.

Source `common-functions.R` in files using cached fits. Prefer
`list_fits()`, `catalog_fits()`, `lazy_fits()`, and `lazy_infos()` over
hard-coded bulk loading. Missing or stale required fits are not passing release
evidence.

## Cache Rules

Standard and certification caches are isolated below
`ROBMA_TEST_FILES_DIR`, or under `tests/testthat/test_files/` by default.
Cache validity includes RoBMA, BayesTools, R, and JAGS versions; the fitting
test hash; cache-affecting RoBMA source; BayesTools' backend fingerprint; and
required attached results.

Regenerate only affected fits with `refresh-standard` after fitting, prior/data
input, native likelihood, or JAGS changes. Do not invalidate fits for unrelated
plotting, summary, documentation, or post-fit changes.

CI restores the standard cache using a key derived from cache-affecting package
sources, fitting tests, R/JAGS versions, and the BayesTools backend fingerprint.
It refreshes missing or stale entries before running the strictly cached
standard profile. Certification caches remain local or release-worker assets.

Relevant controls are `ROBMA_TEST_PROFILE`, `ROBMA_TEST_FILES_DIR`,
`ROBMA_TEST_SKIP_REFIT`, and `ROBMA_TEST_FORCE_REFIT`.

## Correctness Evidence

- Test behavior, not implementation trivia.
- Use analytic identities or independent reference implementations for numerical
  kernels.
- Compare against metafor where estimands match. Save the metafor object in the
  cached fit's `info` during `test-01-*`, then reuse it in post-fit tests.
- Justify tolerances from Monte Carlo or numerical error. Do not use a fixed
  package-wide tolerance merely because it makes a test pass.
- Do not update committed expected results when a test fails. Determine whether
  the implementation or old expectation is wrong and involve the maintainer
  before changing a verified baseline.
- Do not add redundant matrices, samples, fits, or assertions for coverage alone.

The evidence inventory in `tests/testthat/REGRESSION-COVERAGE.md` records the
representative oracle and visual coverage expected for public post-fit methods.
Update it when a method changes tier or loses an oracle.

## Visual Regression

Use `expect_vdiffr_snapshot()` from
`tests/testthat/helper-visuals.R`. Representative, human-reviewed snapshots
belong in the standard profile. Gate redundant galleries with
`skip_if_not_full_visuals()`.

Never auto-accept or auto-update snapshots. Structural `as_data = TRUE` tests
supplement visual snapshots; they do not replace them. Ask the maintainer to
review every intentional visual change.

Set `ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS=TRUE` only during an explicit snapshot
regeneration workflow.

## Maintainer Scenarios

Maintainer analysis scenarios live in `tests/scenarios/` and are excluded from
routine package tests. Use `scenario_fit()`, `scenario_text()`, and
`scenario_plot()` after `scenario_start()`; do not duplicate their cache or
snapshot logic in individual scenarios.

Run one with `Rscript tools/test-scenario.R <name>`. Set
`REGENERATE_SCENARIO_FILES <- TRUE` in the scenario, or pass `--regenerate`, to
refit its models and replace all exercised text and plot baselines. Review every
resulting diff.

Fit caches under `tests/scenarios/cache/` are local and ignored. Text baselines
under `results/` and vdiffr SVGs under `_snaps/` are committed. Missing
baselines may be created locally but must fail on CI.

## Final Verification

Run focused tests after each meaningful change. Before handoff, run the standard
profile when the change scope warrants it. State exactly what ran, what did not,
and why.
