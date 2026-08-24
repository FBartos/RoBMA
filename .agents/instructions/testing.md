# Testing

Use this guide for changes under `tests/testthat/` or the ordinary test-profile
runners under `tools/`. For maintainer analysis scenarios under
`tests/scenarios/` or `tools/test-scenario.R`, use `scenarios.md` instead.

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

For an intentional full run from an interactive session, source
`.dev/user-tests.R` by its absolute path. The script works from any current
working directory and calls `test_tests(refit = TRUE,
stop_on_failure = TRUE)`. `test_tests()` is loaded by the project `.Rprofile`.
Without a filter it refreshes and runs the standard profile, then runs every
certification case with its independent one-hour limit. A regular-expression
`filter` runs only matching standard test files and does not refresh cached
fits unless `refit = TRUE`; filtered runs omit certification. `refit = TRUE`
cleans the standard fit cache before rebuilding; `update = TRUE` permits
missing visual candidates, and `regenerate = TRUE` combines those controls.
Every run ends with snapshot review. Ordinary tests have no per-file timing
baseline, so `update_timings = TRUE` fails clearly; their timing contract is the
standard profile's 15-minute total budget. Interactive calls default to
testthat's progress reporter. Use `reporter = "llm"` for agent-oriented output;
the quiet LLM reporter prints failures and warnings immediately, retains skip
counts, and omits individual skip reports. An unfiltered run is deliberately
much slower than the ordinary standard profile and does not run
`devtools::check()`.

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
`ROBMA_TEST_QUIET_SKIPS=TRUE` retains skip counts while suppressing individual
skip reports in profile-runner output.

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

After sourcing `tests/scenarios/helper-scenarios.R`, call
`review_test_snapshots()` to open testthat's native reviewer for ordinary text
and figure candidates under `tests/testthat/`, followed by per-artifact text and
table candidates under `tests/results/`. Its optional `files` argument has the
same test-name or trailing-slash directory semantics as
`testthat::snapshot_review()`; use `reference_filter` to select reference
groups such as `"interpret"` or `"marginal_means"`. Accept changes only after
maintainer review. Rejected `.new.txt` candidates are deleted and skipped ones
remain available for a later review.

Set `ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS=TRUE` only during an explicit snapshot
regeneration workflow.

## Final Verification

Run focused tests after each meaningful change. Before handoff, run the standard
profile when the change scope warrants it. State exactly what ran, what did not,
and why.
