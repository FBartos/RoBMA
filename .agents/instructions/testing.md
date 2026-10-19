# Testing

Use this guide for changes under `tests/testthat/` or the ordinary test-profile
runners under `tools/`. For maintainer analysis scenarios under
`tests/scenarios/` or `tools/test-scenario.R`, use `scenarios.md` instead.

In BayesToolsVerse, the shared validation guide defines tests, verification,
and scenarios. The runner's `certification` profile is deeper verification,
separate from routine tests. Use the workspace's configured R and private agent
library when available.

## Development Workflow

Always use the LLM reporter:

```r
devtools::test(filter = "topic", reporter = "llm")
```

Run the narrowest relevant test first. Populate affected cached fits only when
needed:

```r
devtools::test(filter = "^01-", reporter = "llm")
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
baseline, so `update_timings = TRUE` fails clearly; the standard profile's
15-minute total limit is an extreme ceiling, not a runtime target. Interactive
calls default to
testthat's progress reporter. Use `reporter = "llm"` for agent-oriented output;
the quiet LLM reporter prints failures and warnings immediately, retains skip
counts, and omits individual skip reports. An unfiltered run is deliberately
much slower than the ordinary standard profile and does not run
`devtools::check()`.

## Test Profiles

- `Rscript tools/test-profile.R standard`: routine unit, integration,
  representative metafor, and human-reviewed visual tests using an already
  valid fit cache; its 15-minute maximum on the reference machine is a ceiling,
  not a target. It never fits models.
- `Rscript tools/test-profile.R refresh-standard`: create only missing or stale
  standard cached fits. Use `--clean` only for an intentional full refresh.
- `Rscript tools/test-profile.R certification --list`: list independently
  runnable expensive numerical and extended-regression cases. Run named cases
  by appending their names; omitting names runs every case in sequence.
- `Rscript tools/test-profile.R release`: refresh standard fits, run the
  standard suite, run all certification cases, then call `devtools::check()`.

Every profile checks the loaded DLL against the native entry points declared
in `src/init.c` before running tests. The complete standard profile also requires
the core oracle and real-estimator cases named by `standard_required_tests()`
to execute with passing expectations. Missing symbols and skipped required cases
are failures; optional dependency and extended-gallery skips remain visible and
permitted. Filtered development runs enforce native availability but do not claim
complete standard evidence.

`tools/profile-iwmde-qcmde.R` prepares the release DLL from its own source
root, disables pkgload compilation and test helpers, and defaults to one native
thread. `--list` reports available and omitted cached-fit cases without running
the estimator. Its output records source/package/DLL identity and labels
elapsed values as including Rprof sampling and allocation profiling. An explicit
`ROBMA_DEBUG_DLL=TRUE` retains the existing DLL and is labelled a debug opt-out.

The complete standard profile gates its test phase at 900 seconds, after
cache validation, including `run_tests()` preparation and required-evidence checks;
the runner separately reports the full profile wall time.
Standard uses the explicit `profile == "standard"` catalog membership, rather
than every `tier == "core"` entry. Certification extends the catalog and
selects its own declared fit subsets; excluded core entries are not standard
evidence. Shared GLMM settings use 1500 samples per chain in standard and 5000
in certification, but `bcg_BMA.glmm` deliberately retains its explicit
three-chain, 1000-sample, 2000-burn-in, 500-adaptation, thin-160 settings.
These profile descriptions do not change any fitting budget.

Every certification case has a hard one-hour limit shared by its cache
preparation and verification phases. Each phase owns its child process tree and
terminates that tree before waiting or draining output at the deadline. The phases execute in separate processes,
so native fitting state cannot leak into post-fit verification. Certification
has no total limit because cases can be selected or rerun independently.
The runners start child R processes with `--vanilla`: they inherit the calling
process's test controls, which a user or site `.Renviron` read again by the
child would otherwise override (for example, redirecting an isolated
`ROBMA_TEST_FILES_DIR` to the default cache).

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
required attached results. A source file is cache-affecting when it executes
while a model is constructed (constructors, data and prior preparation, JAGS
model text, native code of the JAGS module and the R-native routines they
call), while the sampler runs, or while the LOO, WAIC, and marginal-likelihood
results stored with a fit are computed; the list is
`.fit_cache_required_source_files()`, and post-fit methods stay out of it. A
generating test file that reads `test_glmm_fit_settings()` hashes the settings
it returns.

Producer provenance is captured before builders run, once per source checkout
and loaded-product generation. Re-sourcing `common-functions.R` retains that
generation. The profile runner explicitly begins a generation after loading the
product; ordinary devtools tests initialize it on their first helper source.
Saving refuses changed package/build or builder fingerprints before writing any
artifact: reload the product and rerun generation after editing those sources.
This lifecycle cannot retroactively certify a manually loaded namespace whose
sources were already stale before the helper initialized.

The catalog explicitly declares consumed fit/info artifacts for derived caches.
Both known-V IWMDE oracle children depend on the random model's info (including
data and V) and the fixed-random-null fit. Validated loads record exact artifact
bytes and parent provenance in the current generation; memoized objects retain
their own receipt and reload when those inputs change. A save requires those
receipts, refuses changed or conflicting consumption, and persists them in the
child metadata. Validation recursively checks declared parents and required
artifacts, refusing cycles and old child metadata without the dependencies.
Unvalidated loads and receipts from a previous generation cannot authorize a
derived save. Ordinary fits without declared dependencies keep their cache rules.
Present info artifacts that cannot be read fail before consumption is recorded;
truly absent optional info still returns an empty list.

Regenerate only affected fits with `refresh-standard` after fitting, prior/data
input, native likelihood, or JAGS changes. Do not invalidate fits for unrelated
plotting, summary, documentation, or post-fit changes.

CI restores the standard cache using a key derived from cache-affecting package
sources, fitting tests, R/JAGS versions, and the BayesTools backend fingerprint.
It refreshes missing or stale entries before running the strictly cached
standard profile. Certification caches remain local or release-worker assets.

Relevant controls are `ROBMA_TEST_PROFILE`, `ROBMA_TEST_FILES_DIR`,
`ROBMA_TEST_SKIP_REFIT`, `ROBMA_TEST_FORCE_REFIT`, and the case-local
`ROBMA_TEST_ACTIVE_FITS` filter. `tools/test-cache-key.R` clears that inherited
filter because its key describes the whole profile. Explicitly empty
`ROBMA_TEST_PROFILE` remains invalid. Source hashes use binary UTF-8/LF text
and round-trippable numeric deparse; empty text retains its existing semantics.
Cache-key changes require one integrated schema advance before final refresh;
old baseline artifacts are preserved rather than resealed.
`ROBMA_TEST_QUIET_SKIPS=TRUE` retains skip counts while suppressing individual
skip reports in profile-runner output.

Two session caches persist across the tests of one R session: the plans of
`hypothesis()`, `hypothesis_quantities()`, and the marginal-means methods (per
object), and the display-independent qCMDE normalizations (per fit, target, and
control). A reload of the RoBMA or the BayesTools code empties both. A test
that mocks a function they are built with between two requests on one fit calls
`.hypothesis_plan_cache_clear()` and `.iwmde_qcmde_cache_clear()` first and
defers them on exit, so that no plan or normalization built under a mock is
read by another test.

## Correctness Evidence

- Test behavior, not implementation trivia.
- Use analytic identities or independent reference implementations for numerical
  kernels.
- Compare against metafor where estimands match. Save the metafor object in the
  cached fit's `info` during `test-01-*`, then reuse it in post-fit tests.
- Justify tolerances from Monte Carlo or numerical error. Do not use a fixed
  package-wide tolerance merely because it makes a test pass.
- State numerical error explicitly with absolute, relative, log, or peak-scaled
  comparisons where scale matters. Do not rely on testthat edition-specific
  `expect_equal()` normalization for small probabilities or near-zero quantities.
  Keep the current edition unless a separate migration validates its effects.
- Use deterministic call/path guards for algorithmic work limits. Machine-specific
  elapsed-time thresholds belong in focused performance verification, not ordinary
  unit assertions.
- A failing expectation requires diagnosis. Generate a candidate when the
  intended result changed; accept a verified baseline change only after
  maintainer or explicitly delegated review.
- Do not add redundant matrices, samples, fits, or assertions for coverage alone.
- Deterministic collinear replay draws establish structural contracts only;
  they do not provide Monte Carlo, calibration, or posterior-inference evidence.

The evidence inventory in `tests/testthat/REGRESSION-COVERAGE.md` records the
representative oracle and visual coverage expected for public post-fit methods.
Update it when a method changes tier or loses an oracle.

## Visual Regression

Use `expect_vdiffr_snapshot()` from
`tests/testthat/helper-visuals.R`. Representative, human-reviewed snapshots
belong in the standard profile. Gate redundant galleries with
`skip_if_not_full_visuals()`.

Retain candidates for maintainer or explicitly delegated review before
acceptance. Structural `as_data = TRUE` tests supplement visual snapshots;
they do not replace them.

After sourcing `tests/scenarios/helper-scenarios.R`, call
`review_test_snapshots()` to open testthat's native reviewer for ordinary text
and figure candidates under `tests/testthat/`, followed by per-artifact text and
table candidates under `tests/results/`. Its optional `files` argument has the
same test-name or trailing-slash directory semantics as
`testthat::snapshot_review()`; use `reference_filter` to select reference
groups such as `"interpret"` or `"marginal_means"`. Accept changes only after
maintainer or explicitly delegated review. Rejected `.new.txt` candidates are
deleted and skipped ones remain available for a later review.

Set `ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS=TRUE` only during an explicit snapshot
regeneration workflow.

## Final Verification

Run focused tests after each meaningful change. Before handoff, run the standard
profile when the change scope warrants it. State exactly what ran, what did not,
and why.

`tools/verify-ordered-sources.R <manifest.rds>` is focused deeper verification
for the original two-chain, 4000-row, twelve-study ordered fixtures. The RDS
contains a named character vector of existing fit paths (`normal`, `gamma`,
`selected`, `shared`). The manifest selects owner-managed existing fits wherever
they reside; `ROBMA_TEST_FILES_DIR` is printed context only and does not control
these paths. The script never fits, refreshes or writes a cache.
It compares source normalizers and supplier numerators against physical
original-kernel integrals, genuine selected/shared predictors against complete
replay, and twelve production N/D ratio arithmetic proofs against their reported
envelopes. On original rows 1/1990/4000 it also independently integrates both
physical block N and D, and compares their ratio using the complete original
likelihood on replayed T/S coordinates, public group alphas and declared total
priors. Shared allocations are varied once across their consumers. These
selective references carry estimated quadrature errors; they do not establish
global integration accuracy or chain calibration. It loads
the existing development DLL with `compile = FALSE`, records native/dependency
identity, and verifies that input fit hashes remain unchanged. This is
correctness evidence; it does not certify a new native build's performance.
The physical N/D/reference guards name the fixture, original row and quantity;
comparisons use `abs(expm1(logP - logR))` and the reference-scale envelope
`eP * (P/R) + eR + ef * max(1, P/R)`. Reported producer/reference errors are
relative to their own computed values; reference integrals are finite positive
and errors finite nonnegative. Every resulting comparison envelope must remain
at or below the existing 0.025 verification accuracy ceiling. Ratio arithmetic
identity checks remain separate. Set native/BLAS limits in the launching
environment before R starts; the verifier prints those limits and its native
thread option without claiming to reset already initialized BLAS state. Ordered short-row tests explicitly restrict an internal valid
plan population; the public planner's 20-continuous-row minimum stays unchanged.
Missing MCMC error with finite importance ESS remains partial uncertainty.
The ordered-result NA masks also reject malformed available-status all-NA
uncertainty. When the active branch has enough batches in each required fitted
chain, the complete conditioned indicator includes those same active rows.
Ordered fixed-positive singleton step-selection blocks may reduce normalizer
columns only after exact persisted data/design equality and bit-exact current
mean, variance and complete selection-input checks. The Gaussian law and
estimand, original row/power identities and full-observation error budgets are
preserved while their arithmetic is improved. Gaussian evaluation and Normal
completion use direct per-share precision and centered residuals from all
original observations. The actual
fixed-positive singleton product step context must match its declared weights,
bins and source routing before either normalization route uses the direct
observed factor. Source completion replays its Gaussian predictor and retains
the original nuisance-prior constant and both truncation masses. Stationary
polynomial roots, guards and Gamma completion retain their existing formulas.
Conditioned node tests preserve exact original IEEE inputs and compare against
independent declared-model constants at the original output-scale allowance.
Full replay is always related to that target through its actual mean/SD inputs
and independently certified declared-mean enclosures. Captured rounded-input
oracle comparisons additionally run when those replay inputs match; one-ulp
alternate inputs exercise the portable perturbation bound. All original nodes,
allowances and discriminating witnesses remain pinned. Actual unselected and
truncated Normal source contexts retain nuisance priors once and both prior
and posterior support masses. Their physical-integral and independent MP
comparisons include endpoint standardization, the same selected tail as the
shared production mass helper, exp/subtraction conditioning and reference-coordinate rounding.
The named finite-range libm allowance is engineering evidence, not a universal
pnorm accuracy theorem; the tight source case pins the complete allowance below
1e-9. Original conditioned-node allowances remain unchanged. Replay inputs also
obey a-priori 64-operation affine-mean and variance caps before the adaptive
conditioning relation; large deliberate perturbations fail those caps.
Normal source floating diagnostics retain interval conditioning separately
from residual arithmetic; far-tail source cases compare to independent physical
and MP/Arb references without accepting a rounded complementary probability.
Contradictory structural selection metadata
raises the density-method unavailable family; quantitative failures retain
their numerical diagnostics family. Row proofs retain the partition, source,
prior and control identities.
