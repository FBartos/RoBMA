# Maintainer scenarios

Scenario files are deliberately excluded from routine package tests. Name each
file `test-<scenario>.R` and run it with:

```text
Rscript tools/test-scenario.R <scenario>
```

Use `--refit` to refit wrapped models while comparing against the locked
outputs. Use `--update` to reuse valid fit caches and replace exercised text and
plot baselines. `--regenerate` does both. Use `--update-timings` to explicitly
accept current timings, including slower timings. Review all resulting Git
diffs.
Use `scenario_agreement_plot()` for reference-versus-RoBMA difference plots so
their finite-value filtering, agreement band, and axes stay consistent. Set
`reference_label` and `estimate_label` for comparisons that are not
metafor-versus-RoBMA.

The project `.Rprofile` loads `test_scenarios()` for interactive sessions. It
reuses fit caches, suppresses artifact output, and compares without replacing
text or plot baselines. Missing and faster timing measurements are maintained
automatically. Its `filter` is a regular expression matched against scenario
names; `test_scenario()` remains an alias:

```r
test_scenarios()
test_scenarios(filter = "assink")
test_scenarios(filter = "assink|bcg")
test_scenarios(filter = "assink", refit = TRUE, update_timings = TRUE)
review_scenario_snapshots()
review_test_snapshots()
```

`review_scenario_snapshots()` reopens all cached table and figure changes from
the scenarios selected by the latest `test_scenarios()` call. Its optional
`filter` further narrows those scenario names. `review_test_snapshots()` opens
testthat's native reviewer for ordinary snapshots under `tests/testthat/` and
then reviews per-artifact text and table candidates under `tests/results/`.
Pass testthat's `files` selection for ordinary snapshots and
`reference_filter` for reference groups such as `"interpret"` or
`"marginal_means"`.

Minimal scenario:

```r
scenario_start("bcg")

testthat::test_that("BCG analysis remains stable", {

  fit <- scenario_fit("meta-analysis", {
    brma(y = yi, se = sqrt(vi), data = dat)
  })

  scenario_text("summary", {
    summary(fit)
  })

  scenario_plot("forest", {
    forest(fit)
  })
})
```

Sourcing a scenario file directly in an interactive session also reuses fit
caches, but shows artifact output and immediately creates or replaces its text
and SVG baselines. No control variables are needed in scenario files.

Fit caches under `cache/` are local and ignored. Text files under `results/`
and SVG files under `_snaps/` are regression baselines. Per-scenario TSV files
under `timings/` store wall times and peak R-managed memory in GB for every
`scenario_fit()`, `scenario_time()`, `scenario_text()`, and `scenario_plot()`
call. All three baseline types must be committed. Missing text and plot
baselines follow the snapshot workflow; missing performance measurements are
backfilled automatically by direct and managed execution.
Scenario and artifact names may use ASCII letters in either case, numbers,
underscores, hyphens, and internal periods.

Each fit cache records its normalized fitting call and refits automatically when
that call changes. Pass an optional positive `cache_version` to
`scenario_fit()` only when a package-internal fitting change should invalidate
one otherwise unchanged call. A fit cache also records the wall time, peak
R-managed memory, and runtime provenance of the fitting block. Cache hits do not
contribute historical performance measurements or RDS loading performance;
refit a model when its fit performance should be assessed.

`scenario_text()` automatically prints a visible returned value into its locked
output, so summary and table calls do not need an explicit `print()`. Message
conditions emitted while evaluating the expression, such as fitting progress,
are suppressed; warnings and errors remain visible. Changed text is stored as
`<name>.new.txt` and reported with testthat's colored value diff. An interactive
`test_scenarios()` run collects all table and figure mismatches, finishes the
selected scenarios, and then opens testthat's snapshot reviewer for one-by-one
visual Accept, Reject, or Skip decisions. Accepted candidates replace the
baseline, rejected candidates are removed, and skipped candidates remain
cached. Call `review_scenario_snapshots()` to reopen the cached candidates, and
rerun `test_scenarios()` after review to confirm accepted changes. The reviewer
is also opened when a selected scenario stops early after caching a mismatch.
Non-interactive runs retain candidates and fail without prompting.

Use `test_scenarios(update = TRUE)` or CLI `--update` only when a managed run
should replace baselines. `refit = TRUE` / `--refit` refits models without
updating outputs, and `regenerate = TRUE` / `--regenerate` does both. With
updating disabled, changed text is retained as `<name>.new.txt` and direct plot
comparison retains `<name>.new.svg`. `scenario_text()` and `scenario_plot()`
each reset the random seed to 1 before evaluating an artifact.

Timing comparison is deferred until the complete scenario file finishes. One
warning reports every call more than 20% slower than its baseline and an
unweighted mean percentage regression across calls greater than 5%. Per-call
warnings start with the absolute change, followed by the percentage and old and
new wall times, then the timed-part name. Increases greater than two seconds are
highlighted in red. Seconds are shown to one decimal place and percentages as
whole numbers. Calls measured below 0.75 seconds are retained but excluded from
regression warnings. Every successful run automatically adds missing timing
rows and replaces each existing time or memory metric only when it improves.
This also happens while a scenario is executed line by line in an interactive
session.

The `memory_gb` column records the peak reported by R's resettable garbage
collector and excludes external-process memory. A call warns when it exceeds
2 GB and is more than 20% above its baseline, and whenever it exceeds 8 GB
regardless of the baseline. The lowest elapsed time and memory peak are retained
independently. Nested `fit_model`, `fit_loo`, and `fit_marglik` rows have no
separate memory peak; their enclosing `fit` row captures the complete fit.
Existing timing files without `memory_gb` are backfilled as calls are evaluated.

Slower or higher-memory measurements never replace the corresponding baseline
metric automatically. They remain in the ignored
`timings/<scenario>.new.tsv` candidate. Set
`update_timings = TRUE` or use `--update-timings` only to explicitly accept the
current measurements, including an intentional slowdown or memory increase.
Output updates, regeneration, and refitting do not themselves grant that
consent. The old performance baseline is always compared before any missing,
improved, or explicitly accepted measurement is written. The TSV records R
version and platform as review provenance but does not automatically invalidate
measurements.

Base graphics parameters are restored after each drawn evaluation, so settings
such as `par(mar = ...)` do not leak into subsequent plots. Base graphics
overlay mode is reset before and after each drawn block, preventing one scenario
figure from being drawn over the previous figure. Scenario runs report locked
text or SVG artifacts that are no longer referenced by the scenario; they never
delete them automatically.
