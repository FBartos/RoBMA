# Maintainer scenarios

Scenario files are deliberately excluded from routine package tests. Name each
file `test-<scenario>.R` and run it with:

```text
Rscript tools/test-scenario.R <scenario>
```

Use `--refit` to refit wrapped models while comparing against the locked
outputs. Use `--update` to reuse valid fit caches and replace exercised text and
plot, and timing baselines. `--regenerate` does both. Use `--update-timings` to
replace only the timing baseline. Review all resulting Git diffs.
Use `scenario_agreement_plot()` for reference-versus-RoBMA difference plots so
their finite-value filtering, agreement band, and axes stay consistent. Set
`reference_label` and `estimate_label` for comparisons that are not
metafor-versus-RoBMA.

From an interactive session, source `helper-scenarios.R` and run all scenarios
with `test_scenario()`. This reuses fit caches, suppresses artifact output, and
compares without replacing locked baselines. Its `filter` is a regular
expression matched against scenario names:

```r
source("tests/scenarios/helper-scenarios.R")
test_scenario()
test_scenario(filter = "assink")
test_scenario(filter = "assink|bcg")
test_scenario(filter = "assink", refit = TRUE, update_timings = TRUE)
```

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
under `timings/` store wall times for every `scenario_fit()`, `scenario_text()`,
and `scenario_plot()` call. All three baseline types must be committed. Missing
baselines are added by direct interactive execution and rejected or warned
about by managed comparison runs.
Scenario and artifact names may use ASCII letters in either case, numbers,
underscores, hyphens, and internal periods.

Each fit cache records its normalized fitting call and refits automatically when
that call changes. Pass an optional positive `cache_version` to
`scenario_fit()` only when a package-internal fitting change should invalidate
one otherwise unchanged call. A fit cache also records the wall time and runtime
provenance of the fitting block. A cached fit contributes that original fitting
time, not its RDS loading time, so cached and refitted scenario totals remain
comparable. A cache created before timing metadata existed must be refitted once
before its fitting time is available.

`scenario_text()` automatically prints a visible returned value into its locked
output, so summary and table calls do not need an explicit `print()`. Message
conditions emitted while evaluating the expression, such as fitting progress,
are suppressed; warnings and errors remain visible. Changed text is stored as
`<name>.new.txt` and reported with testthat's colored value diff. An interactive
`test_scenario()` run collects all such mismatches, finishes the selected
scenarios, and then opens testthat's snapshot reviewer for one-by-one visual
Accept, Reject, or Skip decisions. Accepted candidates replace the baseline,
rejected candidates are removed, and skipped candidates remain cached. Rerun
`test_scenario()` after review to confirm accepted changes. The reviewer is also
opened when a selected scenario stops early after caching a mismatch.
Non-interactive runs retain candidates and fail without prompting.

Use `test_scenario(update = TRUE)` or CLI `--update` only when a managed run
should replace baselines. `refit = TRUE` / `--refit` refits models without
updating outputs, and `regenerate = TRUE` / `--regenerate` does both. With
updating disabled, changed text is retained as `<name>.new.txt` and direct plot
comparison retains `<name>.new.svg`. `scenario_text()` and `scenario_plot()`
each reset the random seed to 1 before evaluating an artifact.

Timing comparison is deferred until the complete scenario file finishes. One
warning reports every call more than 20% slower than its baseline and an
unweighted mean percentage regression across calls greater than 5%. A
compare-only run writes the ignored candidate
`timings/<scenario>.new.tsv`; it never changes the committed baseline. Output
updates also update timings, while `update_timings = TRUE` / `--update-timings`
updates timings alone. The old timing baseline is always compared first, so a
regression is reported even when the fit cache, output snapshot, and timing
baseline are replaced in the same run. An accepted update stores the measured
time, whether it is faster after an optimization or intentionally slower after
increasing precision; it is not an all-time minimum. `refit = TRUE` alone is
compare-only. The TSV records R version and platform as review provenance but
does not automatically invalidate timings.

Base graphics parameters are restored after each drawn evaluation, so settings
such as `par(mar = ...)` do not leak into subsequent plots. Base graphics
overlay mode is reset before and after each drawn block, preventing one scenario
figure from being drawn over the previous figure. Scenario runs report locked
text or SVG artifacts that are no longer referenced by the scenario; they never
delete them automatically.
