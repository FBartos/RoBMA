# Scenario Tests

Use this guide whenever creating, modifying, regenerating, or reviewing
maintainer scenarios under `tests/scenarios/` or their runner
`tools/test-scenario.R`.

## Purpose

Scenarios are maintainer-facing analysis scripts, not ordinary unit or
validation tests. They compare user-visible RoBMA results with an independent
`metafor` result when the estimands match, or compare RoBMA models and
computations that should agree or differ in a known way. Their text and visual
baselines preserve output after human review and make later regressions visible.

Keep the analysis direct and readable. Do not add direct
`testthat::expect_*()` assertions to scenario files; quantitative assertions
belong in `tests/testthat/`. Do not add multi-step data-wrangling code merely to
assemble a comparison table. A compact `data.frame()` or `cbind.data.frame()`
inside `scenario_text()` is encouraged when it makes several expected
relationships easy to inspect.

## Cached Fits

- Wrap each costly RoBMA fit in `scenario_fit()` with a stable, descriptive
  name. Return a ready-to-analyze object from the block.
- Add `add_loo()` and `add_marglik()` inside the cached block whenever they are
  supported and potentially useful for validation or face-validity checks.
- Reuse cached fits while developing downstream functionality. Restarting R,
  reloading or reinstalling RoBMA, and rerunning only parts of an analysis are
  not reasons to refit.
- Cache validity is based on the stored fitting expression. Changing that
  expression refits automatically; changes to package internals do not. If an
  upstream change may alter the fitted posterior, treat the cache as
  potentially stale and consult the maintainer before regenerating it.
- `cache_version` is optional. Use or increment it only to invalidate one fit
  after an approved package-internal fitting change that leaves the fitting
  expression unchanged.
- Keep fit blocks compact. Scenario files may retain long model calls on one
  line when that is easier to scan on the maintainer's ultrawide monitor.

Fit caches under `tests/scenarios/cache/` are local and ignored. Never replace
them with ordinary package-test caches or duplicate the cache logic in a
scenario.

Each fit cache has companion timing metadata for the successful fitting block,
excluding cache serialization. A cache hit contributes that stored production
time rather than RDS loading time. If timing metadata is unavailable, refit the
cache before accepting a timing baseline; never substitute the cache-loading
time.

## Human-Focused Snapshots

- Use `scenario_text()` for output a user or maintainer would read: summaries,
  estimates, hypotheses, model-fit measures, and compact comparison tables.
- Keep the relevant `metafor` result or alternative RoBMA computation adjacent
  to the locked RoBMA output so the comparison is easy to inspect during
  interactive development.
- Put expected equivalences or differences in short comments beside the values
  being compared. These comments record the intended interpretation; do not
  replace them with direct assertions in the scenario.
- Use `scenario_plot()` for user-facing plots, density-estimator overlays,
  comparisons between simpler and richer models, and vector-valued results
  where visual patterns are more informative than scalar checks.
- Prefer `scenario_agreement_plot()` for pointwise `metafor`-versus-RoBMA
  comparisons such as predictions, BLUPs, random effects, residuals, and
  diagnostics. Its difference scale and reference-relative band are part of the
  intended human review. Set `reference_label` and `estimate_label` when
  comparing two RoBMA results or another reference pair.
- A small local helper is appropriate when it removes repeated plot setup or
  repeated extraction for a family of visual comparisons. Do not introduce a
  general abstraction for one scenario or use helpers to hide the quantities
  being compared.
- `scenario_text()` and `scenario_plot()` set the random seed to 1 internally.
  Do not add artifact-specific seeds unless the maintainer requests a different
  stochastic comparison.
- `scenario_text()` suppresses message conditions emitted while evaluating its
  expression, including fitting progress. Warnings and errors remain visible.

Ordinary package line-width limits do not apply to scenario scripts. Preserve
long calls, aligned comparison columns, compact one-line plot invocations, and
inline relationship comments when they improve navigation and human review. Do
not mechanically reformat an existing scenario to ordinary source-file width.

## Files and Workflow

Call `scenario_start()` before `scenario_fit()`, `scenario_text()`, or
`scenario_plot()`. Use the shared helpers rather than recreating their caching,
printing, graphics, or snapshot behavior.

Run one scenario with:

```text
Rscript tools/test-scenario.R <name>
```

Use `--refit` to refit models and compare them against the locked outputs. Use
`--update` to reuse valid fit caches and replace exercised outputs. The combined
`--regenerate` shortcut performs both actions. Output updates also replace the
timing baseline. Use `--update-timings` when only timings should be replaced.

After sourcing `helper-scenarios.R`, `test_scenario()` runs all scenarios from
an interactive session. Pass a regular-expression `filter` to select scenario
names, for example `test_scenario(filter = "assink|bcg")`. By default, this
reuses fit caches, suppresses artifact output, and compares without replacing
locked baselines. Set `refit = TRUE`, `update = TRUE`, or
`regenerate = TRUE` explicitly to change those behaviors. Set
`update_timings = TRUE` to replace only timing baselines.

Sourcing a scenario file directly in an interactive session defaults to reusing
fit caches, showing artifact output, and immediately creating or replacing text
and SVG baselines. Do not add control variables to scenario files.

Text baselines under `results/` and SVGs under `_snaps/` are committed
human-reviewed artifacts. Missing baselines fail during `test_scenario()`, CLI,
and non-interactive runs unless updating is explicitly enabled. With updating
disabled, changed text is retained as `<name>.new.txt` and reported with
testthat's colored value diff. Interactive `test_scenario()` runs finish all
selected scenarios before listing the text mismatches and opening testthat's
snapshot reviewer for one-by-one Accept, Reject, or Skip decisions. Accepted
candidates replace the baseline, rejected candidates are removed, and skipped
candidates remain cached. Direct plot comparison retains a changed
`<name>.new.svg`. Rerun `test_scenario()` after review to confirm accepted
changes. Run deferred review from the runner's exit path so an early scenario
failure cannot bypass cached candidates. Non-interactive runs never prompt or
accept changes.

Scenario runs report committed text and SVG artifacts that are no longer
referenced by the scenario. Never delete an orphan automatically.

## Performance Timings

`scenario_fit()`, `scenario_text()`, and `scenario_plot()` record wall time in
the tracked `timings/<scenario>.tsv` baseline. Text timing includes expression
evaluation and captured printing. Plot timing includes the canonical assertion
render but excludes any extra interactive preview and file comparison. Fit
timing includes the full cached block, such as `add_loo()` and `add_marglik()`,
but excludes cache reads and writes.

After a complete managed scenario file, issue one warning that lists every call
whose wall time increased by more than 20% and the unweighted mean percentage
change across calls when it is a regression greater than 5%. Cached fits use
their stored production times while text and plot calls use fresh times. Compare
against the old timing baseline before replacing any timing, fit-cache, text, or
plot artifact.

Compare-only runs retain the current measurements in the ignored
`<scenario>.new.tsv`. Direct interactive execution, output updates,
regeneration, and explicit timing updates replace the baseline after comparison.
`refit` alone remains compare-only. An update accepts the current measured time,
including improvements and intentional slowdowns caused by a changed scenario
specification; do not retain an automatic all-time minimum. Record R version and
platform as informational provenance without automatically invalidating the
baseline.

Do not use output updating or regeneration merely to make a scenario pass.
Replace baselines only when the maintainer explicitly requests, approves, or
interactively accepts the change, and review every resulting diff.

## Discrepancies and Scope

If creating or modifying a scenario reveals an unexplained discrepancy, changed
snapshot, unexpected agreement, or unexpected difference, stop and notify the
maintainer before changing anything to accommodate it. Report the scenario and
artifact, the compared models or methods, the observed pattern, and the command
used. Ask how to proceed.

Do not ad hoc change package behavior, scenario inputs, comparison methods,
cached fits, snapshot baselines, or expected relationships. Do not suppress or
normalize away a discrepancy. When a scenario file is already in progress,
make only the requested change; ask before cleaning up, reformatting, expanding,
or correcting adjacent work.
