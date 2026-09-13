# Maintainer Scenarios

Use this guide whenever creating, modifying, regenerating, or reviewing
maintainer scenarios under `tests/scenarios/` or their runner
`tools/test-scenario.R`.

In BayesToolsVerse, the shared validation guide owns scope, runtime criteria,
and the distinction from tests and verification. Use its configured R and
private agent library. The package helpers below own execution and review.

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
- Wrap costly non-fit computations used later in the scenario in
  `scenario_time("name", { ... })`. It returns the computed value unchanged and
  records timing only; it does not cache the value or create a snapshot.
- Put `add_loo()` and `add_marglik()` inside the cached block when the analysis
  needs those supported quantities, so expensive results are reused.
- Reuse cached fits while developing downstream functionality. Restarting R,
  reloading or reinstalling RoBMA, and rerunning only parts of an analysis are
  not reasons to refit.
- Cache validity is based on the stored fitting expression. Changing that
  expression refits automatically; changes to package internals do not. If an
  upstream change may alter the fitted posterior, treat the cache as
  potentially stale. During an authorized current-results refresh, refit the
  affected models when needed, retaining old caches and the reason for
  replacement. Otherwise consult the maintainer before regenerating them.
- `cache_version` is optional. Use or increment it only to invalidate one fit
  after an in-scope package-internal fitting change that leaves the fitting
  expression unchanged.
- Keep fit blocks compact. Scenario files may retain long model calls on one
  line when that is easier to scan on the maintainer's ultrawide monitor.

Fit caches under `tests/scenarios/cache/` are local and ignored. Never replace
them with ordinary package-test caches or duplicate the cache logic in a
scenario.

Each successful uncached fitting block records timing metadata excluding cache
serialization. Calls to unqualified `add_loo()` and `add_marglik()` inside the
block are timed automatically and stored separately from the remaining
model-fitting time; scenario files do not need timing wrappers. A cache hit
contributes no timing measurement: do not compare its historical production
timing with the baseline, update the baseline from it, or substitute the
cache-loading time. Timing completeness checks ignore fit rows served from
cache during the current run.

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
- Keep `ex_m()` and `ex_r()` parameter and statistic vectors on one line, even
  when long. Parameter names become comparison-table labels: omit names that
  merely repeat their selectors, and use names only for different labels or
  repeated selectors that need distinct labels.
- Use `ex_p()` for the pooled-effect column in comparisons with
  `metafor::predict()`.
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

Call `scenario_start()` before `scenario_fit()`, `scenario_time()`,
`scenario_text()`, or `scenario_plot()`. Use the shared helpers rather than
recreating their caching, timing, printing, graphics, or snapshot behavior.

Run one scenario with:

```text
Rscript tools/test-scenario.R <name>
```

Use `--refit` to refit models and compare them against the locked outputs. Use
`--update` to reuse valid fit caches and replace exercised outputs. The combined
`--regenerate` shortcut performs both actions. Use `--update-timings` only to
explicitly accept current timings, including slower timings.

After sourcing `helper-scenarios.R`, `test_scenarios()` runs all scenarios from
an interactive session; the existing singular spelling `test_scenario()` is
retained as an alias. Pass a regular-expression `filter` to select scenario
names, for example `test_scenarios(filter = "assink|bcg")`. By default, this
reuses fit caches, suppresses artifact output, and compares without replacing
text or plot baselines. Missing and faster timing measurements are maintained
automatically. All selected scenarios share one testthat reporter lifecycle, so
interactive runs print one combined progress report and one results block. Set
`refit = TRUE`, `update = TRUE`, or `regenerate = TRUE` explicitly to change fit
or output behavior. Set `update_timings = TRUE` only to accept current timings
even when they are slower.

Call `review_scenario_snapshots()` after a run to review all cached table and
figure candidates from the scenarios selected by the latest
`test_scenario()`. An optional regular-expression `filter` further narrows that
selection. The function uses testthat's one-by-one Accept, Reject, and Skip
reviewer; skipped candidates remain available for another call.

Sourcing a scenario file directly in an interactive session defaults to reusing
fit caches, showing artifact output, creating missing text and SVG baselines,
and comparing existing baselines without replacing them. Pass `update = TRUE`
to `scenario_start()` explicitly when direct execution should replace exercised
outputs. Managed and non-interactive runs neither create nor replace baselines
without an explicit update. Do not add control variables to scenario files.

Text baselines under `results/` and SVGs under `_snaps/` are committed
human-reviewed artifacts. Missing baselines fail during `test_scenario()`, CLI,
and non-interactive runs unless updating is explicitly enabled. With updating
disabled, changed text is retained as `<name>.new.txt` and reported with
testthat's colored value diff. Interactive `test_scenario()` runs finish all
selected scenarios before listing the table and figure mismatches and opening
the same combined reviewer. Accepted candidates replace the baseline, rejected
candidates are removed, and skipped candidates remain cached. Direct plot
comparison retains a changed `<name>.new.svg`. Accepting an already verified
candidate does not itself require rerunning an expensive expression. Rerun
when implementation, inputs, or RNG behavior changed, or reproducibility is
unverified. Run deferred review from the runner's exit
path so an early scenario failure cannot bypass cached candidates.
Non-interactive runs never prompt or accept changes.

Scenario runs report committed text and SVG artifacts that are no longer
referenced by the scenario. Never delete an orphan automatically.

## Performance Timings

`scenario_fit()`, `scenario_time()`, `scenario_text()`, and `scenario_plot()`
record wall time and peak R-managed memory in GB in the tracked
`timings/<scenario>.tsv` baseline. Memory uses R's resettable
garbage-collector peak and excludes external-process memory. A `time` entry
includes only successful expression evaluation. Text timing and memory include
expression evaluation and captured printing. Plot timing and memory include the
canonical assertion render but exclude any extra interactive preview and file
comparison. The `fit` timing and memory include the full cached block but
exclude cache reads and writes.
When the block directly calls unqualified `add_loo()` or `add_marglik()`, the
same measurement also records `fit_model`, `fit_loo`, and `fit_marglik` rows.
`fit_model` is the full block time remaining after the observed post-fit calls.
The redundant total `fit` row remains useful as an end-to-end regression check
but is excluded from the unweighted average when split rows are available.
Nested phase rows do not record separate memory peaks; the enclosing `fit` row
owns the peak-memory measurement.

After a complete managed scenario file, issue one warning that lists every call
whose wall time increased by more than 20% and the unweighted mean percentage
change across calls when it is a regression greater than 5%. Start each per-call
warning with its absolute change, followed by the percentage and old and new wall
times, then the timed-part name. Highlight absolute increases greater than two
seconds in red. Display seconds to one decimal place and percentages as whole
numbers. Only fits evaluated during the current run contribute fit timings;
cached fits contribute none. Text and plot calls use fresh times. Compare
against the old timing baseline before replacing any
timing, fit-cache, text, or plot artifact. Retain calls measured below 0.75
seconds in timing baselines and candidates, but exclude them from both per-call
warnings and the mean regression assessment.

For memory, warn for a call when its peak exceeds 2 GB and is more than 20%
above its baseline. Warn whenever a call exceeds 8 GB even without a baseline
or relative increase. Retain the lowest observed elapsed time and memory peak
independently; `update_timings = TRUE` explicitly accepts both current values.
Timing files without the `memory_gb` column are read with unavailable memory
and backfilled as calls are next evaluated.

Automatically add every available measurement whose timing row is absent from
the baseline. Automatically replace an existing elapsed time or memory peak
only when that metric improves. Apply this maintenance in managed runs and as
successful calls are executed line by line during direct interactive
development. Do not infer or replay fit performance on a cache hit.
Intentionally refit a model whenever its performance should be measured and
tested.

Retain slower current measurements in the ignored `<scenario>.new.tsv` and
leave the faster baseline unchanged. Only explicit `update_timings = TRUE` or
`--update-timings` may accept slower measurements, such as after intentionally
tightening precision. Output updates, regeneration, and refitting do not imply
consent to a slower baseline. Compare against the old baseline before applying
automatic or explicit updates. Record R version and platform as informational
provenance without automatically invalidating the baseline.

Use `plot_scenario_times()` to draw horizontal elapsed-time boxplots from the
committed baselines for one exact scenario or across scenarios. It classifies
individually named calls conservatively, omitting aggregate `fit` rows,
`.new.tsv` candidates, and unmatched compound timings. Selection-model fit
boxes use explicit cond/marg model labels for conditioning specifications. Its `unit` argument sets
seconds, minutes, or hours on the x-axis without changing stored values.

Use `plot_marginal_diagnostics()` for metafor-to-RoBMA or RoBMA-to-RoBMA
diagnostic comparisons. Both RoBMA inputs must use explicit marginal targets;
label two RoBMA fits with `reference_label` and `estimate_label` when their
distinction matters. When either fit inherits from `bselmodel`, this shared
helper compares raw residuals and DFBETAS. Rstandard, hat values, and Cook's
distance remain structurally unavailable for weight-function selection models.
Dedicated selected-model comparisons may also use `rstudent(fit)[["z"]]` for
LOO-PIT with stored estimate-unit LOO. This deletion-conditioned target retains
the original publication event; do not pass `conditioning_depth` to `rstudent()`.

Do not use output updating or regeneration merely to make a scenario pass.
Generate candidates within the requested scenario work. Accept baseline
changes only after maintainer or explicitly delegated review, and review every
resulting diff. Requests to generate or refresh output alone do not delegate
review.

## Discrepancies and Scope

An explicit request to create or refresh scenario snapshots authorizes that
workflow. Generate missing candidates one by one, inspect each result and its
runtime, and investigate obvious errors or suspiciously slow computation as
they arise. Fix clear in-scope implementation errors without asking again for
permission. Validate changed computations with the original workload and RNG
state. Keep candidates pending for maintainer review unless that review was
explicitly delegated.

Preserve the original analysis, model specification, data, comparison methods,
evaluation order, seeds, and parallel settings. Change scenario controls only
within explicit authorization, including authorization already given. When
accuracy requires a different budget, report the proposed public control and
its justification separately from matched-workload performance evidence. Do
not hide budget changes in package code, lower diagnostic requirements, or add
parallel computation merely to improve timings.

Do not suppress or normalize away a discrepancy, or change an expected
relationship merely to make a comparison pass. Investigate against an
independent reference. If a statistical or architectural choice remains
ambiguous, record the artifact, observed pattern, reproduction command, impact,
alternatives, and recommendation in the workspace's
[pending decisions](../../../.agents/decisions.md)
(or a temporary package decision note in a standalone checkout), then request
the needed decision while continuing independent authorized work.
Unrequested scenario expansion, reformatting, and unrelated cleanup remain out
of scope.
