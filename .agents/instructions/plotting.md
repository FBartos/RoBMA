# Plotting

Use this guide for plot APIs, plot data, and base/ggplot renderers.

## Architecture

Follow the architecture of the plot family being changed. For plots supporting
both base graphics and ggplot2, keep statistical computation separate from
rendering:

1. public method validates and dispatches;
2. a data function computes backend-neutral plot data;
3. base and ggplot renderers consume that data.

Do not introduce this layering for a trivial plot that has only one backend, and
do not force unrelated plot families into one abstraction.

When supported, `as_data = TRUE` returns the computed data without drawing.
Both renderers must use the same computed values. Base methods return
`invisible(NULL)`; ggplot methods return the plot object.

Normalize defaults shared by both backends once, in the dispatcher or a focused
options helper. Render semantic layers in the same order in both backends;
backend-specific styling must not change the represented quantities.

Use `.outcome_data_yi()`, `.outcome_data_sei()`,
`.outcome_data_vi()`, and `.outcome_data_weights()` instead of direct fitted
object access.

## Publication-Bias Semantics

User-facing `sampling_bias` means whether the displayed sampling distribution
includes publication bias:

- `sampling_bias = TRUE`: show the distribution users would observe;
- `sampling_bias = FALSE`: show the bias-adjusted distribution.

When delegating to `predict.brma()`, this maps to
`bias_adjusted = !sampling_bias`. Preserve this semantic in funnel and
regression plots. Selection branches must use `.selection_context()`.

## API and Verification

Use `plot_type = "base"` or `"ggplot"` where both backends are available.
Keep released plotting arguments working through the compatibility policy in
root `AGENTS.md`.

Test computed plot data and dispatch separately, then retain representative
human-reviewed visual snapshots for both backends. Structural or render-only
checks do not replace visual regression. See `testing.md`.

Current split implementations live in `R/funnel-*.R`, `R/regplot-*.R`, and
`R/zplot-*.R`. Inspect the current family before adding files or helpers.
