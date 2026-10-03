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

For ordered semantic targets, derive owners from selected catalog quantities'
`term` and `formula_parameter`, including group members, with
`JAGS_parameter_names()`. Intersect them with declared ordered priors before
requesting only those owners from `JAGS_ordered_parameter_spec()`; ordinary
targets do not request ordered metadata. Compare persisted catalog dependencies
with the selected specs' coefficient names. Leave
unspecified whole-factor prior line colors and types to BayesTools' matched
level styles. For selected ordered cells, resolve the scalar posterior color
(`col`, then `color`, then `colour`) and use it for an unspecified prior color
with a dashed line. Explicit prior colors and types, including aliases, win;
`linetype` takes precedence over `lty` in both backends. Preserve ordinary
targets' existing defaults.

For an ordered factor cell, retain the semantic marginal selected by its exact
catalog extraction weights, including its prior and aligned draw metadata.
Displayed contrast columns alone are not primitive ordered source draws; do not
re-evaluate that marginal through `parameter_draws()` on those columns.

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

## Predictive Targets

Use the prediction notation and two-axis contract in
[multivariate-targets.md](multivariate-targets.md).
Regression/forest prediction intervals, normal-model outcome-mode funnel
contours, and default z-plots target marginal new effects or responses.
GLMM funnel contours and regression sampling bands are explicitly descriptive
normal effect-size approximations, not discrete response predictions.
Plug-in funnel heterogeneity uses the within-model RMS SD, with integrated and
conditioned components averaged as variances. Location/bias parameters use their
within-model posterior means. Posterior-predictive funnels still average CDFs
over draws before inversion.
Residual-mode plots inherit the explicitly requested residual conditioning
depth; LOO-PIT remains the deletion-conditioned LOO target. Do not infer
conditioning from whether a plot evaluates fitted rows or an explicit grid.

Selected zplots may reuse checked point-context projections at retained-context
QMC nodes. Preserve the full event at each node, the outer node sequence and
the existing fallback for unsupported covariance geometry. Carry the inner
absolute density and normalization errors through the positive outer weights;
include their effects on MCSE and refinement changes before checking the outer
curve's peak. Standalone point contexts retain their relative-error requirement.

## Credible Bands

Pointwise credible bands of a curve over a grid (the continuous `regplot()`
band and the `zplot()` density bands) are Harrell-Davis quantiles of the
posterior draws, `BayesTools::harrell_davis_quantile()`, so that they are
smooth along the grid. Scalar intervals (categorical `regplot()` levels),
prediction and sampling intervals (posterior predictive mixture quantiles), and
step bands keep their own estimators; document the estimator where a band is
described.

## API and Verification

Use `plot_type = "base"` or `"ggplot"` where both backends are available.
Keep released plotting arguments working for the compatibility periods in
[model-interface.md](model-interface.md).

Test computed plot data and dispatch separately, then retain representative
human-reviewed visual snapshots for both backends. Structural or render-only
checks do not replace visual regression. See `testing.md`.

Current split implementations live in `R/funnel-*.R`, `R/regplot-*.R`, and
`R/zplot-*.R`. Inspect the current family before adding files or helpers.
