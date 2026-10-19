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
level styles. For selected ordered cells with `prior = TRUE`, resolve the scalar posterior color
(`col`, then `color`, then `colour`) and use it for an unspecified prior color
with a dashed line. In this scope, normalize prior color aliases with priority
`col`, `color`, then `colour`. Explicit prior colors and types, including aliases, win;
`linetype` takes precedence over `lty` in both backends. Preserve ordinary
targets' existing defaults.

For an ordered factor cell, retain the semantic marginal selected by its exact
catalog extraction weights, including its prior and aligned draw metadata.
Displayed contrast columns alone are not primitive ordered source draws; do not
re-evaluate that marginal through `parameter_draws()` on those columns.
At an actual scalar backend-source to display-target mapping, label a stored
qCMDE/IWMDE curve with the declared BayesTools plotted leaf and retain its source
parameter as provenance. Explicit precomputed plots require a matching selected
continuous curve; unavailable or wrong-target metadata cannot select KDE.

Ordered qCMDE/IWMDE inference uses the public primitive projection and source
recipes. Both requested methods normalize the same complete conditional source
kernel, with their original row budgets. Supplier/population classification
precedes chart selection and uses declared support and prior provenance.
Fixed and nonzero Bernoulli totals distinguish regular signed interior
contractions before endpoint exponents; uncertified interior limits remain typed
capabilities and a genuinely zero total retains its atom.
Admitted share-supplier endpoints use the allocation group of every nonzero
coefficient, even when their values differ, retaining within-group ratios and
the exact matching declared affine offset. Signed interiors keep their own chart.
An ordered constant-contribution certificate requires one declared total source
slice as well as the existing fixed-allocation/no-other-randomness guards.
Multiple slices retain ordinary uncertainty, including unavailable empirical
zero-variance diagnostics; equal observed contributions do not certify constancy.
Regular offsets with finite first but uncertified fourth reciprocal moments
use a positive block. Normalize source and block rows separately, replace only
the exact offset rows in the common contribution matrix, then aggregate and
compute chain/FPC uncertainty once. Preserve the strict sample alignment gate,
exact offset bits, per-row 0.025 numerical gate and separate atom measure.
Unknown suppliers or missing population precision certificates remain explicit
point capabilities; density curves omit only classified exact unavailable
points and retain their reasons. Keep log/scaled contributions under underflow.
Ordinary IWMDE half-bounded weights fit natural distance on declared row
supports and normalize compact Normal intervals by their actual endpoints.
Transient row-aligned zero certificates admit only valid out-of-proposal-support
draws; preserve their original denominator, chain and conditioning population.
Finite closed-support omega conditioners retain raw endpoints. These fitted
proposals do not certify global second moments or proposal-learning uncertainty.
Lattice qCMDE normalization caches row-resolution estimates from the nested-grid
maximum: the Gaussian physical width divided by that point's chart Jacobian,
or the existing log-density second difference when no kernel is available.
The component is `e / (1 - e)`, with `e = 2 exp(-2 pi^2 (s_z / h_z)^2)` at
the nested spacing, combined by actual ordinate shares. Unknown widths remain
unavailable; nonlattice normalizers contribute zero. The combined stability
metric is an error estimate, not a rigorous bound. Ordered source ratio errors
include this component in their denominator proof and report it without adding
it a second time to the stability gate.
Block capabilities must match the invariant diagonal Gaussian geometry and its
literal Normal/Gamma prior scope. Finite-upper Gamma share inverse suppliers
with no validated route are classified before denominator work. Unselected,
source-state-free blocks may complete literal untruncated Gamma shape 1/2
masses only within the score/conditioning and prospective arithmetic reserve;
retain effective-prior, method and fallback evidence in the row counts. Rejected
completion scales never enter the existing CDF floating-error ledger. Selected
blocks retain the original native selection correction at every source node.
Apply block prior-family eligibility only to resolved block components; raw,
zero and atomic components keep their existing routes. Missing MCMC uncertainty
at some ordered evaluations retains the original selected-row scope with partial
status and reason. Preserve available per-value metrics, reject incomplete
curve availability, and keep affected point precision flags false. Distinguish
these missing metrics from finite log MCSE with absolute representation underflow.

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
Regplot's scalar-intercept pointwise replay fallback matches only
`BayesTools_prediction_new_rows_unavailable` with reason
`posterior_indexed_row_source`. Mapping, source, monitoring and callback errors
stay strict. Its covariance adapter supplies no implicit fitted-row mapping;
an evaluable scale-source callback may support a new grid directly.
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

Ordinary qCMDE/IWMDE rows compose the existing numerical ordered-prior kernel
with the remaining ordinary priors. Cache identities preserve exact inclusion
indicators, including zero; selected components retain their strict validation.
Dependent vector priors use their full joint density and leave unused scalar
focal baselines unavailable. Ordered expression nuisance terms retain the
existing typed ordered-expression refusal for the requested target.

Use `plot_type = "base"` or `"ggplot"` where both backends are available.
Keep released plotting arguments working for the compatibility periods in
[model-interface.md](model-interface.md).
Posterior plots pass explicit limits in display coordinates with
`transformation_settings = TRUE`; prior-facing plots keep their source-coordinate
contract with `FALSE`. Automatic prior windows remain BayesTools-owned.

Test computed plot data and dispatch separately, then retain representative
human-reviewed visual snapshots for both backends. Structural or render-only
checks do not replace visual regression. See `testing.md`.

Current split implementations live in `R/funnel-*.R`, `R/regplot-*.R`, and
`R/zplot-*.R`. Inspect the current family before adding files or helpers.
