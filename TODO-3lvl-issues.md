# TODO: 3-Level Issues

Snapshot date: 2026-04-24

This note captures the current open questions around 3-level normal models.
It is intended as a restart point for a later session.

## Current Status

- The residual / hat-matrix defect was traced to wrong same-data 3-level normal fitted values.
- A custom metafor oracle was added in `tests/testthat/common-functions.R` because `metafor::rstandard(type = "conditional")` does not provide the right oracle for `rma.mv`.
- The currently open question is not only mathematical. It is also semantic:
  - what exactly `study`, `estimate`, and `conditional` should mean in 3-level models,
  - whether `study_ids` should stay as-is or gain a `cluster` alias,
  - whether cluster-level residual diagnostics should be added as a separate concept.
- The `predict.R` change is the main piece still under review.

## Issue 1: Same-data 3-level normal predictions

### Problem

For same-data 3-level normal models, the old logic built predictions in two steps:

1. add the fitted study effect `gamma * tau_between`
2. apply the usual 2-level scalar shrinkage for the estimate-level BLUP

That decomposition is not correct in a 3-level model, because observations inside the same study are coupled through the full study block of the marginal covariance matrix.

The correct same-data conditional predictor is joint, not sequential.

### Current candidate fix

Current branch in `R/predict.R` introduces `multilevel_blup` and routes same-data 3-level normal predictions through:

- `R/predict.R:281` - `multilevel_blup <- NULL`
- `R/predict.R:283` - `.evaluate.brma.multilevel_blup.norm(...)`
- `R/predict.R:295` - use joint study contribution for `type = "study"`
- `R/predict.R:332` - use joint estimate contribution for `type = "estimate"`

The helper is:

- `R/evaluate.R:291` - `.evaluate.brma.multilevel_blup.norm(...)`

It implements the mixed-model identity:

- `b_hat = G M^{-1} (y - X beta)`
- `M = V + diag(tau_within^2) + block(tau_between %*% t(tau_between))`

Relevant lines:

- `R/evaluate.R:318` - compute `M^{-1} (y - X beta)`
- `R/evaluate.R:320` - estimate-level contribution
- `R/evaluate.R:323` - study-level block contribution

### Follow-up decision

Decide whether to keep this `predict.R` hook as the intended design, or refactor the same logic into a more explicit internal API.

### Minimum validation

Add or keep direct oracle tests for:

- `predict(..., type = "study")` on same-data 3-level normal fits
- `predict(..., type = "estimate")` on same-data 3-level normal fits

Do not rely only on residual tests.

## Issue 2: Low-level helper hazard in `evaluate.R`

The old approximation still exists in the low-level helpers:

- `R/evaluate.R:360` - `.evaluate.brma.study_effects(...)`
- `R/evaluate.R:428` - `.evaluate.brma.true_effects.norm(...)`

These helpers are still mathematically fine for:

- non-multilevel models,
- new-data prediction,
- code paths where study and estimate effects are intentionally handled separately.

But they are hazardous for same-data 3-level normal conditional predictions.

### Why this matters

Right now the bug is avoided because `predict.brma()` intercepts same-data 3-level normal predictions first.
If a future caller uses `.evaluate.brma.study_effects()` plus `.evaluate.brma.true_effects.norm()` directly, the old defect can reappear.

### Helper hardening

Choose one of these:

1. add an explicit guard that rejects same-data 3-level normal use of `.evaluate.brma.true_effects.norm()`,
2. refactor the helper layer so the multilevel-aware path is impossible to bypass accidentally,
3. document the restriction very clearly if the API remains split.

## Issue 3: Residual semantics are still ambiguous

Current residual API:

- `R/residuals.R:216` - `estimate` is documented as alias `conditional`
- `R/residuals.R:263` - `conditional` is normalized to `estimate`
- `R/residuals.R:359` - same for `rstudent`
- `R/residuals.R:417` - `conditional` is normalized to `estimate`

### Naming ambiguity

In 3-level models, there are at least two different observation-level conditioning targets:

- `study`: condition on fixed effects plus the study-level random effect
- `estimate`: condition on the full latent estimate-specific effect

So `conditional` is ambiguous in 3-level models.

### Important distinction from metafor

Metafor uses `cluster` for cluster-level multivariate residual diagnostics for `rma.mv`.
That is not the same thing as the current RoBMA `study` residual type.

Relevant current anchors:

- `tests/testthat/test-02-residuals.R:262` - comment that metafor ignores `type = "conditional"` for `rma.mv`
- `tests/testthat/test-02-residuals.R:290` - `# TODO: add clustered residuals`
- `tests/testthat/test-02-funnel.R:141` - `# add cluster-level residuals once implemented`

### API decision

Decide whether to:

1. deprecate `type = "conditional"` for 3-level residuals and require explicit `study` or `estimate`,
2. keep `conditional` as alias to `estimate`, accepting the ambiguity,
3. introduce a new cluster-level diagnostic and reserve `clustered` for that concept.

## Issue 4: Missing cluster-level diagnostics

There is still no actual cluster-level residual / influence diagnostic aligned with metafor's `cluster` concept.

Current state:

- observation-level residuals exist for `marginal`, `study`, and `estimate`
- cluster-level multivariate standardized residuals do not exist yet

The tests already acknowledge this gap:

- `tests/testthat/test-02-residuals.R:290`
- `tests/testthat/test-02-funnel.R:141`

### Missing feature

If metafor alignment is desired, implement cluster-level diagnostics as a distinct feature instead of renaming current observation-level residuals.

## Issue 5: `study_ids` naming vs `cluster`

Public input docs currently say:

- `R/input-data.R:16` - `study_ids` is "used for clustering"

Internally, the codebase consistently uses `study_ids`:

- `R/fit.R`
- `R/input-data.R`
- `R/predict.R`
- `R/evaluate.R`

However, many tests already use a local variable named `cluster` and pass it as `study_ids = cluster`.

### Naming decision

Decide whether to:

1. keep `study_ids` only,
2. add `cluster` as a public alias while keeping `study_ids` internally,
3. fully rename the public API to `cluster` / `cluster_ids`.

If changing names, check docs, examples, input validation, and tests together.

## Issue 6: `cdf` and `pdf` need a semantics review

### `cdf`

`R/cdf.R` currently pulls means from `predict.brma()`:

- `R/cdf.R:249` - `mu_samples <- predict.brma(...)`
- `R/cdf.R:264` - `tau_within_samples <- predict.brma(..., type = "terms.scale")`

That means `cdf` inherits whatever same-data 3-level normal logic `predict.brma()` uses.

So `cdf` does not independently recreate the old bug.

### `pdf`

`R/pdf.R` also routes through `predict.brma()` for the mean:

- `R/pdf.R:456` - `mu_type <- if (is_multilevel) "study" else "terms"`
- `R/pdf.R:457` - `mu_samples <- predict.brma(...)`

So `pdf` also does not independently recreate the old same-data predictor bug.

However, `pdf` feeds LOO / PSIS in `R/loo.R`, so there is a separate modeling question:

- should pointwise multilevel LOO be observation-level conditional on study effects,
- marginal over study effects,
- or cluster-level deletion instead of case deletion?

This is not the same bug, but it still needs an explicit decision.

### LOO decision

Write down the intended target for 3-level `pdf` / `loo` and then verify that:

- `R/pdf.R`
- `R/cdf.R`
- `R/residuals.R`
- `R/loo.R`

all use the same conditioning convention.

## Issue 7: Oracle and tests

Current oracle helper:

- `tests/testthat/common-functions.R:158` - custom multilevel conditional residual oracle for metafor `rma.mv`

Important detail:

- it reconstructs the BLUP-based conditional residual mean using `stats::fitted(model)` plus `metafor::ranef(model, expand = TRUE)`
- it uses the GLS residual variance matrix, not a scalar observation-wise shrinkage approximation

Current residual test anchors:

- `tests/testthat/test-02-residuals.R:264` - use custom oracle for multilevel conditional residuals
- `tests/testthat/test-02-residuals.R:266` - compare `rstandard(..., type = "conditional")`
- `tests/testthat/test-02-residuals.R:267` - compare `predict(..., type = "estimate")` against oracle-derived BLUPs

### Test follow-up

Strengthen tests with direct, narrow oracle checks for:

- same-data `predict(..., type = "study")`
- same-data `predict(..., type = "estimate")`
- any future cluster-level diagnostic

## Suggested Next Session Order

1. Settle terminology first.
2. Re-derive the intended 3-level conditioning targets in writing.
3. Decide whether the `predict.R` `multilevel_blup` path is the intended permanent fix.
4. Add direct oracle tests for `study` and `estimate` predictions.
5. Harden `evaluate.R` so the old approximation cannot be reused accidentally.
6. Review `cdf` / `pdf` / `loo` semantics after the conditioning terminology is fixed.
7. Only then decide whether to add `cluster` as an alias for `study_ids` and whether to add true cluster-level diagnostics.

## Useful Commands

These were the relevant targeted checks during the investigation:

```r
devtools::test(filter = "residuals")
devtools::test(filter = "predict|hatvalues|influence|dfbetas|evaluate")
devtools::test(filter = "funnel|qqnorm|radial")
```

## Out of Scope / Noise

- Unrelated zcurve snapshot diffs should be ignored when working on this topic.
- A `regplot` test failure seen during the investigation came from test harness usage, not from the 3-level predictor math.
