# use_normal Selection-Kernel Pattern

Use this guide when working on selected-normal publication-bias code.

## Purpose

`use_normal` is a logical vector of length S, where S is the number of posterior samples.
It marks posterior rows that come from non-selection bias branches in a product-space fit.

- `use_normal[i] = TRUE`: row can use the ordinary normal kernel.
- `use_normal[i] = FALSE`: row must use the selected-normal kernel.

## Current Pattern

Do not pass `use_normal` around as a loose argument except inside selection helpers.
Build a `.selection_context()` and pass that context to selected-normal functions.

```r
selection_context <- .selection_context(
  object            = object,
  posterior_samples = posterior_samples,
  newdata           = new_data
)

log_lik <- .outcome_pdf.selnorm(
  yi                = yi,
  mu_samples        = mu_samples,
  tau_within        = tau_within,
  sei               = sei,
  selection_context = selection_context
)
```

`.selection_context()` stores:

- `omega`, `alpha`, `phack_kind`, and `kernel_mode` rows
- `bias_indicator`
- `use_normal`
- selection bounds, bins, sign, and telescope settings

For `use_normal` rows, `.selection_context()` sets `kernel_mode` to `SELKERNEL_NORMAL`.

## Where It Is Used

| Function | Normal handling | Selected handling |
|----------|-----------------|-------------------|
| `.outcome_pdf.selnorm()` | native selected-normal kernel with `SELKERNEL_NORMAL` rows | native selected-normal kernel |
| `.outcome_cdf.selnorm()` | native selected-normal CDF with `SELKERNEL_NORMAL` rows | native selected-normal CDF |
| `.outcome_rng.selnorm()` | explicit `.outcome_rng.norm()` fast path for normal rows | `.selection_step_rng_matrix()` |

## Computing use_normal

Use `.extract_use_normal()` in `R/evaluate.R`.
It uses `bias_indicator` for RoBMA mixture priors and falls back to the single-prior model type for `brma`.

```r
use_normal <- .extract_use_normal(
  object            = object,
  posterior_samples = posterior_samples
)
```

Do not infer this from `omega == 1`; that is brittle and misses branch intent.

## Deferred p-hacking kernels

Selected-normal p-hacking branches exist as internal scaffolding only.
They are hidden from users and deferred for later implementation.
Do not expose them, activate them, or treat them as supported selected-normal behavior unless the maintainer explicitly asks for that work.

## Files

- `R/evaluate.R`: `.extract_bias_indicator()`, `.extract_use_normal()`
- `R/selection-mapping.R`: `.selection_context()` and native selected-normal wrappers
- `R/pdf.R`, `R/cdf.R`, `R/rng.R`: selected-normal outcome helpers
- `tests/testthat/test-02-evaluate.R`: `use_normal` tests
