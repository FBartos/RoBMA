# Bias Indicator Extraction for RoBMA

RoBMA models average over multiple publication bias models. The `bias_indicator` 
column in posterior samples tracks which bias model generated each sample.

## Key Functions

Located in `R/evaluate.R`:

- `.extract_bias_indicator()` returns posterior bias-branch indicators, defaulting to `1L` for single-prior models.
- `.extract_use_normal()` returns a logical vector indicating which posterior samples use non-selection bias branches and can run the normal selected-kernel mode.

## How It Works

### For RoBMA (mixture priors)
```r
# bias_indicator column exists, values are 1-indexed bias model numbers
posterior_samples <- coda::as.mcmc(object[["fit"]])
bias_indicator <- as.integer(posterior_samples[, "bias_indicator"])

# Identify which priors use selection kernels
priors_bias <- object[["priors"]][["outcome"]][["bias"]]
wf_indices <- which(sapply(priors_bias, .prior_is_selection_kernel))

# use_normal = TRUE for samples NOT from weightfunction priors
use_normal <- !(bias_indicator %in% wf_indices)
```

### For brma (single prior)
```r
# No bias_indicator column - all samples from same model
is_weightfunction <- .is_priors_weightfunction(priors)
use_normal <- rep(!is_weightfunction, S)
```

## Usage Pattern

Build a `.selection_context()` once per posterior-sample set and pass it to selected-normal helpers:

```r
if (is_weightfunction) {
  selection_context <- .selection_context(
    object            = object,
    posterior_samples = posterior_samples
  )

  log_lik <- .outcome_pdf.selnorm(..., selection_context = selection_context)
}
```

## Why This Approach

- **Authoritative source**: Uses JAGS-set `bias_indicator` rather than floating-point omega comparisons
- **Robust**: Works for both brma (single prior) and RoBMA (mixture priors)
- **Type-safe**: Uses `.prior_is_selection_kernel()` so step weightfunctions and selected-normal kernels stay grouped

## JAGS Behavior

- For non-selection samples in a product-space fit, JAGS stores normal-kernel rows.
- `.selection_context()` sets `kernel_mode = SELKERNEL_NORMAL` for those rows.
- Do not detect normal rows by testing `omega` values.

## Deferred p-hacking kernels

Selected-normal p-hacking branches are scaffolded internally but intentionally deferred.
They are hidden from users and should not be exposed, activated, or documented as supported until the maintainer explicitly asks for that implementation.

## Files

- `R/evaluate.R` - `.extract_bias_indicator()` and `.extract_use_normal()`
- `R/selection-mapping.R` - `.selection_context()` and selected-normal kernel routing
- `tests/testthat/test-02-evaluate.R` - Tests for `.extract_use_normal()`
