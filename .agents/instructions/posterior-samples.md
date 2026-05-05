# Posterior Sample Manipulation

Reference this file when working with posterior samples from brma/RoBMA models.

## When to Use This Guide

Use this guide when working with:

- Extracting or transforming posterior samples from brma/RoBMA fitted models
- Implementing prediction functions that use JAGS posterior samples
- Debugging sign/direction issues in effect size estimates
- Working with selection models, PET, or PEESE bias adjustments
- Understanding the core computation functions (`.evaluate`, `.rng`, `.cdf`, `.pdf`)
- Implementing new model diagnostics or prediction methods

## Effect Direction Flipping in brma Models

### The Problem

When `effect_direction == "negative"` is specified, the JAGS model internally works in "positive space" to simplify selection model and PET/PEESE bias adjustment computations. This creates confusion about when to flip effect signs.

### What JAGS Returns

**Key insight:** The JAGS model uses `-mu` in the likelihood when `effect_direction == "negative"`. This means:

1. **`mu` samples from JAGS are already in the original (negative) scale** - no flipping needed
2. **Study-level effects (`theta`, `gamma`)** are also returned in original scale
3. **Selection binning and selected-normal kernels** carry the effect-direction sign through `.selection_context()`

### Correct Handling by Component

| Component | What JAGS returns | Transformation needed |
|-----------|-------------------|----------------------|
| `mu` (pooled effect) | Original scale | None |
| `tau` (heterogeneity) | Always positive | None |
| `gamma` (study random effects) | Original scale | None |
| `theta` (true study effects) | Original scale | None |
| PET/PEESE bias term | Computed for positive effects | **Subtract** for negative effects |
| Selection model sampling | Selected-normal context carries `sign` | Do not manually flip in callers |
| RNG computation | Normal RNG uses original scale; selection RNG uses `.selection_context()` | Caller should not add extra flips |
| CDF computation | Normal CDF flips internally for negative direction; selection CDF uses `.selection_context()` | Do not double-flip |

### PET/PEESE Bias Adjustment

For bias-unadjusted predictions with PET/PEESE:

- The bias term `PET * sei` inflates effects toward significance
- For **positive** effects: `mu_biased = mu + bias_term`
- For **negative** effects: `mu_biased = mu - bias_term`

```r
# Example from R/evaluate.R
direction <- ifelse(effect_direction == "negative", -1, 1)
mu_samples <- mu_samples + direction * outer(PET_samples, sei_vec)
```

### Selection Model Sampling

When sampling from the selected-normal distribution (`type = "response"`, `bias_adjusted = FALSE`), build a selection context and pass it through. The context contains the effect-direction sign and the `use_normal` rows for product-space models.

```r
# Example from R/predict.R
selection_context <- .selection_context(
  object            = object,
  posterior_samples = posterior_samples,
  newdata           = new_data
)

outcome_samples <- .outcome_rng.selnorm(
  mu_samples        = mu_samples,
  tau_within        = tau_within_samples,
  sei               = outcome_data[["sei"]],
  selection_context = selection_context
)
```

### CDF Computation

When computing cumulative distribution functions, let `R/cdf.R` own direction handling. For ordinary normal models it flips `yi`, `mu`, and `lower.tail`; for selection models it passes the selected-normal context through.

```r
# Non-selection normal branch inside R/cdf.R
if (effect_direction == "negative") {
  mu_samples_cdf <- -mu_samples
  yi_cdf         <- -yi
  lower_tail     <- FALSE
} else {
  mu_samples_cdf <- mu_samples
  yi_cdf         <- yi
  lower_tail     <- TRUE
}
```

## Core Computation Functions

The package uses four families of internal functions for posterior sample manipulation:

| Function Family | Purpose | Key Functions | Effect Direction Handling |
|----------------|---------|---------------|---------------------------|
| `.evaluate.*` | Extract posterior samples from JAGS | `.evaluate.brma.mu()`, `.evaluate.brma.tau()`, `.evaluate.brma.true_effects.norm()` | Returns samples in original scale |
| `.rng.*` | Sample from posterior predictive | `.outcome_rng.norm()`, `.outcome_rng.selnorm()` | Selection direction handled by `.selection_context()` |
| `.cdf.*` | Compute CDFs | `.cdf.brma()`, `.outcome_cdf.norm()`, `.outcome_cdf.selnorm()` | Normal direction handled in `R/cdf.R`; selection direction handled by `.selection_context()` |
| `.pdf.*` | Compute log-likelihoods for LOO-PSIS | `.outcome_pdf.norm()`, `.outcome_pdf.selnorm()` | Selection direction handled by `.selection_context()` |

**Function Organization**: Most families have a main dispatcher (e.g., `.cdf.brma()`) that handles model type detection. Matrix outputs are consistently S x K (samples x observations).

### Common Mistakes to Avoid

1. **Double-flipping `mu`**: JAGS already returns `mu` in original scale; don't flip it again
2. **Flipping `tau`**: Heterogeneity is always positive; never flip it
3. **Ignoring direction for PET/PEESE**: The bias adjustment direction matters
4. **Manually flipping selection samples**: selected-normal helpers already receive sign information through `.selection_context()`
5. **Flipping study effects (`gamma`)**: These are already in original scale
6. **Double-flipping CDF/RNG values**: check whether the branch is ordinary normal or selected-normal before adding direction logic

### Key Files

- `R/evaluate.R` - Posterior sample extraction (`.evaluate.*` functions)
- `R/predict.R` - Prediction logic and RNG calling
- `R/rng.R` - Posterior predictive sampling (`.outcome_rng.*` functions)
- `R/cdf.R` - CDF computation (`.cdf.brma()`, `.outcome_cdf.*` functions)
- `R/pdf.R` - Log-likelihood computation (`.outcome_pdf.*` functions)
- `R/selection-mapping.R` - Selection-kernel context and native selected-normal calls
- `R/fit.R` - JAGS model syntax generation

### Testing Effect Direction Handling

Always test both positive and negative effect directions:

- Fit models with `effect_direction = "positive"` and `effect_direction = "negative"`
- Verify that flipped data produces flipped results
- Compare with metafor reference implementations when available
- See `tests/testthat/test-02-predict.R` for examples
