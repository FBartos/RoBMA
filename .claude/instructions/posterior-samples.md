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
3. **Selection weights and critical values (`crit_yi`)** are computed in "positive" (flipped) space

### Correct Handling by Component

| Component | What JAGS returns | Transformation needed |
|-----------|-------------------|----------------------|
| `mu` (pooled effect) | Original scale | None |
| `tau` (heterogeneity) | Always positive | None |
| `gamma` (study random effects) | Original scale | None |
| `theta` (true study effects) | Original scale | None |
| PET/PEESE bias term | Computed for positive effects | **Subtract** for negative effects |
| Selection model sampling | Weights in positive space | Flip `mu` before sampling, flip result back |
| RNG computation | Computed in original space | **Flip RNG** (- RNG) for negative effects |
| CDF computation | Computed in original space | **Flip CDF** (1 - CDF) for negative effects |

### PET/PEESE Bias Adjustment

For bias-unadjusted predictions with PET/PEESE:

- The bias term `PET * sei` inflates effects toward significance
- For **positive** effects: `mu_biased = mu + bias_term`
- For **negative** effects: `mu_biased = mu - bias_term`

```r
# Example from brma.evaluate.R
direction <- ifelse(effect_direction == "negative", -1, 1)
mu_samples <- mu_samples + direction * outer(PET_samples, sei_vec)
```

### Selection Model Weighted Sampling

When sampling from the weighted distribution (`type = "response"`, `bias_adjusted = FALSE`):

```r
# Example from brma.predict.R
if (effect_direction == "negative") {
  # Flip mu to positive space for weighted sampling
  mu_samples_for_wnorm <- -mu_samples
} else {
  mu_samples_for_wnorm <- mu_samples
}

outcome_samples <- .outcome_rng.wnorm(...)

# Flip samples back to original space
if (effect_direction == "negative") {
  outcome_samples <- -outcome_samples
}
```

### CDF Computation

When computing cumulative distribution functions (e.g., for LOO-PIT residuals):

```r
# Example from brma.residuals.R
cdf_matrix <- .cdf.brma(object)
effect_direction <- .effect_direction(object)
if (effect_direction == "negative") {
  cdf_matrix <- 1 - cdf_matrix
}
```

## Core Computation Functions

The brma package uses four families of internal functions for posterior sample manipulation:

| Function Family | Purpose | Key Functions | Effect Direction Handling |
|----------------|---------|---------------|---------------------------|
| `.evaluate.*` | Extract posterior samples from JAGS | `.evaluate.brma.mu()`, `.evaluate.brma.study_effects()` | Returns samples in original scale |
| `.rng.*` | Sample from posterior predictive | `.outcome_rng.norm()`, `.outcome_rng.wnorm()` | Flipping handled by caller |
| `.cdf.*` | Compute CDFs | `.cdf.brma()`, `.outcome_cdf.norm()` | CDF flipping handled by caller |
| `.pdf.*` | Compute log-likelihoods for LOO-PSIS | `.outcome_pdf.norm()`, `.outcome_pdf.wnorm()` | No flipping needed |

**Function Organization**: Most families have a main dispatcher (e.g., `.cdf.brma()`) that handles model type detection. Matrix outputs are consistently S x K (samples x observations).

### Common Mistakes to Avoid

1. **Double-flipping `mu`**: JAGS already returns `mu` in original scale; don't flip it again
2. **Flipping `tau`**: Heterogeneity is always positive; never flip it
3. **Ignoring direction for PET/PEESE**: The bias adjustment direction matters
4. **Forgetting to flip back**: When sampling in positive space, always flip results back
5. **Flipping study effects (`gamma`)**: These are already in original scale
6. **Forgetting rng/cdf flipping**: For LOO-PIT residuals, always flip CDF values for negative effects

### Key Files

- `R/brma.evaluate.R` - Posterior sample extraction (`.evaluate.*` functions)
- `R/brma.predict.R` - Prediction logic and RNG calling
- `R/brma.rng.R` - Posterior predictive sampling (`.outcome_rng.*` functions)
- `R/brma.cdf.R` - CDF computation (`.cdf.brma()`, `.outcome_cdf.*` functions)
- `R/brma.pdf.R` - Log-likelihood computation (`.outcome_pdf.*` functions)
- `R/fit.R` - JAGS model syntax generation

### Testing Effect Direction Handling

Always test both positive and negative effect directions:

- Fit models with `effect_direction = "positive"` and `effect_direction = "negative"`
- Verify that flipped data produces flipped results
- Compare with metafor reference implementations when available
- See `tests/testthat/test-02-predict.R` for examples
