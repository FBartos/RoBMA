---
name: posterior-sample-manipulation
description: Guide for working with posterior samples from the models (i.e., extracting parameter samples, transforming samples, computing posterior quantities).
---

# Skill: Posterior Sample Manipulation

## When to Use This Skill

Use this skill when working with:
- Extracting or transforming posterior samples from brma/RoBMA fitted models
- Implementing prediction functions that use JAGS posterior samples
- Debugging sign/direction issues in effect size estimates
- Working with selection models, PET, or PEESE bias adjustments

## Effect Direction Flipping in brma Models

### The Problem

When `effect_direction == "negative"` is specified, the JAGS model internally works in "positive space" to simplify selection model and PET/PEESE bias adjustment computations. This creates confusion about when to flip effect signs during prediction and post-processing.

### What JAGS Returns

**Key insight:** The JAGS model uses `-mu` in the likelihood when `effect_direction == "negative"`. This means:

1. **`mu` samples from JAGS are already in the original (negative) scale** - no flipping needed when extracting `mu`
2. **Study-level effects (`theta`, `gamma`)** are also returned in original scale - no additional flipping needed
3. **Selection weights and critical values (`crit_yi`)** are computed in "positive" (flipped) space because the selection process operates on p-values which depend on the absolute effect direction

### Correct Handling by Component

| Component | What JAGS returns | Transformation needed |
|-----------|-------------------|----------------------|
| `mu` (pooled effect) | Original scale | None |
| `tau` (heterogeneity) | Always positive | None |
| `gamma` (study random effects) | Original scale | None |
| `theta` (true study effects) | Original scale | None |
| PET/PEESE bias term | Computed for positive effects | **Subtract** (not add) for negative effects |
| Selection model sampling | Weights computed in positive space | Flip `mu` before sampling, flip result back |

### PET/PEESE Bias Adjustment

For bias-unadjusted predictions with PET/PEESE:
- The bias term `PET * sei` inflates effects toward significance
- For **positive** effects: `mu_biased = mu + bias_term`
- For **negative** effects: `mu_biased = mu - bias_term` (bias pushes toward more negative values)

```r
# Example from brma.evaluate.R
direction <- ifelse(effect_direction == "negative", -1, 1)
mu_samples <- mu_samples + direction * outer(PET_samples, sei_vec)
```

### Selection Model Weighted Sampling

When sampling from the weighted distribution (`type = "response"`, `bias_adjusted = FALSE`):
- Critical values `crit_yi` are computed in positive (flipped) space
- Must flip `mu` to positive before sampling from weighted normal
- Flip the sampled values back to original scale

```r
# Example from brma.predict.R
if (effect_direction == "negative") {
  # Flip mu to positive space for weighted sampling
  mu_samples_for_wnorm <- -mu_samples
} else {
  mu_samples_for_wnorm <- mu_samples
}

outcome_samples <- .outcome_rng.wnorm(
  mu_samples = mu_samples_for_wnorm,
  tau_within = tau_within_samples,
  sei        = outcome_data[["sei"]],
  omega      = omega_samples,
  crit_yi    = fit_data$crit_yi
)

# Flip samples back to original space
if (effect_direction == "negative") {
  outcome_samples <- -outcome_samples
}
```

### Common Mistakes to Avoid

1. **Double-flipping `mu`**: JAGS already returns `mu` in original scale; don't flip it again
2. **Flipping `tau`**: Heterogeneity is always positive; never flip it
3. **Ignoring direction for PET/PEESE**: The bias adjustment direction matters
4. **Forgetting to flip back**: When sampling in positive space, always flip results back
5. **Flipping study effects (`gamma`)**: These are already in original scale

### Key Files

- `R/brma.evaluate.R` - Contains `.evaluate.brma.mu()` and `.evaluate.brma.study_effects()`
- `R/brma.predict.R` - Contains prediction logic and weighted sampling
- `R/fit.R` - Contains JAGS model syntax generation showing how `-mu` is used

### Testing Effect Direction Handling

Always test both positive and negative effect directions:
- Fit models with `effect_direction = "positive"` and `effect_direction = "negative"`
- Verify that flipped data produces flipped results
- Compare with metafor reference implementations when available
- See `tests/testthat/test-02-predict.R` for examples of negative effect tests
