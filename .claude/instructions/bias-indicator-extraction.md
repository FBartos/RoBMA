# Bias Indicator Extraction for RoBMA

RoBMA models average over multiple publication bias models. The `bias_indicator` 
column in posterior samples tracks which bias model generated each sample.

## Key Function: `.extract_use_normal()`

Located in `R/evaluate.R`, this function returns a logical vector indicating which 
posterior samples use non-weightfunction bias models (and can use fast normal path).

## How It Works

### For RoBMA (mixture priors)
```r
# bias_indicator column exists, values are 1-indexed bias model numbers
posterior_samples <- coda::as.mcmc(object[["fit"]])
bias_indicator <- as.integer(posterior_samples[, "bias_indicator"])

# Identify which priors are weightfunctions
priors_bias <- object[["priors"]][["outcome"]][["bias"]]
wf_indices <- which(sapply(priors_bias, BayesTools::is.prior.weightfunction))

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

Call once at top of `.pdf.brma()`, `.cdf.brma()`, `predict.brma()`:

```r
if (is_weightfunction) {
  use_normal <- .extract_use_normal(object)
  
  log_lik <- .outcome_pdf.wnorm(..., use_normal = use_normal)
}
```

## Why This Approach

- **Authoritative source**: Uses JAGS-set `bias_indicator` rather than floating-point omega comparisons
- **Robust**: Works for both brma (single prior) and RoBMA (mixture priors)
- **Type-safe**: Uses `BayesTools::is.prior.weightfunction()` for prior detection

## JAGS Behavior

- For non-weightfunction samples, JAGS sets omega = all 1s
- Weighted normal with omega = 1 mathematically equals normal
- This is why we can always use `.outcome_pdf.wnorm()` when `is_weightfunction = TRUE`

## Files

- `R/evaluate.R` - `.extract_use_normal()` definition
- `tests/testthat/test-02-evaluate.R` - Tests for `.extract_use_normal()`
