---
name: unit-testing
description: Guide for running and writing unit tests in the RoBMA package. Use this when working with tests in tests/testthat/.
---

# Unit Testing in RoBMA

This skill helps you run, debug, and write unit tests for the RoBMA package.

## Test Structure Overview

Tests are organized by prefix number indicating their role:

| Prefix | Purpose | Examples |
|--------|---------|----------|
| `test-00-*` | Input validation, CRAN checks | `test-00-input-data-norm.R` |
| `test-01-*` | **Model fitting** (creates cached fits) | `test-01-brma.norm.R`, `test-01-brma.glmm.R` |
| `test-02-*` | Core method tests using cached fits | `test-02-residuals.R`, `test-02-predict.R` |
| `test-03+` | Higher-level integration tests | `test-04-fit.R`, `test-05-methods.R` |

### Critical Rule: Model Fitting Centralization

**Only `test-01-*` files fit models.** All other tests load pre-fitted models from cache.

## Running Tests

Always use testthat LLM reporting for unit tests. Prefer
`devtools::test(..., reporter = "llm")`; if a wrapper cannot pass `reporter`,
set `AGENT=1` so testthat selects `LlmReporter`.

### Basic Commands

```r
# Run all tests (fits models on first run, caches them)
devtools::test(reporter = "llm")

# Run specific test file
devtools::test(filter = "residuals", reporter = "llm")

# Run only model fitting tests (populates cache)
devtools::test(filter = "01-", reporter = "llm")

# Run tests matching multiple patterns
devtools::test(filter = "predict|residuals", reporter = "llm")
```

### Recommended Workflow

```r
# 1. First run: populate the cache (slow, ~2 min)
devtools::test(filter = "01-", reporter = "llm")

# 2. Iterate on your feature (fast, uses cached fits)
devtools::test(filter = "your-feature", reporter = "llm")

# 3. Final verification before commit
devtools::test(reporter = "llm")
```

## Cache System

### How It Works

- Fitted models are saved to `ROBMA_TEST_FILES_DIR`
- Subdirectories: `fits/`, `margliks/`, `info/`, `temp/`
- Cache persists across R sessions

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `ROBMA_TEST_FILES_DIR` | Cache directory location
| `ROBMA_TEST_SKIP_REFIT` | Skip fitting if cache exists 

### When to Clear Cache

Clear cache when:
- You modify model fitting logic in `R/brma.norm.R`, `R/fit.R`
- You change model definitions in `test-01-*.R` files
- Tests fail unexpectedly due to stale cached fits

```r
# Clear all cached fits
source(testthat::test_path("common-functions.R"))
clean_cached_fits()

## Writing New Tests

### For Tests Using Cached Fits (test-02-* and higher)

```r
# At the top of test file
source(testthat::test_path("common-functions.R"))

# Load pre-fitted models
fits <- list()
info <- list()
for (name in list_fits()) {
  fits[[name]] <- load_fit(name)
  info[[name]] <- load_info(name)
}

# Skip if no fits available
test_that("My test uses cached fits", {
  skip_if_no_fits()
  
  fit_brma    <- fits[["bcg_meta-analysis"]]
  fit_metafor <- info[["bcg_meta-analysis"]][["metafor"]]
  
  # Your test assertions...
  expect_equal(...)
})
```

### For Tests Fitting New Models (test-01-* only)

```r
test_that("Fit my new model type", {
  
  # Skip if already cached
  name <- "my-new-model"
  skip_refit_if_cached(name)
  
  # Fit the model
  fit_brma <- brma(yi ~ 1, sei = sei, data = my_data, ...)
  
  # Also fit metafor reference (optional)
  fit_metafor <- metafor::rma(yi = yi, sei = sei, data = my_data)
  
  # Save to cache
  save_fit(
    name    = name,
    fit     = fit_brma,
    info    = list(metafor = fit_metafor, data = my_data)
  )
  
  # Basic sanity checks

  expect_s3_class(fit_brma, "brma")
})
```

## Helper Functions (common-functions.R)

| Function | Purpose |
|----------|---------|
| `load_fit(name)` | Load cached brma fit |
| `load_info(name)` | Load associated metadata (e.g., metafor reference fit) |
| `load_marglik(name)` | Load cached marginal likelihoods |
| `list_fits()` | List all available cached model names |
| `save_fit(name, fit, info, marglik)` | Save model to cache |
| `skip_if_no_fits()` | Skip test if cache is empty |
| `skip_refit_if_cached(name)` | Skip model fitting if already cached |
| `clean_cached_fits()` | Clear entire cache |
| `clean_cached_fits(name)` | Clear specific model from cache |

## Comparing with metafor

Many tests compare brma results to metafor as a reference implementation:

```r
test_that("My method matches metafor", {
  skip_if_not_installed("metafor")
  skip_if_no_fits()
  
  name        <- "bcg_meta-analysis"
  fit_brma    <- fits[[name]]
  fit_metafor <- info[[name]][["metafor"]]
  
  # Compare results (allow MCMC tolerance)
  brma_result   <- my_method(fit_brma)
  metafor_result <- my_method(fit_metafor)
  
  expect_equal(brma_result, metafor_result, tolerance = 0.10,
               info = "brma should match metafor within MCMC tolerance")
})
```

### Tolerance Guidelines

- Use 0.05 whenever possible.
- Some output in complex models might require 0.10.

## Numerical Regressions

Never update results of commited tests. If tests fail there was either
- an error in the new implementation (fix the implementation, do not skip the test) 
- an error in the previous implementation (notify the user to verify the error)

## Visual Regression Tests (vdiffr)

Some tests use `vdiffr::expect_doppelganger()` for plot comparisons.
Never update visual comparisons. Ask the user to manually check the updates/

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Pre-fitted models not found" | Run `devtools::test(filter = "01-", reporter = "llm")` |
| Tests pass locally, fail on CI | Check if CI has cached fits; may need to run `test-01-*` first |
| Stale cache causing failures | Run `clean_cached_fits()` then retest |
| MCMC tolerance too tight | Increase `tolerance` in `expect_equal()` |
| Snapshot test failures | Ask for human review |
| "JAGS module could not be loaded" | Expected warning on Windows without JAGS; tests still run |

## Test File Naming Conventions

When adding new test files:

1. **Input validation**: `test-00-input-{topic}.R`
2. **Model fitting**: `test-01-{model-type}.R`
3. **Core methods**: `test-02-{method}.R` (residuals, predict, summary, etc.)
4. **Integration**: `test-03+` with sequential numbering

## AI Agent Protocol

1. **Always source common-functions.R** at the top of test files
2. **Always use LLM reporting** when running unit tests: `reporter = "llm"` or `AGENT=1`
3. **Check `list_fits()`** to see available cached models before writing tests
4. **Never fit models** outside `test-01-*` files
5. **Reuse existing cached models** - check `test-01-*.R` for available model types
6. **Use appropriate tolerances** for MCMC-based comparisons
7. **Never modify `GENERATE_REFERENCE_FILES`** flag (maintainer only)
