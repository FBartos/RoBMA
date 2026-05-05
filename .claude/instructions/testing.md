# Unit Testing in RoBMA

Reference this file when working with tests in `tests/testthat/`.

## Test Structure Overview

Tests are organized by prefix number indicating their role:

| Prefix | Purpose | Examples |
|--------|---------|----------|
| `test-00-*` | Input validation, CRAN checks | `test-00-input-data-norm.R` |
| `test-01-*` | **Model fitting** (creates cached fits) | `test-01-brma.norm.R` |
| `test-02-*` | Core method tests using cached fits | `test-02-residuals.R`, `test-02-predict.R` |
| `test-03+` | Higher-level integration tests | `test-03-loo.R` |

### Critical Rule: Model Fitting Centralization

**Only `test-01-*` files fit models.** All other tests load pre-fitted models from cache.

## Running Tests

Always use testthat LLM reporting for unit tests. Prefer
`devtools::test(..., reporter = "llm")`; if a wrapper cannot pass `reporter`,
set `AGENT=1` so testthat selects `LlmReporter`.

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
# 1. First run: populate the cache (slow)
devtools::test(filter = "01-", reporter = "llm")

# 2. Iterate on your feature (fast, uses cached fits)
devtools::test(filter = "your-feature", reporter = "llm")

# 3. Final verification before commit
devtools::test(reporter = "llm")
```

## Cache System

### How It Works

- Fitted models are saved to `ROBMA_TEST_FILES_DIR`
- Default local cache: `tests/testthat/test_files`
- CRAN cache: `tempdir()/RoBMA_test_files`
- Subdirectories: `fits/`, `info/`, `metadata/`, `temp/`
- Cache persists across R sessions by default in local development

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `ROBMA_TEST_FILES_DIR` | Cache directory location |
| `ROBMA_TEST_SKIP_REFIT` | Skip fitting if a valid cache exists; defaults to `TRUE` |
| `ROBMA_TEST_FORCE_REFIT` | Force refitting even if a valid cache exists |
| `ROBMA_TEST_FULL_DIAGNOSTICS` | Run extended redundant residual/influence diagnostics skipped by default |

### Cache Validation

`save_fit()` writes the fit, info object, and metadata. Metadata records the normalized `test-01-*` file hash, the cache-affecting fitting-source hash, and whether the object contains LOO, WAIC, marginal likelihood, and metafor info. The fitting-source hash is intentionally scoped to fitting, prior/data input, cached LOO/marglik, and native likelihood/JAGS code; unrelated summaries, plotting, prediction, documentation, and post-fit methods should not force refits.

`skip_refit_if_cached("brma.norm")` and similar group calls skip a whole `test-01-*` file only when all cataloged fits for that group are valid.

`list_fits(validate = TRUE)` performs metadata-only validation so test setup stays fast. `load_fit()` validates the requested fit and memoizes it within the test file. Use `validate_cached_fit(name, deep = TRUE)` only in dedicated cache-integrity tests because it reads the full fit object.

A cached fit is stale if:

- its `test-01-*` source file changed
- cache-affecting fitting source changed
- metadata is missing
- required LOO/marglik/metafor info is missing
- forbidden fields are present, e.g. marginal likelihood for product-space `RoBMA`, `BMA.norm`, or `BMA.glmm`

Use `list_fits(validate = FALSE)` only when inspecting raw cache files.

Some residual and influence tests are intentionally gated with
`skip_if_not_full_diagnostics()`. They duplicate core metafor/model-family
coverage or have no stable oracle, but remain available by setting
`ROBMA_TEST_FULL_DIAGNOSTICS=TRUE`.

### Fit Catalog

`fit_catalog()` in `common-functions.R` is the source of truth for cached test fits. It records:

- expected class and likelihood family
- source `test-01-*` file
- availability of metafor reference, LOO, WAIC, and marginal likelihood
- test tier and feature tags

Use catalog filters instead of ad hoc hard-coded model lists when possible:

```r
list_fits(has_marglik = TRUE)
list_fits(has_metafor = TRUE, feature = "mods")
catalog_fits(class = "RoBMA")
```

When adding a new cached fit in `test-01-*`, add it to `fit_catalog()` in the same change.

### When to Clear Cache

Clear cache when:

- You modify model fitting logic in `R/brma.norm.R`, `R/fit.R`
- You change model definitions in `test-01-*.R` files
- Tests fail unexpectedly due to stale cached fits

```r
# Clear all cached fits
source(testthat::test_path("common-functions.R"))
clean_cached_fits()
```

## Writing New Tests

### For Tests Using Cached Fits (test-02-* and higher)

```r
# At the top of test file
source(testthat::test_path("common-functions.R"))

# Create lazy fit stores; individual RDS files load on first `[[name]]`.
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

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

  expect_s3_class(fit_brma, "brma")
})
```

## Helper Functions (common-functions.R)

| Function | Purpose |
|----------|---------|
| `load_fit(name)` | Load cached brma fit |
| `load_info(name)` | Load associated metadata (e.g., metafor reference fit) |
| `lazy_fits(names)` | Create a lazy cached-fit store |
| `lazy_infos(names)` | Create a lazy cached-info store |
| `load_fits(names)` | Eagerly load a named fit subset |
| `load_infos(names)` | Eagerly load a named info subset |
| `load_fit_metadata(name)` | Load cache metadata |
| `list_fits()` | List valid cached model names |
| `list_fits(validate = FALSE)` | List raw cached model names |
| `catalog_fits()` | List cataloged model names matching feature filters |
| `fit_catalog()` | Full cached-fit catalog |
| `validate_cached_fit(name)` | Return cache validation problems |
| `save_fit(name, fit, info)` | Save model to cache |
| `skip_if_no_fits()` | Skip test if cache is empty |
| `skip_if_missing_fits(names)` | Skip test if required cached fits are missing/stale |
| `skip_refit_if_cached(name)` | Skip model fitting if already cached |
| `clean_cached_fits()` | Clear entire cache |
| `clean_cached_fits(name)` | Clear specific model from cache |

## Comparing with metafor

Many tests compare brma results to metafor as a reference implementation:

```r
test_that("My method matches metafor", {
  skip_if_no_fits()

  name        <- "bcg_meta-analysis"
  fit_brma    <- fits[[name]]
  fit_metafor <- info[[name]][["metafor"]]

  # Compare results (allow MCMC tolerance)
  brma_result    <- my_method(fit_brma)
  metafor_result <- my_method(fit_metafor)

  expect_equal(brma_result, metafor_result, tolerance = 0.10,
               info = "brma should match metafor within MCMC tolerance")
})
```

Put required packages to the top of the test file (even if a subset of tests could run otherwise).
For example, a heading of the file might start
```r
context("DFBETAS")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)
```

### Tolerance Guidelines

- Use 0.05 whenever possible.
- Some output in complex models might require 0.10.

## Numerical Regressions

Never update results of committed tests. If tests fail there was either:

- an error in the new implementation (fix the implementation, do not skip the test)
- an error in the previous implementation (notify the user to verify the error)

## Visual Regression Tests (vdiffr)

Some tests use `vdiffr::expect_doppelganger()` for plot comparisons.
Never update visual comparisons automatically. Ask the user to manually verify updates.

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Pre-fitted models not found" | Run `devtools::test(filter = "01-", reporter = "llm")` |
| Tests pass locally, fail on CI | Check if CI has cached fits |
| Stale cache causing failures | Run `clean_cached_fits()` then retest |
| MCMC tolerance too tight | Increase `tolerance` in `expect_equal()` |
| Snapshot test failures | Ask for human review |

## Agent Protocol

1. **Always source common-functions.R** at the top of test files
2. **Always use LLM reporting** when running unit tests: `reporter = "llm"` or `AGENT=1`
3. **Check `fit_catalog()` / `list_fits()`** to see available cached models before writing tests
4. **Never fit models** outside `test-01-*` files
5. **Reuse existing cached models** - check `fit_catalog()` and `test-01-*.R` for available model types
6. **Use appropriate tolerances** for MCMC-based comparisons
7. **Never modify `GENERATE_REFERENCE_FILES`** flag (maintainer only)
