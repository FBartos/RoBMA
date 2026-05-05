---
description: 'R unit tests: testing standards and Copilot guidance for interaction with unit tests.'
applyTo: "**/tests/testthat/*.R"
---

# BayesTools Test Guidelines

## Overview

- Model fitting is centralized in `test-01-*`; other tests load cached models
- Tests use `common-functions.R` for shared helpers

## Test Caching (TDD Workflow)

Model fitting is slow. The caching system lets you run the full suite once and reuse fits.

### Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `ROBMA_TEST_FILES_DIR` | Cache directory location | `tests/testthat/test_files` locally, tempdir on CRAN |
| `ROBMA_TEST_SKIP_REFIT` | Skip fitting if a valid cache exists | TRUE |
| `ROBMA_TEST_FORCE_REFIT` | Force refitting even if a valid cache exists | FALSE |
| `ROBMA_TEST_FULL_DIAGNOSTICS` | Run extended redundant residual/influence diagnostics | FALSE |

Cached fits are validated against `fit_catalog()` in `common-functions.R`.
The cache metadata records source file hashes and expected LOO/WAIC/marglik/metafor availability.
Group calls such as `skip_refit_if_cached("brma.norm")` skip only when all cataloged fits for that group are valid.
`list_fits()` uses metadata-only validation by default; use `validate_cached_fit(name, deep = TRUE)` only for explicit cache-integrity checks.
Extended diagnostic tests use `skip_if_not_full_diagnostics()` when they duplicate core coverage or lack a stable oracle; set `ROBMA_TEST_FULL_DIAGNOSTICS=TRUE` to run them.

### Recommended TDD Workflow

Always use testthat LLM reporting for unit tests. Prefer
`devtools::test(..., reporter = "llm")`; if a wrapper cannot pass `reporter`,
set `AGENT=1` so testthat selects `LlmReporter`.

```r
# 1. Run full suite once to verify current code and populate cache (if missing) 
devtools::test(reporter = "llm")

# 2. Iterate on your feature - uses cached fits unless you are modifing model fitting!
devtools::test(filter = "your-feature", reporter = "llm")

# 3. Final verification (disable cache if fit / marglik code or its dependencies changed)
clean_cached_fits()
devtools::test(reporter = "llm")
```

### When to Clear Cache

Clear with `clean_cached_fits()` when 
- You modify model fitting logic.
- Model definitions in `test-01-XXX.R`

## Key Rules

### Model Fitting

- **Only `test-01-*`** fits models and computes LOO/marglik when expected by `fit_catalog()`
- Other tests discover with `list_fits()` and usually use `lazy_fits()` / `lazy_infos()`
- Avoid top-level eager loading of all cached fits unless the test genuinely needs every object
- `RoBMA`, `BMA.norm`, and `BMA.glmm` product-space fits must not cache marginal likelihoods
- Add every new cached fit to `fit_catalog()` in the same change


### Skip Conditions

| Condition | When to Use |
|-----------|-------------|
| `skip_if_no_fits()` | Test loads pre-fitted models |
| `skip_if_missing_fits(names)` | Test requires specific cached models |
| `skip_if_not_full_diagnostics(reason)` | Extended residual/influence diagnostic duplicates core coverage |

### Helper Functions (common-functions.R)

```r
source(testthat::test_path("common-functions.R"))
```

## AI Agent Protocol

1. **Scan `test-01-*` first** - understand available models
2. **Always use LLM reporting** when running unit tests: `reporter = "llm"` or `AGENT=1`
3. **Reuse existing models** - don't create duplicates
4. **Never fit models** outside `test-01-*`
5. **Never modify** `GENERATE_REFERENCE_FILES` flag (maintainer only)

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Pre-fitted models not available" | Run `devtools::test(filter = "01-", reporter = "llm")` |
| Stale cache causing failures | `clean_cached_fits()` then rerun |
| Tests pass locally, fail on CI | Clear cache, run full suite |

````
