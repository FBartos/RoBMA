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
| `ROBMA_TEST_FILES_DIR` | Cache directory location | `../temp/RoBMA_test_files` |
| `ROBMA_TEST_SKIP_REFIT` | Skip fitting if cache exists | TRUE |

### Recommended TDD Workflow

```r
# 1. Run full suite once to verify current code and populate cache (if missing) 
devtools::test()

# 2. Iterate on your feature - uses cached fits unless you are modifing model fitting!
devtools::test(filter = "your-feature")

# 3. Final verification (disable cache if fit / marglik code or its dependencies changed)
clean_cached_fits()
devtools::test()
```

### When to Clear Cache

Clear with `clean_cached_fits()` when 
- You modify model fitting logic.
- Model definitions in `test-01-XXX.R`

## Key Rules

### Model Fitting

- **Only `test-01-*`** fits models and computes marginal likelihoods
- Other tests load with `list_fits(name = )`, `load_info(name = )`, `load_marglik(name = )`


### Skip Conditions

| Condition | When to Use |
|-----------|-------------|
| `skip_if_no_fits()` | Test loads pre-fitted models |

### Helper Functions (common-functions.R)

```r
source(testthat::test_path("common-functions.R"))
```

## AI Agent Protocol

1. **Scan `test-01-*` first** - understand available models
2. **Reuse existing models** - don't create duplicates
3. **Never fit models** outside `test-01-*`
4. **Never modify** `GENERATE_REFERENCE_FILES` flag (maintainer only)

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Pre-fitted models not available" | Run `devtools::test(filter = "01-")` |
| Stale cache causing failures | `clean_cached_fits()` then rerun |
| Tests pass locally, fail on CI | Clear cache, run full suite |

````
