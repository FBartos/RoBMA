# Test File Template Pattern

When creating new `test-02-*.R` files that use pre-fitted models, follow this template:

## Standard Header

```r
context("Feature Name")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)
```

## Test Block Structure

```r
# ============================================================================ #
# Test: Descriptive Test Name
# ============================================================================ #

test_that("Description of what is being tested", {

  name        <- "fit_name"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Section comment for logical grouping
  # --------------------------------------------------

  # Test code here
  expect_equal(...)
})
```

## Key Elements

1. **Context**: Describes the feature being tested
2. **Common functions**: Load shared helpers and set reference directory
3. **Skip conditions**: Skip if fits not available or packages not installed
4. **Lazy fits/info**: Use `lazy_fits()` and `lazy_infos()` unless eager loading is required
5. **Section separators**: Use `# === #` for major sections, `# --- #` for subsections
6. **Info parameter**: Always include `info = "description"` in `expect_*` calls

## File Naming

- `test-01-*.R`: Model fitting tests (generate cached fits)
- `test-02-*.R`: Feature tests using cached fits
- `test-03-*.R`: Integration tests, etc.
