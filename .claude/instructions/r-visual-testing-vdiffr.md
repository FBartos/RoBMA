# Visual Testing with vdiffr for R Packages

**Domain**: R package testing, visualization
**Confidence**: High
**Source**: RoBMA package test patterns (test-02-funnel.R, test-02-regplot.R)

## Pattern

Use `vdiffr::expect_doppelganger()` for visual regression testing of plotting functions. Compare against reference implementations (like metafor) when available.

### Basic Structure

```r
context("Plot function name")

# Load test helpers and cached fits
source(testthat::test_path("common-functions.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()
```

### Side-by-Side Comparison with Reference

When comparing to a reference implementation (e.g., metafor):

```r
test_that("Plot matches metafor structure", {

  name        <- "dataset_model-type"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # Side-by-side visual comparison
  vdiffr::expect_doppelganger("plotname_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::plotfunc(fit_metafor, main = "metafor")
    plotfunc(fit_brma, plot_type = "base", main = "brma")
  })

  # ggplot version standalone
  vdiffr::expect_doppelganger(
    "plotname_brma_ggplot",
    plotfunc(fit_brma, plot_type = "ggplot")
  )
})
```

### Testing with Fixed Parameters

When comparing plots, use identical parameters to ensure fair comparison:

```r
vdiffr::expect_doppelganger("funnel_comparison", function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  metafor::funnel(fit_metafor, main = "metafor", xlim = c(-3, 3), ylim = c(0, 0.8))
  funnel(fit_brma, plot_type = "base", xlim = c(-3, 3), ylim = c(0, 0.8), main = "brma")
})
```

### Testing Customization Options

Test each customization parameter:

```r
test_that("Plot customization works", {

  fit_brma <- fits[["dataset_model"]]

  # Point aesthetics
  vdiffr::expect_doppelganger("plot_custom_points_base", function() {
    plotfunc(fit_brma, plot_type = "base", pch = 19, col = "blue", bg = "lightblue")
  })

  vdiffr::expect_doppelganger(
    "plot_custom_points_ggplot",
    plotfunc(fit_brma, plot_type = "ggplot", pch = 19, col = "blue", bg = "lightblue")
  )

  # Band/region styling
  vdiffr::expect_doppelganger("plot_custom_bands_base", function() {
    plotfunc(fit_brma, plot_type = "base", shade = "lightyellow", col.ci = "lightblue")
  })

  # Axis customization
  vdiffr::expect_doppelganger("plot_custom_labels_base", function() {
    plotfunc(fit_brma, plot_type = "base",
             xlab = "Custom X", ylab = "Custom Y", main = "Custom Title")
  })

  # Axis limits
  vdiffr::expect_doppelganger("plot_custom_limits_base", function() {
    plotfunc(fit_brma, plot_type = "base", xlim = c(0, 60), ylim = c(-3, 1))
  })
})
```

### Testing Interface/Data Return

```r
test_that("Plot has correct interface", {

  fit_brma <- fits[["dataset_model"]]

  # Test as_data = TRUE returns correct structure
  plot_data <- plotfunc(fit_brma, as_data = TRUE)

  expect_true(is.list(plot_data),
    info = "as_data = TRUE should return a list")

  expected_components <- c("points", "xlim", "ylim", "xlab", "ylab")
  expect_true(all(expected_components %in% names(plot_data)),
    info = "plot data should contain expected components")

  # Test points structure
  expect_true(is.data.frame(plot_data$points))
  expect_true(all(c("x", "y") %in% names(plot_data$points)))

  # Test error handling
  expect_error(plotfunc(fit_brma, plot_type = "invalid"))
})
```

### Testing Edge Cases

```r
test_that("Plot handles edge cases", {

  # Model without moderators (if applicable)
  fit_simple <- fits[["simple_model"]]
  expect_error(regplot(fit_simple),
    regexp = "requires.*moderators")

  # Suppress warnings for expected conditions
  vdiffr::expect_doppelganger("plot_with_warnings",
    suppressWarnings(plotfunc(fit_brma, ...))
  )
})
```

### Naming Convention

Use consistent naming for doppelgangers:
- `{plottype}_{variant}_comparison` - side-by-side with reference
- `{plottype}_{variant}_brma_base` - brma base R standalone
- `{plottype}_{variant}_brma_ggplot` - brma ggplot2 standalone
- `{plottype}_custom_{aspect}_base` - customization tests

Examples:
- `funnel_simple_comparison`
- `regplot_continuous_brma_ggplot`
- `funnel_custom_points_base`

## Tips

1. **Use `par()` restoration** - Always restore par settings with `on.exit()`
2. **Match axis limits** - Use identical xlim/ylim for fair comparisons
3. **Skip unavailable fits** - Use `if (is.null(fit)) skip("...")`
4. **Suppress expected warnings** - Wrap with `suppressWarnings()` when needed
5. **Test both backends** - Always test base R and ggplot2 versions

## Running Tests

```r
# Generate/update reference snapshots
devtools::test(filter = "plotname")  # First run creates snapshots

# Accept new snapshots
vdiffr::manage_cases()  # Interactive review

# CI/CD will fail on differences
```

## File Organization

- Test files: `tests/testthat/test-02-{plotname}.R`
- Snapshots: `tests/testthat/_snaps/{testfile}/`
- Numbering (`test-02-*`) controls execution order after fit generation (`test-01-*`)
