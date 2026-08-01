# Visual Testing with vdiffr

Reference this file when adding or updating plot tests in `tests/testthat/`.

## Current Pattern

Use `expect_vdiffr_snapshot()` from `tests/testthat/helper-visuals.R`.
It wraps `vdiffr::expect_doppelganger()` and skips cleanly when snapshots are not available.

```r
context("Plot function name")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)
```

## Side-by-Side Reference Plots

```r
test_that("Plot matches metafor structure", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  expect_vdiffr_snapshot("plotname_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::plotfunc(fit_metafor, main = "metafor")
    plotfunc(fit_brma, plot_type = "base", main = "brma")
  })

  expect_vdiffr_snapshot(
    "plotname_brma_ggplot",
    plotfunc(fit_brma, plot_type = "ggplot")
  )
})
```

## Extended Visual Coverage

Gallery-style or redundant visual cases should be gated:

```r
skip_if_not_full_visuals("This variant duplicates core visual coverage.")
```

Run `Rscript tools/test-profile.R certification` to include those cases.
Use `skip_if_not_full_visuals()` rather than a raw `testthat::skip()` for these
galleries. The helper announces existing conditional snapshots before skipping,
which prevents ordinary test runs from deleting their committed baselines.
The shared cached-fit skip helpers do the same when visual files cannot run
because their required fits are missing or stale.

## Data-Return Tests

Prefer stable structural checks over snapshots when testing `as_data = TRUE`:

```r
plot_data <- plotfunc(fit_brma, as_data = TRUE)
expect_true(is.list(plot_data))
expect_true(all(c("points", "xlim", "ylim", "xlab", "ylab") %in% names(plot_data)))
```

## Running Tests

Always use testthat LLM reporting:

```r
devtools::test(filter = "funnel|regplot|zplot", reporter = "llm")
```

Never auto-update committed snapshots.
Ask the maintainer to review visual differences.

## File Organization

- Plot tests: `tests/testthat/test-02-{plotname}.R` or `test-03-*` for integration-level plots such as `zplot`
- Snapshot helper: `tests/testthat/helper-visuals.R`
- Snapshots: `tests/testthat/_snaps/{testfile}/`
