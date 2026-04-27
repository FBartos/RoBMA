context("Summary")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "summary")

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)


test_that("Model summary", {

  ### default summary for all models
  for (name in names(fits)) {
    test_reference_table(summary(fits[[name]]), paste0("summary-", name, ".txt"), paste0("Summary table for '", name, "' mismatch"))
  }

  ### custom settings
  # custom probability & remove diagnostics
  test_reference_table(summary(fits[["bcg_meta-analysis"]], probs = c(0.01, 0.99), include_MCMC_diagnostics = FALSE),
                       paste0("summary-custom-prob-no-diagnostics.txt"), paste0("Summary table for '", name, "' mismatch"))

  # de-standardize coefficients
  test_reference_table(summary(fits[["bangertdrowns2004_location-scale"]], standardized_coefficients = TRUE),
                       paste0("summary-standardized.txt"), paste0("Summary table for '", name, "' mismatch"))

  # print output
  test_reference_table(fits[["bangertdrowns2004_location-scale"]],
                       paste0("summary-print.txt"), paste0("Summary table for '", name, "' mismatch"))

})
