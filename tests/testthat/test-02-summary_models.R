context("Summary Models")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "summary_models")

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)

summary_model_names <- intersect(
  names(fits),
  catalog_fits(class = c("BMA.norm", "BMA.glmm", "RoBMA"))
)


test_that("summary_models marginal summaries", {

  ### default marginal summary for RoBMA models
  for (name in names(fits)) {
    if (!name %in% summary_model_names) {
      next
    }
    test_reference_table(
      summary_models(fits[[name]], type = "marginal"),
      paste0("summary_models-marginal-", name, ".txt"),
      paste0("Marginal summary_models table for '", name, "' mismatch")
    )
  }
})

test_that("summary_models individual summaries", {

  ### default individual summary for RoBMA models
  for (name in names(fits)) {
    if (!name %in% summary_model_names) {
      next
    }
    test_reference_table(
      summary_models(fits[[name]], type = "individual"),
      paste0("summary_models-individual-", name, ".txt"),
      paste0("Individual summary_models table for '", name, "' mismatch")
    )
  }
})

test_that("summary_models is RoBMA-only", {

  skip_if_missing_fits("bcg_meta-analysis")

  expect_error(
    summary_models(fits[["bcg_meta-analysis"]]),
    "RoBMA objects"
  )
})

test_that("summary_models decodes interaction labels", {

  skip_if_missing_fits("dat.lehmann2018_RoBMA_mods2")

  marginal_output <- capture.output(print(summary_models(
    fits[["dat.lehmann2018_RoBMA_mods2"]],
    type = "marginal"
  )))
  individual_output <- capture.output(print(summary_models(
    fits[["dat.lehmann2018_RoBMA_mods2"]],
    type = "individual"
  )))

  expect_false(any(grepl("__xXx__", marginal_output, fixed = TRUE)))
  expect_false(any(grepl("__xXx__", individual_output, fixed = TRUE)))
  expect_true(any(grepl("Preregistered:Gender", marginal_output, fixed = TRUE)))
  expect_true(any(grepl("Preregistered:Gender", individual_output, fixed = TRUE)))
})
