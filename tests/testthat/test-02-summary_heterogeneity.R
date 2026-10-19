context("Summary heterogeneity")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

REFERENCE_DIR <<- testthat::test_path("..", "results", "summary_heterogeneity")

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(summary_heterogeneity_cases(), function(case) {
  test_that_case("Heterogeneity summary matches metafor", case, {
    expect_summary_heterogeneity_matches_metafor(case)
  })
})

test_that("Heterogeneity for model-averaged 2-level models is summarized", {

  model_averaged_classes <- c("BMA.norm", "BMA.glmm", "RoBMA")
  model_names <- setdiff(
    catalog_fits(class = model_averaged_classes),
    catalog_fits(class = model_averaged_classes, feature = "multilevel")
  )
  skip_if_missing_fits(model_names)

  expected_rows <- c("tau", "tau2", "I2", "H2")

  for (name in model_names) {
    expect_summary_heterogeneity_structure(
      summary_heterogeneity(fits[[name]]),
      expected_rows,
      name
    )
  }
})

test_that("Heterogeneity for model-averaged 3-level models is partitioned", {

  model_averaged_classes <- c("BMA.norm", "BMA.glmm", "RoBMA")
  model_names <- catalog_fits(
    class   = model_averaged_classes,
    feature = "multilevel"
  )
  skip_if_missing_fits(model_names)

  expected_rows <- c(
    "tau",
    "rho",
    "tau [within]",
    "tau [between]",
    "tau2",
    "tau2 [within]",
    "tau2 [between]",
    "I2",
    "I2 [within]",
    "I2 [between]",
    "H2"
  )

  for (name in model_names) {
    heterogeneity <- summary_heterogeneity(fits[[name]])
    expect_summary_heterogeneity_structure(heterogeneity, expected_rows, name)
    expect_equal(
      heterogeneity$estimates["tau2", "Mean"],
      heterogeneity$estimates["tau2 [within]", "Mean"] +
        heterogeneity$estimates["tau2 [between]", "Mean"],
      tolerance = 0.01,
      info      = paste0("tau2 partitions sum to total for '", name, "'")
    )
    expect_equal(
      heterogeneity$estimates["I2", "Mean"],
      heterogeneity$estimates["I2 [within]", "Mean"] +
        heterogeneity$estimates["I2 [between]", "Mean"],
      tolerance = 0.01,
      info      = paste0("I2 partitions sum to total for '", name, "'")
    )
  }
})

test_that("summary_heterogeneity print keeps leading blank line", {

  skip_if_missing_fits("bcg_meta-analysis")

  output <- capture.output(print(summary_heterogeneity(fits[["bcg_meta-analysis"]])))

  expect_equal(output[[1]], "")
})


test_that("ordinary heterogeneity exports retain their values and component", {

  skip_if_missing_fits("bcg_meta-analysis")
  result <- summary_heterogeneity(fits[["bcg_meta-analysis"]])
  frame <- as.data.frame(result)
  raw <- as.data.frame(result[["estimates"]])
  expect_identical(result[["component"]], "heterogeneity")
  expect_identical(frame[["component"]], rep("heterogeneity", nrow(raw)))
  expect_identical(frame[["parameter"]], rownames(raw))
  expect_identical(data.frame(result), frame)
  expect_equal(unname(as.matrix(frame[-c(1L, 2L)])),
    unname(as.matrix(raw)), tolerance = 0)
})
