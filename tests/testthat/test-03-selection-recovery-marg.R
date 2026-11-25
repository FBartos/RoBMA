context("Marginal selection simulation recovery")
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_certification("Selection recovery is reserved for its dedicated certification case.")
skip_if_not_installed("posterior")
skip_if_not_installed("metafor")

test_that("marginal selection recovers simulated fixed and nested parameters", {

  results <- .selection_recovery_run("marg")
  expect_equal(nrow(results), 12L)
})

test_that("marginal selection JAGS posterior matches independent dense integration", {

  .selection_recovery_oracle_run("marg", random = FALSE)
  .selection_recovery_oracle_run("marg", random = TRUE)
})
