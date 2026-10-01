test_that("random-effect prior constructors are exact BayesTools re-exports", {

  constructors <- c(
    "prior_random",
    "random_block",
    "random_variance_allocation",
    "allocation_ref",
    "random_covariance",
    "prior_lkj",
    "random_monitor",
    "random_new_levels",
    "random_sd_source"
  )

  expect_true(all(constructors %in% getNamespaceExports("RoBMA")))
  for (constructor in constructors) {
    expect_identical(
      getExportedValue("RoBMA", constructor),
      getExportedValue("BayesTools", constructor)
    )
  }
})
