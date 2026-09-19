test_that("retained sampling does not enter the scalar marginal selection CDF", {

  # With tau = 0, conditioned sampling, and positive selection weights, the
  # weights cancel and the marginal law is ordinary normal. Folding the
  # sampling variance into the scalar selected CDF incorrectly tilts that law.
  object <- bselmodel(
    yi = c(.2, .5), sei = c(.2, .3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    selection = BayesTools::selection_model(known_sampling_variance = "condition")
  )
  expected <- matrix(c(.3, .7), 1L)
  testthat::local_mocked_bindings(
    .cdf_lik_estimate.brma = function(object) expected,
    .package = "RoBMA"
  )
  expect_equal(unname(.cdf.brma(object, conditioning_depth = "estimate")), expected)
  expect_error(
    .cdf.brma(object, conditioning_depth = "marginal"),
    "Joint-selection CDF evaluation is unavailable at this conditioning depth.",
    fixed = TRUE
  )
})
