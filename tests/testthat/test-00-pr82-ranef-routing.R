test_that("two-level total random effects keep the estimate names", {

  object <- brma(
    yi = c(-.2, .1, .4), sei = c(.2, .2, .2), measure = "OR",
    prior_unit_information_sd = 1, only_priors = TRUE
  )
  samples <- cbind(mu = c(-.2, 0, .2, .4), tau = c(.1, .2, .3, .4))

  estimate <- ranef(object, component = "estimate", .posterior_samples = samples)
  total    <- ranef(object, component = "total", .posterior_samples = samples)
  expect_identical(colnames(total), colnames(estimate))
  expect_false(anyDuplicated(colnames(total)) > 0L)
  expect_equal(unname(as.matrix(total)), unname(as.matrix(estimate)), tolerance = 0)

  # Expanded totals are unchanged.
  expanded <- ranef(object, component = "total", expand = TRUE, .posterior_samples = samples)
  expect_identical(colnames(expanded), colnames(estimate))
})
