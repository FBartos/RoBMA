test_that("fixed selection weights use the complete mixed-prior p-value grid", {

  one <- BayesTools::prior_weightfunction(
    "one-sided", .05, BayesTools::wf_fixed(c(1, .5)))
  two <- BayesTools::prior_weightfunction(
    "two-sided", .05, BayesTools::wf_fixed(c(1, .2)))
  cuts <- c(0, .025, .05, .975, 1)
  expected <- rbind(c(1, 1, .5, .5), c(1, .2, .2, 1), rep(1, 4L))
  observed <- .selection_fixed_omega_by_branch(
    list(one, two, BayesTools::prior_none()), cuts)
  expect_equal(unname(observed), expected, tolerance = 0)
  specification <- .selection_spec(
    list(outcome = list(bias = BayesTools::prior_mixture(list(one, two)))),
    yi = c(-.2, .2), sei = c(.1, .1), effect_direction = "positive"
  )
  expect_equal(unname(specification[["fixed_omega"]]), expected[1:2, ],
               tolerance = 0)
  expect_error(.selection_fixed_omega_branch(one, c(0, 1)),
               "^Selection fixed weights require a compatible global p-value grid\\.$")
})
