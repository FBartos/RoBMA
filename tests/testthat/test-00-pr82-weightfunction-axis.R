test_that("observed weightfunction p-values follow the displayed axis", {

  make_object <- function(prior, direction = "positive") {
    data <- list(outcome = data.frame(yi = c(-2.5, 0, 2.5), sei = 1))
    attr(data, "effect_direction") <- direction
    list(data = data, priors = list(outcome = list(bias = prior)))
  }
  two_sided <- BayesTools::prior_weightfunction(
    "two-sided", 0.05, weights = BayesTools::wf_fixed(c(1, 0.5))
  )
  one_sided <- BayesTools::prior_weightfunction(
    "one-sided", 0.05, weights = BayesTools::wf_fixed(c(1, 0.5))
  )
  expected_two <- c(2 * stats::pnorm(-2.5), 1, 2 * stats::pnorm(-2.5))
  expect_equal(.weightfunction_observed_p_values(make_object(two_sided)),
               expected_two)
  expect_equal(.weightfunction_observed_p_values(make_object(two_sided, "negative")),
               expected_two)
  expect_equal(.weightfunction_observed_p_values(make_object(one_sided, "negative")),
               stats::pnorm(c(-2.5, 0, 2.5)))

  mixture <- BayesTools::prior_mixture(list(two_sided, one_sided))
  expect_equal(.weightfunction_observed_p_values(make_object(mixture)),
               stats::pnorm(c(2.5, 0, -2.5)))
})
