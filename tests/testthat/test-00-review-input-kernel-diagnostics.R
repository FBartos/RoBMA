test_that("kernel-only bias priors do not request step-weight diagnostics", {

  priors <- list(outcome = list(bias = BayesTools::prior_phacking(form = "linear")))
  expect_identical(.convergence_structural_parameters(priors), character())
  priors$outcome$bias <- BayesTools::prior_weightfunction(
    "one-sided", .05, BayesTools::wf_fixed(c(1, .5)))
  expect_identical(.convergence_structural_parameters(priors),
                   c("omega[1]", "omega[0,0.05]"))
})
