test_that("normal LOO tails keep representable upper-tail probability", {

  yi         <- 10
  sei        <- 1
  mu_samples <- matrix(0, nrow = 1, ncol = 1)
  tau_within <- matrix(0, nrow = 1, ncol = 1)

  tails <- .loo_predictive_log_tails_estimate(list(outcome_type = "norm",
    data = structure(list(), outcome_type = "norm"), priors = list(), S = 1L, K = 1L,
    yi = yi, sei = sei, mu = mu_samples, tau_within = tau_within))
  cdf_vals <- exp(tails[["log_upper"]])

  expect_equal(cdf_vals[1, 1], stats::pnorm(yi, lower.tail = FALSE))
  expect_gt(cdf_vals[1, 1], 0)
})


test_that("production step CDF retains its upper-tail normal oracle", {

  yi         <- 3
  sei        <- 1
  mu_samples <- matrix(0, nrow = 1, ncol = 1)
  tau_within <- matrix(0, nrow = 1, ncol = 1)
  prior      <- BayesTools::prior_weightfunction(
    "one-sided", c(.05), BayesTools::wf_fixed(c(1, 1))
  )
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  selection_context <- c(spec, list(
    omega       = matrix(c(1, 1), nrow = 1),
    alpha       = 0,
    phack_kind  = 0L,
    kernel_mode = SELKERNEL_STEP
  ))

  cdf_vals <- .selection_step_cdf_matrix(
    q                 = yi,
    mean              = mu_samples,
    sd                = sqrt(tau_within^2 + sei^2),
    sei               = sei,
    selection_context = selection_context,
    lower.tail        = FALSE
  )

  expect_equal(cdf_vals[1, 1], stats::pnorm(yi, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_gt(cdf_vals[1, 1], 0)
})
