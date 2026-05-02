context("Distribution helpers")
skip_on_cran()

test_that(".rowLogSumExps handles all -Inf rows", {

  x <- matrix(c(-Inf, -Inf, -1000, -1001), nrow = 2, byrow = TRUE)

  expect_equal(.rowLogSumExps(x)[1], -Inf)
  expect_equal(.rowLogSumExps(x)[2], -1000 + log1p(exp(-1)))
})

test_that("native coercion helpers preserve matrix and vector shapes", {

  x <- matrix(1:4, nrow = 2)

  expect_true(is.double(.native_numeric_matrix(x)))
  expect_equal(dim(.native_numeric_matrix(1:4)), c(4L, 1L))
  expect_true(is.double(.native_numeric_vector(1:4)))
  expect_true(is.integer(.native_integer_vector(c(1, 2))))
})

test_that("selected-normal step kernel handles boundary and extreme p-bins", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.001, .025, .05, .50),
    weights = BayesTools::wf_fixed(c(1, .8, .4, .2, .05))
  )
  sei <- c(.05, .10, .20, .40)
  yi  <- stats::qnorm(1 - c(.001, .025, .05, .50)) * sei
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )

  mu <- matrix(c(
    -.05, .00, .05, .10,
     .10, .05, .00, -.05
  ), nrow = 2, byrow = TRUE)
  sigma <- matrix(c(
    .06, .12, .25, .45,
    .08, .15, .30, .50
  ), nrow = 2, byrow = TRUE)
  omega <- matrix(c(
    1, .8, .4, .2, .05,
    1, .1, .9, .3, .01
  ), nrow = 2, byrow = TRUE)

  out <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )

  expect_true(all(is.finite(out)))
  expect_equal(dim(out), c(2L, 4L))
  expect_equal(out, .test_step_reference(yi, mu, sigma, sei, omega, spec),
               tolerance = 1e-12)
})

test_that("selected-normal log-likelihood matrix honors observation weights", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  yi <- c(.05, .15, .30)
  sei <- c(.10, .12, .14)
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  mu      <- matrix(c(.00, .10, .20, .03, .12, .24), nrow = 2, byrow = TRUE)
  sigma   <- matrix(.20, nrow = 2, ncol = 3)
  omega   <- matrix(c(1, .5, .25, 1, .7, .4), nrow = 2, byrow = TRUE)
  weights <- c(1, 2, .25)

  out_unweighted <- .selnorm_kernel_loglik_matrix(
    yi, mu, sigma, sei = sei, omega = omega, selection_spec = spec
  )
  out_weighted <- .selnorm_kernel_loglik_matrix(
    yi, mu, sigma, sei = sei, omega = omega, selection_spec = spec,
    weights = weights
  )

  expect_equal(out_weighted, sweep(out_unweighted, 2, weights, "*"),
               tolerance = 1e-12)
})

test_that("native selected-normal CDF and moments honor normal subdispatch", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.001, .025, .05, .50),
    weights = BayesTools::wf_fixed(c(1, .8, .4, .2, .05))
  )
  yi  <- c(.05, .15, .30, .60)
  sei <- c(.05, .10, .20, .40)
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )

  mu <- matrix(c(
    -.05, .00, .05, .10,
     .10, .05, .00, -.05,
     .20, .10, .05, .00
  ), nrow = 3, byrow = TRUE)
  sigma <- matrix(c(
    .06, .12, .25, .45,
    .08, .15, .30, .50,
    .10, .20, .35, .60
  ), nrow = 3, byrow = TRUE)
  omega <- matrix(c(
    2, 3, 4, 5, 6,
    1, .1, .9, .3, .01,
    1, .8, .4, .2, .05
  ), nrow = 3, byrow = TRUE)
  kernel_mode <- c(SELKERNEL_NORMAL, SELKERNEL_STEP, SELKERNEL_STEP)
  context <- spec
  context[["omega"]]       <- omega
  context[["alpha"]]       <- rep(0, 3)
  context[["phack_kind"]]  <- rep(0L, 3)
  context[["kernel_mode"]] <- kernel_mode
  context[["use_normal"]]  <- kernel_mode == SELKERNEL_NORMAL

  q <- c(-.10, .05, .20, .80)
  cdf_lower <- .selection_step_cdf_matrix(
    q                 = q,
    mean              = mu,
    sd                = sigma,
    sei               = sei,
    selection_context = context,
    lower.tail        = TRUE
  )
  cdf_upper <- .selection_step_cdf_matrix(
    q                 = q,
    mean              = mu,
    sd                = sigma,
    sei               = sei,
    selection_context = context,
    lower.tail        = FALSE
  )
  moments <- .selection_step_moments_matrix(
    mean              = mu,
    sd                = sigma,
    sei               = sei,
    selection_context = context
  )

  expect_equal(
    cdf_lower,
    .test_step_cdf_reference(q, mu, sigma, sei, omega, spec, kernel_mode),
    tolerance = 1e-12
  )
  expect_equal(cdf_lower + cdf_upper, matrix(1, nrow = 3, ncol = 4),
               tolerance = 1e-10)
  expect_equal(
    moments,
    .test_step_moments_reference(mu, sigma, sei, omega, spec, kernel_mode),
    tolerance = 1e-12
  )
  expect_equal(moments[["mean"]][1, ], mu[1, ], tolerance = 1e-12)
  expect_equal(moments[["second"]][1, ], sigma[1, ]^2 + mu[1, ]^2,
               tolerance = 1e-12)
  expect_equal(
    .selection_step_log_norm_matrix(mu, sigma, sei, context)[1, ],
    rep(0, ncol(mu)),
    tolerance = 1e-12
  )
})

test_that("native selected-normal CDF mirrors negative direction and preserves extreme upper tails", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  yi    <- c(.05, .20)
  sei   <- c(.10, .20)
  q     <- c(.10, .50)
  mu    <- matrix(c(.02, .25, .05, .30), nrow = 2, byrow = TRUE)
  sigma <- matrix(c(.15, .30, .18, .35), nrow = 2, byrow = TRUE)
  omega <- matrix(c(1, .5, .25, 1, .7, .4), nrow = 2, byrow = TRUE)

  spec_pos <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  context_pos <- spec_pos
  context_pos[["omega"]]       <- omega
  context_pos[["alpha"]]       <- rep(0, 2)
  context_pos[["phack_kind"]]  <- rep(0L, 2)
  context_pos[["kernel_mode"]] <- rep(SELKERNEL_STEP, 2)

  spec_neg <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = -yi,
    sei              = sei,
    effect_direction = "negative"
  )
  context_neg <- spec_neg
  context_neg[["omega"]]       <- omega
  context_neg[["alpha"]]       <- rep(0, 2)
  context_neg[["phack_kind"]]  <- rep(0L, 2)
  context_neg[["kernel_mode"]] <- rep(SELKERNEL_STEP, 2)

  pos_upper <- .selection_step_cdf_matrix(
    q, mu, sigma, sei, context_pos, lower.tail = FALSE
  )
  neg_lower <- .selection_step_cdf_matrix(
    -q, -mu, sigma, sei, context_neg, lower.tail = TRUE
  )
  expect_equal(neg_lower, pos_upper, tolerance = 1e-12)

  context_extreme <- spec_pos
  context_extreme[["omega"]]       <- matrix(1, nrow = 1, ncol = spec_pos[["n_bins"]])
  context_extreme[["alpha"]]       <- 0
  context_extreme[["phack_kind"]]  <- 0L
  context_extreme[["kernel_mode"]] <- SELKERNEL_STEP

  extreme <- .selection_step_cdf_matrix(
    q                 = 10,
    mean              = matrix(0, nrow = 1, ncol = 1),
    sd                = matrix(1, nrow = 1, ncol = 1),
    sei               = 1,
    selection_context = context_extreme,
    lower.tail        = FALSE
  )
  expect_gt(extreme[1, 1], 0)
  expect_equal(extreme[1, 1], stats::pnorm(10, lower.tail = FALSE),
               tolerance = 1e-20)
})
