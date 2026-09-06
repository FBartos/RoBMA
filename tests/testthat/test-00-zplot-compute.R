test_that("zplot total SD matches replicated-SE evaluation", {

  tau_within <- matrix(
    c(0, .1, 1, 10, .25, .5, 2, 20, .75, 1.5, 3, 30),
    nrow = 4,
    ncol = 3
  )
  sei     <- c(.05, .2, 2)
  sei_mat <- matrix(sei, nrow = nrow(tau_within), ncol = ncol(tau_within),
                    byrow = TRUE)
  expected <- .root_sum_squares(tau_within, sei_mat)

  expect_identical(.zplot_total_sd(tau_within, sei), expected)
  expect_equal(
    .zplot_total_sd(tau_within, sei),
    sqrt(sweep(tau_within^2, 2, sei^2, "+")),
    tolerance = 1e-15
  )
})


test_that("joint z marginals match an independently integrated bivariate law", {

  cutoff <- stats::qnorm(.975)
  sei <- c(.7, 1.3)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), sei,
                          effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2, 1, .2), 2, 2, byrow = TRUE)
  spec$alpha       <- c(0, 0)
  spec$phack_kind  <- c(0L, 0L)
  spec$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_NORMAL)
  spec$use_normal  <- c(FALSE, TRUE)
  rho <- .7
  covariance <- matrix(rep(c(sei[1]^2, prod(sei) * rho, sei[2]^2), each = 2), 2)
  weight <- function(z) ifelse(z > cutoff, 1, .2)
  other_weight <- function(z) {
    .2 + .8 * stats::pnorm(cutoff, rho * z, sqrt(1 - rho^2), lower.tail = FALSE)
  }
  numerator <- function(z) stats::dnorm(z) * weight(z) * other_weight(z)
  normalizer <- stats::integrate(numerator, -Inf, cutoff, rel.tol = 1e-10)$value +
    stats::integrate(numerator, cutoff, Inf, rel.tol = 1e-10)$value
  z <- c(-3, -.5, 0, 1, 2.5, 4)
  control <- set_selection_likelihood_control(points_per_scramble = 256,
    max_points_per_scramble = 16384, relative_tolerance = .002)
  result <- .zplot_joint_block(z, matrix(0, 2, 2),
    covariance, sei, spec, FALSE, control)
  expect_equal(result$density[1, ], numerator(z) / normalizer, tolerance = .002)
  expect_equal(result$density[2, ], stats::dnorm(z), tolerance = 1e-12)
  expect_equal(exp(-result$log_density), c(normalizer, 1), tolerance = .002)
  factor <- .zplot_joint_block(z, matrix(0, 2, 2),
    covariance, sei, spec, FALSE, control,
    factors = list(residual_sd = matrix(rep(sqrt(1 - rho) * sei, each = 2), 2),
                   loading = matrix(rep(sqrt(rho) * sei, each = 2), 2)))
  expect_equal(factor$density[1, ], numerator(z) / normalizer, tolerance = .002)
  probability <- .zplot_joint_block(cutoff, matrix(0, 2, 2),
    covariance, sei, spec, TRUE, control)
  expected_p <- (stats::integrate(numerator, -Inf, -cutoff)$value +
    stats::integrate(numerator, cutoff, Inf)$value) / normalizer
  expect_equal(probability$density[, 1], c(expected_p, .05), tolerance = .002)
  whole_line <- .zplot_joint_block(0, matrix(.3, 2, 2),
    covariance, sei, spec, TRUE, control)
  expect_equal(as.numeric(whole_line$density), c(1, 1), tolerance = 1e-12)
  # The extrapolated density has area 1/C, so its missing mass is 1/C - 1.
  expect_gt(1 / normalizer - 1, 0)
  spec$sign <- -1L
  reflected <- .zplot_joint_block(-z, matrix(0, 2, 2),
    covariance, sei, spec, FALSE, control)
  expect_equal(reflected$density[1, ], numerator(z) / normalizer, tolerance = .002)
})


test_that("approximate z marginals integrate conditional selection normalizers", {

  cutoff <- stats::qnorm(.975)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), 0, 1,
                          effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$use_normal  <- FALSE
  z <- c(-3, 0, 1, 2.5, 4)
  normalizer <- function(u) .2 + .8 * stats::pnorm(cutoff, .3 + .8 * u, 1,
                                                 lower.tail = FALSE)
  reference <- vapply(z, function(value) {
    stats::integrate(function(u) {
      stats::dnorm(u) * stats::dnorm(value, .3 + .8 * u, 1) / normalizer(u)
    }, -Inf, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  result <- .zplot_latent_mixture(z, matrix(.3), matrix(1), matrix(.8), 1,
    spec, FALSE, set_selection_likelihood_control(relative_tolerance = 1e-7))
  expect_equal(as.numeric(result$extrapolated), reference, tolerance = 1e-7)
  expect_equal(as.numeric(result$fitted), reference * ifelse(z > cutoff, 1, .2),
               tolerance = 1e-7)
  threshold <- .zplot_latent_mixture(cutoff, matrix(.3), matrix(1), matrix(.8), 1,
    spec, TRUE, set_selection_likelihood_control(relative_tolerance = 1e-7))
  expected_area <- stats::integrate(function(u) stats::dnorm(u) / normalizer(u),
                                   -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(threshold$weights, expected_area, tolerance = 1e-7)
  expect_equal(result$weights, expected_area, tolerance = 1e-7)
  whole_line <- .zplot_latent_mixture(0, matrix(.3), matrix(1), matrix(.8), 1,
    spec, TRUE, set_selection_likelihood_control(relative_tolerance = 1e-7))
  expect_equal(as.numeric(whole_line$fitted), 1, tolerance = 1e-12)
  expect_equal(whole_line$EDR, 1, tolerance = 1e-12)
})


test_that("factor importance integration agrees with a one-factor reference", {

  cutoff <- stats::qnorm(.975)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .05))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), rep(0, 3), rep(1, 3),
    effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .05), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$use_normal  <- FALSE
  loading <- matrix(0, 3, 3)
  loading[, 1] <- sqrt(.7)
  covariance <- diag(.3, 3) + tcrossprod(loading)
  normalizer <- function(u) .05 + .95 * stats::pnorm(cutoff, sqrt(.7) * u,
    sqrt(.3), lower.tail = FALSE)
  area <- stats::integrate(function(u) stats::dnorm(u) * normalizer(u)^3,
    -Inf, Inf, rel.tol = 1e-10)$value
  z <- c(-1, 0, 2.5, 4)
  reference <- vapply(z, function(value) {
    stats::integrate(function(u) stats::dnorm(u) * normalizer(u)^2 *
      stats::dnorm(value, sqrt(.7) * u, sqrt(.3)), -Inf, Inf,
      rel.tol = 1e-10)$value * ifelse(value > cutoff, 1, .05) / area
  }, numeric(1))
  result <- .zplot_joint_block(z, matrix(0, 1, 3),
    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1), rep(1, 3),
    spec, FALSE, set_selection_likelihood_control(relative_tolerance = .002,
      max_points_per_scramble = 32768),
    factors = list(residual_sd = matrix(sqrt(.3), 1, 3), loading = matrix(loading, 1)))
  expect_equal(as.numeric(result$density), reference, tolerance = .002)
  expect_equal(exp(-result$log_density), area, tolerance = .002)
})


test_that("zplot reports unavailable targets and failed integration explicitly", {

  control <- set_selection_likelihood_control(points_per_scramble = 128,
    max_points_per_scramble = 128, relative_tolerance = 1e-12)
  testthat::local_mocked_bindings(
    .selection_context = function(...) list(),
    .selection_require_step_evaluable = function(...) NULL,
    .package = "RoBMA"
  )
  expect_error(
    .zplot_selection_marginal(NULL, matrix(0), 0, NULL, "estimate", control),
    'Conditional zplot selection targets are unavailable. Use \'conditioning_depth = "marginal"\'.',
    fixed = TRUE
  )

  prior <- BayesTools::prior_weightfunction(
    "one-sided", .025, BayesTools::wf_fixed(c(1, .2))
  )
  spec <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), c(1, 1),
    effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$use_normal  <- FALSE
  failure <- tryCatch(
    .zplot_joint_block(0, matrix(0, 1, 2), matrix(c(1, .7, 1), 1),
      c(1, 1), spec, FALSE, control),
    error = identity
  )
  expect_s3_class(failure, "error")
  message <- conditionMessage(failure)
  observed <- sub("^.*error was ([^ ]+)\\. Increase.*$", "\\1", message)
  expect_gt(as.numeric(observed), control$relative_tolerance)
  expect_identical(message,
    paste0("Zplot marginal integration was rejected by diagnostics: ",
      "relative integration error was ", observed,
      ". Increase 'max_points_per_scramble' ",
      "in 'integration_control = set_selection_likelihood_control()'.")
  )
})


test_that("zplot reuses invariant predictive components", {

  calls <- new.env(parent = emptyenv())
  calls$predictive <- 0L
  calls$density    <- 0L
  calls$paired     <- 0L
  calls$selection  <- FALSE

  predictive <- list(
    mu         = matrix(c(.1, .2, .3, .4), nrow = 2),
    tau_within = matrix(c(.05, .1, .15, .2), nrow = 2),
    sei        = c(.1, .2)
  )
  object <- list(fit = structure(list(), class = "BayesTools_fit"))

  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, nrow = 2, ncol = 1),
    .thin_sample_rows = function(...) NULL,
    .is_PET = function(...) FALSE,
    .is_PEESE = function(...) FALSE,
    .is_weightfunction = function(...) calls$selection,
    .effect_direction = function(...) "positive",
    .zplot_predictive_components = function(..., extrapolate) {
      calls$predictive <- calls$predictive + 1L
      if (extrapolate) {
        stop("Invariant predictive components were recomputed.")
      }
      predictive
    },
    .zplot_selection_context = function(...) {
      if (calls$selection) list(selection = TRUE) else NULL
    },
    .zplot_density_vectorized = function(...) {
      calls$density <- calls$density + 1L
      matrix(1, nrow = 1, ncol = 1)
    },
    .zplot_selnorm_density_pair = function(...) {
      calls$paired <- calls$paired + 1L
      list(
        fitted       = matrix(2, nrow = 1, ncol = 1),
        extrapolated = matrix(3, nrow = 1, ncol = 1)
      )
    },
    .package = "RoBMA"
  )

  ordinary <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 1L)
  expect_identical(calls$density, 1L)
  expect_identical(calls$paired, 0L)
  expect_identical(ordinary$fitted, ordinary$extrapolated)

  calls$predictive <- 0L
  calls$density    <- 0L
  calls$paired     <- 0L
  calls$selection  <- TRUE
  selected <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 1L)
  expect_identical(calls$density, 0L)
  expect_identical(calls$paired, 1L)
  expect_false(identical(selected$fitted, selected$extrapolated))
})


test_that("zplot retains separate PET and PEESE predictive components", {

  calls <- new.env(parent = emptyenv())
  calls$predictive <- 0L
  calls$density    <- 0L

  object <- list(fit = structure(list(), class = "BayesTools_fit"))
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, nrow = 2, ncol = 1),
    .thin_sample_rows = function(...) NULL,
    .is_PET = function(...) TRUE,
    .is_PEESE = function(...) FALSE,
    .is_weightfunction = function(...) FALSE,
    .effect_direction = function(...) "positive",
    .zplot_predictive_components = function(..., extrapolate) {
      calls$predictive <- calls$predictive + 1L
      list(
        mu         = matrix(as.numeric(extrapolate), nrow = 2, ncol = 1),
        tau_within = matrix(.1, nrow = 2, ncol = 1),
        sei        = .2
      )
    },
    .zplot_selection_context = function(...) NULL,
    .zplot_density_vectorized = function(..., mu_samples) {
      calls$density <- calls$density + 1L
      mu_samples
    },
    .package = "RoBMA"
  )

  result <- .zplot_density_pair(object, z_sequence = 0, max_samples = 10)
  expect_identical(calls$predictive, 2L)
  expect_identical(calls$density, 2L)
  expect_equal(result$fitted, matrix(0, nrow = 2, ncol = 1))
  expect_equal(result$extrapolated, matrix(1, nrow = 2, ncol = 1))
})
