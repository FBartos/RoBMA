.postfit_cdf_fixture <- function(weights = c(1, .2), direction = "positive") {

  sei <- c(.7, .9)
  covariance <- matrix(c(1, .4, .4, .9), 2L)
  sign <- if (direction == "positive") 1 else -1
  mean <- matrix(sign * c(.1, -.2), 1L)
  yi <- sign * c(.2, -.1)
  steps <- if (length(weights) == 2L) .025 else c(.025, .17)
  prior <- BayesTools::prior_weightfunction("one-sided", steps, BayesTools::wf_fixed(weights))
  context <- .selection_spec(list(outcome = list(bias = prior)), yi, sei,
    effect_direction = direction, signed_data = FALSE)
  context$omega <- matrix(weights, 1L)
  context$alpha <- 0
  context$phack_kind <- 0L
  context$kernel_mode <- 1L
  context$vector_rule <- 0L
  context$use_normal <- FALSE
  control <- set_selection_likelihood_control(relative_tolerance = 1e-6)
  plan <- c(unclass(control), list(
    quadrature = .selection_joint_cluster_quadrature_rules(c(7L, 15L, 31L, 63L)),
    designs = list(`2` = BayesTools::selection_qmc_design(4L, control$points_per_scramble,
      control$scrambles, control$seed))))
  list(yi = yi, mean = mean, covariance = covariance,
    packed = matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L),
    sei = sei, context = context, plan = plan)
}

.postfit_cdf_original <- function(x) {

  static <- BayesTools::selection_native_static_args(x$context)
  .Call("RoBMA_selnorm_mnorm_step_loglik_batch", as.double(x$yi), x$mean,
    x$packed, x$sei, x$context$omega, static$z_lower, static$z_upper,
    as.integer(x$context$obs_bin), static$sign, static$telescope_probabilities,
    x$context$kernel_mode, as.double(x$plan$designs[["2"]]),
    as.integer(x$plan$points_per_scramble), as.integer(x$plan$scrambles),
    as.double(x$plan$relative_tolerance), TRUE, x$context$vector_rule,
    x$plan$quadrature, PACKAGE = "RoBMA")
}

test_that("postfit CDF routing preserves the independent Gaussian event and fitted plan", {

  # Independent bivariate conditioning integral, split at the first observed
  # weight discontinuity. It does not call a package selection normalizer.
  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  split <- cutoff * .7 - .1
  weighted_second <- function(u) {
    stats::dnorm(u) * (.2 + .8 * stats::pnorm(cutoff * .9,
      -.2 + .4 * u, sqrt(.9 - .4^2), lower.tail = FALSE))
  }
  expected_mass <- .2 * stats::integrate(weighted_second, -Inf, split, rel.tol = 1e-11)$value +
    stats::integrate(weighted_second, split, Inf, rel.tol = 1e-11)$value
  for (direction in c("positive", "negative")) {
    x <- .postfit_cdf_fixture(direction = direction)
    before <- serialize(x$plan, NULL)
    density <- .selection_joint_dense_loglik_block(x$yi, x$mean, x$packed,
      x$sei, x$context, x$plan, 2L, return_normalizer = TRUE)
    mass <- .selection_gaussian_event_mass(x$mean, x$packed, x$sei, x$context, x$plan)
    residual <- forwardsolve(t(chol(x$covariance)), x$yi - as.numeric(x$mean))
    gaussian <- -log(2 * pi) - sum(log(diag(chol(x$covariance)))) - sum(residual^2) / 2
    expect_equal(as.numeric(density$log_normalizer), log(expected_mass), tolerance = 2e-7)
    expect_equal(as.numeric(mass$log_mass), log(expected_mass), tolerance = 2e-7)
    expect_equal(as.numeric(density$log_density), gaussian + 2 * log(.2) - log(expected_mass), tolerance = 2e-7)
    cdf <- attr(density$integration_diagnostics, "cdf_relative_error", exact = TRUE)
    expect_true(length(cdf) == 1L && is.finite(cdf) && cdf > 0)
    expect_true(attr(mass$relative_quadrature_error, "cdf_relative_error", exact = TRUE) > 0)
    total <- density$integration_diagnostics[, "covariance_width"] +
      2 * density$integration_diagnostics[, "quadrature_change"] +
      density$integration_diagnostics[, "tail_bound"] + cdf
    expect_lte(as.numeric(total), x$plan$relative_tolerance)
    expect_lte(as.numeric(mass$relative_quadrature_error), x$plan$relative_tolerance)
    expect_identical(serialize(x$plan, NULL), before)
    expect_null(attr(x$plan$quadrature, "bounded_cdf", exact = TRUE))
  }
})

test_that("bounded postfit calls cannot populate or consume the ordinary cache", {

  old <- .Call("RoBMA_selnorm_cache_control", NULL, 0L, PACKAGE = "RoBMA")
  withr::defer(invisible(.Call("RoBMA_selnorm_cache_control",
    as.double(old$capacity_bytes), 3L, PACKAGE = "RoBMA")))
  invisible(.Call("RoBMA_selnorm_cache_control", as.double(4 * 1024^2), 3L, PACKAGE = "RoBMA"))
  x <- .postfit_cdf_fixture()
  exact <- .postfit_cdf_original(x)
  expect_null(attr(exact$integration_diagnostics, "cdf_relative_error", exact = TRUE))
  before <- .Call("RoBMA_selnorm_cache_control", NULL, 0L, PACKAGE = "RoBMA")
  postfit <- .selection_joint_dense_loglik_block(x$yi, x$mean, x$packed,
    x$sei, x$context, x$plan, 2L, return_normalizer = TRUE)
  after <- .Call("RoBMA_selnorm_cache_control", NULL, 0L, PACKAGE = "RoBMA")
  expect_identical(after, before)
  expect_gt(attr(postfit$integration_diagnostics, "cdf_relative_error", exact = TRUE), 0)
  again <- .postfit_cdf_original(x)
  expect_identical(again$log_density, exact$log_density)
  expect_identical(again$log_normalizer, exact$log_normalizer)
  expect_gt(.Call("RoBMA_selnorm_cache_control", NULL, 0L, PACKAGE = "RoBMA")$exact$hits, before$exact$hits)

  for (weights in list(c(1, 0), c(1, 1e-12), c(1, .2, .2))) {
    x <- .postfit_cdf_fixture(weights)
    exact <- .postfit_cdf_original(x)
    postfit <- .selection_joint_dense_loglik_block(x$yi, x$mean, x$packed,
      x$sei, x$context, x$plan, 2L, return_normalizer = TRUE)
    expect_identical(postfit$log_normalizer, exact$log_normalizer)
    expect_identical(postfit$log_density, exact$log_density)
    expect_identical(attr(postfit$integration_diagnostics, "cdf_relative_error", exact = TRUE), 0)
  }
})

test_that("an eligible bounded envelope retries the ordinary rule before QMC", {

  x <- .postfit_cdf_fixture()
  x$yi <- c(0, 0)
  x$mean <- matrix(c(0, 0), 1L)
  x$sei <- c(1, 1)
  x$packed <- matrix(c(1, .7, 1), 1L)
  x$plan$quadrature <- .selection_joint_cluster_quadrature_rules(c(7L, 15L, 31L))
  invoke <- function(tolerance, bounded) {
    x$plan$relative_tolerance <- tolerance
    if (!bounded) return(.postfit_cdf_original(x))
    .selection_joint_dense_loglik_block(x$yi, x$mean, x$packed,
      x$sei, x$context, x$plan, 2L, return_normalizer = TRUE)
  }
  score <- function(value) {
    diagnostic <- value$integration_diagnostics
    cdf <- attr(diagnostic, "cdf_relative_error", exact = TRUE)
    sum(diagnostic[, c("covariance_width", "tail_bound")]) +
      2 * diagnostic[, "quadrature_change"] + if (is.null(cdf)) 0 else cdf
  }

  # Both pilots must pass the same final rule: the pilot budget is below
  # both coarse scores. Choose an interior tolerance between the measured
  # final scores instead of hard-coding a platform-sensitive boundary.
  coarse <- c(score(invoke(.5, FALSE)), score(invoke(.5, TRUE)))
  pilot_tolerance <- .8 * min(coarse)
  ordinary <- invoke(pilot_tolerance, FALSE)
  bounded <- invoke(pilot_tolerance, TRUE)
  ordinary_score <- score(ordinary)
  bounded_score <- score(bounded)
  log_error <- attr(bounded$integration_diagnostics, "cdf_log_error", exact = TRUE)
  expect_gt(log_error, 0)
  expect_equal(ordinary$integration_diagnostics[, "used_covariance_envelope"][[1L]], 1)
  expect_equal(bounded$integration_diagnostics[, "used_covariance_envelope"][[1L]], 1)
  tolerance <- (ordinary_score + bounded_score) / 2
  expect_gt(tolerance - ordinary_score, 1000 * .Machine$double.eps)
  expect_gt(bounded_score - tolerance, 1000 * .Machine$double.eps)
  expect_lt(tolerance, min(coarse))

  # Eligibility depends on these unchanged weights/K and this CDF-only
  # floor. The final call cannot appear to retry by disabling the table at
  # preflight because its entire approximation budget was already exhausted.
  eligibility_floor <- 3 * exp(log_error) * expm1(2 * log_error) + expm1(log_error)
  expect_gt(tolerance, 100 * eligibility_floor)
  exact <- invoke(tolerance, FALSE)
  retried <- invoke(tolerance, TRUE)
  expect_equal(retried$integration_diagnostics[, "used_covariance_envelope"][[1L]], 1)
  expect_identical(attr(retried$integration_diagnostics, "cdf_relative_error", exact = TRUE), 0)
  expect_identical(retried$log_normalizer, exact$log_normalizer)
  expect_identical(retried$log_density, exact$log_density)

  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  conditional <- function(u) {
    stats::dnorm(u) * (.2 + .8 * stats::pnorm(cutoff,
      .7 * u, sqrt(1 - .7^2), lower.tail = FALSE))
  }
  reference <- .2 * stats::integrate(conditional, -Inf, cutoff, rel.tol = 1e-11)$value +
    stats::integrate(conditional, cutoff, Inf, rel.tol = 1e-11)$value
  expect_lte(abs(expm1(as.numeric(retried$log_normalizer) - log(reference))), 1e-6)
})

test_that("missing native CDF diagnostics cannot silently bypass the postfit contract", {

  x <- .postfit_cdf_fixture()
  native <- .postfit_cdf_original(x)
  guarded <- .selection_joint_dense_loglik_block
  environment(guarded) <- list2env(list(.Call = function(...) native),
    parent = environment(.selection_joint_dense_loglik_block))
  expect_error(guarded(x$yi, x$mean, x$packed, x$sei, x$context, x$plan, 2L, TRUE),
    "Selection CDF integration diagnostics are unavailable.", fixed = TRUE)
})

test_that("dense selection batches preserve independent row results and diagnostics", {

  x <- .postfit_cdf_fixture()
  S <- 19L
  means <- x$mean[rep(1L, S), , drop = FALSE]
  covariance <- x$packed[rep(1L, S), , drop = FALSE]
  omega <- x$context$omega[rep(1L, S), , drop = FALSE]
  kernel <- rep(1L, S)
  rule <- rep(0L, S)
  # Adjacent repeats, each varying-input family, a repeated failure, then the
  # original state again. No fit, altered integration budget, or cache counter.
  means[2:3, 1L] <- means[2:3, 1L] + .11
  covariance[5:6, 2L] <- .1
  omega[8:9, 2L] <- .6
  kernel[11:12] <- 0L
  rule[14:15] <- 1L
  covariance[17:18, 1L] <- -1
  static <- BayesTools::selection_native_static_args(x$context)
  quadrature <- x$plan$quadrature
  attr(quadrature, "bounded_cdf") <- TRUE
  invoke <- function(rows, return_normalizer) {

    .Call("RoBMA_selnorm_mnorm_step_loglik_batch", as.double(x$yi),
      means[rows, , drop = FALSE], covariance[rows, , drop = FALSE], x$sei,
      omega[rows, , drop = FALSE], static$z_lower, static$z_upper,
      as.integer(x$context$obs_bin), static$sign, static$telescope_probabilities,
      kernel[rows], as.double(x$plan$designs[["2"]]),
      as.integer(x$plan$points_per_scramble), as.integer(x$plan$scrambles),
      as.double(x$plan$relative_tolerance), return_normalizer, rule[rows],
      quadrature, PACKAGE = "RoBMA")
  }
  for (return_normalizer in c(FALSE, TRUE)) {
    batch <- invoke(seq_len(S), return_normalizer)
    for (s in seq_len(S)) {
      individual <- invoke(s, return_normalizer)
      row <- batch
      for (name in c("log_density", "relative_mcse", "log_normalizer")) {
        row[[name]] <- row[[name]][s]
      }
      row$integration_diagnostics <- batch$integration_diagnostics[s, , drop = FALSE]
      for (name in c("cdf_relative_error", "cdf_log_error")) {
        attr(row$integration_diagnostics, name) <-
          attr(batch$integration_diagnostics, name, exact = TRUE)[s]
      }
      expect_identical(row, individual)
    }
  }
})
