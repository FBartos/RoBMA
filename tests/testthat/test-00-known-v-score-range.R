.known_v_score_range_values <- function(y, covariance, blocks) {

  for (symbol in c("RoBMA_known_v_covariance_plan_create",
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    "RoBMA_known_v_covariance_plan_conditional_loglik_batch",
    "RoBMA_known_v_covariance_plan_conditional_summary_batch",
    "RoBMA_known_v_covariance_plan_precision_residual_batch")) {
    expect_s3_class(getNativeSymbolInfo(symbol, PACKAGE = "RoBMA"),
      "NativeSymbolInfo")
  }
  K <- length(y)
  plan <- .Call("RoBMA_known_v_covariance_plan_create", as.numeric(y),
    covariance, list(), blocks, PACKAGE = "RoBMA")
  mean <- numeric(K)
  extra <- numeric(K)
  list(
    scalar = .Call("RoBMA_known_v_covariance_plan_conditional_loglik",
      plan, mean, list(), extra, PACKAGE = "RoBMA"),
    batch = .Call("RoBMA_known_v_covariance_plan_conditional_loglik_batch",
      plan, matrix(mean, K), list(list()), matrix(extra, K), PACKAGE = "RoBMA"),
    summary = .Call("RoBMA_known_v_covariance_plan_conditional_summary_batch",
      plan, matrix(mean, K), list(list()), matrix(extra, K), PACKAGE = "RoBMA"),
    precision_residual = .Call("RoBMA_known_v_covariance_plan_precision_residual_batch",
      plan, matrix(mean, K), list(list()), matrix(extra, K), PACKAGE = "RoBMA")
  )
}

test_that("conditional scores standardize finite precision residuals before squaring", {

  covariance <- 1e-200 * matrix(c(1, .5, .5, 1), 2L)
  actual <- .known_v_score_range_values(c(1, 0), covariance, list(1:2))
  residual <- c(1, -.5)
  variance <- rep(.75e-200, 2L)
  reference <- stats::dnorm(residual, sd = sqrt(variance), log = TRUE)
  expect_true(all(is.finite(actual$scalar)))
  expect_lte(max(abs(actual$scalar / reference - 1)), 1e-12)
  expect_equal(as.numeric(actual$batch), actual$scalar, tolerance = 0)
  expect_equal(as.numeric(actual$summary$residual), residual, tolerance = 1e-14)
  expect_lte(max(abs(as.numeric(actual$summary$variance) / variance - 1)), 1e-14)
  expected_precision_residual <- residual / variance
  expect_lte(max(abs(as.numeric(actual$precision_residual) /
    expected_precision_residual - 1)), 1e-14)
})

test_that("singleton scores separate Gaussian constants at huge finite variance", {

  actual <- .known_v_score_range_values(0, matrix(1e308, 1L), list(1L))
  reference <- stats::dnorm(0, sd = 1e154, log = TRUE)
  expect_true(is.finite(actual$scalar))
  expect_equal(actual$scalar, reference, tolerance = 1e-14)
  expect_equal(as.numeric(actual$batch), reference, tolerance = 1e-14)
  expect_identical(as.numeric(actual$summary$residual), 0)
  expect_identical(as.numeric(actual$summary$variance), 1e308)
  expect_identical(as.numeric(actual$precision_residual), 0)
})

test_that("finite Gaussian half squares survive an overflowing unhalved square", {

  z <- 1.5e154
  # Explicit left association keeps the representable half square finite.
  reference <- -.5 * log(2 * pi) - (.5 * z) * z
  expect_true(is.finite(reference))
  expect_equal(reference, stats::dnorm(z, log = TRUE), tolerance = 1e-14)
  actual <- .known_v_score_range_values(z, matrix(1, 1L), list(1L))
  expect_true(is.finite(actual$scalar))
  expect_lte(abs(actual$scalar / reference - 1), 1e-14)
  expect_equal(as.numeric(actual$batch), actual$scalar, tolerance = 0)
  expect_identical(as.numeric(actual$summary$residual), z)
  expect_identical(as.numeric(actual$summary$variance), 1)
  expect_identical(as.numeric(actual$precision_residual), z)

  overflow <- .known_v_score_range_values(2e154, matrix(1, 1L), list(1L))
  expect_identical(overflow$scalar, -Inf)
})

test_that("rank-zero dense covariance scores and ordinary controls preserve their laws", {

  # Factory rank zero cannot enter diagonal_low_rank_conditional_values;
  # this admitted multirow plan exercises its dense covariance fallback.
  actual <- .known_v_score_range_values(c(0, 0), diag(1e308, 2L), list(1:2))
  reference <- rep(stats::dnorm(0, sd = 1e154, log = TRUE), 2L)
  expect_equal(actual$scalar, reference, tolerance = 1e-14)
  expect_equal(as.numeric(actual$summary$variance), rep(1e308, 2L), tolerance = 1e-14)

  covariance <- matrix(c(1, .5, .5, 1), 2L)
  actual <- .known_v_score_range_values(c(1, 0), covariance, list(1:2))
  expect_equal(actual$scalar, stats::dnorm(c(1, -.5), sd = sqrt(.75), log = TRUE),
    tolerance = 1e-14)
  expect_equal(as.numeric(actual$summary$residual), c(1, -.5), tolerance = 1e-14)
  expect_equal(as.numeric(actual$summary$variance), rep(.75, 2L), tolerance = 1e-14)
})


test_that("joint R singleton scores preserve finite Gaussian log densities", {

  cases <- list(c(1e308, 0), c(1e300, 1e155), c(1, 1.5e154), c(2, .5), c(1, 2e154))
  for (mixed in c(FALSE, TRUE)) {
    for (case in cases) {
      variance <- case[[1L]]
      residual <- case[[2L]]
      V <- diag(c(1, 1, variance))
      if (mixed) V[1L, 2L] <- V[2L, 1L] <- .2
      y <- c(0, 0, residual)
      object <- brma.mv(yi = y, V = V, measure = "GEN",
        prior_unit_information_sd = 1,
        prior_effect = BayesTools::prior("point", list(location = 0)),
        only_priors = TRUE, silent = TRUE)
      setup <- .log_lik_evaluated_setup(fit = object$fit, data = object$data,
        priors = object$priors, unit = "estimate", data_hash = NULL,
        mu_samples = matrix(0, 1L, 3L), tau_within_samples = matrix(0, 1L, 3L),
        tau_between_samples = NULL, posterior_samples = NULL)
      context <- list(data = object$data, predictor_cache = new.env(parent = emptyenv()))
      expected <- -log(2 * pi) - if (mixed) .5 * log(.96) else 0
      expected <- expected + stats::dnorm(residual, sd = sqrt(variance), log = TRUE)
      observed <- list(fast = .iwmde_log_lik_known_v_joint_sum_common_shift(context, setup),
        generic = .log_lik_known_v_joint_sum_from_setup(setup))
      for (route in names(observed)) {
        value <- observed[[route]]
        if (is.finite(expected)) {
          expect_length(value, 1L)
          if (length(value) != 1L) next
          expect_true(is.finite(value))
          expect_lte(abs(value / expected - 1), 1e-13)
        } else if (mixed && route == "fast") {
          expect_null(value)
        } else expect_identical(as.numeric(value), -Inf)
      }
    }
  }
})


test_that("joint native Gaussian routes retain representable half squares", {

  group <- function(n, map = rep.int(1L, n), scale = .3) {
    list(type = "group", model_matrix = matrix(1, n, 1L), group_map = map,
      coefficient_structure = "diagonal", coefficient_factor = matrix(scale))
  }
  cases <- list(
    dense = list(V = diag(2L), factors = list(), covariance = diag(2L),
      route = "root_dense_blocks"),
    diagonal = list(V = diag(2L), factors = list(group(2L)),
      covariance = diag(2L) + .09, route = "low_rank_blocks"),
    spectral = list(V = matrix(c(1, .2, .2, 1), 2L), factors = list(group(2L)),
      covariance = matrix(c(1, .2, .2, 1), 2L) + .09, route = "spectral_blocks"),
    block_base = list(V = diag(4L), factors = list(group(4L)),
      covariance = diag(4L) + .09, route = "block_base_blocks")
  )
  cases$block_base$V[1:2, 1:2] <- matrix(c(1, .2, .2, 1), 2L)
  cases$block_base$covariance <- cases$block_base$V + .09
  K <- 12L
  parent <- rep(1:2, each = 6L)
  child <- rep(1:4, each = 3L)
  cases$sparse <- list(V = diag(K), factors = list(group(K, parent), group(K, child, .2)),
    covariance = diag(K) + .09 * outer(parent, parent, "==") +
      .04 * outer(child, child, "=="), route = "sparse_factor_blocks")
  known <- group(4L, rep(1:2, each = 2L))
  known$type <- "known_group"
  known$group_covariance <- matrix(c(1, .2, .2, 1), 2L)
  cases$fixed <- list(V = diag(4L), factors = list(known),
    covariance = diag(4L) + .09 * known$group_covariance[known$group_map, known$group_map],
    route = "fixed_known_group_blocks")
  markov <- list(type = "group", model_matrix = diag(2L), group_map = c(1L, 1L),
    coefficient_structure = "markov", coefficient_factor = .3 * t(chol(matrix(c(1, .2, .2, 1), 2L))),
    coefficient_scale = c(.3, .3), markov_transition = .2,
    markov_innovation_variance = .96)
  cases$markov <- list(V = diag(2L), factors = list(markov),
    covariance = diag(2L) + .09 * matrix(c(1, .2, .2, 1), 2L), route = "markov_blocks")

  for (name in names(cases)) {
    case <- cases[[name]]
    K <- nrow(case$V)
    for (magnitude in c(1.5e154, 1e155, .5)) {
      y <- c(magnitude, rep(0, K - 1L))
      plan <- .Call("RoBMA_known_v_covariance_plan_create", as.double(y),
        case$V, case$factors, list(seq_len(K)), PACKAGE = "RoBMA")
      expect_identical(attr(plan, case$route), 1L)
      state <- lapply(case$factors, function(factor) {
        factor[c("coefficient_factor", "coefficient_scale", "markov_transition",
          "markov_innovation_variance")][names(factor[c("coefficient_factor", "coefficient_scale",
          "markov_transition", "markov_innovation_variance")]) %in% names(factor)]
      })
      observed <- tryCatch(.Call("RoBMA_known_v_covariance_plan_loglik", plan,
        numeric(K), state, numeric(K), PACKAGE = "RoBMA"), error = identity)
      expect_type(observed, "double")
      if (!is.numeric(observed)) next
      root <- chol(case$covariance)
      standardized <- forwardsolve(t(root), y)
      expected <- -.5 * K * log(2 * pi) - sum(log(diag(root))) -
        sum((.5 * standardized) * standardized)
      if (is.finite(expected)) {
        expect_true(is.finite(observed))
        expect_lte(abs(observed / expected - 1), 1e-12)
      } else expect_identical(observed, -Inf)
    }
  }
})


test_that("joint Woodbury residual forms preserve the original PSD-base witness", {

  v <- c(.2, .3)
  y <- c(.1, -.05)
  factor <- list(type = "group", model_matrix = matrix(1, 2L, 1L),
    group_map = c(1L, 1L), coefficient_structure = "diagonal")
  plan <- .Call("RoBMA_known_v_covariance_plan_create", y, tcrossprod(v),
    list(factor), list(1:2), PACKAGE = "RoBMA")
  observed <- .Call("RoBMA_known_v_covariance_plan_loglik", plan, c(0, 0),
    list(list(coefficient_factor = matrix(.1))), c(0, 0), PACKAGE = "RoBMA")
  z <- solve(cbind(v, rep(.1, 2L)), y)
  expected <- -log(2 * pi) - log(.1) - log(abs(v[1L] - v[2L])) - .5 * sum(z^2)
  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("joint shared Gaussian cluster branches retain finite half squares", {

  cases <- list(c(0, 0, 0), c(.2, .2, 0), c(2, 2, 0), c(1e200, 1e-200, 0))
  for (loading in cases) {
    for (magnitude in c(1.5e154, 1e155, .5)) {
      y <- c(0, 1, magnitude)
      observed <- .Call("RoBMA_norm_cluster_analytic_loglik", y, rep(1, 3L),
        matrix(0, 1L, 3L), matrix(0, 1L, 3L), matrix(loading, 1L),
        1:3, 3L, NULL, PACKAGE = "RoBMA")
      if (max(loading) > 1e100) {
        # The tiny second loading changes these log densities by less than
        # 1e-300; evaluate the independent limiting product on its SD scale.
        expected <- stats::dnorm(0, sd = 1e200, log = TRUE) +
          stats::dnorm(1, log = TRUE) + stats::dnorm(magnitude, log = TRUE)
      } else {
        root <- chol(diag(3L) + tcrossprod(loading))
        standardized <- forwardsolve(t(root), y)
        expected <- -1.5 * log(2 * pi) - sum(log(diag(root))) -
          sum((.5 * standardized) * standardized)
      }
      if (is.finite(expected)) {
        expect_true(is.finite(observed))
        expect_lte(abs(as.numeric(observed) / expected - 1), 1e-13)
      } else expect_identical(as.numeric(observed), -Inf)
    }
  }
})


test_that("joint affine grids preserve finite half squares and true overflow", {

  coefficients <- c(0, .2, 1)
  for (magnitude in c(1.5e154, 1e155, .5)) {
    observed <- .iwmde_random_affine_log_likelihood_chunk(
      base_covariance = diag(2L), update_covariance = diag(2L),
      reference_coefficient = 0, coefficients = coefficients,
      means = matrix(0, 1L, 2L), outcome = c(magnitude, 0), blocks = list(1:2),
      plan_cache = new.env(parent = emptyenv()), sampling_covariance = diag(2L))
    expect_true(is.matrix(observed))
    if (!is.matrix(observed)) next
    for (g in seq_along(coefficients)) {
      expected <- sum(stats::dnorm(c(magnitude, 0),
        sd = sqrt(1 + coefficients[[g]]), log = TRUE))
      if (is.finite(expected)) {
        expect_true(is.finite(observed[g, 1L]))
        expect_lte(abs(observed[g, 1L] / expected - 1), 1e-13)
      } else expect_identical(observed[g, 1L], -Inf)
    }
  }
})
