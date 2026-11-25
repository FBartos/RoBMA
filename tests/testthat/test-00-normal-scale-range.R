.normal_scale_range_symbols <- function() {

  for (symbol in c("RoBMA_norm_loglik_row_sum", "RoBMA_norm_predictor_grid_loglik",
    "RoBMA_norm_cluster_analytic_loglik",
    "RoBMA_norm_cluster_analytic_loglik_row_sum")) {
    expect_s3_class(getNativeSymbolInfo(symbol, PACKAGE = "RoBMA"), "NativeSymbolInfo")
  }
}

test_that("normal row sums retain finite extreme scales and likelihood powers", {

  .normal_scale_range_symbols()
  yi <- c(.1, -.2)
  mu <- matrix(c(0, .2, .1, -.1), 2L)
  sei <- c(1, 1)
  for (tau in list(matrix(1e200, 2L, 2L), matrix(c(0, .2, .3, .5), 2L))) {
    # The extreme reference uses the independently representable SD directly.
    sd <- if (all(tau == 1e200)) matrix(1e200, 2L, 2L) else sqrt(tau^2 + 1)
    reference <- stats::dnorm(matrix(yi, 2L, 2L, byrow = TRUE),
      mean = mu, sd = sd, log = TRUE)
    for (weights in list(NULL, c(.5, 1.5))) {
      expected <- rowSums(if (is.null(weights)) reference else
        sweep(reference, 2L, weights, "*"))
      actual <- .Call("RoBMA_norm_loglik_row_sum", yi, mu, tau, sei, weights,
        PACKAGE = "RoBMA")
      expect_true(all(is.finite(actual)))
      expect_equal(actual, expected, tolerance = 1e-14)
      expect_equal(.outcome_pdf_sum.norm(yi, mu, tau, sei, weights), expected,
        tolerance = 1e-14)
      pointwise <- .outcome_pdf.norm(yi, mu, tau, sei)
      expect_equal(pointwise, reference, tolerance = 1e-14)
    }
  }
})

test_that("normal predictor grids keep admitted extreme candidates finite", {

  .normal_scale_range_symbols()
  yi <- c(.1, -.2)
  sei <- c(1, 1)
  mu <- matrix(0, 1L, 2L)
  for (replace_tau in c(FALSE, TRUE)) {
    values <- if (replace_tau) c(1e200, 2e200) else c(0, .2)
    actual <- .Call("RoBMA_norm_predictor_grid_loglik", yi, sei, c(.5, 1.5),
      mu, if (replace_tau) NULL else matrix(1, 1L, 2L),
      if (replace_tau) NULL else matrix(1e200, 1L, 2L), NULL,
      0, values, replace_tau, FALSE, PACKAGE = "RoBMA")
    expect_true(all(actual$valid))
    reference <- vapply(values, function(value) {
      sd <- if (replace_tau) value else 1e200
      location <- if (replace_tau) 0 else value
      sum(c(.5, 1.5) * stats::dnorm(yi, mean = location, sd = sd, log = TRUE))
    }, numeric(1L))
    expect_true(all(is.finite(actual$log_lik)))
    expect_equal(as.numeric(actual$log_lik), reference, tolerance = 1e-14)
    sums <- vapply(values, function(value) {
      .Call("RoBMA_norm_loglik_row_sum", yi,
        matrix(if (replace_tau) 0 else value, 1L, 2L),
        matrix(if (replace_tau) value else 1e200, 1L, 2L),
        sei, c(.5, 1.5), PACKAGE = "RoBMA")
    }, numeric(1L))
    expect_equal(as.numeric(actual$log_lik), sums, tolerance = 1e-14)
  }
  invalid <- .Call("RoBMA_norm_predictor_grid_loglik", yi, sei, NULL, mu, NULL,
    NULL, NULL, 0, c(-1, Inf, NaN), TRUE, FALSE, PACKAGE = "RoBMA")
  expect_identical(as.numeric(invalid$valid), c(0, 0, 0))
  expect_error(.Call("RoBMA_norm_predictor_grid_loglik", yi, sei, NULL, mu,
    NULL, matrix(1, 1L, 2L), NULL, 0, 1, TRUE, FALSE, PACKAGE = "RoBMA"),
    "'tau' must be absent when the grid replaces the scale itself.", fixed = TRUE)
})

test_that("normal cluster likelihood retains finite extreme scales and powers", {

  .normal_scale_range_symbols()
  check <- function(setup, yi, sei, expected) {
    cluster <- .cluster_indices_flatten(setup$cluster)
    native_args <- list(yi, sei, setup$mu, setup$tau_within, setup$tau_between,
      cluster$index, cluster$size, setup$weights)
    matrix_value <- do.call(.Call, c(
      "RoBMA_norm_cluster_analytic_loglik", native_args, PACKAGE = "RoBMA"))
    row_sum <- do.call(.Call, c(
      "RoBMA_norm_cluster_analytic_loglik_row_sum", native_args,
      PACKAGE = "RoBMA"))
    expect_true(all(is.finite(matrix_value)))
    expect_true(all(is.finite(row_sum)))
    expect_equal(matrix_value, expected, tolerance = 1e-14)
    expect_equal(row_sum, rowSums(expected), tolerance = 1e-14)
    expect_equal(.log_lik_cluster_norm_analytic(setup, yi, sei), expected,
      tolerance = 1e-14)
    expect_equal(.log_lik_cluster_norm_analytic_sum(setup, yi, sei),
      rowSums(expected), tolerance = 1e-14)
    expect_equal(.log_lik_cluster_norm_analytic_r(setup, yi, sei), expected,
      tolerance = 1e-14)
  }
  yi <- c(.1, -.2, .3)
  sei <- rep(1, 3L)
  setup <- list(S = 2L, mu = matrix(0, 2L, 3L),
    tau_within = matrix(1e200, 2L, 3L), tau_between = matrix(0, 2L, 3L),
    cluster = list(1:2, 3L), weights = rep(1, 3L))
  # Without a latent loading the cluster density is the product of its rows.
  expected <- matrix(rep(c(sum(stats::dnorm(yi[1:2], sd = 1e200, log = TRUE)),
    stats::dnorm(yi[3L], sd = 1e200, log = TRUE)), each = setup$S), setup$S)
  check(setup, yi, sei, expected)

  setup$tau_within[,] <- .2
  setup$weights <- c(.5, 1, 1.5)
  reference <- setup$weights * stats::dnorm(yi, sd = sqrt(1 + .2^2), log = TRUE)
  expected <- matrix(rep(c(sum(reference[1:2]), reference[3L]), each = setup$S), setup$S)
  check(setup, yi, sei, expected)

  # A latent effect of 1e200 against a sampling SE of 1e-200: the information
  # 1e800 is carried on the log scale; the closed form is
  # -log(2 pi) / 2 - log(sd) - log(1 + (tau / sd)^2) / 2 at a zero residual.
  singleton <- list(S = 1L, mu = matrix(0, 1L, 1L),
    tau_within = matrix(0, 1L, 1L), tau_between = matrix(1e200, 1L, 1L),
    cluster = list(1L), weights = 1)
  expected <- -log(2 * pi) / 2 + 200 * log(10) - (800 * log(10)) / 2
  expect_equal(expected, -461.43595713201381, tolerance = 1e-15)
  check(singleton, 0, 1e-200, matrix(expected, 1L, 1L))
})


test_that("internal normal kernels retain continuous tiny positive total SD", {

  .normal_scale_range_symbols()
  # Internal-only guard: public SE construction may refuse these scales.
  sd <- sqrt(2) * 1e-200
  expected <- stats::dnorm(0, sd = sd, log = TRUE)
  actual <- .Call("RoBMA_norm_loglik_row_sum", 0, matrix(0, 1L),
    matrix(1e-200, 1L), 1e-200, NULL, PACKAGE = "RoBMA")
  expect_true(is.finite(actual))
  expect_equal(actual, expected, tolerance = 1e-14)
  grid <- .Call("RoBMA_norm_predictor_grid_loglik", 0, 1e-200, NULL,
    matrix(0, 1L), NULL, matrix(1e-200, 1L), NULL,
    0, 0, FALSE, FALSE, PACKAGE = "RoBMA")
  expect_true(all(grid$valid))
  expect_true(all(is.finite(grid$log_lik)))
  expect_equal(as.numeric(grid$log_lik), expected, tolerance = 1e-14)
})
