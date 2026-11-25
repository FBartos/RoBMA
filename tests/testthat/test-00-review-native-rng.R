context("Native selection random-number ordering")

test_that("native RNG bin choices retain remote lower-tail mass", {

  skip_if_not(.has_native_selnorm_kernel())
  interval_quantile <- function(lower, upper, u) {
    if (lower > 0) return(stats::qnorm(.test_logsumexp(c(
      log1p(-u) + stats::pnorm(lower, lower.tail = FALSE, log.p = TRUE),
      log(u) + stats::pnorm(upper, lower.tail = FALSE, log.p = TRUE)
    )), lower.tail = FALSE, log.p = TRUE))
    if (upper < 0) return(stats::qnorm(.test_logsumexp(c(
      log1p(-u) + stats::pnorm(lower, log.p = TRUE),
      log(u) + stats::pnorm(upper, log.p = TRUE)
    )), log.p = TRUE))
    stats::qnorm((1 - u) * stats::pnorm(lower) + u * stats::pnorm(upper))
  }
  prior <- BayesTools::prior_weightfunction("one-sided", c(.025, .05),
    BayesTools::wf_fixed(c(1, 0, .5)))
  for (sign in c(1L, -1L)) for (sd in c(.019, .02, .5)) {
    spec <- .selection_spec(list(outcome = list(bias = prior)), 0, 1,
      effect_direction = if (sign == 1L) "positive" else "negative", signed_data = FALSE)
    mean <- sign * 1.8
    lower <- (spec$z_lower - sign * mean) / sd
    upper <- (spec$z_upper - sign * mean) / sd
    log_mass <- log(c(1, 0, .5)) + vapply(1:3, function(bin) {
      .test_interval_log_prob(lower[bin], upper[bin], 0, 1)
    }, numeric(1))
    probability <- exp(log_mass - .test_logsumexp(log_mass))
    withr::local_seed(409L)
    uniforms <- matrix(stats::runif(80L), 2L)
    bins <- vapply(uniforms[1L, ], function(u) which(u <= cumsum(probability))[1L], integer(1))
    expected <- mean + sign * sd * vapply(seq_len(40L), function(j) {
      interval_quantile(lower[bins[j]], upper[bins[j]], uniforms[2L, j])
    }, numeric(1))
    set.seed(409L)
    actual <- .selnorm_kernel_rng_matrix(matrix(mean, 40L), matrix(sd, 40L), 1,
      matrix(rep(c(1, 0, .5), each = 40L), 40L), spec)
    actual_bins <- vapply(as.numeric(actual), function(y) which(sign * y >= spec$z_lower)[1L], integer(1))
    expect_identical(actual_bins, bins)
    expect_equal(as.numeric(actual), expected, tolerance = 1e-11)
  }
})

test_that("kernel RNG consumes the bin uniform before the interval uniform", {

  spec <- .test_step_spec(0, 1)
  withr::local_seed(409L)
  uniforms <- runif(12L)
  expected <- matrix(qnorm(uniforms[seq.int(2L, 12L, by = 2L)]), 6L)
  set.seed(409L)
  actual <- .selnorm_kernel_rng_matrix(matrix(0, 6L), matrix(1, 6L), 1,
    matrix(1, 6L, spec$n_bins), spec, kernel_mode = rep(0L, 6L))
  expect_equal(actual, expected, tolerance = 1e-14)
})
