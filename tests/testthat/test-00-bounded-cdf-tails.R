.bounded_cdf_tail_call <- function(means, omega = c(1, .3), bounded = TRUE,
                                   sign = 1L, denominator = NULL, collect = FALSE,
                                   z = c(-9, -8, 0, 8, 9), cache = NULL, budget = 0) {

  means <- as.matrix(means)
  k <- ncol(means)
  covariance <- diag(1, k) + .25
  nodes <- rep(0, 7L)
  attr(nodes, "bounded_cdf") <- bounded
  if (!is.null(cache)) attr(nodes, "context_rule_cache") <- cache
  if (budget > 0) attr(nodes, "context_rule_cache_budget") <- as.double(budget)
  .Call("RoBMA_selnorm_context_star_zplot", means,
    as.double(covariance[lower.tri(covariance, diag = TRUE)]), rep(1, k),
    omega, c(0, -Inf), c(Inf, 0), sign, TRUE, denominator,
    rep(if (collect) -log(nrow(means)) else -Inf, nrow(means)),
    nodes, rep(-log(length(nodes)), length(nodes)), as.double(z), 0,
    PACKAGE = "RoBMA")
}


test_that("bounded CDF tails saturate only beyond their certified threshold", {

  # At a repeated zero quadrature node every residual row is N(mean, 1).
  # The adjacent doubles below eight retain the ordinary tail evaluation.
  inside <- 8 - .Machine$double.eps * 4
  outside <- 8 + .Machine$double.eps * 8
  for (sign in c(-1L, 1L)) {
    for (omega in list(c(1, .3), c(.3, 1))) {
      for (side in c(-1, 1)) {
        means <- sign * side * c(inside, 8, outside, 40)
        values <- .bounded_cdf_tail_call(cbind(means, means), omega, sign = sign)
        raw <- as.double(values$compact_log_normalizers)
        expect_identical(raw[[2L]], raw[[3L]])
        expect_identical(raw[[2L]], raw[[4L]])
        expect_false(identical(raw[[1L]], raw[[2L]]))
        expect_true(all(attr(values$compact_log_normalizers, "cdf_log_error") > 0))
        ordinary <- .bounded_cdf_tail_call(cbind(means, means), omega, FALSE, sign)
        expect_identical(attr(ordinary$compact_log_normalizers, "cdf_log_error"), rep(0, 4L))
      }
    }
  }
})


test_that("the exact normal-minimum floor retains its ordinary tail reference", {

  # A large relative weight makes the nonzero tail at eight visible above
  # DBL_MIN. Replacing that tail by zero could force a new log fallback, whose
  # addition order differs from the mass-only direct-product evaluator.
  means <- rbind(c(-8, -8), c(-9, 9), c(-40, 40))
  for (sign in c(-1L, 1L)) {
    for (relative in list(c(1e8, 1), c(1, 1e8))) {
      omega <- .Machine$double.xmin * relative
      bounded <- .bounded_cdf_tail_call(sign * means, omega, sign = sign)
      ordinary <- .bounded_cdf_tail_call(sign * means, omega, FALSE, sign)
      expect_identical(as.double(bounded$compact_log_normalizers),
        as.double(ordinary$compact_log_normalizers))
      expect_true(all(attr(bounded$compact_log_normalizers, "cdf_log_error") > 0))
      expect_identical(attr(ordinary$compact_log_normalizers, "cdf_log_error"), rep(0, 3L))
    }
  }
})


test_that("raw context caches from the older CDF primitive are not reused", {

  means <- rbind(c(-8, 8), c(-9, 9))
  initial <- .bounded_cdf_tail_call(means, collect = TRUE, budget = 1024^2)
  cache <- attr(initial, "context_rule_cache")
  expect_type(cache, "list")
  expect_identical(cache[[1L]][[1L]], 2)
  # A stale key must miss before even interpreting its former payload.
  cache[[1L]][[1L]] <- 1
  cache[[2L]][] <- NaN
  actual <- .bounded_cdf_tail_call(means, collect = TRUE, cache = cache)
  expected <- .bounded_cdf_tail_call(means, collect = TRUE)
  expect_null(attr(actual, "context_rule_cache"))
  expect_null(attr(actual, "context_rule_cache_hit"))
  expect_identical(actual, expected)
})


test_that("bounded tail products preserve the independent finite-rule law and ledger", {

  # Absolute weights may make long-double node products underflow or overflow.
  # Those cases must still use the existing complete log-scale fallback.
  configurations <- list(
    list(scale = .Machine$double.xmin, k = 2L),
    list(scale = 1e-300, k = 20L),
    list(scale = 1e300, k = 20L),
    list(scale = 1, k = 2L))
  for (configuration in configurations) {
    k <- configuration$k
    for (sign in c(-1L, 1L)) {
      omega <- configuration$scale * c(2, 1)
      means <- matrix(sign * rep(c(-8, 8, -9, 9), length.out = 2L * k), 2L, k)
      normalizer <- log(omega[[2L]]) + log1p(stats::pnorm(sign * means))
      expected_raw <- rowSums(normalizer)
      mass_only <- .bounded_cdf_tail_call(means, omega, sign = sign)
      raw <- as.double(mass_only$compact_log_normalizers)
      cdf_mass <- attr(mass_only$compact_log_normalizers, "cdf_log_error")
      rounding <- 16 * .Machine$double.eps * pmax(1, abs(expected_raw))
      expect_true(all(is.finite(raw)))
      expect_true(all(abs(raw - expected_raw) <= cdf_mass + rounding))
      expect_identical(mass_only$available, FALSE)
      expect_identical(mass_only$mass, 0)

      # An explicit shared denominator also exercises log-fallback projection.
      actual <- .bounded_cdf_tail_call(means, omega, sign = sign,
        denominator = raw, collect = TRUE)
      expected <- vapply(c(-9, -8, 0, 8, 9), function(point) {
        weight <- omega[[if (sign * point >= 0) 1L else 2L]]
        mean(exp(stats::dnorm(point, means, log = TRUE) + log(weight) - normalizer))
      }, numeric(1L))
      cdf_numerator <- attr(actual$compact_log_normalizers, "cdf_numerator_log_error")
      relative_bound <- expm1(max(cdf_mass + cdf_numerator) + 2 * max(rounding))
      expect_true(actual$available)
      expect_true(all(is.finite(actual$density)))
      expect_true(all(abs(actual$density - expected) <=
        max(expected) * relative_bound + 16 * .Machine$double.eps))
      expect_lte(abs(actual$mass - 1), expm1(2 * max(cdf_mass)) + max(rounding))
      # The existing mass-only and collected log fallbacks have different
      # addition orders. Bound each result against the independent law; the
      # controller's separate shared-D refusal test retains its exact check.
      expect_true(all(abs(as.double(actual$compact_log_normalizers) - expected_raw) <=
        cdf_mass + rounding))
    }
  }
})


test_that("bounded tail routing preserves unsupported weights and invalid inputs", {

  means <- rbind(c(-9, 9), c(-8, 8))
  for (omega in list(c(1, 0), c(1, 1e-12), c(1, 1))) {
    bounded <- .bounded_cdf_tail_call(means, omega)
    ordinary <- .bounded_cdf_tail_call(means, omega, FALSE)
    expect_identical(bounded, ordinary)
  }
  for (invalid in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- means
    bad[[1L]] <- invalid
    for (bounded in c(FALSE, TRUE)) {
      expect_error(.bounded_cdf_tail_call(bad, bounded = bounded), "finite")
    }
  }
})
