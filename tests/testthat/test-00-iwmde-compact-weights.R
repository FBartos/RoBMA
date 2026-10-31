.compact_fixture <- function(lower = 0, upper = Inf, conditioner = TRUE) {
  i <- seq_len(80L)
  u <- seq(-2, 2, length.out = 80L)
  z <- exp(.6 * u + .13 * sin(i))
  values <- if (is.finite(lower)) lower + z else upper - z
  samples <- cbind(T = values)
  priors <- list(T = BayesTools::prior("normal", list(mean = 0, sd = 10),
    truncation = list(lower = lower, upper = upper)))
  if (conditioner) {
    samples <- cbind(samples, U = u)
    priors$U <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  }
  list(context = .iwmde_context_ensure_caches(list(posterior_samples = samples,
    flat_prior_list = priors, indicator_names = character())), values = values,
    z = z, rows = i, support = c(lower, upper))
}
.compact_weight <- function(fixture, rows = fixture$rows, values = fixture$values[rows]) {
  .iwmde_chen_log_weight(fixture$context, "T", list(type = "primitive"),
    rows, values, fixture$rows, fixture$values)
}
.compact_reference <- function(fixture, rows, values) {
  conditioning <- .iwmde_chen_conditioning_matrix(fixture$context, "T", list(type = "primitive"), rows, fixture$rows)
  if (ncol(conditioning$fit)) {
    gaussian <- .iwmde_chen_conditional_gaussian(fixture$z, conditioning$fit, conditioning$eval)
    m <- gaussian$means; s <- gaussian$sd; scale <- gaussian$focal_scale
  } else {
    scale <- max(abs(fixture$z)); m <- rep(mean(fixture$z / scale), length(rows)); s <- sd(fixture$z / scale)
  }
  log_s0 <- pnorm(0, m, s, lower.tail = FALSE, log.p = TRUE)
  a <- qnorm(log_s0 + log(.975), m, s, lower.tail = FALSE, log.p = TRUE)
  b <- qnorm(log_s0 + log(.025), m, s, lower.tail = FALSE, log.p = TRUE)
  z <- if (is.finite(fixture$support[1L])) values - fixture$support[1L] else fixture$support[2L] - values
  r_value <- z / scale
  mass <- vapply(seq_along(rows), function(i) integrate(function(x) dnorm(x, m[i], s),
    a[i], b[i], rel.tol = 1e-11, abs.tol = 0)$value, numeric(1))
  result <- dnorm(r_value, m, s, log = TRUE) - log(mass) - log(scale)
  zero <- r_value < a | r_value > b
  result[zero] <- -Inf
  list(log_weight = result, zero = zero, a = a, b = b, m = m, s = s, scale = scale, mass = mass)
}

testthat::test_that("half-bounded weights use the normalized compact natural-distance conditional", {
  for (bounds in list(c(2, Inf), c(-Inf, 5))) {
    fixture <- .compact_fixture(bounds[1L], bounds[2L])
    rows <- c(80L, 1L, 40L, 7L, 61L)
    values <- fixture$values[rows]
    # An ordinary declared support point far beyond the compact interval has a
    # valid zero proposal; it remains in the same evaluation population.
    values[1L] <- if (is.finite(bounds[1L])) bounds[1L] + 100 else bounds[2L] - 100
    reference <- .compact_reference(fixture, rows, values)
    weight <- .compact_weight(fixture, rows, values)
    testthat::expect_true(any(reference$zero))
    testthat::expect_equal(weight$log_weight, reference$log_weight, tolerance = 2e-12)
    testthat::expect_type(weight$zero_support, "logical")
    testthat::expect_identical(unname(weight$zero_support), reference$zero)
    testthat::expect_identical(names(weight$zero_support), as.character(rows))
    testthat::expect_true(all(reference$a > 0 & reference$a < reference$b))
    testthat::expect_true(all(reference$mass > 0))
    testthat::expect_identical(weight$fallback_count, 0L)
  }
})

testthat::test_that("no-conditioner half-bounded weights use the same compact normal family", {
  fixture <- .compact_fixture(2, Inf, conditioner = FALSE)
  rows <- c(80L, 30L, 1L)
  values <- c(100, fixture$values[30L], 2)
  reference <- .compact_reference(fixture, rows, values)
  weight <- .compact_weight(fixture, rows, values)
  testthat::expect_equal(weight$log_weight, reference$log_weight, tolerance = 2e-12)
  testthat::expect_identical(unname(weight$zero_support), reference$zero)
  testthat::expect_match(weight$method, "compact")
})

testthat::test_that("typed collinearity falls back to the same compact marginal family", {
  fixture <- .compact_fixture(2, Inf)
  fixture$context$posterior_samples <- cbind(fixture$context$posterior_samples,
    V = fixture$context$posterior_samples[, "U"])
  fixture$context$flat_prior_list$V <- BayesTools::prior("normal", list(0, 1))
  fixture$context <- .iwmde_context_ensure_caches(fixture$context, reset = TRUE)
  rows <- c(1L, 30L, 80L)
  values <- fixture$values[rows]
  weight <- .compact_weight(fixture, rows, values)
  marginal <- .compact_fixture(2, Inf, conditioner = FALSE)
  reference <- .compact_reference(marginal, rows, values)
  testthat::expect_equal(weight$log_weight, reference$log_weight, tolerance = 2e-12)
  testthat::expect_identical(weight$fallback_count, 1L)
  testthat::expect_match(weight$fallback_reason, "not positive definite", fixed = TRUE)
  testthat::expect_match(weight$method, "compact.*marginal")
})

testthat::test_that("row-varying translated and reflected half supports retain their physical scale once", {
  i <- seq_len(80L)
  raw <- cbind(a = .1 + (i %% 17L) / 10, b = .2 + (i %% 13L) / 8)
  for (sign in c(1, -1)) {
    samples <- sign * raw
    bound <- if (sign > 0) list(lower = 0, upper = Inf) else list(lower = -Inf, upper = 0)
    prior <- BayesTools::prior("normal", list(0, 1), truncation = bound)
    context <- .iwmde_context_ensure_caches(list(posterior_samples = samples,
      flat_prior_list = list(a = prior, b = prior), indicator_names = character()))
    previous <- NULL
    for (scale in c(1, 2)) {
      spec <- list(type = "linear", weights = c(a = scale, b = scale))
      values <- rowSums(samples) * scale
      supports <- .iwmde_parameter_row_supports(context, "sum", i, spec)
      finite_bound <- sign * abs(raw[, "a"] - raw[, "b"]) * scale
      testthat::expect_equal(supports[, if (sign > 0) 1L else 2L], finite_bound, tolerance = 1e-13)
      weight <- .iwmde_chen_log_weight(context, "sum", spec, i, values, i, values)
      z <- 2 * pmin(raw[, "a"], raw[, "b"]) * scale
      reference_fixture <- .compact_fixture(0, Inf, conditioner = FALSE)
      reference_fixture$z <- z
      reference <- .compact_reference(reference_fixture, i, z)
      testthat::expect_equal(weight$log_weight, reference$log_weight, tolerance = 2e-12)
      testthat::expect_identical(unname(weight$zero_support), reference$zero)
      testthat::expect_identical(names(weight$zero_support), as.character(i))
      testthat::expect_length(weight$partitions, 1L)
      if (!is.null(previous)) testthat::expect_equal(weight$log_weight, previous - log(2), tolerance = 2e-12)
      previous <- weight$log_weight
    }
  }
})

test_that("compact normal construction normalizes actual endpoints for signed means and scales", {

  for (sd in c(.1, 1, 10)) for (ratio in c(-10, 0, 10)) for (scale in c(1e-100, 1, 1e100)) {
    mean <- ratio * sd
    log_s0 <- stats::pnorm(0, mean, sd, lower.tail = FALSE, log.p = TRUE)
    lower <- stats::qnorm(log_s0 + log(.975), mean, sd, lower.tail = FALSE, log.p = TRUE)
    upper <- stats::qnorm(log_s0 + log(.025), mean, sd, lower.tail = FALSE, log.p = TRUE)
    integral <- stats::integrate(function(r) {
      weight <- .iwmde_chen_compact_normal_log_weight(r * scale, mean, sd, scale, "control")
      exp(weight[["log_weight"]] + log(scale))
    }, lower, upper, rel.tol = 1e-11, abs.tol = 0)
    expect_lte(abs(integral[["value"]] - 1), 1e-10 + integral[["abs.error"]])
  }
})

test_that("invalid compact construction never certifies a zero density", {

  inputs <- list(c(NA, 1, 1), c(Inf, 1, 1), c(0, 0, 1), c(0, -1, 1),
    c(0, Inf, 1), c(0, 1, 0), c(0, 1, Inf), c(1e300, 1e-300, 1))
  for (input in inputs) {
    weight <- .iwmde_chen_compact_normal_log_weight(c(0, 1, 2), input[1L], input[2L], input[3L], "control")
    expect_true(all(weight[["log_weight"]] == -Inf))
    expect_false(any(weight[["zero_support"]]))
  }
  weight <- .iwmde_chen_compact_marginal_normal_log_weight(c(0, 1), rep(1, 4L), c(0, Inf))
  expect_true(all(weight[["log_weight"]] == -Inf))
  expect_false(any(weight[["zero_support"]]))
})
