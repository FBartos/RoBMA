probe_state("test-00-selection-predictive-summary.R")
# ============================================================================ #
# test-00-selection-predictive-summary.R
# ============================================================================ #

context("Selected-normal predictive cell summaries")
skip_on_cran()

source(testthat::test_path("helper-selection-kernel.R"))


# Independent base-R reference for one cell: the bin masses in log space from
# R's normal tails, the paired log tails at q, and the truncated-normal moments.
.predictive_cell_reference <- function(q, mean, sd, sei, omega, static) {

  sign <- static[["sign"]]
  m    <- sign * mean
  qs   <- (sign * q - m) / sd
  log_interval <- function(lo, hi) {
    if (!(lo < hi)) return(-Inf)
    if (lo >= 0) {
      a <- stats::pnorm(lo, lower.tail = FALSE, log.p = TRUE)
      b <- stats::pnorm(hi, lower.tail = FALSE, log.p = TRUE)
      return(a + log1p(-exp(b - a)))
    }
    if (hi <= 0) {
      a <- stats::pnorm(hi, log.p = TRUE)
      b <- stats::pnorm(lo, log.p = TRUE)
      return(a + log1p(-exp(b - a)))
    }
    log(stats::pnorm(hi) - stats::pnorm(lo))
  }
  log_sum <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(-Inf)
    max(x) + log(sum(exp(x - max(x))))
  }
  x_phi <- function(x) if (is.finite(x)) x * stats::dnorm(x) else 0

  lower_terms <- upper_terms <- numeric()
  m0 <- c1 <- c2 <- 0
  for (b in seq_along(omega)) {
    if (omega[b] == 0) next
    lo <- (static[["z_lower"]][b] * sei - m) / sd
    hi <- (static[["z_upper"]][b] * sei - m) / sd
    log_w <- log(omega[b])
    lower_terms <- c(lower_terms, log_w + log_interval(lo, min(hi, qs)))
    upper_terms <- c(upper_terms, log_w + log_interval(max(lo, qs), hi))
    p  <- exp(log_interval(lo, hi))
    m0 <- m0 + omega[b] * p
    c1 <- c1 + omega[b] * (stats::dnorm(lo) - stats::dnorm(hi))
    c2 <- c2 + omega[b] * (p + x_phi(lo) - x_phi(hi))
  }
  lower <- log_sum(lower_terms)
  upper <- log_sum(upper_terms)
  total <- log_sum(c(lower, upper))
  tails <- c(lower - total, upper - total)
  if (sign == -1) tails <- rev(tails)
  list(
    log_lower = tails[1L],
    log_upper = tails[2L],
    mean      = sign * (m + sd * c1 / m0),
    variance  = sd^2 * (c2 / m0 - (c1 / m0)^2),
    mass      = m0
  )
}


test_that("selected moments retain bins with underflowing Gaussian probability", {

  skip_if_not(.has_native_selnorm_kernel())

  cuts <- c(stats::pnorm(1, lower.tail = FALSE), .5)
  prior <- prior_weightfunction("one-sided", cuts, wf_fixed(c(1, .5, .25)))
  # In y/se units N(-40, 1) has density phi(40) * exp(-40*y - y^2/2).
  # Scale by the middle bin's weighted boundary density before integration.
  # Below zero, its unconditional moments equal the full Gaussian moments
  # to double precision (the omitted upper tail is smaller than 1e-348).
  tail_integrals <- function(lo, hi) vapply(0:2, function(k) {
    stats::integrate(function(y) y^k * exp(-40*y - y^2/2), lo, hi,
                     abs.tol = 0, rel.tol = 1e-12)$value
  }, numeric(1))
  middle <- tail_integrals(0, 1)
  upper  <- tail_integrals(1, Inf)
  lower  <- c(1, -40, 1601)
  balanced <- exp(log(1e300) + stats::pnorm(-40, log.p = TRUE))
  for (weights in list(c(1, 1, 0), c(1, 1e300, 1e-100), c(1, 1e300, balanced))) {
    moments <- middle + upper * (weights[1] / weights[2]) +
      lower * exp(log(weights[3]) - log(weights[2]) - stats::dnorm(40, log = TRUE))
    ref <- moments / moments[1]
    for (direction in c("positive", "negative")) {
      sign <- if (direction == "positive") 1 else -1
      se <- .01
      spec <- .selection_spec(list(outcome = list(bias = prior)), 0, se,
                              effect_direction = direction, signed_data = FALSE)
      args <- list(mean = matrix(sign * -40 * se), sd = matrix(se), sei = se,
                   omega = matrix(weights, 1L), selection_spec = spec)
      raw <- do.call(.selnorm_kernel_moments_matrix, args)
      predictive <- do.call(.selnorm_kernel_predictive_summary_matrix,
                            c(list(q = sign * se / 40), args))
      expect_lt(abs(raw$second / (ref[3] * se^2) - 1), 1e-8)
      for (value in list(raw, predictive)) {
        expect_lt(abs(value$mean / (sign * se * ref[2]) - 1), 1e-8)
        expect_lt(abs(value$variance / (se^2 * (ref[3] - ref[2]^2)) - 1), 1e-8)
      }
    }
  }
})

test_that("hard selection retains excess and variance across the tail threshold", {

  skip_if_not(.has_native_selnorm_kernel())

  prior <- prior_weightfunction("one-sided", .5, wf_fixed(c(1, 0)))
  spec <- .selection_spec(list(outcome = list(bias = prior)), 0, 1,
                          effect_direction = "positive", signed_data = FALSE)
  for (a in c(35, 36, 37.4, 37.6, 38.5, 40, 54, 1000)) {
    # On s = a*y the density is proportional to exp(-s - s^2/(2*a^2)).
    ints <- vapply(0:2, function(k) stats::integrate(function(s) {
      s^k * exp(-s - .5*(s/a)^2)
    }, 0, Inf, rel.tol = 1e-12)$value, numeric(1))
    ref_mean <- ints[2L] / ints[1L] / a
    ref_var <- (ints[3L] / ints[1L] - (ints[2L] / ints[1L])^2) / a^2
    cell <- .selnorm_kernel_moments_matrix(matrix(-a), matrix(1), 1,
                                            matrix(c(1, 0), 1), spec)
    expect_lt(abs(cell$mean/ref_mean - 1), 1e-8)
    expect_lt(abs(cell$variance/ref_var - 1), 1e-6)
  }

  # Two equally remote supported tails both contribute to the variance.
  spec <- .test_two_sided_step_spec(0, 1)
  static <- BayesTools::selection_native_static_args(spec)
  cutoff <- static$z_lower[1L]
  sd <- cutoff / 40
  w <- numeric(spec$n_bins)
  w[c(1L, length(w))] <- c(1, .5)
  cell <- .selnorm_kernel_moments_matrix(matrix(0), matrix(sd), 1, matrix(w, 1), spec)
  excess <- stats::integrate(function(s) s * exp(-s - .5*(s/40)^2), 0, Inf)$value /
    stats::integrate(function(s) exp(-s - .5*(s/40)^2), 0, Inf)$value / 40
  ref_mean <- sd * (40 + excess) / 3
  ref_var <- sd^2 * (1 + 40 * (40 + excess)) - ref_mean^2
  expect_lt(abs(cell$mean/ref_mean - 1), 1e-8)
  expect_lt(abs(cell$variance/ref_var - 1), 1e-8)

  # Use the actual stored cuts: their magnitudes can differ by one ulp.
  # The peak log ratio factors exactly in outcome units. At 1e7 SDs the
  # conditional excess is below 3e-14 in these units and the relative
  # Mills-factor difference of the two tails is below 1e-15.
  left <- utils::tail(static$z_upper, 1L)
  sd <- cutoff / 1e7
  mu <- 1e-18
  w[length(w)] <- 1
  cell <- .selnorm_kernel_moments_matrix(matrix(mu), matrix(sd), 1, matrix(w, 1), spec)
  midpoint <- (cutoff + left) / 2
  log_ratio <- (cutoff - left) * ((mu - midpoint) / sd^2)
  ref_mean <- midpoint + (cutoff - left) / 2 * tanh(log_ratio / 2)
  expect_lt(abs(cell$mean - ref_mean), 1e-12)
})

test_that("predictive cell summaries match an independent log-space reference", {

  skip_if_not(.has_native_selnorm_kernel())

  sei   <- .2
  grid  <- expand.grid(mean = c(-3, -.4, 0, .25, .6, 3),
                       sd = c(.01, .08, .3, 1.5, 6))
  S     <- nrow(grid)
  # Observations from the center to scores far past the underflow range.
  q     <- c(-12, -.5, -.05, 0, .1, .33, .4, 2, 15)
  K     <- length(q)

  for (two_sided in c(FALSE, TRUE)) for (direction in c("positive", "negative")) {
    spec <- if (two_sided) {
      .test_two_sided_step_spec(q, rep(sei, K), effect_direction = direction)
    } else {
      .test_step_spec(q, rep(sei, K), effect_direction = direction)
    }
    static <- BayesTools::selection_native_static_args(spec)
    B      <- spec[["n_bins"]]
    # Monotone, nonmonotone, zero and single-supported weight rows.
    shapes <- list(
      seq(1, .2, length.out = B),
      rev(seq(.1, 1, length.out = B)),
      replace(rep(1, B), seq(2L, B, by = 2L), 0),
      replace(rep(0, B), B, 1),
      c(1, rep(.5, B - 1L))
    )
    omega <- do.call(rbind, rep(shapes, length.out = S))
    mean  <- matrix(grid$mean, S, K)
    sd    <- matrix(grid$sd, S, K)
    kernel_mode <- rep(SELKERNEL_STEP, S)
    kernel_mode[seq(5L, S, by = 7L)] <- SELKERNEL_NORMAL

    got <- .selnorm_kernel_predictive_summary_matrix(
      q = q, mean = mean, sd = sd, sei = rep(sei, K), omega = omega,
      selection_spec = spec, kernel_mode = kernel_mode
    )
    label <- paste(if (two_sided) "two-sided" else "one-sided", direction)

    for (s in seq_len(S)) for (k in seq_len(K)) {
      info <- paste(label, "row", s, "observation", k)
      if (kernel_mode[s] == SELKERNEL_NORMAL) {
        expected <- list(
          log_lower = stats::pnorm(q[k], mean[s, k], sd[s, k], log.p = TRUE),
          log_upper = stats::pnorm(q[k], mean[s, k], sd[s, k], lower.tail = FALSE, log.p = TRUE),
          mean      = mean[s, k],
          variance  = sd[s, k]^2,
          mass      = 1
        )
      } else {
        expected <- .predictive_cell_reference(
          q[k], mean[s, k], sd[s, k], sei, omega[s, ], static
        )
      }
      for (tail in c("log_lower", "log_upper")) {
        value <- got[[tail]][s, k]
        if (is.infinite(expected[[tail]])) {
          expect_identical(value, expected[[tail]], info = paste(info, tail))
        } else {
          expect_lt(abs(value - expected[[tail]]),
                    1e-11 * max(1, abs(expected[[tail]])), label = paste(info, tail))
        }
      }
      # Moments are compared where the reference mass is representable.
      if (expected[["mass"]] > 1e-250) {
        expect_lt(abs(got$mean[s, k] - expected$mean),
                  1e-10 * max(sd[s, k], abs(expected$mean)), label = paste(info, "mean"))
        expect_lt(abs(got$variance[s, k] - expected$variance),
                  1e-9 * sd[s, k]^2, label = paste(info, "variance"))
      }
    }

    # The existing kernels give the same tails and moments.
    cdf <- .selnorm_kernel_cdf_matrix(
      q = q, mean = mean, sd = sd, sei = rep(sei, K), omega = omega,
      selection_spec = spec, kernel_mode = kernel_mode
    )
    moments <- .selnorm_kernel_moments_matrix(
      mean = mean, sd = sd, sei = rep(sei, K), omega = omega,
      selection_spec = spec, kernel_mode = kernel_mode
    )
    expect_lt(max(abs(exp(got$log_lower) - cdf)), 1e-13, label = label)
    finite <- is.finite(moments$mean) & is.finite(moments$variance)
    expect_true(mean(finite) > .5, label = label)
    expect_lt(max(abs(got$mean - moments$mean)[finite] /
                  pmax(sd, abs(moments$mean))[finite]), 1e-11, label = label)
    expect_lt(max(abs(got$variance - moments$variance)[finite] / sd[finite]^2),
              1e-10, label = label)
  }
})


test_that("selected log tails resolve a positive query distance from a hard cutoff", {

  skip_if_not(.has_native_selnorm_kernel())
  prior <- prior_weightfunction("one-sided", .5, wf_fixed(c(1, 0)))
  for (direction in c("positive", "negative")) {
    sign <- if (direction == "positive") 1 else -1
    for (side in c(-1, 1)) {
      weights <- if (side == 1) c(1, 0) else c(0, 1)
      spec <- .selection_spec(list(outcome = list(bias = prior)), 0, 1,
                              effect_direction = direction, signed_data = FALSE)
      for (mu in c(0, .001, -3)) {
        q <- sign * side * c(0, 1e-308, 1e-17, 1e-12, 1e-6)
        K <- length(q)
        args <- list(mean = matrix(sign*mu, 1L, K), sd = matrix(1, 1L, K),
                     sei = rep(1, K), omega = matrix(weights, 1L), selection_spec = spec)
        got <- do.call(.selnorm_kernel_predictive_summary_matrix, c(list(q = q), args))
        lower <- sign*side > 0
        tail <- drop(if (lower) got$log_lower else got$log_upper)
        expect_identical(tail[1], -Inf)
        # Scale the integral onto [0,1], so a subnormal width appears only
        # in log units. The conditional normalizer is the retained half-line.
        ref <- vapply(abs(q[-1]), function(h) {
          log(h) + log(stats::integrate(function(u) stats::dnorm(side*h*u, mu), 0, 1)$value) -
            stats::pnorm(0, mu, lower.tail = side < 0, log.p = TRUE)
        }, numeric(1))
        expect_lt(max(abs(tail[-1] - ref)), 1e-10)
        cdf <- do.call(.selnorm_kernel_cdf_matrix, c(list(q = q, lower.tail = lower), args))
        expect_lt(max(abs(log(cdf[-1]) - ref)), 1e-10)
      }
    }
  }
})

test_that("a nonzero cut distinguishes its endpoint from the next outcome interval", {

  prior <- prior_weightfunction("one-sided", .025, wf_fixed(c(1, 0)))
  sei <- 2
  spec <- .selection_spec(list(outcome = list(bias = prior)), 0, sei,
                          effect_direction = "positive", signed_data = FALSE)
  cut <- BayesTools::selection_native_static_args(spec)$z_lower[1L] * sei
  q <- c(cut - 2^-40, cut, cut + 2^-40)
  mu <- cut + .01
  out <- .selnorm_kernel_predictive_summary_matrix(q, matrix(mu, 1, 3), matrix(1, 1, 3),
    rep(sei, 3), matrix(c(1, 0), 1), spec)
  ref <- log(q[3] - cut) + stats::dnorm((q[3]+cut)/2, mu, log = TRUE) -
    stats::pnorm(cut, mu, lower.tail = FALSE, log.p = TRUE)
  expect_identical(as.numeric(out$log_lower[1:2]), c(-Inf, -Inf))
  expect_lt(abs(out$log_lower[3] - ref), 1e-12)
})

test_that("a collapsed score interval retains its outcome-scale exponential tail", {

  spec <- .selection_spec(list(outcome = list(bias = prior_weightfunction(
    "one-sided", .5, wf_fixed(c(1, 0))))), 0, 1,
    effect_direction = "positive", signed_data = FALSE)
  h <- 5e-11
  a <- 1e6
  for (mu in c(a, -a)) {
    summary <- .selnorm_kernel_predictive_summary_matrix(h, matrix(mu), matrix(1), 1,
                                                         matrix(c(1,0), 1L), spec)
    # The omitted quadratic term is bounded by h^2/2. In the rare selected
    # tail the additional Mills-ratio correction is at most 1/a^2.
    ref <- if (mu > 0) {
      stats::dnorm(0, mu, log = TRUE) + log(expm1(a*h)/a)
    } else log(-expm1(-a*h))
    # The rare-tail normalization subtracts log masses of order 5e11.
    expect_lt(abs(summary$log_lower[1L] - ref), 3e-4)
  }
})

test_that("predictive cell tails stay finite for widely spaced finite weights", {

  skip_if_not(.has_native_selnorm_kernel())

  # Both tail masses are representable, but their ratio is about 1e350.
  prior <- BayesTools::prior_weightfunction(
    "one-sided", c(.025, .5), BayesTools::wf_fixed(c(1, 1e100, 1e-250))
  )
  # Weights this wide disable probability telescoping, which the cell kernel
  # does not use.
  expect_warning(
    spec <- .selection_spec(list(outcome = list(bias = prior)), yi = -1, sei = 1,
                            effect_direction = "positive", signed_data = FALSE),
    "telescoping disabled"
  )
  args <- list(mean = matrix(0), sd = matrix(1), sei = 1,
               omega = matrix(c(1, 1e100, 1e-250), 1L), selection_spec = spec,
               kernel_mode = SELKERNEL_STEP)
  log_sum  <- function(x) max(x) + log(sum(exp(x - max(x))))
  log_mass <- log_sum(c(log(.025), log(1e100) + log(.475), log(1e-250) + log(.5)))
  expected <- log(1e-250) + stats::pnorm(-1, log.p = TRUE) - log_mass

  summary  <- do.call(.selnorm_kernel_predictive_summary_matrix, c(list(q = -1), args))
  weighted <- do.call(.selnorm_kernel_weighted_summary,
                      c(list(yi = -1, psis_weights = matrix(1)), args))
  expect_equal(summary$log_lower[1L, 1L], expected, tolerance = 1e-12)
  expect_equal(weighted$log_lower, expected, tolerance = 1e-12)
  expect_equal(summary$log_upper[1L, 1L], log1p(-exp(expected)))
})


test_that("predictive cell tails keep heavily weighted bins whose probability underflows", {

  skip_if_not(.has_native_selnorm_kernel())

  # The middle bin's standardized interval [40, 41] underflows in probability
  # space, but its weight makes it dominate the upper tail.
  omega <- c(1, 1e300, 1e-100)
  prior <- BayesTools::prior_weightfunction(
    "one-sided", c(stats::pnorm(1, lower.tail = FALSE), .5), BayesTools::wf_fixed(omega)
  )
  for (direction in c("positive", "negative")) {
    q <- if (direction == "positive") -40 else 40
    spec <- suppressWarnings(.selection_spec(
      list(outcome = list(bias = prior)), yi = q, sei = 1,
      effect_direction = direction, signed_data = FALSE
    ))
    args <- list(mean = matrix(q), sd = matrix(1), sei = 1, omega = matrix(omega, 1L),
                 selection_spec = spec, kernel_mode = SELKERNEL_STEP)
    expected <- .predictive_cell_reference(q, q, 1, 1, omega,
                                           BayesTools::selection_native_static_args(spec))
    summary  <- do.call(.selnorm_kernel_predictive_summary_matrix, c(list(q = q), args))
    weighted <- do.call(.selnorm_kernel_weighted_summary,
                        c(list(yi = q, psis_weights = matrix(1)), args))
    moments  <- do.call(.selnorm_kernel_moments_matrix, args)
    for (tail in c("log_lower", "log_upper")) {
      label <- paste(direction, tail)
      expect_lt(abs(summary[[tail]][1L, 1L] - expected[[tail]]),
                1e-11 * max(1, abs(expected[[tail]])), label = label)
      expect_lt(abs(weighted[[tail]] - expected[[tail]]),
                1e-11 * max(1, abs(expected[[tail]])), label = paste("weighted", label))
    }
    # The heavily weighted bin is far in the signed upper tail.
    expect_lt(expected[[if (direction == "positive") "log_lower" else "log_upper"]], -100)
    expect_identical(summary$mean[1L, 1L], moments$mean[1L, 1L])
    expect_identical(summary$variance[1L, 1L], moments$variance[1L, 1L])
  }
})


test_that("predictive cells outside the fused gate return the base kernels' values", {

  skip_if_not(.has_native_selnorm_kernel())

  # Weights spanning beyond double precision, a query beyond the double range
  # and a subnormal query: each cell takes the unchanged log-tail and
  # partition-moment kernels, including their limits at such inputs.
  cases <- list(
    list(steps = c(stats::pnorm(1, lower.tail = FALSE), .5), omega = c(1, 1e308, 1),
         mean = 500, sd = 1, sei = 1000, q = 500 - 36.8),
    list(steps = c(stats::pnorm(1, lower.tail = FALSE), .5), omega = c(1, 1e300, 1e-100),
         mean = -4, sd = .1, sei = .1, q = 1e308),
    list(steps = .5, omega = c(1e308, 0), mean = 0, sd = 1, sei = 1,
         q = 1.9762625833649862e-323)
  )
  for (case in cases) for (direction in c("positive", "negative")) {
    sign  <- if (direction == "positive") 1 else -1
    prior <- BayesTools::prior_weightfunction(
      "one-sided", case$steps, BayesTools::wf_fixed(replace(pmin(case$omega, 1), 1L, 1))
    )
    spec <- suppressWarnings(.selection_spec(
      list(outcome = list(bias = prior)), yi = sign * case$q, sei = case$sei,
      effect_direction = direction, signed_data = FALSE
    ))
    args <- list(mean = matrix(sign * case$mean), sd = matrix(case$sd), sei = case$sei,
                 omega = matrix(case$omega, 1L), selection_spec = spec,
                 kernel_mode = SELKERNEL_STEP)
    label   <- paste(direction, "q =", case$q)
    summary <- do.call(.selnorm_kernel_predictive_summary_matrix, c(list(q = sign * case$q), args))
    moments <- do.call(.selnorm_kernel_moments_matrix, args)
    cdf     <- do.call(.selnorm_kernel_cdf_matrix, c(list(q = sign * case$q), args))
    expect_identical(summary$mean, moments$mean, label = paste(label, "mean"))
    expect_identical(summary$variance, moments$variance, label = paste(label, "variance"))
    expect_lt(abs(exp(summary$log_lower[1L, 1L]) - cdf[1L, 1L]), 1e-15, label = paste(label, "cdf"))
    # The reducer takes the same cell values: it fails exactly where the base
    # moments are not finite and otherwise returns the matrix tails.
    weighted <- tryCatch(
      do.call(.selnorm_kernel_weighted_summary,
              c(list(yi = sign * case$q, psis_weights = matrix(1)), args)),
      error = identity
    )
    if (is.finite(moments$mean[1L, 1L]) && is.finite(moments$variance[1L, 1L])) {
      expect_false(inherits(weighted, "error"), label = paste(label, "reducer"))
      for (tail in c("log_lower", "log_upper")) {
        expect_equal(weighted[[tail]], summary[[tail]][1L, 1L], tolerance = 1e-14,
                     label = paste(label, "reducer", tail))
      }
    } else {
      expect_s3_class(weighted, "error")
    }
  }
})


test_that("predictive cell moments stay available for bounded weights at tiny masses", {

  skip_if_not(.has_native_selnorm_kernel())

  # The mass sits in bins of weight 1e-290 while the unit-weight bin starts at
  # standardized scores 37.4 to 39, where its probability is subnormal or zero.
  sei    <- 1
  spec   <- .test_step_spec(0, sei)
  static <- BayesTools::selection_native_static_args(spec)
  omega  <- rep(1e-290, spec[["n_bins"]])
  omega[1L] <- 1
  first_cut <- static[["z_lower"]][1L] * sei
  for (score in c(37.4, 37.52, 37.6, 38, 38.5, 39)) {
    mean <- first_cut - score
    got  <- .selnorm_kernel_predictive_summary_matrix(
      q = mean, mean = matrix(mean), sd = matrix(1), sei = sei,
      omega = matrix(omega, 1L), selection_spec = spec, kernel_mode = SELKERNEL_STEP
    )
    expected <- .predictive_cell_reference(mean, mean, 1, sei, omega, static)
    expect_true(is.finite(got$mean[1L, 1L]) && is.finite(got$variance[1L, 1L]),
                label = paste("moments at score", score))
    expect_lt(abs(got$log_lower[1L, 1L] - expected$log_lower), 1e-11,
              label = paste("log_lower at score", score))
  }
})


test_that("predictive cell moments of narrow standardized intervals keep their accuracy", {

  skip_if_not(.has_native_selnorm_kernel())

  # Almost all mass sits in [0, sei] in an interval that is narrow on the
  # scale of a large sd, where shared tail differences cancel.
  prior <- BayesTools::prior_weightfunction(
    "one-sided", c(stats::pnorm(1, lower.tail = FALSE), .5),
    BayesTools::wf_fixed(c(1, 1e100, 0))
  )
  expect_warning(
    spec <- .selection_spec(list(outcome = list(bias = prior)), yi = .5, sei = 1,
                            effect_direction = "positive", signed_data = FALSE),
    "telescoping disabled"
  )
  for (sd in c(1, 10, 100, 1e3, 1e4, 1e5)) {
    args <- list(mean = matrix(0), sd = matrix(sd), sei = 1,
                 omega = matrix(c(1, 1e100, 0), 1L), selection_spec = spec,
                 kernel_mode = SELKERNEL_STEP)
    got <- do.call(.selnorm_kernel_predictive_summary_matrix, c(list(q = .5), args))
    old <- do.call(.selnorm_kernel_moments_matrix, args)
    kernel <- function(x) exp(-.5 * (x / sd)^2)
    mass   <- stats::integrate(kernel, 0, 1, rel.tol = 1e-13)$value
    first  <- stats::integrate(function(x) x * kernel(x), 0, 1, rel.tol = 1e-13)$value / mass
    second <- stats::integrate(function(x) (x - first)^2 * kernel(x), 0, 1,
                               rel.tol = 1e-13)$value / mass
    # Never further from the integrated truncated-normal moments than the
    # partition-moment kernel, and within its accuracy.
    tolerance <- max(1e-9, 2 * abs(old$variance[1L, 1L] - second) / second)
    expect_lt(abs(got$variance[1L, 1L] - second) / second, tolerance,
              label = paste("variance at sd", sd))
    expect_lt(abs(got$mean[1L, 1L] - first), max(1e-9, 2 * abs(old$mean[1L, 1L] - first)),
              label = paste("mean at sd", sd))
  }
})
