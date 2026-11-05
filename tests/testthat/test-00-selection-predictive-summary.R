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
