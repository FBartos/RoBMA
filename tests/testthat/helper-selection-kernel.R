.test_step_reference <- function(yi, mu, sigma, sei, omega, spec,
                                 weights = rep(1, length(yi))) {

  S    <- nrow(mu)
  K    <- ncol(mu)
  sign <- spec[["sign"]]
  out  <- matrix(NA_real_, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      signed_mean <- sign * mu[s, k]
      signed_y    <- sign * yi[k]
      denom       <- 0

      for (b in seq_len(spec[["n_bins"]])) {
        lower <- spec[["z_lower"]][b] * sei[k]
        upper <- spec[["z_upper"]][b] * sei[k]
        prob  <- stats::pnorm(upper, mean = signed_mean, sd = sigma[s, k]) -
          stats::pnorm(lower, mean = signed_mean, sd = sigma[s, k])
        denom <- denom + omega[s, b] * prob
      }

      out[s, k] <- weights[k] * (
        stats::dnorm(signed_y, mean = signed_mean, sd = sigma[s, k], log = TRUE) +
          log(omega[s, spec[["obs_bin"]][k]]) - log(denom)
      )
    }
  }

  return(out)
}

.test_interval_prob <- function(lower, upper, mean, sd) {

  if (!(sd > 0) || lower >= upper) {
    return(0)
  }
  if (is.infinite(lower) && lower < 0) {
    return(stats::pnorm(upper, mean = mean, sd = sd))
  }
  if (is.infinite(upper) && upper > 0) {
    return(stats::pnorm(lower, mean = mean, sd = sd, lower.tail = FALSE))
  }
  return(stats::pnorm(upper, mean = mean, sd = sd) -
    stats::pnorm(lower, mean = mean, sd = sd))
}

.test_logspace_sub <- function(log_a, log_b) {

  if (is.infinite(log_b) && log_b < 0) {
    return(log_a)
  }
  if (log_b > log_a) {
    return(NaN)
  }

  return(log_a + log1p(-exp(log_b - log_a)))
}

.test_interval_log_prob <- function(lower, upper, mean, sd) {

  if (!(sd > 0) || lower >= upper) {
    return(-Inf)
  }
  if (is.infinite(lower) && lower < 0 && is.infinite(upper) && upper > 0) {
    return(0)
  }

  if (lower >= mean) {
    log_s_lower <- if (is.infinite(lower) && lower > 0) {
      -Inf
    } else {
      stats::pnorm(lower, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE)
    }
    log_s_upper <- if (is.infinite(upper) && upper > 0) {
      -Inf
    } else {
      stats::pnorm(upper, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE)
    }
    return(.test_logspace_sub(log_s_lower, log_s_upper))
  }

  if (upper <= mean) {
    log_f_upper <- if (is.infinite(upper) && upper < 0) {
      -Inf
    } else {
      stats::pnorm(upper, mean = mean, sd = sd, lower.tail = TRUE, log.p = TRUE)
    }
    log_f_lower <- if (is.infinite(lower) && lower < 0) {
      -Inf
    } else {
      stats::pnorm(lower, mean = mean, sd = sd, lower.tail = TRUE, log.p = TRUE)
    }
    return(.test_logspace_sub(log_f_upper, log_f_lower))
  }

  prob <- .test_interval_prob(lower, upper, mean, sd)
  return(if (prob > 0) log(prob) else -Inf)
}

.test_logsumexp <- function(x) {

  max_x <- max(x)
  if (!is.finite(max_x)) {
    return(max_x)
  }

  return(max_x + log(sum(exp(x - max_x))))
}

.test_step_cdf_reference <- function(q, mu, sigma, sei, omega, spec,
                                     kernel_mode = rep(SELKERNEL_STEP, nrow(mu)),
                                     lower.tail = TRUE) {

  S <- nrow(mu)
  K <- ncol(mu)
  if (length(q) == 1L) {
    q <- rep(q, K)
  }

  out <- matrix(NA_real_, nrow = S, ncol = K)
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      if (kernel_mode[s] == SELKERNEL_NORMAL) {
        out[s, k] <- stats::pnorm(
          q[k],
          mean       = mu[s, k],
          sd         = sigma[s, k],
          lower.tail = lower.tail
        )
        next
      }

      mean_z     <- spec[["sign"]] * mu[s, k] / sei[k]
      sd_z       <- sigma[s, k] / sei[k]
      q_z        <- spec[["sign"]] * q[k] / sei[k]
      denom      <- 0
      lower_mass <- 0
      upper_mass <- 0

      for (b in seq_len(spec[["n_bins"]])) {
        lower <- spec[["z_lower"]][b]
        upper <- spec[["z_upper"]][b]
        denom <- denom + omega[s, b] *
          .test_interval_prob(lower, upper, mean_z, sd_z)
        lower_mass <- lower_mass + omega[s, b] *
          .test_interval_prob(lower, min(upper, q_z), mean_z, sd_z)
        upper_mass <- upper_mass + omega[s, b] *
          .test_interval_prob(max(lower, q_z), upper, mean_z, sd_z)
      }

      if (spec[["sign"]] == 1L) {
        out[s, k] <- if (lower.tail) lower_mass / denom else upper_mass / denom
      } else {
        out[s, k] <- if (lower.tail) upper_mass / denom else lower_mass / denom
      }
    }
  }

  return(out)
}

.test_step_moments_reference <- function(mu, sigma, sei, omega, spec,
                                         kernel_mode = rep(SELKERNEL_STEP, nrow(mu))) {

  S             <- nrow(mu)
  K             <- ncol(mu)
  moment_mean   <- matrix(NA_real_, nrow = S, ncol = K)
  moment_second <- matrix(NA_real_, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      if (kernel_mode[s] == SELKERNEL_NORMAL) {
        moment_mean[s, k]   <- mu[s, k]
        moment_second[s, k] <- sigma[s, k]^2 + mu[s, k]^2
        next
      }

      signed_mean   <- spec[["sign"]] * mu[s, k]
      denom         <- 0
      signed_first  <- 0
      signed_second <- 0

      for (b in seq_len(spec[["n_bins"]])) {
        lower <- spec[["z_lower"]][b] * sei[k]
        upper <- spec[["z_upper"]][b] * sei[k]
        m0    <- .test_interval_prob(lower, upper, signed_mean, sigma[s, k])
        a     <- (lower - signed_mean) / sigma[s, k]
        b_z   <- (upper - signed_mean) / sigma[s, k]
        pa    <- stats::dnorm(a)
        pb    <- stats::dnorm(b_z)
        a_pa  <- a * pa
        b_pb  <- b_z * pb
        if (!is.finite(a_pa)) {
          a_pa <- 0
        }
        if (!is.finite(b_pb)) {
          b_pb <- 0
        }

        m1 <- signed_mean * m0 + sigma[s, k] * (pa - pb)
        m2 <- (signed_mean^2 + sigma[s, k]^2) * m0 +
          2 * signed_mean * sigma[s, k] * (pa - pb) +
          sigma[s, k]^2 * (a_pa - b_pb)

        denom         <- denom + omega[s, b] * m0
        signed_first  <- signed_first + omega[s, b] * m1
        signed_second <- signed_second + omega[s, b] * m2
      }

      moment_mean[s, k]   <- spec[["sign"]] * signed_first / denom
      moment_second[s, k] <- signed_second / denom
    }
  }

  return(list(
    mean   = moment_mean,
    second = moment_second
  ))
}
