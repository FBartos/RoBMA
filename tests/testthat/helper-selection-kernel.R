.test_step_spec <- function(yi, sei, effect_direction = "positive") {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05, .50),
    weights = BayesTools::wf_fixed(c(1, .7, .35, .2))
  )

  .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = effect_direction,
    signed_data      = FALSE
  )
}

.test_two_sided_step_spec <- function(yi, sei, effect_direction = "positive") {

  prior <- BayesTools::prior_weightfunction(
    side    = "two-sided",
    steps   = c(.05, .10),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )

  .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = effect_direction,
    signed_data      = FALSE
  )
}

.test_step_log_norm_reference <- function(mean, sd, sei, omega, spec) {

  S   <- nrow(mean)
  K   <- ncol(mean)
  out <- matrix(NA_real_, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      mean_z   <- spec[["sign"]] * mean[s, k] / sei[k]
      sd_z     <- sd[s, k] / sei[k]
      log_mass <- vapply(seq_len(spec[["n_bins"]]), function(b) {
        log(omega[s, b]) + .test_interval_log_prob(
          spec[["z_lower"]][b],
          spec[["z_upper"]][b],
          mean_z,
          sd_z
        )
      }, numeric(1))

      out[s, k] <- .test_logsumexp(log_mass)
    }
  }

  return(out)
}

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

.test_interval_prob_vec <- function(lower, upper, mean, sd) {

  lower_tail <- stats::pnorm(upper, mean = mean, sd = sd) -
    stats::pnorm(lower, mean = mean, sd = sd)
  upper_tail <- stats::pnorm(
    lower,
    mean       = mean,
    sd         = sd,
    lower.tail = FALSE
  ) - stats::pnorm(
    upper,
    mean       = mean,
    sd         = sd,
    lower.tail = FALSE
  )
  out <- ifelse(lower >= mean, upper_tail, lower_tail)
  out[lower >= upper] <- 0
  if (any(!is.finite(out)) || any(out < 0) || any(out > 1)) {
    stop("Selection interval probability reference failed.", call. = FALSE)
  }
  return(out)
}

.test_validate_cdf <- function(value, context, operation_count = 1L) {

  if (!is.numeric(value) || any(!is.finite(value))) {
    stop(context, " produced invalid CDF values.", call. = FALSE)
  }

  relative_error <- operation_count * .Machine$double.eps /
    (1 - operation_count * .Machine$double.eps)
  near_zero <- value < 0 & value >= -relative_error
  near_one  <- value > 1 & value <= 1 + relative_error
  value[near_zero] <- 0
  value[near_one]  <- 1

  if (any(value < 0 | value > 1)) {
    stop(context, " produced invalid CDF values.", call. = FALSE)
  }
  return(value)
}

.test_selection_mixture_has_full_support <- function(selection_context,
                                                     selected_rows) {

  if (is.null(selection_context) || any(!selected_rows)) {
    return(TRUE)
  }

  omega <- selection_context[["omega"]][selected_rows, , drop = FALSE]
  return(all(colSums(omega > 0) > 0))
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
