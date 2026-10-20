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
    mean     = moment_mean,
    second   = moment_second,
    variance = moment_second - moment_mean^2
  ))
}


# OpenMP limits and the complete private library list precede namespace import.
# Every numerical comparison executes in the same fresh two-core worker, which
# loads the RoBMA under test: the checkout when the tests run under pkgload,
# the installed namespace of the tests otherwise. The worker's native library
# must be the parent's; a different build would test other code.
.with_kernel_thread_worker <- function(code) {

  withr::local_envvar(c(OMP_NUM_THREADS = "2", OMP_THREAD_LIMIT = "2",
    OMP_DYNAMIC = "FALSE", R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)))
  cluster <- parallel::makePSOCKcluster(1L, rscript_args = "--vanilla")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  namespace_path <- getNamespaceInfo(asNamespace("RoBMA"), "path")
  development    <- requireNamespace("pkgload", quietly = TRUE) &&
    pkgload::is_dev_package("RoBMA")
  parent_dll     <- normalizePath(getLoadedDLLs()[["RoBMA"]][["path"]],
    winslash = "/", mustWork = TRUE)
  load_worker <- function(paths, development, namespace_path) {
    .libPaths(paths)
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    if (development) {
      pkgload::load_all(namespace_path, compile = FALSE, helpers = FALSE, quiet = TRUE)
    } else {
      loadNamespace("RoBMA", lib.loc = dirname(namespace_path))
    }
    RoBMA::RoBMA.options(max_cores = 2L, native_threads = 1L)
    list(
      dll    = normalizePath(getLoadedDLLs()[["RoBMA"]][["path"]],
        winslash = "/", mustWork = TRUE),
      openmp = .Call("RoBMA_native_unwind_probe", "openmp", NULL, 1L, PACKAGE = "RoBMA")
    )
  }
  # Under pkgload this helper descends from the RoBMA namespace; unserializing
  # such a closure would load the installed package in the worker before the
  # build under test, so the function is sent without its defining environment.
  environment(load_worker) <- globalenv()
  loaded <- parallel::clusterCall(cluster, load_worker,
    .libPaths(), development, namespace_path)[[1L]]
  same_build <- if (.Platform$OS.type == "windows") {
    identical(tolower(loaded[["dll"]]), tolower(parent_dll))
  } else identical(loaded[["dll"]], parent_dll)
  if (!same_build) {
    stop("The kernel-thread worker loaded a different RoBMA build than the tests: ",
      "worker library '", loaded[["dll"]], "', tests '", parent_dll, "'.", call. = FALSE)
  }
  testthat::skip_if_not(isTRUE(loaded[["openmp"]]), "This build has no OpenMP support.")

  result <- parallel::clusterCall(cluster, function(code, helper) {
    env <- new.env(parent = asNamespace("RoBMA"))
    for (name in c("expect_identical", "expect_true", "expect_error", "expect_gt")) {
      env[[name]] <- getExportedValue("testthat", name)
    }
    sys.source(helper, envir = env)
    eval(code, envir = env)
    TRUE
  }, substitute(code), normalizePath(testthat::test_path("helper-selection-kernel.R")))[[1L]]
  testthat::expect_true(result)
}


# One synthetic batch with row-wise selection weights, kernel modes and scales.
# Callers choose S from the actual entry point's work formula.
.kernel_thread_inputs <- function(S = 96L, K = 5L) {

  withr::local_seed(20260917)
  yi   <- c(.12, -.05, .31, .22, -.18)[seq_len(K)]
  sei  <- c(.10, .18, .14, .22, .09)[seq_len(K)]
  spec <- .test_step_spec(yi, sei)

  mean <- matrix(stats::rnorm(S * K, mean = .15, sd = .10), nrow = S)
  sd   <- matrix(stats::runif(S * K, .08, .30), nrow = S)
  # A normal row every fourth draw: the selected and unselected branches of the
  # kernels must both be visited inside the parallel region.
  kernel_mode <- rep(SELKERNEL_STEP, S)
  kernel_mode[seq.int(4L, S, by = 4L)] <- SELKERNEL_NORMAL
  omega <- cbind(
    rep(1, S),
    stats::runif(S, .40, .95),
    stats::runif(S, .15, .60),
    stats::runif(S, .05, .40)
  )
  omega <- omega[, seq_len(spec[["n_bins"]]), drop = FALSE]
  omega[, 1L] <- 1

  return(list(
    yi          = yi,
    sei         = sei,
    spec        = spec,
    S           = S,
    K           = K,
    mean        = mean,
    sd          = sd,
    omega       = omega,
    kernel_mode = kernel_mode
  ))
}


.kernel_thread_selection_context <- function(inputs) {

  context <- inputs[["spec"]]
  context[["omega"]]       <- inputs[["omega"]]
  context[["alpha"]]       <- 0
  context[["phack_kind"]]  <- 0L
  context[["kernel_mode"]] <- inputs[["kernel_mode"]]
  context[["vector_rule"]] <- 0L
  context[["use_normal"]]  <- inputs[["kernel_mode"]] == SELKERNEL_NORMAL

  return(BayesTools::selection_context_validate(
    context,
    n_samples = inputs[["S"]]
  ))
}


# Evaluate one entry point at several thread budgets and require the values to
# be identical, not merely close: the row partition must never change a result.
.expect_thread_invariant <- function(label, evaluate,
                                     threads = c(1L, 2L, 8L),
                                     rows = NULL, work_per_row = NULL) {

  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)

  reference <- NULL
  for (count in threads) {
    RoBMA.options(native_threads = count)
    if (count > 1L && !is.null(rows)) {
      schedule <- .Call("RoBMA_selnorm_row_schedule", as.double(rows),
        as.double(work_per_row), PACKAGE = "RoBMA")
      expect_true(schedule[["threads"]] > 1L,
        info = paste(label, "rows", rows, "work per row", work_per_row))
    }
    values <- evaluate()
    if (is.null(reference)) {
      expect_true(
        any(is.finite(unlist(values, use.names = FALSE))),
        info = paste0(label, ": no finite value to compare")
      )
      reference <- values
    } else {
      expect_identical(
        values,
        reference,
        info = paste0(label, " at ", count, " threads")
      )
    }
  }

  return(invisible(reference))
}

.native_selection_jags_fixtures <- function() {

  base <- list(y = c(0, 0), mu = c(0, 0), covariance = c(1, 0, 1),
    sd = c(1, 1), loading = c(.1, .1), sei = c(1, 1), omega = c(1, .2),
    lower = c(0, -1e300), upper = c(1e300, 0), bins = c(1, 1),
    sign = 1, telescope = 1, mode = 0, points = 4, maximum = 4,
    scrambles = 2, tolerance = .005, rule = 0,
    nodes = rep(0, 9), log_weights = rep(0, 9), orders = c(1, 3, 5))
  cases <- list(
    dense = list(distribution = "dselnorm_mnorm_step", arguments = c(
      "mu", "covariance", "sei", "omega", "lower", "upper", "bins",
      "sign", "telescope", "mode", "qmc", "points", "scrambles", "tolerance",
      "rule", "nodes", "log_weights", "orders"),
      extra = list(qmc = rep(.5, 32)), bad = c("mode", "sign", "telescope", "bins", "points")),
    cluster = list(distribution = "dselnorm_cluster_step", arguments = c(
      "mu", "sd", "loading", "sei", "omega", "lower", "upper", "bins",
      "sign", "telescope", "mode", "nodes", "log_weights", "orders", "qmc",
      "points", "maximum", "scrambles", "tolerance", "rule"),
      extra = list(qmc = rep(.5, 16)), bad = c("mode", "sign", "telescope", "bins")),
    factor = list(distribution = "dselnorm_factor_step", arguments = c(
      "mu", "sd", "factor_loading", "sei", "omega", "lower", "upper", "bins",
      "sign", "telescope", "mode", "nodes", "log_weights", "orders", "qmc",
      "points", "maximum", "scrambles", "tolerance", "rule"),
      extra = list(factor_loading = rep(.1, 4), qmc = rep(.5, 32)),
      bad = c("mode", "bins", "points")),
    conditioned = list(distribution = "dselnorm_sampling_conditioned", arguments = c(
      "mu", "covariance", "residual", "loading", "rank", "mu", "mu", "sei", "omega",
      "lower", "upper", "bins", "sign", "mode", "telescope", "rule", "groups", "qmc",
      "points", "maximum", "scrambles", "tolerance", "nodes", "log_weights", "orders",
      "nodes", "log_weights", "orders"),
      extra = list(residual = c(0, 0), rank = 1, groups = c(1, 2), qmc = rep(.5, 32)),
      bad = c("rank", "bins", "groups", "rule", "scrambles"))
  )
  list(base = base, cases = cases)
}
