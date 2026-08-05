# ============================================================================ #
# Adaptive quadrature for ordinary GLMM likelihoods
# ============================================================================ #

.glmm_aghq_rule_cache <- new.env(parent = emptyenv())


# Return cached standard-normal Gauss-Hermite rules for adaptive refinement.
.glmm_aghq_rules <- function(orders = c(5L, 9L, 13L, 17L, 25L, 33L, 41L, 49L)) {

  orders <- as.integer(orders)
  key    <- paste(orders, collapse = ",")

  if (exists(key, envir = .glmm_aghq_rule_cache, inherits = FALSE)) {
    return(get(key, envir = .glmm_aghq_rule_cache, inherits = FALSE))
  }

  rules <- lapply(orders, .gauss_hermite_nodes)
  out   <- list(
    nodes       = lapply(rules, `[[`, "nodes"),
    log_weights = lapply(rules, `[[`, "log_weights")
  )

  assign(key, out, envir = .glmm_aghq_rule_cache)
  return(out)
}


# Check whether the native specialized quadrature is available for a family.
.has_native_glmm_aghq <- function(outcome_type, row_sum = FALSE) {

  suffix <- if (row_sum) "_row_sum" else ""
  symbol <- switch(
    outcome_type,
    "bin"  = paste0("RoBMA_glmm_binom_aghq", suffix),
    "pois" = paste0("RoBMA_glmm_pois_aghq", suffix),
    stop("Unknown GLMM outcome type.", call. = FALSE)
  )

  return(is.loaded(symbol, PACKAGE = "RoBMA"))
}


# Extract an untruncated nuisance-prior specification supported by AGHQ.
.glmm_aghq_prior_spec <- function(prior, outcome_type) {

  if (outcome_type == "bin") {
    if (!identical(prior[["distribution"]], "beta") ||
        !identical(as.numeric(prior[["truncation"]][["lower"]]), 0) ||
        !identical(as.numeric(prior[["truncation"]][["upper"]]), 1)) {
      return(NULL)
    }

    alpha <- as.numeric(prior[["parameters"]][["alpha"]])
    beta  <- as.numeric(prior[["parameters"]][["beta"]])
    if (length(alpha) != 1L || length(beta) != 1L ||
        !is.finite(alpha) || !is.finite(beta) || alpha <= 0 || beta <= 0) {
      return(NULL)
    }

    return(c(alpha = alpha, beta = beta))
  }

  if (outcome_type == "pois") {
    if (!identical(prior[["distribution"]], "normal") ||
        !is.infinite(prior[["truncation"]][["lower"]]) ||
        prior[["truncation"]][["lower"]] >= 0 ||
        !is.infinite(prior[["truncation"]][["upper"]]) ||
        prior[["truncation"]][["upper"]] <= 0) {
      return(NULL)
    }

    mean <- as.numeric(prior[["parameters"]][["mean"]])
    sd   <- as.numeric(prior[["parameters"]][["sd"]])
    if (length(mean) != 1L || length(sd) != 1L ||
        !is.finite(mean) || !is.finite(sd) || sd <= 0) {
      return(NULL)
    }

    return(c(mean = mean, sd = sd))
  }

  stop("Unknown GLMM outcome type.", call. = FALSE)
}


# Require the certified default marginal-likelihood route.
.glmm_aghq_require_default <- function(prior_spec, outcome_type,
                                       row_sum = FALSE) {

  if (is.null(prior_spec)) {
    prior_label <- if (identical(outcome_type, "bin")) {
      "an untruncated beta 'prior_pi'"
    } else {
      "an untruncated normal 'prior_phi'"
    }
    stop(
      "Default ordinary GLMM marginal log-likelihood evaluation requires ",
      prior_label, ". Truncated, point, and other nuisance priors are not ",
      "supported by the certified adaptive quadrature.",
      call. = FALSE
    )
  }
  if (!.has_native_glmm_aghq(outcome_type, row_sum = row_sum)) {
    stop(
      "The native adaptive ordinary-GLMM quadrature is unavailable.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Common native AGHQ controls. Optional arguments support focused numerical tests.
.glmm_aghq_control <- function(orders = NULL, tolerance = 1e-6,
                               consecutive = 2L, mode_tolerance = 1e-12) {

  if (is.null(orders)) {
    orders <- c(5L, 9L, 13L, 17L, 25L, 33L, 41L, 49L)
  }
  BayesTools::check_int(orders, "orders", lower = 3, check_length = 0)
  BayesTools::check_real(tolerance, "tolerance", lower = 0,
                         allow_bound = FALSE, check_length = 1)
  BayesTools::check_int(consecutive, "consecutive", lower = 1,
                        check_length = 1)
  BayesTools::check_real(mode_tolerance, "mode_tolerance", lower = 0,
                         allow_bound = FALSE, check_length = 1)

  orders <- as.integer(orders)
  if (is.unsorted(orders, strictly = TRUE) || any(orders %% 2L == 0L)) {
    stop("'orders' must be strictly increasing odd integers.", call. = FALSE)
  }
  if (length(orders) < consecutive + 1L) {
    stop("'orders' must contain enough rules for the convergence requirement.",
         call. = FALSE)
  }

  return(list(
    rules          = .glmm_aghq_rules(orders),
    tolerance      = as.numeric(tolerance),
    consecutive    = as.integer(consecutive),
    mode_tolerance = as.numeric(mode_tolerance)
  ))
}


# Independent theta and nuisance refinement for the prior-CDF fallback.
.glmm_grid_control <- function(
    theta_orders    = c(15L, 21L, 29L, 39L, 49L),
    nuisance_orders = c(30L, 50L, 75L, 100L, 150L, 200L, 300L, 450L,
                        600L, 800L, 1000L, 1200L),
    tolerance       = 1e-6,
    consecutive     = 2L) {

  BayesTools::check_int(theta_orders, "theta_orders", lower = 1,
                        check_length = 0)
  BayesTools::check_int(nuisance_orders, "nuisance_orders", lower = 1,
                        check_length = 0)
  BayesTools::check_real(tolerance, "tolerance", lower = 0,
                         allow_bound = FALSE, check_length = 1)
  BayesTools::check_int(consecutive, "consecutive", lower = 1,
                        check_length = 1)

  theta_orders    <- as.integer(theta_orders)
  nuisance_orders <- as.integer(nuisance_orders)
  consecutive     <- as.integer(consecutive)
  if (is.unsorted(theta_orders, strictly = TRUE) ||
      is.unsorted(nuisance_orders, strictly = TRUE)) {
    stop("Quadrature orders must be strictly increasing.", call. = FALSE)
  }
  if (length(theta_orders) < consecutive + 1L ||
      length(nuisance_orders) < consecutive + 1L) {
    stop("Each quadrature ladder must contain enough refinement rules.",
         call. = FALSE)
  }

  return(list(
    theta_orders    = theta_orders,
    nuisance_orders = nuisance_orders,
    tolerance       = as.numeric(tolerance),
    consecutive     = consecutive
  ))
}


.glmm_grid_change <- function(current, previous) {

  same  <- current == previous
  delta <- abs(current - previous)
  delta[!is.na(same) & same] <- 0
  if (anyNA(delta)) {
    return(Inf)
  }

  return(max(delta))
}


.glmm_grid_refine <- function(evaluate, outcome_type, observation, control,
                              theta_active = TRUE, nuisance_active = TRUE) {

  theta_orders <- if (theta_active) {
    control[["theta_orders"]]
  } else {
    1L
  }
  nuisance_orders <- if (nuisance_active) {
    control[["nuisance_orders"]]
  } else {
    1L
  }
  consecutive     <- control[["consecutive"]]
  theta_index     <- if (theta_active) consecutive + 1L else 1L
  nuisance_index  <- if (nuisance_active) consecutive + 1L else 1L
  cache            <- new.env(parent = emptyenv())

  value_at <- function(theta_index, nuisance_index) {

    key <- paste(theta_index, nuisance_index, sep = ":")
    if (!exists(key, envir = cache, inherits = FALSE)) {
      value <- evaluate(
        theta_orders[theta_index],
        nuisance_orders[nuisance_index]
      )
      assign(key, as.numeric(value), envir = cache)
    }
    return(get(key, envir = cache, inherits = FALSE))
  }

  repeat {
    current <- value_at(theta_index, nuisance_index)
    theta_changes <- if (theta_active) {
      vapply(seq_len(consecutive), function(offset) {
        .glmm_grid_change(
          value_at(theta_index - offset + 1L, nuisance_index),
          value_at(theta_index - offset, nuisance_index)
        )
      }, numeric(1))
    } else {
      0
    }
    nuisance_changes <- if (nuisance_active) {
      vapply(seq_len(consecutive), function(offset) {
        .glmm_grid_change(
          value_at(theta_index, nuisance_index - offset + 1L),
          value_at(theta_index, nuisance_index - offset)
        )
      }, numeric(1))
    } else {
      0
    }
    theta_ok    <- all(theta_changes <= control[["tolerance"]])
    nuisance_ok <- all(nuisance_changes <= control[["tolerance"]])

    if (theta_ok && nuisance_ok) {
      return(list(
        value            = current,
        theta_order      = theta_orders[theta_index],
        nuisance_order   = nuisance_orders[nuisance_index],
        max_change       = max(theta_changes, nuisance_changes)
      ))
    }

    if (theta_active && !theta_ok) {
      theta_index <- theta_index + 1L
    }
    if (nuisance_active && !nuisance_ok) {
      nuisance_index <- nuisance_index + 1L
    }
    if (theta_index > length(theta_orders) ||
        nuisance_index > length(nuisance_orders)) {
      label <- if (identical(outcome_type, "bin")) "Binomial" else "Poisson"
      stop(
        label, " prior-CDF quadrature failed to converge at observation ",
        observation, ". Last theta changes: ",
        paste(format(theta_changes, digits = 4), collapse = ", "),
        "; last nuisance changes: ",
        paste(format(nuisance_changes, digits = 4), collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
}


.glmm_aghq_numerical_failure <- function(condition, outcome_type) {

  label <- if (identical(outcome_type, "bin")) "Binomial" else "Poisson"
  pattern <- paste0("^", label, " AGHQ (mode failed|failed to converge)")
  return(grepl(pattern, conditionMessage(condition), perl = TRUE))
}


.glmm_aghq_dispatch <- function(S, K, outcome_type, evaluate_aghq = NULL,
                                evaluate_grid, row_sum = FALSE) {

  value               <- matrix(NA_real_, nrow = S, ncol = K)
  max_order           <- 0L
  max_change          <- 0
  max_mode_iterations <- 0L
  order_counts        <- integer(0)
  exact_count         <- 0L
  grid_columns        <- integer(0)
  grid_theta_order    <- 0L
  grid_nuisance_order <- 0L
  grid_max_change     <- 0

  for (k in seq_len(K)) {
    aghq <- if (is.null(evaluate_aghq)) {
      NULL
    } else {
      tryCatch(
        evaluate_aghq(k),
        error = function(condition) {
          if (!.glmm_aghq_numerical_failure(condition, outcome_type)) {
            stop(condition)
          }
          return(NULL)
        }
      )
    }

    if (is.null(aghq)) {
      grid                    <- evaluate_grid(k)
      value[, k]              <- grid[["value"]]
      grid_columns            <- c(grid_columns, k)
      grid_theta_order        <- max(grid_theta_order, grid[["theta_order"]])
      grid_nuisance_order     <- max(grid_nuisance_order, grid[["nuisance_order"]])
      grid_max_change         <- max(grid_max_change, grid[["max_change"]])
      next
    }

    value[, k]          <- aghq[["value"]][, 1L]
    max_order           <- max(max_order, aghq[["max_order"]])
    max_change          <- max(max_change, aghq[["max_change"]])
    max_mode_iterations <- max(max_mode_iterations,
                               aghq[["max_mode_iterations"]])
    if (length(order_counts) == 0L) {
      order_counts <- aghq[["order_counts"]]
    } else {
      order_counts <- order_counts + aghq[["order_counts"]]
    }
    exact_count <- exact_count + aghq[["exact_count"]]
  }

  if (row_sum) {
    value <- rowSums(value)
  }
  attr(value, "glmm_aghq_diagnostics") <- list(
    max_order            = max_order,
    max_change           = max_change,
    max_mode_iterations  = max_mode_iterations,
    order_counts         = order_counts,
    exact_count          = exact_count,
    grid_columns         = grid_columns,
    grid_max_theta_order = grid_theta_order,
    grid_max_nuisance_order = grid_nuisance_order,
    grid_max_change      = grid_max_change
  )
  return(value)
}


# Attach compact convergence diagnostics to a likelihood result.
.glmm_aghq_value <- function(out) {

  value <- out[["value"]]
  attr(value, "glmm_aghq_diagnostics") <- out[c(
    "max_order",
    "max_change",
    "max_mode_iterations",
    "order_counts",
    "exact_count"
  )]
  return(value)
}


# Native binomial AGHQ dispatch. Returns values and convergence diagnostics.
.glmm_binom_aghq <- function(ai, ci, n1i, n2i, mu_samples, tau_within,
                             weights, prior_spec, row_sum = FALSE,
                             control = .glmm_aghq_control()) {

  symbol <- if (row_sum) {
    "RoBMA_glmm_binom_aghq_row_sum"
  } else {
    "RoBMA_glmm_binom_aghq"
  }
  weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
    .glmm_likelihood_weights(weights, ncol(mu_samples))
  )

  return(.Call(
    symbol,
    .native_integer_vector(ai),
    .native_integer_vector(ci),
    .native_integer_vector(n1i),
    .native_integer_vector(n2i),
    .native_numeric_matrix(mu_samples),
    .native_numeric_matrix(tau_within),
    weights_arg,
    as.numeric(prior_spec[["alpha"]]),
    as.numeric(prior_spec[["beta"]]),
    control[["rules"]][["nodes"]],
    control[["rules"]][["log_weights"]],
    control[["tolerance"]],
    control[["consecutive"]],
    control[["mode_tolerance"]],
    PACKAGE = "RoBMA"
  ))
}


# Native Poisson AGHQ dispatch. Returns values and convergence diagnostics.
.glmm_pois_aghq <- function(x1i, x2i, t1i, t2i, mu_samples, tau_within,
                            weights, prior_spec, row_sum = FALSE,
                            control = .glmm_aghq_control()) {

  symbol <- if (row_sum) {
    "RoBMA_glmm_pois_aghq_row_sum"
  } else {
    "RoBMA_glmm_pois_aghq"
  }
  weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
    .glmm_likelihood_weights(weights, ncol(mu_samples))
  )

  return(.Call(
    symbol,
    .native_integer_vector(x1i),
    .native_integer_vector(x2i),
    .native_numeric_vector(t1i),
    .native_numeric_vector(t2i),
    .native_numeric_matrix(mu_samples),
    .native_numeric_matrix(tau_within),
    weights_arg,
    as.numeric(prior_spec[["mean"]]),
    as.numeric(prior_spec[["sd"]]),
    control[["rules"]][["nodes"]],
    control[["rules"]][["log_weights"]],
    control[["tolerance"]],
    control[["consecutive"]],
    control[["mode_tolerance"]],
    PACKAGE = "RoBMA"
  ))
}


.glmm_binom_marginal <- function(
    ai, ci, n1i, n2i, mu_samples, tau_within, weights, prior_pi,
    prior_spec, row_sum = FALSE, control = .glmm_aghq_control(),
    grid_control = .glmm_grid_control()) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  evaluate_aghq <- NULL
  if (!is.null(prior_spec)) {
    .glmm_aghq_require_default(prior_spec, "bin")
    evaluate_aghq <- function(k) {
      .glmm_binom_aghq(
        ai          = ai[k],
        ci          = ci[k],
        n1i         = n1i[k],
        n2i         = n2i[k],
        mu_samples  = mu_samples[, k, drop = FALSE],
        tau_within  = tau_within[, k, drop = FALSE],
        weights     = if (is.null(weights)) NULL else weights[k],
        prior_spec  = prior_spec,
        control     = control
      )
    }
  }
  evaluate_grid <- function(k) {
    .glmm_grid_refine(
      evaluate = function(n_theta, n_nuisance) {
        .outcome_pdf.binom(
          ai          = ai[k],
          ci          = ci[k],
          n1i         = n1i[k],
          n2i         = n2i[k],
          mu_samples  = mu_samples[, k, drop = FALSE],
          tau_within  = tau_within[, k, drop = FALSE],
          prior_pi    = prior_pi,
          weights     = if (is.null(weights)) NULL else weights[k],
          n_theta     = n_theta,
          n_pi        = n_nuisance
        )
      },
      outcome_type    = "bin",
      observation     = k,
      control         = grid_control,
      theta_active    = any(tau_within[, k] != 0),
      nuisance_active = !BayesTools::is.prior.point(prior_pi)
    )
  }

  return(.glmm_aghq_dispatch(
    S             = S,
    K             = K,
    outcome_type  = "bin",
    evaluate_aghq = evaluate_aghq,
    evaluate_grid = evaluate_grid,
    row_sum       = row_sum
  ))
}


.glmm_pois_marginal <- function(
    x1i, x2i, t1i, t2i, mu_samples, tau_within, weights, prior_phi,
    prior_spec, row_sum = FALSE, control = .glmm_aghq_control(),
    grid_control = .glmm_grid_control()) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  evaluate_aghq <- NULL
  if (!is.null(prior_spec)) {
    .glmm_aghq_require_default(prior_spec, "pois")
    evaluate_aghq <- function(k) {
      .glmm_pois_aghq(
        x1i        = x1i[k],
        x2i        = x2i[k],
        t1i        = t1i[k],
        t2i        = t2i[k],
        mu_samples = mu_samples[, k, drop = FALSE],
        tau_within = tau_within[, k, drop = FALSE],
        weights    = if (is.null(weights)) NULL else weights[k],
        prior_spec = prior_spec,
        control    = control
      )
    }
  }
  evaluate_grid <- function(k) {
    .glmm_grid_refine(
      evaluate = function(n_theta, n_nuisance) {
        .outcome_pdf.pois(
          x1i        = x1i[k],
          x2i        = x2i[k],
          t1i        = t1i[k],
          t2i        = t2i[k],
          mu_samples = mu_samples[, k, drop = FALSE],
          tau_within = tau_within[, k, drop = FALSE],
          prior_phi  = prior_phi,
          weights    = if (is.null(weights)) NULL else weights[k],
          n_theta    = n_theta,
          n_phi      = n_nuisance
        )
      },
      outcome_type    = "pois",
      observation     = k,
      control         = grid_control,
      theta_active    = any(tau_within[, k] != 0),
      nuisance_active = !BayesTools::is.prior.point(prior_phi)
    )
  }

  return(.glmm_aghq_dispatch(
    S             = S,
    K             = K,
    outcome_type  = "pois",
    evaluate_aghq = evaluate_aghq,
    evaluate_grid = evaluate_grid,
    row_sum       = row_sum
  ))
}
