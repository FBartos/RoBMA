# ============================================================================ #
# Adaptive quadrature for ordinary GLMM likelihoods
# ============================================================================ #

.glmm_aghq_rule_cache <- new.env(parent = emptyenv())


# Return cached standard-normal Gauss-Hermite rules for adaptive refinement.
.glmm_aghq_rules <- function(orders = c(5L, 9L, 13L, 17L, 25L, 33L, 41L, 49L, 65L)) {

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
    orders <- c(5L, 9L, 13L, 17L, 25L, 33L, 41L, 49L, 65L)
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
