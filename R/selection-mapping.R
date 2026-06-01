# ============================================================================ #
# Selection-Model Mapping Helpers
# ============================================================================ #
#
# Selected-normal kernels use BayesTools p-bin order directly. R code prepares
# compact kernel metadata; likelihood probability loops stay in native code.
#
# ============================================================================ #


SELKERNEL_NORMAL           <- 0L
SELKERNEL_STEP             <- 1L
SELKERNEL_PHACK_POWER      <- 2L
SELKERNEL_STEP_PHACK_POWER <- 3L


.is_prior_weightfunction <- function(x) {

  return(!is.null(x) && BayesTools::is.prior.weightfunction(x))
}

.is_prior_phacking <- function(x) {

  return(!is.null(x) && BayesTools::is_prior_phacking(x))
}

.is_prior_bias_kernel <- function(x) {

  return(!is.null(x) && BayesTools::is_prior_bias(x))
}

.prior_is_selection_kernel <- function(prior) {

  return(
    .is_prior_weightfunction(prior) ||
      .is_prior_phacking(prior) ||
      .is_prior_bias_kernel(prior)
  )
}

.prior_has_phacking <- function(prior) {

  if (is.null(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(any(vapply(as.list(prior), .prior_has_phacking, logical(1))))
  }
  if (.is_prior_weightfunction(prior)) {
    return(FALSE)
  }
  if (!.prior_is_selection_kernel(prior)) {
    return(FALSE)
  }

  return(.selection_backend_has_phacking(.selection_backend_spec(list(prior))))
}

.prior_has_selection <- function(prior) {

  if (is.null(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(any(vapply(as.list(prior), .prior_has_selection, logical(1))))
  }
  if (.is_prior_weightfunction(prior)) {
    return(TRUE)
  }
  if (!.prior_is_selection_kernel(prior)) {
    return(FALSE)
  }

  backend <- .selection_backend_spec(list(prior))
  return(backend[["mode"]] %in% c("step", "step_phack_power"))
}

.selection_stop_phacking_deferred <- function() {

  stop(
    "P-hacking selection kernels are not enabled in this RoBMA release.",
    call. = FALSE
  )
}

.selection_bias_priors <- function(priors) {

  priors_bias <- priors[["outcome"]][["bias"]]
  if (is.null(priors_bias)) {
    return(list())
  }
  if (BayesTools::is.prior.mixture(priors_bias)) {
    out <- as.list(priors_bias)
    class(out) <- "list"
    return(out)
  }
  return(list(priors_bias))
}

.selection_backend_spec <- function(priors_bias) {

  return(BayesTools::selection_backend_spec(priors_bias, backend = "jags"))
}

.selection_backend_has_phacking <- function(backend) {

  return(
    backend[["mode"]] %in% c("phack_power", "step_phack_power")
  )
}

.selection_kernel_breaks <- function(priors) {

  if (is.null(priors)) {
    return(NULL)
  }

  priors_list <- if (BayesTools::is.prior.mixture(priors)) {
    as.list(priors)
  } else if (BayesTools::is.prior(priors)) {
    list(priors)
  } else {
    priors
  }

  keep <- vapply(priors_list, .prior_has_selection, logical(1))
  if (!any(keep)) {
    return(NULL)
  }

  backend <- .selection_backend_spec(priors_list[keep])
  breaks  <- backend[["step"]][["breaks"]]
  if (is.null(breaks)) {
    breaks <- c(0, 1)
  }
  return(.selection_assert_p_cuts(breaks))
}

.selection_kernel_mode <- function(mode) {

  switch(
    mode,
    "normal"            = SELKERNEL_NORMAL,
    "step"              = SELKERNEL_STEP,
    "phack_power"       = SELKERNEL_PHACK_POWER,
    "step_phack_power"  = SELKERNEL_STEP_PHACK_POWER,
    stop("Unsupported selection kernel mode.", call. = FALSE)
  )
}

.selection_branch_kernel_modes <- function(backend) {

  branch_type <- backend[["branch_type"]]
  if (is.null(branch_type)) {
    return(.selection_kernel_mode(backend[["mode"]]))
  }

  kernel_mode <- rep(SELKERNEL_NORMAL, length(branch_type))
  has_step  <- branch_type %in% c("weightfunction", "combined")
  has_phack <- branch_type %in% c("phack", "phacking", "combined")

  if (backend[["mode"]] == "step") {
    kernel_mode[has_step] <- SELKERNEL_STEP
  } else if (backend[["mode"]] == "phack_power") {
    kernel_mode[has_phack] <- SELKERNEL_PHACK_POWER
  } else if (backend[["mode"]] == "step_phack_power") {
    kernel_mode[has_step & !has_phack]  <- SELKERNEL_STEP
    kernel_mode[has_phack & !has_step]  <- SELKERNEL_PHACK_POWER
    kernel_mode[has_step & has_phack]   <- SELKERNEL_STEP_PHACK_POWER
  }

  return(as.integer(kernel_mode))
}

.selection_jags_kernel_mode_expression <- function(branch_kernel_mode) {

  if (length(branch_kernel_mode) == 1L) {
    return(as.character(branch_kernel_mode[1L]))
  }

  terms <- paste0(
    branch_kernel_mode,
    " * equals(bias_indicator, ",
    seq_along(branch_kernel_mode),
    ")"
  )
  return(paste(terms, collapse = " + "))
}

.selection_obs_bin <- function(yi, sei, p_cuts, sign) {

  p_value <- stats::pnorm(sign * yi / sei, lower.tail = FALSE)
  return(.selection_p_bin(p_value, p_cuts))
}

.selection_step_bin_from_z <- function(z, p_cuts) {

  p_value <- stats::pnorm(z, lower.tail = FALSE)
  return(.selection_p_bin(p_value, p_cuts))
}

.selection_p_bin <- function(p_value, p_cuts) {

  p_cuts <- .selection_assert_p_cuts(p_cuts)

  for (cut in p_cuts) {
    close <- !is.na(p_value) & abs(p_value - cut) <= 1e-12
    p_value[close] <- cut
  }
  bin     <- findInterval(p_value, p_cuts, rightmost.closed = TRUE, left.open = TRUE)
  bin     <- pmin(pmax(bin, 1L), length(p_cuts) - 1L)
  return(as.integer(bin))
}

.selection_segment_midpoint <- function(lower, upper) {

  if (is.infinite(lower) && lower < 0) {
    return(upper - 1)
  }
  if (is.infinite(upper) && upper > 0) {
    return(lower + 1)
  }
  return((lower + upper) / 2)
}

.selection_segments <- function(p_cuts, phack_z_source = c(0, 0),
                                phack_z_dest = c(0, 0), has_phacking = FALSE) {

  p_cuts <- .selection_assert_p_cuts(p_cuts)

  z_lower <- stats::qnorm(1 - p_cuts[-1])
  z_upper <- stats::qnorm(1 - p_cuts[-length(p_cuts)])
  bounds  <- c(-Inf, Inf, z_lower[is.finite(z_lower)], z_upper[is.finite(z_upper)])

  if (has_phacking) {
    bounds <- c(bounds, phack_z_source, phack_z_dest)
  }

  bounds <- sort(unique(bounds))
  n_segments <- length(bounds) - 1L
  step_bin   <- integer(n_segments)
  region     <- integer(n_segments)

  for (i in seq_len(n_segments)) {
    mid <- .selection_segment_midpoint(bounds[i], bounds[i + 1L])
    step_bin[i] <- .selection_step_bin_from_z(mid, p_cuts)

    if (has_phacking) {
      if (mid >= phack_z_source[1] && mid <= phack_z_source[2]) {
        region[i] <- 1L
      } else if (mid > phack_z_dest[1] && mid <= phack_z_dest[2]) {
        region[i] <- 2L
      }
    }
  }

  return(list(
    bounds       = bounds,
    step_bin     = step_bin,
    phack_region = region
  ))
}

.selection_jags_bounds <- function(x) {

  x[is.infinite(x) & x < 0] <- -1e300
  x[is.infinite(x) & x > 0] <-  1e300
  return(x)
}

.selection_assert_p_cuts <- function(p_cuts) {

  if (!is.numeric(p_cuts) || length(p_cuts) < 2L ||
      anyNA(p_cuts) || any(!is.finite(p_cuts))) {
    stop("Selection p-value cut breaks must be finite numeric values.", call. = FALSE)
  }
  if (abs(p_cuts[1L]) > 1e-12 || abs(p_cuts[length(p_cuts)] - 1) > 1e-12) {
    stop("Selection p-value cut breaks must include endpoints 0 and 1.", call. = FALSE)
  }
  if (any(diff(p_cuts) <= 0)) {
    stop("Selection p-value cut breaks must be strictly increasing.", call. = FALSE)
  }

  p_cuts[1L]             <- 0
  p_cuts[length(p_cuts)] <- 1
  return(as.numeric(p_cuts))
}

.selection_fixed_omega_by_branch <- function(priors_bias, p_cuts) {

  n_bins <- length(p_cuts) - 1L
  if (n_bins <= 0L) {
    return(matrix(numeric(), nrow = length(priors_bias), ncol = 0L))
  }

  rows <- lapply(priors_bias, function(prior) {
    .selection_fixed_omega_branch(prior, p_cuts)
  })
  out <- do.call(rbind, rows)
  colnames(out) <- paste0("omega[", seq_len(n_bins), "]")

  return(out)
}

.selection_fixed_omega_branch <- function(prior, p_cuts) {

  n_bins <- length(p_cuts) - 1L
  selection <- NULL

  if (BayesTools::is.prior.weightfunction(prior)) {
    selection <- prior
  } else if (BayesTools::is_prior_bias(prior) &&
             !is.null(prior[["selection"]])) {
    selection <- prior[["selection"]]
  }

  if (is.null(selection)) {
    return(rep(1, n_bins))
  }

  weights <- selection[["weights"]]
  if (is.null(weights) || !identical(weights[["type"]], "fixed")) {
    return(rep(NA_real_, n_bins))
  }

  mapping <- BayesTools::weightfunctions_mapping(
    prior_list = list(selection),
    one_sided  = TRUE
  )[[1L]]

  return(as.numeric(weights[["omega"]][mapping]))
}

.selection_prior_quantile <- function(prior, probability) {

  distribution <- prior[["distribution"]]
  parameters   <- prior[["parameters"]]
  truncation   <- prior[["truncation"]]
  lower        <- truncation[["lower"]]
  upper        <- truncation[["upper"]]

  q <- switch(
    distribution,
    "beta" = stats::qbeta(
      probability,
      shape1 = parameters[["alpha"]],
      shape2 = parameters[["beta"]]
    ),
    "normal" = stats::qnorm(
      probability,
      mean = parameters[["mean"]],
      sd   = parameters[["sd"]]
    ),
    Inf
  )

  if (is.finite(lower)) {
    q <- max(q, lower)
  }
  if (is.finite(upper)) {
    q <- min(q, upper)
  }

  return(q)
}

.selection_weightfunction_omega_bound <- function(prior, probability = .999) {

  weights <- prior[["weights"]]
  if (is.null(weights[["type"]])) {
    return(Inf)
  }

  if (identical(weights[["type"]], "fixed")) {
    return(max(weights[["omega"]], 1, na.rm = TRUE))
  }
  if (identical(weights[["type"]], "cumulative")) {
    return(1)
  }
  if (!identical(weights[["type"]], "independent")) {
    return(Inf)
  }

  prior <- weights[["prior"]]
  if (identical(weights[["scale"]], "omega")) {
    upper <- prior[["truncation"]][["upper"]]
    if (is.finite(upper)) {
      return(max(upper, 1))
    }
    return(max(.selection_prior_quantile(prior, probability), 1))
  }
  if (identical(weights[["scale"]], "log_omega")) {
    q <- .selection_prior_quantile(prior, probability)
    if (!is.finite(q) || q >= log(.Machine$double.xmax)) {
      return(Inf)
    }
    return(max(exp(q), 1))
  }

  return(Inf)
}

.selection_telescope_probability_check <- function(priors_bias, z_lower, z_upper) {

  z_threshold     <- 8
  omega_threshold <- 100
  finite_z        <- c(z_lower, z_upper)
  finite_z        <- finite_z[is.finite(finite_z)]
  reasons         <- character()

  if (any(abs(finite_z) > z_threshold)) {
    reasons <- c(reasons, "extreme selection cutpoints")
  }

  step_priors <- priors_bias[vapply(priors_bias, .prior_has_selection, logical(1))]
  omega_bound <- vapply(
    step_priors,
    .selection_weightfunction_omega_bound,
    numeric(1)
  )
  if (any(!is.finite(omega_bound) | omega_bound > omega_threshold)) {
    reasons <- c(reasons, "wide selection-weight priors")
  }

  if (length(reasons) > 0L) {
    warning(
      "Selection-model probability telescoping disabled: ",
      paste(unique(reasons), collapse = " and "),
      ".",
      call. = FALSE
    )
    return(FALSE)
  }

  return(TRUE)
}

.selection_spec <- function(priors, yi, sei, effect_direction, signed_data = FALSE) {

  priors_bias <- .selection_bias_priors(priors)
  if (!any(vapply(priors_bias, .prior_is_selection_kernel, logical(1)))) {
    return(NULL)
  }
  if (any(vapply(priors_bias, .prior_has_phacking, logical(1)))) {
    .selection_stop_phacking_deferred()
  }

  backend <- .selection_backend_spec(priors_bias)
  if (is.null(backend)) {
    return(NULL)
  }

  branch_kernel_mode   <- .selection_branch_kernel_modes(backend)
  jags_use_step_switch <- backend[["mode"]] == "step" &&
    any(branch_kernel_mode == SELKERNEL_NORMAL) &&
    any(branch_kernel_mode == SELKERNEL_STEP)

  p_cuts      <- backend[["step"]][["breaks"]]
  if (is.null(p_cuts)) {
    p_cuts <- c(0, 1)
  }
  p_cuts      <- .selection_assert_p_cuts(p_cuts)
  fixed_omega <- .selection_fixed_omega_by_branch(priors_bias, p_cuts)
  n_bins      <- length(p_cuts) - 1L
  has_step    <- backend[["mode"]] %in% c("step", "step_phack_power")
  has_phack   <- backend[["mode"]] %in% c("phack_power", "step_phack_power")
  phack       <- backend[["phacking"]]
  sign        <- if (signed_data || effect_direction != "negative") 1L else -1L
  z_lower     <- stats::qnorm(1 - p_cuts[-1])
  z_upper     <- stats::qnorm(1 - p_cuts[-length(p_cuts)])
  telescope_probabilities <- has_step &&
    .selection_telescope_probability_check(priors_bias, z_lower, z_upper)
  phack_q_values <- if (has_phack) .selection_phack_q_values(phack[["q"]]) else 1L
  phack_q        <- .selection_phack_static_q(phack_q_values)
  phack_source <- if (has_phack) .selection_phack_static_geometry(phack[["z_source"]], "source") else c(0, 0)
  phack_dest   <- if (has_phack) .selection_phack_static_geometry(phack[["z_destination"]], "destination") else c(0, 0)
  segments     <- .selection_segments(
    p_cuts          = p_cuts,
    phack_z_source  = phack_source,
    phack_z_dest    = phack_dest,
    has_phacking    = has_phack
  )

  jags_omega      <- if (has_step || has_phack) backend[["jags_omega"]] else "sel_omega"
  jags_alpha      <- if (has_phack) backend[["jags_alpha"]] else "sel_phack_alpha"
  jags_phack_kind <- if (has_phack) backend[["jags_phack_kind"]] else "sel_phack_kind"

  jags_data <- list(
    sel_z_lower             = .selection_jags_bounds(z_lower),
    sel_z_upper             = .selection_jags_bounds(z_upper),
    sel_obs_bin             = .selection_obs_bin(yi, sei, p_cuts, sign),
    sel_sign                = sign,
    sel_telescope_probabilities = as.integer(telescope_probabilities)
  )

  if (!identical(backend[["mode"]], "step")) {
    jags_data <- c(jags_data, list(
      sel_kernel_mode          = .selection_kernel_mode(backend[["mode"]]),
      phack_z_source           = phack_source,
      phack_z_dest             = phack_dest,
      sel_segment_bounds       = .selection_jags_bounds(segments[["bounds"]]),
      sel_segment_step_bin     = segments[["step_bin"]],
      sel_segment_phack_region = segments[["phack_region"]]
    ))

    if (!has_step && identical(jags_omega, "sel_omega")) {
      jags_data[[jags_omega]] <- rep(1, n_bins)
    }
    if (!has_phack) {
      jags_data[[jags_alpha]] <- 0
    } else {
      jags_data[["phack_z_source"]] <- NULL
      jags_data[["phack_z_dest"]]   <- NULL
    }
    if (identical(jags_phack_kind, "sel_phack_kind")) {
      jags_data[[jags_phack_kind]] <- if (has_phack) phack_q else 0L
    }
  }
  return(list(
    mode            = backend[["mode"]],
    kernel_mode     = .selection_kernel_mode(backend[["mode"]]),
    p_rule          = "signed_one_sided",
    p_cuts          = p_cuts,
    z_lower         = z_lower,
    z_upper         = z_upper,
    obs_bin         = jags_data[["sel_obs_bin"]],
    sign            = sign,
    n_bins          = n_bins,
    telescope_probabilities = telescope_probabilities,
    has_step        = has_step,
    has_phack       = has_phack,
    phack_q         = phack_q,
    phack_q_values  = phack_q_values,
    mixed_phack_q   = length(phack_q_values) > 1L,
    phack_z_source  = phack_source,
    phack_z_dest    = phack_dest,
    segments        = segments,
    branch_kernel_mode    = branch_kernel_mode,
    fixed_omega           = fixed_omega,
    jags_use_step_switch  = jags_use_step_switch,
    jags_kernel_mode      = if (jags_use_step_switch) "sel_kernel_mode_active" else "sel_kernel_mode",
    jags_kernel_mode_expr = .selection_jags_kernel_mode_expression(branch_kernel_mode),
    jags_omega      = jags_omega,
    jags_alpha      = jags_alpha,
    jags_pi_null    = backend[["jags_pi_null"]],
    jags_beta_null  = backend[["jags_beta_null"]],
    jags_phack_kind = jags_phack_kind,
    jags_phack_z_source = backend[["jags_phack_z_source"]],
    jags_phack_z_dest   = backend[["jags_phack_z_dest"]],
    jags_code       = paste(
      c(backend[["prior_code"]], backend[["transform_code"]]),
      collapse = "\n"
    ),
    backend_data    = backend[["data"]],
    jags_data       = jags_data
  ))
}

.selection_phack_static_q <- function(q) {

  q <- as.integer(q)
  if (length(q) == 0L) {
    return(1L)
  }
  if (length(unique(q)) == 1L) {
    return(q[1L])
  }

  # The active JAGS node phack_kind selects the form for fitted mixtures. The
  # native R wrapper still needs one fallback q for inactive/default calls.
  return(1L)
}

.selection_phack_q_values <- function(q) {

  q <- sort(unique(as.integer(q)))
  q <- q[!is.na(q)]
  if (length(q) == 0L) {
    q <- 1L
  }
  return(q)
}

.selection_phack_static_geometry <- function(z, label) {

  if (is.null(z)) {
    return(c(0, 0))
  }

  if (is.matrix(z)) {
    if (ncol(z) != 2L) {
      stop("P-hacking ", label, " geometry must have two z cut points.", call. = FALSE)
    }
    unique_rows <- unique(z)
    if (nrow(unique_rows) > 1L) {
      stop(
        "RoBMA currently requires p-hacking mixture components to share ",
        label, " geometry.",
        call. = FALSE
      )
    }
    return(as.numeric(unique_rows[1L, ]))
  }

  z <- as.numeric(z)
  if (length(z) != 2L) {
    stop("P-hacking ", label, " geometry must have two z cut points.", call. = FALSE)
  }

  return(z)
}

.selection_reset_native_cache <- function(selection_context) {

  if (!is.null(selection_context)) {
    selection_context[["native_cache"]] <- new.env(parent = emptyenv())
  }

  return(selection_context)
}

.has_native_selnorm_kernel <- local({

  cache   <- FALSE
  symbols <- c(
    "RoBMA_selnorm_kernel_loglik_matrix",
    "RoBMA_selnorm_kernel_log_norm_matrix",
    "RoBMA_selnorm_kernel_cdf_matrix",
    "RoBMA_selnorm_kernel_moments_matrix",
    "RoBMA_selnorm_kernel_rng_matrix",
    "RoBMA_selnorm_kernel_weighted_summary"
  )

  function() {

    if (isTRUE(cache)) {
      return(TRUE)
    }

    cache <<- all(vapply(symbols, is.loaded, logical(1), PACKAGE = "RoBMA"))
    return(cache)
  }
})

.has_native_selnorm_log_norm_delta <- function() {

  return(is.loaded(
    "RoBMA_selnorm_kernel_log_norm_delta_grid",
    PACKAGE = "RoBMA"
  ))
}

.has_native_selnorm_loglik_row_sum <- function() {

  return(is.loaded(
    "RoBMA_selnorm_kernel_loglik_row_sum",
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_loglik_matrix <- function(yi, mu_num, sigma_num,
                                          mu_norm = mu_num,
                                          sigma_norm = sigma_num,
                                          sei, omega, selection_spec,
                                          alpha = NULL,
                                          phack_kind = NULL,
                                          kernel_mode = NULL,
                                          weights = NULL) {

  S <- nrow(mu_num)
  K <- ncol(mu_num)

  if (is.null(weights)) {
    weights <- rep(1, K)
  }
  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]

  if (!.has_native_selnorm_kernel()) {
    stop("The selected-normal native kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_loglik_matrix",
    .native_numeric_vector(yi),
    .native_numeric_matrix(mu_num),
    .native_numeric_matrix(sigma_num),
    .native_numeric_matrix(mu_norm),
    .native_numeric_matrix(sigma_norm),
    .native_numeric_vector(sei),
    .native_numeric_vector(weights),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    .native_integer_vector(selection_spec[["obs_bin"]]),
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_loglik_row_sum <- function(yi, mu_num, sigma_num,
                                           mu_norm = mu_num,
                                           sigma_norm = sigma_num,
                                           sei, omega, selection_spec,
                                           alpha = NULL,
                                           phack_kind = NULL,
                                           kernel_mode = NULL,
                                           weights = NULL) {

  S <- nrow(mu_num)
  K <- ncol(mu_num)

  if (is.null(weights)) {
    weights <- rep(1, K)
  }
  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]

  if (!.has_native_selnorm_loglik_row_sum()) {
    return(rowSums(.selnorm_kernel_loglik_matrix(
      yi             = yi,
      mu_num         = mu_num,
      sigma_num      = sigma_num,
      mu_norm        = mu_norm,
      sigma_norm     = sigma_norm,
      sei            = sei,
      omega          = omega,
      selection_spec = selection_spec,
      alpha          = alpha,
      phack_kind     = phack_kind,
      kernel_mode    = kernel_mode,
      weights        = weights
    )))
  }

  return(.Call(
    "RoBMA_selnorm_kernel_loglik_row_sum",
    .native_numeric_vector(yi),
    .native_numeric_matrix(mu_num),
    .native_numeric_matrix(sigma_num),
    .native_numeric_matrix(mu_norm),
    .native_numeric_matrix(sigma_norm),
    .native_numeric_vector(sei),
    .native_numeric_vector(weights),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    .native_integer_vector(selection_spec[["obs_bin"]]),
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_log_norm_delta_grid <- function(mean, sd, basis,
                                                current_log_norm,
                                                current, values, sei,
                                                weights = NULL, omega,
                                                selection_spec,
                                                alpha = NULL,
                                                phack_kind = NULL,
                                                kernel_mode = NULL) {

  S <- nrow(mean)
  K <- ncol(mean)

  if (is.null(weights)) {
    weights <- rep(1, K)
  }

  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]

  .selection_require_native_no_active_phack(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    S           = S,
    caller      = ".selnorm_kernel_log_norm_delta_grid()"
  )

  if (!.has_native_selnorm_log_norm_delta()) {
    stop("The selected-normal log-normalizer delta kernel is not loaded.",
         call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_log_norm_delta_grid",
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_matrix(basis),
    .native_numeric_matrix(current_log_norm),
    .native_numeric_vector(current),
    .native_numeric_vector(values),
    .native_numeric_vector(sei),
    .native_numeric_vector(weights),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_log_norm_matrix <- function(mean, sd, sei, omega,
                                            selection_spec, alpha = NULL,
                                            phack_kind = NULL,
                                            kernel_mode = NULL) {

  S <- nrow(mean)

  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]

  if (!.has_native_selnorm_kernel()) {
    stop("The selected-normal native kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_log_norm_matrix",
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_cdf_matrix <- function(q, mean, sd, sei, omega,
                                       selection_spec, alpha = NULL,
                                       phack_kind = NULL,
                                       kernel_mode = NULL,
                                       lower.tail = TRUE) {

  S <- nrow(mean)

  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]
  .selection_require_native_no_active_phack(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    S           = S,
    caller      = ".selnorm_kernel_cdf_matrix()"
  )

  if (!.has_native_selnorm_kernel()) {
    stop("The selected-normal native kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_cdf_matrix",
    .native_numeric_vector(q),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    as.logical(lower.tail),
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_moments_matrix <- function(mean, sd, sei, omega,
                                           selection_spec, alpha = NULL,
                                           phack_kind = NULL,
                                           kernel_mode = NULL) {

  S <- nrow(mean)

  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]
  .selection_require_native_no_active_phack(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    S           = S,
    caller      = ".selnorm_kernel_moments_matrix()"
  )

  if (!.has_native_selnorm_kernel()) {
    stop("The selected-normal native kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_moments_matrix",
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_rng_matrix <- function(mean, sd, sei, omega,
                                       selection_spec, alpha = NULL,
                                       phack_kind = NULL,
                                       kernel_mode = NULL) {

  S <- nrow(mean)

  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]
  .selection_require_native_no_active_phack(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    S           = S,
    caller      = ".selnorm_kernel_rng_matrix()"
  )

  if (!.has_native_selnorm_kernel()) {
    stop("The selected-normal native kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_rng_matrix",
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.selnorm_kernel_weighted_summary <- function(yi, mean, sd, sei, psis_weights,
                                             omega, selection_spec,
                                             alpha = NULL, phack_kind = NULL,
                                             kernel_mode = NULL) {

  S <- nrow(mean)

  native_args   <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_spec,
    S              = S,
    alpha          = alpha,
    phack_kind     = phack_kind,
    kernel_mode    = kernel_mode
  )
  native_static <- native_args[["static"]]
  alpha         <- native_args[["alpha"]]
  phack_kind    <- native_args[["phack_kind"]]
  kernel_mode   <- native_args[["kernel_mode"]]
  .selection_require_native_no_active_phack(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    S           = S,
    caller      = ".selnorm_kernel_weighted_summary()"
  )

  if (!.has_native_selnorm_kernel()) {
    stop("The selected-normal native kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_kernel_weighted_summary",
    .native_numeric_vector(yi),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(psis_weights),
    .native_numeric_matrix(omega),
    .native_numeric_vector(alpha),
    .native_integer_vector(phack_kind),
    .native_integer_vector(kernel_mode),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}

.extract_selection_omega_samples <- function(posterior_samples, selection_spec) {

  omega_name <- selection_spec[["jags_omega"]]
  if (is.null(omega_name) || !nzchar(omega_name)) {
    stop("Selection specification is missing 'jags_omega'.", call. = FALSE)
  }

  omega <- BayesTools::JAGS_indexed_parameter_matrix(
    samples   = posterior_samples,
    parameter = omega_name
  )
  if (!is.null(omega)) {
    if (ncol(omega) != selection_spec[["n_bins"]]) {
      stop(
        "Expected ", selection_spec[["n_bins"]], " posterior ",
        omega_name, " column(s), found ", ncol(omega), ".",
        call. = FALSE
      )
    }
    storage.mode(omega) <- "double"
    return(omega)
  }

  omega <- .extract_selection_fixed_omega_samples(
    posterior_samples = posterior_samples,
    selection_spec    = selection_spec
  )
  if (!is.null(omega)) {
    return(omega)
  }

  fixed_omega <- selection_spec[["jags_data"]][[omega_name]]
  if (is.null(fixed_omega) && identical(omega_name, "sel_omega")) {
    fixed_omega <- rep(1, selection_spec[["n_bins"]])
  }
  if (!is.null(fixed_omega)) {
    fixed_omega <- as.numeric(fixed_omega)
    if (length(fixed_omega) != selection_spec[["n_bins"]] ||
        any(!is.finite(fixed_omega))) {
      stop("Invalid fixed selection weights in 'jags_data'.", call. = FALSE)
    }
    omega <- matrix(
      fixed_omega,
      nrow  = nrow(posterior_samples),
      ncol  = selection_spec[["n_bins"]],
      byrow = TRUE
    )
    colnames(omega) <- paste0(omega_name, "[", seq_len(selection_spec[["n_bins"]]), "]")
    storage.mode(omega) <- "double"
    return(omega)
  }

  stop(
    "Missing posterior selection-weight columns for '", omega_name, "'.",
    call. = FALSE
  )
}

.extract_selection_fixed_omega_samples <- function(posterior_samples,
                                                   selection_spec) {

  fixed_omega <- selection_spec[["fixed_omega"]]
  if (is.null(fixed_omega)) {
    return(NULL)
  }

  fixed_omega <- as.matrix(fixed_omega)
  if (ncol(fixed_omega) != selection_spec[["n_bins"]]) {
    return(NULL)
  }

  if (nrow(fixed_omega) == 1L) {
    omega <- matrix(
      fixed_omega[1L, ],
      nrow  = nrow(posterior_samples),
      ncol  = selection_spec[["n_bins"]],
      byrow = TRUE
    )
  } else {
    if (!"bias_indicator" %in% colnames(posterior_samples)) {
      return(NULL)
    }
    indicator <- suppressWarnings(as.integer(round(posterior_samples[, "bias_indicator"])))
    if (any(!is.finite(indicator)) ||
        any(indicator < 1L | indicator > nrow(fixed_omega))) {
      return(NULL)
    }
    omega <- fixed_omega[indicator, , drop = FALSE]
  }

  if (any(!is.finite(omega))) {
    return(NULL)
  }

  colnames(omega) <- paste0(
    selection_spec[["jags_omega"]],
    "[",
    seq_len(selection_spec[["n_bins"]]),
    "]"
  )
  storage.mode(omega) <- "double"

  return(omega)
}

.extract_selection_alpha_samples <- function(posterior_samples, selection_spec) {

  if (!selection_spec[["has_phack"]]) {
    return(rep(0, nrow(posterior_samples)))
  }

  alpha_name <- selection_spec[["jags_alpha"]]
  if (!is.null(alpha_name) && alpha_name %in% colnames(posterior_samples)) {
    return(as.numeric(posterior_samples[, alpha_name]))
  }

  fixed_alpha <- selection_spec[["jags_data"]][[alpha_name]]
  if (!is.null(fixed_alpha)) {
    return(rep(as.numeric(fixed_alpha)[1L], nrow(posterior_samples)))
  }

  stop("Missing posterior p-hacking alpha column for '", alpha_name, "'.",
       call. = FALSE)
}

.extract_selection_phack_kind <- function(posterior_samples, selection_spec) {

  if (!selection_spec[["has_phack"]]) {
    return(rep(0L, nrow(posterior_samples)))
  }
  kind_name <- selection_spec[["jags_phack_kind"]]
  if (!is.null(kind_name) && kind_name %in% colnames(posterior_samples)) {
    return(as.integer(posterior_samples[, kind_name]))
  }
  fixed_kind <- selection_spec[["jags_data"]][[kind_name]]
  if (!is.null(fixed_kind)) {
    return(rep(as.integer(fixed_kind)[1L], nrow(posterior_samples)))
  }
  if (isTRUE(selection_spec[["mixed_phack_q"]])) {
    stop(
      "Missing posterior '", kind_name,
      "' column for mixed linear/quadratic p-hacking forms.",
      call. = FALSE
    )
  }
  return(rep(selection_spec[["phack_q"]], nrow(posterior_samples)))
}

.selection_row_routing <- function(priors, posterior_samples,
                                   selection_spec = NULL) {

  S           <- nrow(posterior_samples)
  priors_bias <- priors[["outcome"]][["bias"]]
  prior_list  <- .selection_bias_priors(priors)

  if ("bias_indicator" %in% colnames(posterior_samples)) {
    bias_indicator <- .extract_posterior_indicator(
      posterior_samples = posterior_samples,
      parameter         = "bias",
      prior             = prior_list,
      column            = "bias_indicator"
    )
  } else if (BayesTools::is.prior.mixture(priors_bias)) {
    stop("Missing posterior model indicator: 'bias_indicator'.",
         call. = FALSE)
  } else {
    bias_indicator <- rep(1L, S)
  }

  if (!.is_priors_weightfunction(priors)) {
    use_normal <- rep(TRUE, S)
  } else {
    selection_indices <- which(vapply(prior_list, .prior_is_selection_kernel, logical(1)))
    use_normal        <- !(bias_indicator %in% selection_indices)
  }

  kernel_mode <- rep(
    if (is.null(selection_spec)) SELKERNEL_NORMAL else selection_spec[["kernel_mode"]],
    S
  )
  kernel_mode[use_normal] <- SELKERNEL_NORMAL

  return(list(
    bias_indicator = bias_indicator,
    use_normal     = use_normal,
    kernel_mode    = kernel_mode
  ))
}

.selection_context_from_parts <- function(fit, data, priors, posterior_samples,
                                           effect_direction,
                                           newdata = NULL) {

  if (!.is_priors_weightfunction(priors)) {
    return(NULL)
  }

  posterior_samples <- .get_posterior_samples(fit, posterior_samples)
  outcome_data      <- if (is.null(newdata)) {
    data[["outcome"]]
  } else if (is.list(newdata) && !is.null(newdata[["outcome"]])) {
    newdata[["outcome"]]
  } else {
    newdata
  }
  yi                <- outcome_data[["yi"]]
  sei               <- outcome_data[["sei"]]

  selection_spec <- .selection_spec(
    priors           = priors,
    yi               = yi,
    sei              = sei,
    effect_direction = effect_direction,
    signed_data      = FALSE
  )

  if (is.null(selection_spec)) {
    return(NULL)
  }

  routing <- .selection_row_routing(
    priors               = priors,
    posterior_samples    = posterior_samples,
    selection_spec       = selection_spec
  )

  selection_context <- selection_spec
  selection_context[["family"]]         <- selection_spec[["mode"]]
  selection_context[["yi"]]             <- yi
  selection_context[["sei"]]            <- sei
  selection_context[["omega"]]          <- .extract_selection_omega_samples(posterior_samples, selection_spec)
  selection_context[["alpha"]]          <- .extract_selection_alpha_samples(posterior_samples, selection_spec)
  selection_context[["phack_kind"]]     <- .extract_selection_phack_kind(posterior_samples, selection_spec)
  selection_context[["kernel_mode"]]    <- routing[["kernel_mode"]]
  selection_context[["bias_indicator"]] <- routing[["bias_indicator"]]
  selection_context[["use_normal"]]     <- routing[["use_normal"]]

  selection_context <- BayesTools::selection_context_validate(
    context           = selection_context,
    n_samples         = nrow(posterior_samples),
    required          = c("omega", "alpha", "phack_kind", "kernel_mode",
                          "bias_indicator", "use_normal")
  )

  return(.selection_reset_native_cache(selection_context))
}

.selection_context <- function(object, posterior_samples = NULL, newdata = NULL) {

  return(.selection_context_from_parts(
    fit                  = object[["fit"]],
    data                 = object[["data"]],
    priors               = object[["priors"]],
    posterior_samples    = posterior_samples,
    effect_direction     = .effect_direction(object),
    newdata              = newdata
  ))
}

.selection_active_phack <- function(selection_context) {

  S           <- nrow(selection_context[["omega"]])
  alpha       <- selection_context[["alpha"]]
  phack_kind  <- selection_context[["phack_kind"]]
  kernel_mode <- selection_context[["kernel_mode"]]

  if (length(alpha) == 1L) {
    alpha <- rep(alpha, S)
  }
  if (length(phack_kind) == 1L) {
    phack_kind <- rep(phack_kind, S)
  }
  if (length(kernel_mode) == 1L) {
    kernel_mode <- rep(kernel_mode, S)
  }

  return(
    kernel_mode %in% c(SELKERNEL_PHACK_POWER, SELKERNEL_STEP_PHACK_POWER) &
      phack_kind > 0L &
      alpha > 0
  )
}

.selection_active_phack_inputs <- function(alpha, phack_kind, kernel_mode, S) {

  alpha <- BayesTools::selection_row_arg(
    x    = alpha,
    n    = S,
    name = "alpha"
  )
  phack_kind <- BayesTools::selection_row_arg(
    x    = phack_kind,
    n    = S,
    name = "phack_kind"
  )
  kernel_mode <- BayesTools::selection_row_arg(
    x    = kernel_mode,
    n    = S,
    name = "kernel_mode"
  )

  return(
    kernel_mode %in% c(SELKERNEL_PHACK_POWER, SELKERNEL_STEP_PHACK_POWER) &
      phack_kind > 0L &
      alpha > 0
  )
}

.selection_require_native_no_active_phack <- function(alpha, phack_kind,
                                                      kernel_mode, S, caller) {

  if (any(.selection_active_phack_inputs(alpha, phack_kind, kernel_mode, S))) {
    stop(
      caller, " does not yet support active p-hacking selection kernels.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}

.selection_require_step_evaluable <- function(selection_context, caller) {

  if (any(.selection_active_phack(selection_context))) {
    stop(
      caller, " does not yet support active p-hacking selection kernels.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}

.selection_interval_prob_vec <- function(lower, upper, mean, sd) {

  out <- stats::pnorm(upper, mean = mean, sd = sd) -
    stats::pnorm(lower, mean = mean, sd = sd)
  out[lower >= upper] <- 0
  out[out < 0 & out > -1e-12] <- 0
  return(out)
}

.selection_step_log_norm_matrix <- function(mean, sd, sei, selection_context) {

  .selection_require_step_evaluable(selection_context, ".selection_step_log_norm_matrix()")

  return(.selnorm_kernel_log_norm_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]]
  ))
}

.selection_step_cdf_matrix <- function(q, mean, sd, sei, selection_context,
                                       lower.tail = TRUE) {

  .selection_require_step_evaluable(selection_context, ".selection_step_cdf_matrix()")

  return(.selnorm_kernel_cdf_matrix(
    q              = q,
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]],
    lower.tail     = lower.tail
  ))
}

.selection_step_logpdf_matrix <- function(y, mean, sd, sei, selection_context) {

  .selection_require_step_evaluable(selection_context, ".selection_step_logpdf_matrix()")

  mean <- as.matrix(mean)
  sd   <- as.matrix(sd)

  if (length(y) == 1L) {
    y <- rep(y, ncol(mean))
  }

  selection_context[["obs_bin"]] <- .selection_obs_bin(
    y,
    sei,
    selection_context[["p_cuts"]],
    selection_context[["sign"]]
  )

  return(.selnorm_kernel_loglik_matrix(
    yi             = y,
    mu_num         = mean,
    sigma_num      = sd,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]]
  ))
}

.selection_step_moments_matrix <- function(mean, sd, sei, selection_context) {

  .selection_require_step_evaluable(selection_context, ".selection_step_moments_matrix()")

  return(.selnorm_kernel_moments_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]]
  ))
}

.selection_step_rng_matrix <- function(mean, sd, sei, selection_context) {

  .selection_require_step_evaluable(selection_context, ".selection_step_rng_matrix()")

  mean <- as.matrix(mean)
  sd   <- as.matrix(sd)

  return(.selnorm_kernel_rng_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]]
  ))
}

.selection_step_weighted_summary <- function(yi, mean, sd, sei, psis_weights,
                                             selection_context) {

  .selection_require_step_evaluable(selection_context, ".selection_step_weighted_summary()")

  mean         <- as.matrix(mean)
  sd           <- as.matrix(sd)
  psis_weights <- as.matrix(psis_weights)

  return(.selnorm_kernel_weighted_summary(
    yi             = yi,
    mean           = mean,
    sd             = sd,
    sei            = sei,
    psis_weights   = psis_weights,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]]
  ))
}
