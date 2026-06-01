# ============================================================================ #
# Normal Cluster-Unit Log-Likelihoods
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_analytic
# ---------------------------------------------------------------------------- #
#
# Analytic cluster-unit likelihood for unweighted normal multilevel models
# without selection. The normal cluster effect is integrated by the block
# covariance.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_analytic <- function(setup, yi, vi) {

  cluster_indices <- setup[["cluster"]]
  S               <- setup[["S"]]
  log_lik         <- matrix(NA_real_, nrow = S, ncol = length(cluster_indices))

  for (g in seq_along(cluster_indices)) {
    idx       <- cluster_indices[[g]]
    diagonal  <- sweep(setup[["tau_within"]][, idx, drop = FALSE]^2, 2L, vi[idx], "+")
    rank_one  <- setup[["tau_between"]][, idx, drop = FALSE]
    residual  <- matrix(yi[idx], nrow = S, ncol = length(idx), byrow = TRUE) -
      setup[["mu"]][, idx, drop = FALSE]
    valid     <- rowSums(!is.finite(diagonal) | diagonal <= 0) == 0L

    if (any(valid)) {
      inv_diag       <- 1 / diagonal[valid, , drop = FALSE]
      rank_valid     <- rank_one[valid, , drop = FALSE]
      residual_valid <- residual[valid, , drop = FALSE]
      denom          <- 1 + rowSums(rank_valid^2 * inv_diag)
      valid_rows     <- which(valid)
      valid_denom    <- is.finite(denom) & denom > .Machine$double.eps

      if (any(valid_denom)) {
        rows      <- valid_rows[valid_denom]
        inv_use   <- inv_diag[valid_denom, , drop = FALSE]
        rank_use  <- rank_valid[valid_denom, , drop = FALSE]
        res_use   <- residual_valid[valid_denom, , drop = FALSE]
        denom_use <- denom[valid_denom]
        log_det   <- rowSums(log(diagonal[rows, , drop = FALSE])) + log(denom_use)
        quad      <- rowSums(res_use^2 * inv_use) -
          rowSums(rank_use * res_use * inv_use)^2 / denom_use

        log_lik[rows, g] <- -0.5 * (length(idx) * log(2 * pi) + log_det + quad)
      }
    }

    fallback <- which(!is.finite(log_lik[, g]))
    for (s in fallback) {
      log_lik[s, g] <- .log_dmvnorm_diag_rank_one(
        x        = yi[idx],
        mean     = setup[["mu"]][s, idx],
        diagonal = vi[idx] + setup[["tau_within"]][s, idx]^2,
        rank_one = setup[["tau_between"]][s, idx]
      )
    }
  }

  return(log_lik)
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm.brma
# ---------------------------------------------------------------------------- #
#
# Analytic cluster-unit likelihood for unweighted normal multilevel models
# without selection. The normal cluster effect is integrated by the block
# covariance.
#
# @param object brma object.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm.brma <- function(object) {

  setup <- .log_lik_cluster_setup.brma(object)
  yi    <- .outcome_data_yi(object)
  vi    <- .outcome_data_vi(object)

  if (setup[["effect_direction"]] == "negative") {
    setup[["mu"]] <- -setup[["mu"]]
    yi            <- -yi
  }

  log_lik <- .log_lik_cluster_norm_analytic(setup = setup, yi = yi, vi = vi)

  return(.add_cluster_log_lik_metadata(log_lik, setup[["cluster"]], setup[["data_hash"]]))
}



# ---------------------------------------------------------------------------- #
# .log_dmvnorm_diag_rank_one
# ---------------------------------------------------------------------------- #
#
# Log-density for N(mean, diag(diagonal) + rank_one %*% t(rank_one)).
#
# ---------------------------------------------------------------------------- #
.log_dmvnorm_diag_rank_one <- function(x, mean, diagonal, rank_one) {

  residual <- x - mean

  if (!all(is.finite(diagonal)) || any(diagonal <= 0)) {
    covariance <- diag(diagonal, nrow = length(diagonal), ncol = length(diagonal)) +
      tcrossprod(rank_one)
    return(mvtnorm::dmvnorm(
      x     = x,
      mean  = mean,
      sigma = covariance,
      log   = TRUE
    ))
  }

  inv_diag <- 1 / diagonal
  inv_rank <- rank_one * inv_diag
  denom    <- 1 + sum(rank_one * inv_rank)

  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    covariance <- diag(diagonal, nrow = length(diagonal), ncol = length(diagonal)) +
      tcrossprod(rank_one)
    return(mvtnorm::dmvnorm(
      x     = x,
      mean  = mean,
      sigma = covariance,
      log   = TRUE
    ))
  }

  log_det <- sum(log(diagonal)) + log(denom)
  quad    <- sum(residual^2 * inv_diag) -
    sum(rank_one * residual * inv_diag)^2 / denom

  return(-0.5 * (length(x) * log(2 * pi) + log_det + quad))
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_quadrature.brma
# ---------------------------------------------------------------------------- #
#
# Gamma-quadrature cluster-unit likelihood for selected-normal models and
# data-weighted normal models. Conditional on gamma, per-estimate likelihood
# contributions factorize.
#
# @param object brma object.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_quadrature.brma <- function(object) {

  setup             <- .log_lik_cluster_setup.brma(object)
  is_weightfunction <- .is_weightfunction(object)
  yi                <- .outcome_data_yi(object)
  sei               <- .outcome_data_sei(object)

  if (setup[["effect_direction"]] == "negative" && !is_weightfunction) {
    setup[["mu"]] <- -setup[["mu"]]
    yi            <- -yi
  }

  if (is_weightfunction) {
    posterior_samples <- setup[["posterior_samples"]]
    selection_context <- .selection_context(
      object            = object,
      posterior_samples = posterior_samples
    )
  } else {
    selection_context <- NULL
  }

  log_lik <- .log_lik_cluster_norm_quadrature(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = is_weightfunction,
    selection_context = selection_context
  )

  return(.add_cluster_log_lik_metadata(log_lik, setup[["cluster"]], setup[["data_hash"]]))
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_quadrature
# ---------------------------------------------------------------------------- #
#
# Gamma-quadrature cluster-unit likelihood for selected-normal models and
# data-weighted normal models. Conditional on gamma, per-estimate likelihood
# contributions factorize.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_quadrature <- function(setup, yi, sei,
                                             is_weightfunction,
                                             selection_context = NULL,
                                             n_gamma = .get_cluster_likelihood_n_gamma()) {

  if (.has_native_norm_cluster_quadrature(selection = is_weightfunction)) {
    return(.log_lik_cluster_norm_quadrature_native(
      setup             = setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = is_weightfunction,
      selection_context = selection_context,
      n_gamma           = n_gamma
    ))
  }

  return(.log_lik_cluster_norm_quadrature_r(
    setup             = setup,
    yi                = yi,
    sei               = sei,
    is_weightfunction = is_weightfunction,
    selection_context = selection_context,
    n_gamma           = n_gamma
  ))
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_quadrature_native
# ---------------------------------------------------------------------------- #
#
# Native batched implementation for normal and selected-normal cluster
# likelihoods.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_quadrature_native <- function(setup, yi, sei,
                                                    is_weightfunction,
                                                    selection_context = NULL,
                                                    n_gamma = .get_cluster_likelihood_n_gamma()) {

  gh      <- .gauss_hermite_nodes(n_gamma)
  cluster <- .cluster_indices_flatten(setup[["cluster"]])
  weights <- if (is.null(setup[["weights"]])) {
    NULL
  } else {
    .native_numeric_vector(setup[["weights"]])
  }

  if (is_weightfunction) {
    if (is.null(selection_context)) {
      stop("Selection context is required for selected-normal cluster likelihood.",
           call. = FALSE)
    }
    native_static <- BayesTools::selection_native_static_args(selection_context)

    return(.Call(
      "RoBMA_selnorm_cluster_loglik",
      .native_numeric_vector(yi),
      .native_numeric_vector(sei),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh[["nodes"]]),
      .native_numeric_vector(gh[["log_weights"]]),
      .native_numeric_matrix(selection_context[["omega"]]),
      .native_numeric_vector(selection_context[["alpha"]]),
      .native_integer_vector(selection_context[["phack_kind"]]),
      .native_integer_vector(selection_context[["kernel_mode"]]),
      native_static[["z_lower"]],
      native_static[["z_upper"]],
      .native_integer_vector(selection_context[["obs_bin"]]),
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

  return(.Call(
    "RoBMA_norm_cluster_loglik",
    .native_numeric_vector(yi),
    .native_numeric_vector(sei),
    .native_numeric_matrix(setup[["mu"]]),
    .native_numeric_matrix(setup[["tau_within"]]),
    .native_numeric_matrix(setup[["tau_between"]]),
    cluster[["index"]],
    cluster[["size"]],
    weights,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    PACKAGE = "RoBMA"
  ))
}


.log_lik_cluster_norm_quadrature_sum <- function(setup, yi, sei,
                                                 is_weightfunction,
                                                 selection_context = NULL,
                                                 n_gamma = .get_cluster_likelihood_n_gamma()) {

  if (!.has_native_norm_loglik_row_sum(
    selection = is_weightfunction,
    cluster   = TRUE
  )) {
    return(rowSums(.log_lik_cluster_norm_quadrature(
      setup             = setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = is_weightfunction,
      selection_context = selection_context,
      n_gamma           = n_gamma
    )))
  }

  gh      <- .gauss_hermite_nodes(n_gamma)
  cluster <- .cluster_indices_flatten(setup[["cluster"]])
  weights <- if (is.null(setup[["weights"]])) {
    NULL
  } else {
    .native_numeric_vector(setup[["weights"]])
  }

  if (is_weightfunction) {
    if (is.null(selection_context)) {
      stop("Selection context is required for selected-normal cluster likelihood.",
           call. = FALSE)
    }
    native_static <- BayesTools::selection_native_static_args(selection_context)

    return(.Call(
      "RoBMA_selnorm_cluster_loglik_row_sum",
      .native_numeric_vector(yi),
      .native_numeric_vector(sei),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh[["nodes"]]),
      .native_numeric_vector(gh[["log_weights"]]),
      .native_numeric_matrix(selection_context[["omega"]]),
      .native_numeric_vector(selection_context[["alpha"]]),
      .native_integer_vector(selection_context[["phack_kind"]]),
      .native_integer_vector(selection_context[["kernel_mode"]]),
      native_static[["z_lower"]],
      native_static[["z_upper"]],
      .native_integer_vector(selection_context[["obs_bin"]]),
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

  return(.Call(
    "RoBMA_norm_cluster_loglik_row_sum",
    .native_numeric_vector(yi),
    .native_numeric_vector(sei),
    .native_numeric_matrix(setup[["mu"]]),
    .native_numeric_matrix(setup[["tau_within"]]),
    .native_numeric_matrix(setup[["tau_between"]]),
    cluster[["index"]],
    cluster[["size"]],
    weights,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    PACKAGE = "RoBMA"
  ))
}


.log_lik_cluster_selnorm_location_grid <- function(setup, yi, sei, basis,
                                                   current, values,
                                                   selection_context,
                                                   n_gamma = .get_cluster_likelihood_n_gamma()) {

  if (!.has_native_selnorm_cluster_location_grid()) {
    return(NULL)
  }
  if (is.null(selection_context)) {
    return(NULL)
  }

  gh      <- .gauss_hermite_nodes(n_gamma)
  cluster <- .cluster_indices_flatten(setup[["cluster"]])
  native_static <- BayesTools::selection_native_static_args(selection_context)
  weights <- if (is.null(setup[["weights"]])) {
    NULL
  } else {
    .native_numeric_vector(setup[["weights"]])
  }

  return(.Call(
    "RoBMA_selnorm_cluster_location_grid",
    .native_numeric_vector(yi),
    .native_numeric_vector(sei),
    .native_numeric_matrix(setup[["mu"]]),
    .native_numeric_matrix(basis),
    .native_numeric_vector(current),
    .native_numeric_vector(values),
    .native_numeric_matrix(setup[["tau_within"]]),
    .native_numeric_matrix(setup[["tau_between"]]),
    cluster[["index"]],
    cluster[["size"]],
    weights,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    .native_integer_vector(selection_context[["obs_bin"]]),
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



# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_quadrature_r
# ---------------------------------------------------------------------------- #
#
# R-composed reference implementation for normal and selected-normal cluster
# likelihoods.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_quadrature_r <- function(setup, yi, sei,
                                               is_weightfunction,
                                               selection_context = NULL,
                                               n_gamma = .get_cluster_likelihood_n_gamma()) {

  gh              <- .gauss_hermite_nodes(n_gamma)
  cluster_indices <- setup[["cluster"]]
  S               <- setup[["S"]]
  G               <- length(cluster_indices)
  log_lik         <- matrix(NA_real_, nrow = S, ncol = G)
  max_cells       <- 2e6

  for (g in seq_along(cluster_indices)) {
    idx                 <- cluster_indices[[g]]
    m                   <- length(idx)
    log_terms           <- matrix(NA_real_, nrow = S, ncol = length(gh[["nodes"]]))
    node_step           <- max(1L, floor(max_cells / max(S * m, 1L)))
    node_step           <- min(node_step, length(gh[["nodes"]]))
    selection_context_g <- if (is_weightfunction) {
      BayesTools::selection_context_subset_observations(
        context = selection_context,
        idx     = idx
      )
    } else {
      NULL
    }

    for (start in seq(1L, length(gh[["nodes"]]), by = node_step)) {
      node_rows <- start:min(start + node_step - 1L, length(gh[["nodes"]]))
      B         <- length(node_rows)
      row_index  <- rep(seq_len(S), times = B)
      node_index <- rep(seq_len(B), each = S)

      mu_node <- setup[["mu"]][row_index, idx, drop = FALSE] +
        gh[["nodes"]][node_rows][node_index] *
        setup[["tau_between"]][row_index, idx, drop = FALSE]

      point_log_lik <- if (is_weightfunction) {
        context_idx <- BayesTools::selection_context_subset_rows(
          context = selection_context_g,
          rows    = row_index
        )

        .outcome_pdf.selnorm(
          yi                = yi[idx],
          mu_samples        = mu_node,
          tau_within        = setup[["tau_within"]][row_index, idx, drop = FALSE],
          sei               = sei[idx],
          selection_context = context_idx
        )
      } else {
        .outcome_pdf.norm(
          yi         = yi[idx],
          mu_samples = mu_node,
          tau_within = setup[["tau_within"]][row_index, idx, drop = FALSE],
          sei        = sei[idx]
        )
      }
      point_log_lik <- .apply_log_lik_weights(point_log_lik, setup[["weights"]][idx])
      log_terms[, node_rows] <- matrix(
        rowSums(point_log_lik),
        nrow = S,
        ncol = B
      ) + matrix(
        gh[["log_weights"]][node_rows],
        nrow  = S,
        ncol  = B,
        byrow = TRUE
      )
    }

    log_lik[, g] <- .rowLogSumExps(log_terms)
  }

  return(log_lik)
}
