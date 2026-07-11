# ============================================================================ #
# Known-V Estimate-Unit Log-Likelihoods
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .known_v_estimate_target_uses_backend
# ---------------------------------------------------------------------------- #
#
# Estimate-wise diagnostics use the known-V backend whenever the target cannot
# be represented by independent scalar likelihood factors after adding latent
# sampling factors to the mean. Only correlated known-V implies Schur
# conditioning on neighboring observed estimates.
#
# ---------------------------------------------------------------------------- #
.known_v_estimate_target_uses_backend <- function(data) {

  .is_data_known_v(data) &&
    (
      .data_known_v_correlated(data) ||
        .is_data_random(data) ||
        .data_has_marginalized_random_effects(data)
    )
}


.known_v_estimate_target_uses_schur_conditioning <- function(data) {

  .is_data_known_v(data) && .data_known_v_correlated(data)
}


# ---------------------------------------------------------------------------- #
# .known_v_dependency_blocks
# ---------------------------------------------------------------------------- #
#
# Dependency blocks used by known-V Schur conditionals and target metadata.
#
# ---------------------------------------------------------------------------- #
.known_v_dependency_blocks <- function(data, K) {

  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    return(as.list(seq_len(K)))
  }

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V block metadata is missing and cannot be reconstructed.",
         call. = FALSE)
  }
  blocks        <- .known_v_blocks(known_V)
  block_indices <- lapply(blocks, `[[`, "index")

  .known_v_validate_dependency_blocks(block_indices, K)
  return(block_indices)
}


.known_v_validate_dependency_blocks <- function(block_indices, K) {

  if (!is.list(block_indices) || length(block_indices) == 0L) {
    stop("Known-V block metadata must be a non-empty list.", call. = FALSE)
  }

  flat <- unlist(block_indices, use.names = FALSE)
  if (!is.numeric(flat) && !is.integer(flat)) {
    stop("Known-V block metadata must contain row indices.", call. = FALSE)
  }
  if (anyNA(flat) || any(!is.finite(flat))) {
    stop("Known-V block metadata must partition the fitted rows.",
         call. = FALSE)
  }
  if (any(flat != as.integer(flat))) {
    stop("Known-V block metadata must contain integer row indices.",
         call. = FALSE)
  }
  flat <- as.integer(flat)
  if (length(flat) != K || !identical(sort(flat), seq_len(K))) {
    stop("Known-V block metadata must partition the fitted rows.",
         call. = FALSE)
  }

  return(invisible(TRUE))
}


# ---------------------------------------------------------------------------- #
# .log_lik_known_v_estimate_target_from_setup
# ---------------------------------------------------------------------------- #
#
# Compute one log-density per estimate and posterior row for the known-V
# estimate target. Correlated dependency blocks use Schur-complement
# conditionals; singleton blocks reduce to the usual scalar normal density.
#
# ---------------------------------------------------------------------------- #
.log_lik_known_v_estimate_target_from_setup <- function(setup) {

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "Known-V estimate-unit log-likelihood target is only available ",
      "for normal outcome models.",
      call. = FALSE
    )
  }
  if (isTRUE(setup[["is_weightfunction"]])) {
    stop(
      "Known-V estimate-unit log-likelihood target is not available ",
      "for selection models.",
      call. = FALSE
    )
  }
  if (!is.null(setup[["weights"]])) {
    stop(
      "Known-V estimate-unit log-likelihood target is not available ",
      "for weighted likelihoods.",
      call. = FALSE
    )
  }

  data       <- setup[["data"]]
  known_V    <- .data_known_v_data(data)
  K          <- setup[["K"]]
  S          <- setup[["S"]]
  yi         <- setup[["yi"]]
  mu_samples <- setup[["mu"]]

  if (identical(setup[["effect_direction"]], "negative")) {
    yi         <- -yi
    mu_samples <- -mu_samples
  }

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance metadata is inconsistent with the outcome data.",
         call. = FALSE)
  }
  block_data        <- .known_v_blocks(known_V)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  .known_v_validate_dependency_blocks(block_indices, K)

  extra_variance <- .known_v_extra_variance_from_setup(setup)

  log_lik <- matrix(NA_real_, nrow = S, ncol = K)
  for (s in seq_len(S)) {
    for (block in seq_along(block_indices)) {
      idx        <- block_indices[[block]]
      covariance <- block_covariances[[block]] +
        diag(extra_variance[s, idx], nrow = length(idx))
      log_lik[s, idx] <- .log_lik_known_v_component_conditional(
        yi         = yi[idx],
        mu         = mu_samples[s, idx],
        covariance = covariance
      )
    }
  }

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .log_lik_known_v_joint_sum_from_setup
# ---------------------------------------------------------------------------- #
#
# Full observed-data known-V log-likelihood for one posterior/evaluated row.
#
# ---------------------------------------------------------------------------- #
.log_lik_known_v_joint_sum_from_setup <- function(setup) {

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "Known-V joint log-likelihood is only available for normal outcome models.",
      call. = FALSE
    )
  }
  if (isTRUE(setup[["is_weightfunction"]])) {
    stop(
      "Known-V joint log-likelihood is not available for selection models.",
      call. = FALSE
    )
  }
  if (!is.null(setup[["weights"]])) {
    stop(
      "Known-V joint log-likelihood is not available for weighted likelihoods.",
      call. = FALSE
    )
  }

  data       <- setup[["data"]]
  known_V    <- .data_known_v_data(data)
  K          <- setup[["K"]]
  S          <- setup[["S"]]
  yi         <- setup[["yi"]]
  mu_samples <- setup[["mu"]]

  if (identical(setup[["effect_direction"]], "negative")) {
    yi         <- -yi
    mu_samples <- -mu_samples
  }

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance metadata is inconsistent with the outcome data.",
         call. = FALSE)
  }

  block_data        <- .known_v_blocks(known_V)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  .known_v_validate_dependency_blocks(block_indices, K)
  extra_variance    <- .known_v_extra_variance_from_setup(setup)
  log_lik           <- numeric(S)

  for (s in seq_len(S)) {
    for (block in seq_along(block_indices)) {
      idx        <- block_indices[[block]]
      covariance <- block_covariances[[block]] +
        diag(extra_variance[s, idx], nrow = length(idx))
      log_lik[[s]] <- log_lik[[s]] + .marglik_mvn_log_density(
        y          = yi[idx],
        mean       = mu_samples[s, idx],
        covariance = covariance,
        context    = "joint likelihood"
      )
    }
  }

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .cdf_known_v_estimate_target_from_setup
# ---------------------------------------------------------------------------- #
#
# CDF values matching `.log_lik_known_v_estimate_target_from_setup()`.
#
# ---------------------------------------------------------------------------- #
.cdf_known_v_estimate_target_from_setup <- function(setup) {

  summary <- .known_v_estimate_target_summary_from_setup(setup)

  return(summary[["cdf"]])
}


# ---------------------------------------------------------------------------- #
# .known_v_estimate_target_summary_from_setup
# ---------------------------------------------------------------------------- #
#
# Normal CDF and first two moments for each known-V estimate target. Correlated
# dependency blocks use Schur-complement conditionals.
#
# ---------------------------------------------------------------------------- #
.known_v_estimate_target_summary_from_setup <- function(setup) {

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "Known-V estimate-unit CDF target is only available for normal ",
      "outcome models.",
      call. = FALSE
    )
  }
  if (isTRUE(setup[["is_weightfunction"]])) {
    stop(
      "Known-V estimate-unit CDF target is not available for selection ",
      "models.",
      call. = FALSE
    )
  }
  if (!is.null(setup[["weights"]])) {
    stop(
      "Known-V estimate-unit CDF target is not available for weighted ",
      "likelihoods.",
      call. = FALSE
    )
  }

  data       <- setup[["data"]]
  known_V    <- .data_known_v_data(data)
  K          <- setup[["K"]]
  S          <- setup[["S"]]
  yi         <- setup[["yi"]]
  mu_samples <- setup[["mu"]]
  lower_tail <- TRUE

  if (identical(setup[["effect_direction"]], "negative")) {
    yi         <- -yi
    mu_samples <- -mu_samples
    lower_tail <- FALSE
  }

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance metadata is inconsistent with the outcome data.",
         call. = FALSE)
  }

  block_data        <- .known_v_blocks(known_V)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  .known_v_validate_dependency_blocks(block_indices, K)
  extra_variance    <- .known_v_extra_variance_from_setup(setup)

  cdf      <- matrix(NA_real_, nrow = S, ncol = K)
  mean     <- matrix(NA_real_, nrow = S, ncol = K)
  variance <- matrix(NA_real_, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    for (block in seq_along(block_indices)) {
      idx        <- block_indices[[block]]
      covariance <- block_covariances[[block]] +
        diag(extra_variance[s, idx], nrow = length(idx))
      distribution <- .known_v_component_conditional_distribution(
        yi         = yi[idx],
        mu         = mu_samples[s, idx],
        covariance = covariance
      )
      cdf[s, idx] <- stats::pnorm(
        distribution[["residual"]],
        mean       = 0,
        sd         = sqrt(distribution[["variance"]]),
        lower.tail = lower_tail
      )
      mean[s, idx]     <- distribution[["mean"]]
      variance[s, idx] <- distribution[["variance"]]
    }
  }

  if (identical(setup[["effect_direction"]], "negative")) {
    mean <- -mean
  }

  second <- variance + mean^2

  return(list(
    cdf      = cdf,
    mean     = mean,
    variance = variance,
    second   = second
  ))
}


# ---------------------------------------------------------------------------- #
# .known_v_extra_variance_from_setup
# ---------------------------------------------------------------------------- #
#
# Diagonal variance beyond supplied V. For random-formula known-V models, sampled
# random effects are conditioned in the mean and only marginalized estimate-level
# random effects contribute here.
#
# ---------------------------------------------------------------------------- #
.known_v_extra_variance_from_setup <- function(setup) {

  data              <- setup[["data"]]
  posterior_samples <- setup[["posterior_samples"]]
  K                 <- setup[["K"]]
  S                 <- setup[["S"]]

  extra_variance <- if (.is_data_random(data)) {
    .evaluate_marginalized_random_variance(
      data              = data,
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = setup[["marginalized_random_source_samples"]]
    )
  } else {
    setup[["tau_within"]]^2
  }
  extra_variance <- as.matrix(extra_variance)

  if (nrow(extra_variance) != S || ncol(extra_variance) != K) {
    stop("Known-V diagonal variance contributions have inconsistent dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(extra_variance)) || any(extra_variance < 0)) {
    stop("Known-V diagonal variance contributions must be non-negative.",
         call. = FALSE)
  }

  return(extra_variance)
}


.known_v_marginalized_random_source_samples <- function(fit, data, priors,
                                                        posterior_samples) {

  if (!.is_data_known_v(data) ||
      !.data_has_marginalized_random_effects(data) ||
      !.is_data_scale(data)) {
    return(NULL)
  }

  object <- list(
    fit    = fit,
    data   = data,
    priors = priors
  )

  .predict_known_v_newdata_marginalized_source_samples(
    object            = object,
    data              = data,
    posterior_samples = posterior_samples
  )
}


.known_v_marginalized_random_source_samples_from_tau <- function(
    data, tau_within_samples) {

  if (!.is_data_known_v(data) ||
      !.data_has_marginalized_random_effects(data) ||
      !.is_data_scale(data)) {
    return(NULL)
  }

  source_names <- .known_v_marginalized_random_row_source_names(data)
  if (length(source_names) != 1L) {
    return(NULL)
  }

  stats::setNames(
    list(as.matrix(tau_within_samples)),
    source_names
  )
}


.known_v_marginalized_random_row_source_names <- function(data) {

  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return(character())
  }

  out <- character()
  for (term in terms) {
    sources <- .predict_known_v_marginalized_sd_sources(term)
    for (source in sources) {
      if (identical(source[["shape"]], "row")) {
        out <- c(out, source[["name"]])
      }
    }
  }

  unique(out[!is.na(out) & nzchar(out)])
}


.known_v_extra_sd_from_setup <- function(setup) {

  sqrt(pmax(.known_v_extra_variance_from_setup(setup), 0))
}


# ---------------------------------------------------------------------------- #
# .known_v_estimate_blup_from_setup
# ---------------------------------------------------------------------------- #
#
# Posterior means of estimate-level true effects for known-V diagnostics.
# Sampled random effects are already included in setup[["mu"]]; diagonal
# marginalized heterogeneity enters through `.known_v_extra_variance_from_setup()`.
#
# ---------------------------------------------------------------------------- #
.known_v_estimate_blup_from_setup <- function(setup) {

  extra_variance <- .known_v_extra_variance_from_setup(setup)

  .evaluate.brma.known_v_blup.norm(
    mu_samples = setup[["mu"]],
    tau_within = sqrt(pmax(extra_variance, 0)),
    yi         = setup[["yi"]],
    known_V    = .data_known_v_data(setup[["data"]])
  )
}


# ---------------------------------------------------------------------------- #
# .log_lik_known_v_component_conditional
# ---------------------------------------------------------------------------- #
#
# Conditional normal log-densities for all coordinates in one covariance block.
# For a joint normal distribution, Var(y_i | y_-i) is 1 / precision[i, i] and
# E(y_i | y_-i) shifts by the precision-weighted joint residual.
#
# ---------------------------------------------------------------------------- #
.log_lik_known_v_component_conditional <- function(yi, mu, covariance) {

  distribution <- .known_v_component_conditional_distribution(
    yi         = yi,
    mu         = mu,
    covariance = covariance
  )

  -0.5 * (
    log(2 * pi * distribution[["variance"]]) +
      distribution[["residual"]]^2 / distribution[["variance"]]
  )
}


# ---------------------------------------------------------------------------- #
# .known_v_component_conditional_distribution
# ---------------------------------------------------------------------------- #
#
# Conditional normal distribution for all coordinates in one block.
#
# ---------------------------------------------------------------------------- #
.known_v_component_conditional_distribution <- function(yi, mu, covariance) {

  size <- length(yi)
  if (size == 1L) {
    return(list(
      mean     = mu,
      variance = covariance[1L, 1L],
      residual = yi - mu
    ))
  }

  chol_covariance <- .known_v_chol_covariance(
    covariance = covariance,
    context    = "conditional"
  )

  precision          <- chol2inv(chol_covariance)
  precision_diagonal <- diag(precision)
  if (any(!is.finite(precision_diagonal)) || any(precision_diagonal <= 0)) {
    stop(
      "Known-V conditional precision has non-positive diagonal entries.",
      call. = FALSE
    )
  }

  residual <- yi - mu
  conditional_residual <- as.vector(precision %*% residual) /
    precision_diagonal
  conditional_variance <- 1 / precision_diagonal
  conditional_mean     <- yi - conditional_residual

  return(list(
    mean     = conditional_mean,
    variance = conditional_variance,
    residual = conditional_residual
  ))
}


.known_v_chol_covariance <- function(covariance, context) {

  factorization <- .covariance_factorization(covariance)
  chol_covariance <- .covariance_cholesky(factorization)
  if (is.null(chol_covariance)) {
    .known_v_stop_non_positive_definite_covariance(covariance, context)
  }

  chol_covariance
}


.known_v_stop_non_positive_definite_covariance <- function(covariance, context) {

  factorization <- tryCatch(
    .covariance_factorization(covariance),
    error = function(e) NULL
  )

  if (!is.null(factorization)) {
    if (.covariance_is_positive_semidefinite(factorization)) {
      stop(
        "Known-V ", context, " covariance is positive semidefinite, ",
        "not positive definite; this Cholesky-based target is not ",
        "available for singular known-V blocks unless posterior extra ",
        "variance makes the block positive definite.",
        call. = FALSE
      )
    }
  }

  stop(
    "Known-V ", context, " covariance is not positive definite.",
    call. = FALSE
  )
}
