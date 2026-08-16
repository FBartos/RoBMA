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

  block_data <- .known_v_dependency_block_data(data, K)
  return(lapply(block_data, `[[`, "index"))
}


# Return validated row indices and covariances for each known-V block.
.known_v_dependency_block_data <- function(data, K) {

  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    stop("Known-V covariance metadata is unavailable.", call. = FALSE)
  }

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V block metadata is missing and cannot be reconstructed.",
         call. = FALSE)
  }
  block_data    <- .known_v_blocks(known_V)
  block_indices <- lapply(block_data, `[[`, "index")

  .known_v_validate_dependency_blocks(block_indices, K)
  return(block_data)
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


.known_v_plan_block_data <- function(block_data, selected_blocks) {

  selected_data    <- block_data[selected_blocks]
  selected_indices <- lapply(selected_data, `[[`, "index")
  global_indices   <- as.integer(unlist(selected_indices, use.names = FALSE))
  block_sizes      <- lengths(selected_indices)
  local_blocks     <- vector("list", length(block_sizes))
  covariance       <- matrix(
    0,
    nrow = length(global_indices),
    ncol = length(global_indices)
  )
  rank_one_factors <- lapply(selected_data, function(block) {
    .covariance_exact_rank_one_factor(block[["covariance"]])
  })
  is_rank_one      <- !vapply(rank_one_factors, is.null, logical(1L))
  if (any(is_rank_one) && !all(is_rank_one)) {
    stop("Known-V covariance plan blocks mix incompatible representations.",
         call. = FALSE)
  }

  offset <- 0L
  for (block in seq_along(selected_data)) {
    local_index <- offset + seq_len(block_sizes[[block]])
    local_blocks[[block]] <- local_index
    if (!is_rank_one[[block]]) {
      covariance[local_index, local_index] <-
        selected_data[[block]][["covariance"]]
    }
    offset <- offset + block_sizes[[block]]
  }

  if (all(is_rank_one)) {
    # Treat an exact rank-one sampling block as a fixed low-rank factor. This
    # is an algebraic representation of supplied V, not an added random term;
    # it avoids losing row-specific sub-ULP diagonal variance in V + diag(d).
    factor_plan  <- list(
      type                  = "group",
      model_matrix          = matrix(unlist(
        rank_one_factors,
        use.names = FALSE
      ), ncol = 1L),
      group_map             = rep.int(seq_along(block_sizes), block_sizes),
      coefficient_structure = "diagonal"
    )
    factor_state  <- list(coefficient_factor = matrix(1, 1L, 1L))
    factor_plans  <- list(factor_plan)
    factor_states <- list(factor_state)
  } else {
    factor_plans  <- list()
    factor_states <- list()
  }

  return(list(
    global_indices  = global_indices,
    block_indices   = local_blocks,
    covariance      = covariance,
    factor_plans    = factor_plans,
    factor_states   = factor_states
  ))
}


.known_v_plan_block_groups <- function(block_data, selected_blocks) {

  is_rank_one <- vapply(selected_blocks, function(block) {
    !is.null(.covariance_exact_rank_one_factor(
      block_data[[block]][["covariance"]]
    ))
  }, logical(1L))
  groups <- list(
    selected_blocks[!is_rank_one],
    selected_blocks[is_rank_one]
  )

  return(groups[lengths(groups) > 0L])
}


.known_v_rethrow_conditional_error <- function(
    error, block_data, selected_blocks, extra_variance) {

  for (s in seq_len(nrow(extra_variance))) {
    for (block in block_data[selected_blocks]) {
      idx <- block[["index"]]
      covariance <- block[["covariance"]] +
        diag(extra_variance[s, idx], nrow = length(idx))
      .known_v_chol_covariance(
        covariance = covariance,
        context    = "conditional"
      )
    }
  }
  stop(conditionMessage(error), call. = FALSE)
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
  block_data     <- .known_v_dependency_block_data(data, K)
  block_sizes    <- lengths(lapply(block_data, `[[`, "index"))
  singleton      <- which(block_sizes == 1L)
  dependent      <- which(block_sizes > 1L)
  extra_variance <- .known_v_extra_variance_from_setup(setup)
  log_lik        <- matrix(NA_real_, nrow = S, ncol = K)

  if (length(singleton) > 0L) {
    singleton_indices <- vapply(
      singleton,
      function(block) block_data[[block]][["index"]][[1L]],
      integer(1L)
    )
    sampling_variance <- vapply(
      singleton,
      function(block) block_data[[block]][["covariance"]][1L, 1L],
      numeric(1L)
    )
    variance <- sweep(
      extra_variance[, singleton_indices, drop = FALSE],
      MARGIN = 2L,
      STATS  = sampling_variance,
      FUN    = "+"
    )
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      invalid <- which(
        !is.finite(variance) | variance <= 0,
        arr.ind = TRUE
      )[1L, ]
      .known_v_chol_covariance(
        covariance = matrix(variance[invalid[[1L]], invalid[[2L]]], 1L, 1L),
        context    = "conditional"
      )
    }
    residual <- sweep(
      mu_samples[, singleton_indices, drop = FALSE],
      MARGIN = 2L,
      STATS  = yi[singleton_indices],
      FUN    = "-"
    )
    log_lik[, singleton_indices] <- -0.5 * (
      log(2 * pi * variance) + residual^2 / variance
    )
  }

  if (length(dependent) > 0L) {
    for (blocks in .known_v_plan_block_groups(block_data, dependent)) {
      plan_data <- .known_v_plan_block_data(block_data, blocks)
      idx       <- plan_data[["global_indices"]]
      states    <- rep(list(plan_data[["factor_states"]]), S)
      log_lik[, idx] <- tryCatch(
        .marglik_covariance_plan_conditional_loglik_batch(
          cache                    = NULL,
          y                        = as.double(yi[idx]),
          means                    = mu_samples[, idx, drop = FALSE],
          sampling_covariance      = plan_data[["covariance"]],
          random_covariance_plans  = plan_data[["factor_plans"]],
          random_covariance_states = states,
          block_indices            = plan_data[["block_indices"]],
          extra_variances          = extra_variance[, idx, drop = FALSE]
        ),
        error = function(error) {
          .known_v_rethrow_conditional_error(
            error           = error,
            block_data      = block_data,
            selected_blocks = blocks,
            extra_variance  = extra_variance
          )
        }
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

  block_data     <- .known_v_dependency_block_data(data, K)
  block_sizes    <- lengths(lapply(block_data, `[[`, "index"))
  singleton      <- which(block_sizes == 1L)
  dependent      <- which(block_sizes > 1L)
  extra_variance <- .known_v_extra_variance_from_setup(setup)
  log_lik        <- numeric(S)

  if (length(singleton) > 0L) {
    singleton_indices <- vapply(
      singleton,
      function(block) block_data[[block]][["index"]][[1L]],
      integer(1L)
    )
    sampling_variance <- vapply(
      singleton,
      function(block) block_data[[block]][["covariance"]][1L, 1L],
      numeric(1L)
    )
    variance <- sweep(
      extra_variance[, singleton_indices, drop = FALSE],
      MARGIN = 2L,
      STATS  = sampling_variance,
      FUN    = "+"
    )
    residual <- sweep(
      mu_samples[, singleton_indices, drop = FALSE],
      MARGIN = 2L,
      STATS  = yi[singleton_indices],
      FUN    = "-"
    )
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      invalid <- which(
        !is.finite(variance) | variance <= 0,
        arr.ind = TRUE
      )[1L, ]
      .known_v_chol_covariance(
        covariance = matrix(variance[invalid[[1L]], invalid[[2L]]], 1L, 1L),
        context    = "joint likelihood"
      )
    }
    log_lik <- rowSums(-0.5 * (
      log(2 * pi * variance) + residual^2 / variance
    ))
  }

  if (length(dependent) > 0L) {
    for (blocks in .known_v_plan_block_groups(block_data, dependent)) {
      plan_data <- .known_v_plan_block_data(block_data, blocks)
      idx       <- plan_data[["global_indices"]]
      states    <- rep(list(plan_data[["factor_states"]]), S)
      log_lik <- log_lik + tryCatch(
        .marglik_covariance_plan_loglik_batch(
          cache                    = NULL,
          y                        = as.double(yi[idx]),
          means                    = mu_samples[, idx, drop = FALSE],
          sampling_covariance      = plan_data[["covariance"]],
          random_covariance_plans  = plan_data[["factor_plans"]],
          random_covariance_states = states,
          block_indices            = plan_data[["block_indices"]],
          extra_variances          = extra_variance[, idx, drop = FALSE]
        ),
        error = function(error) {
          for (s in seq_len(S)) {
            for (block in block_data[blocks]) {
              block_idx <- block[["index"]]
              covariance <- block[["covariance"]] +
                diag(extra_variance[s, block_idx], nrow = length(block_idx))
              .known_v_chol_covariance(
                covariance = covariance,
                context    = "joint likelihood"
              )
            }
          }
          stop(conditionMessage(error), call. = FALSE)
        }
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

  summary <- .known_v_estimate_target_summary_from_setup(
    setup      = setup,
    components = "cdf"
  )

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
.known_v_estimate_target_summary_from_setup <- function(
    setup,
    components = c("cdf", "log_lower", "log_upper", "mean", "variance")) {

  components <- match.arg(
    components,
    c("cdf", "log_lower", "log_upper", "mean", "variance"),
    several.ok = TRUE
  )

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

  block_data     <- .known_v_dependency_block_data(data, K)
  block_sizes    <- lengths(lapply(block_data, `[[`, "index"))
  singleton      <- which(block_sizes == 1L)
  dependent      <- which(block_sizes > 1L)
  extra_variance <- .known_v_extra_variance_from_setup(setup)
  residual       <- matrix(NA_real_, nrow = S, ncol = K)
  variance       <- matrix(NA_real_, nrow = S, ncol = K)

  if (length(singleton) > 0L) {
    singleton_indices <- vapply(
      singleton,
      function(block) block_data[[block]][["index"]][[1L]],
      integer(1L)
    )
    sampling_variance <- vapply(
      singleton,
      function(block) block_data[[block]][["covariance"]][1L, 1L],
      numeric(1L)
    )
    residual[, singleton_indices] <- -sweep(
      mu_samples[, singleton_indices, drop = FALSE],
      MARGIN = 2L,
      STATS  = yi[singleton_indices],
      FUN    = "-"
    )
    singleton_variance <- sweep(
      extra_variance[, singleton_indices, drop = FALSE],
      MARGIN = 2L,
      STATS  = sampling_variance,
      FUN    = "+"
    )
    if (any(!is.finite(singleton_variance)) ||
        any(singleton_variance <= 0)) {
      invalid <- which(
        !is.finite(singleton_variance) | singleton_variance <= 0,
        arr.ind = TRUE
      )[1L, ]
      .known_v_chol_covariance(
        covariance = matrix(
          singleton_variance[invalid[[1L]], invalid[[2L]]],
          1L,
          1L
        ),
        context    = "conditional"
      )
    }
    variance[, singleton_indices] <- singleton_variance
  }

  if (length(dependent) > 0L) {
    for (blocks in .known_v_plan_block_groups(block_data, dependent)) {
      plan_data   <- .known_v_plan_block_data(block_data, blocks)
      idx         <- plan_data[["global_indices"]]
      states      <- rep(list(plan_data[["factor_states"]]), S)
      conditional <- tryCatch(
        .marglik_covariance_plan_conditional_summary_batch(
          cache                    = NULL,
          y                        = as.double(yi[idx]),
          means                    = mu_samples[, idx, drop = FALSE],
          sampling_covariance      = plan_data[["covariance"]],
          random_covariance_plans  = plan_data[["factor_plans"]],
          random_covariance_states = states,
          block_indices            = plan_data[["block_indices"]],
          extra_variances          = extra_variance[, idx, drop = FALSE]
        ),
        error = function(error) {
          .known_v_rethrow_conditional_error(
            error           = error,
            block_data      = block_data,
            selected_blocks = blocks,
            extra_variance  = extra_variance
          )
        }
      )
      residual[, idx] <- conditional[["residual"]]
      variance[, idx] <- conditional[["variance"]]
    }
  }

  sd  <- sqrt(variance)
  out <- list()

  if ("cdf" %in% components) {
    out[["cdf"]] <- stats::pnorm(
      residual,
      sd         = sd,
      lower.tail = lower_tail
    )
  }
  if ("log_lower" %in% components) {
    out[["log_lower"]] <- stats::pnorm(
      residual,
      sd         = sd,
      lower.tail = lower_tail,
      log.p      = TRUE
    )
  }
  if ("log_upper" %in% components) {
    out[["log_upper"]] <- stats::pnorm(
      residual,
      sd         = sd,
      lower.tail = !lower_tail,
      log.p      = TRUE
    )
  }
  if ("mean" %in% components) {
    mean <- matrix(yi, nrow = S, ncol = K, byrow = TRUE) - residual
    if (identical(setup[["effect_direction"]], "negative")) {
      mean <- -mean
    }
    out[["mean"]] <- mean
  }
  if ("variance" %in% components) {
    out[["variance"]] <- variance
  }

  return(out)
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

  sqrt(.known_v_extra_variance_from_setup(setup))
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
    tau_within = sqrt(extra_variance),
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


# Conditional moments for diag(diagonal) + rank_one rank_one'.
.known_v_diagonal_rank_one_conditional <- function(yi, mu, diagonal,
                                                   rank_one) {

  size <- length(yi)
  if (length(mu) != size || length(diagonal) != size ||
      length(rank_one) != size || anyNA(diagonal) ||
      any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop(
      "Known-V conditional covariance is positive semidefinite, not positive ",
      "definite; positive diagonal extra variance is required.",
      call. = FALSE
    )
  }

  residual             <- yi - mu
  conditional_residual <- numeric(size)
  conditional_variance <- numeric(size)
  for (i in seq_len(size)) {
    other       <- seq_len(size) != i
    denominator <- 1 + sum(rank_one[other]^2 / diagonal[other])
    adjustment  <- rank_one[[i]] *
      sum(rank_one[other] * residual[other] / diagonal[other]) /
      denominator
    conditional_residual[[i]] <- residual[[i]] - adjustment
    conditional_variance[[i]] <- diagonal[[i]] +
      rank_one[[i]]^2 / denominator
  }

  if (any(!is.finite(conditional_residual)) ||
      any(!is.finite(conditional_variance)) ||
      any(conditional_variance <= 0)) {
    stop("Known-V conditional covariance could not be evaluated faithfully.",
         call. = FALSE)
  }

  return(list(
    mean     = yi - conditional_residual,
    variance = conditional_variance,
    residual = conditional_residual
  ))
}


# Joint density without materializing a diagonal-plus-rank-one covariance.
.known_v_diagonal_rank_one_log_density <- function(y, mean, diagonal, rank_one,
                                                   context) {

  if (anyNA(diagonal) || any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop(
      "Known-V ", context, " covariance is positive semidefinite, not ",
      "positive definite; positive diagonal extra variance is required.",
      call. = FALSE
    )
  }

  .log_dmvnorm_diag_rank_one(
    x        = y,
    mean     = mean,
    diagonal = diagonal,
    rank_one = rank_one
  )
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
