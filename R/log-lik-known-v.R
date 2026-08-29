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


# ---------------------------------------------------------------------------- #
# .estimate_normal_target_uses_covariance_backend
# ---------------------------------------------------------------------------- #
#
# Normal estimate-deletion targets integrate Gaussian local effects. Known-V
# random-formula models use BayesTools' covariance-factor metadata; the
# specialized multilevel interface uses the same backend with its random
# intercept represented as a row-scaled grouped factor.
#
# ---------------------------------------------------------------------------- #
.estimate_normal_target_uses_covariance_backend <- function(data, priors) {

  if (!identical(.data_outcome_type(data), "norm") || .is_data_weights(data)) {
    return(FALSE)
  }

  if (.is_priors_weightfunction(priors)) {
    return(.is_data_exact_selection(data))
  }

  if (.known_v_estimate_target_uses_backend(data)) {
    return(TRUE)
  }

  .is_data_multilevel(data)
}


.estimate_normal_covariance_target_plan_from_setup <- function(setup) {

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "Gaussian covariance estimate targets are only available for normal ",
      "outcome models.",
      call. = FALSE
    )
  }
  if (isTRUE(setup[["is_weightfunction"]]) &&
      !.is_data_exact_selection(setup[["data"]])) {
    stop(
      "Gaussian covariance estimate targets are not available for selection ",
      "models.",
      call. = FALSE
    )
  }
  if (!is.null(setup[["weights"]])) {
    stop(
      "Gaussian covariance estimate targets are not available for weighted ",
      "likelihoods.",
      call. = FALSE
    )
  }

  K         <- setup[["K"]]
  S         <- setup[["S"]]
  mu_random <- setup[["mu_random"]]
  if (is.null(mu_random)) {
    mu_random <- matrix(0, nrow = S, ncol = K)
  }
  means <- setup[["mu"]] - mu_random

  if (.is_data_known_v(setup[["data"]])) {
    object <- list(
      fit    = setup[["fit"]],
      data   = setup[["data"]],
      priors = setup[["priors"]]
    )
    plan <- .known_v_marginal_factor_plan(
      object            = object,
      posterior_samples = setup[["posterior_samples"]],
      known_V           = .data_known_v_data(setup[["data"]]),
      extra_variances   = .known_v_extra_variance_from_setup(setup)
    )
  } else if (isTRUE(setup[["is_multilevel"]])) {
    group_map <- integer(K)
    for (group in seq_along(setup[["cluster"]])) {
      group_map[setup[["cluster"]][[group]]] <- group
    }
    if (any(group_map == 0L)) {
      stop("Cluster metadata must assign every estimate to one cluster.",
           call. = FALSE)
    }

    factor_plan <- list(
      type                  = "row_group",
      model_matrix          = matrix(1, nrow = K, ncol = 1L),
      group_map             = group_map,
      coefficient_structure = "diagonal"
    )
    factor_states <- lapply(seq_len(S), function(draw) {
      list(list(
        coefficient_factor = matrix(1, nrow = 1L, ncol = 1L),
        row_scale           = as.double(setup[["tau_between"]][draw, ])
      ))
    })
    plan <- list(
      sampling_covariance      = diag(setup[["sei"]]^2, nrow = K, ncol = K),
      random_covariance_plans  = list(factor_plan),
      random_covariance_states = factor_states,
      block_indices            = unname(setup[["cluster"]]),
      extra_variances          = setup[["tau_within"]]^2
    )
  } else if (.is_data_exact_selection(setup[["data"]])) {
    plan <- list(
      sampling_covariance      = diag(setup[["sei"]]^2, nrow = K, ncol = K),
      random_covariance_plans  = list(),
      random_covariance_states = rep(list(list()), S),
      block_indices            = as.list(seq_len(K)),
      extra_variances          = setup[["tau_within"]]^2
    )
  } else {
    stop("Gaussian covariance estimate target metadata are unavailable.",
         call. = FALSE)
  }

  y          <- setup[["yi"]]
  lower_tail <- TRUE
  if (identical(setup[["effect_direction"]], "negative")) {
    y          <- -y
    means      <- -means
    lower_tail <- FALSE
  }

  c(
    list(
      y          = as.double(y),
      means      = means,
      lower_tail = lower_tail
    ),
    plan[c(
      "sampling_covariance",
      "random_covariance_plans",
      "random_covariance_states",
      "block_indices",
      "extra_variances"
    )]
  )
}


# ---------------------------------------------------------------------------- #
# .log_lik_normal_covariance_estimate_target_from_setup
# ---------------------------------------------------------------------------- #
#
# Compute one log-density per estimate and posterior row for a Gaussian
# covariance target. Dependency blocks use Schur-complement conditionals;
# singleton blocks reduce to the usual scalar normal density.
#
# ---------------------------------------------------------------------------- #
.log_lik_normal_covariance_estimate_target_from_setup <- function(
    setup, add_dependency_metadata = FALSE) {

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "Gaussian covariance estimate target is only available ",
      "for normal outcome models.",
      call. = FALSE
    )
  }
  if (isTRUE(setup[["is_weightfunction"]]) &&
      !.is_data_exact_selection(setup[["data"]])) {
    stop(
      "Gaussian covariance estimate target is not available ",
      "for selection models.",
      call. = FALSE
    )
  }
  if (!is.null(setup[["weights"]])) {
    stop(
      "Gaussian covariance estimate target is not available ",
      "for weighted likelihoods.",
      call. = FALSE
    )
  }

  plan <- .estimate_normal_covariance_target_plan_from_setup(setup)

  if (.is_data_exact_selection(setup[["data"]])) {
    conditional <- .marglik_covariance_plan_conditional_summary_batch(
      cache                    = NULL,
      y                        = plan[["y"]],
      means                    = plan[["means"]],
      sampling_covariance      = plan[["sampling_covariance"]],
      random_covariance_plans  = plan[["random_covariance_plans"]],
      random_covariance_states = plan[["random_covariance_states"]],
      block_indices            = plan[["block_indices"]],
      extra_variances          = plan[["extra_variances"]]
    )
    conditional_mean <- matrix(
      plan[["y"]],
      nrow  = setup[["S"]],
      ncol  = setup[["K"]],
      byrow = TRUE
    ) - conditional[["residual"]]
    conditional_sd <- sqrt(conditional[["variance"]])
    selection_context <- .selection_exact_signed_context(
      setup     = setup,
      signed_yi = plan[["y"]]
    )
    log_lik <- .selnorm_kernel_loglik_matrix(
      yi             = plan[["y"]],
      mu_num         = conditional_mean,
      sigma_num      = conditional_sd,
      mu_norm        = conditional_mean,
      sigma_norm     = conditional_sd,
      sei            = setup[["selection_sei"]],
      omega          = selection_context[["omega"]],
      selection_spec = selection_context,
      alpha          = selection_context[["alpha"]],
      phack_kind     = selection_context[["phack_kind"]],
      kernel_mode    = selection_context[["kernel_mode"]]
    )
    if (isTRUE(add_dependency_metadata)) {
      attr(log_lik, "RoBMA_dependency_blocks") <- plan[["block_indices"]]
    }
    return(log_lik)
  }

  log_lik <- .marglik_covariance_plan_conditional_loglik_batch(
    cache                    = NULL,
    y                        = plan[["y"]],
    means                    = plan[["means"]],
    sampling_covariance      = plan[["sampling_covariance"]],
    random_covariance_plans  = plan[["random_covariance_plans"]],
    random_covariance_states = plan[["random_covariance_states"]],
    block_indices            = plan[["block_indices"]],
    extra_variances          = plan[["extra_variances"]]
  )
  if (isTRUE(add_dependency_metadata)) {
    attr(log_lik, "RoBMA_dependency_blocks") <- plan[["block_indices"]]
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
  if (isTRUE(setup[["is_weightfunction"]]) &&
      !.is_data_exact_selection(setup[["data"]])) {
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

  if (.is_data_exact_selection(setup[["data"]])) {
    return(.selection_exact_joint_loglik_from_setup(setup))
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
# .cdf_normal_covariance_estimate_target_from_setup
# ---------------------------------------------------------------------------- #
#
# CDF values matching
# `.log_lik_normal_covariance_estimate_target_from_setup()`.
#
# ---------------------------------------------------------------------------- #
.cdf_normal_covariance_estimate_target_from_setup <- function(setup) {

  summary <- .normal_covariance_estimate_target_summary_from_setup(
    setup      = setup,
    components = "cdf"
  )

  return(summary[["cdf"]])
}


# ---------------------------------------------------------------------------- #
# .normal_covariance_estimate_target_summary_from_setup
# ---------------------------------------------------------------------------- #
#
# Normal CDF and first two moments for each covariance estimate target.
# Dependency blocks use Schur-complement conditionals.
#
# ---------------------------------------------------------------------------- #
.normal_covariance_estimate_target_summary_from_setup <- function(
    setup,
    components = c("cdf", "log_lower", "log_upper", "mean", "variance")) {

  components <- match.arg(
    components,
    c("cdf", "log_lower", "log_upper", "mean", "variance"),
    several.ok = TRUE
  )

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "Gaussian covariance estimate target is only available for normal ",
      "outcome models.",
      call. = FALSE
    )
  }
  if (isTRUE(setup[["is_weightfunction"]]) &&
      !.is_data_exact_selection(setup[["data"]])) {
    stop(
      "Gaussian covariance estimate target is not available for selection ",
      "models.",
      call. = FALSE
    )
  }
  if (!is.null(setup[["weights"]])) {
    stop(
      "Gaussian covariance estimate target is not available for weighted ",
      "likelihoods.",
      call. = FALSE
    )
  }

  plan <- .estimate_normal_covariance_target_plan_from_setup(setup)
  conditional <- .marglik_covariance_plan_conditional_summary_batch(
    cache                    = NULL,
    y                        = plan[["y"]],
    means                    = plan[["means"]],
    sampling_covariance      = plan[["sampling_covariance"]],
    random_covariance_plans  = plan[["random_covariance_plans"]],
    random_covariance_states = plan[["random_covariance_states"]],
    block_indices            = plan[["block_indices"]],
    extra_variances          = plan[["extra_variances"]]
  )
  residual   <- conditional[["residual"]]
  variance   <- conditional[["variance"]]
  lower_tail <- plan[["lower_tail"]]
  yi         <- plan[["y"]]
  S          <- setup[["S"]]
  K          <- setup[["K"]]

  sd  <- sqrt(variance)
  out <- list()

  if (.is_data_exact_selection(setup[["data"]])) {
    conditional_mean <- matrix(
      yi,
      nrow  = S,
      ncol  = K,
      byrow = TRUE
    ) - residual
    selection_context <- .selection_exact_signed_context(
      setup     = setup,
      signed_yi = yi
    )
    if (any(c("cdf", "log_lower", "log_upper") %in% components)) {
      selected_lower <- .selection_step_cdf_matrix(
        q                 = yi,
        mean              = conditional_mean,
        sd                = sd,
        sei               = setup[["selection_sei"]],
        selection_context = selection_context,
        lower.tail        = lower_tail
      )
      selected_upper <- .selection_step_cdf_matrix(
        q                 = yi,
        mean              = conditional_mean,
        sd                = sd,
        sei               = setup[["selection_sei"]],
        selection_context = selection_context,
        lower.tail        = !lower_tail
      )
      if ("cdf" %in% components) {
        out[["cdf"]] <- selected_lower
      }
      if ("log_lower" %in% components) {
        out[["log_lower"]] <- log(selected_lower)
      }
      if ("log_upper" %in% components) {
        out[["log_upper"]] <- log(selected_upper)
      }
    }
    if (any(c("mean", "variance") %in% components)) {
      moments <- .selection_step_moments_matrix(
        mean              = conditional_mean,
        sd                = sd,
        sei               = setup[["selection_sei"]],
        selection_context = selection_context
      )
      if ("mean" %in% components) {
        selected_mean <- moments[["mean"]]
        if (identical(setup[["effect_direction"]], "negative")) {
          selected_mean <- -selected_mean
        }
        out[["mean"]] <- selected_mean
      }
      if ("variance" %in% components) {
        selected_variance <- moments[["second"]] - moments[["mean"]]^2
        if (any(!is.finite(selected_variance)) ||
            any(selected_variance < 0)) {
          stop("Selected-normal predictive variance is invalid.",
               call. = FALSE)
        }
        out[["variance"]] <- selected_variance
      }
    }
    return(out)
  }

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
# Diagonal variance beyond supplied V. The covariance factor plan integrates
# sampled random-effect blocks; already-marginalized row effects contribute here.
#
# ---------------------------------------------------------------------------- #
.known_v_extra_variance_from_setup <- function(setup) {

  data              <- setup[["data"]]
  posterior_samples <- setup[["posterior_samples"]]
  K                 <- setup[["K"]]
  S                 <- setup[["S"]]

  extra_variance <- if (.is_data_exact_selection(data) &&
                        .is_data_random(data)) {
    matrix(0, nrow = S, ncol = K)
  } else if (.is_data_random(data)) {
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
