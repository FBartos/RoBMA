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

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V block metadata is missing and cannot be reconstructed.",
         call. = FALSE)
  }
  if (identical(.known_v_storage(known_V), "factor")) {
    block_indices <- known_V[["block_indices"]]
    .known_v_validate_dependency_blocks(block_indices, K)
    return(block_indices)
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
    return(.is_data_joint_selection(data))
  }

  if (.known_v_estimate_target_uses_backend(data)) {
    return(TRUE)
  }

  .is_data_multilevel(data)
}


.estimate_normal_covariance_target_location_from_setup <- function(setup) {

  mu_random <- setup[["mu_random"]]
  if (is.null(mu_random)) {
    mu_random <- matrix(0, nrow = setup[["S"]], ncol = setup[["K"]])
  }
  # Selection setup includes only retained context. Its integrated sources
  # already belong to the conditional covariance and must not remove H here.
  means <- if (.is_data_joint_selection(setup[["data"]])) {
    setup[["mu"]]
  } else {
    setup[["mu"]] - mu_random
  }
  y     <- setup[["yi"]]

  if (identical(setup[["effect_direction"]], "negative")) {
    y     <- -y
    means <- -means
  }

  list(
    y          = as.double(y),
    means      = means,
    lower_tail = !identical(setup[["effect_direction"]], "negative")
  )
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
      !.is_data_joint_selection(setup[["data"]])) {
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

  K        <- setup[["K"]]
  S        <- setup[["S"]]
  location <- .estimate_normal_covariance_target_location_from_setup(setup)

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
  } else if (.setup_uses_joint_selection_likelihood(setup)) {
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

  c(
    location,
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
.selection_joint_conditional_summary_from_setup <- function(setup) {

  location <- .estimate_normal_covariance_target_location_from_setup(setup)
  data     <- setup[["data"]]
  plan     <- .data_selection_execution_plan(data)
  S        <- setup[["S"]]
  K        <- setup[["K"]]
  random_factors <- list(
    factor_plans  = list(),
    factor_states = rep(list(list()), S)
  )
  extra_variances <- if (.selection_integrates_estimate(data)) {
    setup[["tau_within"]]^2
  } else {
    matrix(0, S, K)
  }
  if (.is_data_random(data)) {
    extra_variances <- matrix(0, S, K)
    blocks <- plan[["random_covariance"]][["term_names"]]
    if (length(blocks) > 0L) {
      # Only the resolved integrated sources enter covariance. In particular,
      # NULL would request every term, including retained contextual effects.
      random_factors <- .brma_mv_random_effects_marginal_factor_plan(
        object            = list(
          fit = setup[["fit"]], data = data, priors = setup[["priors"]]
        ),
        posterior_samples = setup[["posterior_samples"]],
        blocks            = blocks,
        row_blocks        = plan[["row_blocks"]]
      )
    }
  } else if (isTRUE(setup[["is_multilevel"]]) &&
             !.selection_retains_other_random(data)) {
    group_map <- integer(K)
    for (group in seq_along(setup[["cluster"]])) {
      group_map[setup[["cluster"]][[group]]] <- group
    }
    if (any(group_map == 0L)) {
      stop("Cluster metadata must assign every estimate to one cluster.",
           call. = FALSE)
    }
    random_factors[["factor_plans"]] <- list(list(
      type                  = "row_group",
      model_matrix          = matrix(1, K, 1L),
      group_map             = group_map,
      coefficient_structure = "diagonal"
    ))
    random_factors[["factor_states"]] <- lapply(seq_len(S), function(draw) {
      list(list(
        coefficient_factor = matrix(1, 1L, 1L),
        row_scale           = as.double(setup[["tau_between"]][draw, ])
      ))
    })
  }

  # Keep the native Markov accelerator for interior states. At a valid
  # zero-innovation boundary, the same authoritative root remains usable by
  # the general native factor route. Invalid metadata must still fail.
  for (i in seq_along(random_factors[["factor_plans"]])) {
    if (identical(random_factors[["factor_plans"]][[i]][[
        "coefficient_structure"
      ]], "markov")) {
      n_columns <- ncol(random_factors[["factor_plans"]][[i]][["model_matrix"]])
      boundary <- vapply(random_factors[["factor_states"]], function(state) {
        markov <- .marglik_validate_random_covariance_markov_state(
          state[[i]], n_columns, allow_zero_innovation = TRUE
        )
        any(markov[["markov_innovation_variance"]] == 0)
      }, logical(1L))
      if (any(boundary)) {
        random_factors[["factor_plans"]][[i]][["coefficient_structure"]] <- "dense"
      }
    }
  }
  native_supported <- all(vapply(
    random_factors[["factor_plans"]],
    function(factor) {
      identical(factor[["type"]], "dense") ||
        (isTRUE(factor[["type"]] %in% c("group", "row_group", "known_group")) &&
         isTRUE(factor[["coefficient_structure"]] %in% c("diagonal", "dense", "markov")))
    },
    logical(1L)
  ))

  if (native_supported) {
    sampling <- plan[["sampling"]]
    if (identical(sampling[["representation"]], "dense")) {
      sampling_covariance <- sampling[["covariance"]]
    } else if (identical(sampling[["representation"]], "diagonal_factor")) {
      # Blocks a recovered factor does not represent keep their supplied
      # entries in the base; the loading adds only the recovered blocks.
      sampling_covariance <- .selection_joint_sampling_dense_base(sampling)
      loading <- sampling[["loading"]]
      rank    <- ncol(loading)
      if (rank > 0L) {
        random_factors[["factor_plans"]] <- c(
          random_factors[["factor_plans"]],
          list(list(
            type                  = "group",
            model_matrix          = loading,
            group_map             = rep(1L, K),
            coefficient_structure = "diagonal"
          ))
        )
        sampling_state <- list(coefficient_factor = diag(1, rank, rank))
        random_factors[["factor_states"]] <- lapply(
          random_factors[["factor_states"]],
          function(state) c(state, list(sampling_state))
        )
      }
    } else {
      stop("Selection sampling-covariance metadata are invalid.",
           call. = FALSE)
    }
    summary <- .marglik_covariance_plan_conditional_summary_batch(
      cache                    = NULL,
      y                        = location[["y"]],
      means                    = location[["means"]],
      sampling_covariance      = sampling_covariance,
      random_covariance_plans  = random_factors[["factor_plans"]],
      random_covariance_states = random_factors[["factor_states"]],
      block_indices            = plan[["row_blocks"]],
      extra_variances          = extra_variances
    )
    variances <- summary[["variance"]]
    means <- matrix(location[["y"]], S, K, byrow = TRUE) - summary[["residual"]]
  } else {
    # Unsupported factor contracts retain the existing general covariance
    # evaluator. Numerical errors in an admitted native route are not retried.
    factors <- .selection_joint_random_factor_samples(setup)
    covariance <- if (is.null(factors)) {
      .selection_joint_random_covariance_samples(setup)
    } else {
      NULL
    }
    means <- variances <- matrix(NA_real_, S, K)
    for (block in seq_along(plan[["row_blocks"]])) {
      index <- plan[["row_blocks"]][[block]]
      k <- length(index)
      lower <- .selection_joint_covariance_lower(
        setup                     = setup,
        block_index               = block,
        random_covariance_samples = covariance,
        random_factor_samples     = factors
      )
      pairs <- .selection_joint_lower_pairs(plan, seq_len(k))
      for (draw in seq_len(S)) {
        sigma <- matrix(0, k, k)
        sigma[cbind(pairs[["row_1"]], pairs[["row_2"]])] <- lower[draw, ]
        sigma[cbind(pairs[["row_2"]], pairs[["row_1"]])] <- lower[draw, ]
        factor <- tryCatch(chol(sigma), error = function(e) NULL)
        if (is.null(factor)) {
          stop("Selection deletion covariance must be positive definite.",
               call. = FALSE)
        }
        precision <- chol2inv(factor)
        variance <- 1 / diag(precision)
        residual <- as.vector(precision %*%
          (location[["y"]][index] - location[["means"]][draw, index]))
        variances[draw, index] <- variance
        means[draw, index] <- location[["y"]][index] - variance * residual
      }
    }
  }
  list(
    y                = location[["y"]],
    means            = means,
    variance         = variances,
    lower_tail       = location[["lower_tail"]],
    selection_context = .selection_joint_signed_context(
      setup = setup, signed_yi = location[["y"]]
    ),
    groups           = .data_selection_model(data)[["groups"]][["row_blocks"]],
    dependency_blocks = plan[["row_blocks"]]
  )
}


# Fix all other observed outcomes in the original publication event. For best
# selection this caps candidate p-values at the best retained p-value, including
# the mirrored z partition of a two-sided weight function.
.selection_joint_deleted_row_context <- function(context, y, sei, groups, row) {

  group <- which(vapply(groups, function(index) row %in% index, logical(1L)))
  if (length(group) != 1L) {
    stop("Selection publication groups must partition the observed rows.",
         call. = FALSE)
  }
  retained <- setdiff(groups[[group]], row)
  S <- nrow(context[["omega"]])
  original_rule <- rep_len(context[["vector_rule"]], S)
  focal_bin <- rep.int(context[["obs_bin"]][row], S)
  for (s in which(original_rule != 0L)) {
    focal_bin[s] <- .selection_joint_best_bin(y[row], sei[row], context, original_rule[s])
  }
  context <- .selection_joint_condition_event_context(
    context, matrix(y[retained], S, length(retained), byrow = TRUE), sei[retained]
  )
  context[["vector_rule"]] <- rep.int(0L, S)
  context <- BayesTools::selection_context_subset_observations(context, row)
  context[["obs_bin_by_sample"]] <- focal_bin
  context[["row_fields"]] <- unique(c(context[["row_fields"]], "obs_bin_by_sample"))
  context
}


# Conditional coordinate moments also admit a singular Gaussian source. The
# existing covariance policy owns its rank; no variance repair is performed.
.selection_deleted_gaussian_coordinates <- function(covariance, values) {

  K <- nrow(covariance)
  out <- list(mean = numeric(K), variance = numeric(K))
  if (all(covariance == 0)) return(out)
  factor <- tryCatch(chol(covariance), error = function(e) NULL)
  if (!is.null(factor)) {
    precision <- chol2inv(factor)
    out[["variance"]] <- 1 / diag(precision)
    out[["mean"]] <- values - out[["variance"]] * as.vector(precision %*% values)
    return(out)
  }
  for (row in seq_len(K)) {
    conditional <- .selection_deleted_gaussian_block(covariance, values, row)
    out[["mean"]][row] <- conditional[["mean"]]
    out[["variance"]][row] <- conditional[["covariance"]][1L, 1L]
  }
  out
}


# Delete the focal sampling error along with its observed outcome. Retain only
# e_-i, the other observed rows, and the model's declared random context. Given
# these, the integrated random coordinate and the sampling coordinate are two
# independent scalar Gaussians before the original selection correction.
.selection_conditioned_sampling_replicate_setup <- function(setup, state, draw, n) {

  rows <- rep(draw, n)
  current <- setup
  current[["S"]] <- n
  for (name in c("posterior_samples", "tau_within", "tau_between")) {
    current[[name]] <- setup[[name]][rows, , drop = FALSE]
  }
  if (!is.null(state[["candidate_factors"]])) {
    factors <- state[["candidate_factors"]]
    factors[["diagonal"]] <- factors[["diagonal"]][rows, , drop = FALSE]
    factors[["loadings"]] <- lapply(factors[["loadings"]], function(loading) {
      loading[rows, , , drop = FALSE]
    })
    current[["selection_conditioned_factors"]] <- factors
  }
  current
}


.selection_conditioned_sampling_independent_targets <- function(setup, state, selection, components) {

  plan <- .data_selection_execution_plan(setup[["data"]])
  independent_blocks <- all(lengths(plan[["row_blocks"]]) == 1L)
  if (!independent_blocks &&
      (length(plan[["factor_ranks"]]) != length(plan[["row_blocks"]]) ||
       any(plan[["factor_ranks"]] != 0L) || any(selection[["vector_rule"]] != 0L))) return(NULL)
  S <- setup[["S"]]
  K <- setup[["K"]]
  if (is.null(state)) {
    # A singleton deletion integrates its own sampling context. For ordinary
    # univariate sources the other contexts do not enter this scalar law, so
    # neither Gaussian auxiliary reconstruction nor S x K x K arrays are used.
    if (!independent_blocks || .is_data_random(setup[["data"]]) ||
        .is_data_multilevel(setup[["data"]]) || any(selection[["omega"]] <= 0) ||
        any(!selection[["kernel_mode"]] %in% c(SELKERNEL_NORMAL, SELKERNEL_STEP)) ||
        any(selection[["vector_rule"]] != 0L)) return(NULL)
    sampling_variances <- diag(plan[["sampling_covariance"]])
    tau_variances <- setup[["tau_within"]]^2
    if (nrow(tau_variances) == 1L && S > 1L) {
      tau_variances <- tau_variances[rep(1L, S), , drop = FALSE]
    }
    total_variances <- sweep(tau_variances, 2L, sampling_variances, "+")
    integrated_variances <- if (.selection_integrates_estimate(setup[["data"]])) {
      tau_variances
    } else matrix(0, S, K)
    baseline_mu <- setup[["mu"]]
    if (identical(components, "log_density")) {
      direction <- if (identical(setup[["effect_direction"]], "negative")) -1 else 1
      quadrature <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
      return(list(log_density = .Call(
        "RoBMA_selnorm_sampling_deletion_loglik_batch",
        as.numeric(direction * setup[["yi"]]), direction * baseline_mu,
        as.numeric(sampling_variances), integrated_variances, total_variances,
        as.numeric(setup[["selection_sei"]]), selection[["omega"]],
        as.numeric(selection[["z_lower"]]), as.numeric(selection[["z_upper"]]),
        as.integer(selection[["obs_bin"]]), as.integer(selection[["kernel_mode"]]),
        isTRUE(selection[["telescope_probabilities"]]),
        as.numeric(quadrature[["nodes"]]), as.numeric(quadrature[["log_weights"]]),
        as.integer(quadrature[["orders"]]), as.numeric(plan[["relative_tolerance"]]),
        PACKAGE = "RoBMA"
      )))
    }
  } else {
    sampling_variances <- diag(state[["sampling_covariance"]])
    baseline_mu <- state[["baseline_mu"]]
    integrated_variances <- .selection_covariance_diagonal(state[["integrated_covariance"]])
    total_variances <- .selection_covariance_diagonal(state[["total_covariance"]])
  }
  sampling_means <- matrix(0, S, K)
  if (!independent_blocks) {
    # Independent candidate effects make every other product factor cancel.
    # Sampling errors may still correlate: condition e_i on the retained e_-i.
    factor <- tryCatch(chol(state[["sampling_covariance"]]), error = function(e) NULL)
    if (is.null(factor)) return(NULL)
    precision <- chol2inv(factor)
    sampling_variances <- 1 / diag(precision)
    sampling_means <- state[["e"]] - sweep(state[["e"]] %*% precision,
      2L, sampling_variances, "*")
  }
  out <- stats::setNames(lapply(components, function(component) matrix(NA_real_, S, K)), components)
  direction <- if (identical(setup[["effect_direction"]], "negative")) -1 else 1
  y <- direction * setup[["yi"]]
  equal_weights <- apply(selection[["omega"]], 1L,
                         function(weights) all(weights == weights[1L]))
  for (row in seq_len(K)) {
    variance <- integrated_variances[, row]
    ordinary <- variance == 0 | selection[["kernel_mode"]] == SELKERNEL_NORMAL |
      equal_weights
    # When selection cancels in a dependent block, retain the ordinary joint
    # Gaussian deletion calculation below, which also integrates random context.
    normal_rows <- if (independent_blocks) ordinary else rep(FALSE, S)
    ordinary_mean <- direction * setup[["mu"]][normal_rows, row]
    ordinary_variance <- total_variances[normal_rows, row]
    ordinary_sd <- sqrt(ordinary_variance)
    if ("log_density" %in% components) out[["log_density"]][normal_rows, row] <- stats::dnorm(y[row], ordinary_mean, ordinary_sd, log = TRUE)
    if ("cdf" %in% components) out[["cdf"]][normal_rows, row] <- stats::pnorm(y[row], ordinary_mean, ordinary_sd, lower.tail = direction == 1)
    if ("log_lower" %in% components) out[["log_lower"]][normal_rows, row] <- stats::pnorm(y[row], ordinary_mean, ordinary_sd, lower.tail = direction == 1, log.p = TRUE)
    if ("log_upper" %in% components) out[["log_upper"]][normal_rows, row] <- stats::pnorm(y[row], ordinary_mean, ordinary_sd, lower.tail = direction != 1, log.p = TRUE)
    if ("mean" %in% components) out[["mean"]][normal_rows, row] <- direction * ordinary_mean
    if ("variance" %in% components) out[["variance"]][normal_rows, row] <- ordinary_variance
    active <- which(!ordinary)
    if (!length(active)) next
    sampling_variance <- sampling_variances[[row]]
    mean <- direction * (baseline_mu[, row] + sampling_means[, row])
    total <- sampling_variance + variance
    conditional_mean <- mean + sampling_variance / total * (y[row] - mean)
    conditional_sd <- sqrt(sampling_variance * variance / total)
    context <- BayesTools::selection_context_subset_observations(selection, row)
    log_weight <- .selection_joint_log_weight(matrix(y[row], S, 1L),
      setup[["selection_sei"]][row], context)
    previous <- NULL
    # Integrate inverse acceptance under E | Y instead of multiplying a
    # narrow Gaussian likelihood by a broad sampling-error quadrature rule.
    for (order in SELNORM_CLUSTER_QUADRATURE_ORDERS) {
      rule <- .gauss_hermite_nodes(order)
      n <- length(active)
      contexts <- BayesTools::selection_context_subset_rows(context, rep(active, order))
      sd <- matrix(rep(sqrt(variance[active]), order), n * order, 1L)
      sei <- setup[["selection_sei"]][row]
      current <- list()
      if ("log_density" %in% components) {
        means <- outer(conditional_mean[active], rep(1, order)) +
          outer(conditional_sd[active], rule[["nodes"]])
        log_mass <- .selection_step_log_norm_matrix(matrix(as.vector(means), n * order, 1L), sd, sei, contexts)
        terms <- sweep(-matrix(log_mass, n, order), 2L, rule[["log_weights"]], "+")
        current[["log_density"]] <- stats::dnorm(y[row], mean[active], sqrt(total[active]), log = TRUE) +
          log_weight[active] + .rowLogSumExps(terms)
      }
      if (any(components != "log_density")) {
        means <- outer(mean[active], rep(1, order)) +
          outer(rep(sqrt(sampling_variance), n), rule[["nodes"]])
        means <- matrix(as.vector(means), n * order, 1L)
        tails <- intersect(components, c("cdf", "log_lower", "log_upper"))
        tail_values <- list()
        if (length(tails)) log_mass <- .selection_step_log_norm_matrix(means, sd, sei, contexts)
        for (tail in tails) {
          lower_tail <- if (tail == "log_upper") direction != 1 else direction == 1
          key <- as.character(lower_tail)
          value <- tail_values[[key]]
          if (is.null(value)) {
            partial <- .selection_gaussian_event_mass(means, sd^2, sei, contexts, plan,
              lower = if (lower_tail) NULL else y[row], upper = if (lower_tail) y[row] else NULL)[["log_mass"]]
            terms <- sweep(matrix(partial - log_mass, n, order), 2L, rule[["log_weights"]], "+")
            value <- .rowLogSumExps(terms)
            tail_values[[key]] <- value
          }
          current[[tail]] <- if (tail == "cdf") exp(value) else value
        }
        if (any(c("mean", "variance") %in% components)) {
          moments <- .selection_step_moments_matrix(means, sd, sei, contexts)
          selected_mean <- as.vector(matrix(moments[["mean"]], n, order) %*% rule[["weights"]])
          if ("mean" %in% components) current[["mean"]] <- direction * selected_mean
          if ("variance" %in% components) current[["variance"]] <-
            as.vector(matrix(moments[["second"]], n, order) %*% rule[["weights"]]) - selected_mean^2
        }
      }
      if (!is.null(previous)) {
        error <- numeric(n)
        valid <- rep(TRUE, n)
        for (component in components) {
          value <- current[[component]]
          change <- if (component %in% c("log_density", "log_lower", "log_upper")) {
            abs(expm1(value - previous[[component]]))
          } else {
            scale <- if (component == "mean") sqrt(total[active]) else
              if (component == "variance") total[active] else abs(value)
            abs(value - previous[[component]]) / scale
          }
          same <- !is.na(value) & !is.na(previous[[component]]) & value == previous[[component]]
          change[same] <- 0
          error <- pmax(error, change)
          valid <- valid & !is.na(value) & value != Inf
          if (component == "variance") valid <- valid & value >= 0
        }
        accepted <- which(valid & is.finite(error) & error <= plan[["relative_tolerance"]])
        if (length(accepted)) {
          draws <- active[accepted]
          for (component in components) out[[component]][draws, row] <- current[[component]][accepted]
        }
        failed <- setdiff(seq_len(n), accepted)
        if (!length(failed)) break
        active <- active[failed]
        current <- lapply(current, function(value) value[failed])
      }
      previous <- current
    }
  }
  out
}


.selection_conditioned_sampling_estimate_targets <- function(setup, components) {

  S <- setup[["S"]]
  K <- setup[["K"]]
  direction <- if (identical(setup[["effect_direction"]], "negative")) -1 else 1
  y <- direction * setup[["yi"]]
  plan <- .data_selection_execution_plan(setup[["data"]])
  selection <- .selection_joint_signed_context(setup, y)
  out <- .selection_conditioned_sampling_independent_targets(setup, NULL, selection, components)
  state <- NULL
  if (is.null(out)) {
    state <- .selection_conditioned_sampling_state(setup)
    out <- .selection_conditioned_sampling_independent_targets(setup, state, selection, components)
  }
  if (is.null(out)) {
    out <- stats::setNames(lapply(components, function(component) matrix(NA_real_, S, K)),
                          components)
  }
  if (!any(vapply(out, anyNA, logical(1L)))) {
    attr(out, "dependency_blocks") <- plan[["row_blocks"]]
    return(out)
  }
  if (is.null(state)) state <- .selection_conditioned_sampling_state(setup)
  baseline <- direction * state[["baseline_mu"]]
  errors <- direction * state[["e"]]
  V <- state[["sampling_covariance"]]
  groups <- .data_selection_model(setup[["data"]])[["groups"]][["row_blocks"]]
  tolerance <- plan[["relative_tolerance"]]
  subdivisions <- as.integer(max(1, floor(plan[["max_points_per_scramble"]] / 21)))
  remedy <- paste0("Increase 'max_points_per_scramble' in ",
                   "'selection_control = set_selection_likelihood_control()'.")
  ordinary_covariance <- state[["total_covariance"]]
  for (draw in seq_len(S)) {
    if (!any(vapply(out, function(value) anyNA(value[draw, ]), logical(1L)))) next
    covariance <- .selection_covariance_draw(state[["integrated_covariance"]], draw)
    ordinary <- all(covariance == 0) || selection[["kernel_mode"]][draw] == SELKERNEL_NORMAL ||
      all(selection[["omega"]][draw, ] == selection[["omega"]][draw, 1L])
    if (ordinary) {
      fixed <- direction * setup[["mu"]][draw, ]
      gaussian <- .selection_deleted_gaussian_coordinates(
        .selection_covariance_draw(ordinary_covariance, draw), y - fixed
      )
      means <- fixed + gaussian[["mean"]]
      sd <- sqrt(gaussian[["variance"]])
      pending <- Reduce(`|`, lapply(out, function(value) is.na(value[draw, ])))
      if (any(sd[pending] == 0)) {
        stop("Selection estimate deletion is unavailable because the deleted outcome is determined by retained outcomes. Use a larger deletion unit.",
             call. = FALSE)
      }
      for (component in components) {
        rows <- which(is.na(out[[component]][draw, ]))
        out[[component]][draw, rows] <- switch(component,
          log_density = stats::dnorm(y[rows], means[rows], sd[rows], log = TRUE),
          cdf = stats::pnorm(y[rows], means[rows], sd[rows], lower.tail = direction == 1),
          log_lower = stats::pnorm(y[rows], means[rows], sd[rows], lower.tail = direction == 1, log.p = TRUE),
          log_upper = stats::pnorm(y[rows], means[rows], sd[rows], lower.tail = direction != 1, log.p = TRUE),
          mean = direction * means[rows],
          variance = gaussian[["variance"]][rows]
        )
      }
      next
    }
    sampling <- .selection_deleted_gaussian_coordinates(V, errors[draw, ])
    random <- .selection_deleted_gaussian_coordinates(
      covariance, y - baseline[draw, ] - errors[draw, ]
    )
    for (row in seq_len(K)) {
      if (!any(vapply(out, function(value) is.na(value[draw, row]), logical(1L)))) next
      v_e <- sampling[["variance"]][row]
      v_u <- random[["variance"]][row]
      variance <- v_e + v_u
      if (variance == 0) {
        stop("Selection estimate deletion is unavailable because the deleted outcome is determined by retained sampling errors and outcomes. Use 'known_sampling_variance = \"integrate\"' when fitting to obtain sampling-marginal deletion scores.",
             call. = FALSE)
      }
      m_e <- sampling[["mean"]][row]
      location <- baseline[draw, row] + random[["mean"]][row]
      mean <- location + m_e
      context <- .selection_joint_deleted_row_context(selection, y, setup[["selection_sei"]], groups, row)
      context <- BayesTools::selection_context_subset_rows(context, draw)
      sei <- setup[["selection_sei"]][row]
      cuts <- sort(unique(c(context[["z_lower"]], context[["z_upper"]]) * sei))
      cuts <- cuts[is.finite(cuts)]
      cache <- new.env(parent = emptyenv())
      log_normalizer <- function(e) {

        key <- paste(format(e, digits = 17), collapse = "/")
        if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
        n <- length(e)
        current <- .selection_conditioned_sampling_replicate_setup(setup, state, draw, n)
        means <- matrix(baseline[draw, ] + errors[draw, ], n, K, byrow = TRUE)
        means[, row] <- baseline[draw, row] + e
        result <- .selection_conditioned_sampling_normalizer(current, direction * means)
        mass <- result[["log_mass"]]
        if (any(!is.finite(mass))) {
          stop("Selection estimate deletion is unavailable because a retained sampling context has zero selection mass.",
               call. = FALSE)
        }
        assign(key, mass, envir = cache)
        mass
      }
      anchor <- log_normalizer(errors[draw, row])
      integrate_probabilities <- function(fun, boundaries, moment_scale = NULL) {

        parts <- length(boundaries) - 1L
        value <- error <- 0
        for (part in seq_len(parts)) {
          integral <- stats::integrate(
            fun,
            lower = boundaries[part], upper = boundaries[part + 1L],
            subdivisions = subdivisions, rel.tol = tolerance / (2 * parts),
            abs.tol = if (is.null(moment_scale)) 0 else tolerance * moment_scale / (2 * parts),
            stop.on.error = FALSE
          )
          value <- value + integral[["value"]]
          error <- error + integral[["abs.error"]]
          if (!identical(integral[["message"]], "OK") || !is.finite(value)) {
            stop(paste0("Selection estimate deletion was rejected by diagnostics: scalar integration reported ",
                        integral[["message"]], ". ", remedy), call. = FALSE)
          }
        }
        error_scale <- if (is.null(moment_scale)) abs(value) else moment_scale
        relative_error <- if (error_scale == 0 && error == 0) 0 else error / error_scale
        if (!is.finite(relative_error) || relative_error > tolerance) {
          stop(paste0("Selection estimate deletion was rejected by diagnostics: relative integration error was ",
                      format(relative_error, digits = 4), ". ", remedy), call. = FALSE)
        }
        value
      }
      integrate_gaussian <- function(fun, mean_e, variance_e, moment_scale = NULL) {

        if (variance_e == 0) return(as.numeric(fun(mean_e)))
        sd_e <- sqrt(variance_e)
        z <- (cuts - location - mean_e) / sd_e
        lower_boundaries <- sort(unique(c(0, .5, stats::pnorm(z[z < 0]))))
        upper_boundaries <- sort(unique(c(0, .5, stats::pnorm(z[z > 0], lower.tail = FALSE))))
        # Work in each tail's own probabilities. A narrow interval next to
        # CDF = 1 otherwise rounds interior quadrature points to exactly one,
        # sending infinite locations into the Gaussian selection kernel.
        half_scale <- if (is.null(moment_scale)) NULL else moment_scale / 2
        integrate_probabilities(function(p) fun(mean_e + sd_e * stats::qnorm(p)),
          lower_boundaries, half_scale) +
          integrate_probabilities(function(p) fun(mean_e + sd_e * stats::qnorm(p, lower.tail = FALSE)),
            upper_boundaries, half_scale)
      }
      weighted_mass <- function(e, lower = NULL, upper = NULL, moment = NULL) {

        n <- length(e)
        contexts <- BayesTools::selection_context_subset_rows(context, rep(1L, n))
        means <- matrix(location + e, n, 1L)
        mass <- if (v_u == 0) {
          value <- .selection_joint_log_weight(means, sei, contexts)
          outside <- rep(FALSE, n)
          if (!is.null(lower)) outside <- outside | as.numeric(means) < lower
          if (!is.null(upper)) outside <- outside | as.numeric(means) > upper
          value[outside] <- -Inf
          value
        } else .selection_gaussian_event_mass(
          means, matrix(v_u, n, 1L), sei, contexts, plan,
          lower = lower, upper = upper
        )[["log_mass"]]
        value <- exp(mass + anchor - log_normalizer(e))
        if (!is.null(moment)) {
          moments <- if (v_u == 0) {
            list(mean = means, second = means^2)
          } else .selection_step_moments_matrix(means, matrix(sqrt(v_u), n, 1L), sei, contexts)
          value <- value * as.numeric(moments[[moment]])
        }
        as.numeric(value)
      }
      singleton_block <- which(vapply(plan[["row_blocks"]], function(index) {
        length(index) == 1L && index == row
      }, logical(1L)))
      normalizer <- if (length(singleton_block)) {
        # The complete independent event integrates to one under the original
        # sampling-error law. Only the numerical anchor remains in this scale.
        exp(state[["log_normalizer"]][draw, singleton_block])
      } else integrate_gaussian(weighted_mass, m_e, v_e)
      if (!is.finite(normalizer) || normalizer <= 0) {
        stop("Selection estimate deletion is unavailable because its conditional selection event cannot be normalized.",
             call. = FALSE)
      }
      if ("log_density" %in% components) {
        conditional_mean <- m_e + v_e / variance * (y[row] - mean)
        conditional_variance <- v_e * v_u / variance
        reciprocal <- integrate_gaussian(function(e) exp(anchor - log_normalizer(e)),
                                          conditional_mean, conditional_variance)
        weight <- context[["omega"]][1L, context[["obs_bin_by_sample"]][1L]]
        out[["log_density"]][draw, row] <- stats::dnorm(y[row], mean, sqrt(variance), log = TRUE) +
          log(weight) + log(reciprocal) - log(normalizer)
      }
      for (tail in intersect(components, c("cdf", "log_lower", "log_upper"))) {
        lower_tail <- if (tail == "log_upper") direction != 1 else direction == 1
        probability <- integrate_gaussian(function(e) weighted_mass(e,
          lower = if (lower_tail) NULL else y[row], upper = if (lower_tail) y[row] else NULL),
          m_e, v_e) / normalizer
        log_probability <- log(probability)
        if (probability == 0) {
          # A positive tail may be too small for ordinary probability-space
          # mixing. Factor out its Gaussian tail and integrate the bounded
          # selection density ratio under that truncated Gaussian instead.
          log_base <- stats::pnorm(y[row], mean, sqrt(variance),
                                    lower.tail = lower_tail, log.p = TRUE)
          cut_probabilities <- stats::pnorm(cuts, mean, sqrt(variance),
                                            lower.tail = lower_tail, log.p = TRUE)
          boundaries <- sort(unique(c(0, 1,
            exp(cut_probabilities[cut_probabilities < log_base] - log_base))))
          ratio <- integrate_probabilities(function(p) {
            proposed <- stats::qnorm(log_base + log(p), mean, sqrt(variance),
                                      lower.tail = lower_tail, log.p = TRUE)
            contexts <- BayesTools::selection_context_subset_rows(context, rep(1L, length(p)))
            log_weights <- .selection_joint_log_weight(matrix(proposed, length(p), 1L), sei, contexts)
            vapply(seq_along(proposed), function(i) {
              conditional_mean <- m_e + v_e / variance * (proposed[i] - mean)
              reciprocal <- integrate_gaussian(function(e) exp(anchor - log_normalizer(e)),
                conditional_mean, v_e * v_u / variance)
              exp(log_weights[i] + log(reciprocal) - log(normalizer))
            }, numeric(1L))
          }, boundaries)
          log_probability <- log_base + log(ratio)
        }
        out[[tail]][draw, row] <- if (tail == "cdf") exp(log_probability) else log_probability
      }
      if (any(c("mean", "variance") %in% components)) {
        # Scale moment error by the Gaussian SD/second moment, so a zero
        # selected mean remains a valid result rather than a relative-error
        # singularity. Density and event-mass criteria remain relative to mass.
        selected_mean <- integrate_gaussian(function(e) weighted_mass(e, moment = "mean"),
          m_e, v_e, moment_scale = normalizer * sqrt(variance)) / normalizer
        if ("mean" %in% components) out[["mean"]][draw, row] <- direction * selected_mean
        if ("variance" %in% components) {
          second <- integrate_gaussian(function(e) weighted_mass(e, moment = "second"),
            m_e, v_e, moment_scale = normalizer * (variance + mean^2)) / normalizer
          selected_variance <- second - selected_mean^2
          if (!is.finite(selected_variance) || selected_variance < 0) {
            stop("Selected-normal predictive variance is invalid.", call. = FALSE)
          }
          out[["variance"]][draw, row] <- selected_variance
        }
      }
    }
  }
  attr(out, "dependency_blocks") <- plan[["row_blocks"]]
  out
}


.selection_joint_estimate_targets <- function(
    setup, components = c("log_density", "cdf", "log_lower", "log_upper", "mean", "variance")) {

  if (.selection_retains_sampling(setup[["data"]])) {
    return(.selection_conditioned_sampling_estimate_targets(setup, components))
  }

  conditional <- .selection_joint_conditional_summary_from_setup(setup)
  out <- stats::setNames(lapply(components, function(x) {
    matrix(NA_real_, setup[["S"]], setup[["K"]])
  }), components)
  model <- .data_selection_model(setup[["data"]])
  branches <- model[["branches"]][model[["active_branches"]]]
  product <- length(branches) > 0L && all(vapply(branches, function(branch) {
    identical(branch[["weight_rule"]], "product")
  }, logical(1L)))
  if (product) {
    context <- conditional[["selection_context"]]
    mean <- conditional[["means"]]
    variance <- conditional[["variance"]]
    sd <- sqrt(variance)
    sei <- setup[["selection_sei"]]
    y <- conditional[["y"]]
    if ("log_density" %in% components) {
      out[["log_density"]] <- .selection_joint_singleton_loglik_matrix(
        yi = y, means = mean, variances = variance, sei = sei,
        selection_context = context
      )
      # Preserve the existing zero-event error; an observed zero weight alone
      # remains a valid -Inf score. Ordinary finite rows need no second mass.
      if (any(!is.finite(out[["log_density"]]))) {
        normalizer <- .selection_step_log_norm_matrix(mean, sd, sei, context)
        if (any(!is.finite(normalizer))) {
          stop("Selected row deletion is unavailable because its conditional selection event cannot be normalized.", call. = FALSE)
        }
      }
    }
    if ("cdf" %in% components) {
      out[["cdf"]] <- .selection_step_cdf_matrix(
        q = y, mean = mean, sd = sd, sei = sei,
        selection_context = context, lower.tail = conditional[["lower_tail"]]
      )
    }
    if (any(c("mean", "variance") %in% components)) {
      moments <- .selection_step_moments_matrix(mean, sd, sei, context)
      if ("mean" %in% components) {
        out[["mean"]] <- if (conditional[["lower_tail"]]) {
          moments[["mean"]]
        } else {
          -moments[["mean"]]
        }
      }
      if ("variance" %in% components) {
        variance <- moments[["second"]] - moments[["mean"]]^2
        if (any(!is.finite(variance)) || any(variance < 0)) {
          stop("Selected-normal predictive variance is invalid.", call. = FALSE)
        }
        out[["variance"]] <- variance
      }
    }
    components <- intersect(components, c("log_lower", "log_upper"))
  }
  for (row in if (length(components)) seq_len(setup[["K"]]) else integer()) {
    context <- .selection_joint_deleted_row_context(
      context = conditional[["selection_context"]],
      y       = conditional[["y"]],
      sei     = setup[["selection_sei"]],
      groups  = conditional[["groups"]],
      row     = row
    )
    mean <- conditional[["means"]][, row, drop = FALSE]
    sd <- sqrt(conditional[["variance"]][, row, drop = FALSE])
    sei <- setup[["selection_sei"]][row]
    y <- conditional[["y"]][row]
    if (any(c("log_density", "log_lower", "log_upper") %in% components)) {
      normalizer <- .selection_gaussian_event_mass(
        mean, sd^2, sei, context, .data_selection_execution_plan(setup[["data"]])
      )[["log_mass"]]
      if (any(!is.finite(normalizer))) {
        stop("Selected row deletion is unavailable because its conditional selection event cannot be normalized.", call. = FALSE)
      }
    }
    if ("log_density" %in% components) {
      log_weight <- numeric(setup[["S"]])
      selected <- which(context[["kernel_mode"]] != SELKERNEL_NORMAL)
      log_weight[selected] <- log(context[["omega"]][cbind(
        selected, context[["obs_bin_by_sample"]][selected]
      )])
      out[["log_density"]][, row] <- stats::dnorm(y, mean, sd, log = TRUE) +
        log_weight - normalizer
    }
    for (tail in intersect(components, c("cdf", "log_lower", "log_upper"))) {
      lower_tail <- if (tail == "log_upper") !conditional[["lower_tail"]] else
        conditional[["lower_tail"]]
      out[[tail]][, row] <- if (tail == "cdf") {
        .selection_step_cdf_matrix(
          q = y, mean = mean, sd = sd, sei = sei,
          selection_context = context, lower.tail = lower_tail
        )
      } else {
        .selection_gaussian_event_mass(
          mean, sd^2, sei, context, .data_selection_execution_plan(setup[["data"]]),
          lower = if (lower_tail) NULL else y,
          upper = if (lower_tail) y else NULL
        )[["log_mass"]] - normalizer
      }
    }
    if (any(c("mean", "variance") %in% components)) {
      moments <- .selection_step_moments_matrix(
        mean = mean, sd = sd, sei = sei, selection_context = context
      )
      if ("mean" %in% components) {
        out[["mean"]][, row] <- if (conditional[["lower_tail"]]) {
          moments[["mean"]]
        } else {
          -moments[["mean"]]
        }
      }
      if ("variance" %in% components) {
        variance <- moments[["second"]] - moments[["mean"]]^2
        if (any(!is.finite(variance)) || any(variance < 0)) {
          stop("Selected-normal predictive variance is invalid.", call. = FALSE)
        }
        out[["variance"]][, row] <- variance
      }
    }
  }
  attr(out, "dependency_blocks") <- conditional[["dependency_blocks"]]
  out
}


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
      !.is_data_joint_selection(setup[["data"]])) {
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

  if (.setup_uses_joint_selection_likelihood(setup)) {
    selected <- .selection_joint_estimate_targets(setup, "log_density")
    log_lik <- selected[["log_density"]]
    if (isTRUE(add_dependency_metadata)) {
      attr(log_lik, "RoBMA_dependency_blocks") <- attr(selected, "dependency_blocks")
    }
    return(log_lik)
  }
  plan <- .estimate_normal_covariance_target_plan_from_setup(setup)
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
      !.is_data_joint_selection(setup[["data"]])) {
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

  if (.setup_uses_joint_selection_likelihood(setup)) {
    return(.selection_joint_loglik_from_setup(setup))
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
      !.is_data_joint_selection(setup[["data"]])) {
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

  if (.setup_uses_joint_selection_likelihood(setup)) {
    return(.selection_joint_estimate_targets(setup, components))
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

  extra_variance <- if (.is_data_joint_selection(data) &&
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
