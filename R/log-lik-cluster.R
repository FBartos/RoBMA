# ============================================================================ #
# Cluster-Unit Log-Likelihood Dispatch
# ============================================================================ #

.log_lik_cluster_from_setup <- function(setup) {

  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  cluster_setup     <- setup

  if (.setup_uses_joint_selection_likelihood(setup) &&
      .selection_retains_sampling(setup[["data"]])) {
    return(.selection_conditioned_sampling_deletion_loglik(setup, setup[["cluster"]]))
  }

  if (.setup_uses_joint_selection_likelihood(setup)) {
    plan <- .data_selection_execution_plan(setup[["data"]])
    blocks <- plan[["row_blocks"]]
    cluster_blocks <- unname(setup[["cluster"]])
    membership <- integer(setup[["K"]])
    for (cluster in seq_along(cluster_blocks)) membership[cluster_blocks[[cluster]]] <- cluster
    if (all(vapply(blocks, function(rows) length(unique(membership[rows])) == 1L,
                   logical(1L)))) {
      # Sum complete independent event contributions, never deletion row scores.
      log_lik <- .selection_joint_block_loglik_from_setup(setup)
      out <- matrix(0, setup[["S"]], length(cluster_blocks))
      for (block in seq_along(blocks)) {
        cluster <- membership[blocks[[block]][1L]]
        out[, cluster] <- out[, cluster] + log_lik[, block]
      }
      return(out)
    }
    return(.selection_joint_deletion_loglik_from_setup(setup, cluster_blocks))
  }

  if (outcome_type == "norm" && !is_weightfunction &&
      is.null(setup[["weights"]])) {
    if (setup[["effect_direction"]] == "negative") {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    return(.log_lik_cluster_norm_analytic(
      setup = cluster_setup,
      yi    = yi,
      vi    = sei^2
    ))
  }

  if (outcome_type == "norm") {
    if (setup[["effect_direction"]] == "negative" && !is_weightfunction) {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    selection_context <- if (is_weightfunction) {
      .selection_context_from_parts(
        fit                  = setup[["fit"]],
        data                 = setup[["data"]],
        priors               = setup[["priors"]],
        posterior_samples    = setup[["posterior_samples"]],
        effect_direction     = setup[["effect_direction"]]
      )
    } else {
      NULL
    }

    return(.log_lik_cluster_norm_quadrature(
      setup             = cluster_setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = is_weightfunction,
      selection_context = selection_context
    ))
  }

  if (outcome_type %in% c("bin", "pois")) {
    return(.log_lik_cluster_glmm(
      setup        = cluster_setup,
      data         = setup[["data"]],
      priors       = setup[["priors"]],
      outcome_type = outcome_type
    ))
  }

  stop("Unsupported outcome type for cluster-unit log-likelihood.",
       call. = FALSE)
}


.selection_joint_deletion_loglik_from_setup <- function(setup, units) {

  if (.selection_retains_sampling(setup[["data"]])) {
    return(.selection_conditioned_sampling_deletion_loglik(setup, units))
  }

  plan <- .data_selection_execution_plan(setup[["data"]])
  location <- .estimate_normal_covariance_target_location_from_setup(setup)
  selection <- .selection_joint_signed_context(setup, location[["y"]])
  factors <- .selection_joint_random_factor_samples(setup)
  random_covariance <- if (is.null(factors)) .selection_joint_random_covariance_samples(setup) else NULL
  full_log_lik <- NULL
  out <- matrix(0, setup[["S"]], length(units))
  for (block in seq_along(plan[["row_blocks"]])) {
    index <- plan[["row_blocks"]][[block]]
    covariance <- .selection_joint_covariance_lower(setup, block, random_covariance, factors)
    mean <- location[["means"]][, index, drop = FALSE]
    y <- location[["y"]][index]
    sei <- setup[["selection_sei"]][index]
    context <- BayesTools::selection_context_subset_observations(selection, index)
    full <- NULL
    for (unit in seq_along(units)) {
      deleted <- which(index %in% units[[unit]])
      if (!length(deleted)) next
      retained <- setdiff(seq_along(index), deleted)
      if (!length(retained)) {
        if (is.null(full_log_lik)) full_log_lik <- .selection_joint_block_loglik_from_setup(setup)
        out[, unit] <- out[, unit] + full_log_lik[, block]
        next
      }
      if (is.null(full)) {
        full <- .selection_joint_event_numerator(
          mean, covariance, sei, context, plan, seq_along(index), y
        )[["log_numerator"]]
        if (any(!is.finite(full))) {
          stop("Selection deletion is unavailable because the observed outcomes have zero or invalid event density.",
               call. = FALSE)
        }
      }
      partial <- .selection_joint_checked_event(
        compute = function(control, rows = NULL) {
          if (is.null(rows)) rows <- seq_len(setup[["S"]])
          .selection_joint_event_numerator(
            mean[rows, , drop = FALSE], covariance[rows, , drop = FALSE], sei,
            BayesTools::selection_context_subset_rows(context, rows), control,
            retained, y[retained]
          )
        },
        execution_plan = plan, subject = "Selection deletion integral",
        remedy = paste0("Increase 'max_points_per_scramble' in ",
                        "'selection_control = set_selection_likelihood_control()' when fitting.")
      )
      if (any(!is.finite(partial[["log_numerator"]]))) {
        stop("Selection deletion is unavailable because the retained outcomes have zero or invalid event density.", call. = FALSE)
      }
      out[, unit] <- out[, unit] + full - partial[["log_numerator"]]
    }
  }
  out
}


.selection_deleted_gaussian_block <- function(covariance, values, deleted) {

  retained <- setdiff(seq_along(values), deleted)
  dimension <- length(deleted)
  if (!length(retained)) {
    return(list(mean = numeric(dimension), covariance = covariance[deleted, deleted, drop = FALSE]))
  }
  if (all(covariance == 0)) {
    return(list(mean = numeric(dimension), covariance = matrix(0, dimension, dimension)))
  }
  kept <- covariance[retained, retained, drop = FALSE]
  cross <- covariance[retained, deleted, drop = FALSE]
  decomposition <- .covariance_factorization(kept)
  if (!.covariance_is_positive_semidefinite(decomposition)) {
    stop("Selection deletion covariance must be positive semidefinite.", call. = FALSE)
  }
  factor <- .covariance_cholesky(decomposition)
  if (!is.null(factor)) {
    coefficient <- backsolve(factor, forwardsolve(t(factor), cross))
  } else {
    positive <- decomposition[["spectral_values"]] > 0
    vectors <- decomposition[["eigenvectors"]][, positive, drop = FALSE]
    coefficient <- if (any(positive)) {
      vectors %*% (crossprod(vectors, cross) / decomposition[["spectral_values"]][positive])
    } else matrix(0, length(retained), dimension)
  }
  full <- .covariance_factorization(covariance)
  root <- .covariance_sampling_factor(full)
  if (is.null(root)) {
    stop("Selection deletion covariance must be positive semidefinite.", call. = FALSE)
  }
  residual_root <- root[, deleted, drop = FALSE] - root[, retained, drop = FALSE] %*% coefficient
  conditional_covariance <- crossprod(residual_root)
  if (!.covariance_is_numerically_positive_definite(full)) {
    kept_rank <- sum(decomposition[["spectral_values"]] > 0)
    if (kept_rank == sum(full[["spectral_values"]] > 0)) conditional_covariance[,] <- 0
  }
  list(mean = as.vector(crossprod(coefficient, values[retained])),
       covariance = conditional_covariance)
}


.selection_conditioned_sampling_deletion_loglik <- function(setup, units) {

  if (all(lengths(units) == 1L)) {
    scores <- .selection_conditioned_sampling_estimate_targets(setup, "log_density")[["log_density"]]
    return(scores[, unlist(units, use.names = FALSE), drop = FALSE])
  }
  state <- .selection_conditioned_sampling_state(setup)
  S <- setup[["S"]]
  K <- setup[["K"]]
  direction <- if (identical(setup[["effect_direction"]], "negative")) -1 else 1
  y <- direction * setup[["yi"]]
  baseline <- direction * state[["baseline_mu"]]
  errors <- direction * state[["e"]]
  V <- state[["sampling_covariance"]]
  plan <- .data_selection_execution_plan(setup[["data"]])
  selection <- .selection_joint_signed_context(setup, y)
  groups <- .data_selection_model(setup[["data"]])[["groups"]][["row_blocks"]]
  designs <- new.env(parent = emptyenv())
  out <- matrix(NA_real_, S, length(units))
  ordinary_covariance <- state[["total_covariance"]]
  for (draw in seq_len(S)) {
    covariance <- .selection_covariance_draw(state[["integrated_covariance"]], draw)
    ordinary <- all(covariance == 0) || selection[["kernel_mode"]][draw] == SELKERNEL_NORMAL ||
      all(selection[["omega"]][draw, ] == selection[["omega"]][draw, 1L])
    context <- BayesTools::selection_context_subset_rows(selection, draw)
    weight_groups <- if (context[["vector_rule"]] == 0L) list(seq_len(K)) else groups
    log_weights <- function(values) {

      weights <- numeric(nrow(values))
      for (group in weight_groups) {
        group_context <- BayesTools::selection_context_subset_observations(context, group)
        group_context <- BayesTools::selection_context_subset_rows(group_context, rep(1L, nrow(values)))
        weights <- weights + .selection_joint_log_weight(values[, group, drop = FALSE],
          setup[["selection_sei"]][group], group_context)
      }
      weights
    }
    observed_weight <- log_weights(matrix(y, 1L))
    for (unit in seq_along(units)) {
      deleted <- units[[unit]]
      dimension <- length(deleted)
      if (ordinary) {
        fixed <- direction * setup[["mu"]][draw, ]
        gaussian <- .selection_deleted_gaussian_block(.selection_covariance_draw(ordinary_covariance, draw),
                                                      y - fixed, deleted)
        mean <- fixed[deleted] + gaussian[["mean"]]
        total <- gaussian[["covariance"]]
      } else {
        sampling <- .selection_deleted_gaussian_block(V, errors[draw, ], deleted)
        random <- .selection_deleted_gaussian_block(covariance,
          y - baseline[draw, ] - errors[draw, ], deleted)
        mean <- baseline[draw, deleted] + sampling[["mean"]] + random[["mean"]]
        total <- sampling[["covariance"]] + random[["covariance"]]
      }
      total_factor <- .covariance_cholesky(.covariance_factorization(total))
      if (is.null(total_factor)) {
        stop("Selection deletion is unavailable because the deleted outcome vector has no density conditional on retained sampling errors and outcomes. Use 'known_sampling_variance = \"integrate\"' when fitting to obtain sampling-marginal deletion scores.",
             call. = FALSE)
      }
      residual <- y[deleted] - mean
      log_base <- sum(stats::dnorm(forwardsolve(t(total_factor), residual), log = TRUE)) -
        sum(log(diag(total_factor)))
      if (ordinary) {
        out[draw, unit] <- log_base
        next
      }
      sampling_root <- .covariance_sampling_factor(.covariance_factorization(sampling[["covariance"]]))
      random_root <- .covariance_sampling_factor(.covariance_factorization(random[["covariance"]]))
      if (is.null(sampling_root) || is.null(random_root)) {
        stop("Selection deletion covariance must be positive semidefinite.", call. = FALSE)
      }
      gain <- sampling[["covariance"]] %*% chol2inv(total_factor)
      posterior_mean <- sampling[["mean"]] + as.vector(gain %*% residual)
      # A conditional-draw root avoids subtracting nearly equal covariances.
      posterior_root <- rbind(sampling_root %*% t(diag(dimension) - gain),
                              -random_root %*% t(gain))
      affected <- vapply(plan[["row_blocks"]], function(index) any(index %in% deleted), logical(1L))
      complete_blocks <- all(vapply(plan[["row_blocks"]][affected], function(index) {
        all(index %in% deleted)
      }, logical(1L)))
      log_denominator <- NULL
      if (complete_blocks) {
        # Whole independent events have unit integrated selected probability.
        # Retained outside events contribute only their fixed W/A ratio.
        outside <- setdiff(seq_len(K), deleted)
        outside_weight <- 0
        for (group in weight_groups) {
          rows <- intersect(group, outside)
          if (length(rows)) outside_weight <- outside_weight + .selection_joint_log_weight(
            matrix(y[rows], 1L), setup[["selection_sei"]][rows],
            BayesTools::selection_context_subset_observations(context, rows))
        }
        log_denominator <- outside_weight - sum(state[["log_normalizer"]][draw, !affected])
      }
      normalizer <- function(error_values) {

        n <- nrow(error_values)
        current <- .selection_conditioned_sampling_replicate_setup(setup, state, draw, n)
        means <- matrix(baseline[draw, ] + errors[draw, ], n, K, byrow = TRUE)
        means[, deleted] <- sweep(error_values, 2L, baseline[draw, deleted], "+")
        mass <- .selection_conditioned_sampling_normalizer(current, direction * means)[["log_mass"]]
        if (any(!is.finite(mass))) {
          stop("Selection deletion is unavailable because a retained sampling context has zero selection mass.",
               call. = FALSE)
        }
        mass
      }
      points <- plan[["points_per_scramble"]]
      scrambles <- plan[["scrambles"]]
      previous <- NULL
      repeat {
        key <- paste(dimension, points, sep = "/")
        if (!exists(key, envir = designs, inherits = FALSE)) {
          uniforms <- BayesTools::selection_qmc_design(2L * dimension, points, scrambles, plan[["seed"]])
          assign(key, matrix(stats::qnorm(uniforms), points * scrambles, 2L * dimension), envir = designs)
        }
        normals <- get(key, envir = designs, inherits = FALSE)
        denominator <- if (is.null(log_denominator)) {
          e <- sweep(normals[, seq_len(dimension), drop = FALSE] %*% sampling_root,
                     2L, sampling[["mean"]], "+")
          u <- sweep(normals[, dimension + seq_len(dimension), drop = FALSE] %*% random_root,
                     2L, random[["mean"]], "+")
          proposed <- matrix(y, nrow(e), K, byrow = TRUE)
          proposed[, deleted] <- sweep(e + u, 2L, baseline[draw, deleted], "+")
          log_weights(proposed) - normalizer(e)
        } else rep(log_denominator, points * scrambles)
        e_given_y <- sweep(normals %*% posterior_root, 2L, posterior_mean, "+")
        numerator <- -normalizer(e_given_y)
        denominator_scale <- max(denominator)
        numerator_scale <- max(numerator)
        denominator <- rowMeans(matrix(exp(denominator - denominator_scale), scrambles, points))
        numerator <- rowMeans(matrix(exp(numerator - numerator_scale), scrambles, points))
        log_ratio <- log(mean(numerator)) - log(mean(denominator)) + numerator_scale - denominator_scale
        relative_error <- stats::sd(numerator / mean(numerator) - denominator / mean(denominator)) /
          sqrt(scrambles)
        if (!is.null(previous)) relative_error <- max(relative_error, abs(expm1(log_ratio - previous)))
        if (is.finite(log_ratio) && is.finite(relative_error) &&
            relative_error <= plan[["relative_tolerance"]]) {
          out[draw, unit] <- log_base + observed_weight + log_ratio
          break
        }
        if (points >= plan[["max_points_per_scramble"]]) {
          stop(paste0("Selection deletion was rejected by diagnostics: relative integration error was ",
                      format(relative_error, digits = 4), ". Increase 'max_points_per_scramble' in ",
                      "'selection_control = set_selection_likelihood_control()'."), call. = FALSE)
        }
        previous <- log_ratio
        points <- min(2L * points, plan[["max_points_per_scramble"]])
      }
    }
  }
  out
}


.log_lik_cluster_sum_from_setup <- function(setup) {

  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  cluster_setup     <- setup

  if (.setup_uses_joint_selection_likelihood(setup)) {
    return(.selection_joint_loglik_from_setup(setup))
  }

  if (outcome_type == "norm" && !is_weightfunction &&
      is.null(setup[["weights"]])) {
    if (setup[["effect_direction"]] == "negative") {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    return(.log_lik_cluster_norm_analytic_sum(
      setup = cluster_setup,
      yi    = yi,
      vi    = sei^2
    ))
  }

  if (outcome_type == "norm") {
    if (setup[["effect_direction"]] == "negative" && !is_weightfunction) {
      cluster_setup[["mu"]] <- -cluster_setup[["mu"]]
      yi                    <- -yi
    }

    selection_context <- if (is_weightfunction) {
      .selection_context_from_parts(
        fit                  = setup[["fit"]],
        data                 = setup[["data"]],
        priors               = setup[["priors"]],
        posterior_samples    = setup[["posterior_samples"]],
        effect_direction     = setup[["effect_direction"]]
      )
    } else {
      NULL
    }

    return(.log_lik_cluster_norm_quadrature_sum(
      setup             = cluster_setup,
      yi                = yi,
      sei               = sei,
      is_weightfunction = is_weightfunction,
      selection_context = selection_context
    ))
  }

  if (outcome_type %in% c("bin", "pois")) {
    return(.log_lik_cluster_glmm_sum(
      setup        = cluster_setup,
      data         = setup[["data"]],
      priors       = setup[["priors"]],
      outcome_type = outcome_type
    ))
  }

  stop("Unsupported outcome type for cluster-unit log-likelihood.",
       call. = FALSE)
}



# ---------------------------------------------------------------------------- #
# .log_lik_cluster.brma
# ---------------------------------------------------------------------------- #
#
# Cluster-unit likelihood for multilevel models.
#
# Each column is the joint held-out-cluster log-likelihood contribution.
#
# @param object brma object.
#
# @return S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster.brma <- function(object) {

  return(.log_lik_from_posterior_samples(
    fit                  = object[["fit"]],
    posterior_samples    = .get_posterior_samples(object[["fit"]]),
    data                 = object[["data"]],
    priors               = object[["priors"]],
    unit                 = "cluster",
    add_metadata         = TRUE,
    data_hash            = .get_outcome_hash(object)
  ))
}



# ---------------------------------------------------------------------------- #
# .add_cluster_log_lik_metadata
# ---------------------------------------------------------------------------- #
#
# @param log_lik         S x G log-likelihood matrix.
# @param cluster_indices named list of cluster index vectors.
# @param data_hash       character; hash of the outcome target.
#
# @return log-likelihood matrix with names and metadata.
#
# ---------------------------------------------------------------------------- #
.add_cluster_log_lik_metadata <- function(log_lik, cluster_indices, data_hash) {

  cluster_labels    <- names(cluster_indices)
  colnames(log_lik) <- paste0("log_lik_cluster[", cluster_labels, "]")
  attr(log_lik, "RoBMA_target") <- list(
    unit             = "cluster",
    retained_context = "remaining_data",
    target           = "cluster_joint",
    n                = length(cluster_indices),
    targets          = cluster_labels,
    cluster_partition = unname(cluster_indices),
    data_hash        = data_hash
  )

  return(log_lik)
}
