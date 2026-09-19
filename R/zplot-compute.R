# ============================================================================ #
# Internal Helper Functions
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# .zplot_fun.brma
# ---------------------------------------------------------------------------- #
#
# Core computation function for zplot estimates and densities.
#
# This function operates in two modes:
# 1. EDR mode (z_threshold set): Computes Expected Discovery Rate
# 2. Density mode (z_sequence set): Computes density over z-values
#
# The extrapolate parameter controls bias adjustment:
# - extrapolate=TRUE:  Removes publication bias (PET/PEESE/weights) to estimate
#                      "true" power distribution
# - extrapolate=FALSE: Includes all bias adjustments for fitted distribution
#
# @param object      brma object with fit and priors
# @param z_threshold z-value threshold for significance (EDR mode)
# @param z_sequence  vector of z-values for density evaluation (density mode)
# @param max_samples maximum posterior samples for computation
# @param extrapolate whether to remove bias adjustments
#
# @return EDR mode:     list with EDR (S-vector) and weights (S-vector)
#         Density mode: S x length(z_sequence) matrix of densities
#
.zplot_fun.brma <- function(object, z_threshold = NULL, z_sequence = NULL,
                             max_samples = 10000, extrapolate = FALSE,
                             conditioning_depth = "marginal",
                             integration_control = set_selection_likelihood_control(),
                             parallel = FALSE, cores = min(4, RoBMA.get_option("max_cores"))) {

  max_samples <- .normalize_max_samples(max_samples, "max_samples")
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)

  ### extract model info
  is_weightfunction <- .is_weightfunction(object)
  effect_direction  <- .effect_direction(object)

  ### 1. Thin posterior samples before derived quantities
  posterior_samples <- .get_posterior_samples(object[["fit"]])
  selected_ind      <- .thin_sample_rows(nrow(posterior_samples), max_samples)
  if (!is.null(selected_ind)) {
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }

  if (.zplot_requires_selection_marginal(object, conditioning_depth)) {
    result <- if (parallel && !extrapolate && is.null(z_threshold) &&
        inherits(object, "brma.mv") && identical(conditioning_depth, "marginal") &&
        cores > 1L && nrow(posterior_samples) >= cores) {
      .zplot_selection_marginal_parallel(object, posterior_samples, z_sequence,
        conditioning_depth, integration_control, cores)
    } else {
      # Serial worker rows inherit the fitted model's parallel setup; the
      # PSOCK branch above disables native threading in its workers.
      .native_threads_configure(.resolve_native_threads(object))
      on.exit(.native_threads_configure(1L), add = TRUE)
      # Only the fitted curve leaves this call when no extrapolated curve and
      # no threshold summary is requested, so the selection routes can skip the
      # inverse weights and the extrapolated mixture they would discard.
      .zplot_selection_marginal(object, posterior_samples, z_sequence, z_threshold,
        conditioning_depth, integration_control, extrapolate_only = extrapolate,
        fitted_only = !extrapolate && is.null(z_threshold))
    }
    if (!is.null(z_threshold)) {
      return(list(
        EDR     = if (extrapolate) result$EDR else result$fitted[, 1L],
        weights = result$weights
      ))
    }
    return(result[[if (extrapolate) "extrapolated" else "fitted"]])
  }

  # The vectorized route's kernels are row-parallel too, so it takes the same
  # thread budget as the selection-marginal route above and restores one thread
  # on exit.
  .native_threads_configure(.resolve_native_threads(object))
  on.exit(.native_threads_configure(1L), add = TRUE)

  predictive <- .zplot_predictive_components(
    object             = object,
    posterior_samples  = posterior_samples,
    extrapolate        = extrapolate,
    conditioning_depth = conditioning_depth
  )
  mu_samples <- predictive[["mu"]]
  tau_within <- predictive[["tau_within"]]
  sei        <- predictive[["sei"]]

  ### 2. Prepare for Computation
  selection <- .zplot_selection_context(
    object            = object,
    posterior_samples = posterior_samples,
    is_weightfunction = is_weightfunction
  )


  ### 3. Dispatch: EDR (Threshold) vs Density (Sequence)

  if (!is.null(z_threshold)) {
    return(.zplot_threshold_vectorized(
      z_threshold      = z_threshold,
      mu_samples       = mu_samples,
      tau_within       = tau_within,
      sei              = sei,
      selection        = selection,
      extrapolate      = extrapolate,
      effect_direction = effect_direction
    ))
  }

  if (!is.null(z_sequence)) {
    return(.zplot_density_vectorized(
      z_sequence       = z_sequence,
      mu_samples       = mu_samples,
      tau_within       = tau_within,
      sei              = sei,
      selection        = selection,
      extrapolate      = extrapolate,
      effect_direction = effect_direction
    ))
  }

  return(NULL)
}


.zplot_predictive_components <- function(object, posterior_samples,
                                         extrapolate,
                                         conditioning_depth = "marginal",
                                         predictive_heterogeneity = NULL) {

  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
  .check_unit_conditioning_depth(
    object             = object,
    unit               = "estimate",
    conditioning_depth = conditioning_depth,
    caller             = ".zplot_predictive_components()"
  )

  predict_type <- switch(conditioning_depth,
    "marginal" = "terms",
    "cluster"  = "cluster",
    "estimate" = "blup"
  )
  mu_samples <- predict.brma(
    object             = object,
    newdata            = NULL,
    type               = predict_type,
    bias_adjusted      = extrapolate,
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  mu_samples <- as.matrix(mu_samples)

  if (is.null(predictive_heterogeneity)) {
    predictive_heterogeneity <- .zplot_predictive_heterogeneity(
      object             = object,
      posterior_samples  = posterior_samples,
      conditioning_depth = conditioning_depth
    )
  }

  return(list(
    mu         = mu_samples,
    tau_within = predictive_heterogeneity,
    sei        = .outcome_data_sei(object)
  ))
}


.zplot_requires_selection_marginal <- function(object, conditioning_depth) {

  .is_weightfunction(object) &&
    (!identical(conditioning_depth, "marginal") ||
       .is_multilevel(object) || inherits(object, "brma.mv") ||
       .zplot_vector_selection_target(object))
}


.zplot_selection_marginal <- function(
    object, posterior_samples, z_sequence, z_threshold,
    conditioning_depth, integration_control, extrapolate_only = FALSE,
    fitted_only = FALSE) {

  control   <- .check_selection_likelihood_control(integration_control)
  selection <- .selection_context(object, posterior_samples = posterior_samples)
  .selection_require_step_evaluable(selection, "zplot()")
  if (conditioning_depth != "marginal") {
    stop(
      "Conditional zplot selection targets are unavailable. ",
      "Use 'conditioning_depth = \"marginal\"'.", call. = FALSE
    )
  }
  predictive <- .zplot_predictive_components(object, posterior_samples, FALSE)
  predictive$mu_extrapolated <- if (.is_PET(object) || .is_PEESE(object)) {
    .zplot_predictive_components(object, posterior_samples, TRUE,
      predictive_heterogeneity = predictive$tau_within)$mu
  } else {
    predictive$mu
  }
  probability <- !is.null(z_threshold)
  z <- if (probability) z_threshold else z_sequence

  if (.is_data_joint_selection(object[["data"]])) {
    return(.zplot_joint_marginal(
      object, posterior_samples, predictive, selection, z, probability, control,
      extrapolate_only, fitted_only
    ))
  }

  K <- length(predictive$sei)
  if (.is_data_random(object[["data"]]) || .is_data_known_v(object[["data"]])) {
    known_V <- if (.is_data_known_v(object[["data"]])) {
      .data_known_v_data(object[["data"]])
    } else {
      .known_v_newdata_prepare(diag(predictive$sei^2, K), K)
    }
    plan <- .known_v_marginal_factor_plan(object, posterior_samples, known_V)
    if (any(plan$sampling_covariance[row(plan$sampling_covariance) !=
                                    col(plan$sampling_covariance)] != 0)) {
      stop("The conditional selection zplot target is unavailable without a conditional sampling factorization.",
           call. = FALSE)
    }
    sd <- sqrt(sweep(plan$extra_variances, 2L, diag(plan$sampling_covariance), "+"))
    latent_sd <- sqrt(plan$latent_variances)
  } else {
    scales    <- .zplot_tau_samples(object, posterior_samples)
    sd        <- .zplot_total_sd(scales$tau_within, predictive$sei)
    latent_sd <- scales$tau_between
  }
  .zplot_latent_mixture(
    z, predictive$mu, sd, latent_sd, predictive$sei, selection,
    probability, control, predictive$mu_extrapolated, fitted_only
  )
}


.zplot_vector_selection_target <- function(object) {

  data <- object[["data"]]
  model <- .data_selection_model(data)
  if (is.null(model)) return(FALSE)
  branches <- model[["branches"]][model[["active_branches"]]]
  # The sampling axis also applies to ordinary independent errors. Integrating
  # those errors retains the established univariate extrapolation convention.
  dependent_sampling <- .is_data_known_v(data) &&
    length(.known_v_correlated_blocks(.data_known_v_data(data))) > 0L
  integrated_context <-
    (isTRUE(model[["applicability"]][["other_random_effects"]]) &&
       identical(model[["other_random_effects"]], "integrate")) ||
    (dependent_sampling && identical(model[["known_sampling_variance"]], "integrate"))
  inherits(object, "brma.mv") || .selection_retains_estimate(data) ||
    .selection_retains_sampling(data) ||
    integrated_context || any(vapply(branches, function(branch) {
    identical(branch[["weight_rule"]], "best")
  }, logical(1)))
}


# Place one chunk's rows in the assembled marginal. A chunk carries only the
# quantities its own request produced: a fitted-only chunk has no extrapolated
# curve, no inverse weights and no EDR, and the assembled result then has none
# either.
.zplot_marginal_chunk_assign <- function(result, chunk, rows, probability) {

  for (name in c("fitted", "extrapolated")) {
    if (is.null(chunk[[name]])) result[[name]] <- NULL else
      result[[name]][rows, ] <- chunk[[name]]
  }
  for (name in c("weights", if (probability) "EDR")) {
    if (is.null(chunk[[name]])) result[[name]] <- NULL else
      result[[name]][rows] <- chunk[[name]]
  }

  result
}


.zplot_joint_marginal <- function(
    object, posterior_samples, predictive, selection, z, probability, control,
    extrapolate_only = FALSE, fitted_only = FALSE) {

  S <- nrow(posterior_samples)
  K <- length(predictive[["sei"]])
  data <- object[["data"]]
  vector_target <- .zplot_vector_selection_target(object)
  if (vector_target && extrapolate_only) {
    gaussian <- .zplot_gaussian_marginal_reference(object, posterior_samples, predictive)
    sd <- sqrt(gaussian[["variance"]])
    if (probability) {
      q <- matrix(z * gaussian[["sei"]], S, K, byrow = TRUE)
      reference <- matrix(rowMeans(stats::pnorm(q, gaussian[["mu"]], sd, lower.tail = FALSE) +
        stats::pnorm(-q, gaussian[["mu"]], sd)), S, 1L)
    } else {
      reference <- .zplot_normal_density_matrix(z, gaussian[["mu"]], sd, gaussian[["sei"]])
    }
    return(list(fitted = NULL, extrapolated = reference, weights = rep(1, S),
      EDR = if (probability) reference[, 1L] else NULL))
  }
  if (vector_target) {
    diagonal <- .zplot_diagonal_product_marginal(
      object, posterior_samples, predictive, selection, z, probability, control
    )
    if (!is.null(diagonal)) return(diagonal)
    streamed <- .zplot_joint_factor_marginal(
      object, posterior_samples, predictive, selection, z, probability, control
    )
    if (!is.null(streamed)) return(streamed)
  }
  chunks <- .known_v_covariance_chunk_indices(
    S = S, K = K, max_bytes = .known_v_covariance_max_bytes() / 4
  )
  if (length(chunks) > 1L) {
    result <- list(fitted = matrix(0, S, length(z)),
      extrapolated = matrix(0, S, length(z)), weights = numeric(S),
      EDR = if (probability) numeric(S) else NULL)
    for (rows in chunks) {
      chunk_predictive <- predictive
      for (name in c("mu", "mu_extrapolated", "tau_within")) {
        if (!is.null(predictive[[name]])) {
          chunk_predictive[[name]] <- predictive[[name]][rows, , drop = FALSE]
        }
      }
      chunk <- .zplot_joint_marginal(object, posterior_samples[rows, , drop = FALSE],
        chunk_predictive, BayesTools::selection_context_subset_rows(selection, rows),
        z, probability, control, extrapolate_only, fitted_only)
      result <- .zplot_marginal_chunk_assign(result, chunk, rows, probability)
    }
    return(result)
  }
  scales <- if (.is_data_random(data)) {
    list(tau_within = matrix(0, S, K), tau_between = matrix(0, S, K))
  } else {
    .zplot_tau_samples(object, posterior_samples)
  }
  parts <- .predict_joint_selection_gaussian_parts(
    object = object, data = data, posterior_samples = posterior_samples,
    fixed_mu = predictive[["mu"]], within = scales[["tau_within"]],
    between = scales[["tau_between"]], draw_context = FALSE
  )
  variance <- .block_covariance_diag_matrix(parts[["covariance"]])
  context_variance <- .block_covariance_diag_matrix(parts[["context_covariance"]])
  normal <- function(mean) {
    sd <- sqrt(variance + context_variance)
    if (!probability) return(.zplot_normal_density_matrix(z, mean, sd, predictive[["sei"]]))
    q <- matrix(z * predictive[["sei"]], S, K, byrow = TRUE)
    matrix(rowMeans(stats::pnorm(q, mean, sd, lower.tail = FALSE) +
      stats::pnorm(-q, mean, sd)), S, 1L)
  }
  reference <- normal(predictive[["mu_extrapolated"]])
  result <- list(fitted = NULL, extrapolated = reference, weights = rep(1, S),
    EDR = if (probability) reference[, 1L] else NULL)

  diagonal_kernel <- .block_covariance_diagonal_draws(parts[["covariance"]])
  zero_context <- .block_covariance_zero_draws(parts[["context_covariance"]])
  rules <- rep_len(selection[["vector_rule"]], S)
  if (!vector_target && all(rules == 0L) && all(diagonal_kernel)) {
    return(.zplot_latent_mixture(
      z, predictive[["mu"]], sqrt(variance), sqrt(context_variance),
      predictive[["sei"]], selection, probability, control,
      predictive[["mu_extrapolated"]], fitted_only
    ))
  }
  if (!vector_target && all(rules == 0L) && all(zero_context)) {
    return(.zplot_integrated_marginal(object, posterior_samples, predictive,
      selection, z, probability, control, extrapolate_only, fitted_only))
  }

  result[["fitted"]] <- normal(predictive[["mu"]])
  active <- which(!selection[["use_normal"]])
  # With no refreshed source, positive selection weights cancel at each
  # retained realization. Its marginal projection is the ordinary Gaussian.
  deterministic <- .block_covariance_zero_draws(parts[["covariance"]])
  active <- setdiff(active, which(deterministic))
  scalar <- active[rules[active] == 0L & diagonal_kernel[active] &
    rowSums(variance[active, , drop = FALSE] <= 0) == 0L]
  if (length(scalar)) {
    result[["fitted"]][scalar, ] <- .zplot_latent_mixture(
      z, predictive[["mu"]][scalar, , drop = FALSE],
      sqrt(variance[scalar, , drop = FALSE]), sqrt(context_variance[scalar, , drop = FALSE]),
      predictive[["sei"]], BayesTools::selection_context_subset_rows(selection, scalar),
      probability, control, fitted_only = TRUE
    )[["fitted"]]
  }
  active <- setdiff(active, scalar)
  integrated <- active[rules[active] == 0L & zero_context[active]]
  if (length(integrated)) {
    subset_predictive <- predictive
    for (name in c("mu", "mu_extrapolated", "tau_within")) {
      subset_predictive[[name]] <- predictive[[name]][integrated, , drop = FALSE]
    }
    result[["fitted"]][integrated, ] <- .zplot_integrated_marginal(
      object, posterior_samples[integrated, , drop = FALSE], subset_predictive,
      BayesTools::selection_context_subset_rows(selection, integrated),
      z, probability, control, fitted_only = TRUE
    )[["fitted"]]
  }
  active <- setdiff(active, integrated)
  if (length(active)) {
    execution_plan <- .selection_joint_execution_plan_with_control(
      .data_selection_execution_plan(data), control
    )
    setup <- list(fit = object[["fit"]], data = data, priors = object[["priors"]],
      posterior_samples = posterior_samples[active, , drop = FALSE], S = length(active), K = K,
      tau_within = scales[["tau_within"]][active, , drop = FALSE],
      tau_between = scales[["tau_between"]][active, , drop = FALSE],
      is_multilevel = .is_multilevel(object))
    if (.selection_retains_sampling(data)) {
      random_factors <- .selection_conditioned_sampling_factors(setup)
      block_factors <- lapply(seq_along(execution_plan[["row_blocks"]]), function(block) {
        rows <- execution_plan[["row_blocks"]][[block]]
        loading <- random_factors[["loadings"]][[block]]
        list(diagonal = random_factors[["diagonal"]][, rows, drop = FALSE],
             residual_sd = sqrt(random_factors[["diagonal"]][, rows, drop = FALSE]),
             loading = matrix(loading, nrow = length(active)), rank = dim(loading)[[3L]])
      })
    } else {
      random_factors <- .selection_joint_random_factor_samples(setup)
      block_factors <- lapply(seq_along(execution_plan[["row_blocks"]]), function(block) {
        if (!execution_plan[["block_methods"]][[block]] %in% c("factor", "rank_one")) return(NULL)
        .selection_joint_factor_block_samples(setup, block, random_factors)
      })
    }
    result[["fitted"]][active, ] <- .zplot_full_event_context_mixture(
      z = z, mean = parts[["means"]][active, , drop = FALSE],
      covariance = .block_covariance_draws(parts[["covariance"]], active),
      context_covariance = .block_covariance_draws(
        parts[["context_covariance"]], active),
      sei = predictive[["sei"]],
      selection = BayesTools::selection_context_subset_rows(selection, active),
      probability = probability, control = control, execution_plan = execution_plan,
      block_factors = block_factors,
      sampling_factor_blocks = execution_plan[["sampling_factor_blocks"]],
      random_covariance = .block_covariance_draws(
        parts[["random_covariance"]], active),
      publication_groups = .data_selection_model(data)[["groups"]][["group_index"]]
    )
  }
  result
}


# For diagonal integrated covariance and product weights, conditional selection
# factorizes by outcome. Each displayed marginal therefore needs only its own
# retained Gaussian variance, even when retained sources correlate outcomes.
# Resolve diagonality from the compiled source graph before constructing cubes.
.zplot_diagonal_product_marginal <- function(
    object, posterior_samples, predictive, selection, z, probability, control) {

  if (!inherits(object, "brma.mv") || is.null(selection) ||
      any(selection[["vector_rule"]] != 0L)) return(NULL)
  S <- nrow(posterior_samples)
  K <- length(predictive[["sei"]])
  data <- object[["data"]]
  model <- .data_selection_model(data)
  if (is.null(model)) return(NULL)
  retained_sampling <- .selection_retains_sampling(data)
  known_V <- .data_known_v_data(data)
  if (!retained_sampling && !is.null(known_V) &&
      length(.known_v_correlated_blocks(known_V)) > 0L) return(NULL)
  variance <- context_variance <- matrix(0, S, K)
  if (.is_data_random(data)) {
    sources <- model[["sources"]][["random"]]
    retained <- vapply(sources, function(source) isTRUE(source[["retained"]]), logical(1L))
    source_names <- vapply(sources, `[[`, character(1L), "name")
    integrated_names <- source_names[!retained]
    if (length(integrated_names)) {
      design <- .fitted_formula_design(object, "mu", required = TRUE)
      dependency <- BayesTools::random_effects_dependency_matrix(
        random_effects = design[["random_effects"]], n_rows = K, blocks = integrated_names
      )
      if (any(dependency[row(dependency) != col(dependency)])) return(NULL)
    }
    components <- .brma_mv_heterogeneity_components(object, posterior_samples,
      include_known_group_covariance = TRUE)
    if (!all(source_names %in% names(components))) {
      stop("Diagonal selection source variances are inconsistent.", call. = FALSE)
    }
    for (index in seq_along(sources)) {
      sd <- .expand_brma_mv_heterogeneity_samples(components[[source_names[index]]], S, K)
      if (any(!is.finite(sd)) || any(sd < 0)) {
        stop("Diagonal selection source standard deviations are invalid.", call. = FALSE)
      }
      if (retained[index]) context_variance <- context_variance + sd^2 else
        variance <- variance + sd^2
    }
  }
  sampling_variance <- if (is.null(known_V)) predictive[["sei"]]^2 else .known_v_diagonal(known_V)
  if (retained_sampling) context_variance <- sweep(context_variance, 2L, sampling_variance, "+") else
    variance <- sweep(variance, 2L, sampling_variance, "+")
  total <- variance + context_variance
  if (any(!is.finite(total)) || any(total <= 0)) return(NULL)
  active <- which(!selection[["use_normal"]])
  positive <- variance > 0
  zero_rows <- active[rowSums(!positive[active, , drop = FALSE]) > 0L]
  # A zero candidate coordinate cancels only where its retained event has
  # positive weight. Preserve the existing full-event handling otherwise.
  if (length(zero_rows) && any(selection[["omega"]][zero_rows, , drop = FALSE] <= 0)) return(NULL)
  normal <- function(mean, rows = seq_len(S), columns = seq_len(K)) {
    sd <- sqrt(total[rows, columns, drop = FALSE])
    sei <- predictive[["sei"]][columns]
    mean <- mean[rows, columns, drop = FALSE]
    if (!probability) return(.zplot_normal_density_matrix(z, mean, sd, sei))
    threshold <- matrix(z * sei, length(rows), length(columns), byrow = TRUE)
    matrix(rowMeans(stats::pnorm(threshold, mean, sd, lower.tail = FALSE) +
      stats::pnorm(-threshold, mean, sd)), length(rows), 1L)
  }
  reference <- normal(predictive[["mu_extrapolated"]])
  fitted <- if (identical(predictive[["mu"]], predictive[["mu_extrapolated"]])) reference else
    normal(predictive[["mu"]])
  if (length(active)) {
    groups <- split(active, apply(positive[active, , drop = FALSE], 1L, paste0, collapse = ""))
    for (rows in groups) {
      columns <- which(positive[rows[1L], ])
      value <- matrix(0, length(rows), length(z))
      if (length(columns)) {
        context <- BayesTools::selection_context_subset_observations(
          BayesTools::selection_context_subset_rows(selection, rows), columns)
        value <- .zplot_latent_mixture(z, predictive[["mu"]][rows, columns, drop = FALSE],
          sqrt(variance[rows, columns, drop = FALSE]),
          sqrt(context_variance[rows, columns, drop = FALSE]),
          predictive[["sei"]][columns], context, probability, control,
          fitted_only = TRUE)[["fitted"]] * (length(columns) / K)
      }
      cancelled <- which(!positive[rows[1L], ])
      if (length(cancelled)) {
        value <- value + normal(predictive[["mu"]], rows, cancelled) * (length(cancelled) / K)
      }
      fitted[rows, ] <- value
    }
  }
  list(fitted = fitted, extrapolated = reference, weights = rep(1, S),
    EDR = if (probability) reference[, 1L] else NULL)
}


# Average conditional selected projections over newly realized retained
# contexts. Source roles and normalization units have already been compiled.
.zplot_full_event_context_mixture <- function(
    z, mean, covariance, context_covariance, sei, selection, probability,
    control, execution_plan, block_factors = NULL, sampling_factor_blocks = NULL,
    random_covariance = NULL, publication_groups = NULL) {

  S <- nrow(mean)
  K <- ncol(mean)
  output <- matrix(0, S, length(z))
  designs <- new.env(parent = emptyenv())
  if (is.null(publication_groups)) publication_groups <- rep.int(1L, K)
  for (block in seq_along(execution_plan[["row_blocks"]])) {
    observations <- execution_plan[["row_blocks"]][[block]]
    k <- length(observations)
    for (draw in seq_len(S)) {
      sigma <- .block_covariance_sub(covariance, draw, observations)
      latent <- .block_covariance_sub(context_covariance, draw, observations)
      context <- BayesTools::selection_context_subset_observations(
        BayesTools::selection_context_subset_rows(selection, draw), observations)
      rank_one <- if (!is.null(random_covariance)) {
        .selection_joint_declared_rank_one_loading(sampling_factor_blocks[[block]],
          .block_covariance_sub(random_covariance, draw, observations))
      } else NULL
      factors <- block_factors[[block]]
      if (!is.null(factors)) {
        for (name in c("residual_sd", "loading")) {
          if (!is.null(factors[[name]])) factors[[name]] <- factors[[name]][draw, , drop = FALSE]
        }
        if (identical(execution_plan[["statistical_target"]], "whole_sampling_error_selection")) {
          factors$diagonal <- factors$diagonal[draw, , drop = FALSE]
        }
      }
      current <- .zplot_full_event_context_draw(z, mean[draw, observations], sigma, latent,
        sei[observations], context, probability, control, execution_plan,
        factors = factors, rank_one = rank_one, publication_groups = publication_groups[observations],
        designs = designs)
      output[draw, ] <- output[draw, ] + current[1L, ] * k / K
    }
  }
  output
}


# Evaluate one original full-event context cell. Both covariance-cube and
# streamed callers use this same fallback and its unchanged diagnostics.
.zplot_full_event_context_draw <- function(
    z, mean, sigma, latent, sei, context, probability, control, execution_plan,
    factors = NULL, rank_one = NULL, publication_groups = rep.int(1L, length(sei)),
    designs = new.env(parent = emptyenv())) {

  k <- length(sei)
  factor_projection <- !is.null(factors) && identical(
    execution_plan[["statistical_target"]], "whole_sampling_error_selection")
  zero_variance <- if (factor_projection) which(diag(sigma) == 0) else integer()
  unchanged <- numeric(length(z))
  if (length(zero_variance)) {
    sd <- sqrt(diag(latent)[zero_variance])
    if (any(sd <= 0)) {
      stop("Zplot marginal density is unavailable because an outcome has no continuous Gaussian variation.",
           call. = FALSE)
    }
    unchanged <- vapply(z, function(value) {
      q <- value * sei[zero_variance]
      mu <- mean[zero_variance]
      if (probability) {
        sum(stats::pnorm(-q, mu, sd) + stats::pnorm(q, mu, sd, lower.tail = FALSE)) / k
      } else {
        sum(sei[zero_variance] * stats::dnorm(q, mu, sd)) / k
      }
    }, numeric(1L))
  }
  if (!is.null(factors) && identical(factors[["rank"]], 1L) &&
      all(factors[["residual_sd"]][1L, ] == 0) &&
      all(factors[["loading"]][1L, ] != 0)) {
    rank_one <- as.numeric(factors[["loading"]][1L, ])
  }
  project_direct <- function(means, diagnostics = FALSE) {

    local_factors <- factors
    if (!is.null(local_factors)) {
      local_factors$residual_sd <- local_factors$residual_sd[rep.int(1L, nrow(means)), , drop = FALSE]
      local_factors$loading <- local_factors$loading[rep.int(1L, nrow(means)), , drop = FALSE]
      if (factor_projection) {
        local_factors$diagonal <- local_factors$diagonal[rep.int(1L, nrow(means)), , drop = FALSE]
      }
    }
    context_rows <- BayesTools::selection_context_subset_rows(context, rep.int(1L, nrow(means)))
    if (factor_projection) {
      dimensions <- max(2L * k, local_factors[["rank"]] + k)
      key <- paste("projection", dimensions, control[["max_points_per_scramble"]], sep = "/")
      if (!exists(key, designs, inherits = FALSE)) {
        assign(key, BayesTools::selection_qmc_design(
          dimensions = dimensions, points = control[["max_points_per_scramble"]],
          scrambles = control[["scrambles"]], seed = control[["seed"]]), designs)
      }
      return(.selection_factor_projection(
        means, local_factors[["diagonal"]], local_factors[["loading"]], sei,
        context_rows, execution_plan, z, probability, publication_groups,
        get(key, designs))[["density"]])
    }
    packed <- matrix(sigma[lower.tri(sigma, diag = TRUE)],
      nrow(means), k * (k + 1L) / 2L, byrow = TRUE)
    if (all(context_rows[["vector_rule"]] == 0L) && is.null(rank_one)) {
      projected <- .zplot_joint_block(z, means, packed, sei,
        context_rows, probability, control, designs, local_factors)
      density <- projected[["density"]]
      if (diagnostics) {
        relative_error <- projected[["relative_mcse"]]
        if (!is.numeric(relative_error) || length(relative_error) != nrow(means) ||
            any(!is.finite(relative_error)) || any(relative_error < 0)) {
          stop("Zplot fallback integration diagnostics are unavailable.", call. = FALSE)
        }
        # This is the generic kernel's checked MCSE/refinement estimate,
        # including its normalizer check, rather than a certified bound.
        attr(density, "integration_error") <- list(
          absolute = relative_error * apply(density, 1L, max), mass = relative_error)
      }
      return(density)
    }
    .zplot_full_event_projection(
      z, means, packed[1L, , drop = FALSE], sei, context_rows, probability,
      execution_plan, rank_one_loading = if (!is.null(rank_one))
        matrix(rank_one, nrow(means), k, byrow = TRUE) else NULL
    )
  }
  point_projection <- !probability && !factor_projection && k > 1L &&
    is.null(rank_one) && context[["vector_rule"]] == 0L &&
    context[["kernel_mode"]] == SELKERNEL_STEP && all(context[["omega"]] > 0)
  project <- function(means, absolute_tolerance = NULL, diagnostics = FALSE,
                      groups = seq_len(nrow(means))) {

    # The outer QMC nodes are fixed retained realizations. Each therefore has
    # a point-context full-event density, including the same original dense
    # sampling covariance. Reuse its checked factor projection before the
    # generic inner QMC calculation, without changing the outer nodes or gates.
    if (!point_projection) {
      return(project_direct(means))
    }
    output <- matrix(0, nrow(means), length(z))
    absolute_error <- mass_error <- numeric(nrow(means))
    pending <- logical(nrow(means))
    point_factor <- matrix(0, 1L, k)
    # Nodes differ only in their context mean, so they are projected together:
    # one certified projection of a batch's equal-weight average replaces one
    # per node, whose setup and native calls dominated the cost. A batch never
    # mixes scrambles, so each scramble's node sum and the spread between
    # scrambles the outer error rests on are unchanged; every node of a batch
    # carries the batch average and its bound.
    batches <- unlist(lapply(split(seq_len(nrow(means)), groups), function(rows) {
      split(rows, ceiling(seq_along(rows) / .zplot_context_batch_size()))
    }), recursive = FALSE, use.names = FALSE)
    for (rows in batches) {
      batch <- means[rows, , drop = FALSE]
      projected <- .zplot_context_projection(z, colMeans(batch), sigma,
        point_factor, sei, context, control, absolute_tolerance,
        context_means = batch)
      if (is.null(projected)) {
        pending[rows] <- TRUE
      } else {
        density <- .zplot_context_projection_density(
          projected, z, control, absolute_tolerance)
        output[rows, ] <- matrix(density[1L, ], length(rows), length(z), byrow = TRUE)
        absolute_error[rows] <- projected[["integration_error"]][["absolute"]]
        mass_error[rows] <- projected[["mass_error"]]
      }
    }
    pending <- which(pending)
    if (length(pending)) {
      direct <- project_direct(means[pending, , drop = FALSE],
        diagnostics = diagnostics)
      output[pending, ] <- direct
      if (diagnostics) {
        errors <- attr(direct, "integration_error", exact = TRUE)
        absolute_error[pending] <- errors[["absolute"]]
        mass_error[pending] <- errors[["mass"]]
      }
    }
    if (diagnostics) {
      attr(output, "integration_error") <- list(absolute = absolute_error, mass = mass_error)
    }
    output
  }
  if (length(zero_variance) == k) {
    current <- matrix(0, 1L, length(z))
  } else if (all(latent == 0)) {
    current <- project(matrix(mean, 1L, k))
  } else {
    factor <- .covariance_sampling_factor(.covariance_factorization(latent))
    if (is.null(factor)) {
      stop("Zplot retained-context covariance must be positive semidefinite.", call. = FALSE)
    }
    points <- control[["points_per_scramble"]]
    scrambles <- control[["scrambles"]]
    point_absolute_tolerance <- NULL
    if (point_projection) {
      reference <- .zplot_normal_density_matrix(z, matrix(mean, 1L, k),
        matrix(sqrt(diag(sigma) + diag(latent)), 1L, k), sei)
      # This allocates inner work only. Acceptance below uses the combined
      # errors and the computed selected curve, not this Gaussian reference.
      allocation <- control[["relative_tolerance"]] * max(reference) / 16
      if (is.finite(allocation) && allocation > 0) point_absolute_tolerance <- allocation
    }
    used <- 0L
    sums <- matrix(0, scrambles, length(z))
    error_sums <- matrix(0, scrambles, 2L)
    previous <- NULL
    previous_absolute_error <- 0
    repeat {
      key <- paste("context", k, points, sep = "/")
      if (!exists(key, designs, inherits = FALSE)) {
        assign(key, BayesTools::selection_qmc_design(
          dimensions = k, points = points, scrambles = scrambles,
          seed = control[["seed"]]), designs)
      }
      uniforms <- get(key, designs)[, seq.int(used + 1L, points), , drop = FALSE]
      contexts <- matrix(stats::qnorm(uniforms), (points - used) * scrambles, k) %*% factor
      means <- sweep(contexts, 2L, mean, "+")
      values <- project(means, point_absolute_tolerance, diagnostics = point_projection,
        groups = rep(seq_len(scrambles), points - used))
      sums <- sums + rowsum(values, rep(seq_len(scrambles), points - used), reorder = FALSE)
      numerical <- attr(values, "integration_error", exact = TRUE)
      if (!is.null(numerical)) {
        error_sums <- error_sums + rowsum(cbind(numerical[["absolute"]], numerical[["mass"]]),
          rep(seq_len(scrambles), points - used), reorder = FALSE)
      }
      estimates <- sums / points
      current <- matrix(colMeans(estimates), 1L, length(z))
      peak <- max(current)
      # Positive averaging carries the inner bounds and checked error estimates
      # through the QMC weights. Their perturbation of the scramble means can
      # also affect MCSE; use the corresponding Euclidean norm allowance.
      scramble_errors <- error_sums[, 1L] / points
      numerical_error <- mean(scramble_errors)
      mass_error <- mean(error_sums[, 2L] / points)
      mcse_error <- sqrt(sum(scramble_errors^2) / (scrambles * (scrambles - 1L)))
      error <- max(apply(estimates, 2L, stats::sd)) / sqrt(scrambles)
      error <- error + mcse_error
      if (!is.null(previous)) {
        change <- max(abs(current - previous)) + numerical_error + previous_absolute_error
        error <- max(error, change)
      }
      error <- error + numerical_error
      peak_lower <- peak - numerical_error
      if (peak_lower > 0) error <- error / peak_lower else if (numerical_error > 0) error <- Inf
      if (all(is.finite(current)) && is.finite(error) &&
          error <= control[["relative_tolerance"]] && mass_error <= control[["relative_tolerance"]]) break
      if (points >= control[["max_points_per_scramble"]]) {
        if (mass_error > control[["relative_tolerance"]]) {
          stop("Zplot retained-context integration was rejected by diagnostics: normalization error was ",
            format(mass_error, digits = 4), ". Inspect the inner selection integration diagnostics.", call. = FALSE)
        }
        stop("Zplot retained-context integration was rejected by diagnostics: relative integration error was ",
          format(error, digits = 4), ". Increase 'max_points_per_scramble' in ",
          "'integration_control = set_selection_likelihood_control()'.", call. = FALSE)
      }
      previous <- current
      previous_absolute_error <- numerical_error
      used <- points
      points <- min(2L * points, control[["max_points_per_scramble"]])
    }
  }
  matrix(current[1L, ] + unchanged, 1L, length(z))
}


# Keep the original event in both the denominator and every partial numerator.
.zplot_full_event_projection <- function(
    z, mean, covariance_lower, sei, selection, probability, execution_plan,
    rank_one_loading = NULL, diagonal_best = TRUE) {

  S <- nrow(mean)
  K <- ncol(mean)
  remedy <- paste0("Increase 'max_points_per_scramble' in ",
    "'integration_control = set_selection_likelihood_control()'.")
  evaluate <- function(observed = integer(), values = numeric(), lower = NULL, upper = NULL) {

    .selection_joint_checked_event(function(plan, rows) {
      if (is.null(rows)) rows <- seq_len(S)
      context <- BayesTools::selection_context_subset_rows(selection, rows)
      packed <- if (is.null(rank_one_loading)) {
        if (nrow(covariance_lower) == 1L) covariance_lower else
          covariance_lower[rows, , drop = FALSE]
      } else NULL
      loading <- if (!is.null(rank_one_loading)) rank_one_loading[rows, , drop = FALSE] else NULL
      if (!length(observed)) {
        if (!is.null(packed) && nrow(packed) == 1L) {
          packed <- packed[rep.int(1L, length(rows)), , drop = FALSE]
        }
        return(.selection_gaussian_event_mass(mean[rows, , drop = FALSE],
          packed, sei, context, plan, lower = lower, upper = upper,
          rank_one_loading = loading))
      }
      .selection_joint_event_numerator(mean[rows, , drop = FALSE],
        packed, sei, context, plan, observed = observed, values = values,
        rank_one_loading = loading)
    }, execution_plan, "Zplot selection event integration", remedy)
  }
  normalizer <- evaluate()[["log_mass"]]
  if (any(!is.finite(normalizer))) {
    stop("Zplot selection normalizers must be finite and positive.", call. = FALSE)
  }
  if (diagonal_best && is.null(rank_one_loading) && nrow(covariance_lower) == 1L &&
      all(selection[["vector_rule"]] != 0L)) {
    pairs <- which(lower.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
    diagonal <- pairs[, 1L] == pairs[, 2L]
    if (all(covariance_lower[1L, !diagonal] == 0) &&
        all(covariance_lower[1L, diagonal] > 0)) {
      fast <- .zplot_diagonal_best_projection(z, mean, covariance_lower[1L, diagonal],
        sei, selection, probability, execution_plan, normalizer)
      if (!is.null(fast)) {
        if (length(fast[["boundary_columns"]])) {
          columns <- fast[["boundary_columns"]]
          fast[["density"]][, columns] <- .zplot_full_event_projection(
            z[columns], mean, covariance_lower, sei, selection, probability,
            execution_plan, diagonal_best = FALSE)
        }
        return(fast[["density"]])
      }
    }
  }
  result <- matrix(0, S, length(z))
  for (observation in seq_len(K)) {
    for (point in seq_along(z)) {
      q <- z[[point]] * sei[[observation]]
      if (probability) {
        lower <- rep(-Inf, K)
        upper <- rep(Inf, K)
        upper[[observation]] <- -q
        left <- evaluate(upper = upper)[["log_mass"]]
        upper[[observation]] <- Inf
        lower[[observation]] <- q
        right <- evaluate(lower = lower)[["log_mass"]]
        value <- exp(left - normalizer) + exp(right - normalizer)
      } else {
        numerator <- evaluate(observation, q)[["log_numerator"]]
        value <- sei[[observation]] * exp(numerator - normalizer)
      }
      result[, point] <- result[, point] + value / K
    }
  }
  result
}


# With independent conditional errors, the other estimates' event weight is
# constant within a focal selection bin. Reuse those weights across the grid.
.zplot_diagonal_best_projection <- function(
    z, mean, variance, sei, selection, probability, execution_plan, log_normalizer) {

  S <- nrow(mean)
  K <- ncol(mean)
  bins <- ncol(selection[["omega"]])
  midpoint <- vapply(seq_len(bins), function(bin) {
    .selection_segment_midpoint(selection[["z_lower"]][[bin]], selection[["z_upper"]][[bin]])
  }, numeric(1L))
  output <- matrix(0, S, length(z))
  boundary_columns <- integer()
  for (observation in seq_len(K)) {
    others <- setdiff(seq_len(K), observation)
    weights <- matrix(0, S, bins)
    for (bin in seq_len(bins)) {
      observed <- matrix(selection[["sign"]] * midpoint[[bin]] * sei[[observation]], 1L, 1L)
      context <- .selection_joint_condition_event_context(selection, observed, sei[[observation]])
      if (!length(others)) {
        weights[, bin] <- exp(.selection_joint_log_weight(
          matrix(observed, S, 1L), sei[[observation]], selection))
        next
      }
      sigma <- diag(variance[others], length(others))
      lower <- matrix(sigma[lower.tri(sigma, diag = TRUE)], S,
        length(others) * (length(others) + 1L) / 2L, byrow = TRUE)
      mass <- .selection_joint_checked_event(function(plan, rows) {
        if (is.null(rows)) rows <- seq_len(S)
        .selection_gaussian_event_mass(mean[rows, others, drop = FALSE],
          lower[rows, , drop = FALSE], sei[others],
          BayesTools::selection_context_subset_rows(context, rows), plan)
      }, execution_plan, "Zplot selection event integration", paste0(
        "Increase 'max_points_per_scramble' in ",
        "'integration_control = set_selection_likelihood_control()'."))
      weights[, bin] <- exp(mass[["log_mass"]])
      if (any(!is.finite(weights[, bin])) ||
          any(weights[, bin] == 0 & is.finite(mass[["log_mass"]]))) return(NULL)
    }
    context <- BayesTools::selection_context_subset_observations(selection, observation)
    context[["omega"]] <- weights
    context[["vector_rule"]] <- integer(S)
    sd <- matrix(sqrt(variance[[observation]]), S, 1L)
    mu <- mean[, observation, drop = FALSE]
    local_normalizer <- .selection_step_log_norm_matrix(mu, sd, sei[[observation]], context)[, 1L]
    correction <- exp(local_normalizer - log_normalizer)
    if (any(!is.finite(correction))) return(NULL)
    value <- if (probability) {
      matrix(.zplot_selnorm_threshold_summary(z, mu, sd, sei[[observation]],
        context, FALSE)[["EDR"]], S, length(z))
    } else {
      .zplot_selnorm_density_matrix(z, mu, sd, sei[[observation]], context, FALSE)
    }
    if (any(!is.finite(value))) return(NULL)
    output <- output + value * correction / K
    if (!probability) {
      native_bins <- .selection_step_bin_from_z(selection[["sign"]] * z, selection[["p_cuts"]])
      for (rule in unique(selection[["vector_rule"]][!selection[["use_normal"]]])) {
        represented <- vapply(native_bins, function(bin) {
          .selection_joint_best_bin(selection[["sign"]] * midpoint[[bin]] * sei[[observation]],
            sei[[observation]], selection, rule)
        }, integer(1L))
        actual <- vapply(z, function(value) {
          .selection_joint_best_bin(value * sei[[observation]], sei[[observation]], selection, rule)
        }, integer(1L))
        boundary_columns <- union(boundary_columns, which(represented != actual))
      }
    }
  }
  list(density = output, boundary_columns = boundary_columns)
}


.selection_joint_declared_rank_one_loading <- function(sampling, added_covariance) {

  if (is.null(sampling) || !identical(sampling[["rank"]], 1L) ||
      any(sampling[["diagonal"]] != 0) || any(added_covariance != 0)) return(NULL)
  as.numeric(sampling[["loading"]][, 1L])
}


# Integrated product-event projection with the legacy extrapolation helper.
# Public vector targets request only the fitted projection from this path.
.zplot_integrated_marginal <- function(
    object, posterior_samples, predictive, selection, z, probability, control,
    extrapolate_only = FALSE, fitted_only = FALSE) {

  S <- nrow(posterior_samples)
  K <- length(predictive$sei)
  scales <- if (.is_data_random(object[["data"]])) {
    list(tau_within = matrix(0, S, K), tau_between = matrix(0, S, K))
  } else {
    .zplot_tau_samples(object, posterior_samples)
  }
  fitted       <- if (extrapolate_only) NULL else matrix(0, S, length(z))
  extrapolated <- matrix(0, S, length(z))
  weights      <- numeric(S)
  designs      <- new.env(parent = emptyenv())
  fitted_plan <- .data_selection_execution_plan(object[["data"]])
  execution_plan <- fitted_plan
  if (extrapolate_only) {
    execution_plan <- .selection_joint_execution_plan_with_control(
      fitted_plan, control
    )
  }
  factor_setup <- list(
    fit = object[["fit"]], data = object[["data"]], priors = object[["priors"]],
    posterior_samples = posterior_samples, S = S, K = K,
    tau_within = scales$tau_within, tau_between = scales$tau_between,
    is_multilevel = .is_multilevel(object)
  )
  random_factors <- .selection_joint_random_factor_samples(factor_setup)
  random_covariance <- if (is.null(random_factors)) {
    .selection_joint_random_covariance_samples(factor_setup)
  } else {
    NULL
  }
  for (block in seq_along(fitted_plan$row_blocks)) {
    rows <- fitted_plan$row_blocks[[block]]
    k <- length(rows)
    factors <- if (fitted_plan$block_methods[block] %in% c("factor", "rank_one")) {
      .selection_joint_factor_block_samples(factor_setup, block, random_factors)
    } else {
      NULL
    }
    if (!extrapolate_only && !is.null(factors) &&
        (factors$rank == 0L || factors$rank > k || any(factors$residual_sd <= 0))) {
      factors <- NULL
    }
    pairs <- .selection_joint_lower_pairs(fitted_plan, seq_len(k))
    covariance <- .selection_joint_covariance_lower(
      factor_setup, block, random_covariance, random_factors
    )
    result <- if (extrapolate_only) {
      context <- BayesTools::selection_context_subset_observations(selection, rows)
      .zplot_normalizer_block(
        .outcome_data_yi(object)[rows], predictive$mu[, rows, drop = FALSE],
        covariance, predictive$sei[rows], context, control,
        execution_plan, block, designs, factors
      )
    } else {
      .zplot_joint_block(
        z, predictive$mu[, rows, drop = FALSE], covariance,
        predictive$sei[rows], selection, probability, control, designs, factors
      )
    }
    if (fitted_only) {
      fitted <- fitted + result$density * (k / K)
      next
    }
    sd <- sqrt(covariance[, pairs$row_1 == pairs$row_2, drop = FALSE])
    normal <- if (probability) {
      matrix(0, S, 1L)
    } else {
      .zplot_normal_density_matrix(z, predictive$mu_extrapolated[, rows, drop = FALSE], sd,
                                   predictive$sei[rows])
    }
    if (probability) {
      q <- matrix(z * predictive$sei[rows], S, k, byrow = TRUE)
      normal[, 1L] <- rowMeans(
        stats::pnorm(q, predictive$mu_extrapolated[, rows, drop = FALSE], sd, lower.tail = FALSE) +
          stats::pnorm(-q, predictive$mu_extrapolated[, rows, drop = FALSE], sd)
      )
    }
    inverse <- exp(result$log_density)
    if (!extrapolate_only) fitted <- fitted + result$density * (k / K)
    extrapolated <- extrapolated + normal * inverse * (k / K)
    weights      <- weights + inverse * (k / K)
  }
  if (fitted_only) return(list(fitted = fitted))
  list(fitted = fitted, extrapolated = extrapolated, weights = weights,
       EDR = if (probability) extrapolated[, 1L] / weights else NULL)
}


# Refine only rows rejected by the existing Gaussian-event MCSE criterion.
# compute(plan, rows) returns row vectors, with NULL selecting the initial batch.
.selection_joint_checked_event <- function(compute, execution_plan, subject, remedy) {

  plan <- execution_plan
  current <- compute(plan, NULL)
  S <- length(current[["relative_mcse"]])
  result <- current
  active <- seq_len(S)
  repeat {
    if (!is.list(current) || !length(current) || !S ||
        any(vapply(current, function(value) {
          !is.numeric(value) || !is.null(dim(value)) || length(value) != length(active)
        }, logical(1)))) {
      stop("Gaussian selection event results are inconsistent.", call. = FALSE)
    }
    for (name in names(current)) result[[name]][active] <- current[[name]]
    quadrature_error <- current[["relative_quadrature_error"]]
    failed_quadrature <- which(!is.finite(quadrature_error) |
      quadrature_error > plan[["relative_tolerance"]])
    if (length(failed_quadrature)) {
      stop(subject, " was rejected by diagnostics: relative quadrature error was ",
        format(quadrature_error[failed_quadrature[[1L]]], digits = 4),
        ". Inspect the integration diagnostics.", call. = FALSE)
    }
    error <- current[["relative_mcse"]]
    failed <- which(!is.finite(error) | error > plan[["relative_tolerance"]])
    if (!length(failed)) return(result)
    if (plan[["points_per_scramble"]] >= plan[["max_points_per_scramble"]]) {
      stop(subject, " was rejected by diagnostics: relative integration MCSE was ",
        format(error[failed[[1L]]], digits = 4), ". ", remedy, call. = FALSE)
    }
    active <- active[failed]
    plan[["points_per_scramble"]] <- min(
      2L * plan[["points_per_scramble"]], plan[["max_points_per_scramble"]]
    )
    design_names <- names(plan[["designs"]])
    if (length(design_names)) {
      plan[["designs"]] <- plan[["designs"]][!grepl("^[0-9]+$", design_names)]
    }
    current <- compute(plan, active)
  }
}


# EDR and extrapolated curves need the joint normalizer, but no selected
# threshold or density projection. Reuse the likelihood's certified routes.
.zplot_normalizer_block <- function(
    yi, mean, covariance_lower, sei, selection, control, execution_plan,
    block_index, designs, factors) {

  S      <- nrow(mean)
  K      <- ncol(mean)
  method <- execution_plan$block_methods[[block_index]]
  active <- seq_len(S)
  points <- control$points_per_scramble
  result <- list(log_density = numeric(S), relative_mcse = numeric(S))
  previous <- NULL
  repeat {
    context <- BayesTools::selection_context_subset_rows(selection, active)
    current <- if (method == "rank_one") {
      .selection_joint_cluster_loglik_block(
        yi, mean, factors$residual_sd, factors$loading, sei, context,
        execution_plan, return_normalizer = TRUE
      )
    } else if (method == "factor") {
      .selection_joint_factor_loglik_block(
        yi, mean, factors$residual_sd, factors$loading, sei, context,
        execution_plan, block_index, return_normalizer = TRUE
      )
    } else {
      key <- paste(K, points, sep = "/")
      if (K > 1L && !exists(key, designs, inherits = FALSE)) {
        design <- if (points == execution_plan$points_per_scramble) {
          execution_plan$designs[[as.character(K)]]
        } else {
          BayesTools::selection_qmc_design(
            dimensions = 2L * K, points = points,
            scrambles = control$scrambles, seed = control$seed
          )
        }
        assign(key, design, designs)
      }
      current_plan <- execution_plan
      current_plan$points_per_scramble <- points
      current_plan$designs[[as.character(K)]] <- if (K > 1L) {
        get(key, designs)
      } else numeric()
      .selection_joint_dense_loglik_block(
        yi, mean[active, , drop = FALSE],
        covariance_lower[active, , drop = FALSE], sei, context,
        current_plan, K, return_normalizer = TRUE
      )
    }
    error <- current$relative_mcse
    if (method %in% c("rank_one", "factor")) {
      error <- pmax(error, current$relative_change)
    }
    if (!is.null(previous)) {
      error <- pmax(error, abs(expm1(previous - current$log_normalizer)))
    }
    result$log_density[active]  <- -current$log_normalizer
    result$relative_mcse[active] <- error
    failed <- which(!is.finite(error) | error > control$relative_tolerance |
                      !is.finite(current$log_normalizer))
    if (!length(failed)) return(result)
    if (method %in% c("rank_one", "factor") ||
        points >= control$max_points_per_scramble) {
      stop(
        "Zplot marginal integration was rejected by diagnostics: relative integration error was ",
        format(error[failed[1L]], digits = 4),
        ". Increase 'max_points_per_scramble' in 'integration_control = set_selection_likelihood_control()'.",
        call. = FALSE
      )
    }
    active   <- active[failed]
    previous <- current$log_normalizer[failed]
    points   <- min(2L * points, control$max_points_per_scramble)
  }
}


# Depth of the forest the native nested rule integrates, or NA when the
# supports are not a forest: every pair must be nested or disjoint. A chain is
# the special case of one child per level. The rule expands one axis per level,
# so the depth, not the rank, decides what one rule costs.
.selection_factor_support_forest_depth <- function(support) {

  rank <- ncol(support)
  if (rank <= 1L) {
    return(0L)
  }
  ordered <- support[, order(colSums(support), decreasing = TRUE), drop = FALSE]
  depth <- integer(rank)
  for (outer in seq_len(rank)) {
    parent <- 0L
    for (inner in seq_len(outer - 1L)) {
      if (!any(ordered[, inner] & ordered[, outer])) {
        next
      }
      if (!all(ordered[, outer] <= ordered[, inner])) {
        return(NA_integer_)
      }
      parent <- inner
    }
    depth[[outer]] <- if (parent == 0L) 0L else depth[[parent]] + 1L
  }

  max(depth)
}


.zplot_joint_block <- function(
    z, mean, covariance_lower, sei, selection, probability, control,
    designs = new.env(parent = emptyenv()), factors = NULL) {

  S      <- nrow(mean)
  K      <- ncol(mean)
  points <- control$points_per_scramble
  active <- seq_len(S)
  result <- list(log_density = numeric(S), relative_mcse = numeric(S),
                 density = matrix(0, S, length(z)))
  previous <- NULL
  rank <- if (is.null(factors)) 0L else ncol(factors$loading) %/% K
  support <- factors$loading_support
  nested  <- FALSE
  if (rank > 0L && !is.null(support)) {
    if (!is.logical(support) || anyNA(support) ||
        !identical(dim(support), c(K, rank))) {
      stop("Selection factor loading supports are invalid.", call. = FALSE)
    }
    if (rank <= K &&
        identical(dim(factors$residual_sd), c(S, K)) &&
        identical(dim(factors$loading), c(S, K * rank)) &&
        is.numeric(factors$residual_sd) && is.numeric(factors$loading) &&
        all(is.finite(factors$residual_sd) & factors$residual_sd > 0) &&
        all(is.finite(factors$loading))) {
      loading_active <- matrix(vapply(seq_len(rank), function(column) {

        rowSums(factors$loading[, (column - 1L) * K + seq_len(K), drop = FALSE] != 0) > 0L
      }, logical(S)), S, rank)
      active_rank <- rowSums(loading_active)
      if (any(active_rank > 0L & active_rank < rank)) {
        groups <- split(seq_len(S), apply(loading_active, 1L, paste0, collapse = ""))
        for (rows in groups) {
          keep <- which(loading_active[rows[1L], ])
          # Entirely zero states already use the native diagonal calculation.
          if (!length(keep)) keep <- seq_len(rank)
          columns <- unlist(lapply(keep, function(column) {

            (column - 1L) * K + seq_len(K)
          }), use.names = FALSE)
          group_factors <- list(
            residual_sd = factors$residual_sd[rows, , drop = FALSE],
            loading = factors$loading[rows, columns, drop = FALSE],
            loading_support = support[, keep, drop = FALSE]
          )
          current <- .zplot_joint_block(
            z, mean[rows, , drop = FALSE], covariance_lower[rows, , drop = FALSE],
            sei, BayesTools::selection_context_subset_rows(selection, rows),
            probability, control, designs, group_factors
          )
          result$log_density[rows] <- current$log_density
          result$relative_mcse[rows] <- current$relative_mcse
          result$density[rows, ] <- current$density
        }
        return(result)
      }
    }
  }
  forest_depth <- if (rank >= 3L && !is.null(support)) {
    .selection_factor_support_forest_depth(support)
  } else NA_integer_
  nested <- rank >= 3L && !is.na(forest_depth)
  quadrature <- rank %in% c(1L, 2L) || nested
  orders <- if (rank == 1L) SELNORM_CLUSTER_QUADRATURE_ORDERS else
    SELNORM_FACTOR_QUADRATURE_ORDERS
  # Same node budget the kernels apply: a forest support pays
  # `order^(depth + 1)`, a plain tensor `order^rank`.
  cost_exponent <- if (nested) forest_depth + 1L else rank
  if (quadrature) {
    orders <- orders[.selnorm_factor_rule_affordable(orders, cost_exponent)]
    quadrature <- length(orders) > 0L
  }
  rule_index <- 1L
  repeat {
    rule <- if (quadrature) .gauss_hermite_nodes(orders[rule_index]) else NULL
    if (quadrature) points <- length(rule$nodes)
    key <- paste(K, points, sep = "/")
    if (!quadrature && !exists(key, designs, inherits = FALSE)) {
      assign(key, BayesTools::selection_qmc_design(
        dimensions = 2L * K, points = points,
        scrambles = control$scrambles, seed = control$seed
      ), designs)
    }
    context <- BayesTools::selection_context_subset_rows(selection, active)
    static  <- BayesTools::selection_native_static_args(context)
    current <- .Call(
      "RoBMA_selnorm_mnorm_zplot_batch",
      as.numeric(mean[1L, ]), .native_numeric_matrix(mean[active, , drop = FALSE]),
      .native_numeric_matrix(covariance_lower[active, , drop = FALSE]),
      as.numeric(sei), .native_numeric_matrix(context$omega),
      static$z_lower, static$z_upper, rep.int(1L, K), static$sign,
      static$telescope_probabilities, .native_integer_vector(context$kernel_mode),
      if (quadrature) numeric() else as.numeric(get(key, designs)),
      as.integer(points), as.integer(if (quadrature) 2L else control$scrambles),
      as.numeric(control$relative_tolerance), as.numeric(z), as.logical(probability),
      if (is.null(factors)) NULL else c(list(
        .native_numeric_matrix(factors$residual_sd[active, , drop = FALSE]),
        .native_numeric_matrix(factors$loading[active, , drop = FALSE])
      ), if (quadrature) list(rule$nodes, rule$log_weights, support)),
      .native_integer_vector(context[["vector_rule"]]),
      PACKAGE = "RoBMA"
    )
    error <- current$relative_mcse
    error[rowSums(!is.finite(current$density)) > 0L] <- Inf
    if (!is.null(previous)) {
      peak <- apply(current$density, 1L, max)
      change <- apply(abs(current$density - previous$density), 1L, max)
      positive <- is.finite(peak) & peak > 0
      change[positive] <- change[positive] / peak[positive]
      normalizer_change <- abs(expm1(current$log_density - previous$log_density))
      change[!is.finite(change)] <- Inf
      normalizer_change[!is.finite(normalizer_change)] <- Inf
      error <- pmax(error, change, normalizer_change)
    }
    accepted_error <- if (quadrature) {
      if (is.null(previous)) rep(Inf, length(active)) else pmax(error, previous$error)
    } else {
      error
    }
    result$log_density[active]  <- current$log_density
    result$relative_mcse[active] <- error
    result$density[active, ]    <- current$density
    failed <- which(!is.finite(accepted_error) | accepted_error > control$relative_tolerance |
                      !is.finite(current$log_density))
    if (!length(failed)) return(result)
    if (!quadrature && points >= control$max_points_per_scramble) {
      stop(
        "Zplot marginal integration was rejected by diagnostics: relative integration error was ",
        format(error[failed[1L]], digits = 4),
        ". Increase 'max_points_per_scramble' in 'integration_control = set_selection_likelihood_control()'.",
        call. = FALSE
      )
    }
    active <- active[failed]
    previous <- list(log_density = current$log_density[failed],
                     density = current$density[failed, , drop = FALSE],
                     error = if (is.null(previous)) rep(Inf, length(failed)) else error[failed])
    if (quadrature) {
      rule_index <- rule_index + 1L
      if (rule_index > length(orders)) {
        quadrature <- FALSE
        previous   <- NULL
        points     <- control$points_per_scramble
      }
    } else {
      points <- min(2L * points, control$max_points_per_scramble)
    }
  }
}


# Conditional selection normalizes before integrating the Gaussian latent
# effect. Only its one-dimensional marginal is needed for this scalar display.
.zplot_latent_mixture <- function(
    z, mean, sd, latent_sd, sei, selection, probability, control,
    mean_extrapolated = mean, fitted_only = FALSE) {

  # Resolve the default before refinement subsets its source argument.
  force(mean_extrapolated)
  S <- nrow(mean)
  normal_rows <- which(selection$use_normal)
  if (length(normal_rows)) {
    total_sd <- .root_sum_squares(sd, latent_sd)
    normal_result <- function(location) {

      if (!probability) return(.zplot_normal_density_matrix(z, location, total_sd, sei))
      q <- matrix(z * sei, S, length(sei), byrow = TRUE)
      matrix(rowMeans(stats::pnorm(q, location, total_sd, lower.tail = FALSE) +
                        stats::pnorm(-q, location, total_sd)), S, 1L)
    }
    fitted       <- normal_result(mean)
    extrapolated <- normal_result(mean_extrapolated)
    weights      <- rep(1, S)
    active       <- which(!selection$use_normal)
    if (length(active)) {
      selected <- .zplot_latent_mixture(z, mean[active, , drop = FALSE],
        sd[active, , drop = FALSE], latent_sd[active, , drop = FALSE], sei,
        BayesTools::selection_context_subset_rows(selection, active),
        probability, control, mean_extrapolated[active, , drop = FALSE], fitted_only)
      fitted[active, ]       <- selected$fitted
      extrapolated[active, ] <- selected$extrapolated
      weights[active] <- selected$weights
    }
    return(list(fitted = fitted, extrapolated = extrapolated, weights = weights,
                EDR = if (probability && !fitted_only) extrapolated[, 1L] / weights else NULL))
  }
  previous <- previous_rule <- NULL
  active   <- seq_len(S)
  result <- list(fitted = matrix(0, S, length(z)),
                 extrapolated = matrix(0, S, length(z)), weights = numeric(S))
  for (order in SELNORM_CLUSTER_QUADRATURE_ORDERS) {
    rule <- .gauss_hermite_nodes(order)
    fitted       <- matrix(0, S, length(z))
    extrapolated <- matrix(0, S, length(z))
    weights      <- numeric(S)
    if (!probability && !is.null(previous_rule)) {
      attr(rule, "density_control") <- list(as.double(control$relative_tolerance), previous_rule)
      # Condition the retained Gaussian source on the evaluated outcome before
      # quadrature. Integrating its original population density can miss the
      # arbitrarily narrow sampling kernel when the integrated SD is small.
      if (fitted_only) {
        fitted <- .zplot_selnorm_density_matrix(z, mean, sd, sei, selection,
          FALSE, latent_sd = latent_sd, quadrature = rule)
      } else {
        pair <- .zplot_selnorm_density_pair(z, mean, sd, sei, selection,
          latent_sd = latent_sd, quadrature = rule)
        fitted       <- pair$fitted
        extrapolated <- pair$extrapolated
        if (!identical(mean, mean_extrapolated)) {
          extrapolated <- .zplot_selnorm_density_matrix(z, mean_extrapolated,
            sd, sei, selection, TRUE, latent_sd = latent_sd, quadrature = rule)
        }
      }
    }
    if (probability || !fitted_only) {
      for (j in seq_along(rule$nodes)) {
        conditional_mean <- mean + latent_sd * rule$nodes[j]
        conditional_extrapolated <- mean_extrapolated + latent_sd * rule$nodes[j]
        if (probability && fitted_only) {
          fitted[, 1L] <- fitted[, 1L] + rule$weights[j] *
            .zplot_selnorm_threshold_summary(z, conditional_mean, sd, sei,
              selection, FALSE)$EDR
        } else if (probability) {
          fit <- .zplot_selnorm_threshold_summary(z, conditional_mean, sd, sei, selection, FALSE)
          ext <- .zplot_selnorm_threshold_summary(z, conditional_extrapolated, sd, sei, selection, TRUE)
          fitted[, 1L] <- fitted[, 1L] + rule$weights[j] * fit$EDR
          extrapolated[, 1L] <- extrapolated[, 1L] + rule$weights[j] * ext$EDR * ext$weights
          weights <- weights + rule$weights[j] * ext$weights
        } else {
          weights <- weights + rule$weights[j] * rowMeans(
            .zplot_inverse_selection_weights(conditional_mean, sd, sei, selection)
          )
        }
      }
    }
    current <- list(fitted = fitted, extrapolated = extrapolated,
                    weights = matrix(weights, S, 1L))
    if (!is.null(previous)) {
      error <- Reduce(pmax, lapply(names(current), function(component) {
        value <- current[[component]]
        if (!probability && component != "weights" && !(fitted_only && component == "extrapolated")) {
          error <- attr(value, "relative_integration_error")
          if (!is.numeric(error) || length(error) != S || anyNA(error) || any(error < 0)) {
            stop("Zplot density integration diagnostics are unavailable.", call. = FALSE)
          }
          return(error)
        }
        peak  <- apply(abs(value), 1L, max)
        change <- apply(abs(value - previous[[component]]), 1L, max)
        positive <- is.finite(peak) & peak > 0
        change[positive] <- change[positive] / peak[positive]
        change[!is.finite(change)] <- Inf
        change
      }))
      accepted <- which(is.finite(error) & error <= control$relative_tolerance)
      result$fitted[active[accepted], ]       <- fitted[accepted, , drop = FALSE]
      result$extrapolated[active[accepted], ] <- extrapolated[accepted, , drop = FALSE]
      result$weights[active[accepted]]       <- weights[accepted]
      failed <- setdiff(seq_len(S), accepted)
      if (!length(failed)) {
        return(c(result, list(
          EDR = if (probability && !fitted_only) result$extrapolated[, 1L] / result$weights else NULL
        )))
      }
      # Each draw has its own convergence criterion; retain completed rows.
      active            <- active[failed]
      S                 <- length(active)
      mean              <- mean[failed, , drop = FALSE]
      mean_extrapolated  <- mean_extrapolated[failed, , drop = FALSE]
      sd                <- sd[failed, , drop = FALSE]
      latent_sd         <- latent_sd[failed, , drop = FALSE]
      selection <- BayesTools::selection_context_subset_rows(selection, failed)
      current   <- lapply(current, function(value) value[failed, , drop = FALSE])
    }
    previous <- current
    attr(rule, "density_control") <- NULL
    previous_rule <- rule
  }
  stop(
    "Zplot marginal integration was rejected by diagnostics: relative integration error was ",
    format(max(error), digits = 4), ". Inspect the fitted selection and heterogeneity parameters.",
    call. = FALSE
  )
}


.zplot_stored_conditioning_depth <- function(object) {

  conditioning_depth <- object[["zplot"]][["data"]][["conditioning_depth"]]
  if (is.null(conditioning_depth)) {
    conditioning_depth <- if (.is_multilevel(object)) "cluster" else "marginal"
  }

  return(.normalize_conditioning_depth(conditioning_depth))
}


# The pre-selection Gaussian reference marginalizes every random source.
# Its scalar row projections need only the covariance diagonal; source
# retention affects the selected law, not this reference covariance.
.zplot_gaussian_marginal_reference <- function(object, posterior_samples,
                                               predictive) {

  S <- nrow(posterior_samples)
  K <- length(predictive[["sei"]])
  mean <- predictive[["mu_extrapolated"]]
  heterogeneity <- predictive[["tau_within"]]
  if (is.null(heterogeneity)) {
    heterogeneity <- .zplot_predictive_heterogeneity(
      object, posterior_samples, conditioning_depth = "marginal"
    )
  }
  heterogeneity <- .expand_brma_mv_heterogeneity_samples(heterogeneity, S, K)
  data <- object[["data"]]
  sampling_variance <- if (.is_data_known_v(data)) {
    .known_v_diagonal(.data_known_v_data(data))
  } else {
    predictive[["sei"]]^2
  }
  variance <- sweep(heterogeneity^2, 2L, sampling_variance, "+")
  if (!is.matrix(mean) || !identical(dim(mean), c(S, K)) ||
      length(sampling_variance) != K || any(!is.finite(mean)) ||
      any(!is.finite(variance)) || any(variance < 0)) {
    stop("Gaussian zplot marginal reference parameters are invalid.", call. = FALSE)
  }
  list(mu = mean, variance = variance, sei = predictive[["sei"]])
}


.zplot_predictive_heterogeneity <- function(object, posterior_samples,
                                            conditioning_depth) {

  if (conditioning_depth == "estimate") {
    conditional_variance <- .zplot_estimate_conditional_variance(
      object            = object,
      posterior_samples = posterior_samples
    )
    return(sqrt(conditional_variance))
  }

  if (inherits(object, "brma.mv") && conditioning_depth == "marginal") {
    components <- .brma_mv_heterogeneity_components(
      object                         = object,
      posterior_samples              = posterior_samples,
      include_known_group_covariance = TRUE
    )
    return(.total_brma_mv_heterogeneity_samples(components))
  }

  tau_result <- .zplot_tau_samples(
    object            = object,
    posterior_samples = posterior_samples
  )

  if (conditioning_depth == "cluster") {
    return(tau_result[["tau_within"]])
  }

  return(tau_result[["tau_total"]])
}


.zplot_tau_samples <- function(object, posterior_samples) {

  K <- length(.outcome_data_sei(object))

  return(.evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) {
      .create_fit_formula_list(data = object[["data"]], "scale")
    } else {
      NULL
    },
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = .is_multilevel(object),
    K                 = K,
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  ))
}


# Posterior conditional variance of the fitted latent true effects. The
# corresponding conditional means are returned by predict(type = "blup").
.zplot_estimate_conditional_variance <- function(object, posterior_samples) {

  if (inherits(object, "brma.mv") && .is_random(object)) {
    return(.zplot_mv_random_conditional_variance(
      object            = object,
      posterior_samples = posterior_samples
    ))
  }

  tau_result <- .zplot_tau_samples(
    object            = object,
    posterior_samples = posterior_samples
  )
  tau_within <- tau_result[["tau_within"]]
  fit_vi     <- .outcome_data_sei(object)^2 / .outcome_data_weights(object)

  if (.is_multilevel(object)) {
    if (.is_weightfunction(object)) {
      return(.zplot_independent_conditional_variance(tau_within, fit_vi))
    }

    return(.zplot_multilevel_conditional_variance(
      tau_within  = tau_within,
      tau_between = tau_result[["tau_between"]],
      vi          = fit_vi,
      cluster     = object[["data"]][["outcome"]][["cluster"]]
    ))
  }

  if (.is_data_known_v(object[["data"]])) {
    return(.zplot_known_v_conditional_variance(
      tau_within = tau_within,
      known_V    = .data_known_v_data(object[["data"]])
    ))
  }

  return(.zplot_independent_conditional_variance(tau_within, fit_vi))
}


.zplot_independent_conditional_variance <- function(tau_within, vi) {

  tau2       <- tau_within^2
  vi_samples <- matrix(vi, nrow = nrow(tau_within), ncol = ncol(tau_within),
                       byrow = TRUE)
  denominator <- tau2 + vi_samples

  if (any(!is.finite(denominator)) || any(denominator <= 0)) {
    stop(
      "Cannot evaluate the estimate-depth conditional variance from the fitted ",
      "heterogeneity and sampling variances.",
      call. = FALSE
    )
  }

  return(tau2 * vi_samples / denominator)
}


.zplot_multilevel_conditional_variance <- function(tau_within, tau_between,
                                                   vi, cluster) {

  S             <- nrow(tau_within)
  K             <- ncol(tau_within)
  block_indices <- .get_multilevel_block_indices(cluster)
  out           <- matrix(0, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    for (idx in block_indices) {
      within2     <- tau_within[s, idx]^2
      denominator <- within2 + vi[idx]
      if (any(!is.finite(denominator)) || any(denominator <= 0)) {
        stop(
          "Cannot evaluate the multilevel estimate-depth conditional variance.",
          call. = FALSE
        )
      }

      within_variance <- within2 * vi[idx] / denominator
      gamma_variance  <- 1 / (
        1 + sum(tau_between[s, idx]^2 / denominator)
      )
      gamma_loading  <- tau_between[s, idx] * vi[idx] / denominator

      out[s, idx] <- within_variance + gamma_loading^2 * gamma_variance
    }
  }

  return(out)
}


.zplot_known_v_conditional_variance <- function(tau_within, known_V) {

  S <- nrow(tau_within)
  K <- ncol(tau_within)
  if (.known_v_nrow(known_V) != K) {
    stop(
      "Known-V covariance dimensions do not match the estimate-depth target.",
      call. = FALSE
    )
  }

  block_data <- .known_v_blocks(known_V)
  .known_v_validate_dependency_blocks(
    lapply(block_data, `[[`, "index"),
    K
  )

  out <- matrix(0, nrow = S, ncol = K)
  for (block in block_data) {
    idx     <- block[["index"]]
    V_block <- block[["covariance"]]

    if (length(idx) == 1L) {
      tau2       <- tau_within[, idx]^2
      denominator <- tau2 + V_block[1L, 1L]
      if (any(!is.finite(denominator)) || any(denominator <= 0)) {
        stop(
          "Cannot evaluate a known-V estimate-depth conditional variance block.",
          call. = FALSE
        )
      }
      out[, idx] <- tau2 * V_block[1L, 1L] / denominator
      next
    }

    for (s in seq_len(S)) {
      latent_covariance <- diag(
        tau_within[s, idx]^2,
        nrow = length(idx),
        ncol = length(idx)
      )
      out[s, idx] <- .zplot_gaussian_conditional_variance(
        latent_covariance   = latent_covariance,
        sampling_covariance = V_block
      )
    }
  }

  return(out)
}


.zplot_mv_random_conditional_variance <- function(object, posterior_samples,
                                                  max_bytes = NULL) {

  known_V <- .data_known_v_data(object[["data"]])
  K       <- nrow(object[["data"]][["outcome"]])
  S       <- nrow(posterior_samples)
  if (is.null(known_V) || .known_v_nrow(known_V) != K) {
    stop(
      "Random-formula brma.mv estimate-depth variance requires matching known-V metadata.",
      call. = FALSE
    )
  }

  sampling_covariance <- .known_v_materialize(known_V)
  out                 <- matrix(0, nrow = S, ncol = K)
  chunks              <- .known_v_covariance_chunk_indices(
    S         = S,
    K         = K,
    max_bytes = max_bytes
  )

  for (rows in chunks) {
    random_vcov <- .brma_mv_random_effects_marginal_vcov(
      object            = object,
      posterior_samples = posterior_samples[rows, , drop = FALSE],
      diagonal_only     = FALSE,
      data              = object[["data"]],
      new_levels        = "error"
    )
    covariance_samples <- random_vcov[["samples"]]
    expected_dim       <- c(length(rows), K, K)
    if (!is.numeric(covariance_samples) ||
        !identical(dim(covariance_samples), expected_dim) ||
        any(!is.finite(covariance_samples))) {
      stop(
        "Random-effect covariance samples have inconsistent dimensions.",
        call. = FALSE
      )
    }

    for (draw_i in seq_along(rows)) {
      out[rows[draw_i], ] <- .zplot_gaussian_conditional_variance(
        latent_covariance   = matrix(
          covariance_samples[draw_i, , ],
          nrow = K,
          ncol = K
        ),
        sampling_covariance = sampling_covariance
      )
    }
  }

  return(out)
}


# Matheron's conditional Gaussian covariance, evaluated as a sum of squares.
.zplot_gaussian_conditional_variance <- function(latent_covariance,
                                                 sampling_covariance) {

  if (!is.matrix(latent_covariance) ||
      !is.matrix(sampling_covariance) ||
      !identical(dim(latent_covariance), dim(sampling_covariance))) {
    stop(
      "Latent and sampling covariance matrices must have matching dimensions.",
      call. = FALSE
    )
  }

  marginal_covariance <- latent_covariance + sampling_covariance
  chol_marginal <- .covariance_cholesky(
    .covariance_factorization(marginal_covariance)
  )
  if (is.null(chol_marginal)) {
    stop(
      "The estimate-depth Gaussian conditional variance is unavailable because positive definiteness of the marginal covariance cannot be resolved at working precision.",
      call. = FALSE
    )
  }

  latent_factor <- .covariance_sampling_factor(
    .covariance_factorization(latent_covariance)
  )
  sampling_factor <- .covariance_sampling_factor(
    .covariance_factorization(sampling_covariance)
  )
  if (is.null(latent_factor) || is.null(sampling_factor)) {
    stop("Latent and sampling covariance matrices must be positive semidefinite.",
         call. = FALSE)
  }

  # For T = Q + V, the conditional covariance is
  # (T^-1 V)' Q (T^-1 V) + (T^-1 Q)' V (T^-1 Q).
  # Solve both terms directly: subtracting Q T^-1 Q from Q can erase an
  # entire small conditional variance even when the result is nonnegative.
  solved_sampling <- backsolve(
    chol_marginal,
    forwardsolve(t(chol_marginal), sampling_covariance)
  )
  solved_latent <- backsolve(
    chol_marginal,
    forwardsolve(t(chol_marginal), latent_covariance)
  )
  conditional_variance <- colSums((latent_factor %*% solved_sampling)^2) +
    colSums((sampling_factor %*% solved_latent)^2)

  # A PSD source with a zero diagonal has an exactly zero row and column.
  # Either a deterministic latent coordinate or a noiseless observation
  # therefore has exactly zero conditional variance in that coordinate.
  deterministic <- diag(latent_covariance) == 0 | diag(sampling_covariance) == 0
  conditional_variance[deterministic] <- 0
  if (any(!is.finite(conditional_variance))) {
    stop("Cannot compute the estimate-depth conditional latent variance.",
         call. = FALSE)
  }
  unname(conditional_variance)
}


# ---------------------------------------------------------------------------- #
# .zplot_threshold_vectorized
# ---------------------------------------------------------------------------- #
#
# Vectorized EDR computation for normal and selected-normal rows.
#
# ---------------------------------------------------------------------------- #
.zplot_total_sd <- function(tau_within, sei) {

  total_sd <- vapply(
    seq_len(ncol(tau_within)),
    function(i) .root_sum_squares(tau_within[, i], sei[[i]]),
    numeric(nrow(tau_within))
  )
  total_sd <- matrix(
    total_sd,
    nrow     = nrow(tau_within),
    ncol     = ncol(tau_within),
    dimnames = dimnames(tau_within)
  )

  return(total_sd)
}

.zplot_threshold_vectorized <- function(z_threshold, mu_samples, tau_within,
                                         sei, selection, extrapolate,
                                         effect_direction) {

  S        <- nrow(mu_samples)
  K        <- ncol(mu_samples)
  total_sd <- .zplot_total_sd(tau_within, sei)

  if (!is.null(selection) && .has_native_zplot_threshold()) {
    return(.zplot_selnorm_threshold_summary(
      z_threshold       = z_threshold,
      mean              = mu_samples,
      sd                = total_sd,
      sei               = sei,
      selection_context = selection,
      extrapolate       = extrapolate
    ))
  }

  q_upper   <- z_threshold * sei
  q_lower   <- -z_threshold * sei

  thresholds <- stats::pnorm(
    matrix(q_upper, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = FALSE
  ) + stats::pnorm(
    matrix(q_lower, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = TRUE
  )
  weights       <- .zplot_inverse_selection_weights(
    mean      = mu_samples,
    sd        = total_sd,
    sei       = sei,
    selection = selection
  )
  weighted_rows <- if (is.null(selection)) integer(0) else which(!selection[["use_normal"]])
  if (length(weighted_rows) > 0) {
    selection_weight <- BayesTools::selection_context_subset_rows(
      context = selection,
      rows    = weighted_rows
    )
    mean_weight      <- mu_samples[weighted_rows, , drop = FALSE]
    sd_weight        <- total_sd[weighted_rows, , drop = FALSE]

    if (!extrapolate) {
      prob_upper <- .selection_step_cdf_matrix(
        q                 = q_upper,
        mean              = mean_weight,
        sd                = sd_weight,
        sei               = sei,
        selection_context = selection_weight,
        lower.tail        = FALSE
      )
      prob_lower <- .selection_step_cdf_matrix(
        q                 = q_lower,
        mean              = mean_weight,
        sd                = sd_weight,
        sei               = sei,
        selection_context = selection_weight,
        lower.tail        = TRUE
      )

      thresholds[weighted_rows, ] <- prob_upper + prob_lower
    }
  }

  if (!is.null(selection) && extrapolate) {
    EDR <- rowSums(thresholds * weights) / rowSums(weights)
  } else {
    EDR <- rowMeans(thresholds)
  }

  return(list(
    EDR     = EDR,
    weights = rowMeans(weights)
  ))
}


# ---------------------------------------------------------------------------- #
# .zplot_selnorm_threshold_summary
# ---------------------------------------------------------------------------- #
#
# Native EDR and inverse-selection-weight reductions for zplot thresholds.
#
# ---------------------------------------------------------------------------- #
.has_native_zplot_threshold <- function() {

  return(is.loaded(
    "RoBMA_selnorm_zcurve_threshold_summary",
    PACKAGE = "RoBMA"
  ))
}

.zplot_selnorm_threshold_summary <- function(z_threshold, mean, sd, sei,
                                             selection_context, extrapolate) {

  .selection_require_step_evaluable(selection_context, ".zplot_threshold_vectorized()")
  native_static <- BayesTools::selection_native_static_args(selection_context)

  return(.Call(
    "RoBMA_selnorm_zcurve_threshold_summary",
    .native_numeric_vector(z_threshold),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    as.logical(extrapolate),
    native_static[["telescope_probabilities"]],
    PACKAGE = "RoBMA"
  ))
}


# ---------------------------------------------------------------------------- #
# .zplot_density_vectorized
# ---------------------------------------------------------------------------- #
#
# Vectorized z-density computation for normal and selected-normal rows.
#
# ---------------------------------------------------------------------------- #
.has_native_zplot_density <- function(selection = FALSE) {

  symbols <- "RoBMA_zcurve_normal_density_matrix"
  if (selection) {
    symbols <- c(symbols, "RoBMA_selnorm_zcurve_density_matrix")
  }

  return(all(vapply(symbols, is.loaded, logical(1), PACKAGE = "RoBMA")))
}

.zplot_normal_density_matrix <- function(z_sequence, mean, sd, sei) {

  if (!.has_native_zplot_density(selection = FALSE)) {
    stop("The native zplot density kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_zcurve_normal_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    PACKAGE = "RoBMA"
  ))
}

.zplot_selnorm_density_matrix <- function(z_sequence, mean, sd, sei,
                                           selection_context, extrapolate,
                                           latent_sd = NULL, quadrature = NULL) {

  .selection_require_step_evaluable(selection_context, ".zplot_density_vectorized()")

  if (!.has_native_zplot_density(selection = TRUE)) {
    stop("The native selected-normal zplot density kernel is not loaded.", call. = FALSE)
  }
  native_static <- BayesTools::selection_native_static_args(selection_context)

  return(.Call(
    "RoBMA_selnorm_zcurve_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    as.logical(extrapolate),
    native_static[["telescope_probabilities"]],
    if (is.null(latent_sd)) NULL else .native_numeric_matrix(latent_sd),
    quadrature,
    PACKAGE = "RoBMA"
  ))
}

.zplot_selnorm_density_pair <- function(z_sequence, mean, sd, sei,
                                         selection_context,
                                         latent_sd = NULL, quadrature = NULL) {

  .selection_require_step_evaluable(selection_context, ".zplot_density_pair()")

  if (!.has_native_zplot_density(selection = TRUE)) {
    stop("The native selected-normal zplot density kernel is not loaded.", call. = FALSE)
  }
  native_static <- BayesTools::selection_native_static_args(selection_context)

  return(.Call(
    "RoBMA_selnorm_zcurve_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    native_static[["sign"]],
    native_static[["phack_q"]],
    native_static[["phack_z_source"]],
    native_static[["phack_z_dest"]],
    native_static[["segment_bounds"]],
    native_static[["segment_step_bin"]],
    native_static[["segment_phack_region"]],
    c(FALSE, TRUE),
    native_static[["telescope_probabilities"]],
    if (is.null(latent_sd)) NULL else .native_numeric_matrix(latent_sd),
    quadrature,
    PACKAGE = "RoBMA"
  ))
}

.zplot_density_pair <- function(object, z_sequence, max_samples,
                                conditioning_depth = "marginal",
                                integration_control = set_selection_likelihood_control(),
                                parallel = FALSE, cores = min(4, RoBMA.get_option("max_cores"))) {

  posterior_samples <- .get_posterior_samples(object[["fit"]])
  selected_ind      <- .thin_sample_rows(nrow(posterior_samples), max_samples)
  if (!is.null(selected_ind)) {
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }

  if (.zplot_requires_selection_marginal(object, conditioning_depth)) {
    if (parallel && inherits(object, "brma.mv") &&
        identical(conditioning_depth, "marginal") && cores > 1L && nrow(posterior_samples) >= cores) {
      return(.zplot_selection_marginal_parallel(object, posterior_samples, z_sequence,
        conditioning_depth, integration_control, cores))
    }
    return(.zplot_selection_marginal(object, posterior_samples, z_sequence, NULL,
      conditioning_depth, integration_control))
  }

  predictive_fit <- .zplot_predictive_components(
    object             = object,
    posterior_samples  = posterior_samples,
    extrapolate        = FALSE,
    conditioning_depth = conditioning_depth
  )
  has_location_bias <- .is_PET(object) || .is_PEESE(object)
  if (has_location_bias) {
    predictive_extrapolated <- .zplot_predictive_components(
      object                      = object,
      posterior_samples           = posterior_samples,
      extrapolate                 = TRUE,
      conditioning_depth          = conditioning_depth,
      predictive_heterogeneity    = predictive_fit[["tau_within"]]
    )
  } else {
    predictive_extrapolated <- predictive_fit
  }
  selection <- .zplot_selection_context(
    object            = object,
    posterior_samples = posterior_samples,
    is_weightfunction = .is_weightfunction(object)
  )

  same_predictive <- identical(predictive_fit, predictive_extrapolated)
  if (same_predictive && !is.null(selection)) {
    total_sd <- .zplot_total_sd(
      predictive_fit[["tau_within"]],
      predictive_fit[["sei"]]
    )

    return(.zplot_selnorm_density_pair(
      z_sequence        = z_sequence,
      mean              = predictive_fit[["mu"]],
      sd                = total_sd,
      sei               = predictive_fit[["sei"]],
      selection_context = selection
    ))
  }

  if (same_predictive && is.null(selection)) {
    density <- .zplot_density_vectorized(
      z_sequence       = z_sequence,
      mu_samples       = predictive_fit[["mu"]],
      tau_within       = predictive_fit[["tau_within"]],
      sei              = predictive_fit[["sei"]],
      selection        = NULL,
      extrapolate      = FALSE,
      effect_direction = .effect_direction(object)
    )
    return(list(fitted = density, extrapolated = density))
  }

  return(list(
    fitted = .zplot_density_vectorized(
      z_sequence       = z_sequence,
      mu_samples       = predictive_fit[["mu"]],
      tau_within       = predictive_fit[["tau_within"]],
      sei              = predictive_fit[["sei"]],
      selection        = selection,
      extrapolate      = FALSE,
      effect_direction = .effect_direction(object)
    ),
    extrapolated = .zplot_density_vectorized(
      z_sequence       = z_sequence,
      mu_samples       = predictive_extrapolated[["mu"]],
      tau_within       = predictive_extrapolated[["tau_within"]],
      sei              = predictive_extrapolated[["sei"]],
      selection        = selection,
      extrapolate      = TRUE,
      effect_direction = .effect_direction(object)
    )
  ))
}

.zplot_density_vectorized <- function(z_sequence, mu_samples, tau_within,
                                       sei, selection, extrapolate,
                                       effect_direction) {

  total_sd <- .zplot_total_sd(tau_within, sei)

  if (is.null(selection)) {
    return(.zplot_normal_density_matrix(
      z_sequence = z_sequence,
      mean       = mu_samples,
      sd         = total_sd,
      sei        = sei
    ))
  }

  density <- .zplot_selnorm_density_matrix(
    z_sequence        = z_sequence,
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    selection_context = selection,
    extrapolate       = extrapolate
  )

  return(density)
}


# ---------------------------------------------------------------------------- #
# .zplot_inverse_selection_weights
# ---------------------------------------------------------------------------- #
#
# Reference-bin-normalized inverse-weight extrapolation factors for the
# released univariate diagnostic, not identified attempted-study counts.
# A count interpretation requires additional reporting/stopping assumptions.
# Normal/no-bias branches have weight one.
#
# ---------------------------------------------------------------------------- #
.zplot_inverse_selection_weights <- function(mean, sd, sei, selection) {

  weights <- matrix(1, nrow = nrow(mean), ncol = ncol(mean))

  if (is.null(selection)) {
    return(weights)
  }

  weighted_rows <- which(!selection[["use_normal"]])
  if (length(weighted_rows) == 0L) {
    return(weights)
  }

  selection_weight <- BayesTools::selection_context_subset_rows(
    context = selection,
    rows    = weighted_rows
  )
  log_norm <- .selection_step_log_norm_matrix(
    mean              = mean[weighted_rows, , drop = FALSE],
    sd                = sd[weighted_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_weight
  )

  weights[weighted_rows, ] <- exp(-log_norm)
  return(weights)
}


# ---------------------------------------------------------------------------- #
# .zplot_selection_context
# ---------------------------------------------------------------------------- #
#
# Prepare posterior-row selection metadata for branch-aware zplot evaluation.
#
# ---------------------------------------------------------------------------- #
.zplot_selection_context <- function(object, posterior_samples,
                                     is_weightfunction) {

  if (!is_weightfunction) {
    return(NULL)
  }

  return(.selection_context(
    object            = object,
    posterior_samples = posterior_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .zplot_selection_args
# ---------------------------------------------------------------------------- #
#
# Extract active omega and cutpoints for one posterior row and estimate.
#
# ---------------------------------------------------------------------------- #
.zplot_selection_args <- function(selection, row, estimate, n = 1L) {

  .selection_require_step_evaluable(selection, ".zplot_selection_args()")

  omega <- selection[["omega"]][row, , drop = FALSE]

  if (n > 1L) {
    omega <- matrix(as.numeric(omega), nrow = n, ncol = ncol(omega), byrow = TRUE)
  }

  return(list(
    omega   = omega,
    crit_yi = stats::qnorm(rev(selection[["p_cuts"]])[-c(1L, length(selection[["p_cuts"]]))],
                           lower.tail = FALSE) *
      selection[["sei"]][estimate]
  ))
}


# ============================================================================ #
# Graphical Helper Functions
# ============================================================================ #
#
# These functions extract and set default graphical parameters for zplot
# plotting components. They handle the translation between base graphics
# and ggplot2 parameter naming conventions.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .get_dots_hist_zplot
# ---------------------------------------------------------------------------- #
#
# Extracts histogram graphical parameters with defaults.
#
# @param dots      list of user-supplied graphical parameters
# @param plot_type "base" or "ggplot"
# @param max_density maximum density value for setting ylim
#
# @return list of graphical parameters appropriate for plot_type
#
