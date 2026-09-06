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
                             integration_control = set_selection_likelihood_control()) {

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
    result <- .zplot_selection_marginal(
      object, posterior_samples, z_sequence, z_threshold,
      conditioning_depth, integration_control
    )
    if (!is.null(z_threshold)) {
      return(list(
        EDR     = if (extrapolate) result$EDR else result$fitted[, 1L],
        weights = result$weights
      ))
    }
    return(result[[if (extrapolate) "extrapolated" else "fitted"]])
  }

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
       .is_multilevel(object) || inherits(object, "brma.mv"))
}


.zplot_selection_marginal <- function(
    object, posterior_samples, z_sequence, z_threshold,
    conditioning_depth, integration_control) {

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

  if (.is_data_exact_selection(object[["data"]])) {
    return(.zplot_exact_marginal(
      object, posterior_samples, predictive, selection, z, probability, control
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
      stop("The approximate selection zplot target is unavailable without a conditional sampling factorization.",
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
    probability, control, predictive$mu_extrapolated
  )
}


# Both curves use the same Gaussian block and the same joint normalizer.
# Its inverse is the area of the extrapolated curve per observed estimate.
.zplot_exact_marginal <- function(
    object, posterior_samples, predictive, selection, z, probability, control) {

  S <- nrow(posterior_samples)
  K <- length(predictive$sei)
  scales <- if (.is_data_random(object[["data"]])) {
    list(tau_within = matrix(0, S, K), tau_between = matrix(0, S, K))
  } else {
    .zplot_tau_samples(object, posterior_samples)
  }
  fitted       <- matrix(0, S, length(z))
  extrapolated <- matrix(0, S, length(z))
  weights      <- numeric(S)
  designs      <- new.env(parent = emptyenv())
  exact_plan <- .data_exact_selection_setup(object[["data"]])
  factor_setup <- list(
    fit = object[["fit"]], data = object[["data"]], priors = object[["priors"]],
    posterior_samples = posterior_samples, S = S, K = K,
    tau_within = scales$tau_within, tau_between = scales$tau_between,
    is_multilevel = .is_multilevel(object)
  )
  random_factors <- .selection_exact_random_factor_samples(factor_setup)
  random_covariance <- if (is.null(random_factors)) {
    .selection_exact_random_covariance_samples(factor_setup)
  } else {
    NULL
  }
  for (block in seq_along(exact_plan$row_blocks)) {
    rows <- exact_plan$row_blocks[[block]]
    k <- length(rows)
    factors <- if (exact_plan$block_methods[block] %in% c("factor", "cluster")) {
      .selection_exact_factor_block_samples(factor_setup, block, random_factors)
    } else {
      NULL
    }
    if (!is.null(factors) &&
        (factors$rank == 0L || factors$rank > k || any(factors$residual_sd <= 0))) {
      factors <- NULL
    }
    pairs <- .selection_exact_lower_pairs(exact_plan, seq_len(k))
    covariance <- .selection_exact_covariance_lower(
      factor_setup, block, random_covariance, random_factors
    )
    result <- .zplot_joint_block(
      z, predictive$mu[, rows, drop = FALSE], covariance,
      predictive$sei[rows], selection, probability, control, designs, factors
    )
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
    fitted       <- fitted + result$density * (k / K)
    extrapolated <- extrapolated + normal * inverse * (k / K)
    weights      <- weights + inverse * (k / K)
  }
  list(fitted = fitted, extrapolated = extrapolated, weights = weights,
       EDR = if (probability) extrapolated[, 1L] / weights else NULL)
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
  quadrature <- rank %in% c(1L, 2L)
  orders <- if (rank == 1L) SELNORM_CLUSTER_QUADRATURE_ORDERS else
    SELNORM_FACTOR_QUADRATURE_ORDERS[["2"]]
  rule_index <- 1L
  repeat {
    rule <- if (quadrature) .gauss_hermite_nodes(orders[rule_index]) else NULL
    if (quadrature) points <- length(rule$nodes)^rank
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
      ), if (quadrature) list(rule$nodes, rule$log_weights)),
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


# Approximate selection normalizes before integrating the Gaussian latent
# effect. Only its one-dimensional marginal is needed for this scalar display.
.zplot_latent_mixture <- function(
    z, mean, sd, latent_sd, sei, selection, probability, control,
    mean_extrapolated = mean) {

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
        probability, control, mean_extrapolated[active, , drop = FALSE])
      fitted[active, ]       <- selected$fitted
      extrapolated[active, ] <- selected$extrapolated
      weights[active] <- selected$weights
    }
    return(list(fitted = fitted, extrapolated = extrapolated, weights = weights,
                EDR = if (probability) extrapolated[, 1L] / weights else NULL))
  }
  previous <- NULL
  for (order in SELNORM_CLUSTER_QUADRATURE_ORDERS) {
    rule <- .gauss_hermite_nodes(order)
    fitted       <- matrix(0, S, length(z))
    extrapolated <- matrix(0, S, length(z))
    weights      <- numeric(S)
    for (j in seq_along(rule$nodes)) {
      conditional_mean <- mean + latent_sd * rule$nodes[j]
      conditional_extrapolated <- mean_extrapolated + latent_sd * rule$nodes[j]
      if (probability) {
        fit <- .zplot_selnorm_threshold_summary(z, conditional_mean, sd, sei, selection, FALSE)
        ext <- .zplot_selnorm_threshold_summary(z, conditional_extrapolated, sd, sei, selection, TRUE)
        fitted[, 1L] <- fitted[, 1L] + rule$weights[j] * fit$EDR
        extrapolated[, 1L] <- extrapolated[, 1L] + rule$weights[j] * ext$EDR * ext$weights
        weights <- weights + rule$weights[j] * ext$weights
      } else {
        pair <- .zplot_selnorm_density_pair(z, conditional_mean, sd, sei, selection)
        weights <- weights + rule$weights[j] * rowMeans(
          .zplot_inverse_selection_weights(conditional_mean, sd, sei, selection)
        )
        if (!identical(mean, mean_extrapolated)) {
          pair$extrapolated <- .zplot_selnorm_density_matrix(
            z, conditional_extrapolated, sd, sei, selection, TRUE
          )
        }
        fitted       <- fitted + rule$weights[j] * pair$fitted
        extrapolated <- extrapolated + rule$weights[j] * pair$extrapolated
      }
    }
    current <- list(fitted = fitted, extrapolated = extrapolated,
                    weights = matrix(weights, S, 1L))
    if (!is.null(previous)) {
      error <- max(vapply(names(current), function(component) {
        value <- current[[component]]
        peak  <- apply(abs(value), 1L, max)
        change <- apply(abs(value - previous[[component]]), 1L, max)
        positive <- is.finite(peak) & peak > 0
        change[positive] <- change[positive] / peak[positive]
        change[!is.finite(change)] <- Inf
        max(change)
      }, numeric(1)))
      if (is.finite(error) && error <= control$relative_tolerance) {
        return(list(fitted = fitted, extrapolated = extrapolated, weights = weights,
                    EDR = if (probability) extrapolated[, 1L] / weights else NULL))
      }
    }
    previous <- current
  }
  stop(
    "Zplot marginal integration was rejected by diagnostics: relative quadrature change was ",
    format(error, digits = 4), ". Inspect the fitted selection and heterogeneity parameters.",
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


# Diagonal of Q - Q (Q + V)^-1 Q. The covariance factorization policy treats
# only eigensolver artifacts inside its backward-error envelope as null-space
# values; materially indefinite conditional covariances fail.
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
  chol_marginal <- tryCatch(
    chol(marginal_covariance),
    error = function(e) NULL
  )
  if (is.null(chol_marginal)) {
    stop(
      "Cannot solve the estimate-depth marginal covariance; it is not positive definite.",
      call. = FALSE
    )
  }

  solved <- backsolve(
    chol_marginal,
    forwardsolve(t(chol_marginal), latent_covariance)
  )
  reduction           <- latent_covariance %*% solved
  conditional_variance <- diag(latent_covariance) - diag(reduction)
  if (all(is.finite(conditional_variance)) &&
      all(conditional_variance >= 0)) {
    return(unname(conditional_variance))
  }

  conditional_covariance <- latent_covariance - reduction
  conditional_factorization <- .covariance_factorization(
    conditional_covariance
  )
  if (!.covariance_is_positive_semidefinite(conditional_factorization)) {
    stop(
      "The estimate-depth conditional latent covariance is not positive semidefinite.",
      call. = FALSE
    )
  }

  conditional_factor <- .covariance_sampling_factor(
    conditional_factorization
  )
  if (is.null(conditional_factor)) {
    stop(
      "Cannot factor the estimate-depth conditional latent covariance.",
      call. = FALSE
    )
  }

  return(colSums(conditional_factor^2))
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
                                           selection_context, extrapolate) {

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
    PACKAGE = "RoBMA"
  ))
}

.zplot_selnorm_density_pair <- function(z_sequence, mean, sd, sei,
                                         selection_context) {

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
    PACKAGE = "RoBMA"
  ))
}

.zplot_density_pair <- function(object, z_sequence, max_samples,
                                conditioning_depth = "marginal",
                                integration_control = set_selection_likelihood_control()) {

  posterior_samples <- .get_posterior_samples(object[["fit"]])
  selected_ind      <- .thin_sample_rows(nrow(posterior_samples), max_samples)
  if (!is.null(selected_ind)) {
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }

  if (.zplot_requires_selection_marginal(object, conditioning_depth)) {
    return(.zplot_selection_marginal(
      object, posterior_samples, z_sequence, NULL,
      conditioning_depth, integration_control
    ))
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
# Expected attempted studies represented by each observed study under the
# selection model. Normal/no-bias branches have weight one.
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
