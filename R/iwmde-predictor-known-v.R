# ============================================================================ #
# IWMDE Known-V Predictor Kernels
# ============================================================================ #

.iwmde_log_lik_known_v_joint_sum_from_samples <- function(
    context, posterior_samples, active_setup, unit, data_hash = NULL) {

  conditioned_random_effects <- .iwmde_conditioned_random_effects_from_latent(
    context           = context,
    posterior_samples = posterior_samples,
    unit              = unit
  )
  setup <- .log_lik_posterior_setup(
    fit                        = context[["object"]][["fit"]],
    posterior_samples          = posterior_samples,
    data                       = context[["data"]],
    priors                     = active_setup[["priors"]],
    unit                       = unit,
    data_hash                  = data_hash,
    conditioned_random_effects = conditioned_random_effects
  )
  fast <- .iwmde_log_lik_known_v_joint_sum_common_shift(
    context = context,
    setup   = setup
  )
  if (!is.null(fast)) {
    return(fast)
  }

  return(.log_lik_known_v_joint_sum_from_setup(setup))
}


.iwmde_log_lik_known_v_joint_sum_from_evaluated_predictors <- function(
    context, active_setup, mu_samples, tau_within_samples,
    tau_between_samples = NULL, posterior_samples = NULL, unit = "estimate",
    data_hash = NULL, random_effects_conditioning = "none") {

  setup <- .log_lik_evaluated_setup(
    fit                         = context[["object"]][["fit"]],
    data                        = context[["data"]],
    priors                      = active_setup[["priors"]],
    unit                        = unit,
    data_hash                   = data_hash,
    mu_samples                  = mu_samples,
    tau_within_samples          = tau_within_samples,
    tau_between_samples         = tau_between_samples,
    posterior_samples           = posterior_samples,
    random_effects_conditioning = random_effects_conditioning
  )
  fast <- .iwmde_log_lik_known_v_joint_sum_common_shift(
    context = context,
    setup   = setup
  )
  if (!is.null(fast)) {
    return(fast)
  }

  return(.log_lik_known_v_joint_sum_from_setup(setup))
}


.iwmde_conditioned_random_effects_from_latent <- function(
    context, posterior_samples, unit) {

  data <- context[["data"]]
  if (!identical(unit, "estimate") || !.is_data_random(data)) {
    return(NULL)
  }

  plan <- .iwmde_scalar_random_effect_plan(context)
  if (is.null(plan)) {
    return(NULL)
  }

  S   <- nrow(posterior_samples)
  K   <- nrow(data[["outcome"]])
  out <- NULL

  for (block in plan) {
    contribution <- .iwmde_scalar_random_effect_from_plan(
      block             = block,
      posterior_samples = posterior_samples,
      K                 = K
    )
    if (is.null(contribution)) {
      return(NULL)
    }
    out <- if (is.null(out)) contribution else out + contribution
  }

  if (is.null(out)) {
    out <- matrix(0, nrow = S, ncol = K)
  }
  return(out)
}


.iwmde_scalar_random_effect_plan <- function(context) {

  cache <- context[["predictor_cache"]]
  key   <- "scalar_random_effect_plan"
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  formula_design <- .fitted_formula_design(
    context[["object"]],
    "mu",
    required = TRUE
  )
  terms <- .formula_design_random_effects_by_mode(formula_design, "sampled")
  K     <- nrow(context[["data"]][["outcome"]])
  plan  <- lapply(terms, .iwmde_scalar_random_effect_block_plan, K = K)
  if (any(vapply(plan, is.null, logical(1)))) {
    plan <- NULL
  }

  assign(key, plan, envir = cache)
  return(plan)
}


.iwmde_scalar_random_effect_block_plan <- function(term, K) {

  if (!identical(term[["compile_mode"]], "sampled") ||
      !identical(term[["parameterization_resolved"]], "noncentered") ||
      !identical(term[["n_columns"]], 1L) ||
      !is.null(term[["latent_layout"]]) ||
      .random_effect_term_has_known_group_covariance(term)) {
    return(NULL)
  }

  model_matrix <- term[["model_matrix"]]
  group_map    <- term[["group_map"]]
  n_groups     <- term[["n_groups"]]
  sd_name      <- term[["sd_parameter_names"]]
  stem         <- term[["parameter_stem"]]
  if (!is.numeric(model_matrix) || !identical(dim(model_matrix), c(K, 1L)) ||
      any(!is.finite(model_matrix)) || !is.numeric(group_map) ||
      length(group_map) != K || anyNA(group_map) ||
      any(!is.finite(group_map)) || any(group_map != as.integer(group_map)) ||
      !is.numeric(n_groups) || length(n_groups) != 1L || is.na(n_groups) ||
      !is.finite(n_groups) || n_groups != as.integer(n_groups) ||
      n_groups < 1L || length(term[["group_levels"]]) != n_groups ||
      any(group_map < 1L) ||
      any(group_map > n_groups) || !is.character(sd_name) ||
      length(sd_name) != 1L || is.na(sd_name) || !nzchar(sd_name) ||
      !is.character(stem) || length(stem) != 1L || is.na(stem) ||
      !nzchar(stem)) {
    return(NULL)
  }

  z_names <- paste0(stem, "_xRE_Zx[", seq_len(n_groups), ",1]")
  return(list(
    model_values      = model_matrix[, 1L],
    unit_model        = all(model_matrix[, 1L] == 1),
    group_map         = as.integer(group_map),
    sd_name           = sd_name,
    z_names           = z_names,
    coefficient_names = paste0(
      stem,
      "_xRE_COEFx[",
      seq_len(n_groups),
      ",1]"
    )
  ))
}


.iwmde_scalar_random_effect_from_plan <- function(block, posterior_samples, K) {

  if (length(block[["model_values"]]) != K ||
      length(block[["group_map"]]) != K) {
    return(NULL)
  }

  sd_name           <- block[["sd_name"]]
  z_names           <- block[["z_names"]]
  coefficient_names <- block[["coefficient_names"]]
  if (!sd_name %in% colnames(posterior_samples) ||
      !all(z_names %in% colnames(posterior_samples)) ||
      any(coefficient_names %in% colnames(posterior_samples))) {
    return(NULL)
  }

  sd_samples <- posterior_samples[, sd_name]
  z_samples  <- posterior_samples[, z_names, drop = FALSE]
  if (any(!is.finite(sd_samples)) || any(sd_samples < 0) ||
      any(!is.finite(z_samples))) {
    return(NULL)
  }

  group_effects <- z_samples * sd_samples
  contribution  <- group_effects[, block[["group_map"]], drop = FALSE]
  if (!isTRUE(block[["unit_model"]])) {
    contribution <- sweep(
      contribution,
      2L,
      block[["model_values"]],
      "*"
    )
  }
  dimnames(contribution) <- NULL

  return(contribution)
}


.iwmde_log_lik_known_v_joint_sum_common_shift <- function(context, setup) {

  if (!identical(setup[["outcome_type"]], "norm") ||
      isTRUE(setup[["is_weightfunction"]]) ||
      !is.null(setup[["weights"]])) {
    return(NULL)
  }

  known_V <- .data_known_v_data(context[["data"]])
  K       <- setup[["K"]]
  S       <- setup[["S"]]
  mu      <- setup[["mu"]]
  if (.known_v_nrow(known_V) != K || !identical(dim(mu), c(S, K))) {
    return(NULL)
  }

  block_data        <- .known_v_dependency_block_data(context[["data"]], K)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  if (any(vapply(block_covariances, function(covariance) {
    !is.null(.covariance_exact_rank_one_factor(covariance)) &&
      nrow(covariance) > 1L
  }, logical(1)))) {
    return(NULL)
  }

  extra_variance <- .known_v_extra_variance_from_setup(setup)
  yi             <- setup[["yi"]]
  if (all(lengths(block_indices) == 1L)) {
    sampling_variance <- numeric(K)
    for (block in seq_along(block_indices)) {
      sampling_variance[block_indices[[block]]] <-
        block_covariances[[block]][1L, 1L]
    }
    variance <- sweep(extra_variance, 2L, sampling_variance, "+")
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      return(NULL)
    }
    residual <- matrix(yi, nrow = S, ncol = K, byrow = TRUE) - mu

    return(rowSums(-.5 * (
      log(2 * pi * variance) + residual^2 / variance
    )))
  }

  common_shift <- vapply(block_indices, function(index) {
    extra <- extra_variance[, index, drop = FALSE]

    length(index) == 1L || all(extra == extra[, 1L])
  }, logical(1))
  if (!all(common_shift)) {
    return(NULL)
  }

  spectral <- .iwmde_known_v_location_spectral_blocks(
    context    = context,
    block_data = block_data
  )
  log_lik <- numeric(S)
  for (block in seq_along(block_indices)) {
    index      <- block_indices[[block]]
    covariance <- block_covariances[[block]]
    residual <- matrix(
      yi[index],
      nrow  = S,
      ncol  = length(index),
      byrow = TRUE
    ) - mu[, index, drop = FALSE]
    shift <- extra_variance[, index[[1L]]]

    if (length(index) == 1L) {
      variance <- covariance[1L, 1L] + shift
      if (any(!is.finite(variance)) || any(variance <= 0)) {
        return(NULL)
      }
      log_lik <- log_lik - .5 * (
        log(2 * pi * variance) + residual[, 1L]^2 / variance
      )
      next
    }

    spectral_values <- spectral[[block]][["values"]]
    invariant_shift <- S == 1L || all(shift == shift[[1L]])
    denominator <- if (invariant_shift) {
      spectral_values + shift[[1L]]
    } else {
      outer(shift, spectral_values, "+")
    }
    if (any(!is.finite(denominator)) || any(denominator <= 0)) {
      return(NULL)
    }
    residual_rotated <- residual %*% spectral[[block]][["vectors"]]
    if (invariant_shift) {
      quadratic <- rowSums(sweep(
        residual_rotated^2,
        2L,
        denominator,
        "/"
      ))
      log_lik <- log_lik - .5 * (
        length(index) * log(2 * pi) + sum(log(denominator)) + quadratic
      )
    } else {
      log_lik <- log_lik - .5 * (
        length(index) * log(2 * pi) +
          rowSums(log(denominator) + residual_rotated^2 / denominator)
      )
    }
  }

  if (any(!is.finite(log_lik))) {
    return(NULL)
  }

  return(log_lik)
}


# Evaluate known-V tau candidates while reusing each block factorization across
# posterior rows.
.iwmde_log_q_grid_normal_known_v_tau_group <- function(
    context, parameter, values, row_states, replacement, setup, basis) {

  data <- context[["data"]]
  if (!identical(parameter, "tau") ||
      !.iwmde_uses_known_v_joint_likelihood(context) ||
      .is_data_random(data) ||
      .is_data_multilevel(data) ||
      isTRUE(setup[["is_weightfunction"]]) ||
      !is.null(setup[["weights"]]) ||
      !identical(basis[["scale_update"]], "tau") ||
      !is.null(basis[["mu_basis"]]) ||
      !is.null(basis[["log_tau_basis"]]) ||
      isTRUE(basis[["formula_mu"]]) ||
      isTRUE(basis[["formula_logtau"]])) {
    return(NULL)
  }

  known_V <- .data_known_v_data(data)
  G       <- length(values)
  S       <- length(row_states)
  K       <- setup[["K"]]
  if (.known_v_nrow(known_V) != K ||
      !identical(dim(setup[["mu"]]), c(S, K))) {
    return(NULL)
  }

  log_prior <- .iwmde_predictor_log_prior(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  if (is.null(log_prior) || length(log_prior) != G * S) {
    return(NULL)
  }

  block_data        <- .known_v_dependency_block_data(data, K)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  if (any(vapply(
    block_covariances,
    function(covariance) {
      !is.null(.covariance_exact_rank_one_factor(covariance)) &&
        nrow(covariance) > 1L
    },
    logical(1)
  ))) {
    return(NULL)
  }

  residual <- matrix(
    setup[["yi"]],
    nrow  = S,
    ncol  = K,
    byrow = TRUE
  ) - setup[["mu"]]
  log_lik <- matrix(-Inf, nrow = G, ncol = S)
  valid   <- is.finite(values) & values >= 0

  for (g in which(valid)) {
    candidate_log_lik <- numeric(S)
    for (block in seq_along(block_indices)) {
      idx        <- block_indices[[block]]
      covariance <- block_covariances[[block]] +
        diag(values[[g]]^2, nrow = length(idx))
      chol_covariance <- .known_v_chol_covariance(
        covariance = covariance,
        context    = "joint likelihood"
      )
      residual_whitened <- backsolve(
        chol_covariance,
        t(residual[, idx, drop = FALSE]),
        transpose = TRUE
      )
      candidate_log_lik <- candidate_log_lik - .5 * (
        length(idx) * log(2 * pi) +
          2 * sum(log(diag(chol_covariance))) +
          colSums(residual_whitened^2)
      )
    }
    log_lik[g, ] <- candidate_log_lik
  }

  return(log_lik + matrix(log_prior, nrow = G, ncol = S))
}


.iwmde_normal_location_likelihood_change_known_v <- function(context, setup,
                                                              yi, mu,
                                                              mu_basis) {

  if (isTRUE(setup[["is_weightfunction"]]) ||
      !is.null(setup[["weights"]])) {
    return(NULL)
  }

  known_V <- .data_known_v_data(context[["data"]])
  K       <- setup[["K"]]
  S       <- setup[["S"]]
  if (.known_v_nrow(known_V) != K ||
      !identical(dim(mu), c(S, K)) ||
      !identical(dim(mu_basis), c(S, K))) {
    return(NULL)
  }

  block_data        <- .known_v_dependency_block_data(context[["data"]], K)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  if (any(vapply(
    block_covariances,
    function(covariance) {
      !is.null(.covariance_exact_rank_one_factor(covariance)) &&
        nrow(covariance) > 1L
    },
    logical(1)
  ))) {
    return(NULL)
  }

  extra_variance <- .known_v_extra_variance_from_setup(setup)
  structured <- .iwmde_normal_location_likelihood_change_known_v_structured(
    context           = context,
    yi                = yi,
    mu                = mu,
    mu_basis          = mu_basis,
    block_data        = block_data,
    extra_variance    = extra_variance
  )
  if (!is.null(structured)) {
    return(structured)
  }

  return(.iwmde_normal_location_likelihood_change_known_v_cholesky(
    yi             = yi,
    mu             = mu,
    mu_basis       = mu_basis,
    block_data     = block_data,
    extra_variance = extra_variance
  ))
}


.iwmde_normal_location_likelihood_change_known_v_structured <- function(
    context, yi, mu, mu_basis, block_data, extra_variance) {

  S                 <- nrow(mu)
  K                 <- ncol(mu)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")

  if (all(lengths(block_indices) == 1L)) {
    sampling_variance <- numeric(K)
    for (block in seq_along(block_indices)) {
      sampling_variance[block_indices[[block]]] <-
        block_covariances[[block]][1L, 1L]
    }
    variance <- sweep(extra_variance, 2L, sampling_variance, "+")
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      return(NULL)
    }

    residual  <- matrix(yi, nrow = S, ncol = K, byrow = TRUE) - mu
    linear    <- rowSums(residual * mu_basis / variance)
    quadratic <- rowSums(mu_basis^2 / variance)

    return(.iwmde_normal_location_likelihood_change_result(
      linear    = linear,
      quadratic = quadratic
    ))
  }

  common_shift <- vapply(block_indices, function(index) {
    extra <- extra_variance[, index, drop = FALSE]

    length(index) == 1L || all(extra == extra[, 1L])
  }, logical(1))
  if (!all(common_shift)) {
    return(NULL)
  }

  spectral <- .iwmde_known_v_location_spectral_blocks(
    context    = context,
    block_data = block_data
  )
  linear    <- numeric(S)
  quadratic <- numeric(S)
  for (block in seq_along(block_indices)) {
    index      <- block_indices[[block]]
    covariance <- block_covariances[[block]]
    residual <- matrix(
      yi[index],
      nrow  = S,
      ncol  = length(index),
      byrow = TRUE
    ) - mu[, index, drop = FALSE]
    basis <- mu_basis[, index, drop = FALSE]
    shift <- extra_variance[, index[[1L]]]

    if (length(index) == 1L) {
      variance <- covariance[1L, 1L] + shift
      if (any(!is.finite(variance)) || any(variance <= 0)) {
        return(NULL)
      }
      linear    <- linear + residual[, 1L] * basis[, 1L] / variance
      quadratic <- quadratic + basis[, 1L]^2 / variance
      next
    }

    denominator <- sweep(
      matrix(
        spectral[[block]][["values"]],
        nrow  = S,
        ncol  = length(index),
        byrow = TRUE
      ),
      1L,
      shift,
      "+"
    )
    if (any(!is.finite(denominator)) || any(denominator <= 0)) {
      return(NULL)
    }
    residual_rotated <- residual %*% spectral[[block]][["vectors"]]
    basis_rotated    <- basis %*% spectral[[block]][["vectors"]]
    linear <- linear + rowSums(
      residual_rotated * basis_rotated / denominator
    )
    quadratic <- quadratic + rowSums(
      basis_rotated^2 / denominator
    )
  }

  return(.iwmde_normal_location_likelihood_change_result(
    linear    = linear,
    quadratic = quadratic
  ))
}


.iwmde_known_v_location_spectral_blocks <- function(context, block_data) {

  key <- .iwmde_predictor_cache_key(
    prefix      = "known_v_location_spectral_blocks",
    block_index = lapply(block_data, `[[`, "index"),
    covariance  = lapply(block_data, `[[`, "covariance")
  )
  cache <- context[["predictor_cache"]]
  if (is.environment(cache) &&
      exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  spectral <- lapply(block_data, function(block) {
    if (length(block[["index"]]) == 1L) {
      return(NULL)
    }

    eigen(block[["covariance"]], symmetric = TRUE)
  })
  if (is.environment(cache)) {
    assign(key, spectral, envir = cache)
  }

  return(spectral)
}


.iwmde_normal_location_likelihood_change_known_v_cholesky <- function(
    yi, mu, mu_basis, block_data, extra_variance) {

  S                 <- nrow(mu)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  linear    <- numeric(S)
  quadratic <- numeric(S)

  for (block in seq_along(block_indices)) {
    idx        <- block_indices[[block]]
    covariance <- block_covariances[[block]]
    residual   <- matrix(yi[idx], nrow = S, ncol = length(idx), byrow = TRUE) -
      mu[, idx, drop = FALSE]
    basis      <- mu_basis[, idx, drop = FALSE]
    extra      <- extra_variance[, idx, drop = FALSE]

    invariant_extra <- S == 1L || all(extra == matrix(
      extra[1L, ],
      nrow  = S,
      ncol  = length(idx),
      byrow = TRUE
    ))
    if (invariant_extra) {
      covariance <- covariance + diag(extra[1L, ], nrow = length(idx))
      chol_covariance <- .known_v_chol_covariance(
        covariance = covariance,
        context    = "joint likelihood"
      )
      residual_whitened <- backsolve(
        chol_covariance,
        t(residual),
        transpose = TRUE
      )
      basis_whitened <- backsolve(
        chol_covariance,
        t(basis),
        transpose = TRUE
      )
      linear    <- linear + colSums(residual_whitened * basis_whitened)
      quadratic <- quadratic + colSums(basis_whitened^2)
      next
    }

    for (s in seq_len(S)) {
      row_covariance <- covariance +
        diag(extra[s, ], nrow = length(idx))
      chol_covariance <- .known_v_chol_covariance(
        covariance = row_covariance,
        context    = "joint likelihood"
      )
      residual_whitened <- backsolve(
        chol_covariance,
        residual[s, ],
        transpose = TRUE
      )
      basis_whitened <- backsolve(
        chol_covariance,
        basis[s, ],
        transpose = TRUE
      )
      linear[[s]]    <- linear[[s]] +
        sum(residual_whitened * basis_whitened)
      quadratic[[s]] <- quadratic[[s]] + sum(basis_whitened^2)
    }
  }

  return(.iwmde_normal_location_likelihood_change_result(
    linear    = linear,
    quadratic = quadratic
  ))
}


.iwmde_normal_location_likelihood_change_result <- function(linear,
                                                             quadratic) {

  if (any(!is.finite(linear)) || any(!is.finite(quadratic)) ||
      any(quadratic < 0)) {
    return(NULL)
  }

  return(list(
    linear    = linear,
    quadratic = quadratic
  ))
}
