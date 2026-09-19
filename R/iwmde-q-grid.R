# ============================================================================ #
# IWMDE Q-Grid Dispatch and GLMM Batches
# ============================================================================ #

.iwmde_log_q_grid <- function(context, parameter, values, row_states,
                              replacement) {

  transformed <- vapply(row_states, function(state) {
    !is.null(state[["conditioning_transform"]])
  }, logical(1L))
  if (any(transformed)) {
    out <- matrix(NA_real_, length(values), length(row_states))
    out[, transformed] <- .iwmde_log_q_grid_retained_location(context, parameter,
      values, row_states[transformed], replacement)
    if (any(!transformed)) {
      ordinary <- .iwmde_log_q_grid(context, parameter, values,
        row_states[!transformed], replacement)
      out[, !transformed] <- ordinary
      change <- attr(ordinary, "max_quadrature_relative_change", exact = TRUE)
      if (!is.null(change)) attr(out, "max_quadrature_relative_change") <- change
    }
    return(out)
  }

  if (.is_data_joint_selection(context[["data"]])) {
    # An all-singleton selection plan has no dependent block to integrate, so
    # the batched predictor route evaluates the same joint density as the
    # generic candidate route below: the closed-form normal location change and
    # the native normalizer delta grid replace one full likelihood per
    # (value, row) candidate. Retained-location rows have already returned
    # above, and the generic route stays the fallback whenever the predictor
    # batch declines.
    plan <- .data_selection_execution_plan(context[["data"]])
    if (length(plan[["block_methods"]]) > 0L &&
        all(plan[["block_methods"]] == "singleton")) {
      predictor <- .iwmde_log_q_grid_predictor_batch(
        context     = context,
        parameter   = parameter,
        values      = values,
        row_states  = row_states,
        replacement = replacement
      )
      if (is.matrix(predictor) &&
          nrow(predictor) == length(values) &&
          ncol(predictor) == length(row_states)) {
        return(predictor)
      }
    }

    modes <- unique(vapply(row_states, `[[`, character(1L), "likelihood_mode"))
    if (length(modes) == 1L) {
      result <- .iwmde_log_q_grid_from_samples(
        context, parameter, values, row_states, replacement,
        likelihood_mode = modes,
        log_lik_fun = function(samples, active_setup, batch) {
          affine <- .iwmde_joint_affine_log_likelihood(context, parameter, values,
            samples, active_setup, batch, replacement)
          if (!is.null(affine)) return(affine)
          rank_one <- .iwmde_joint_rank_one_variance_log_likelihood(
            context, parameter, values, samples, active_setup, batch, replacement)
          if (!is.null(rank_one)) return(rank_one)
          .iwmde_log_lik_from_posterior_samples_sum_active_branch(
            context, samples, active_setup, unit = "estimate"
          )
        }
      )
      if (is.matrix(result)) return(result)
    }
    return(.iwmde_log_q_grid_scalar(context, parameter, values, row_states, replacement))
  }

  out <- .iwmde_log_q_grid_predictor_batch(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )

  if (is.matrix(out) &&
      nrow(out) == length(values) &&
      ncol(out) == length(row_states)) {
    return(out)
  }

  out <- .iwmde_log_q_grid_glmm_conditional_batch(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )

  if (is.matrix(out) &&
      nrow(out) == length(values) &&
      ncol(out) == length(row_states)) {
    return(out)
  }

  out <- .iwmde_log_q_grid_batch(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )

  if (is.matrix(out) &&
      nrow(out) == length(values) &&
      ncol(out) == length(row_states)) {
    return(out)
  }

  return(.iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  ))
}


.iwmde_log_q_grid_glmm_conditional_batch <- function(context, parameter, values,
                                                     row_states, replacement) {

  if (!.data_outcome_type(context[["data"]]) %in% c("bin", "pois")) {
    return(NULL)
  }

  return(.iwmde_log_q_grid_from_samples(
    context         = context,
    parameter       = parameter,
    values          = values,
    row_states      = row_states,
    replacement     = replacement,
    likelihood_mode = "conditional",
    log_lik_fun     = function(samples, active_setup, batch) {
      .iwmde_glmm_conditional_log_likelihood_samples(
        context      = context,
        samples      = samples,
        active_setup = active_setup
      )
    }
  ))
}


.iwmde_glmm_conditional_log_likelihood_samples <- function(context, samples,
                                                           active_setup) {

  data         <- context[["data"]]
  outcome_type <- .data_outcome_type(data)
  if (!outcome_type %in% c("bin", "pois")) {
    return(NULL)
  }

  K             <- nrow(data[["outcome"]])
  is_mods       <- .is_data_mods(data)
  is_scale      <- .is_data_scale(data)
  is_multilevel <- .is_data_multilevel(data)

  mu_samples <- .evaluate.brma.mu(
    fit               = context[["object"]][["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = active_setup[["priors"]][["mods"]],
    is_mods           = is_mods,
    is_PET            = active_setup[["is_PET"]],
    is_PEESE          = active_setup[["is_PEESE"]],
    effect_direction  = .data_effect_direction(data),
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = samples,
    priors            = active_setup[["priors"]]
  )
  tau_result <- .evaluate.brma.tau(
    fit               = context[["object"]][["fit"]],
    scale_data        = data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors      = active_setup[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = samples,
    fixed_tau         = .fixed_tau_prior_value(active_setup[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(active_setup[["priors"]])
  )

  return(.log_lik_glmm_conditional_sum_from_evaluated_predictors(
    fit                 = context[["object"]][["fit"]],
    data                = data,
    posterior_samples   = samples,
    mu_samples          = mu_samples,
    tau_within_samples  = tau_result[["tau_within"]],
    tau_between_samples = tau_result[["tau_between"]]
  ))
}


.iwmde_log_q_grid_predictor_batch <- function(context, parameter, values,
                                              row_states, replacement) {

  if (length(values) == 0L || length(row_states) == 0L) {
    return(matrix(-Inf, nrow = length(values), ncol = length(row_states)))
  }

  likelihood_modes <- vapply(row_states, function(state) {
    state[["likelihood_mode"]]
  }, character(1))
  if (!all(likelihood_modes == "marginal")) {
    return(NULL)
  }

  out    <- matrix(-Inf, nrow = length(values), ncol = length(row_states))
  groups <- .iwmde_row_state_groups(context, row_states)
  quadrature_change <- NA_real_

  for (state_cols in groups) {
    group_states <- row_states[state_cols]
    group_out    <- .iwmde_log_q_grid_predictor_group(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = group_states,
      replacement = replacement
    )

    if (!is.matrix(group_out)) {
      return(NULL)
    }

    out[, state_cols] <- group_out
    group_quadrature_change <- attr(
      group_out,
      "max_quadrature_relative_change",
      exact = TRUE
    )
    if (length(group_quadrature_change) == 1L &&
        is.finite(group_quadrature_change)) {
      quadrature_change <- max(
        quadrature_change,
        group_quadrature_change,
        na.rm = TRUE
      )
    }
  }

  if (is.finite(quadrature_change)) {
    attr(out, "max_quadrature_relative_change") <- quadrature_change
  }

  return(out)
}


.iwmde_log_q_grid_predictor_group <- function(context, parameter, values,
                                              row_states, replacement) {

  active_setup <- row_states[[1L]][["active_setup"]]
  unit         <- if (.is_data_multilevel(context[["data"]])) {
    "cluster"
  } else {
    "estimate"
  }
  random_scalar_out <-
    .iwmde_log_q_grid_known_v_random_group_iid(
      context      = context,
      parameter    = parameter,
      values       = values,
      row_states   = row_states,
      replacement  = replacement,
      active_setup = active_setup
    )
  if (is.matrix(random_scalar_out)) {
    return(random_scalar_out)
  }
  random_affine_out <- .iwmde_log_q_grid_known_v_random_affine(
    context      = context,
    parameter    = parameter,
    values       = values,
    row_states   = row_states,
    replacement  = replacement,
    active_setup = active_setup
  )
  if (is.matrix(random_affine_out)) {
    return(random_affine_out)
  }
  random_factor_out <- .iwmde_log_q_grid_known_v_random_factor(
    context      = context,
    parameter    = parameter,
    values       = values,
    row_states   = row_states,
    replacement  = replacement,
    active_setup = active_setup
  )
  if (is.matrix(random_factor_out)) {
    return(random_factor_out)
  }
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = unit
  )
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )
  if (is.null(basis)) {
    return(NULL)
  }
  normal_out <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  if (is.matrix(normal_out)) {
    return(normal_out)
  }
  if (.iwmde_uses_known_v_random_marginal_likelihood(
      context, priors = active_setup[["priors"]])) {
    return(NULL)
  }

  native <- .iwmde_predictor_normal_grid_log_lik(
    context      = context,
    active_setup = active_setup,
    setup        = setup,
    basis        = basis,
    values       = values,
    unit         = unit
  )
  if (!is.null(native)) {
    # Prior ordinates never depended on the candidate matrices for the batches
    # this route serves: no formula replacement means no replacement samples.
    log_prior <- .iwmde_predictor_log_prior(
      context             = context,
      parameter           = parameter,
      values              = values,
      row_states          = row_states,
      replacement         = replacement,
      replacement_samples = NULL
    )
    if (is.null(log_prior) || length(log_prior) != length(native[["log_lik"]])) {
      return(NULL)
    }
    log_q <- native[["log_lik"]] + log_prior
    log_q[!native[["valid"]]] <- -Inf

    return(matrix(log_q, nrow = length(values), ncol = length(row_states)))
  }

  candidates <- .iwmde_predictor_candidates(
    context     = context,
    active_setup = active_setup,
    setup       = setup,
    basis       = basis,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  posterior_samples <- candidates[["posterior_samples"]]
  log_prior <- .iwmde_predictor_log_prior(
    context             = context,
    parameter           = parameter,
    values              = values,
    row_states          = row_states,
    replacement         = replacement,
    replacement_samples = candidates[["replacement_samples"]]
  )
  if (is.null(log_prior) || length(log_prior) != nrow(candidates[["mu"]])) {
    return(NULL)
  }

  random_factor_samples <- NULL
  if (.setup_uses_joint_selection_likelihood(setup) &&
      is.null(candidates[["replacement_samples"]])) {
    # These candidates only expand existing posterior rows. Predictor changes
    # have already been evaluated, so their random-effect factors are invariant.
    random_factor_samples <- .selection_joint_random_factor_samples(setup)
    if (!is.null(random_factor_samples)) {
      rows <- candidates[["row_index"]]
      random_factor_samples[["diagonal"]] <-
        random_factor_samples[["diagonal"]][rows, , drop = FALSE]
      random_factor_samples[["loadings"]] <- lapply(
        random_factor_samples[["loadings"]],
        function(loading) loading[rows, , , drop = FALSE]
      )
    }
  }

  log_lik <- .iwmde_log_lik_from_evaluated_predictors_sum_active_branch(
    context               = context,
    active_setup          = active_setup,
    mu_samples            = candidates[["mu"]],
    tau_within_samples    = candidates[["tau_within"]],
    tau_between_samples   = candidates[["tau_between"]],
    posterior_samples     = posterior_samples,
    unit                  = unit,
    random_factor_samples = random_factor_samples
  )
  log_q <- log_lik + log_prior
  log_q[!candidates[["valid"]]] <- -Inf

  return(matrix(log_q, nrow = length(values), ncol = length(row_states)))
}


# Reuse a declared fixed-coefficient direction in the existing full joint setup.
# Candidate posterior rows still own all priors, selection and covariance inputs.
.iwmde_joint_affine_log_likelihood <- function(context, parameter, values,
    samples, active_setup, batch, replacement) {

  data <- context[["data"]]
  tracked <- FALSE
  on.exit({
    if (is.environment(context[["normalizer_grid"]]) && !tracked) context[["normalizer_grid"]]$untracked <- TRUE
  }, add = TRUE)
  if (!replacement[["type"]] %in% c("linear", "scalar") ||
      !.is_data_joint_selection(data) || !.is_data_known_v(data) ||
      .data_outcome_type(data) != "norm" || .selection_retains_sampling(data)) return(NULL)
  plan <- .data_selection_execution_plan(data)
  # Rebuilding only the changed predictor rows is a property of the mean
  # sweep, not of how a block normalizer is evaluated, so every certified
  # block route qualifies.
  if (!any(plan[["block_methods"]] %in% c("dense", "rank_one", "factor")) ||
      any(!plan[["block_methods"]] %in%
          c("dense", "rank_one", "factor", "singleton"))) return(NULL)
  # The static chart resolver rejects mean translations and dynamic formula
  # multipliers. Verify the full dependency contract at this evaluation seam.
  design <- .fitted_formula_design(context[["object"]], "mu", required = TRUE)
  if (any(vapply(design[["random_effects"]], function(term) {
    !is.null(term[["mean_translation"]])
  }, logical(1L)))) return(NULL)
  states <- batch[["row_states"]]
  columns <- if (identical(replacement[["type"]], "scalar")) {
    parameter
  } else {
    unique(unlist(lapply(states, function(state) {
      .iwmde_linear_replacement_state(context, state, replacement)[["active_columns"]]
    }), use.names = FALSE))
  }
  if (!length(columns)) return(NULL)
  dependencies <- BayesTools::JAGS_formula_coordinate_dependencies(
    context[["object"]][["fit"]], columns)
  if (any(dependencies[["formula_parameter"]] != "mu") ||
      any(dependencies[["dependency_type"]] != "coefficient")) return(NULL)
  states <- batch[["row_states"]]
  baseline <- .iwmde_predictor_setup(context, states, active_setup, "estimate")
  basis <- .iwmde_predictor_update_basis(context, parameter, states, replacement, baseline)
  if (is.null(basis) || !identical(basis[["scale_update"]], "none") ||
      isTRUE(basis[["formula_mu"]]) || !is.null(basis[["log_tau_basis"]]) ||
      !is.matrix(basis[["mu_basis"]]) ||
      !identical(dim(basis[["mu_basis"]]), dim(baseline[["mu"]]))) return(NULL)
  candidates <- batch[["candidates"]]
  positions <- batch[["valid_positions"]]
  state_index <- candidates[["state_index"]][positions]
  grid_index <- candidates[["grid_index"]][positions]
  delta <- values[grid_index] - basis[["current"]][state_index]
  mu <- baseline[["mu"]][state_index, , drop = FALSE]
  update <- basis[["mu_basis"]][state_index, , drop = FALSE] * delta
  if (any(!is.finite(update))) return(NULL)
  # Exactly unchanged rows retain the baseline bytes. This is algebraic zero
  # work, not a tolerance or a rounded-parameter lookup.
  changing <- update != 0
  mu[changing] <- mu[changing] + update[changing]
  if (any(!is.finite(mu))) return(NULL)
  setup <- .iwmde_affine_candidate_setup(context, baseline, samples,
    active_setup, state_index)
  if (!identical(dim(mu), dim(setup[["mu"]]))) {
    stop("Affine joint predictors do not match the candidate rows.", call. = FALSE)
  }
  setup[["mu"]] <- mu
  setup[["selection_static"]] <- .iwmde_selection_joint_static(context)
  if (is.environment(context[["normalizer_grid"]])) {
    setup[["normalizer_grid"]] <- list(shared = context[["normalizer_grid"]],
      rows = vapply(states, `[[`, integer(1L), "row_index"),
      current = basis[["current"]], basis = basis[["mu_basis"]],
      mean = baseline[["mu"]], state_index = state_index, values = values[grid_index])
  }
  result <- .log_lik_estimate_sum_from_setup(setup)
  tracked <- TRUE
  result
}


# The candidate rows of an affine mean sweep repeat one posterior row per state
# in every coordinate the sweep does not write, and the dependency guard above
# admits only mu coefficient coordinates. Every field of the likelihood setup
# other than the location is therefore the state's own field, so the setup is
# built once for the S states - it already is, as the cached predictor setup the
# affine direction is read from - and its rows are repeated for the candidates.
# The candidates keep their own posterior rows, so a consumer that reads a swept
# coordinate still sees the candidate's value.
.iwmde_affine_candidate_setup <- function(context, state_setup, samples,
                                          active_setup, state_index) {

  setup <- state_setup
  for (field in c("mu", "mu_random", "tau_total", "tau_within", "tau_between")) {
    value <- setup[[field]]
    if (is.matrix(value)) {
      setup[[field]] <- value[state_index, , drop = FALSE]
    }
  }
  if (!is.null(setup[["rho"]])) {
    setup[["rho"]] <- setup[["rho"]][state_index]
  }
  sources <- setup[["marginalized_random_source_samples"]]
  if (is.list(sources) && length(sources) > 0L) {
    setup[["marginalized_random_source_samples"]] <- lapply(sources, function(value) {
      if (is.matrix(value)) value[state_index, , drop = FALSE] else value
    })
  }
  setup[["S"]] <- length(state_index)
  setup[["posterior_samples"]] <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  setup[["state_rows"]] <- .iwmde_affine_state_rows(state_setup, state_index)

  return(setup)
}


# The state repetition a candidate batch carries, or NULL when the batch does
# not repeat whole states and every construction has to be evaluated per row.
.iwmde_affine_state_rows <- function(state_setup, state_index) {

  S <- state_setup[["S"]]
  if (!is.integer(state_index) || length(state_index) == 0L ||
      anyNA(state_index) || !is.numeric(S) || length(S) != 1L ||
      min(state_index) < 1L || max(state_index) > S) {
    return(NULL)
  }

  return(list(setup = state_setup, index = state_index))
}


# One joint-selection static per density-line context. The migrated execution
# plan and the caches its block constructions key by evaluated state count are
# functions of the fitted data alone, so a density line resolves them once for
# all of its replacement chunks instead of once per chunk.
.iwmde_selection_joint_static <- function(context) {

  cache <- context[["predictor_cache"]]
  key   <- "selection_joint_static"
  if (is.environment(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  static <- .selection_joint_static(context[["data"]])
  if (is.environment(cache)) {
    assign(key, static, envir = cache)
  }

  return(static)
}
