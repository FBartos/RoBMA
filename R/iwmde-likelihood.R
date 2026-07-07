# ============================================================================ #
# IWMDE Likelihood and Prior Ordinates
# ============================================================================ #

.iwmde_log_lik_posterior_setup_active_branch <- function(context,
                                                         posterior_samples,
                                                         active_setup,
                                                         unit,
                                                         data_hash = NULL) {

  posterior_samples <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = posterior_samples,
    active_setup = active_setup
  )

  return(.log_lik_posterior_setup(
    fit                  = context[["object"]][["fit"]],
    posterior_samples    = posterior_samples,
    data                 = context[["data"]],
    priors               = active_setup[["priors"]],
    unit                 = unit,
    data_hash            = data_hash
  ))
}


.iwmde_log_lik_from_posterior_samples_active_branch <- function(context,
                                                                posterior_samples,
                                                                active_setup,
                                                                unit,
                                                                add_metadata = FALSE,
                                                                data_hash = NULL) {

  posterior_samples <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = posterior_samples,
    active_setup = active_setup
  )

  return(.log_lik_from_posterior_samples(
    fit                  = context[["object"]][["fit"]],
    posterior_samples    = posterior_samples,
    data                 = context[["data"]],
    priors               = active_setup[["priors"]],
    unit                 = unit,
    add_metadata         = add_metadata,
    data_hash            = data_hash
  ))
}


.iwmde_log_lik_from_posterior_samples_sum_active_branch <- function(context,
                                                                    posterior_samples,
                                                                    active_setup,
                                                                    unit,
                                                                    data_hash = NULL) {

  posterior_samples <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = posterior_samples,
    active_setup = active_setup
  )

  if (.iwmde_uses_known_v_joint_likelihood(context)) {
    return(.log_lik_known_v_joint_sum_from_posterior_samples(
      fit               = context[["object"]][["fit"]],
      posterior_samples = posterior_samples,
      data              = context[["data"]],
      priors            = active_setup[["priors"]],
      unit              = unit,
      data_hash         = data_hash
    ))
  }

  return(.log_lik_from_posterior_samples_sum(
    fit                  = context[["object"]][["fit"]],
    posterior_samples    = posterior_samples,
    data                 = context[["data"]],
    priors               = active_setup[["priors"]],
    unit                 = unit,
    data_hash            = data_hash
  ))
}


.iwmde_log_lik_from_evaluated_predictors_sum_active_branch <- function(context,
                                                                       active_setup,
                                                                       mu_samples,
                                                                       tau_within_samples,
                                                                       tau_between_samples = NULL,
                                                                       posterior_samples = NULL,
                                                                       unit = "estimate",
                                                                       data_hash = NULL) {

  if (!is.null(posterior_samples)) {
    posterior_samples <- .iwmde_likelihood_posterior_samples(
      context      = context,
      samples      = posterior_samples,
      active_setup = active_setup
    )
  }

  if (.iwmde_uses_known_v_joint_likelihood(context)) {
    return(.log_lik_known_v_joint_sum_from_evaluated_predictors(
      fit                 = context[["object"]][["fit"]],
      data                = context[["data"]],
      priors              = active_setup[["priors"]],
      mu_samples          = mu_samples,
      tau_within_samples  = tau_within_samples,
      tau_between_samples = tau_between_samples,
      posterior_samples   = posterior_samples,
      unit                = unit,
      data_hash           = data_hash
    ))
  }

  return(.log_lik_from_evaluated_predictors_sum(
    fit                  = context[["object"]][["fit"]],
    data                 = context[["data"]],
    priors               = active_setup[["priors"]],
    mu_samples           = mu_samples,
    tau_within_samples   = tau_within_samples,
    tau_between_samples  = tau_between_samples,
    posterior_samples    = posterior_samples,
    unit                 = unit,
    data_hash            = data_hash
  ))
}


.iwmde_uses_known_v_joint_likelihood <- function(context) {

  .known_v_estimate_target_uses_backend(context[["data"]])
}


.iwmde_selection_context_active_branch <- function(context, active_setup,
                                                   posterior_samples,
                                                   newdata = NULL) {

  posterior_samples <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = posterior_samples,
    active_setup = active_setup
  )

  return(.selection_context_from_parts(
    fit                  = context[["object"]][["fit"]],
    data                 = context[["data"]],
    priors               = active_setup[["priors"]],
    posterior_samples    = posterior_samples,
    effect_direction     = .data_effect_direction(context[["data"]]),
    newdata              = newdata
  ))
}


.iwmde_log_q_state <- function(context, row, active_setup, parameters,
                               prior_list, likelihood_mode) {

  likelihood_row <- if (identical(likelihood_mode, "marginal") &&
                        !.iwmde_marginal_likelihood_requires_row(context) &&
                        !("gamma" %in% names(parameters)) &&
                        !("theta" %in% names(parameters)) &&
                        !("pi" %in% names(parameters)) &&
                        !("phi" %in% names(parameters))) {
    NULL
  } else {
    row
  }
  log_lik <- .iwmde_log_likelihood_parameters(
    context         = context,
    parameters      = parameters,
    active_setup    = active_setup,
    likelihood_mode = likelihood_mode,
    row             = likelihood_row
  )
  if (!is.finite(log_lik)) {
    return(-Inf)
  }

  log_prior <- .iwmde_log_prior_row(row, prior_list)
  out       <- log_lik + log_prior
  if (!is.finite(out)) {
    return(-Inf)
  }

  return(out)
}


.iwmde_log_likelihood_parameters <- function(context, parameters,
                                             active_setup,
                                             likelihood_mode = "conditional",
                                             row = NULL) {

  if (identical(likelihood_mode, "marginal")) {
    return(.iwmde_log_likelihood_parameters_marginal(
      context      = context,
      parameters   = parameters,
      active_setup = active_setup,
      row          = row
    ))
  }

  log_lik <- .log_posterior(
    parameters        = parameters,
    data              = active_setup[["fit_data"]],
    is_mods           = .is_data_mods(context[["data"]]),
    is_scale          = .is_data_scale(context[["data"]]),
    is_random         = .is_data_random(context[["data"]]),
    is_multilevel     = .is_data_multilevel(context[["data"]]),
    is_weights        = .is_data_weights(context[["data"]]),
    is_known_v        = .is_data_known_v(context[["data"]]),
    known_v_parameterization = .data_known_v_parameterization(context[["data"]]),
    model_data        = context[["data"]],
    is_PET            = active_setup[["is_PET"]],
    is_PEESE          = active_setup[["is_PEESE"]],
    is_weightfunction = active_setup[["is_weightfunction"]],
    effect_direction  = .data_effect_direction(context[["data"]]),
    outcome_type      = .data_outcome_type(context[["data"]])
  )
  if (length(log_lik) != 1L || !is.finite(log_lik)) {
    return(-Inf)
  }

  return(as.numeric(log_lik))
}


.iwmde_log_likelihood_parameters_marginal <- function(context, parameters,
                                                      active_setup,
                                                      row = NULL) {

  log_lik <- if (!is.null(row)) {
    .iwmde_log_likelihood_row_marginal(
      context      = context,
      row          = row,
      active_setup = active_setup
    )
  } else {
    .iwmde_log_likelihood_parameters_marginal_internal(
      context      = context,
      parameters   = parameters,
      active_setup = active_setup
    )
  }
  if (length(log_lik) != 1L || !is.finite(log_lik)) {
    return(-Inf)
  }

  return(as.numeric(log_lik))
}


.iwmde_log_likelihood_row_marginal <- function(context, row, active_setup) {

  samples <- matrix(
    as.numeric(row[colnames(context[["posterior_samples"]])]),
    nrow = 1L
  )
  colnames(samples) <- colnames(context[["posterior_samples"]])
  unit <- if (.is_data_multilevel(context[["data"]])) {
    "cluster"
  } else {
    "estimate"
  }

  log_lik <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(
    context           = context,
    posterior_samples = samples,
    active_setup      = active_setup,
    unit              = unit
  )

  return(log_lik)
}


.iwmde_log_likelihood_parameters_marginal_internal <- function(context,
                                                               parameters,
                                                               active_setup) {

  data         <- context[["data"]]
  outcome_type <- .data_outcome_type(data)
  K            <- nrow(data[["outcome"]])

  mu_samples <- .marglik_get_mu_samples(
    parameters       = parameters,
    is_mods          = .is_data_mods(data),
    is_PET           = active_setup[["is_PET"]],
    is_PEESE         = active_setup[["is_PEESE"]],
    effect_direction = .data_effect_direction(data),
    sei              = if (outcome_type == "norm") data[["outcome"]][["sei"]] else NULL,
    K                = K
  )
  tau_result <- .marglik_get_tau_samples(
    parameters    = parameters,
    is_scale      = .is_data_scale(data),
    is_multilevel = .is_data_multilevel(data),
    K             = K,
    fixed_tau     = .fixed_tau_prior_value(active_setup[["priors"]])
  )

  if (!outcome_type %in% c("norm", "bin", "pois")) {
    return(-Inf)
  }

  posterior_samples <- .iwmde_parameters_posterior_samples(parameters)
  unit <- if (.is_data_multilevel(data)) {
    "cluster"
  } else {
    "estimate"
  }

  return(.iwmde_log_lik_from_evaluated_predictors_sum_active_branch(
    context             = context,
    active_setup        = active_setup,
    mu_samples          = mu_samples,
    tau_within_samples  = tau_result[["tau_within"]],
    tau_between_samples = tau_result[["tau_between"]],
    posterior_samples   = posterior_samples,
    unit                = unit
  ))
}


.iwmde_parameters_posterior_samples <- function(parameters) {

  .parameters_as_sample_matrix(parameters)
}


.iwmde_log_prior_row <- function(row, prior_list) {

  if (length(prior_list) == 0L) {
    return(0)
  }

  row <- .resolve_fixed_prior_row(
    row        = row,
    prior_list = prior_list
  )
  log_prior <- BayesTools::JAGS_marglik_priors(row, prior_list)
  if (length(log_prior) != 1L || !is.finite(log_prior)) {
    return(-Inf)
  }

  return(as.numeric(log_prior))
}


.iwmde_focal_log_prior <- function(prior, value, parameter = NULL) {

  out <- .iwmde_focal_log_prior_values(
    prior     = prior,
    values    = value,
    parameter = parameter
  )
  if (length(out) != 1L || is.na(out[[1L]])) {
    return(NA_real_)
  }

  return(as.numeric(out[[1L]]))
}


.iwmde_focal_log_prior_values <- function(prior, values, parameter = NULL) {

  if (is.null(prior) || BayesTools::is.prior.none(prior)) {
    return(rep(0, length(values)))
  }

  fast <- .iwmde_focal_log_prior_values_fast(
    prior     = prior,
    values    = values,
    parameter = parameter
  )
  if (!is.null(fast)) {
    return(fast)
  }

  return(vapply(values, function(value) {
    out <- tryCatch(
      BayesTools::lpdf(prior, value),
      error = function(e) NA_real_
    )
    if (length(out) != 1L || is.na(out[[1L]])) {
      return(NA_real_)
    }
    return(as.numeric(out[[1L]]))
  }, numeric(1)))
}


.iwmde_focal_log_prior_values_fast <- function(prior, values,
                                               parameter = NULL) {

  if (inherits(prior, "prior.simple") &&
      identical(prior[["distribution"]], "normal")) {
    parameters <- prior[["parameters"]]
    out <- stats::dnorm(
      values,
      mean = parameters[["mean"]],
      sd   = parameters[["sd"]],
      log  = TRUE
    )
    truncation <- prior[["truncation"]]
    lower      <- if (!is.null(truncation[["lower"]])) {
      truncation[["lower"]]
    } else {
      -Inf
    }
    upper      <- if (!is.null(truncation[["upper"]])) {
      truncation[["upper"]]
    } else {
      Inf
    }
    log_normalizer <- .iwmde_log_normal_interval_prob(
      lower = lower,
      upper = upper,
      mean  = parameters[["mean"]],
      sd    = parameters[["sd"]]
    )
    if (!is.finite(log_normalizer)) {
      return(NULL)
    }
    out <- out - log_normalizer
    out[values < lower | values > upper] <- -Inf
    return(out)
  }

  if (inherits(prior, "prior.vector") &&
      identical(prior[["distribution"]], "mnormal") &&
      .iwmde_prior_has_default_truncation(prior)) {
    parameters <- prior[["parameters"]]
    mean       <- .iwmde_prior_coordinate_parameter(
      x         = parameters[["mean"]],
      parameter = parameter
    )
    sd         <- .iwmde_prior_coordinate_parameter(
      x         = parameters[["sd"]],
      parameter = parameter
    )
    if (length(mean) == 1L && length(sd) == 1L &&
        is.finite(mean) && is.finite(sd) && sd > 0) {
      return(stats::dnorm(
        values,
        mean = mean,
        sd   = sd,
        log  = TRUE
      ))
    }
  }

  return(NULL)
}


.iwmde_can_use_focal_prior_delta <- function(prior) {

  if (is.null(prior) || BayesTools::is.prior.none(prior)) {
    return(TRUE)
  }
  if (.iwmde_prior_is_invgamma(prior)) {
    return(FALSE)
  }

  return(inherits(prior, "prior.simple"))
}


.iwmde_log_normal_interval_prob <- function(lower, upper, mean, sd) {

  if (!is.finite(sd) || sd <= 0 || lower >= upper) {
    return(-Inf)
  }

  log_cdf_upper <- stats::pnorm(upper, mean = mean, sd = sd, log.p = TRUE)
  log_cdf_lower <- stats::pnorm(lower, mean = mean, sd = sd, log.p = TRUE)
  log_left      <- .iwmde_logspace_subtract(log_cdf_upper, log_cdf_lower)

  log_tail_lower <- stats::pnorm(
    lower, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE
  )
  log_tail_upper <- stats::pnorm(
    upper, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE
  )
  log_right <- .iwmde_logspace_subtract(log_tail_lower, log_tail_upper)
  candidates <- c(log_left, log_right)
  candidates <- candidates[is.finite(candidates)]
  if (length(candidates) == 0L) {
    return(-Inf)
  }

  return(max(candidates))
}


.iwmde_logspace_subtract <- function(log_x, log_y) {

  if (!is.finite(log_x)) {
    return(-Inf)
  }
  if (!is.finite(log_y)) {
    return(log_x)
  }
  if (log_y > log_x) {
    return(NaN)
  }

  return(log_x + log1p(-exp(log_y - log_x)))
}


.iwmde_prior_has_default_truncation <- function(prior) {

  truncation <- prior[["truncation"]]
  lower      <- truncation[["lower"]]
  upper      <- truncation[["upper"]]
  if (is.null(lower)) {
    lower <- -Inf
  }
  if (is.null(upper)) {
    upper <- Inf
  }

  return(identical(as.numeric(lower), -Inf) &&
    identical(as.numeric(upper), Inf))
}


.iwmde_prior_coordinate_parameter <- function(x, parameter = NULL) {

  x <- as.numeric(x)
  if (length(x) <= 1L || is.null(parameter)) {
    return(x[[1L]])
  }

  index <- .iwmde_parameter_index(parameter)
  if (is.na(index) || index < 1L || index > length(x)) {
    return(NA_real_)
  }

  return(x[[index]])
}


.iwmde_parameter_index <- function(parameter) {

  indexed <- regexec("^.+\\[([0-9]+)\\]$", parameter)
  match   <- regmatches(parameter, indexed)[[1]]
  if (length(match) != 2L) {
    return(NA_integer_)
  }

  return(as.integer(match[[2L]]))
}
