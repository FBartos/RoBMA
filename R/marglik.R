# ============================================================================ #
# marglik.R
# ============================================================================ #
#
# This file implements marginal likelihood computation for brma class objects
# via bridge sampling. The marginal likelihood is used for Bayesian model
# comparison via Bayes factors.
#
# The implementation reuses the evaluate and pdf helper functions from
# evaluate.R and pdf.R by converting single-sample parameters
# to 1-row matrices.
#
# ============================================================================ #


### RoBMA 4.0.0


# ---------------------------------------------------------------------------- #
# add_marglik generic and method
# ---------------------------------------------------------------------------- #

#' @export
add_marglik <- function(object, ...) UseMethod("add_marglik")


#' @title Add Marginal Likelihood to brma Objects
#'
#' @description Compute the marginal likelihood of a brma model using
#' bridge sampling and store the result in the object.
#'
#' @param object a brma model object.
#' @param ... additional arguments (currently not used).
#'
#' @details
#' The marginal likelihood is computed using the \code{bridgesampling} package
#' via the \code{BayesTools::JAGS_bridgesampling} wrapper. The result is stored
#' in the object and can be extracted using \code{\link{bridge_sampler.brma}}.
#'
#' Product-space model-averaging objects (\code{BMA.norm}, \code{BMA.glmm},
#' and \code{RoBMA}) do not expose a bridge-sampling marginal likelihood;
#' use predictive comparison methods such as \code{\link{loo.brma}} instead.
#' For \code{brma.mv()} known-\code{V} objects, bridge sampling evaluates the
#' joint likelihood corresponding to the fitted known-\code{V} backend, not the
#' conditional estimate-wise target used by LOO/WAIC diagnostics.
#'
#' @return The brma object with the marginal likelihood result stored in
#' \code{object[["marglik"]]}.
#'
#' @seealso \code{\link{bridge_sampler.brma}}, \code{\link{logml.brma}},
#' \code{\link{bf.brma}}, \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit <- add_marglik(fit)
#'
#'   bridge <- bridge_sampler(fit)
#'   print(bridge)
#'
#'   logml(fit)
#' }
#' }
#'
#' @aliases add_marglik
#' @export
add_marglik.brma <- function(object, ...) {

  .check_marglik_available(object, "add_marglik()")
  marglik <- .marglik(object)
  if (inherits(marglik, "error")) {
    stop(conditionMessage(marglik), call. = FALSE)
  }
  marglik <- .brma_mv_attach_marglik_target_metadata(marglik, object)
  object[["marglik"]] <- marglik
  return(object)
}


# ---------------------------------------------------------------------------- #
# .marglik internal function
# ---------------------------------------------------------------------------- #

# Compute marginal likelihood for brma objects using bridge sampling.
# This is an internal function called by add_marglik.brma.
#
# @param object a brma model object.
#
# @return An object of class "bridge" as returned by bridgesampling::bridge_sampler.
#
# @keywords internal
.marglik <- function(object) {

  data   <- object[["data"]]
  priors <- object[["priors"]]
  fit    <- object[["fit"]]

  .check_marglik_available(object, ".marglik()")

  ### create model base
  fit_priors       <- .create_fit_priors(data = data, priors = priors)
  fit_data         <- .create_fit_data(data = data, priors = priors)
  fit_data         <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = priors,
    effect_direction = .data_effect_direction(data)
  )
  fit_formula_args <- .create_jags_formula_args(data = data, priors = priors)
  bridge_sd_source_spec <- .marglik_bridge_sd_source_spec(
    add_parameters = fit_formula_args[["add_parameters"]],
    fit            = fit,
    K              = nrow(data[["outcome"]])
  )

  ### compute marginal likelihood
  marglik <- BayesTools::JAGS_bridgesampling(
    fit                = fit,
    log_posterior      = .log_posterior,
    data               = fit_data,
    prior_list         = fit_priors,
    formula_list       = .optional_jags_list(fit_formula_args[["formula_list"]]),
    formula_data_list  = .optional_jags_list(fit_formula_args[["formula_data_list"]]),
    formula_prior_list = .optional_jags_list(fit_formula_args[["formula_prior_list"]]),
    formula_scale_list = .optional_jags_list(fit_formula_args[["formula_scale_list"]]),
    formula_random_prior_list           = .optional_jags_list(fit_formula_args[["formula_random_prior_list"]]),
    formula_random_effects_compile_list = .optional_jags_list(fit_formula_args[["formula_random_effects_compile_list"]]),
    add_parameters                      = .optional_jags_character(bridge_sd_source_spec[["parameters"]]),
    add_bounds                          = bridge_sd_source_spec[["bounds"]],
    bridge_context                      = .marglik_needs_bridge_context(data),
    # additional arguments passed to .log_posterior via ...
    is_mods                  = .is_data_mods(data),
    is_scale                 = .is_data_scale(data),
    is_random                = .is_data_random(data),
    is_multilevel            = .is_data_multilevel(data),
    is_weights               = .is_data_weights(data),
    is_known_v               = .is_data_known_v(data),
    known_v_parameterization = .data_known_v_parameterization(data),
    model_data               = data,
    is_PET                   = .is_priors_PET(priors),
    is_PEESE                 = .is_priors_PEESE(priors),
    is_weightfunction        = .is_priors_weightfunction(priors),
    fixed_tau                = .fixed_tau_prior_value(priors),
    fixed_zero_tau           = .has_fixed_zero_tau_prior(priors),
    effect_direction         = .data_effect_direction(data),
    outcome_type             = .data_outcome_type(data)
  )

  return(marglik)
}


.marglik_needs_bridge_context <- function(data) {

  .is_data_known_v(data) &&
    .is_data_random(data) &&
    .data_has_marginalized_random_effects(data)
}


.check_marglik_available <- function(object, caller) {

  if (inherits(object, "RoBMA")) {
    stop(
      "Marginal likelihood is not available for product-space ",
      "model-averaging objects (BMA.norm, BMA.glmm, RoBMA).",
      call. = FALSE
    )
  }
  if (.is_random(object) &&
      !(inherits(object, "brma.mv") && .is_data_known_v(object[["data"]]))) {
    .check_random_formula_postfit_deferred(object, caller)
  }

  # Public constructors reject p-hacking/composed selection kernels earlier.
  if (.marglik_has_composed_bias(object[["priors"]][["outcome"]][["bias"]])) {
    stop(
      "Marginal likelihood is not available for combined ",
      "prior_bias(selection, phacking) models yet.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.marglik_has_composed_bias <- function(prior) {

  if (is.null(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(any(vapply(as.list(prior), .marglik_has_composed_bias, logical(1))))
  }
  if (!.is_prior_bias_kernel(prior)) {
    return(FALSE)
  }

  return(.prior_has_phacking(prior))
}


.marglik_add_selection_bridge_data <- function(fit_data, priors, effect_direction) {

  if (!.is_priors_weightfunction(priors) || is.null(fit_data[["yi"]])) {
    return(fit_data)
  }

  selection_spec <- .selection_spec(
    priors           = priors,
    yi               = fit_data[["yi"]],
    sei              = fit_data[["sei"]],
    effect_direction = effect_direction,
    signed_data      = TRUE
  )

  if (is.null(selection_spec)) {
    return(fit_data)
  }

  fit_data[["sel_kernel_mode"]]           <- selection_spec[["kernel_mode"]]
  fit_data[["sel_segment_bounds"]]        <- .selection_jags_bounds(selection_spec[["segments"]][["bounds"]])
  fit_data[["sel_segment_step_bin"]]      <- selection_spec[["segments"]][["step_bin"]]
  fit_data[["sel_segment_phack_region"]]  <- selection_spec[["segments"]][["phack_region"]]
  fit_data[["sel_phack_q"]]               <- selection_spec[["phack_q"]]
  fit_data[["sel_phack_z_source"]]        <- selection_spec[["phack_z_source"]]
  fit_data[["sel_phack_z_dest"]]          <- selection_spec[["phack_z_dest"]]

  return(fit_data)
}


.marglik_bridge_sd_source_spec <- function(add_parameters, fit, K) {

  posterior_names <- colnames(suppressWarnings(as.matrix(coda::as.mcmc(fit))))
  parameters      <- character()
  lb              <- numeric()
  ub              <- numeric()

  add_parameter <- function(parameter, lower, upper) {

    if (parameter %in% parameters) {
      return(invisible(NULL))
    }

    parameters <<- c(parameters, parameter)
    lb <<- c(lb, stats::setNames(lower, parameter))
    ub <<- c(ub, stats::setNames(upper, parameter))

    invisible(NULL)
  }

  for (parameter in add_parameters) {
    if (parameter %in% posterior_names) {
      add_parameter(parameter, 0, Inf)
      next
    }

    indexed_parameter <- paste0(parameter, "[", seq_len(K), "]")
    if (all(indexed_parameter %in% posterior_names)) {
      for (indexed in indexed_parameter) {
        add_parameter(indexed, 0, Inf)
      }
      next
    }

    add_parameter(parameter, 0, Inf)
  }

  if (length(parameters) == 0L) {
    return(list(
      parameters = character(),
      bounds     = NULL
    ))
  }

  return(list(
    parameters = parameters,
    bounds     = list(
      lb = lb[parameters],
      ub = ub[parameters]
    )
  ))
}


# ---------------------------------------------------------------------------- #
# .log_posterior
# ---------------------------------------------------------------------------- #
#
# Compute the log-likelihood of the data given a single posterior sample.
# This function is called by BayesTools::JAGS_bridgesampling for each
# posterior sample during bridge sampling.
#
# The function converts single-sample parameters to 1-row matrices to reuse
# the existing evaluate and pdf helper functions from evaluate.R and pdf.R.
# (in contrast to .pdf.brma, this function does not add likelihood contributions from 
#  estimate-level effects / baserate/ lograte for GLMMs since those are automatically 
#  handled via JAGS_bridgesampling function -- i.e., all parameters generated directly 
#  by BayesTools are automatically evaluated in the likelihood calculation, the only 
#  part missing is the final p(data | parameters) term) 
#
# @param parameters       named list of parameter values (preprocessed by BayesTools)
#                         - mu: scalar (no mods) or vector of length K (with mods)
#                         - tau: scalar (no scale regression)
#                         - log_tau: vector of length K (scale regression)
#                         - rho: scalar (if multilevel)
#                         - gamma: vector of length n_clusters (if multilevel)
#                         - PET: scalar (if PET model)
#                         - PEESE: scalar (if PEESE model)
#                         - omega: vector of weights (if weightfunction)
#                         - pi: vector of baserates (if binomial)
#                         - phi: vector of log-rates (if Poisson)
#                         - theta: vector of random effects (if GLMM)
# @param data             list containing fit_data (from .create_fit_data)
# @param is_mods          logical; whether model has moderators
# @param is_scale         logical; whether model has scale regression
# @param is_multilevel    logical; whether model is multilevel
# @param is_weights       logical; whether model uses weights (currently unused)
# @param is_PET           logical; whether model includes PET
# @param is_PEESE         logical; whether model includes PEESE
# @param is_weightfunction logical; whether model includes selection weights
# @param effect_direction character; "positive" or "negative"
# @param outcome_type     character; "norm", "bin", or "pois"
#
# @return scalar; sum of log-likelihoods across all observations
#
# ---------------------------------------------------------------------------- #
.log_posterior <- function(
    parameters, data,
    is_mods, is_scale, is_multilevel, is_weights,
    is_known_v, is_PET, is_PEESE, is_weightfunction, effect_direction,
    outcome_type, is_random = FALSE, known_v_parameterization = "latent",
    model_data = NULL, bridge_context = NULL, fixed_tau = NULL,
    fixed_zero_tau = FALSE) {

  ### extract number of observations
  K <- data[["K"]]

  ### compute mu samples as 1 x K matrix
  # BayesTools evaluates the formula and returns:
  # - scalar mu (no mods) -> replicate to K columns

  # - vector mu of length K (with mods) -> use directly
  mu_samples <- .marglik_get_mu_samples(
    parameters       = parameters,
    is_mods          = is_mods,
    is_PET           = is_PET,
    is_PEESE         = is_PEESE,
    effect_direction = effect_direction,
    sei              = if (outcome_type == "norm") data[["sei"]] else NULL,
    K                = K
  )

  ### compute tau samples as 1 x K matrix
  # BayesTools returns:
  # - scalar tau (no scale regression) -> replicate to K columns
  # - vector log_tau of length K (scale regression) -> exponentiate, use directly
  tau_result <- if (is_random) {
    .marglik_get_zero_tau_samples(K)
  } else {
    .marglik_get_tau_samples(
      parameters      = parameters,
      is_scale        = is_scale,
      is_multilevel   = is_multilevel,
      K               = K,
      fixed_tau       = fixed_tau,
      fixed_zero_tau  = fixed_zero_tau
    )
  }
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  ### add cluster-level (gamma) contribution for multilevel models
  if (is_multilevel) {
    cluster_contribution <- .marglik_get_cluster_effects(
      parameters  = parameters,
      tau_between = tau_between_samples,
      cluster     = data[["cluster"]]
    )
    mu_samples <- mu_samples + cluster_contribution
  }

  if (is_known_v) {
    sampling_dependency <- .marglik_get_sampling_dependency(
      parameters       = parameters,
      data             = data,
      effect_direction = effect_direction,
      K                = K
    )
    mu_samples <- mu_samples + sampling_dependency
  }

  ### dispatch to appropriate log-likelihood computation based on outcome type
  if (outcome_type == "norm") {

    # for selection models, PET, and PEESE, the outcome is in "positive" space
    # (yi flipped for negative effect direction in .create_fit_data)
    # so we need to flip mu_samples to match
    if (effect_direction == "negative") {
      mu_samples <- -mu_samples
    }

    if (is_known_v) {

      log_lik <- .marglik_known_v_norm_log_lik(
        parameters               = parameters,
        data                     = data,
        model_data               = model_data,
        bridge_context           = bridge_context,
        mu_samples               = mu_samples,
        tau_within_samples       = tau_within_samples,
        is_random                = is_random,
        is_weightfunction        = is_weightfunction,
        known_v_parameterization = known_v_parameterization,
        effect_direction         = effect_direction,
        K                        = K
      )

    } else if (is_weightfunction) {

      selection_context <- .marglik_selection_context(parameters, data)
      log_lik <- .outcome_pdf.selnorm(
        yi                = data[["yi"]],
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = data[["sei"]],
        selection_sei     = data[["sei"]],
        selection_context = selection_context
      )

    } else {

      # standard normal likelihood
      log_lik <- .outcome_pdf.norm(
        yi         = data[["yi"]],
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = data[["sei"]]
      )

    }

  } else if (outcome_type == "bin") {

    # add GLMM random effects (theta * tau_within) for binomial models
    theta_contribution <- .marglik_get_theta_samples(
      parameters = parameters,
      tau_within = tau_within_samples,
      K          = K
    )
    mu_samples <- mu_samples + theta_contribution

    # get baserate samples
    logit_baserate <- .marglik_get_baserate_samples(
      parameters = parameters,
      K          = K
    )

    log_lik <- .outcome_pdf.binom_conditional(
      ai             = data[["ai"]],
      ci             = data[["ci"]],
      n1i            = data[["n1i"]],
      n2i            = data[["n2i"]],
      mu_samples     = mu_samples,
      logit_baserate = logit_baserate
    )

  } else if (outcome_type == "pois") {

    # add GLMM random effects (theta * tau_within) for Poisson models
    theta_contribution <- .marglik_get_theta_samples(
      parameters = parameters,
      tau_within = tau_within_samples,
      K          = K
    )
    mu_samples <- mu_samples + theta_contribution

    # get log-rate samples
    log_phi <- .marglik_get_lograte_samples(
      parameters = parameters,
      K          = K
    )

    log_lik <- .outcome_pdf.pois_conditional(
      x1i        = data[["x1i"]],
      x2i        = data[["x2i"]],
      t1i        = data[["t1i"]],
      t2i        = data[["t2i"]],
      mu_samples = mu_samples,
      log_phi    = log_phi
    )

  }

  if (is_weights) {
    log_lik <- .apply_log_lik_weights(log_lik, data[["weight"]])
  }

  ### return sum of log-likelihoods (scalar)
  return(sum(log_lik))
}


# ---------------------------------------------------------------------------- #
# Helper functions for .log_posterior
# ---------------------------------------------------------------------------- #
#
# These functions convert single-sample parameters from BayesTools to 1-row
# matrices compatible with the evaluate and pdf helper functions.
# ---------------------------------------------------------------------------- #


#' @keywords internal
.marglik_get_mu_samples <- function(parameters, is_mods, is_PET, is_PEESE,
                                    effect_direction, sei, K) {

  # BayesTools returns:
  # - scalar mu (no mods)
  #   => replicate to 1 x K matrix
  # - vector mu of length K (with mods, formula already evaluated)
  #   => vector of length K -> 1 x K matrix
  mu_samples <- matrix(parameters[["mu"]], nrow = 1, ncol = K)

  # add PET adjustment (PET * sei)
  # direction multiplier: +1 for positive, -1 for negative
  if (is_PET) {
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * parameters[["PET"]] * sei
  }

  # add PEESE adjustment (PEESE * sei^2)
  if (is_PEESE) {
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * parameters[["PEESE"]] * sei^2
  }

  return(mu_samples)
}


#' @keywords internal
.marglik_get_tau_samples <- function(parameters, is_scale, is_multilevel, K,
                                     fixed_tau = NULL,
                                     fixed_zero_tau = FALSE) {

  # BayesTools returns:
  # - scalar tau (no scale regression)
  # - vector log_tau of length K (scale regression, need to exponentiate)
  if (is_scale) {
    # log_tau vector -> exponentiate and convert to 1 x K matrix
    tau_samples <- matrix(exp(parameters[["log_tau"]]), nrow = 1, ncol = K)
  } else {
    missing_tau <- .resolve_missing_tau_value(
      if (!is.null(fixed_tau)) fixed_tau else fixed_zero_tau
    )
    if (!is.null(missing_tau) && is.null(parameters[["tau"]])) {
      tau_samples <- matrix(missing_tau, nrow = 1, ncol = K)
    } else {
      # scalar tau -> replicate to 1 x K matrix
      tau_samples <- matrix(parameters[["tau"]], nrow = 1, ncol = K)
    }
  }

  # split tau into within/between components for multilevel models
  if (is_multilevel) {
    # extract rho (proportion of variance at cluster-level)
    rho <- parameters[["rho"]]

    # clamp rho to [0, 1] to handle numerical precision issues
    rho <- min(max(rho, 0), 1)

    # tau_within  = tau * sqrt(1 - rho)  (estimate-level heterogeneity)
    # tau_between = tau * sqrt(rho)      (cluster-level heterogeneity)
    tau_within_samples  <- tau_samples * sqrt(1 - rho)
    tau_between_samples <- tau_samples * sqrt(rho)

  } else {

    # non-multilevel: all heterogeneity is at estimate-level
    tau_within_samples  <- tau_samples
    tau_between_samples <- matrix(0, nrow = 1, ncol = K)

  }

  return(list(
    tau_total   = tau_samples,
    tau_within  = tau_within_samples,
    tau_between = tau_between_samples
  ))
}


.marglik_get_zero_tau_samples <- function(K) {

  tau_samples <- matrix(0, nrow = 1L, ncol = K)

  return(list(
    tau_total   = tau_samples,
    tau_within  = tau_samples,
    tau_between = tau_samples
  ))
}


#' @keywords internal
.marglik_get_cluster_effects <- function(parameters, tau_between, cluster) {

  K <- ncol(tau_between)

  # extract gamma samples (vector of length n_clusters)
  # BayesTools returns gamma as a vector
  gamma <- parameters[["gamma"]]

  # gamma[cluster] maps each observation to its cluster-level effect
  # multiply by tau_between to get the contribution
  cluster_contribution <- matrix(gamma[cluster] * tau_between, nrow = 1, ncol = K)

  return(cluster_contribution)
}


.marglik_selection_context <- function(parameters, data) {

  n_bins <- length(data[["sel_z_lower"]])
  phack_z_source <- if (!is.null(data[["phack_z_source"]])) {
    data[["phack_z_source"]]
  } else if (!is.null(data[["sel_phack_z_source"]])) {
    data[["sel_phack_z_source"]]
  } else {
    c(0, 0)
  }
  phack_z_dest <- if (!is.null(data[["phack_z_dest"]])) {
    data[["phack_z_dest"]]
  } else if (!is.null(data[["sel_phack_z_dest"]])) {
    data[["sel_phack_z_dest"]]
  } else {
    c(0, 0)
  }
  omega <- if (!is.null(parameters[["omega"]])) {
    matrix(parameters[["omega"]], nrow = 1)
  } else if (!is.null(data[["sel_omega"]])) {
    matrix(data[["sel_omega"]], nrow = 1)
  } else {
    matrix(1, nrow = 1, ncol = n_bins)
  }

  alpha <- if (!is.null(parameters[["alpha"]])) parameters[["alpha"]] else {
    if (!is.null(data[["sel_phack_alpha"]])) data[["sel_phack_alpha"]] else 0
  }
  phack_kind <- if (!is.null(parameters[["phack_kind"]])) {
    as.integer(parameters[["phack_kind"]])
  } else {
    if (!is.null(data[["sel_phack_kind"]])) {
      as.integer(data[["sel_phack_kind"]])
    } else if (!is.null(data[["sel_phack_q"]]) &&
               data[["sel_kernel_mode"]] %in% c(SELKERNEL_PHACK_POWER, SELKERNEL_STEP_PHACK_POWER)) {
      as.integer(data[["sel_phack_q"]])
    } else {
      0L
    }
  }

  selection_context <- list(
    kernel_mode    = data[["sel_kernel_mode"]],
    z_lower        = data[["sel_z_lower"]],
    z_upper        = data[["sel_z_upper"]],
    obs_bin        = data[["sel_obs_bin"]],
    sign           = data[["sel_sign"]],
    n_bins         = n_bins,
    has_step       = data[["sel_kernel_mode"]] %in% c(SELKERNEL_STEP, SELKERNEL_STEP_PHACK_POWER),
    has_phack      = data[["sel_kernel_mode"]] %in% c(SELKERNEL_PHACK_POWER, SELKERNEL_STEP_PHACK_POWER),
    phack_q        = if (phack_kind > 0L) phack_kind else 1L,
    phack_z_source = phack_z_source,
    phack_z_dest   = phack_z_dest,
    segments       = list(
      bounds       = data[["sel_segment_bounds"]],
      step_bin     = data[["sel_segment_step_bin"]],
      phack_region = data[["sel_segment_phack_region"]]
    ),
    omega          = omega,
    alpha          = alpha,
    phack_kind     = phack_kind
  )

  return(.selection_reset_native_cache(selection_context))
}


.marglik_get_sampling_dependency <- function(parameters, data,
                                             effect_direction, K) {

  if (is.null(data[["sampling_rank"]]) || data[["sampling_rank"]] == 0L) {
    return(matrix(0, nrow = 1, ncol = K))
  }

  sampling_z <- parameters[["sampling_z"]]
  if (is.null(sampling_z)) {
    sampling_z_names <- paste0("sampling_z[", seq_len(data[["sampling_rank"]]), "]")
    if (!all(sampling_z_names %in% names(parameters))) {
      stop("Missing known-V latent sampling factors.", call. = FALSE)
    }
    sampling_z <- unlist(parameters[sampling_z_names], use.names = FALSE)
  }

  sampling_dependency <- matrix(
    as.vector(data[["sampling_B"]] %*% sampling_z),
    nrow = 1,
    ncol = K
  )

  if (effect_direction == "negative") {
    sampling_dependency <- -sampling_dependency
  }

  return(sampling_dependency)
}


.marglik_known_v_norm_log_lik <- function(parameters, data, model_data,
                                          bridge_context,
                                          mu_samples, tau_within_samples,
                                          is_random, is_weightfunction,
                                          known_v_parameterization,
                                          effect_direction, K) {

  extra_variance <- .marglik_known_v_extra_variance(
    parameters         = parameters,
    model_data         = model_data,
    bridge_context     = bridge_context,
    tau_within_samples = tau_within_samples,
    is_random          = is_random,
    K                  = K
  )
  tau_within <- sqrt(pmax(extra_variance, 0))

  if (known_v_parameterization == "latent") {
    if (is_weightfunction) {
      selection_context <- .marglik_selection_context(parameters, data)
      return(.outcome_pdf.selnorm(
        yi                = data[["yi"]],
        mu_samples        = mu_samples,
        tau_within        = tau_within,
        sei               = sqrt(data[["sampling_var"]]),
        selection_sei     = data[["sei"]],
        selection_context = selection_context
      ))
    }

    return(.outcome_pdf.norm(
      yi         = data[["yi"]],
      mu_samples = mu_samples,
      tau_within = tau_within,
      sei        = sqrt(data[["sampling_var"]])
    ))
  }

  if (is_weightfunction) {
    stop(
      "Selection-model bridge likelihoods for known-V data require ",
      "known_v_parameterization = 'latent'.",
      call. = FALSE
    )
  }
  if (is.null(model_data) || !.is_data_known_v(model_data)) {
    stop("Known-V bridge likelihood requires known-V model metadata.",
         call. = FALSE)
  }

  if (known_v_parameterization == "whitened") {
    return(.marglik_known_v_whitened_log_lik(
      data           = data,
      mu_samples     = mu_samples,
      extra_variance = extra_variance,
      K              = K
    ))
  }
  if (known_v_parameterization == "block_mvn") {
    return(.marglik_known_v_block_mvn_log_lik_sum(
      model_data       = model_data,
      mu_samples       = mu_samples,
      extra_variance   = extra_variance,
      effect_direction = effect_direction,
      K                = K
    ))
  }

  stop("Unknown known-V parameterization: ", known_v_parameterization,
       call. = FALSE)
}


.marglik_known_v_extra_variance <- function(parameters, model_data,
                                            bridge_context,
                                            tau_within_samples, is_random, K) {

  extra_variance <- if (is_random) {
    if (is.null(model_data)) {
      stop(
        "Random-formula known-V bridge likelihood requires model metadata.",
        call. = FALSE
      )
    }
    .evaluate_marginalized_random_variance(
      data              = model_data,
      posterior_samples = .marglik_bridge_posterior_samples(
        parameters     = parameters,
        bridge_context = bridge_context
      ),
      K                 = K
    )
  } else {
    tau_within_samples^2
  }
  extra_variance <- as.matrix(extra_variance)

  if (nrow(extra_variance) != 1L || ncol(extra_variance) != K) {
    stop(
      "Known-V bridge diagonal variance contributions have inconsistent dimensions.",
      call. = FALSE
    )
  }
  if (any(!is.finite(extra_variance)) || any(extra_variance < 0)) {
    stop(
      "Known-V bridge diagonal variance contributions must be non-negative.",
      call. = FALSE
    )
  }

  return(extra_variance)
}


.marglik_known_v_whitened_log_lik <- function(data, mu_samples,
                                              extra_variance, K) {

  whitening_mu <- as.vector(data[["whitening_matrix"]] %*%
    as.numeric(mu_samples[1L, ]))
  variance     <- data[["whitening_var"]] + as.numeric(extra_variance[1L, ])

  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Known-V whitened bridge variances must be positive.", call. = FALSE)
  }

  matrix(
    stats::dnorm(
      x    = data[["whitening_y"]],
      mean = whitening_mu,
      sd   = sqrt(variance),
      log  = TRUE
    ),
    nrow = 1L,
    ncol = K
  )
}


.marglik_known_v_block_mvn_log_lik_sum <- function(model_data, mu_samples,
                                                   extra_variance,
                                                   effect_direction, K) {

  known_V       <- .data_known_v_data(model_data)
  yi            <- model_data[["outcome"]][["yi"]]
  block_indices <- known_V[["block_indices"]]
  log_lik       <- 0

  if (effect_direction == "negative") {
    yi <- -yi
  }
  if (is.null(block_indices) || length(block_indices) == 0L) {
    block_indices <- as.list(seq_len(K))
  }

  for (idx in block_indices) {
    covariance <- known_V[["V"]][idx, idx, drop = FALSE] +
      diag(as.numeric(extra_variance[1L, idx]), nrow = length(idx))
    log_lik <- log_lik + .marglik_mvn_log_density(
      y          = yi[idx],
      mean       = as.numeric(mu_samples[1L, idx]),
      covariance = covariance
    )
  }

  return(log_lik)
}


.marglik_mvn_log_density <- function(y, mean, covariance) {

  chol_covariance <- tryCatch(
    chol(covariance),
    error = function(e) NULL
  )
  if (is.null(chol_covariance)) {
    stop("Known-V bridge covariance is not positive definite.", call. = FALSE)
  }

  residual <- y - mean
  z        <- backsolve(chol_covariance, residual, transpose = TRUE)
  size     <- length(y)

  -0.5 * (
    size * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      sum(z^2)
  )
}


.parameters_as_sample_matrix <- function(parameters) {

  values <- numeric()
  for (name in names(parameters)) {
    value <- parameters[[name]]
    if (length(value) == 0L) {
      next
    }
    value <- as.numeric(value)
    names(value) <- if (length(value) == 1L) {
      name
    } else {
      paste0(name, "[", seq_along(value), "]")
    }
    values <- c(values, value)
  }

  posterior_row <- attr(parameters, "posterior_samples", exact = TRUE)
  if (!is.null(posterior_row)) {
    posterior_row <- as.matrix(posterior_row)
    if (nrow(posterior_row) != 1L || is.null(colnames(posterior_row))) {
      stop("Bridge posterior row metadata have invalid dimensions.",
           call. = FALSE)
    }
    extra <- setdiff(colnames(posterior_row), names(values))
    if (length(extra) > 0L) {
      extra_values <- as.numeric(posterior_row[1L, extra, drop = TRUE])
      names(extra_values) <- extra
      values <- c(values, extra_values)
    }
  }

  out <- matrix(values, nrow = 1L)
  colnames(out) <- names(values)
  return(out)
}


.marglik_bridge_posterior_samples <- function(parameters, bridge_context = NULL) {

  posterior_samples <- .parameters_as_sample_matrix(parameters)

  if (!inherits(bridge_context, "BayesTools_bridge_context")) {
    return(posterior_samples)
  }

  nodes <- bridge_context[["nodes"]]
  if (length(nodes) == 0L) {
    return(posterior_samples)
  }
  if (!is.numeric(nodes) || is.null(names(nodes))) {
    stop("Bridge context nodes must be a named numeric vector.",
         call. = FALSE)
  }

  nodes <- nodes[!is.na(names(nodes)) & nzchar(names(nodes))]
  extra <- setdiff(names(nodes), colnames(posterior_samples))
  if (length(extra) == 0L) {
    return(posterior_samples)
  }

  extra_samples <- matrix(
    as.numeric(nodes[extra]),
    nrow     = 1L,
    dimnames = list(NULL, extra)
  )
  posterior_samples <- cbind(posterior_samples, extra_samples)

  return(posterior_samples)
}


#' @keywords internal
.marglik_get_theta_samples <- function(parameters, tau_within, K) {

  # theta is a vector of length K (estimate-level random effects for GLMM)
  theta <- parameters[["theta"]]

  # theta contribution = theta[k] * tau_within[k]
  theta_contribution <- matrix(theta * tau_within, nrow = 1, ncol = K)

  return(theta_contribution)
}


#' @keywords internal
.marglik_get_baserate_samples <- function(parameters, K) {

  # pi is a vector of length K (baserates for binomial models)
  pi <- parameters[["pi"]]

  # return logit(pi) as 1 x K matrix
  logit_baserate <- matrix(.logit(pi), nrow = 1, ncol = K)

  return(logit_baserate)
}


#' @keywords internal
.marglik_get_lograte_samples <- function(parameters, K) {

  # phi is a vector of length K (log-rates for Poisson models)
  phi <- parameters[["phi"]]

  # return log(phi) as 1 x K matrix
  # Note: phi is already on log scale in JAGS model
  log_phi <- matrix(phi, nrow = 1, ncol = K)

  return(log_phi)
}
