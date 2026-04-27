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
#' @return The brma object with the marginal likelihood result stored in
#' \code{object[["marglik"]]}.
#'
#' @seealso \code{\link{bridge_sampler.brma}}, \code{\link{logml.brma}},
#' \code{\link{bf.brma}}, \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' fit <- brma(y = d, se = se, data = Bem2011)
#'
#' # Compute and store marginal likelihood
#' fit <- add_marglik(fit)
#'
#' # Extract bridge sampling object
#' bridge <- bridge_sampler(fit)
#' print(bridge)
#'
#' # Get log marginal likelihood
#' logml(fit)
#' }
#'
#' @export
add_marglik.brma <- function(object, ...) {
  marglik <- .marglik(object)
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

  ### create arguments to be passed to BayesTools::JAGS_bridgesampling
  fit_formula_list        <- list()
  fit_formula_data_list   <- list()
  fit_formula_prior_list  <- list()
  fit_formula_scale_list  <- list()

  ### create model base
  fit_priors <- .create_fit_priors(data = data, priors = priors)
  fit_data   <- .create_fit_data(data = data, priors = priors)

  ### add effect regressions
  if (.is_data_mods(data)) {
    fit_formula_list[["mu"]]       <- .create_fit_formula_list(data = data, parameter = "mods")
    fit_formula_data_list[["mu"]]  <- .create_fit_formula_data_list(data = data, parameter = "mods")
    fit_formula_prior_list[["mu"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "mods")
    fit_formula_scale_list[["mu"]] <- .data_standardize_continuous_predictors(data)
  }

  ### add heterogeneity regressions
  if (.is_data_scale(data)) {
    fit_formula_list[["log_tau"]]       <- .create_fit_formula_list(data = data, parameter = "scale")
    fit_formula_data_list[["log_tau"]]  <- .create_fit_formula_data_list(data = data, parameter = "scale")
    fit_formula_prior_list[["log_tau"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "scale")
    fit_formula_scale_list[["log_tau"]] <- .data_standardize_continuous_predictors(data)
  }

  ### compute marginal likelihood
  marglik <- BayesTools::JAGS_bridgesampling(
    fit                = fit,
    log_posterior      = .log_posterior,
    data               = fit_data,
    prior_list         = fit_priors,
    formula_list       = if (length(fit_formula_list)       > 0) fit_formula_list       else NULL,
    formula_data_list  = if (length(fit_formula_data_list)  > 0) fit_formula_data_list  else NULL,
    formula_prior_list = if (length(fit_formula_prior_list) > 0) fit_formula_prior_list else NULL,
    formula_scale      = if (length(fit_formula_scale_list) > 0) fit_formula_scale_list else NULL,
    # additional arguments passed to .log_posterior via ...
    is_mods            = .is_data_mods(data),
    is_scale           = .is_data_scale(data),
    is_multilevel      = .is_data_multilevel(data),
    is_weights         = .is_data_weights(data),
    is_PET             = .is_priors_PET(priors),
    is_PEESE           = .is_priors_PEESE(priors),
    is_weightfunction  = .is_priors_weightfunction(priors),
    effect_direction   = .data_effect_direction(data),
    outcome_type       = .data_outcome_type(data)
  )

  return(marglik)
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
    is_PET, is_PEESE, is_weightfunction, effect_direction,
    outcome_type) {

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
  tau_result <- .marglik_get_tau_samples(
    parameters    = parameters,
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  ### add cluster-level (gamma) contribution for multilevel models
  if (is_multilevel) {
    cluster_contribution <- .marglik_get_cluster_effects(
      parameters       = parameters,
      tau_between      = tau_between_samples,
      cluster          = data[["cluster"]],
      effect_direction = effect_direction
    )
    mu_samples <- mu_samples + cluster_contribution
  }

  ### dispatch to appropriate log-likelihood computation based on outcome type
  if (outcome_type == "norm") {

    # for selection models, PET, and PEESE, the outcome is in "positive" space
    # (yi flipped for negative effect direction in .create_fit_data)
    # so we need to flip mu_samples to match
    if (effect_direction == "negative") {
      mu_samples <- -mu_samples
    }

    if (is_weightfunction) {

      # weighted normal likelihood (selection models)
      log_lik <- .outcome_pdf.wnorm(
        yi         = data[["yi"]],
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = data[["sei"]],
        omega      = .marglik_get_omega_samples(parameters),
        crit_yi    = data[["crit_yi"]],
        bias_indicator      = .marglik_get_bias_indicator(parameters, data),
        crit_yi_mapping     = data[["crit_yi_mapping"]],
        crit_yi_mapping_max = data[["crit_yi_mapping_max"]]
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
.marglik_get_tau_samples <- function(parameters, is_scale, is_multilevel, K) {

  # BayesTools returns:
  # - scalar tau (no scale regression)
  # - vector log_tau of length K (scale regression, need to exponentiate)
  if (is_scale) {
    # log_tau vector -> exponentiate and convert to 1 x K matrix
    tau_samples <- matrix(exp(parameters[["log_tau"]]), nrow = 1, ncol = K)
  } else {
    # scalar tau -> replicate to 1 x K matrix
    tau_samples <- matrix(parameters[["tau"]], nrow = 1, ncol = K)
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


#' @keywords internal
.marglik_get_cluster_effects <- function(parameters, tau_between, cluster,
                                         effect_direction) {

  K <- ncol(tau_between)

  # extract gamma samples (vector of length n_clusters)
  # BayesTools returns gamma as a vector
  gamma <- parameters[["gamma"]]

  # gamma[cluster] maps each observation to its cluster-level effect
  # multiply by tau_between to get the contribution
  cluster_contribution <- matrix(gamma[cluster] * tau_between, nrow = 1, ncol = K)

  return(cluster_contribution)
}


#' @keywords internal
.marglik_get_omega_samples <- function(parameters) {

  # omega is a vector of weights from BayesTools
  # convert to 1 x W matrix for compatibility with .dwnorm_fast.ss.matrix
  omega <- parameters[["omega"]]
  return(matrix(omega, nrow = 1))
}

.marglik_get_bias_indicator <- function(parameters, data) {
  if (!is.null(parameters[["bias_indicator"]])) {
    return(as.integer(parameters[["bias_indicator"]]))
  }
  if (!is.null(data[["bias_indicator"]])) {
    return(as.integer(data[["bias_indicator"]]))
  }
  return(1L)
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
