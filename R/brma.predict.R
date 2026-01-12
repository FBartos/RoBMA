#' @title Predict From brma Object
#'
#' @description \code{predict.brma} predicts values
#'
#' @inheritParams summary.brma
#' @param newdata a data.frame (if prediction for a meta-regression is performed) or
#' a list named list with the effect size measure and variability metrics (if prediction
#' for a meta-analysis is performed) for new studies. Note that the input has to corresponds
#' to the format and naming that was used to estimate the original fit. Defaults to
#' \code{NULL} which corresponds to prediction for the observed data.
#' @param type type of prediction to be performed. Defaults to \code{"response"} which
#' produces predictions for the observed effect size estimates. Alternatives are
#' \code{"terms"} which produces the mean effect size estimate at the given predictors
#' levels (not accounting for the random-effects) and \code{"effect"} which predicts the
#' distribution of the true study effects at the given predictors levels
#' (i.e., incorporating heterogeneity into \code{"terms"}).
#' @param bias_adjusted whether sampling of new values should adjust for publication bias
#' (note that selection models do not affect the mean paramater when \code{"terms"}
#' (equal mean parameter under normal vs. weighted likelihood equals different expectation).
#'
#' @details
#' Note that in contrast to \link[metafor]{predict}, the \code{type = "response"} produces
#' predictions for the new effect size estimates (instead of the true study effects).
#' To obtain results corresponding to the metafor's predict function, use the
#' \code{type = "terms"} to obtain the mean effect size estimate in its credible interval
#' and \code{type = "effect"} to obtain the distribution of the true study effects (i.e.,
#' prediction interval).
#'
#'
#' @examples \dontrun{
#' }
#'
#' @return \code{pooled_effect} returns a list of tables of class 'BayesTools_table'.
#' @seealso [true_effects()], [residuals.RoBMA()]
#' @export
predict.brma <- function(object, newdata,
                         type = "terms",
                         probs = c(.025, .975),
                         bias_adjusted = TRUE,
                         as_samples = FALSE,
                         ...){

  # some options checked inside BayesTools table directly
  BayesTools::check_char(type, "type", allow_values = c("response", "terms", "terms.scale", "effect"))
  BayesTools::check_bool(as_samples, "as_samples")
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")

  ### types of predictions
  # terms:       fixed effects terms for the overall effect (mu) / incorporating mods if present
  # terms.scale: fixed effects terms for the overall heterogeneity (tau) / incorporating scale if present
  # effect:      incorporating between-study heterogeneity into terms to obtain the true study effects
  #              (via empirical Bayes for random effects necessary, in case of new_data, new random effect is sampled)
  # response:    incorporating between-study heterogeneity and sampling variability
  #              (via marginalized random-effects)

  ### dispatch between prediction on the current data vs. new data
  if (missing(newdata) || is.null(newdata)) {

    # an existing data are used
    same_data <- TRUE
    new_data  <- object[["data"]]

  } else {

    # prepare newdata using the same settings as the original fit
    same_data <- FALSE
    new_data  <- .prepare_newdata(
      object                       = object,
      newdata                      = newdata,
      type                         = type,
      incorporate_publication_bias = incorporate_publication_bias
    )

  }

  ### extract priors and structural information about the model
  priors            <- object[["priors"]]
  is_mods           <- .is_mods(object)
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)

  ### extract outcome data and fit data for convenience
  outcome_data      <- new_data[["outcome"]]
  fit_data          <- .create_fit_data(data = new_data, priors = priors)

  # outcome dimensions
  K  <- nrow(outcome_data)

  ### obtain tau samples using helper function
  # returns list(tau_within, tau_between) - all S x K matrices
  # see .evaluate.brma.tau() in brma.evaluate.R for details
  tau_result          <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = new_data[["scale"]],
    scale_formula = if (is_scale) attr(object[["data"]][["scale"]], "formula") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]
  S                   <- nrow(tau_within_samples)  # number of posterior samples

  ### return only tau samples if type = "terms.scale" is selected
  if (type == "terms.scale") {
    # reconstruct total tau from components: tau = sqrt(tau_within^2 + tau_between^2)
    tau_samples <- sqrt(tau_within_samples^2 + tau_between_samples^2)

    # rename samples
    colnames(tau_samples) <- paste0("tau[", seq_len(K), "]")

    if (as_samples) {
      return(tau_samples)
    } else {
      outcome_table <- BayesTools::ensemble_estimates_table(
        samples    = asplit(tau_samples, 2),
        parameters = colnames(tau_samples),
        probs      = probs,
        title      = "Scale Term Posterior Prediction:"
      )
      out <- list(
        summary = outcome_table,
        data    = new_data
      )
      class(out) <- "brma.predict"
      return(out)
    }
  }

  ### get the base mu samples using helper function
  # returns S x K matrix of location samples
  # see .evaluate.brma.mu() in brma.evaluate.R for details
  mu_samples <- .evaluate.brma.mu(
    fit                           = object[["fit"]],
    outcome_data                  = outcome_data,
    mods_data                     = new_data[["mods"]],
    mods_formula                  = if (is_mods) attr(object[["data"]][["mods"]], "formula") else NULL,
    mods_priors                   = priors[["mods"]],
    is_mods                       = is_mods,
    is_PET                        = is_PET,
    is_PEESE                      = is_PEESE,
    effect_direction              = effect_direction,
    incorporate_publication_bias  = bias_adjusted,
    K                             = K
  )

  ### return only mu samples if type = "terms" is selected
  # terms incorporate fixed effects only (i.e., random effects are not incorporated)
  if (type == "terms") {
    # rename samples
    colnames(mu_samples) <- paste0("mu[", seq_len(K), "]")

    if (as_samples) {
      return(mu_samples)
    } else {
      outcome_table <- BayesTools::ensemble_estimates_table(
        samples    = asplit(mu_samples, 2),
        parameters = colnames(mu_samples),
        probs      = probs,
        title      = "Location Term Posterior Prediction:"
      )
      out <- list(
        summary = outcome_table,
        data    = new_data
      )
      class(out) <- "brma.predict"
      return(out)
    }
  }

  ### include 3-level (study-level) random-effects using helper function
  # returns contribution matrix gamma[study_id] * tau_between for multilevel models
  # see .evaluate.brma.study_effects() in brma.evaluate.R for details
  if (is_multilevel) {
    study_contribution <- .evaluate.brma.study_effects(
      fit              = object[["fit"]],
      tau_between      = tau_between_samples,
      study_ids        = fit_data[["study_ids"]],
      same_data        = same_data,
      effect_direction = effect_direction
    )
    mu_samples <- mu_samples + study_contribution
  }

  ### create true-effect prediction using helper function
  # computes empirical Bayes (BLUP) estimates of true study effects
  # see .evaluate.brma.true_effects() in brma.evaluate.R for details
  if (type == "effect") {

    true_effects_samples <- .evaluate.brma.true_effects(
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      yi         = outcome_data[["yi"]],
      sei        = outcome_data[["sei"]]
    )

    # rename samples
    colnames(true_effects_samples) <- paste0("theta[", seq_len(K), "]")

    if (as_samples) {
      return(true_effects_samples)
    } else {
      outcome_table <- BayesTools::ensemble_estimates_table(
        samples    = asplit(true_effects_samples, 2),
        parameters = colnames(true_effects_samples),
        probs      = probs,
        title      = "True Effect Posterior Prediction:"
      )
      out <- list(
        summary = outcome_table,
        data    = new_data
      )
      class(out) <- "brma.predict"
      return(out)
    }

  }

  ### create observed effects prediction
  if (type == "response") {

    # different model types have different output structures
    if (is.element(outcome_type, c("bin", "pois"))) {

      # the estimate-level effects are not marginalized for GLMMs
      # include the sampled random effects or marginalize over them for new data
      # see .evaluate.brma.theta.glmm() in brma.evaluate.R for details
      theta_contribution <- .evaluate.brma.theta.glmm(
        fit        = object[["fit"]],
        tau_within = tau_within_samples,
        same_data  = same_data,
        K          = K
      )
      mu_samples <- mu_samples + theta_contribution

      if (outcome_type == "bin") {

        ### incorporate with base-rate using helper function
        # .evaluate.brma.baserate() returns S x K matrix of logit(pi) samples
        logit_baserate <- .evaluate.brma.baserate(fit = object[["fit"]], K = K)

        ### sample outcome using RNG helper
        outcome_samples <- .outcome_rng.binom(
          mu_samples     = mu_samples,
          logit_baserate = logit_baserate,
          n1i            = outcome_data[["n1i"]],
          n2i            = outcome_data[["n2i"]]
        )

      } else if (outcome_type == "pois") {

        ### incorporate with log-rate using helper function
        # .evaluate.brma.lograte() returns S x K matrix of log(phi) samples
        log_phi <- .evaluate.brma.lograte(fit = object[["fit"]], K = K)

        ### sample outcome using RNG helper
        outcome_samples <- .outcome_rng.pois(
          mu_samples = mu_samples,
          log_phi    = log_phi,
          t1i        = outcome_data[["t1i"]],
          t2i        = outcome_data[["t2i"]]
        )

      }


    } else if (outcome_type == "norm") {

      # normal outcome models dispatch between:
      # - bias_adjusted = TRUE: sample from unweighted normal (as if no bias)
      # - bias_adjusted = FALSE with weightfunction: sample from weighted normal
      # - bias_adjusted = FALSE without weightfunction: sample from unweighted normal

      if (bias_adjusted || !is_weightfunction) {

        # sample from unweighted normal distribution
        # y[s,k] ~ N(mu[s,k], sqrt(tau_within[s,k]^2 + sei[k]^2))
        outcome_samples <- .outcome_rng.norm(
          mu_samples = mu_samples,
          tau_within = tau_within_samples,
          sei        = outcome_data[["sei"]]
        )

      } else {

        # sample from weighted normal distribution for selection models
        # extract omega samples for weight function
        posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
        omega_samples     <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

        outcome_samples <- .outcome_rng.wnorm(
          mu_samples = mu_samples,
          tau_within = tau_within_samples,
          sei        = outcome_data[["sei"]],
          omega      = omega_samples,
          crit_yi    = fit_data$crit_yi
        )

      }

      # rename samples
      colnames(outcome_samples) <- paste0("yi[", seq_len(K), "]")
    }

    if (as_samples) {
      return(outcome_samples)
    } else {
      outcome_table <- BayesTools::ensemble_estimates_table(
        samples    = asplit(outcome_samples, 2),
        parameters = colnames(outcome_samples),
        probs      = probs,
        title      = "Observations Posterior Prediction:"
      )
      out <- list(
        summary = outcome_table,
        data    = new_data
      )
      class(out) <- "brma.predict"
      return(out)
    }
  }
}



#' @export
summary.brma.predict <- function(object, ...) {
  print(object, ...)
}

#' @export
print.brma.predict <- function(x, ...) {

  cat("\n")
  print(x[["summary"]])
  cat("\n")

  return(invisible(x))
}


# ============================================================================ #
# Outcome RNG Functions
# ============================================================================ #
#
# These are prediction-specific RNG functions for sampling from the posterior
# predictive distribution of observed effect sizes. They are kept separate from
# the .evaluate.brma.* functions because they involve random number generation
# and are specific to the prediction task.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .outcome_rng.norm
# ---------------------------------------------------------------------------- #
#
# Sample observed effect sizes from normal posterior predictive distribution.
#
# For normal outcome models, the observed effect y_i is:
#   y_i ~ N(mu_i, tau_within^2 + se_i^2)
#
# This function samples from this distribution for each posterior sample.
#
# @param mu_samples       S x K matrix of location samples (with study effects if multilevel)
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
#
# @return S x K matrix of sampled observed effect sizes
#
# ---------------------------------------------------------------------------- #
.outcome_rng.norm <- function(mu_samples, tau_within, sei) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate sei across samples: matrix(sei, S, K, byrow = TRUE)
  # each row contains the same se values
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute total SD: sqrt(tau^2 + se^2)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # sample from N(mu, total_sd) for each cell
  # matrix(rnorm(S*K), S, K) * total_sd + mu is vectorized sampling
  # equivalent: for each s,k: y[s,k] ~ N(mu[s,k], total_sd[s,k])
  response_samples <- mu_samples + matrix(stats::rnorm(S * K), nrow = S, ncol = K) * total_sd

  return(response_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.wnorm
# ---------------------------------------------------------------------------- #
#
# Sample observed effect sizes from weighted normal distribution (selection models).
#
# For selection models, the observed effect y_i follows a weighted normal:
#   y_i ~ f(y) * omega(y) / Z
# where f(y) is normal density and omega(y) is the selection weight function.
#
# This uses the fast spike-and-slab weighted normal sampler.
#
# @param mu_samples       S x K matrix of location samples
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param omega            S x W matrix of omega (weight) samples
# @param crit_yi          W x K matrix of critical values for each observation
#
# @return S x K matrix of sampled observed effect sizes from weighted distribution
#
# ---------------------------------------------------------------------------- #
.outcome_rng.wnorm <- function(mu_samples, tau_within, sei, omega, crit_yi) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # pre-compute total SD for each observation
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # sample from weighted normal for each observation
  # (loop is necessary due to .rwnorm_fast.ss interface)
  response_samples <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    response_samples[, k] <- .rwnorm_fast.ss(
      mean   = mu_samples[, k],
      sd     = total_sd[, k],
      omega  = omega,
      crit_x = crit_yi[, k]
    )
  }

  return(response_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.binom
# ---------------------------------------------------------------------------- #
#
# Sample observed counts from binomial posterior predictive distribution.
#
# For binomial outcome models (log-odds ratio), the observed counts are:
#   ai ~ Binom(n1i, p1) where logit(p1) = logit(pi) + 0.5 * mu (treatment)
#   ci ~ Binom(n2i, p2) where logit(p2) = logit(pi) - 0.5 * mu (control)
#
# The 0.5 factor arises from the log-OR parameterization where mu represents
# the log-odds ratio between treatment and control groups.
#
# @param mu_samples       S x K matrix of log-odds ratio samples
# @param logit_baserate   S x K matrix of logit(pi) base-rate samples
# @param n1i              numeric vector of length K; treatment group sizes
# @param n2i              numeric vector of length K; control group sizes
#
# @return S x (2*K) matrix with columns interleaved: ai[1], ci[1], ai[2], ci[2], ...
#
# ---------------------------------------------------------------------------- #
.outcome_rng.binom <- function(mu_samples, logit_baserate, n1i, n2i) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # compute logit probabilities for each group
  # group 1: logit(p1) = logit(pi) + 0.5 * mu (treatment/exposed)
  # group 2: logit(p2) = logit(pi) - 0.5 * mu (control/unexposed)
  logit_p1 <- logit_baserate + 0.5 * mu_samples
  logit_p2 <- logit_baserate - 0.5 * mu_samples

  # sample from binomial for each group
  outcome_samples_ai <- matrix(NA_integer_, nrow = S, ncol = K)
  outcome_samples_ci <- matrix(NA_integer_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    outcome_samples_ai[, k] <- stats::rbinom(n = S, size = n1i[k], prob = .inv_logit(logit_p1[, k]))
    outcome_samples_ci[, k] <- stats::rbinom(n = S, size = n2i[k], prob = .inv_logit(logit_p2[, k]))
  }

  # name columns
  colnames(outcome_samples_ai) <- paste0("ai[", seq_len(K), "]")
  colnames(outcome_samples_ci) <- paste0("ci[", seq_len(K), "]")

  # merge samples with interleaved columns: ai[1], ci[1], ai[2], ci[2], ...
  outcome_samples <- matrix(NA_integer_, nrow = S, ncol = 2 * K)
  outcome_samples[, seq(1, 2 * K, by = 2)] <- outcome_samples_ai
  outcome_samples[, seq(2, 2 * K, by = 2)] <- outcome_samples_ci
  colnames(outcome_samples) <- as.vector(rbind(
    colnames(outcome_samples_ai),
    colnames(outcome_samples_ci)
  ))

  return(outcome_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.pois
# ---------------------------------------------------------------------------- #
#
# Sample observed counts from Poisson posterior predictive distribution.
#
# For Poisson outcome models (log incidence rate ratio), the observed counts are:
#   x1i ~ Pois(t1i * exp(phi + 0.5 * mu))  (treatment/exposed)
#   x2i ~ Pois(t2i * exp(phi - 0.5 * mu))  (control/unexposed)
#
# where phi is the log baseline rate and t1i, t2i are exposure times.
#
# @param mu_samples       S x K matrix of log incidence rate ratio samples
# @param log_phi          S x K matrix of log baseline rate samples
# @param t1i              numeric vector of length K; treatment exposure times
# @param t2i              numeric vector of length K; control exposure times
#
# @return S x (2*K) matrix with columns interleaved: x1i[1], x2i[1], x1i[2], x2i[2], ...
#
# ---------------------------------------------------------------------------- #
.outcome_rng.pois <- function(mu_samples, log_phi, t1i, t2i) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate exposure times to S x K matrices for vectorized ops
  log_t1i <- matrix(log(t1i), nrow = S, ncol = K, byrow = TRUE)
  log_t2i <- matrix(log(t2i), nrow = S, ncol = K, byrow = TRUE)

  # compute log rates for each group
  # group 1: log(lambda1) = phi + 0.5 * mu + log(t1i)
  # group 2: log(lambda2) = phi - 0.5 * mu + log(t2i)
  log_lambda1 <- log_phi + 0.5 * mu_samples + log_t1i
  log_lambda2 <- log_phi - 0.5 * mu_samples + log_t2i

  # sample from Poisson for each group
  outcome_samples_x1i <- matrix(NA_integer_, nrow = S, ncol = K)
  outcome_samples_x2i <- matrix(NA_integer_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    outcome_samples_x1i[, k] <- stats::rpois(n = S, lambda = exp(log_lambda1[, k]))
    outcome_samples_x2i[, k] <- stats::rpois(n = S, lambda = exp(log_lambda2[, k]))
  }

  # name columns
  colnames(outcome_samples_x1i) <- paste0("x1i[", seq_len(K), "]")
  colnames(outcome_samples_x2i) <- paste0("x2i[", seq_len(K), "]")

  # merge samples with interleaved columns: x1i[1], x2i[1], x1i[2], x2i[2], ...
  outcome_samples <- matrix(NA_integer_, nrow = S, ncol = 2 * K)
  outcome_samples[, seq(1, 2 * K, by = 2)] <- outcome_samples_x1i
  outcome_samples[, seq(2, 2 * K, by = 2)] <- outcome_samples_x2i
  colnames(outcome_samples) <- as.vector(rbind(
    colnames(outcome_samples_x1i),
    colnames(outcome_samples_x2i)
  ))

  return(outcome_samples)
}