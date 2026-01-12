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

  ### extract posterior samples (and obtain conditional indicator)
  posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
  priors            <- object[["priors"]]

  ### extract structural information about the model
  is_mods           <- .is_mods(object)
  is_multilevel     <- .is_multilevel(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)

  ### extract outcome data and fit data for convenience
  outcome_data     <- new_data[["outcome"]]
  effect_direction <- .effect_direction(object)
  fit_data         <- .create_fit_data(data = new_data, priors = priors)

  # outcome dimensions
  K  <- nrow(outcome_data)
  S  <- nrow(posterior_samples)

  ### obtain mu and tau samples (following the same structure as the JAGS syntax created by .create_model_syntax)
  # JAGS constructs tau_estimate as:
  # tau                                   (if not multilevel and not scale)
  # tau = exp(log_tau)                    (if scale)
  # tau_within = tau * sqrt(rho)          (if multilevel)
  # tau_between = tau * sqrt(1-rho)
  #
  # JAGS constructs mu_estimate as:
  # mu +
  #  + gamma[study_ids[i]] * tau_between  (if multilevel)
  #  + PET * sei[i]                       (if PET & incorporate_publication_bias)
  #  + PEESE * sei[i]^2                   (if PEESE & incorporate_publication_bias)
  # - mu and gamma are flipped if effect_direction is "negative"
  #   (data are flipped in .create_fit_data; to keep PET and PEESE coefficients positive)
  # - multilevel effects are added later in case only terms are predicted (also require tau_between samples)

  ### obtain tau samples
  if (.is_scale(object)) {

    # from effect size regression for models with mods
    log_tau_samples  <- t(BayesTools::JAGS_evaluate_formula(
      fit         = object[["fit"]],
      formula     = attr(object[["data"]][["scale"]], "formula"),
      parameter   = "log_tau",
      data        = new_data[["scale"]],
      prior_list  = priors[["scale"]]
    ))
    tau_samples      <- exp(log_tau_samples)

  } else {

    # from the overall effect for the remaining models
    tau_samples <- matrix(posterior_samples[,"tau"], ncol = K, nrow = S)

  }

  ### split tau samples into between (study-level) and within (estimate-level) components
  if (is_multilevel) {
    # multilevel models have explicit split into the components via rho
    rho_samples <- posterior_samples[,"rho"]
    # deal with computer precision errors from JAGS
    rho_samples[rho_samples > 1] <- 1
    rho_samples[rho_samples < 0] <- 0
    # tau_within  = tau * sqrt(rho)
    # tau_between = tau * sqrt(1-rho)
    tau_within_samples  <- tau_samples * sqrt(rho_samples)
    tau_between_samples <- tau_samples * sqrt(1 - rho_samples)
  } else {
    # simple random-effects model only model estimate-level heterogeneity
    # (specify tau_between_samples for simplification of following code)
    tau_within_samples  <- tau_samples
    tau_between_samples <- matrix(0, ncol = K, nrow = S)
  }

  ### return only tau samples if type = "terms.scale" is selected
  if (type == "terms.scale") {
    # rename samples
    colnames(tau_samples) <- paste0("tau[", seq_len(K), "]")

    if (as_samples) {
      return(mu_samples)
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

  ### get the base mu samples
  if (is_mods) {
    # from effect size regression for models with mods
    mu_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit         = object[["fit"]],
      formula     = attr(object[["data"]][["mods"]], "formula"),
      parameter   = "mu",
      data        = new_data[["mods"]],
      prior_list  = priors[["mods"]]
    ))
  } else {
    # from the overall effect for the remaining models
    mu_samples <- matrix(posterior_samples[,"mu"], ncol = K, nrow = S)
  }

  ### apply effect direction flipping (to match JAGS model)
  # JAGS uses: ifelse(effect_direction == "negative", "- mu", "mu")
  # data are flipped within `.create_fit_data()`, so we flip mu samples to match
  if (effect_direction == "negative") {
    mu_samples <- - mu_samples
  }

  ### add PET adjustment (+ PET * sei[i])
  if (is_PET && !bias_adjusted) {
    PET_samples <- posterior_samples[,"PET"]
    for (i in seq_len(K)) {
      mu_samples[,i] <- mu_samples[,i] + PET_samples * outcome_data[["sei"]][i]
    }
  }

  ### add PEESE adjustment (+ PEESE * sei[i]^2)
  if (is_PEESE && !bias_adjusted) {
    PEESE_samples <- posterior_samples[,"PEESE"]
    for (i in seq_len(K)) {
      mu_samples[,i] <- mu_samples[,i] + PEESE_samples * outcome_data[["sei"]][i]^2
    }
  }

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

  ### include 3-level (study-level) random-effects
  # study-level effects are always sampled -- not marginalized
  if (is_multilevel) {
    if (same_data) {
      # when same data are used we apply the existing random effects samples
      gamma_samples <- posterior_samples[,paste0("gamma[", 1:max(fit_data[["study_ids"]]),"]")]
      for (i in seq_len(K)) {
        mu_samples[,i] <- mu_samples[,i] + ifelse(effect_direction == "negative", -1, 1) * gamma_samples[,fit_data[["study_ids"]][i]] * tau_between_samples[,i]
      }
    } else {
      # when new data are used we marginalize over the distribution of random effects (they are scaled normal)
      for (i in seq_len(K)) {
        mu_samples[,i] <- mu_samples[,i] + ifelse(effect_direction == "negative", -1, 1) * stats::rnorm(S, 0, 1) * tau_between_samples[,i]
      }
    }
  }

  ### create true-effect prediction (i.e., the latent estimate-level random-effects)
  if (type == "effect") {
    # the estimate-level random-effects are always marginalized
    # as such, we need to use empirical Bayes to obtain the true study effects
    # (corresponds to BLUPs in metafor)
    # (the selection weights cancel out and the same formula applies for all models)

    effect_sizes    <- matrix(outcome_data[["yi"]],  ncol = K, nrow = S)
    standard_errors <- matrix(outcome_data[["sei"]], ncol = K, nrow = S)

    # get the shrinkage matrix
    lambda <- tau_within_samples^2 / (tau_within_samples^2 + standard_errors^2)

    # compute BLUPs
    true_effects_samples <- lambda * effect_sizes + (1 - lambda) * mu_samples

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

      # the estimate-level estimates are not marginalized for GLMMs
      # include the sampled random effects or marginalized over them for new data
      if (same_data) {
        # when same data are used we apply the existing random effects samples
        # (no need to check direction -- GLMMs always set to positive)
        for (i in seq_len(K)) {
          mu_samples[,i] <- mu_samples[,i] + posterior_samples[[paste0("theta[",i,"]")]] * tau_within_samples[,i]
        }
      } else {
        # when new data are used we marginalize over the distribution of random effects (they are scaled normal)
        # (no need to check direction -- GLMMs always set to positive)
        for (i in seq_len(K)) {
          mu_samples[,i] <- mu_samples[,i] + stats::rnorm(S, 0, 1) * tau_within_samples[,i]
        }
      }

      if (outcome_type == "bin") {

        ### incorporate with base-rate
        logit_p1_samples <- matrix(NA, ncol = K, nrow = S)
        logit_p2_samples <- matrix(NA, ncol = K, nrow = S)

        for (i in seq_len(K)) {
          logit_p1_samples[,i] <- .logit(posterior_samples[[paste0("pi[", i, "]")]]) + 0.5 * mu_samples[,i]
          logit_p2_samples[,i] <- .logit(posterior_samples[[paste0("pi[", i, "]")]]) - 0.5 * mu_samples[,i]
        }

        ### sample outcome (observed number of successes in each group)
        outcome_samples_ai <- matrix(NA, ncol = K, nrow = S)
        outcome_samples_ci <- matrix(NA, ncol = K, nrow = S)

        for (i in seq_len(K)) {
          outcome_samples_ai[,i] <- stats::rbinom(n = S, size = outcome_data[["n1i"]][i], prob = .inv_logit(logit_p1_samples[,i]))
          outcome_samples_ci[,i] <- stats::rbinom(n = S, size = outcome_data[["n2i"]][i], prob = .inv_logit(logit_p2_samples[,i]))
        }

        # rename samples
        colnames(outcome_samples_ai) <- paste0("ai[", seq_len(K), "]")
        colnames(outcome_samples_ci) <- paste0("ci[", seq_len(K), "]")

        # merge samples
        outcome_samples <- matrix(NA, ncol = 2*K, nrow = S)
        outcome_samples[,(seq_len(K)) * 2 - 1] <- outcome_samples_ai
        outcome_samples[,(seq_len(K)) * 2    ] <- outcome_samples_ci
        colnames(outcome_samples)[,(seq_len(K)) * 2 - 1] <- colnames(outcome_samples_ai)
        colnames(outcome_samples)[,(seq_len(K)) * 2    ] <- colnames(outcome_samples_ci)

      } else if (outcome_type == "pois") {

        ### incorporate with log-rate
        log_r1_samples <- matrix(NA, ncol = K, nrow = S)
        log_r2_samples <- matrix(NA, ncol = K, nrow = S)

        for (i in seq_len(K)) {
          log_r1_samples[,i] <- posterior_samples[[paste0("phi[", i, "]")]] + 0.5 * mu_samples[,i] + log(outcome_data[["t1i"]][i])
          log_r2_samples[,i] <- posterior_samples[[paste0("phi[", i, "]")]] - 0.5 * mu_samples[,i] + log(outcome_data[["t2i"]][i])
        }

        ### sample outcome (observed number of events in each group)
        outcome_samples_x1i <- matrix(NA, ncol = K, nrow = S)
        outcome_samples_x2i <- matrix(NA, ncol = K, nrow = S)

        for (i in seq_len(K)) {
          outcome_samples_x1i[,i] <- stats::rpois(n = S, lambda = exp(log_r1_samples[,i]))
          outcome_samples_x2i[,i] <- stats::rpois(n = S, lambda = exp(log_r2_samples[,i]))
        }

        # rename samples
        colnames(outcome_samples_x1i) <- paste0("x1i[", seq_len(K), "]")
        colnames(outcome_samples_x2i) <- paste0("x2i[", seq_len(K), "]")

        # merge samples
        outcome_samples <- matrix(NA, ncol = 2*K, nrow = S)
        outcome_samples[,(seq_len(K)) * 2 - 1] <- outcome_samples_x1i
        outcome_samples[,(seq_len(K)) * 2    ] <- outcome_samples_x2i
        colnames(outcome_samples)[,(seq_len(K)) * 2 - 1] <- colnames(outcome_samples_x1i)
        colnames(outcome_samples)[,(seq_len(K)) * 2    ] <- colnames(outcome_samples_x2i)

      }


    } else if (outcome_type == "norm") {

      # outcome samples holder (observed effect sizes)
      outcome_samples <- matrix(NA, ncol = K, nrow = S)

      # normal outcome models need to dispatch between
      # normal and selection models
      # (bias_adjusted results return the samples if no bias was present)

      if (bias_adjusted) {
        # sample from the outcome distribution directly
        # (i.e., as if publication bias did not exist)
        for (i in seq_len(K)) {
          outcome_samples[,i] <- stats::rnorm(
            n    = S,
            mean = mu_samples[,i],
            sd   = sqrt(tau_within_samples[,i]^2 + outcome_data[["sei"]][i]^2)
          )
        }

      } else {

        if (is_weightfunction) {
          # sample from the weighted distribution for selection models
          for (i in seq_len(K)) {
            outcome_samples[,i] <- .rwnorm_fast.ss(
              mean   = mu_samples[,i],
              sd     = sqrt(tau_within_samples[,i]^2 + outcome_data[["se"]][i]^2),
              omega  = posterior_samples[, grep("omega", colnames(posterior_samples)),drop = FALSE],
              crit_x = fit_data$crit_yi[, i]
            )
          }
        } else {
          # sample from normal distirbution for remaining models
          for (i in seq_len(K)) {
            outcome_samples[,i] <- stats::rnorm(
              n    = S,
              mean = mu_samples[,i],
              sd   = sqrt(tau_within_samples[,i]^2 + outcome_data[["sei"]][i]^2)
            )
          }
        }
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
summary.brma.predict <- function(x, ...) {
  print(summary.brma.predict(x, ...))
}

#' @export
print.brma.predict <- function(x, ...) {

  cat("\n")
  print(x[["summary"]])
  cat("\n")

  return(invisible(x))
}

