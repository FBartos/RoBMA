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
#' @param incorporate_publication_bias whether sampling of new values should incorporate
#' the estimated publication bias (note that selection models do not affect the mean paramater
#' when \code{"terms"} (equal mean parameter under normal vs. weighted likelihood equals different
#' expectation).
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
                         incorporate_publication_bias = TRUE,
                         as_samples = FALSE,
                         ...){

  # some options checked inside BayesTools table directly
  BayesTools::check_char(type, "type", allow_values = c("response", "terms", "terms.scale", "effect"))
  BayesTools::check_bool(as_samples, "as_samples")

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
  is_mods       <- .is_mods(object)
  is_multilevel <- .is_multilevel(object)
  is_PET        <- .is_PET(object)
  is_PEESE      <- .is_PEESE(object)

  ### extract outcome data for convenience
  outcome_data     <- new_data[["outcome"]]
  effect_direction <- .effect_direction(object)
  K                <- nrow(outcome_data)

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
      prior_list  = priors
    ))
    tau_samples      <- exp(log_tau_samples)

  } else {

    # from the overall effect for the remaining models
    tau_samples <- matrix(posterior_samples[,"tau"], ncol = K, nrow = nrow(posterior_samples))

  }

  # splt tau samples into between and within study components if multilevel
  if (is_multilevel) {
    rho_samples <- posterior_samples[,"rho"]
    # deal with computer precision errors from JAGS
    rho_samples[rho_samples>1] <- 1
    rho_samples[rho_samples<0] <- 0
    # tau_within  = tau * sqrt(rho)
    # tau_between = tau * sqrt(1-rho)
    tau_within_samples  <- tau_samples * sqrt(rho_samples)
    tau_between_samples <- tau_samples * sqrt(1 - rho_samples)
  }

   ### return only tau samples if type = "terms.scale" is selected
  if (type == "terms.scale") {
    # TODO: implement summary table for tau samples only....
    return(tau_samples)
  }

  ### get the base mu samples
  if (is_mods) {
    # from effect size regression for models with mods
    mu_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit         = object[["fit"]],
      formula     = .create_fit_formula_list(data = new_data, parameter = "mods"),
      parameter   = "mu",
      data        = .create_fit_formula_data_list(data = new_data, parameter = "mods"),
      prior_list  = .create_fit_formula_prior_list(priors = priors, parameter = "mods")
    ))
  } else {
    # from the overall effect for the remaining models
    mu_samples <- matrix(posterior_samples[,"mu"], ncol = K, nrow = nrow(posterior_samples))
  }

  ### apply effect direction flipping (to match JAGS model)
  # JAGS uses: ifelse(effect_direction == "negative", "- mu", "mu")
  # data are flipped within `.create_fit_data()`, so we flip mu samples to match
  if (effect_direction == "negative") {
    mu_samples <- - mu_samples
  }

  ### add PET adjustment (+ PET * sei[i])
  if (is_PET && incorporate_publication_bias) {
    PET_samples <- posterior_samples[,"PET"]
    for (i in seq_len(K)) {
      mu_samples[,i] <- mu_samples[,i] + PET_samples * outcome_data[i, "se"]
    }
  }

  ### add PEESE adjustment (+ PEESE * sei[i]^2)
  if (is_PEESE && incorporate_publication_bias) {
    PEESE_samples <- posterior_samples[,"PEESE"]
    for (i in seq_len(K)) {
      mu_samples[,i] <- mu_samples[,i] + PEESE_samples * outcome_data[i, "se"]^2
    }
  }

  ### return only mu samples if type = "terms" is selected
  if (type == "terms") {
    # TODO: implement summary table for mu samples only....
    return(mu_samples)
  }

  ### include random-effects if type = `effect` is selected



  # create response prediction if required
  if(type %in% c("response", "effect")){

    if(!incorporate_publication_bias || inherits(object, "NoBMA.reg") || inherits(object, "BiBMA.reg") || (length(priors$bias) == 1 && is.prior.none(priors$bias[[1]]))){

      # predicting responses without selection models does not require incorporating the between-study random-effects
      # (the marginalized and non-marginalized parameterization are equivalent)
      tau_samples <- posterior_samples[,"tau"]
      tau_samples <- .scale(tau_samples, object$add_info[["output_scale"]], model_scale)

      # sample the effects / observed studies
      outcome_samples <- matrix(NA, nrow = nrow(mu_samples), ncol = ncol(mu_samples))
      if (type == "effect"){
        for(i in seq_len(ncol(mu_samples))){
          outcome_samples[,i] <- stats::rnorm(nrow(mu_samples), mu_samples[,i], tau_samples)
        }
      }else if (type == "response"){
        for(i in seq_len(ncol(mu_samples))){
          outcome_samples[,i] <- stats::rnorm(nrow(mu_samples), mu_samples[,i], sqrt(tau_samples^2 + outcome_data[i, "se"]^2))
        }
      }


    }else{

      # required for study ids / crit_x values in selection models
      fit_data <- .fit_data_ss(
        data             = outcome_data,
        priors           = priors,
        effect_direction = effect_direction,
        prior_scale      = object$add_info[["prior_scale"]],
        weighted         = FALSE,
        weighted_type    = FALSE,
        multivariate     = if (same_data) .is_multivariate(object) else FALSE
      )

      # predicting response requires incorporating the between-study random effects if selection models are present
      # (we use approximate selection likelihood which samples the true study effects instead of marginalizing them)
      if(.is_multivariate(object)){

        tau_samples <- posterior_samples[,"tau"]
        rho_samples <- posterior_samples[,"rho"]
        # deal with computer precision errors from JAGS
        rho_samples[rho_samples>1] <- 1
        rho_samples[rho_samples<0] <- 0
        # tau_within  = tau * sqrt(rho)
        # tau_between = tau * sqrt(1-rho)
        tau_within_samples  <- tau_samples * sqrt(rho_samples)
        tau_between_samples <- tau_samples * sqrt(1-rho_samples)
        gamma_samples       <- posterior_samples[,grep("gamma", colnames(posterior_samples)),drop = FALSE]

        tau_between_samples <- .scale(tau_between_samples, object$add_info[["output_scale"]], model_scale)
        tau_within_samples  <- .scale(tau_within_samples,  object$add_info[["output_scale"]], model_scale)

        # incorporate within study heterogeneity into the predictor
        # either estimated for prediction on the same data or integrated over for new data
        if(same_data){
          for(i in seq_len(K)){
            mu_samples[,i] <- mu_samples[,i] + gamma_samples[,fit_data$study_ids[i]] * tau_within_samples
          }
        }else{
          for(i in seq_len(K)){
            mu_samples[,i] <- mu_samples[,i] + stats::rnorm(nrow(mu_samples)) * tau_within_samples
          }
        }


        # tau_between samples work as tau for the final sampling step
        tau_samples <- tau_between_samples

      }else{

        tau_samples  <- posterior_samples[,"tau"]
        tau_samples  <- .scale(tau_samples,  object$add_info[["output_scale"]], model_scale)

      }

      outcome_samples <- matrix(NA, nrow = nrow(mu_samples), ncol = ncol(mu_samples))

      # selection models are sampled separately for increased efficiency
      bias_indicator           <- posterior_samples[,"bias_indicator"]
      weightfunction_indicator <- bias_indicator %in% which(sapply(priors[["bias"]], is.prior.weightfunction))

      # sample the effects / observed studies
      if(type == "effect"){

        for(i in seq_len(ncol(mu_samples))){

          # sample normal models/PET/PEESE
          if(any(!weightfunction_indicator)){
            outcome_samples[!weightfunction_indicator,i] <- stats::rnorm(
              n    = sum(!weightfunction_indicator),
              mean = mu_samples[!weightfunction_indicator,i],
              sd   = tau_samples[!weightfunction_indicator]
            )
          }

          # sample selection models (.rwnorm_true_fast.ss returns the implied random effects for given selection)
          if(any(weightfunction_indicator)){
            outcome_samples[weightfunction_indicator,i] <- .rwnorm_true_fast.ss(
              mean   = mu_samples[weightfunction_indicator,i],
              tau    = tau_samples[weightfunction_indicator],
              se     = outcome_data[i, "se"],
              omega  = posterior_samples[weightfunction_indicator, grep("omega", colnames(posterior_samples)),drop = FALSE],
              crit_x = fit_data$crit_y[, i]
            )
          }
        }

      }else if(type == "response"){

        for(i in seq_len(ncol(mu_samples))){

          # sample normal models/PET/PEESE
          if(any(!weightfunction_indicator)){
            outcome_samples[!weightfunction_indicator,i] <- stats::rnorm(
              n    = sum(!weightfunction_indicator),
              mean = mu_samples[!weightfunction_indicator,i],
              sd   = sqrt(tau_samples[!weightfunction_indicator]^2 + outcome_data[i, "se"]^2)
            )
          }

          # sample selection models
          if(any(weightfunction_indicator)){
            outcome_samples[weightfunction_indicator,i] <- .rwnorm_fast.ss(
              mean   = mu_samples[weightfunction_indicator,i],
              sd     = sqrt(tau_samples[weightfunction_indicator]^2 + outcome_data[i, "se"]^2),
              omega  = posterior_samples[weightfunction_indicator, grep("omega", colnames(posterior_samples)),drop = FALSE],
              crit_x = fit_data$crit_y[, i]
            )
          }
        }

      }

    }

  }else if(type == "terms"){
    # terms only returns the mean for the prediction
    outcome_samples <- mu_samples
  }

  # flip outcome samples back to the original direction
  if (effect_direction == "negative") {
    outcome_samples <- -outcome_samples
  }

  # select conditional estimates
  if(conditional){
    outcome_samples_conditional <- outcome_samples[mu_indicator %in% which(!mu_is_null),, drop=FALSE]
    outcome_samples_conditional <- lapply(1:ncol(outcome_samples_conditional), function(i) {
      .transform_mu(outcome_samples_conditional[,i], from = model_scale, to = output_scale)
    })
    names(outcome_samples_conditional) <- switch(
      type,
      "terms"    = sapply(seq_along(outcome_samples_conditional), function(x) paste0("mu[", x, "]")),
      "effect"   = sapply(seq_along(outcome_samples_conditional), function(x) paste0("theta[", x, "]")),
      "response" = sapply(seq_along(outcome_samples_conditional), function(x) paste0("estimate[", x, "]"))
    )
  }

  # transform the effect sizes (and name the matrix)
  outcome_samples <- lapply(1:ncol(outcome_samples), function(i) {
    .transform_mu(outcome_samples[,i], from = model_scale, to = output_scale)
  })
  names(outcome_samples) <- switch(
    type,
    "terms"    = sapply(seq_along(outcome_samples), function(x) paste0("mu[", x, "]")),
    "effect"   = sapply(seq_along(outcome_samples), function(x) paste0("theta[", x, "]")),
    "response" = sapply(seq_along(outcome_samples), function(x) paste0("estimate[", x, "]"))
  )

  # return only samples if requested
  if(as_samples){
    if(conditional){
      return(outcome_samples_conditional)
    }else{
      return(outcome_samples)
    }
  }

  # obtain estimates tables
  estimates <- BayesTools::ensemble_estimates_table(
    samples    = outcome_samples,
    parameters = names(outcome_samples),
    probs      = probs,
    title      = "Posterior predictions:",
    footnotes  = c(.scale_note_simple(object$add_info[["prior_scale"]], output_scale))
  )

  if(conditional){
    estimates_conditional <- BayesTools::ensemble_estimates_table(
      samples    = outcome_samples_conditional,
      parameters = names(outcome_samples_conditional),
      probs      = probs,
      title      = "Conditional posterior predictions:",
      footnotes  = c(.scale_note_simple(object$add_info[["prior_scale"]], output_scale))
    )
  }

  # create the output object
  output <- list(
    call       = object[["call"]],
    title      = .object_title(object),
    estimates  = estimates,
    footnotes  = c(.scale_note_simple(object$add_info[["prior_scale"]], output_scale))
  )

  if(conditional){
    output$estimates_conditional <- estimates_conditional
  }

  class(output) <- "summary.RoBMA"
  attr(output, "type") <- "ensemble"

  return(output)
}





