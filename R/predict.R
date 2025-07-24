#' @title Predict method for Robust Bayesian Meta-Analysis Fits
#'
#' @description \code{predict.RoBMA} predicts values based on the RoBMA model.
#' Only available for normal-normal models estimated using the spike-and-slab
#' algorithm (i.e., \code{algorithm = "ss"}).
#'
#' @inheritParams summary.RoBMA
#' @inheritParams pooled_effect
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
#' The conditional estimate is calculated conditional on the presence of the effect
#' (in meta-analysis) or the intercept (in meta-regression).
#'
#' @examples \dontrun{
#' require(metafor)
#' dat <- escalc(measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
#'               data = dat.bcg)
#'
#' # fit meta-regression
#' robma_dat <- data.frame(
#'   logOR  = dat$yi,
#'   se     = sqrt(dat$vi),
#'   ablat  = dat$ablat,
#'   alloc  = dat$alloc
#' )
#'
#' fit <- NoBMA.reg(~ ablat + alloc, data = robma_dat,
#'                  seed = 1, algorithm = "ss", parallel = TRUE)
#'
#' # prediction for the mean effect, prediction interval, and posterior predictive
#' mean_effect          <- predict(fit, type = "terms")
#' prediction_interval  <- predict(fit, type = "effect")
#' posterior_predictive <- predict(fit, type = "response")
#'
#' # visualize the estimates vs predictions
#' plot(NA, type = "n", xlim = c(-3, 3), ylim = c(0, nrow(robma_dat) + 1),
#'      xlab = "logOR", ylab = "Observation", las = 1)
#' points(robma_dat$logOR, seq_along(robma_dat$logOR), cex = 2, pch = 16)
#' points(mean_effect$estimates$Mean, seq_along(robma_dat$logOR) + 0.2,
#'        cex = 1.5, pch = 16, col = "blue")
#' sapply(seq_along(robma_dat$logOR), function(i){
#'   lines(c(robma_dat$logOR[i] - 1.96 * robma_dat$se[i],
#'           robma_dat$logOR[i] + 1.96 * robma_dat$se[i]),
#'         c(i, i), lwd = 2)
#'   lines(c(mean_effect$estimates[i,"0.025"],
#'           mean_effect$estimates[i,"0.975"]),
#'         c(i + 0.2, i + 0.2), lwd = 2, col = "blue")
#'   lines(c(prediction_interval$estimates[i,"0.025"],
#'           prediction_interval$estimates[i,"0.975"]),
#'         c(i + 0.3, i + 0.3), lwd = 2, lty = 2, col = "blue")
#'   lines(c(posterior_predictive$estimates[i,"0.025"],
#'           posterior_predictive$estimates[i,"0.975"]),
#'         c(i + 0.4, i + 0.4), lwd = 2, lty = 3, col = "blue")
#' })
#' legend("bottomright", col = c("black", rep("blue", 3)),
#'        lwd = 2, lty = c(1, 1, 2, 3), bty = "n",
#'        legend = c("Observed + CI", "Predicted + CI",
#'                   "Prediction Int.", "Sampling Int.")
#'       )
#'
#'  # prediction across a lattitude
#'  fit2 <- NoBMA.reg(~ ablat, data = robma_dat,
#'                    seed = 1, algorithm = "ss", parallel = TRUE)
#'
#'  new_df <- data.frame(
#'    logOR  = 0,   # only relevant for "response" (not plotted here)
#'    se     = 0.1, # only relevant for "response" (not plotted here)
#'    ablat  = 10:60
#'  )
#'  new_mean_effect          <- predict(fit2, newdata = new_df, type = "terms")
#'  new_prediction_interval  <- predict(fit2, newdata = new_df, type = "effect")
#'
#'  # create bubble plot
#'  plot(robma_dat$ablat, robma_dat$logOR,  ylim = c(-2, 1),
#'       xlab = "Latitude", ylab = "logOR", las = 1, cex = 0.3/robma_dat$se, pch = 16)
#'  polygon(c(10:60, rev(10:60)),
#'          c(new_mean_effect$estimates[,"0.025"],
#'            rev(new_mean_effect$estimates[,"0.975"])),
#'          col = rgb(0, 0, 1, alpha = 0.2), border = NA)
#'  polygon(c(10:60, rev(10:60)),
#'          c(new_prediction_interval$estimates[,"0.025"],
#'            rev(new_prediction_interval$estimates[,"0.975"])),
#'          col = rgb(0, 0, 1, alpha = 0.2), border = NA)
#' }
#'
#' @return \code{pooled_effect} returns a list of tables of class 'BayesTools_table'.
#' @export
predict.RoBMA <- function(object, newdata = NULL, type = "response",
                          conditional = FALSE, output_scale = NULL, probs = c(.025, .975),
                          incorporate_publication_bias = TRUE, as_samples = FALSE, ...){

  # some options checked inside BayesTools table directly
  BayesTools::check_char(type, "type", allow_values = c("response", "terms", "effect"))
  if((is.RoBMA.reg(object) || is.NoBMA.reg(object) || is.BiBMA.reg(object)) && !(is.null(newdata) || is.data.frame(newdata)))
    stop("The 'newdata' argument must be a data frame or NULL when performing prediction for meta-regression models.", call. = FALSE)
  if(!(is.RoBMA.reg(object) || is.NoBMA.reg(object) || is.BiBMA.reg(object)) && !(is.null(newdata) || is.list(newdata)))
    stop("The 'newdata' argument must be a list or NULL when performing prediction for meta-analytic models.", call. = FALSE)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(incorporate_publication_bias, "incorporate_publication_bias")
  BayesTools::check_bool(as_samples, "as_samples")
  dots <- list(...)
  if(object[["add_info"]][["algorithm"]] != "ss")
    stop("Predictions can only be computed for spike and slab models.")
  if(inherits(object, "BiBMA") || inherits(object, "BiBMA.reg"))
    stop("The true effects can only be computed for normal-normal (NoBMA / RoBMA) models.")


  # get the model fitting scale
  if (is.BiBMA(object)) {
    model_scale <- "logOR"
  } else {
    model_scale <- object$add_info[["effect_measure"]]
  }
  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale) != "y"){
    stop("Models estimated using the general effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale)
  }

  # predict for the current data if no data are specified
  if(is.null(newdata)){

    # an existing data are used
    same_data <- TRUE
    # - when predicting effects/outcomes for multilevel outcomes, use estimated gamma levels

    # dispatch between meta-regression / meta-analysis input
    if(inherits(object, "RoBMA.reg") || inherits(object, "NoBMA.reg") || inherits(object, "BiBMA.reg")){
      newdata.predictors <- do.call(cbind.data.frame, object$data[["predictors"]])
      newdata.outcome    <- object$data[["outcome"]]
    }else{
      newdata.outcome <- object$data
    }

  }else{

    # a new data are used
    same_data <- FALSE
    # - when predicting effects/outcomes for multilevel outcomes, average across possible tau_within levels via sampling

    # dispatch between meta-regression / meta-analysis input
    if(inherits(object, "RoBMA.reg") || inherits(object, "NoBMA.reg") || inherits(object, "BiBMA.reg")){

      RoBMA.options(check_scaling = FALSE)
      newdata <- .combine_data.reg(
        formula                = object[["formula"]],
        data                   = newdata,
        standardize_predictors = FALSE,
        transformation         = .transformation_invar(model_scale),
        study_names            = NULL,
        study_ids              = NULL
      )
      RoBMA.options(check_scaling = TRUE)

      # manually standardize predictors to match the original scaling
      if(object$add_info[["standardize_predictors"]]){
        variables_info <- attr(object$data[["predictors"]], "variables_info")
        for(i in seq_along(variables_info)){
          if(variables_info[[i]][["type"]] == "continuous"){
            newdata[["predictors"]][[i]] <- (newdata[["predictors"]][[i]] - variables_info[[i]][["mean"]]) / variables_info[[i]][["sd"]]
          }
        }
      }

      newdata.predictors <- do.call(cbind.data.frame, newdata[["predictors"]])
      newdata.outcome    <- newdata[["outcome"]]

    }else{

      newdata.outcome <- combine_data(d = newdata[["d"]], r = newdata[["r"]], z = newdata[["z"]], logOR = newdata[["logOR"]],
                                      OR = newdata[["OR"]], t = newdata[["t"]], y = newdata[["y"]], se = newdata[["se"]],
                                      v = newdata[["v"]], n = newdata[["n"]], lCI = newdata[["lCI"]], uCI = newdata[["uCI"]],
                                      study_names = NULL, study_ids = NULL, weight = NULL, data = NULL, transformation = .transformation_invar(model_scale))

    }
  }


  # extract posterior samples (and obtain conditional indicator)
  posterior_samples <- suppressWarnings(coda::as.mcmc(object[["model"]][["fit"]]))
  priors            <- object[["priors"]]

  # obtain the (study-specific) mu estimate
  # meta-regression and meta-analysis separately
  if(inherits(object, "RoBMA.reg") || inherits(object, "NoBMA.reg") || inherits(object, "BiBMA.reg")){

    mu_samples  <- t(BayesTools::JAGS_evaluate_formula(
      fit         = object$model$fit,
      formula     = object$formula,
      parameter   = "mu",
      data        = newdata.predictors,
      prior_list  = attr(object$model$fit, "prior_list")
    ))
    mu_is_null   <- attr(object[["model"]]$priors$terms[["intercept"]], "components") == "null"
    mu_indicator <- posterior_samples[,"mu_intercept_indicator"]

  }else{

    mu_samples   <- matrix(posterior_samples[,"mu"],  ncol = nrow(newdata.outcome), nrow = nrow(posterior_samples))
    mu_is_null   <- attr(object[["model"]]$priors$mu, "components") == "null"
    mu_indicator <- posterior_samples[,"mu_indicator"]

  }

  # transform the samples to the model fitting scale (the data are at the model fitting scale)
  mu_samples  <- .scale(mu_samples,  object$add_info[["output_scale"]], model_scale)

  # add conditional warnings
  if(conditional){
    if(sum(mu_indicator %in% which(!mu_is_null)) < 100)
      warning(gettextf("There is only a very small number of posterior samples (%1s) assuming presence of the effect. The resulting estimates are not reliable.",
                       sum(mu_indicator %in% which(!mu_is_null))), call. = FALSE, immediate. = TRUE)
    if(sum(mu_indicator %in% which(!mu_is_null)) <= 2)
      stop("Less or equal to 2 posterior samples assuming presence of the effects. The estimates could not be computed.")
  }

  # add PET/PEESE adjustment
  priors_bias  <- priors[["bias"]]
  if(incorporate_publication_bias){
    if(any(sapply(priors_bias, is.prior.PET))){
      PET_samples <- posterior_samples[,"PET"]
      # PET is scale invariant (no-scaling needed)
    }else{
      PET_samples <- rep(0, nrow(posterior_samples))
    }
    if(any(sapply(priors_bias, is.prior.PEESE))){
      PEESE_samples <- posterior_samples[,"PEESE"]
      # PEESE scales with the inverse
      PEESE_samples <- .scale(PEESE_samples, model_scale, object$add_info[["output_scale"]])
    }else{
      PEESE_samples <- rep(0, nrow(posterior_samples))
    }

    for(i in seq_len(ncol(mu_samples))){
      mu_samples[,i] <- mu_samples[,i] + PET_samples * newdata.outcome[i,"se"] + PEESE_samples * newdata.outcome[i,"se"]^2
    }
  }


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
          outcome_samples[,i] <- stats::rnorm(nrow(mu_samples), mu_samples[,i], sqrt(tau_samples^2 + newdata.outcome[i,"se"]^2))
        }
      }


    }else{

      # required for study ids / crit_x values in selection models
      fit_data <- .fit_data_ss(
        data             = newdata.outcome,
        priors           = priors,
        effect_direction = object$add_info[["effect_direction"]],
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
          for(i in seq_len(nrow(newdata.outcome))){
            mu_samples[,i] <- mu_samples[,i] + gamma_samples[,fit_data$study_ids[i]] * tau_within_samples
          }
        }else{
          for(i in seq_len(nrow(newdata.outcome))){
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

          # sample selection models (.rwnorm_predict_true_fast returns the implied random effects for given selection)
          if(any(weightfunction_indicator)){
            outcome_samples[weightfunction_indicator,i] <- .rwnorm_predict_true_fast(
              mean   = mu_samples[weightfunction_indicator,i],
              tau    = tau_samples[weightfunction_indicator],
              se     = newdata.outcome[i,"se"],
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
              sd   = sqrt(tau_samples[!weightfunction_indicator]^2 + newdata.outcome[i,"se"]^2)
            )
          }

          # sample selection models
          if(any(weightfunction_indicator)){
            outcome_samples[weightfunction_indicator,i] <- .rwnorm_predict_fast(
              mean   = mu_samples[weightfunction_indicator,i],
              sd     = sqrt(tau_samples[weightfunction_indicator]^2 + newdata.outcome[i,"se"]^2),
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
    footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale))
  )

  if(conditional){
    estimates_conditional <- BayesTools::ensemble_estimates_table(
      samples    = outcome_samples_conditional,
      parameters = names(outcome_samples_conditional),
      probs      = probs,
      title      = "Conditional posterior predictions:",
      footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale))
    )
  }

  # create the output object
  output <- list(
    call       = object[["call"]],
    title      = .object_title(object),
    estimates  = estimates,
    footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale))
  )

  if(conditional){
    output$estimates_conditional <- estimates_conditional
  }

  class(output) <- "summary.RoBMA"
  attr(output, "type") <- "ensemble"

  return(output)
}
