#' @title Compute pooled effect size
#'
#' @description \code{pooled_effect} computes the pooled effect size
#' for a fitted RoBMA.reg and BiBMA.reg object. Only available for models
#' estimated using the spike-and-slab algorithm (i.e., \code{algorithm = "ss"}).
#'
#' @inheritParams summary.RoBMA
#' @param as_samples whether posterior samples instead of a summary table should
#' be returned. Defaults to \code{FALSE}.
#'
#' @details
#' The meta-regression specification results in the intercept corresponding
#' to the adjusted effect estimate (i.e., adjusting for the effect of moderators).
#' In case of moderators inbalance, the adjusted effect estimate might not be
#' representative of the sample of studies. The pooled effect size function averages
#' the effect size estimate across the moderators proportionately to the
#' moderators levels observed in the data set. Note that there is no Bayes factor
#' test for the presence of the pooled effect (the summary function provides the
#' adjusted effect and the test for the presence of the adjusted effect).
#'
#' The conditional estimate is calculated conditional on the presence of the adjusted
#' effect (i.e., the intercept).
#'
#'
#' @return \code{pooled_effect} returns a list of tables of class 'BayesTools_table'.
#' @seealso [adjusted_effect()]
#' @export
pooled_effect <- function(object, conditional = FALSE, output_scale = NULL, probs = c(.025, .975), ...) {
  return(.compute_effect(object, conditional = conditional, output_scale = output_scale, probs = probs, type = "pooled", ...))
}


#' @title Compute adjusted effect size
#'
#' @description \code{adjusted_effect} computes the adjusted effect size
#' for a fitted RoBMA.reg and BiBMA.reg object. Only available for models
#' estimated using the spike-and-slab algorithm (i.e., \code{algorithm = "ss"}).
#'
#' @inheritParams summary.RoBMA
#' @inheritParams pooled_effect
#'
#' @details
#' Non-default meta-regression specification (i.e., using treatment contrasts for
#' predictors) might results in the intercept corresponding
#' to the effect estimate in the baseline group. (i.e., adjusting for the effect of moderators).
#' The adjusted effect size function averages the effect size estimate across the moderators
#' levels. Note that there is no Bayes factor test for the presence of the adjusted effect
#' (the summary function provides the effect estimate in the baseline group and the test for the
#' presence of the effect in the baseline group if a treatment contrasts are specified).
#'
#' The conditional estimate is calculated conditional on the presence of the baseline group
#' effect (i.e., the intercept).
#'
#'
#' @return \code{pooled_effect} returns a list of tables of class 'BayesTools_table'.
#' @seealso [pooled_effect()]
#' @export
adjusted_effect <- function(object, conditional = FALSE, output_scale = NULL, probs = c(.025, .975), ...) {
  return(.compute_effect(object, conditional = conditional, output_scale = output_scale, probs = probs, type = "adjusted", ...))
}

.compute_effect <- function(object, conditional = FALSE, output_scale = NULL, probs = c(.025, .975), type = "pooled", as_samples = FALSE) {

  .check_is_any_RoBMA_object(object)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_char(type, "type", allow_values = c("pooled", "adjusted"))
  BayesTools::check_bool(as_samples, "as_samples")
  if(!(inherits(object, "RoBMA.reg") || inherits(object, "NoBMA.reg") || inherits(object, "BiBMA.reg")))
    stop("The pooled effect size can only be computed for regression models.")
  if(object[["add_info"]][["algorithm"]] != "ss")
    stop("The pooled effect size can only be computed for spike and slab models.")

  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale) != "y"){
    stop("Models estimated using the generall effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale)
  }

  # get posterior samples
  posterior_samples <- suppressWarnings(coda::as.mcmc(object[["model"]][["fit"]]))

  ### compute pooled/adjusted effect
  if (type == "pooled") {
    # use sample proportions
    design_data <- do.call(cbind.data.frame, object[["data"]][["predictors"]])
  } else if (type == "adjusted") {
    # use equivalent proportions
    predictors  <- object[["data"]][["predictors"]]
    design_data <- lapply(names(predictors), function(predictor) {
      if (attr(predictors,"variables_info")[[predictor]][["type"]] == "continuous") {
        return(mean(predictors[[predictor]]))
      } else if (attr(predictors,"variables_info")[[predictor]][["type"]] == "factor") {
        return(unique(predictors[[predictor]]))
      }
    })
    names(design_data) <- names(predictors)
    design_data <- data.frame(expand.grid(design_data))
  }

  ### compute adjusted effect
  pooled_estimate <- BayesTools::JAGS_evaluate_formula(
    fit         = posterior_samples,
    formula     = object[["formula"]],
    parameter   = "mu",
    data        = design_data,
    prior_list  = attr(object[["model"]][["fit"]], "prior_list")
  )

  # average across the design matrix multiplied by the samples
  pooled_estimate <- rowMeans(t(pooled_estimate))

  ### compute prediction interval
  # simulate predictions for PI
  if(is.null(object$add_info[["seed"]])){
    set.seed(1)
  }else{
    set.seed(object$add_info[["seed"]])
  }

  if(BayesTools::is.prior.point(object[["model"]]$priors[["tau"]])){
    tau <- object[["model"]]$priors[["tau"]]$parameters[["location"]]
  }else{
    tau <- posterior_samples[,"tau"]
  }

  # these needs to be simulated as the posteriors can be non-normal
  predictions <- stats::rnorm(length(pooled_estimate), mean = pooled_estimate, sd = tau)

  # compute conditional estimates
  mu_is_null   <- attr(object[["model"]]$priors$terms[["intercept"]], "components") == "null"
  mu_indicator <- posterior_samples[,"mu_intercept_indicator"]

  pooled_estimate_conditional <- pooled_estimate[mu_indicator %in% which(!mu_is_null)]
  predictions_conditional     <- predictions[mu_indicator %in% which(!mu_is_null)]

  # return samples if requested
  if (as_samples){
    return(list(
      estimate    = .transform_mu(pooled_estimate, from = object$add_info[["prior_scale"]], to = output_scale),
      predictions = .transform_mu(predictions,     from = object$add_info[["prior_scale"]], to = output_scale),
      estimate_conditional    = .transform_mu(pooled_estimate_conditional, from = object$add_info[["prior_scale"]], to = output_scale),
      predictions_conditional = .transform_mu(predictions_conditional,     from = object$add_info[["prior_scale"]], to = output_scale)
    ))
  }

  # obtain estimates tables
  estimates <- BayesTools::ensemble_estimates_table(
    samples    = list(
      estimate  = .transform_mu(pooled_estimate, from = object$add_info[["prior_scale"]], to = output_scale),
      PI        = .transform_mu(predictions,     from = object$add_info[["prior_scale"]], to = output_scale)
    ),
    parameters = c("estimate","PI"),
    probs      = probs,
    title      = if (type == "pooled") "Model-averaged pooled effect estimate:" else "Model-averaged adjusted effect estimate:",
    warnings   = .collect_errors_and_warnings(object)
  )

  estimates_conditional <- BayesTools::ensemble_estimates_table(
    samples    = list(
      estimate  = .transform_mu(pooled_estimate_conditional, from = object$add_info[["prior_scale"]], to = output_scale),
      PI        = .transform_mu(predictions_conditional,     from = object$add_info[["prior_scale"]], to = output_scale)
    ),
    parameters = c("estimate","PI"),
    probs      = probs,
    title      = if (type == "pooled") "Conditional pooled effect estimate:" else "Conditional adjusted effect estimate:",
    warnings   = .collect_errors_and_warnings(object),
    footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale))
  )

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


#' @title Compute estimated true effect sizes
#'
#' @description \code{true_effects} computes the estimated true effect size
#' for a fitted RoBMA object. These estimates correspond to the frequentist
#' "Best Linear Unbiased Predictions (BLUPs)". Only available for normal-normal models
#' estimated using the spike-and-slab algorithm (i.e., \code{algorithm = "ss"}).
#'
#' @inheritParams summary.RoBMA
#' @inheritParams pooled_effect
#'
#' @details
#' The conditional estimate is calculated conditional on the presence of the effect
#' (in meta-analysis) or the intercept (in meta-regression).
#'
#' @return \code{pooled_effect} returns a list of tables of class 'BayesTools_table'.
#' @export
true_effects <- function(object, conditional = FALSE, output_scale = NULL, probs = c(.025, .975), as_samples = FALSE){

  .check_is_any_RoBMA_object(object)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(as_samples, "as_samples")
  if(inherits(object, "BiBMA") || inherits(object, "BiBMA.reg"))
    stop("The true effects can only be computed for normal-normal (NoBMA / RoBMA) models.")
  if(object[["add_info"]][["algorithm"]] != "ss")
    stop("The true effects can only be computed for spike and slab models.")

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

  # extract posterior samples (and obtain conditional indicator)
  posterior_samples <- suppressWarnings(coda::as.mcmc(object[["model"]][["fit"]]))

  # obtain the (study-specific) mu estimate
  # meta-regression and meta-analysis separately
  if(inherits(object, "RoBMA.reg") || inherits(object, "NoBMA.reg") || inherits(object, "BiBMA.reg")){

    mu_samples  <- BayesTools::JAGS_evaluate_formula(
      fit         = object$model$fit,
      formula     = object$formula,
      parameter   = "mu",
      data        = do.call(cbind.data.frame, object$data$predictors),
      prior_list  = attr(object$model$fit, "prior_list")
    )
    tau_samples  <- matrix(posterior_samples[,"tau"], byrow = TRUE, nrow = nrow(object$data$outcome), ncol = nrow(posterior_samples))

    mu_is_null   <- attr(object[["model"]]$priors$terms[["intercept"]], "components") == "null"
    mu_indicator <- posterior_samples[,"mu_intercept_indicator"]

    effect_size    <- matrix(object$data$outcome[,"y"],  nrow = nrow(object$data$outcome), ncol = nrow(posterior_samples))
    standard_error <- matrix(object$data$outcome[,"se"], nrow = nrow(object$data$outcome), ncol = nrow(posterior_samples))

  }else{

    mu_samples   <- matrix(posterior_samples[,"mu"],  byrow = TRUE, nrow = nrow(object$data), ncol = nrow(posterior_samples))
    tau_samples  <- matrix(posterior_samples[,"tau"], byrow = TRUE, nrow = nrow(object$data), ncol = nrow(posterior_samples))

    mu_is_null   <- attr(object[["model"]]$priors$mu, "components") == "null"
    mu_indicator <- posterior_samples[,"mu_indicator"]

    effect_size    <- matrix(object$data[,"y"],  nrow = nrow(object$data), ncol = nrow(posterior_samples))
    standard_error <- matrix(object$data[,"se"], nrow = nrow(object$data), ncol = nrow(posterior_samples))

  }

  # transform the samples to the model fitting scale (the data are at the model fitting scale)
  mu_samples  <- .scale(mu_samples,  object$add_info[["output_scale"]], model_scale)
  tau_samples <- .scale(tau_samples, object$add_info[["output_scale"]], model_scale)

  # add conditional warnings
  if(conditional){
    if(sum(mu_indicator %in% which(!mu_is_null)) < 100)
      warning(gettextf("There is only a very small number of posterior samples (%1s) assuming presence of the effect. The resulting estimates are not reliable.",
                       sum(mu_indicator %in% which(!mu_is_null))), call. = FALSE, immediate. = TRUE)
    if(sum(mu_indicator %in% which(!mu_is_null)) <= 2)
      stop("Less or equal to 2 posterior samples assuming presence of the effects. The estimates could not be computed.")
  }

  # get the shrinkage matrix
  lambda <- tau_samples^2 / (tau_samples^2 + standard_error^2)

  # get the blups matrix
  true_effects_samples <- lambda * effect_size + (1 - lambda) * mu_samples

  # select conditional estimates
  if(conditional){
    true_effects_samples_conditional <- true_effects_samples[,mu_indicator %in% which(!mu_is_null), drop=FALSE]
    true_effects_samples_conditional <- lapply(1:nrow(true_effects_samples_conditional), function(i) {
      .transform_mu(true_effects_samples_conditional[i,], from = model_scale, to = output_scale)
    })
    names(true_effects_samples_conditional) <- sapply(seq_along(true_effects_samples_conditional), function(x) paste0("theta[", x, "]"))
  }

  # transform the effect sizes (and name the matrix)
  true_effects_samples <- lapply(1:nrow(true_effects_samples), function(i) {
    .transform_mu(true_effects_samples[i,], from = model_scale, to = output_scale)
  })
  names(true_effects_samples) <- sapply(seq_along(true_effects_samples), function(x) paste0("theta[", x, "]"))

  # return samples if requested
  if (as_samples){
    if(conditional){
      return(do.call(cbind, true_effects_samples_conditional))
    }else{
      return(do.call(cbind, true_effects_samples))
    }
  }

  # obtain estimates tables
  estimates <- BayesTools::ensemble_estimates_table(
    samples    = true_effects_samples,
    parameters = names(true_effects_samples),
    probs      = probs,
    title      = "True effect estimates:",
    footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale))
  )

  if(conditional){
    estimates_conditional <- BayesTools::ensemble_estimates_table(
      samples    = true_effects_samples_conditional,
      parameters = names(true_effects_samples_conditional),
      probs      = probs,
      title      = "Conditional true effect estimates:",
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
