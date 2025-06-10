#' @title Compute pooled effect size
#'
#' @description \code{pooled_effect} computes the pooled effect size
#' for a fitted RoBMA.reg and BiBMA.reg object.
#'
#' @inheritParams summary.RoBMA
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
#' for a fitted RoBMA.reg and BiBMA.reg object.
#'
#' @inheritParams summary.RoBMA
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

.compute_effect <- function(object, conditional = FALSE, output_scale = NULL, probs = c(.025, .975), type = "pooled", ...) {

  .check_is_any_RoBMA_object(object)
  if(.is_model_regression(object))
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

  dots <- list(...)

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
  if (!is.null(dots[["as_samples"]]) && isTRUE(dots[["as_samples"]])){
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
