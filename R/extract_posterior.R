#' @title Extract Posterior Samples from a RoBMA Model
#'
#' @description Extracts posterior samples for a specified parameter from a RoBMA model object.
#'
#' @param x a fitted RoBMA object
#' @param parameter a parameter for which posterior samples
#' should be extracted. Defaults to \code{"mu"} (for the effect size).
#' The additional options  are \code{"tau"} (for the heterogeneity),
#' \code{"weightfunction"} (for the estimated weightfunction),
#' or \code{"PET"} and \code{"PEESE"} (for the PET-PEESE coefficients).
#' @param conditional whether conditional estimates should be
#' extracted. Defaults to \code{FALSE} which extracts the model-averaged
#' estimates. Note that both \code{"weightfunction"} and
#' \code{"PET-PEESE"} are always ignoring the other type of
#' publication bias adjustment.
#' @param output_scale transform the effect sizes and the meta-analytic
#' effect size estimate to a different scale. Defaults to \code{NULL}
#' which returns the same scale as the model was estimated on.
#' @param ... additional arguments passed to the method.
#'
#' @return A matrix containing the posterior samples for the specified parameter.
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a fitted RoBMA model:
#' posterior_mu <- extract_posterior(fit, parameter = "mu")
#' }
#'
#' @export
extract_posterior <- function(x, parameter = "mu", conditional = FALSE, output_scale = NULL, ...) {

  dots <- list(...)

  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  # all remaining input is checked in the .extract_posterior.RoBMA function

  ### transform the samples
  # check the scales
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else if(x$add_info[["output_scale"]] == "y" & .transformation_var(output_scale, estimation = FALSE) != "y"){
    stop("Models estimated using the generall effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale, estimation = FALSE)
  }

  # transform the estimates if needed
  if(x$add_info[["output_scale"]] != output_scale){
    x <- .transform_posterior(x, x$add_info[["output_scale"]], output_scale)
  }

  samples <- .extract_posterior.RoBMA(x, parameter, conditional)

  # get the name of the parameter
  parameter_samples <- .check_and_transform_parameter_name(parameter, x$add_info[["predictors"]])
  if (parameter_samples == "PETPEESE") {
    samples <- samples[names(samples)[names(samples) %in% c("mu", "mu_intercept", "PET", "PEESE")]]
    samples <- do.call(cbind, samples)
  } else {
    samples <- samples[[parameter_samples]]
    attributes(samples)[!names(attributes(samples)) %in% c("dim", "dimnames")] <- NULL
  }

  # return metadata if requested
  if (!is.null(dots[["metadata"]]) && isTRUE(dots[["metadata"]]))
    return(samples)

  # transform to a matrix
  if (is.null(dim(samples))) {
    samples <- matrix(samples, ncol = 1)
    colnames(samples) <- parameter
  } else if (is.null(colnames(samples))) {
    parameter <- as.matrix(parameter)
  }

  return(samples)
}

.extract_posterior.RoBMA <- function(x, parameter, conditional) {

  # apply version changes to RoBMA object
  x <- .update_object(x)

  # check whether plotting is possible
  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  # check settings
  BayesTools::check_char(parameter, "parameter")
  BayesTools::check_bool(conditional, "conditional")

  # get the name of the parameter
  parameter         <- .check_and_transform_parameter_name(parameter, x$add_info[["predictors"]])
  parameter_samples <- .check_and_transform_parameter_samples_name(parameter, x$add_info[["predictors"]])

  # choose the samples
  if(conditional && parameter == "PETPEESE"){

    # get model-averaged posterior across PET and PEESE parameters
    models <- x[["models"]]

    effect <- sapply(x[["models"]], function(model)!.is_component_null(model[["priors"]], "effect"))
    PET    <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
    PEESE  <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

    if(any(PET) & any(PEESE)){
      is_conditional  <- effect & (PET | PEESE)
      if(sum(is_conditional) == 0)
        stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of PET-PEESE publication bias adjustment. Please, verify that you specified at least one model assuming the presence of PET-PEESE publication bias adjustment.")
      parameters      <- c("mu", "PET", "PEESE")
      parameters_null <- c("mu" = list(!is_conditional), "PET" = list(!is_conditional), "PEESE" = list(!is_conditional))
    }else if(any(PET)){
      is_conditional  <- effect & PET
      if(sum(is_conditional) == 0)
        stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of PET-PEESE publication bias adjustment. Please, verify that you specified at least one model assuming the presence of PET-PEESE publication bias adjustment.")
      parameters      <- c("mu", "PET")
      parameters_null <- c("mu" = list(!is_conditional), "PET"   = list(!is_conditional))
    }else if(any(PEESE)){
      is_conditional  <- effect & PEESE
      if(sum(is_conditional) == 0)
        stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of PET-PEESE publication bias adjustment. Please, verify that you specified at least one model assuming the presence of PET-PEESE publication bias adjustment.")
      parameters      <- c("mu", "PEESE")
      parameters_null <- c("mu" = list(!is_conditional), "PEESE" = list(!is_conditional))
    }else{
      stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of PET-PEESE publication bias adjustment. Please, verify that you specified at least one model assuming the presence of PET-PEESE publication bias adjustment.")
    }

    samples <- BayesTools::mix_posteriors(
      model_list   = models,
      parameters   = parameters,
      is_null_list = parameters_null,
      seed         = x$add_info[["seed"]],
      conditional  = TRUE
    )

  }else{
    if(conditional & parameter %in% c("mu", "tau", "rho", "PET", "PEESE", "PETPEESE", "omega")){
      samples <- x[["RoBMA"]][["posteriors_conditional"]]
    }else if(conditional & parameter %in% x$add_info[["predictors"]]){
      samples <- x[["RoBMA"]][["posteriors_predictors_conditional"]]
    }else if(!conditional &  parameter %in% c("mu", "tau", "rho", "PET", "PEESE", "PETPEESE", "omega")){
      samples <- x[["RoBMA"]][["posteriors"]]
    }else if(!conditional &  parameter %in% x$add_info[["predictors"]]){
      samples <- x[["RoBMA"]][["posteriors_predictors"]]
    }
  }


  if(parameter %in% c("mu", "tau", "rho", "omega", "PET", "PEESE")){
    if(conditional && is.null(samples[[parameter]])){
      switch(
        parameter,
        "mu"    = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect."),
        "tau"   = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the heterogeneity. Please, verify that you specified at least one model assuming the presence of the heterogeneity."),
        "PET"   = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the PET models. Please, verify that you specified at least one model assuming the presence of the PET models."),
        "PEESE" = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the PEESE models. Please, verify that you specified at least one model assuming the presence of the PEESE models."),
        "rho"   = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the three-level structure. Please, verify that you specified at least one model assuming the presence of the three-level structure."),
        "omega" = stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of selection models publication bias adjustment. Please, verify that you specified at least one model assuming the presence of selection models publication bias adjustment.")
      )
    }else if(is.null(samples[[parameter]])){
      switch(
        parameter,
        "mu"    = stop("The ensemble does not contain any posterior samples model-averaged across the effect. Please, verify that you specified at least one model for the effect."),
        "tau"   = stop("The ensemble does not contain any posterior samples model-averaged across the heterogeneity. Please, verify that you specified at least one model for the heterogeneity."),
        "PET"   = stop("The ensemble does not contain any posterior samples model-averaged across the PET. Please, verify that you specified at least one model for the PET."),
        "PEESE" = stop("The ensemble does not contain any posterior samples model-averaged across the PEESE. Please, verify that you specified at least one model for the PEESE."),
        "rho"   = stop("The ensemble does not contain any posterior samples model-averaged across the three-level structure. Please, verify that you specified at least one model for the three-level structure."),
        "omega" = stop("The ensemble does not contain any posterior samples model-averaged across the selection models publication bias adjustment. Please, verify that you specified at least one selection models publication bias adjustment.")
      )
    }
  }else if(parameter %in% "PETPEESE"){
    # checking for mu since it's the common parameter for PET-PEESE
    if(conditional && (is.null(samples[["mu"]]) || is.null(samples[["PET"]]) && is.null(samples[["PEESE"]]))){
      stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of PET-PEESE publication bias adjustment. Please, verify that you specified at least one model assuming the presence of PET-PEESE publication bias adjustment.")
    }else if(is.null(samples[["mu"]]) || is.null(samples[["PET"]]) && is.null(samples[["PEESE"]])){
      stop("The ensemble does not contain any posterior samples model-averaged across the PET-PEESE publication bias adjustment. Please, verify that you specified at least one PET-PEESE publication bias adjustment.")
    }
  }else if(parameter %in% x$add_info[["predictors"]]){
    if(conditional && is.null(samples[[parameter_samples]])){
      stop(sprintf("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the '%1$s' predictor. Please, verify that you specified at least one model assuming the presence of '%1$s' predictor.", parameter))
    }else if(is.null(samples[[parameter_samples]])){
      stop(sprintf("The ensemble does not contain any posterior samples model-averaged across the '%1$s' predictor. Please, verify that you specified at least one model containing the '%1$s' predictor.", parameter))
    }
  }

  return(samples)
}

.check_and_transform_parameter_name         <- function(parameter, predictors) {
  # deal with bad parameter names for PET-PEESE, weightfunction
  if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
    parameter         <- "omega"
  }else if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) == "petpeese"){
    parameter         <- "PETPEESE"
  }else if(parameter %in% c("mu", "tau", "rho", "PET", "PEESE")){
    parameter         <- parameter
  }else if(length(predictors) > 0 && parameter %in% predictors){
    parameter         <- parameter
  }else{
    if(length(predictors) > 0){
      stop(paste0("The passed parameter does not correspond to any of main model parameter ('mu', 'tau', 'omega', 'PET', 'PEESE') or any of the specified predictors: ", paste0("'", predictors, "'", collapse = ", "), ". See '?plot.RoBMA' for more details."))
    }else{
      stop(paste0("The passed parameter does not correspond to any of main model parameter ('mu', 'tau', 'omega', 'PET', 'PEESE'). See '?plot.RoBMA' for more details."))
    }
  }

  return(parameter)
}
.check_and_transform_parameter_samples_name <- function(parameter, predictors) {
  # deal with bad parameter names for PET-PEESE, weightfunction
  if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
    parameter_samples <- "omega"
  }else if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) == "petpeese"){
    parameter_samples <- "PETPEESE"
  }else if(parameter %in% c("mu", "tau", "rho", "PET", "PEESE")){
    parameter_samples <- parameter
  }else if(length(predictors) > 0 && parameter %in% predictors){
    parameter_samples <- .BayesTools_parameter_name(parameter)
  }else{
    if(length(predictors) > 0){
      stop(paste0("The passed parameter does not correspond to any of main model parameter ('mu', 'tau', 'omega', 'PET', 'PEESE') or any of the specified predictors: ", paste0("'", predictors, "'", collapse = ", "), ". See '?plot.RoBMA' for more details."))
    }else{
      stop(paste0("The passed parameter does not correspond to any of main model parameter ('mu', 'tau', 'omega', 'PET', 'PEESE'). See '?plot.RoBMA' for more details."))
    }
  }

  return(parameter_samples)
}

