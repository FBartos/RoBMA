#' @title Check fitted RoBMA object for errors and warnings
#'
#' @description Checks fitted RoBMA object
#' for warnings and errors and prints them to the
#' console.
#'
#' @param fit a fitted RoBMA object.
#'
#'
#' @return \code{check_RoBMA} returns a vector of error and
#' warning messages. \code{check_RoBMA_convergence} returns
#' a logical vector indicating whether the models have
#' converged.
#'
#' @name check_RoBMA
#' @aliases check_RoBMA check_RoBMA_convergence
#' @export check_RoBMA
#' @export check_RoBMA_convergence

#' @rdname check_RoBMA
check_RoBMA <- function(fit){
  .check_is_any_RoBMA_object(fit)
  .print_errors_and_warnings(fit, max_print = Inf)
}
#' @rdname check_RoBMA
check_RoBMA_convergence <- function(fit){
  .check_is_any_RoBMA_object(fit)
  return(.get_model_convergence(fit))
}


.update_model_checks       <- function(model, convergence_checks){

  fit     <- model[["fit"]]
  marglik <- model[["marglik"]]

  # skip models with an error or emty models
  if(inherits(fit, "error") || inherits(fit, "null_model")){
    return(model)
  }

  errors    <- NULL
  warnings  <- NULL

  # re-evaluate convergence checks with the new settigns
  check_fit     <- BayesTools::JAGS_check_convergence(
    fit          = fit,
    prior_list   = attr(fit, "prior_list"),
    max_Rhat     = convergence_checks[["max_Rhat"]],
    min_ESS      = convergence_checks[["min_ESS"]],
    max_error    = convergence_checks[["max_error"]],
    max_SD_error = convergence_checks[["max_SD_error"]]
  )
  warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
  if(convergence_checks[["remove_failed"]] && !check_fit){
    converged <- FALSE
  }else{
    converged <- TRUE
  }

  # check the marginal likelihood and keep the errors
  if(is.na(marglik$logml)){

    errors         <- c(errors, attr(marglik, "errors"))
    converged      <- FALSE

  }else{

    # forward warnings if present
    warnings <- c(warnings, attr(marglik, "warnings"))

  }

  model$errors    <- errors
  model$warnings  <- warnings
  model$converged <- converged

  return(model)
}
.is_model_constant         <- function(priors){
  # checks whether there is at least one non-nill prior
  if(is.null(priors[["terms"]])){
    # in simple models
    return(all(sapply(priors, function(prior) is.prior.point(prior) | is.prior.none(prior))) & is.null(priors[["omega"]]))
  }else{
    # in regression models
    non_terms <- all(sapply(priors[names(priors) != "terms"], function(prior) is.prior.point(prior) | is.prior.none(prior))) & is.null(priors[["omega"]])
    terms     <- all(sapply(priors[["terms"]], function(prior) is.prior.point(prior)))
    return(non_terms && terms)
  }
  return()
}
.remove_model_posteriors   <- function(object){
  for(i in seq_along(object[["models"]])){
    if(inherits(object$models[[i]][["fit"]], "runjags")){
      object$models[[i]]$fit[["mcmc"]] <- NULL
    }
  }
  return(object)
}
.remove_model_margliks     <- function(object){
  for(i in seq_along(object[["models"]])){
    if(inherits(object$models[[i]][["marglik"]], "bridge")){
      object$models[[i]]$marglik[["q11"]] <- NULL
      object$models[[i]]$marglik[["q12"]] <- NULL
      object$models[[i]]$marglik[["q21"]] <- NULL
      object$models[[i]]$marglik[["q22"]] <- NULL
    }
  }
  return(object)
}
.print_errors_and_warnings <- function(object, max_print = 5){

  errors_and_warnings <- .collect_errors_and_warnings(object, max_print = max_print)

  for(i in seq_along(errors_and_warnings))
    warning(errors_and_warnings[i], immediate. = TRUE, call. = FALSE)

  return(invisible(errors_and_warnings))
}
.shorten_warnings    <- function(warnings, n_warnings = 5){
  if(is.null(warnings)){
    return(NULL)
  }else if(length(warnings) <= n_warnings){
    return(warnings)
  }else{
    return(c(warnings[1:n_warnings], paste0("There were another ", length(warnings) - n_warnings - 1, " warnings. To see all warnings call 'check_RoBMA(fit)'.")))
  }
}
.shorten_errors      <- function(errors, n_errors = 5){
  if(is.null(errors)){
    return(NULL)
  }else if(length(errors) <= n_errors){
    return(errors)
  }else{
    return(c(errors[1:n_errors], paste0("There were another ", length(errors) - n_errors - 1, " errors. To see all errors call 'check_RoBMA(fit)'.")))
  }
}
.convergence_warning <- function(object){
  if(any(!.get_model_convergence(object))){
    return(paste0(sum(!.get_model_convergence(object)), ifelse(sum(!.get_model_convergence(object)) == 1, " model", " models"), " failed to converge."))
  }else{
    return(NULL)
  }
}
.collect_errors_and_warnings <- function(object, max_print = 5){

  short_warnings <- .shorten_warnings(object$add_info[["warnings"]], max_print)
  short_errors   <- .shorten_errors(object$add_info[["errors"]],     max_print)
  conv_warning   <- .convergence_warning(object)

  return(c(short_warnings, short_errors, conv_warning))
}
.get_model_convergence       <- function(object){
  return(sapply(object[["models"]], function(model) if(is.null(model[["converged"]])) FALSE else model[["converged"]]))
}
.get_model_warnings          <- function(object){
  return(unlist(sapply(seq_along(object[["models"]]), function(i){
    if(!is.null(object[["models"]][[i]][["warnings"]])){
      paste0("Model (", i, "): ", object[["models"]][[i]][["warnings"]])
    }
  })))
}
.get_model_errors            <- function(object){
  return(unlist(sapply(seq_along(object[["models"]]), function(i){
    if(!is.null(object[["models"]][[i]][["errors"]])){
      paste0("Model (", i, "): ", object[["models"]][[i]][["errors"]])
    }
  })))
}
.note_omega                  <- function(object){

  bias_distributions <- sapply(object[["priors"]][["bias"]], function(prior) prior[["distribution"]])

  if(length(bias_distributions) == 0){
    return(NULL)
  }

  is_one.sided       <- grepl("one.sided", bias_distributions)
  is_two.sided       <- grepl("two.sided", bias_distributions)

  if(any(is_one.sided) & any(is_two.sided)){
    return("(Estimated publication weights omega correspond to one-sided p-values.)")
  }else if(any(is_one.sided)){
    return("(Estimated publication weights omega correspond to one-sided p-values.)")
  }else if(any(is_two.sided)){
    return("(Estimated publication weights omega correspond to two-sided p-values.)")
  }else{
    return(NULL)
  }

}
.is_component_null           <- function(priors, component){
  if(component == "effect"){
    if(is.null(priors[["terms"]])){
      if(is.null(priors[["mu"]])){
        return(TRUE)
      }else{
        return(priors[["mu"]][["is_null"]])
      }
    }else{
      return(priors[["terms"]][["intercept"]][["is_null"]])
    }
  }else if(component == "heterogeneity"){
    if(is.null(priors[["tau"]])){
      return(TRUE)
    }else{
      return(priors[["tau"]][["is_null"]])
    }
  }else if(component == "baseline"){
    if(is.null(priors[["pi"]])){
      return(TRUE)
    }else{
      return(priors[["pi"]][["is_null"]])
    }
  }else if(component == "hierarchical"){
    if((is.prior.point(priors[["tau"]]) && priors[["tau"]]$parameters[["location"]] == 0) || is.null(priors[["rho"]])){
      return(TRUE)
    }else{
      return(priors[["rho"]][["is_null"]])
    }
  }else if(component == "bias"){
    if(!is.null(priors[["omega"]])){
      return(priors[["omega"]][["is_null"]])
    }else if(!is.null(priors[["PET"]])){
      return(priors[["PET"]][["is_null"]])
    }else if(!is.null(priors[["PEESE"]])){
      return(priors[["PEESE"]][["is_null"]])
    }else{
      return(TRUE)
    }
  }
}
.multivariate_warning        <- function(){
  warning("You are about to estimate multivariate models. Note that this is an extremely computationaly expensive experimental feature.", immediate. = TRUE, call. = FALSE)
}
.weighted_warning            <- function(){
  warning("You are about to estimate weighted models. Note that this is an experimental feature.", immediate. = TRUE, call. = FALSE)
}
.update_object               <- function(object){

  # no package version number saved prior to 2.4
  if(!all("version" %in% names(object[["add_info"]]))){

    # 2.1 -> 2.2
    if(is.null(object[["formula"]]) && is.null(attr(object$data, "all_independent"))){
      attr(object$data, "all_independent") <- TRUE
    }

    object[["add_info"]][["version"]] <- list(c(2,2,0))
  }

  # 2.2 -> 2.3
  if(.object_version_older(object, "2.2.0")){

    if(is.null(object[["formula"]]) && !is.null(attr(object$data, "all_independent")) && is.null(attr(object$data, "weighted"))){
      attr(object$data, "weighted") <- FALSE
    }

    object[["add_info"]][["version"]] <- list(c(2,3,0))
  }

  # 2.3 -> 2.4
  if(.object_version_older(object, "2.3.0")){

    if(!all(c("predictors", "predictors_test", "standardize_predictors") %in% names(object[["add_info"]]))){
      object[["add_info"]][["predictors"]]             <- NULL
      object[["add_info"]][["predictors_test"]]        <- NULL
      object[["add_info"]][["standardize_predictors"]] <- NULL
      object[["add_info"]] <- object[["add_info"]][c(
        "model_type",
        "predictors",
        "predictors_test",
        "prior_scale",
        "output_scale",
        "effect_measure",
        "effect_direction",
        "standardize_predictors",
        "seed",
        "save",
        "warnings",
        "errors"
      )]
    }
    if(!all("version" %in% names(object[["add_info"]]))){
      object[["add_info"]][["version"]] <- utils::packageVersion("RoBMA")
    }
    if(!all(c("inference_predictors", "inference_predictors_conditional", "posteriors_predictors", "posteriors_predictors_conditional") %in% names(object[["RoBMA"]]))){
      object[["RoBMA"]] <- list(
        "inference"                         = object[["RoBMA"]][["inference"]],
        "inference_conditional"             = object[["RoBMA"]][["inference_conditional"]],
        "inference_predictors"              = NULL,
        "inference_predictors_conditional"  = NULL,
        "posteriors"                        = object[["RoBMA"]][["posteriors"]],
        "posteriors_conditional"            = object[["RoBMA"]][["posteriors_conditional"]],
        "posteriors_predictors"             = NULL,
        "posteriors_predictors_conditional" = NULL
      )
    }
    names(object[["priors"]])[names(object[["priors"]]) == "rho"] <- "hierarchical"
    object[["add_info"]][["version"]] <- list(c(2,4,0))
  }

  return(object)
}
.object_version_older        <- function(object, version){

  object  <- unlist(object[["add_info"]][["version"]])
  current <- as.numeric(unlist(strsplit(version, ".", fixed = TRUE)))

  if(length(object) < 3 | length(current) < 3){
    return(object[1] <= current[1] && object[2] <= current[2])
  }else{
    return(object[1] <= current[1] && object[2] <= current[2] && object[3] <= current[3])
  }
}
.check_is_any_RoBMA_object   <- function(x){
  if(!(is.RoBMA(x) || is.RoBMA.reg(x) || is.NoBMA(x) || is.NoBMA.reg(x) || is.BiBMA(x))){
    stop("The object is not a model fitted with the RoBMA package.", call. = FALSE)
  }
}
.get_one_sided_cuts          <- function(prior){
  steps <- prior[["parameters"]][["steps"]]
  if (grepl("two.sided", prior[["distribution"]])) {
    return(c(1-rev(steps/2), steps/2))
  } else {
    return(steps)
  }
}



# object type checks
.is_multivariate <- function(object){
  if(.is_regression(object)){
    return(!attr(object[["data"]][["outcome"]], "all_independent"))
  }else{
    return(!attr(object[["data"]], "all_independent"))
  }
}
.is_weighted     <- function(object){
  if(.is_regression(object)){
    return(attr(object[["data"]][["outcome"]], "weighted"))
  }else{
    return(attr(object[["data"]], "weighted"))
  }
}
.is_regression   <- function(object){
  return(!is.null(object[["formula"]]))
}

# parameter naming functions
.BayesTools_parameter_name   <- function(parameter){
  return(BayesTools::JAGS_parameter_names(parameter, formula_parameter = "mu"))
}
.output_parameter_names      <- function(parameter){
  return(BayesTools::format_parameter_names(parameter, formula_parameters = "mu", formula_prefix = FALSE))
}
.reserved_words              <- function() c("intercept", "Intercept", "terms", "mu", "tau", "theta", "omega", "rho", "eta", "PET", "PEESE", "pi", "gamma",
                                             "weightfunction", "weigthfunction", "PET-PEESE", "PETPEESE",
                                             "d", "t", "r", "z", "y", "logOR", "OR", "lCI", "uCI", "v", "se", "n", "weight", "x1", "x2", "n1", "n2")
