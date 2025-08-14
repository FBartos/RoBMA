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

# Internal function for re-evaluating model convergence with updated criteria:
# Purpose: Allows post-hoc changes to convergence diagnostics without refitting
# 
# Process:
# 1. Extracts MCMC fit and marginal likelihood from model object
# 2. Skips evaluation for error/null models (automatically non-converged)
# 3. Re-runs JAGS convergence diagnostics with new thresholds
# 4. For bridge sampling: validates marginal likelihood computation
# 5. Updates model convergence status and error/warning messages
#
# Key behaviors:
# - Bridge sampling requires both MCMC convergence AND valid marginal likelihood
# - Spike-and-slab (ss) only requires MCMC convergence 
# - Failed models can be kept or removed based on remove_failed setting
# - Aggregates warnings from MCMC fitting, convergence checks, and marginal likelihood
#
# This enables users to adjust convergence criteria after model fitting
# without requiring computationally expensive re-fitting
.update_model_checks       <- function(model, convergence_checks, algorithm = "bridge"){

  fit     <- model[["fit"]]
  marglik <- model[["marglik"]]

  # Skip models that failed to fit or are null models (automatically non-converged)
  if(inherits(fit, "error") || inherits(fit, "null_model")){
    return(model)
  }

  errors    <- NULL
  warnings  <- NULL

  # Re-evaluate MCMC convergence diagnostics with updated criteria
  # Uses BayesTools to check Rhat, ESS, and Monte Carlo error thresholds
  check_fit     <- BayesTools::JAGS_check_convergence(
    fit          = fit,
    prior_list   = attr(fit, "prior_list"),
    max_Rhat     = convergence_checks[["max_Rhat"]],      # Gelman-Rubin statistic threshold
    min_ESS      = convergence_checks[["min_ESS"]],       # Effective sample size minimum
    max_error    = convergence_checks[["max_error"]],     # Monte Carlo standard error maximum
    max_SD_error = convergence_checks[["max_SD_error"]]   # SD of MC error maximum
  )
  
  # Aggregate warnings from fit and convergence check
  warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
  
  # Determine convergence based on MCMC diagnostics and user preferences
  if(convergence_checks[["remove_failed"]] && !check_fit){
    converged <- FALSE
  }else{
    converged <- TRUE
  }

  # Bridge sampling algorithm requires additional marginal likelihood validation
  # Spike-and-slab (ss) algorithm doesn't use marginal likelihoods
  if(algorithm == "bridge"){
    # Check that marginal likelihood computation succeeded
    if(is.na(marglik$logml)){
      errors         <- c(errors, attr(marglik, "errors"))
      converged      <- FALSE  # Failed marginal likelihood = non-converged model
    }else{
      # Forward any marginal likelihood warnings (but don't fail the model)
      warnings <- c(warnings, attr(marglik, "warnings"))
    }
  }

  # Update model object with new convergence assessment
  model$errors    <- errors
  model$warnings  <- warnings
  model$converged <- converged

  return(model)
}
# Internal function to determine if a model has only constant (point) priors:
# Purpose: Identifies models that don't require MCMC sampling because all parameters are fixed
#
# Logic:
# - Simple models: all priors must be point/none AND no publication bias (omega) priors
# - Regression models: all non-term priors must be point/none AND all term priors must be point 
#   AND no publication bias priors
#
# Key concepts:
# - Point priors: fix parameter to specific value (no uncertainty)
# - None priors: parameter not included in model
# - If all parameters are fixed, no MCMC sampling needed (analytical solution)
# - Publication bias (omega) priors require sampling even if other priors are point
#
# Returns: TRUE if model is constant (no sampling needed), FALSE otherwise
.is_model_constant         <- function(priors){
  
  if(is.null(priors[["terms"]])){
    # Simple (non-regression) models: check all components except terms
    # Must have no omega priors (publication bias requires sampling)
    return(all(sapply(priors, function(prior) is.prior.point(prior) | is.prior.none(prior))) && is.null(priors[["omega"]]))
  }else{
    # Regression models: separate checks for non-term and term priors
    # Non-term priors (effect, heterogeneity, etc.): can be point or none, no omega allowed
    non_terms <- all(sapply(priors[names(priors) != "terms"], function(prior) is.prior.point(prior) | is.prior.none(prior))) && is.null(priors[["omega"]])
    
    # Term priors (regression coefficients): must all be point priors (not none)
    # Regression coefficients can't be "none" - they must be either estimated or fixed
    terms     <- all(sapply(priors[["terms"]], function(prior) is.prior.point(prior)))
    
    # Both conditions must be satisfied for model to be constant
    return(non_terms && terms)
  }
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
# Internal function to aggregate and format errors/warnings for user display:
# Purpose: Collects errors, warnings, and convergence issues across all models
# 
# Process:
# 1. Shortens long warning/error lists to manageable size (max_print limit)
# 2. Adds convergence failure summary if any models failed
# 3. Returns formatted message vector for user feedback
#
# Design rationale:
# - Users shouldn't be overwhelmed with hundreds of similar warnings
# - Convergence failures are critical and always reported
# - Full details available via check_RoBMA() for debugging
.collect_errors_and_warnings <- function(object, max_print = 5){

  short_warnings <- .shorten_warnings(object$add_info[["warnings"]], max_print)
  short_errors   <- .shorten_errors(object$add_info[["errors"]],     max_print)
  conv_warning   <- .convergence_warning(object)

  return(c(short_warnings, short_errors, conv_warning))
}

# Internal function to extract convergence status across different algorithms:
# Purpose: Returns logical vector indicating which models converged successfully
#
# Algorithm differences:
# - Bridge sampling: multiple models stored in object[["models"]], each with convergence status  
# - Spike-and-slab: single model stored in object[["model"]] with single convergence status
#
# Returns: logical vector (bridge) or single logical value (ss) indicating convergence
.get_model_convergence       <- function(object){
  if(object$add_info[["algorithm"]] == "bridge"){
    # Bridge sampling: check convergence status for each model in the ensemble
    # Returns logical vector with length = number of models
    return(sapply(object[["models"]], function(model) if(is.null(model[["converged"]])) FALSE else model[["converged"]]))
  }else if(object$add_info[["algorithm"]] == "ss"){
    # Spike-and-slab: single model convergence status
    # Returns single logical value
    return(if(is.null(object[["model"]][["converged"]])) FALSE else object[["model"]][["converged"]])
  }
}

# Internal function to extract warnings from models with algorithm-specific formatting:
# Purpose: Collects all warning messages across fitted models for user feedback
# 
# Formatting differences:
# - Bridge sampling: prefixes warnings with "Model (i):" to identify source model
# - Spike-and-slab: returns warnings directly (only one model)
.get_model_warnings          <- function(object){
  if(object$add_info[["algorithm"]] == "bridge"){
    # Bridge sampling: format warnings with model index for identification
    return(unlist(sapply(seq_along(object[["models"]]), function(i){
      if(!is.null(object[["models"]][[i]][["warnings"]])){
        paste0("Model (", i, "): ", object[["models"]][[i]][["warnings"]])
      }
    })))
  }else if(object$add_info[["algorithm"]] == "ss"){
    # Spike-and-slab: direct warning passthrough
    return(object[["model"]][["warnings"]])
  }
}

# Internal function to extract errors from models with algorithm-specific formatting:
# Purpose: Collects all error messages across fitted models for debugging
# Same formatting logic as warnings but for errors
.get_model_errors            <- function(object){
  if(object$add_info[["algorithm"]] == "bridge"){
    return(unlist(sapply(seq_along(object[["models"]]), function(i){
      if(!is.null(object[["models"]][[i]][["errors"]])){
        paste0("Model (", i, "): ", object[["models"]][[i]][["errors"]])
      }
    })))
  }else if(object$add_info[["algorithm"]] == "ss"){
    return(object[["model"]][["errors"]])
  }
}
.note_omega                  <- function(object){

  bias_distributions <- sapply(object[["priors"]][["bias"]], function(prior) prior[["distribution"]])

  if(length(bias_distributions) == 0){
    return(NULL)
  }

  if(object[["add_info"]][["algorithm"]] == "ss"){
    return("(Estimated publication weights omega correspond to one-sided p-values.)")
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
# Internal function to determine if a model component is specified as "null":
# Purpose: Checks if effect, heterogeneity, bias, or hierarchical components are absent
#
# Component-specific logic:
# - Effect: checks mu (simple models) or intercept (regression models)
# - Heterogeneity: checks tau parameter 
# - Bias: checks omega (publication bias) parameter
# - Hierarchical: checks rho (between-study correlation) parameter
#
# Key concepts:
# - NULL prior: component not included in model at all
# - is_null flag: component included but set to test null hypothesis (typically point prior at 0)
# - Simple vs regression models have different parameter structures
#
# Returns: TRUE if component represents null hypothesis, FALSE if alternative hypothesis
.is_component_null           <- function(priors, component){
  if(component == "effect"){
    if(is.null(priors[["terms"]])){
      # Simple models: check mu parameter directly
      if(is.null(priors[["mu"]])){
        return(TRUE)  # No effect parameter = null effect
      }else{
        return(priors[["mu"]][["is_null"]])  # Check if effect is set to null (e.g., point prior at 0)
      }
    }else{
      # Regression models: effect is represented by intercept term
      return(priors[["terms"]][["intercept"]][["is_null"]])
    }
  }else if(component == "heterogeneity"){
    # Heterogeneity: tau parameter controls between-study variance
    if(is.null(priors[["tau"]])){
      return(TRUE)  # No tau = homogeneous effects (null heterogeneity)
    }else{
      return(priors[["tau"]][["is_null"]])  # Check if tau set to null (typically point prior at 0)
    }
  }else if(component == "baseline"){
    # Baseline: pi parameter for baseline event rates in BiBMA models
    if(is.null(priors[["pi"]])){
      return(TRUE)  # No baseline parameter specified
    }else{
      return(priors[["pi"]][["is_null"]])  # Check if baseline is null
    }
  }else if(component == "hierarchical"){
    # Hierarchical: rho parameter for between-study correlations in multilevel models
    # Special logic: hierarchical models require both tau > 0 AND rho parameter
    if((is.prior.point(priors[["tau"]]) && priors[["tau"]]$parameters[["location"]] == 0) || is.null(priors[["rho"]])){
      return(TRUE)  # No hierarchical structure if tau=0 or no rho parameter
    }else{
      return(priors[["rho"]][["is_null"]])  # Check if rho set to null
    }
  }else if(component == "bias"){
    # Publication bias: can be modeled via omega (weight functions), PET, or PEESE
    # Check for any bias parameter in order of preference
    if(!is.null(priors[["omega"]])){
      return(priors[["omega"]][["is_null"]])  # Weight functions
    }else if(!is.null(priors[["PET"]])){
      return(priors[["PET"]][["is_null"]])    # PET models
    }else if(!is.null(priors[["PEESE"]])){
      return(priors[["PEESE"]][["is_null"]])  # PEESE models
    }else{
      return(TRUE)  # No bias parameters = null bias
    }
  }
}
.is.prior_null               <- function(prior){
  if(is.null(prior)){
    return()
  }else if(is.prior.mixture(prior)){
    attr(prior, "components") == "null"
  }else if(is.prior.spike_and_slab(prior)){
    return(c(TRUE, FALSE))
  }else{
    return(prior[["is_null"]])
  }
}
.multivariate_warning        <- function(){
  warning("You are about to estimate multivariate models. Note that this is an extremely computationaly expensive experimental feature.", immediate. = TRUE, call. = FALSE)
}
.weighted_warning            <- function(){
  warning("You are about to estimate weighted models. Note that this is an experimental feature.", immediate. = TRUE, call. = FALSE)
}
.update_object               <- function(object){

  .check_is_any_RoBMA_object(object)

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

  # 3.2 -> 3.3
  if(.object_version_older(object, "3.3.0")){
    object[["add_info"]][["algorithm"]] <- "bridge"
    object[["add_info"]][["version"]]   <- list(c(2,4,0))
  }

  return(object)
}
.object_version_older        <- function(object, version){

  object  <- unlist(object[["add_info"]][["version"]])
  current <- as.numeric(unlist(strsplit(version, ".", fixed = TRUE)))

  # deal with potential missing trailing zeroes
  if(length(object) == 2){
    object <- c(object, 0)
  }
  if(length(current) == 2){
    current <- c(current, 0)
  }

  if(length(object) != 3)
    stop("The current version number is not in the correct format.")
  if(length(current) != 3)
    stop("The version number to compare with is not in the correct format.")

  # transform by 10^3 to compare the versions
  object  <- object[1] * 1000^2 + object[2] * 1000 + object[3]
  current <- current[1] * 1000^2 + current[2] * 1000 + current[3]

  return(object < current)
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
.get_K                       <- function(object){
  if(.is_regression(object)){
    ids <- object$data$outcome[["study_ids"]]
  }else{
    ids <- object$data[["study_ids"]]
  }
  ids_NA <- is.na(ids)
  K      <- length(unique(ids[!ids_NA])) + sum(ids_NA)
  return(K)
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
.BayesTools_parameter_name    <- function(parameter){
  return(BayesTools::JAGS_parameter_names(parameter, formula_parameter = "mu"))
}
.BayesTools_make_column_names <- function(table, column = 1){
  rownames(table) <- table[,column]
  attr(table, "rownames") <- TRUE
  table <- BayesTools::remove_column(table, column)
  return(table)
}
.output_parameter_names       <- function(parameter){
  return(BayesTools::format_parameter_names(parameter, formula_parameters = "mu", formula_prefix = FALSE))
}
# Internal function returning reserved keywords that cannot be used as predictor names:
# Purpose: Prevents naming conflicts between user predictors and internal parameters
#
# Categories of reserved words:
# 1. JAGS parameter names: mu, tau, theta, omega, rho, eta, PET, PEESE, pi, gamma
# 2. Internal component names: intercept, terms, component_*
# 3. Model type identifiers: weightfunction, PET-PEESE, etc.
# 4. Data column names: d, t, r, z, y, logOR, OR, lCI, uCI, v, se, n, weight, x1, x2, n1, n2
# 5. Component identifiers: component_effect, component_heterogeneity, etc.
#
# Using these names as predictors would cause conflicts in JAGS model code generation
# or internal data processing, leading to unpredictable behavior
.reserved_words               <- function() c("intercept", "Intercept", "terms", "mu", "tau", "theta", "omega", "rho", "eta", "PET", "PEESE", "pi", "gamma",
                                              "weightfunction", "weigthfunction", "PET-PEESE", "PETPEESE",
                                              "d", "t", "r", "z", "y", "logOR", "OR", "lCI", "uCI", "v", "se", "n", "weight", "x1", "x2", "n1", "n2",
                                              "component_effect", "component_heterogeneity", "component_bias", "component_hierarchical", "component_baseline")
