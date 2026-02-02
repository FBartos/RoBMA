#' @title Options for the RoBMA package
#'
#' @description A placeholder object and functions for the RoBMA package.
#' (adapted from the runjags R package).
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below.
#'
#' @details
#' The available options are:
#' \describe{
#'   \item{\code{max_cores}}{number of cores to use for parallel computing (default to all available cores - 1)}
#'   \item{\code{check_scaling}}{whether to check scaling of predictors (default \code{TRUE})}
#'   \item{\code{silent}}{whether to suppress output (default \code{FALSE})}
#'   \item{\code{autocompute.loo}}{whether to automatically compute LOO (default \code{FALSE})}
#'   \item{\code{autocompute.waic}}{whether to automatically compute WAIC (default \code{FALSE})}
#'   \item{\code{autocompute.marglik}}{whether to automatically compute marginal likelihood (default \code{FALSE})}
#'   \item{\code{default_UISD.effect}}{default scaling of the unit information standard deviation for the effect size parameter (default \code{0.5})}
#'   \item{\code{default_UISD.heterogeneity}}{default scaling of the unit information standard deviation for the heterogeneity parameter (default \code{0.25})}
#'   \item{\code{default_UISD.mods}}{default scaling of the unit information standard deviation for the moderators (default \code{0.25})}
#'   \item{\code{default_UISD.scale}}{default scaling of the unit information standard deviation for the scale parameter (default \code{0.5})}
#'   \item{\code{default_informed_priors.mods}}{default scaling of informed priors for moderators (default \code{0.5})}
#'   \item{\code{default_informed_priors.scale}}{default scaling of informed priors for the scale parameter (default \code{0.5})}
#'   \item{\code{default_bias_weightfunction.alpha}}{default alpha for the weightfunction (default \code{1})}
#'   \item{\code{default_bias_PET.scale}}{default scale for the PET (default \code{1})}
#'   \item{\code{default_bias_PEESE.scale}}{default scale for the PEESE (default \code{5})}
#' }
#'
#' @return The current value of all available RoBMA options (after applying any
#' changes specified) is returned invisibly as a named list.
#'
#' @export RoBMA.options
#' @export RoBMA.get_option
#' @name RoBMA_options
#' @aliases RoBMA_options RoBMA.options RoBMA.get_option
NULL


#' @rdname RoBMA_options
RoBMA.options    <- function(...){

	opts <- list(...)

	for(i in seq_along(opts)){

	  if(!names(opts)[i] %in% names(RoBMA.private))
	    stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))

	  assign(names(opts)[i], opts[[i]] , envir = RoBMA.private)
	}

	return(invisible(RoBMA.private$options))
}

#' @rdname RoBMA_options
RoBMA.get_option <- function(name){

	if(length(name)!=1)
	  stop("Only 1 option can be retrieved at a time")

	if(!name %in% names(RoBMA.private))
	  stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))

	# Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
	return(eval(RoBMA.private[[name]]))
}

# export the function directly to suppress import warnings
.runjags__findjags <- function() runjags::findjags()

# adapted from the runjags package version 2.2.0
RoBMA.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("defaultoptions",  list(
  jagspath = expression(.runjags__findjags()),
  envir    = RoBMA.private))

assign("options",         RoBMA.private$defaultoptions,   envir = RoBMA.private)
assign("RoBMA_version",   "notset",                       envir = RoBMA.private)
assign("min_jags_major",  4,                              envir = RoBMA.private)
assign("max_jags_major",  4,                              envir = RoBMA.private)
assign("max_cores",       parallel::detectCores(logical = TRUE) - 1,  envir = RoBMA.private)
assign("check_scaling",   TRUE,                                       envir = RoBMA.private) # TODO: remove?

assign("silent", FALSE,  envir = RoBMA.private)

assign("autocompute.loo",     FALSE, envir = RoBMA.private)
assign("autocompute.waic",    FALSE, envir = RoBMA.private)
assign("autocompute.marglik", FALSE, envir = RoBMA.private)


### default scaling of unit information prior for different parameters
assign("default_UISD.effect",          1/2,  envir = RoBMA.private)
assign("default_UISD.heterogeneity",   1/4,  envir = RoBMA.private)
assign("default_UISD.mods",            1/4,  envir = RoBMA.private)
assign("default_UISD.scale",           1/2,  envir = RoBMA.private)

### default scaling of informed priors for moderators
assign("default_informed_priors.mods",  1/2,  envir = RoBMA.private)
assign("default_informed_priors.scale", 1/2,  envir = RoBMA.private)

### default setting for the publication bias priors
assign("default_bias_weightfunction.alpha",  1,  envir = RoBMA.private)
assign("default_bias_PET.scale",             1,  envir = RoBMA.private)
assign("default_bias_PEESE.scale",           5,  envir = RoBMA.private)


# check and fix number of threads (sometimes bugs out during installation)
.check_max_cores <- function(){
  if(RoBMA.private$max_cores > parallel::detectCores(logical = TRUE) - 1){
    RoBMA.options(max_cores = parallel::detectCores(logical = TRUE) - 1)
  }
}


# TODO: find a better home for those functions
#' @title Control MCMC fitting process
#'
#' @description \code{set_autofit_control} and \code{set_convergence_checks} control settings for the
#' autofit process of the MCMC JAGS sampler and specify termination criteria and convergence checks.
#'
#' @param max_Rhat maximum value of the R-hat diagnostic.
#' Defaults to \code{1.05}.
#' @param min_ESS minimum estimated sample size.
#' Defaults to \code{500}.
#' @param max_error maximum value of the MCMC error.
#' Defaults to \code{NULL}. Be aware that PEESE publication bias
#' adjustment can have estimates on different scale than the rest of
#' the output, resulting in relatively large max MCMC error.
#' @param max_SD_error maximum value of the proportion of MCMC error
#' of the estimated SD of the parameter.
#' Defaults to \code{NULL}.
#' @param max_time list with the time and unit specifying the maximum
#' autofitting process per model. Passed to \link[base]{difftime} function
#' (possible units are \code{"secs", "mins", "hours", "days", "weeks", "years"}).
#' Defaults to \code{list(time = 60, unit = "mins")}.
#' @param sample_extend number of samples to extend the fitting process if
#' the criteria are not satisfied.
#' Defaults to \code{1000}.
#' @param remove_failed whether models not satisfying the convergence checks should
#' be removed from the inference. Defaults to \code{FALSE} - only a warning is raised.
#' @param balance_probability whether prior model probability should be balanced
#' across the combinations of models with the same H0/H1 for effect / heterogeneity / bias
#' in the case of non-convergence. Defaults to \code{TRUE}.
#' @param restarts number of times new initial values should be generated in case a
#' model fails to initialize. Defaults to \code{10}.
#' @param max_extend number of times after which the automatic fitting function is stopped.
#' Defaults to \code{10}.
#'
#' @details
#' \code{set_autofit_control} controls the automatic model fitting process, determining when to
#' stop iterating based on convergence criteria. \code{set_convergence_checks} sets thresholds
#' for determining whether a model has converged adequately.
#'
#' The autofit control manages computational resources by setting maximum time limits and
#' determining how many additional samples to draw if convergence criteria are not met.
#' The convergence checks determine the quality standards that fitted models must meet.
#'
#' @examples
#' # Set custom autofit control with shorter time limit
#' autofit_ctrl <- set_autofit_control(
#'   max_time = list(time = 30, unit = "mins"),
#'   sample_extend = 2000
#' )
#'
#' # Set custom convergence checks with stricter criteria
#' conv_checks <- set_convergence_checks(
#'   max_Rhat = 1.01,
#'   min_ESS = 1000
#' )
#'
#' \dontrun{
#' # Use in RoBMA function
#' fit <- RoBMA(d = c(0.5, 0.3, 0.1),
#'              se = c(0.2, 0.15, 0.1),
#'              autofit_control = autofit_ctrl,
#'              convergence_checks = conv_checks)
#' }
#'
#' @return \code{set_autofit_control} returns a list of autofit control settings
#' and \code{set_convergence_checks} returns a list of convergence checks settings.
#'
#' @export set_autofit_control
#' @export set_convergence_checks
#' @name RoBMA_control
#' @aliases set_autofit_control set_convergence_checks
#'
#' @seealso [RoBMA()], [update.RoBMA()]
NULL

#' @rdname RoBMA_control
set_autofit_control     <- function(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10){

  autofit_settings <- list(
    max_Rhat      = max_Rhat,
    min_ESS       = min_ESS,
    max_error     = max_error,
    max_SD_error  = max_SD_error,
    max_time      = max_time,
    sample_extend = sample_extend,
    restarts      = restarts,
    max_extend    = max_extend
  )
  autofit_settings <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_settings, call = "Checking 'autofit_control':\n\t")

  return(autofit_settings)
}
#' @rdname RoBMA_control
set_convergence_checks  <- function(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL, remove_failed = FALSE, balance_probability = TRUE){

  convergence_checks <- list(
    max_Rhat            = max_Rhat,
    min_ESS             = min_ESS,
    max_error           = max_error,
    max_SD_error        = max_SD_error,
    remove_failed       = remove_failed,
    balance_probability = balance_probability
  )
  # allows NULL arguments so it can be used in this way too
  convergence_checks <- .check_and_list_convergence_checks(convergence_checks)

  return(convergence_checks)
}


# Internal control parameter update functions:
# Purpose: Merge user-specified parameters with existing control settings
# Used for updating model fitting, convergence checking, and autofit parameters
# LLM Note: These functions implement parameter inheritance for model updates

# Update MCMC fitting control parameters:
# Handles chains, adaptation, burnin, sampling parameters for JAGS
.update_fit_control     <- function(old_fit_control, chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed){

  if(is.null(chains)){
    chains <- old_fit_control[["chains"]]
  }
  if(is.null(adapt)){
    adapt  <- old_fit_control[["adapt"]]
  }
  if(is.null(burnin)){
    burnin <- old_fit_control[["burnin"]]
  }
  if(is.null(sample)){
    sample <- old_fit_control[["sample"]]
  }
  if(is.null(thin)){
    thin  <- old_fit_control[["thin"]]
  }
  if(is.null(autofit)){
    autofit  <- old_fit_control[["autofit"]]
  }
  if(is.null(parallel)){
    parallel <- old_fit_control[["parallel"]]
  }
  if(is.null(silent)){
    silent <- old_fit_control[["silent"]]
  }
  if(is.null(seed)){
    seed   <- old_fit_control[["seed"]]
  }

  new_fit_control <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)

  return(new_fit_control)
}
.update_autofit_control <- function(old_autofit_control, autofit_control){

  if(!is.null(autofit_control[["max_Rhat"]])){
    max_Rhat <- autofit_control[["max_Rhat"]]
  }else{
    max_Rhat <- old_autofit_control[["max_Rhat"]]
  }
  if(!is.null(autofit_control[["min_ESS"]])){
    min_ESS <- autofit_control[["min_ESS"]]
  }else{
    min_ESS <- old_autofit_control[["min_ESS"]]
  }
  if(!is.null(autofit_control[["max_error"]])){
    max_error <- autofit_control[["max_error"]]
  }else{
    max_error <- old_autofit_control[["max_error"]]
  }
  if(!is.null(autofit_control[["max_SD_error"]])){
    max_SD_error <- autofit_control[["max_SD_error"]]
  }else{
    max_SD_error <- old_autofit_control[["max_SD_error"]]
  }
  if(!is.null(autofit_control[["max_time"]])){
    max_time <- autofit_control[["max_time"]]
  }else{
    max_time <- old_autofit_control[["max_time"]]
  }
  if(!is.null(autofit_control[["sample_extend"]])){
    sample_extend <- autofit_control[["sample_extend"]]
  }else{
    sample_extend <- old_autofit_control[["sample_extend"]]
  }
  if(!is.null(autofit_control[["restarts"]])){
    restarts <- autofit_control[["restarts"]]
  }else{
    restarts <- old_autofit_control[["restarts"]]
  }

  new_autofit_control <- set_autofit_control(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, max_time = max_time, sample_extend = sample_extend, restarts = restarts)
  new_autofit_control <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = new_autofit_control)

  return(new_autofit_control)
}
.update_convergence_checks <- function(old_convergence_checks, convergence_checks){

  if(!is.null(convergence_checks[["max_Rhat"]])){
    max_Rhat <- convergence_checks[["max_Rhat"]]
  }else{
    max_Rhat <- old_convergence_checks[["max_Rhat"]]
  }
  if(!is.null(convergence_checks[["min_ESS"]])){
    min_ESS <- convergence_checks[["min_ESS"]]
  }else{
    min_ESS <- old_convergence_checks[["min_ESS"]]
  }
  if(!is.null(convergence_checks[["max_error"]])){
    max_error <- convergence_checks[["max_error"]]
  }else{
    max_error <- old_convergence_checks[["max_error"]]
  }
  if(!is.null(convergence_checks[["max_SD_error"]])){
    max_SD_error <- convergence_checks[["max_SD_error"]]
  }else{
    max_SD_error <- old_convergence_checks[["max_SD_error"]]
  }
  if(!is.null(convergence_checks[["remove_failed"]])){
    remove_failed <- convergence_checks[["remove_failed"]]
  }else{
    remove_failed <- old_convergence_checks[["remove_failed"]]
  }
  if(!is.null(convergence_checks[["balance_probability"]])){
    balance_probability <- convergence_checks[["balance_probability"]]
  }else{
    balance_probability <- old_convergence_checks[["balance_probability"]]
  }

  new_convergence_checks <- set_convergence_checks(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, remove_failed = remove_failed, balance_probability = balance_probability)
  new_convergence_checks <- .check_and_list_convergence_checks(new_convergence_checks)
}


.check_and_list_convergence_checks <- function(convergence_checks){

  remove_failed       <- convergence_checks[["remove_failed"]]
  balance_probability <- convergence_checks[["balance_probability"]]
  convergence_checks["remove_failed"]       <- NULL
  convergence_checks["balance_probability"] <- NULL
  convergence_checks <- BayesTools::JAGS_check_and_list_autofit_settings(convergence_checks, skip_sample_extend = TRUE, call = "Checking 'convergence_checks':\n\t")

  BayesTools::check_bool(remove_failed,       "remove_failed",       call = "Checking 'convergence_checks':\n\t")
  BayesTools::check_bool(balance_probability, "balance_probability", call = "Checking 'convergence_checks':\n\t")
  convergence_checks[["remove_failed"]]       <- remove_failed
  convergence_checks[["balance_probability"]] <- balance_probability
  return(convergence_checks)
}


### TODO: save helpers
.get_one_sided_cuts          <- function(prior){
  steps <- prior[["parameters"]][["steps"]]
  if (grepl("two.sided", prior[["distribution"]])) {
    return(c(1-rev(steps/2), steps/2))
  } else {
    return(steps)
  }
}
