#' @title Options for the RoBMA package
#'
#' @description Inspect or change package-level defaults used by RoBMA.
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
#'   \item{\code{cluster_likelihood.n_gamma}}{number of Gauss-Hermite nodes used for cluster-unit log-likelihoods (default \code{15})}
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
#' @return `RoBMA.options()` invisibly returns a named list with all current
#' options after applying any changes. `RoBMA.get_option()` returns the current
#' value of the requested option.
#'
#' @export RoBMA.options
#' @export RoBMA.get_option
#' @name RoBMA_options
#' @aliases RoBMA_options RoBMA.options RoBMA.get_option
NULL


#' @rdname RoBMA_options
RoBMA.options    <- function(...) {

  opts <- list(...)

  if (length(opts) > 0 && (is.null(names(opts)) || any(names(opts) == ""))) {
    stop("All options must be named.", call. = FALSE)
  }

  option_names <- RoBMA.private[["option_names"]]

  if (length(opts) > 0 && anyDuplicated(names(opts))) {
    stop("Option names must be unique.", call. = FALSE)
  }

  for (i in seq_along(opts)) {

    if (!names(opts)[i] %in% option_names) {
      stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep = ""), call. = FALSE)
    }

    value <- .RoBMA_validate_option(names(opts)[i], opts[[i]])
    assign(names(opts)[i], value, envir = RoBMA.private)
  }

  return(invisible(.RoBMA_current_options()))
}

#' @rdname RoBMA_options
RoBMA.get_option <- function(name) {

  if (!is.character(name) || length(name) != 1 || is.na(name)) {
    stop("Only 1 option can be retrieved at a time", call. = FALSE)
  }

  option_names <- RoBMA.private[["option_names"]]

  if (!name %in% option_names) {
    stop(paste("Unmatched or ambiguous option '", name, "'", sep = ""), call. = FALSE)
  }

  return(RoBMA.private[[name]])
}

# export the function directly to suppress import warnings
.runjags__findjags <- function() runjags::findjags()

# adapted from the runjags package version 2.2.0
RoBMA.private <- new.env()

assign("defaultoptions",  list(
  jagspath = expression(.runjags__findjags())
), envir = RoBMA.private)

assign("options",         RoBMA.private$defaultoptions,   envir = RoBMA.private)
assign("RoBMA_version",   "notset",                       envir = RoBMA.private)
assign("min_jags_major",  4,                              envir = RoBMA.private)
assign("max_jags_major",  4,                              envir = RoBMA.private)

.RoBMA_default_max_cores <- function() {

  cores <- parallel::detectCores(logical = TRUE)
  if (length(cores) != 1 || is.na(cores) || !is.finite(cores)) {
    return(1L)
  }

  return(max(1L, as.integer(cores) - 1L))
}

.RoBMA_check_option_bool <- function(value, name) {

  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop(paste0("Option '", name, "' must be TRUE or FALSE."), call. = FALSE)
  }

  return(value)
}

.RoBMA_check_option_int <- function(value, name, lower = -Inf) {

  if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
      !is.finite(value) || value != as.integer(value) || value < lower) {
    stop(paste0("Option '", name, "' must be an integer >= ", lower, "."), call. = FALSE)
  }

  return(as.integer(value))
}

.RoBMA_check_option_real <- function(value, name, lower = -Inf, allow_bound = TRUE) {

  if (!is.numeric(value) || length(value) != 1 || is.na(value) || !is.finite(value)) {
    stop(paste0("Option '", name, "' must be a finite number."), call. = FALSE)
  }

  below_bound <- if (allow_bound) value < lower else value <= lower
  if (below_bound) {
    bound_text <- if (allow_bound) " >= " else " > "
    stop(paste0("Option '", name, "' must be", bound_text, lower, "."), call. = FALSE)
  }

  return(as.numeric(value))
}

.RoBMA_option_schema <- list(
  "max_cores" = list(
    default  = .RoBMA_default_max_cores,
    validate = function(value, name) .RoBMA_check_option_int(value, name, lower = 1L)
  ),
  "check_scaling" = list(
    default  = TRUE,
    validate = .RoBMA_check_option_bool
  ),
  "silent" = list(
    default  = FALSE,
    validate = .RoBMA_check_option_bool
  ),
  "autocompute.loo" = list(
    default  = FALSE,
    validate = .RoBMA_check_option_bool
  ),
  "autocompute.waic" = list(
    default  = FALSE,
    validate = .RoBMA_check_option_bool
  ),
  "autocompute.marglik" = list(
    default  = FALSE,
    validate = .RoBMA_check_option_bool
  ),
  "cluster_likelihood.n_gamma" = list(
    default  = 15L,
    validate = function(value, name) .RoBMA_check_option_int(value, name, lower = 3L)
  ),
  "default_UISD.effect" = list(
    default  = 1/2,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_UISD.heterogeneity" = list(
    default  = 1/4,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_UISD.mods" = list(
    default  = 1/4,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_UISD.scale" = list(
    default  = 1/2,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_informed_priors.mods" = list(
    default  = 1/2,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_informed_priors.scale" = list(
    default  = 1/2,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_bias_weightfunction.alpha" = list(
    default  = 1,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_bias_PET.scale" = list(
    default  = 1,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  ),
  "default_bias_PEESE.scale" = list(
    default  = 5,
    validate = function(value, name) .RoBMA_check_option_real(value, name, lower = 0, allow_bound = FALSE)
  )
)

assign("option_names", names(.RoBMA_option_schema), envir = RoBMA.private)

for (option_name in names(.RoBMA_option_schema)) {

  default <- .RoBMA_option_schema[[option_name]][["default"]]
  value   <- if (is.function(default)) default() else default
  assign(option_name, value, envir = RoBMA.private)
}

.RoBMA_validate_option <- function(name, value) {

  validator <- .RoBMA_option_schema[[name]][["validate"]]
  return(validator(value, name))
}


.RoBMA_current_options <- function() {

  option_names <- RoBMA.private[["option_names"]]
  out          <- mget(option_names, envir = RoBMA.private, inherits = FALSE)

  return(out)
}


# check and fix number of threads (sometimes bugs out during installation)
.check_max_cores <- function() {

  max_cores <- .RoBMA_default_max_cores()
  if (RoBMA.private$max_cores > max_cores || RoBMA.private$max_cores < 1) {
    RoBMA.options(max_cores = max_cores)
  }
}


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
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- RoBMA(
#'     yi                 = yi,
#'     vi                 = vi,
#'     data               = dat.lehmann2018,
#'     measure            = "SMD",
#'     autofit_control    = autofit_ctrl,
#'     convergence_checks = conv_checks
#'   )
#' }
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
#' @seealso [RoBMA()]
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
  if(is.null(cores)){
    cores <- old_fit_control[["cores"]]
  }
  if(is.null(silent)){
    silent <- old_fit_control[["silent"]]
  }
  if(is.null(seed)){
    seed   <- old_fit_control[["seed"]]
  }

  new_fit_control <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = cores, silent = silent, seed = seed)

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
  if(!is.null(autofit_control[["max_extend"]])){
    max_extend <- autofit_control[["max_extend"]]
  }else{
    max_extend <- old_autofit_control[["max_extend"]]
  }

  new_autofit_control <- set_autofit_control(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, max_time = max_time, sample_extend = sample_extend, restarts = restarts, max_extend = max_extend)
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

  return(new_convergence_checks)
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

# Apply shared metafor-style defaults for observed-study point markers.
.plot_point_style_defaults <- function(dots) {

  if (is.null(dots[["pch"]]))  dots[["pch"]]  <- 21
  if (is.null(dots[["col"]]))  dots[["col"]]  <- "black"
  if (is.null(dots[["bg"]]))   dots[["bg"]]   <- "#A6A6A6"
  if (is.null(dots[["cex"]]))  dots[["cex"]]  <- 1
  if (is.null(dots[["size"]])) dots[["size"]] <- 2

  return(dots)
}
