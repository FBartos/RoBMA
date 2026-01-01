### RoBMA 4.0.0

### creates basic RoBMA object with
# - fit_control
# - autofit_control
# - convergence_checks
#
# @param dots additional arguments originally passed as ...
# @param chains Number of MCMC chains
# @param adapt Number of adaptation iterations
# @param burnin Number of burnin iterations
# @param sample Number of sampling iterations
# @param thin Thinning interval
# @param autofit Whether to use autofit
# @param parallel Whether to run chains in parallel
# @param silent Whether to suppress output
# @param seed Random seed
# @param autofit_control Autofit control settings
# @param convergence_checks Convergence check settings
#
# @return A list containing fit_control, autofit_control, convergence_checks
.createObject <- function(
    dots, class,
    chains, adapt, burnin, sample, thin,
    autofit, parallel, silent, seed,
    autofit_control, convergence_checks) {

  object       <- NULL

  ### input global settings if unspecified
  if (missing(silent)) {
    silent <- RoBMA.get_option("silent")
  }

  ### check and store MCMC settings
  object$fit_control <- BayesTools::JAGS_check_and_list_fit_settings(
    chains = chains, adapt = adapt, burnin = burnin, sample = sample,
    thin = thin, autofit = autofit, parallel = parallel, cores = chains,
    silent = silent, seed = seed
  )
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks = convergence_checks)


  ### include JASP indicators for progress bars
  if (!is.null(dots[["is_JASP"]])) {
    object[["is_JASP"]]        <- dots[["is_JASP"]]
    object[["is_JASP_prefix"]] <- dots[["is_JASP_prefix"]]
  }


  ### add class
  class(object) <- class

  return(object)
}


### object tools options
# add simple summary and model coefficients to the object
# (this differ from more customizable user facing summary function)
.object_summary      <- function(object) {

  # provide a simple summary
  estimates <- BayesTools::JAGS_estimates_table(
    fit               = object[["fit"]],
    transform_factors = TRUE,
    transform_scaled  = TRUE,
    remove_parameters = c(
      "theta", # remove random-effects (estimate-level)
      "gamma", # remove random-effects (study-level)
      "pi",    # remove baserate for OR models
      "phi"    # remove lograte for IRR models
    )
  )

  return(estimates)
}
.object_coefficients <- function(object) {

  estimates        <- object[["summary"]][,"Mean"]
  names(estimates) <- rownames(object[["summary"]])

  return(estimates)
}
