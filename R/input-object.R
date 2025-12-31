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
# @param standardize_continuous_predictors Whether to standardize continuous predictors
#
# @return A list containing fit_control, autofit_control, convergence_checks
.createObject <- function(
    dots, class,
    chains, adapt, burnin, sample, thin,
    autofit, parallel, silent, seed,
    autofit_control, convergence_checks,
    standardize_continuous_predictors) {

  object       <- NULL

  ### check and store MCMC settings
  object$fit_control <- BayesTools::JAGS_check_and_list_fit_settings(
    chains = chains, adapt = adapt, burnin = burnin, sample = sample,
    thin = thin, autofit = autofit, parallel = parallel, cores = chains,
    silent = silent, seed = seed
  )
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks = convergence_checks)

  ### dealt with additional settings
  # automatic scaling of continuous predictors: used directly in BayesTools::JAGS_fit
  BayesTools::check_bool(standardize_continuous_predictors, "standardize_continuous_predictors", allow_NA = FALSE)
  object$standardize_continuous_predictors <- standardize_continuous_predictors

  ### include JASP indicators for progress bars
  if (!is.null(dots[["is_JASP"]])) {
    object[["is_JASP"]]        <- dots[["is_JASP"]]
    object[["is_JASP_prefix"]] <- dots[["is_JASP_prefix"]]
  }


  ### add class
  class(object) <- c("brma", class)

  return(object)
}
