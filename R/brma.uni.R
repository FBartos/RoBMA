
#' @export
brma.uni <- function(
    # input specification
    yi, vi, sei, ni, weights,
    mods, scale, study_ids,
    data, slab, subset,
    measure = "GEN",

    # prior specification
    rescale_continuous_predictors = TRUE,
    contrast_factor_predictors    = "treatment",
    priors = NULL, rescale_priors = 1,

    # MCMC fitting settings
    sample = 5000, burnin = 2000, adapt = 500,
    chains = 3, thin = 1, parallel = FALSE,
    autofit = FALSE, autofit_control = set_autofit_control(),
    convergence_checks = set_convergence_checks(),

    # additional settings
    seed = NULL, silent = TRUE, ...

    ){

  ### create the output object
  time.start   <- proc.time()
  dots         <- list(...)
  object       <- NULL
  object$call  <- match.call()

  ### check and store MCMC settings
  object$fit_control        <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)

  ## check and store the data
  object$data <- .check_and_list_data(.call  = object$call, .envir = parent.frame())
  if (isTRUE(dots[["only_data"]]))
    return(object[["data"]])

  ## check and store priors
  object$priors <- .check_and_list_priors.brma(
    rescale_continuous_predictors = rescale_continuous_predictors,
    contrast_factor_predictors    = contrast_factor_predictors,
    priors = priors, rescale_priors = rescale_priors,
    data = data, measure = measure
  )
  if (isTRUE(dots[["only_priors"]]))
    return(object[["priors"]])

  ## fit the model
  object$fit <- .fit.brma(object)



  object <- list(
    info      = info
  )


}




# ----


# -----
.check_measure <- function(measure, allow_GEN = TRUE) {

  BayesTools::check_char(measure, allow_values = c("SMD", "ZCOR", "RR", "OR", "HR", "GEN"))
  if (!allow_GEN && measure == "GEN")
    stop("'GEN' measure is not available for a given operation")

  return()
}

### prior scale functions
# returns the sd of unit information
# based on Chapter 2.4 in Spiegelhalter, Abrams, and Myles 2004 and Chapter 1 in Grieve 2022
# as reported in Table 2 of Pawel, S., & Held, L. (2025). Closed-form power and sample size calculations for Bayes factors. The American Statistician 9(3), 330-344, 10.1080/00031305.2025.2467919
.get_default_unit_information_sd <- function(measure) {

  .check_measure(measure, allow_GEN = FALSE)

  return(switch(
    measure,
    "SMD"  = sqrt(2),
    "ZCOR" = 1,
    "RR"   = sqrt(4),
    "OR"   = sqrt(4),
    "HR"   = sqrt(4)
  ))
}

# computes the unit information based on data
# based on Eq. 6 of Röver, C., Bender, R., Dias, S., Schmid, C. H., Schmidli, H., Sturtz, S., ... & Friede, T. (2021). On weakly informative prior distributions for the heterogeneity parameter in Bayesian random‐effects meta‐analysis. Research Synthesis Methods, 12(4), 448-474. 10.1002/jrsm.1475
.compute_unit_information_sd     <- function(sei, ni) {

  UISD <- sqrt(sum(ni) / sum(sei^(-2)))

  return(UISD)
}
