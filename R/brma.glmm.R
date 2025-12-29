### RoBMA 4.0.0

#' @export
brma.glmm <- function(
  # input specification
  ai, bi, ci, di, n1i, n2i, weights,
  mods, scale, study_ids,
  data, slab, subset,

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale,
  prior_heterogeneity_allocation,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd, rescale_priors = 1,
  prior_informed_field, prior_informed_subfield,

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
  object       <- .createObject(
    dots = dots, class = "brma.glmm",
    # MCMC and fitting settings
    chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
    autofit = autofit, parallel = parallel, silent = silent, seed = seed,
    autofit_control = autofit_control, convergence_checks = convergence_checks,
    # additional options
    standardize_continuous_predictors = standardize_continuous_predictors
  )

  ### check and store the data
  object$data <- .check_and_list_data(.call = match.call(), .envir = parent.frame(), class = "glmm")
  if (isTRUE(dots[["only_data"]]))
    return(object)

}
