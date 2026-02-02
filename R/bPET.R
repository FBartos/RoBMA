#' @title Bayesian Precision-Effect Test (PET) Model
#'
#' @description Function for fitting random-effects, meta-regression, multilevel,
#' and location-scale meta-analytic PET models.
#'
#' @inheritParams data_input
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#' @export
bPET <- function(
  # input specification
  yi, vi, sei, weights, ni,
  mods, scale, study_ids,
  data, slab, subset,
  measure = "GEN",

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale,
  prior_heterogeneity_allocation, prior_bias,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd, rescale_priors = 1,
  prior_informed_field, prior_informed_subfield,
  effect_direction = "detect",

  # MCMC fitting settings
  sample = 5000, burnin = 2000, adapt = 500,
  chains = 3, thin = 1, parallel = FALSE,
  autofit = FALSE, autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),

  # additional settings
  seed = NULL, silent, ...
) {

  ### create the output object
  time.start   <- proc.time()
  dots         <- list(...)
  object       <- .createObject(
    dots = dots, class = c("bPET", "brma"),
    # MCMC and fitting settings
    chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
    autofit = autofit, parallel = parallel, silent = silent, seed = seed,
    autofit_control = autofit_control, convergence_checks = convergence_checks
  )

  ### check and store the data
  object$data <- .check_and_list_data(
    .call = match.call(), .envir = parent.frame(), class = "norm",
    set_contrast_factor_predictors = set_contrast_factor_predictors,
    standardize_continuous_predictors = standardize_continuous_predictors,
    effect_direction = effect_direction, measure = measure)
  if (isTRUE(dots[["only_data"]]))
    return(object)

  ### check and store priors
  object$priors <- .check_and_list_priors.brma(
    prior_effect = prior_effect, prior_heterogeneity = prior_heterogeneity,
    prior_mods = prior_mods, prior_scale = prior_scale,
    prior_heterogeneity_allocation = prior_heterogeneity_allocation,
    prior_bias = prior_bias,
    rescale_priors                    = rescale_priors,
    prior_unit_information_sd         = prior_unit_information_sd,
    prior_informed_field              = prior_informed_field,
    prior_informed_subfield           = prior_informed_subfield,
    data = object[["data"]], bias_type = "PET")
  if (isTRUE(dots[["only_priors"]]))
    return(object)

  ### fit the model
  object$fit <- .fit(object)

  ### store simple summary & coefficients
  object$summary       <- .object_summary(object)
  object$coefficients  <- .object_coefficients(object)

  ### autocompute
  if (RoBMA.get_option("autocompute.loo")) {
    object <- add_loo(object)
  }
  if (RoBMA.get_option("autocompute.waic")) {
    object <- add_waic(object)
  }
  if (RoBMA.get_option("autocompute.marglik")) {
    object <- add_marglik(object)
  }

  return(object)
}
