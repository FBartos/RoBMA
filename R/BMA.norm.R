#' @title Bayesian Model-Averaged Meta-Analysis
#'
#' @description Convenience wrapper for fitting Bayesian model-averaged meta-analytic
#' models without publication bias adjustment. This function calls \code{\link{RoBMA}}
#' internally without the \code{model_type} argument, omitting selection models and
#' PET-PEESE bias adjustment.
#'
#' @inheritParams data_input
#' @inheritParams RoBMA_prior_specification
#' @inheritParams fitting_specification
#'
#' @details
#' \code{BMA.norm} (and its alias \code{BMA}) provides a simplified interface for
#' Bayesian model-averaged meta-analysis when publication bias adjustment is not needed.
#' It is equivalent to calling \code{\link{RoBMA}} without specifying \code{prior_bias},
#' \code{prior_bias_null}, or \code{model_type}.
#'
#' For publication bias adjusted meta-analysis, use \code{\link{RoBMA}} directly.
#'
#' @return \code{BMA.norm} returns an object of class \code{"BMA.norm"}.
#'
#' @seealso [RoBMA()] [brma.norm()] [summary.brma()] [plot.brma()]
#'
#' @export
BMA <- BMA.norm <- function(
  # input specification
  yi, vi, sei, weights, ni,
  mods, scale, study_ids,
  data, slab, subset,
  measure = "GEN",

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale,
  prior_heterogeneity_allocation,
  prior_effect_null, prior_heterogeneity_null, prior_mods_null,
  prior_scale_null, prior_heterogeneity_allocation_null,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "meandif",
  prior_unit_information_sd, rescale_priors = 1,
  prior_informed_field, prior_informed_subfield,

  # MCMC fitting settings
  sample = 5000, burnin = 2000, adapt = 500,
  chains = 3, thin = 1, parallel = FALSE,
  autofit = FALSE, autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),

  # additional settings
  seed = NULL, silent = TRUE, ...
) {

  ### create the output object
  time.start   <- proc.time()
  dots         <- list(...)
  object       <- .createObject(
    dots = dots, class = c("BMA.norm", "RoBMA", "brma"),
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
    measure = measure
  )
  if (isTRUE(dots[["only_data"]]))
    return(object)

  ### check and store priors (using RoBMA-style mixture priors, without bias)
  object$priors <- .check_and_list_priors.RoBMA(
    prior_effect = prior_effect, prior_heterogeneity = prior_heterogeneity,
    prior_mods = prior_mods, prior_scale = prior_scale,
    prior_heterogeneity_allocation = prior_heterogeneity_allocation,
    prior_effect_null = prior_effect_null, prior_heterogeneity_null = prior_heterogeneity_null,
    prior_mods_null = prior_mods_null, prior_scale_null = prior_scale_null,
    prior_heterogeneity_allocation_null = prior_heterogeneity_allocation_null,
    rescale_priors                    = rescale_priors,
    prior_unit_information_sd         = prior_unit_information_sd,
    prior_informed_field              = prior_informed_field,
    prior_informed_subfield           = prior_informed_subfield,
    data = object[["data"]]
  )
  if (isTRUE(dots[["only_priors"]]))
    return(object)

  ### fit the model
  object$fit <- .fit(object)

  ### store simple summary & coefficients
  object$summary       <- .object_summary(object)
  object$coefficients  <- .object_coefficients(object)

  return(object)
}