#' @title Bayesian Model-Averaged Meta-Analysis
#'
#' @description Fits Bayesian model-averaged meta-analytic models without
#' publication-bias adjustment. `BMA()` is an alias for `BMA.norm()`.
#'
#' @inheritParams data_input
#' @inheritParams RoBMA_prior_specification
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#'
#' @details
#' \code{BMA.norm} (and its alias \code{BMA}) provides a simplified interface for
#' Bayesian model-averaged meta-analysis when publication bias adjustment is not needed.
#' It uses the same product-space mixture-prior machinery as `RoBMA()` but
#' omits selection models and PET-PEESE bias adjustment.
#'
#' For publication bias adjusted meta-analysis, use \code{\link{RoBMA}} directly.
#'
#' @return A fitted object of class `c("BMA.norm", "RoBMA", "brma")`.
#' The object contains checked `data`, checked mixture `priors`, the JAGS
#' `fit`, cached `summary`, and cached `coefficients`.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- BMA(
#'     yi      = yi,
#'     vi      = vi,
#'     mods    = ~ Preregistered,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary(fit)
#' }
#' }
#'
#' @seealso [RoBMA()] [brma()] [summary.brma()] [plot.brma()]
#'
#' @aliases BMA.norm
#' @export BMA.norm
#' @export
BMA <- BMA.norm <- function(
  # input specification
  yi, vi, sei, weights, ni,
  mods, scale, cluster,
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
  dots         <- list(...)
  dots         <- .validate_constructor_dots(dots, caller = "BMA.norm()")
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
    return(.set_only_priors_class(object))

  ### fit the model
  object$fit <- .fit(object)

  ### store simple summary & coefficients
  object$summary       <- .object_summary(object)
  object$coefficients  <- .object_coefficients(object)
  object               <- .autocompute_brma(object)

  return(object)
}
