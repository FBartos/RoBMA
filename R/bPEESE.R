#' @title Bayesian Precision-Effect Estimate with Standard Errors (PEESE) Model
#'
#' @description Function for fitting random-effects, meta-regression, multilevel,
#' and location-scale meta-analytic PEESE models.
#'
#' @inheritParams data_input
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#' @param prior_bias PEESE regression-coefficient prior created by
#' \code{prior_PEESE()}. If omitted or \code{NULL}, the default is a Cauchy
#' prior centered at 0, truncated to positive values, with scale from
#' \code{RoBMA.get_option("default_bias_PEESE.scale")} after UISD adjustment.
#'
#' @return A fitted object of class `c("bPEESE", "brma")` containing a single
#' PEESE publication-bias model fit.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- bPEESE(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary(fit)
#'   funnel(fit)
#' }
#' }
#'
#' @seealso [publication_bias_prior_specification], [RoBMA()], [bPET()],
#' [bselmodel()], [summary.brma()], [funnel.brma()]
#' @export
bPEESE <- function(
    # input specification
  yi, vi, sei, weights, ni,
  mods, scale, cluster,
  data, slab, subset,
  measure,

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
  dots            <- list(...)
  missing_measure <- missing(measure)
  if (missing_measure && !isTRUE(dots[["only_data"]])) {
    .stop_missing_measure("bPEESE()")
  }
  if (missing_measure) {
    measure <- "GEN"
  }
  dots            <- .validate_constructor_dots(dots, caller = "bPEESE()")
  object          <- .createObject(
    dots = dots, class = c("bPEESE", "brma"),
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
    data = object[["data"]], bias_type = "PEESE")
  if (isTRUE(dots[["only_priors"]]))
    return(.set_only_priors_class(object))

  ### fit the model
  object$fit <- .fit(object)
  .stop_fit_errors(object$fit)

  ### store simple summary & coefficients
  object$summary       <- .object_summary(object)
  object$coefficients  <- .object_coefficients(object)

  object               <- .autocompute_brma(object)

  return(object)
}
