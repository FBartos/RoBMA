#' @title Bayesian Selection Model
#'
#' @description Function for fitting random-effects, meta-regression, multilevel,
#' and location-scale meta-analytic selection models.
#'
#' @inheritParams data_input
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#' @param steps numeric vector of one-sided p-value cut points for the
#' selection model. If omitted, the default is `0.025`, yielding intervals
#' `[0, .025]` and `(.025, 1]`.
#'
#' @return A fitted object of class `c("bselmodel", "brma")` containing a
#' single Bayesian selection model fit.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- bselmodel(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     steps   = 0.025,
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary(fit)
#'   funnel(fit)
#' }
#' }
#'
#' @seealso [RoBMA()], [bPET()], [bPEESE()], [summary.brma()],
#' [funnel.brma()]
#' @export
bselmodel <- function(
    # input specification
  yi, vi, sei, weights, ni,
  mods, scale, cluster,
  data, slab, subset,
  measure = "GEN",

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale,
  prior_heterogeneity_allocation, prior_bias,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "treatment",
  prior_unit_information_sd, rescale_priors = 1,
  prior_informed_field, prior_informed_subfield,
  effect_direction = "detect", steps,

  # MCMC fitting settings
  sample = 5000, burnin = 2000, adapt = 500,
  chains = 3, thin = 1, parallel = FALSE,
  autofit = FALSE, autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),

  # additional settings
  seed = NULL, silent, ...
) {

  ### create the output object
  dots         <- list(...)
  .check_unused_dots(
    dots    = dots,
    allowed = c("only_data", "only_priors", "is_JASP", "is_JASP_prefix"),
    caller  = "bselmodel()"
  )
  object       <- .createObject(
    dots = dots, class = c("bselmodel", "brma"),
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
    data = object[["data"]], bias_type = "selmodel", steps = steps)
  if (isTRUE(dots[["only_priors"]]))
    return(.set_only_priors_class(object))

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
