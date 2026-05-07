#' @title Bayesian Selection Model
#'
#' @description Function for fitting random-effects, meta-regression, multilevel,
#' and location-scale meta-analytic selection models.
#'
#' @inheritParams data_input
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#' @param prior_bias selection-model bias prior, usually created by
#' \code{prior_weightfunction()}. If omitted or \code{NULL}, a default
#' one-sided weightfunction prior is constructed from \code{steps}.
#' @param steps numeric vector of one-sided p-value cut points for the
#' default selection model. If `prior_bias` is supplied, the prior carries its
#' own side, steps, and weights. If omitted, the default is `0.025`, yielding
#' intervals `[0, .025]` and `(.025, 1]`.
#'
#' @details
#' `bselmodel()` is a normal/effect-size selection-model constructor. Custom
#' `prior_bias` can be a weightfunction prior or a supported BayesTools
#' selection-kernel prior; p-hacking kernels are not supported in active RoBMA.
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
#' @seealso [publication_bias_prior_specification], [RoBMA()], [bPET()],
#' [bPEESE()], [summary.brma()], [funnel.brma()]
#' @export
bselmodel <- function(
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
  dots            <- list(...)
  missing_measure <- missing(measure)
  if (missing_measure && !isTRUE(dots[["only_data"]])) {
    .stop_missing_measure("bselmodel()")
  }
  if (missing_measure) {
    measure <- "GEN"
  }
  dots            <- .validate_constructor_dots(dots, caller = "bselmodel()")
  object          <- .createObject(
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

  object               <- .autocompute_brma(object)

  return(object)
}
