#' @title Robust Bayesian Model-Averaged Meta-Analysis
#'
#' @description Fits a robust Bayesian model-averaged meta-analysis. The
#' default ensemble averages across models with and without an effect,
#' heterogeneity, and publication-bias adjustment.
#'
#' @inheritParams data_input
#' @inheritParams RoBMA_prior_specification
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#'
#' @details
#' `RoBMA()` uses product-space Bayesian model averaging. Inclusion Bayes
#' factors and model-averaged estimates are obtained from mixture priors for
#' effect, heterogeneity, moderators, scale regression, and publication-bias
#' components.
#'
#' By default, `model_type = "PSMA"` includes selection-model weight functions
#' together with PET and PEESE publication-bias adjustments. Use `BMA()` for
#' model averaging without publication-bias adjustment, or `brma()` for fitting
#' a single meta-analytic model.
#'
#' @return A fitted object of class `c("RoBMA", "brma")`. The object contains
#' checked `data`, checked `priors`, the JAGS `fit`, cached `summary`, and
#' cached `coefficients`. It can be passed to `summary()`, `plot()`,
#' `predict()`, `funnel()`, `add_loo()`, and related methods.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- RoBMA(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary(fit)
#'   plot(fit)
#' }
#' }
#'
#' @seealso [BMA()], [brma()], [bselmodel()], [bPET()], [bPEESE()],
#' [summary.brma()], [plot.brma()]
#' @export
RoBMA <- function(
  # input specification
  yi, vi, sei, weights, ni,
  mods, scale, cluster,
  data, slab, subset,
  measure = "GEN", effect_direction = "detect",

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale, prior_heterogeneity_allocation, prior_bias,
  prior_effect_null, prior_heterogeneity_null, prior_mods_null, prior_scale_null, prior_heterogeneity_allocation_null, prior_bias_null,
  standardize_continuous_predictors = TRUE,
  set_contrast_factor_predictors = "meandif",
  prior_unit_information_sd, rescale_priors = 1,
  prior_informed_field, prior_informed_subfield,
  model_type = "PSMA",

  # MCMC fitting settings
  sample = 5000, burnin = 2000, adapt = 500,
  chains = 3, thin = 1, parallel = FALSE,
  autofit = FALSE, autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),

  # additional settings
  seed = NULL, silent = TRUE, ...) {

  ### create the output object
  dots         <- list(...)
  .check_unused_dots(
    dots    = dots,
    allowed = c("only_data", "only_priors", "is_JASP", "is_JASP_prefix"),
    caller  = "RoBMA()"
  )
  object       <- .createObject(
    dots = dots, class = c("RoBMA", "brma"),
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
  # checks and store the base priors
  object$priors <- .check_and_list_priors.RoBMA(
    prior_effect = prior_effect, prior_heterogeneity = prior_heterogeneity,
    prior_mods = prior_mods, prior_scale = prior_scale,
    prior_heterogeneity_allocation = prior_heterogeneity_allocation, prior_bias = prior_bias,
    prior_effect_null = prior_effect_null, prior_heterogeneity_null = prior_heterogeneity_null,
    prior_mods_null = prior_mods_null, prior_scale_null = prior_scale_null,
    prior_heterogeneity_allocation_null = prior_heterogeneity_allocation_null, prior_bias_null = prior_bias_null,
    rescale_priors                    = rescale_priors,
    prior_unit_information_sd         = prior_unit_information_sd,
    prior_informed_field              = prior_informed_field,
    prior_informed_subfield           = prior_informed_subfield,
    data = object[["data"]], model_type = model_type)
  if (isTRUE(dots[["only_priors"]]))
    return(.set_only_priors_class(object))

  ### fit the model
  object$fit <- .fit(object)

  ### store simple summary & coefficients
  object$summary       <- .object_summary(object)
  object$coefficients  <- .object_coefficients(object)

  return(object)
}
