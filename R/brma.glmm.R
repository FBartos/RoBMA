#' @title Bayesian Generalized Meta-Analysis
#'
#' @description Function for fitting random-effects, meta-regression, multilevel,
#' and location-scale meta-analytic models directly to either binary or count data.
#'
#' @inheritParams data_input
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#' @details
#' Model for odds ratios (`measure = "OR"`) corresponds to Model 4 described in
#' \insertCite{jackson2018comparison;textual}{RoBMA}.
#' `logit(pi[i])` is the study-specific midpoint of the two arm logits.
#' `prior_baserate` defines the estimate-specific prior distribution on `pi[i]`.
#'
#' Model for incidence rate ratios (`measure = "IRR"`) corresponds to
#' \insertCite{bagos2009mixed;textual}{RoBMA}.
#' `phi[i]` is the study-specific midpoint of the two arm log incidence rates.
#' `prior_lograte` defines the estimate-specific prior distribution on `phi[i]`.
#' If unspecified, a unit-information prior is based on the data and used
#' independently for each estimate.
#'
#' When `weights` are supplied, they are treated as likelihood weights on the
#' paired two-arm study contribution.
#'
#' @return A fitted object of class `c("brma.glmm", "brma")`. The object
#' contains checked `data`, checked `priors`, the JAGS `fit`, cached `summary`,
#' and cached `coefficients`. If the corresponding package options are enabled,
#' it can also contain cached LOO, WAIC, or marginal likelihood results.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'
#'   fit <- brma.glmm(
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
#'     mods    = ~ alloc,
#'     data    = dat.bcg,
#'     measure = "OR",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary(fit)
#' }
#' }
#'
#' @seealso [brma()], [BMA.glmm()], [summary.brma()], [predict.brma()]
#' @export
brma.glmm <- function(
  # input specification
  ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, weights,
  mods, scale, cluster,
  data, slab, subset,
  measure = "OR",

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale,
  prior_heterogeneity_allocation, prior_baserate, prior_lograte,
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
  seed = NULL, silent, ...
){

  ### create the output object
  dots         <- list(...)
  dots         <- .validate_constructor_dots(dots, caller = "brma.glmm()")
  object       <- .createObject(
    dots = dots, class = c("brma.glmm", "brma"),
    # MCMC and fitting settings
    chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
    autofit = autofit, parallel = parallel, silent = silent, seed = seed,
    autofit_control = autofit_control, convergence_checks = convergence_checks
  )

  ### check and store the data
  object$data <- .check_and_list_data(
    .call = match.call(), .envir = parent.frame(), class = "glmm",
    set_contrast_factor_predictors = set_contrast_factor_predictors,
    measure = measure, standardize_continuous_predictors = standardize_continuous_predictors)
  if (isTRUE(dots[["only_data"]]))
    return(object)

  ### check and store priors
  object$priors <- .check_and_list_priors.brma(
    prior_effect = prior_effect, prior_heterogeneity = prior_heterogeneity,
    prior_mods = prior_mods, prior_scale = prior_scale,
    prior_heterogeneity_allocation = prior_heterogeneity_allocation,
    prior_baserate = prior_baserate, prior_lograte = prior_lograte,
    rescale_priors                    = rescale_priors,
    prior_unit_information_sd         = prior_unit_information_sd,
    prior_informed_field              = prior_informed_field,
    prior_informed_subfield           = prior_informed_subfield,
    data = object[["data"]])
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
