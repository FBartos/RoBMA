#' @title Bayesian Model-Averaged Generalized Meta-Analysis
#'
#' @description Function for fitting Bayesian model-averaged meta-analytic models
#' directly to binary or count data using a generalized linear mixed model (GLMM)
#' framework. Unlike \code{\link{RoBMA}}, this function does not adjust for
#' publication bias, as weight function and regression-based bias adjustment
#' methods are not available for GLMM outcomes.
#'
#' @inheritParams data_input
#' @inheritParams RoBMA_prior_specification
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#'
#' @details
#' \code{BMA.glmm} combines the data input style of \code{\link{brma.glmm}} with
#' the mixture prior specification of \code{\link{RoBMA}} for Bayesian model-averaging.
#'
#' Model for odds ratios (\code{measure = "OR"}) uses a binomial-normal model
#' as described in \insertCite{jackson2018comparison;textual}{RoBMA}.
#'
#' Model for incidence rate ratios (\code{measure = "IRR"}) uses a Poisson-normal
#' model as described in \insertCite{bagos2009mixed;textual}{RoBMA}.
#'
#' When \code{weights} are supplied, they are treated as likelihood weights on
#' the paired two-arm study contribution.
#'
#' @return A fitted object of class
#' `c("BMA.glmm", "RoBMA", "brma.glmm", "brma")`. The object contains checked
#' `data`, checked mixture `priors`, the JAGS `fit`, cached `summary`, and
#' cached `coefficients`.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'
#'   fit <- BMA.glmm(
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
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
#' @seealso [brma.glmm()] [RoBMA()] [summary.brma()] [plot.brma()]
#' @export
BMA.glmm <- function(
  # input specification
  ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, weights,
  mods, scale, cluster,
  data, slab, subset,
  measure = "OR",

  # prior specification
  prior_effect, prior_heterogeneity, prior_mods, prior_scale,
  prior_heterogeneity_allocation, prior_baserate, prior_lograte,
  prior_effect_null, prior_heterogeneity_null, prior_mods_null,
  prior_scale_null, prior_heterogeneity_allocation_null,
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
) {

  ### create the output object
  dots         <- list(...)
  dots         <- .validate_constructor_dots(dots, caller = "BMA.glmm()")
  object       <- .createObject(
    dots = dots, class = c("BMA.glmm", "RoBMA", "brma.glmm", "brma"),
    # MCMC and fitting settings
    chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
    autofit = autofit, parallel = parallel, silent = silent, seed = seed,
    autofit_control = autofit_control, convergence_checks = convergence_checks
  )

  ### check and store the data
  object$data <- .check_and_list_data(
    .call = match.call(), .envir = parent.frame(), class = "glmm",
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
    prior_baserate = prior_baserate, prior_lograte = prior_lograte,
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
