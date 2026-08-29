#' @title Bayesian Multivariate and Multilevel PET Model
#'
#' @description
#' Fits Bayesian precision-effect test (PET) models with known sampling
#' covariance and BayesTools random-effect formulas.
#'
#' @inheritParams brma.mv
#' @inheritParams bPET
#'
#' @details
#' `bPET.mv()` combines the PET publication-bias regression of [bPET()] with
#' the known-covariance and random-formula model of [brma.mv()]. The PET
#' predictor is the marginal sampling standard error `sqrt(diag(V))`; known
#' off-diagonal sampling covariance remains part of the joint likelihood.
#'
#' The supported known-`V` backends, random-effect marginalization rules,
#' predictive targets, and diagnostic semantics are the same as for
#' [brma.mv()]. In particular, `random` replaces the specialized `cluster`
#' argument used by [bPET()].
#'
#' @return A fitted object of class
#' `c("bPET.mv", "bPET", "brma.mv", "brma.norm", "brma")`.
#'
#' @examples \dontrun{
#' V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
#' fit <- bPET.mv(
#'   yi      = c(0.10, 0.20),
#'   V       = V,
#'   measure = "GEN",
#'   seed    = 1,
#'   silent  = TRUE
#' )
#' summary(fit)
#' }
#'
#' @seealso [bPET()], [bPEESE.mv()], [brma.mv()], [summary.brma()]
#'
#' @export
bPET.mv <- function(
    # input specification
    yi, V, ni,
    mods, scale, random,
    R = NULL, Rscale = "cor",
    data, slab, subset,
    measure,

    # prior specification
    prior_effect, prior_heterogeneity, prior_mods, prior_scale, prior_bias,
    standardize_continuous_predictors = TRUE,
    set_contrast_factor_predictors = "treatment",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,
    effect_direction = "detect",

    # MCMC fitting settings
    known_v_parameterization = "auto",
    known_v_residual_fraction = 0.10,
    marginalize_estimate_level = TRUE,
    sample = 5000, burnin = 2000, adapt = 500,
    chains = 3, thin = 1, parallel = FALSE,
    autofit = FALSE, autofit_control = set_autofit_control(),
    convergence_checks = set_convergence_checks(),

    # additional settings
    seed = NULL, silent, ...,
    vi = NULL, sei = NULL) {

  .fit_bias_regression_mv(
    matched_call_unevaluated            = match.call(expand.dots = FALSE),
    matched_call                        = match.call(),
    envir                               = parent.frame(),
    caller                              = "bPET.mv()",
    class                               = c(
      "bPET.mv", "bPET", "brma.mv", "brma.norm", "brma"
    ),
    bias_type                           = "PET",
    known_v_residual_fraction_specified = !missing(
      known_v_residual_fraction
    ),
    R                                   = R,
    Rscale                              = Rscale,
    measure                             = measure,
    prior_effect                        = prior_effect,
    prior_heterogeneity                 = prior_heterogeneity,
    prior_mods                          = prior_mods,
    prior_scale                         = prior_scale,
    prior_bias                          = prior_bias,
    standardize_continuous_predictors   = standardize_continuous_predictors,
    set_contrast_factor_predictors      = set_contrast_factor_predictors,
    prior_unit_information_sd           = prior_unit_information_sd,
    rescale_priors                      = rescale_priors,
    prior_informed_field                = prior_informed_field,
    prior_informed_subfield             = prior_informed_subfield,
    effect_direction                    = effect_direction,
    known_v_parameterization            = known_v_parameterization,
    known_v_residual_fraction           = known_v_residual_fraction,
    marginalize_estimate_level          = marginalize_estimate_level,
    sample                              = sample,
    burnin                              = burnin,
    adapt                               = adapt,
    chains                              = chains,
    thin                                = thin,
    parallel                            = parallel,
    autofit                             = autofit,
    autofit_control                     = autofit_control,
    convergence_checks                  = convergence_checks,
    seed                                = seed,
    silent                              = silent,
    dots                                = list(...)
  )
}


#' @title Bayesian Multivariate and Multilevel PEESE Model
#'
#' @description
#' Fits Bayesian precision-effect estimate with standard errors (PEESE) models
#' with known sampling covariance and BayesTools random-effect formulas.
#'
#' @inheritParams brma.mv
#' @inheritParams bPEESE
#'
#' @details
#' `bPEESE.mv()` combines the PEESE publication-bias regression of [bPEESE()]
#' with the known-covariance and random-formula model of [brma.mv()]. The PEESE
#' predictor is the marginal sampling variance `diag(V)`; known off-diagonal
#' sampling covariance remains part of the joint likelihood.
#'
#' The supported known-`V` backends, random-effect marginalization rules,
#' predictive targets, and diagnostic semantics are the same as for
#' [brma.mv()]. In particular, `random` replaces the specialized `cluster`
#' argument used by [bPEESE()].
#'
#' @return A fitted object of class
#' `c("bPEESE.mv", "bPEESE", "brma.mv", "brma.norm", "brma")`.
#'
#' @examples \dontrun{
#' V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
#' fit <- bPEESE.mv(
#'   yi      = c(0.10, 0.20),
#'   V       = V,
#'   measure = "GEN",
#'   seed    = 1,
#'   silent  = TRUE
#' )
#' summary(fit)
#' }
#'
#' @seealso [bPEESE()], [bPET.mv()], [brma.mv()], [summary.brma()]
#'
#' @export
bPEESE.mv <- function(
    # input specification
    yi, V, ni,
    mods, scale, random,
    R = NULL, Rscale = "cor",
    data, slab, subset,
    measure,

    # prior specification
    prior_effect, prior_heterogeneity, prior_mods, prior_scale, prior_bias,
    standardize_continuous_predictors = TRUE,
    set_contrast_factor_predictors = "treatment",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,
    effect_direction = "detect",

    # MCMC fitting settings
    known_v_parameterization = "auto",
    known_v_residual_fraction = 0.10,
    marginalize_estimate_level = TRUE,
    sample = 5000, burnin = 2000, adapt = 500,
    chains = 3, thin = 1, parallel = FALSE,
    autofit = FALSE, autofit_control = set_autofit_control(),
    convergence_checks = set_convergence_checks(),

    # additional settings
    seed = NULL, silent, ...,
    vi = NULL, sei = NULL) {

  .fit_bias_regression_mv(
    matched_call_unevaluated            = match.call(expand.dots = FALSE),
    matched_call                        = match.call(),
    envir                               = parent.frame(),
    caller                              = "bPEESE.mv()",
    class                               = c(
      "bPEESE.mv", "bPEESE", "brma.mv", "brma.norm", "brma"
    ),
    bias_type                           = "PEESE",
    known_v_residual_fraction_specified = !missing(
      known_v_residual_fraction
    ),
    R                                   = R,
    Rscale                              = Rscale,
    measure                             = measure,
    prior_effect                        = prior_effect,
    prior_heterogeneity                 = prior_heterogeneity,
    prior_mods                          = prior_mods,
    prior_scale                         = prior_scale,
    prior_bias                          = prior_bias,
    standardize_continuous_predictors   = standardize_continuous_predictors,
    set_contrast_factor_predictors      = set_contrast_factor_predictors,
    prior_unit_information_sd           = prior_unit_information_sd,
    rescale_priors                      = rescale_priors,
    prior_informed_field                = prior_informed_field,
    prior_informed_subfield             = prior_informed_subfield,
    effect_direction                    = effect_direction,
    known_v_parameterization            = known_v_parameterization,
    known_v_residual_fraction           = known_v_residual_fraction,
    marginalize_estimate_level          = marginalize_estimate_level,
    sample                              = sample,
    burnin                              = burnin,
    adapt                               = adapt,
    chains                              = chains,
    thin                                = thin,
    parallel                            = parallel,
    autofit                             = autofit,
    autofit_control                     = autofit_control,
    convergence_checks                  = convergence_checks,
    seed                                = seed,
    silent                              = silent,
    dots                                = list(...)
  )
}


.fit_bias_regression_mv <- function(
    matched_call_unevaluated, matched_call, envir, caller, class, bias_type,
    known_v_residual_fraction_specified,
    R, Rscale, measure,
    prior_effect, prior_heterogeneity, prior_mods, prior_scale, prior_bias,
    standardize_continuous_predictors,
    set_contrast_factor_predictors,
    prior_unit_information_sd, rescale_priors,
    prior_informed_field, prior_informed_subfield, effect_direction,
    known_v_parameterization, known_v_residual_fraction,
    marginalize_estimate_level,
    sample, burnin, adapt, chains, thin, parallel,
    autofit, autofit_control, convergence_checks,
    seed, silent, dots) {

  .brma_mv_reject_unsupported_dots(
    matched_call_unevaluated,
    caller = caller
  )

  missing_measure <- missing(measure)
  if (missing_measure && !isTRUE(dots[["only_data"]])) {
    .stop_missing_measure(caller)
  }
  if (missing_measure) {
    measure <- "GEN"
  }
  known_v_parameterization <- match.arg(
    known_v_parameterization,
    c("auto", "latent", "whitened", "block_mvn")
  )
  random_group_covariance <- .brma_mv_make_group_covariance(
    R      = R,
    Rscale = Rscale
  )
  dots <- .validate_constructor_dots(
    dots   = dots,
    caller = caller
  )

  object <- .createObject(
    dots  = dots,
    class = class,
    chains = chains, adapt = adapt, burnin = burnin, sample = sample,
    thin = thin, autofit = autofit, parallel = parallel, silent = silent,
    seed = seed, autofit_control = autofit_control,
    convergence_checks = convergence_checks
  )

  object[["data"]] <- .check_and_list_data(
    .call                              = matched_call,
    .envir                             = envir,
    class                              = "mv",
    set_contrast_factor_predictors    = set_contrast_factor_predictors,
    standardize_continuous_predictors = standardize_continuous_predictors,
    effect_direction                  = effect_direction,
    measure                           = measure,
    random_group_covariance            = random_group_covariance,
    known_v_parameterization            = known_v_parameterization,
    known_v_residual_fraction           = known_v_residual_fraction,
    known_v_residual_fraction_specified =
      known_v_residual_fraction_specified
  )
  if (isTRUE(dots[["only_data"]])) {
    return(object)
  }

  object[["priors"]] <- .check_and_list_priors.brma(
    prior_effect              = prior_effect,
    prior_heterogeneity       = prior_heterogeneity,
    prior_mods                = prior_mods,
    prior_scale               = prior_scale,
    prior_bias                = prior_bias,
    rescale_priors            = rescale_priors,
    prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field      = prior_informed_field,
    prior_informed_subfield   = prior_informed_subfield,
    data                      = object[["data"]],
    bias_type                 = bias_type
  )

  .finalize_brma_mv_object(
    object                     = object,
    marginalize_estimate_level = marginalize_estimate_level,
    only_priors                = isTRUE(dots[["only_priors"]])
  )
}
