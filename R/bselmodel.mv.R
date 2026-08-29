#' @title Bayesian Multivariate and Multilevel Selection Model
#'
#' @description
#' Fits estimate-level Bayesian selection models with known sampling covariance
#' and BayesTools random-effect formulas.
#'
#' @inheritParams brma.mv
#' @inheritParams bselmodel
#'
#' @details
#' `bselmodel.mv()` combines the selection-model interface of [bselmodel()] with
#' the known-covariance and random-formula model of [brma.mv()]. Selection is
#' always defined at the estimate level, so the selection weight is the product
#' of the row-specific weights.
#'
#' With `selection_likelihood = "exact"`, all Gaussian random effects are
#' analytically marginalized and each connected sampling/random covariance
#' block is fitted with its joint selected-Gaussian density. This target is
#' independent of `known_v_parameterization`; that argument remains relevant to
#' the approximate likelihood and downstream known-`V` machinery. With
#' `selection_likelihood = "approximate"`, the existing row-wise selected-normal
#' likelihood is used. The automatic known-`V` backend resolves to `"latent"`
#' for this approximate target; explicitly requesting `"whitened"` or
#' `"block_mvn"` is rejected because those transformations do not preserve
#' estimate-level selection events.
#'
#' `marginalize_estimate_level` applies only to the approximate likelihood. The
#' exact likelihood necessarily marginalizes every Gaussian random-effect block
#' into the joint observation covariance.
#'
#' @return A fitted object of class
#' `c("bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma")`.
#'
#' @examples \dontrun{
#' V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
#' fit <- bselmodel.mv(
#'   yi      = c(0.10, 0.20),
#'   V       = V,
#'   measure = "GEN",
#'   steps   = 0.025,
#'   seed    = 1,
#'   silent  = TRUE
#' )
#' summary(fit)
#' }
#'
#' @seealso [bselmodel()], [brma.mv()], [set_selection_likelihood_control()],
#'   [summary.brma()]
#'
#' @export
bselmodel.mv <- function(
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
    effect_direction = "detect", steps,

    # selection likelihood
    selection_likelihood = c("exact", "approximate"),
    selection_control = set_selection_likelihood_control(),

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

  matched_call_unevaluated            <- match.call(expand.dots = FALSE)
  .brma_mv_reject_unsupported_dots(
    matched_call_unevaluated,
    caller = "bselmodel.mv()"
  )

  dots                                <- list(...)
  missing_measure                     <- missing(measure)
  known_v_residual_fraction_specified <- !missing(known_v_residual_fraction)
  selection_likelihood                <- match.arg(selection_likelihood)
  known_v_parameterization <- match.arg(
    known_v_parameterization,
    c("auto", "latent", "whitened", "block_mvn")
  )
  if (identical(selection_likelihood, "approximate")) {
    if (identical(known_v_parameterization, "auto")) {
      known_v_parameterization <- "latent"
    } else if (!identical(known_v_parameterization, "latent")) {
      stop(
        "'selection_likelihood = \"approximate\"' requires ",
        "'known_v_parameterization = \"latent\"'.",
        call. = FALSE
      )
    }
  }

  if (missing_measure && !isTRUE(dots[["only_data"]])) {
    .stop_missing_measure("bselmodel.mv()")
  }
  if (missing_measure) {
    measure <- "GEN"
  }
  random_group_covariance <- .brma_mv_make_group_covariance(
    R      = R,
    Rscale = Rscale
  )
  dots <- .validate_constructor_dots(
    dots   = dots,
    caller = "bselmodel.mv()"
  )

  object <- .createObject(
    dots  = dots,
    class = c("bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma"),
    chains = chains, adapt = adapt, burnin = burnin, sample = sample,
    thin = thin, autofit = autofit, parallel = parallel, silent = silent,
    seed = seed, autofit_control = autofit_control,
    convergence_checks = convergence_checks
  )

  matched_call <- match.call()
  matched_call[["known_v_parameterization"]] <- known_v_parameterization
  object[["data"]] <- .check_and_list_data(
    .call                              = matched_call,
    .envir                             = parent.frame(),
    class                              = "mv",
    set_contrast_factor_predictors    = set_contrast_factor_predictors,
    standardize_continuous_predictors = standardize_continuous_predictors,
    effect_direction                  = effect_direction,
    measure                           = measure,
    random_group_covariance            = random_group_covariance,
    known_v_parameterization            = known_v_parameterization,
    known_v_residual_fraction           = known_v_residual_fraction,
    known_v_residual_fraction_specified = known_v_residual_fraction_specified
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
    bias_type                 = "selmodel",
    steps                     = steps
  )

  .finalize_bselmodel_object(
    object                     = object,
    selection_likelihood       = selection_likelihood,
    selection_control          = selection_control,
    only_priors                = isTRUE(dots[["only_priors"]]),
    marginalize_estimate_level = marginalize_estimate_level
  )
}
