#' @title Robust Bayesian Model-Averaged Multivariate Meta-Analysis
#'
#' @description Fits robust Bayesian model-averaged multivariate and
#' multilevel meta-analytic models with known sampling covariance and formula
#' random effects.
#'
#' @inheritParams brma.mv
#' @inheritParams RoBMA
#' @inheritParams bselmodel.mv
#'
#' @details
#' `RoBMA.mv()` combines the publication-bias product space of [RoBMA()] with
#' the known-covariance, formula, random-effect, and prediction machinery of
#' [brma.mv()]. Selection models, PET, PEESE, and the unadjusted branch are
#' averaged together according to `model_type` or custom `prior_bias` and
#' `prior_bias_null` specifications. Selection is always defined at the
#' estimate level. PET uses `sqrt(diag(V))` and PEESE uses `diag(V)` as their
#' bias predictors while the likelihood retains the full known sampling
#' covariance.
#'
#' Random structures supplied through `random` use the same independent
#' component gates as [BMA.mv()]. A Dirichlet prior allocates the slab total
#' variance across multiple top-level components and each allocation is
#' multiplied by its own inclusion indicator without renormalizing the
#' remaining allocations. The default prior inclusion probability is 0.5 for
#' every component.
#'
#' With `selection_likelihood = "exact"`, all Gaussian random effects are
#' analytically marginalized and each connected sampling/random covariance
#' block uses its joint selected-Gaussian density in selection-model branches.
#' Non-selection branches in the product space retain their ordinary Gaussian,
#' PET, or PEESE likelihood contribution. With
#' `selection_likelihood = "approximate"`, selection-model branches use the
#' row-wise selected-normal likelihood conditional on sampled random effects.
#' Correlated known sampling covariance then requires the latent known-`V`
#' representation. `known_v_parameterization = "auto"` is changed to that
#' representation only when the fitted bias mixture actually contains a
#' selection model; explicitly requesting `"whitened"` or `"block_mvn"` is
#' rejected for this target.
#'
#' `marginalize_estimate_level` applies to models without an exact selection
#' branch and to the approximate selection likelihood. Exact selection
#' necessarily marginalizes every Gaussian random-effect block into the joint
#' observation covariance.
#'
#' Product-space marginal likelihood and bridge-sampling methods are not
#' available. Predictive comparison through [loo.brma()] and [waic.brma()]
#' remains available. If `V` is singular, every allowed product-space branch
#' must contain structural variance sufficient to regularize its null space.
#'
#' @return A fitted object of class
#' `c("RoBMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")`.
#'
#' @examples \dontrun{
#' data("dat.assink2016", package = "metadat")
#' V <- metafor::vcalc(
#'   vi,
#'   cluster = study,
#'   type = deltype,
#'   obs = esid,
#'   rho = c(0.7, 0.5),
#'   data = dat.assink2016
#' )
#' fit <- RoBMA.mv(
#'   yi = yi,
#'   V = V,
#'   mods = ~ deltype,
#'   random = ~ 1 | study / esid,
#'   data = dat.assink2016,
#'   measure = "SMD",
#'   seed = 1,
#'   silent = TRUE
#' )
#' summary(fit)
#' summary_models(fit)
#' }
#'
#' @seealso [RoBMA()], [BMA.mv()], [brma.mv()], [bselmodel.mv()], [bPET.mv()],
#'   [bPEESE.mv()], [set_selection_likelihood_control()], [summary.brma()],
#'   [summary_models()]
#'
#' @export
RoBMA.mv <- function(
    # input specification
    yi, V, ni,
    mods, scale, random,
    R = NULL, Rscale = "cor",
    data, slab, subset,
    measure, effect_direction = "detect",

    # prior specification
    prior_effect, prior_heterogeneity, prior_mods, prior_scale, prior_bias,
    prior_effect_null, prior_heterogeneity_null,
    prior_mods_null, prior_scale_null, prior_bias_null,
    standardize_continuous_predictors = TRUE,
    set_contrast_factor_predictors = "treatment",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,
    model_type = "PSMA",

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

  selection_likelihood <- match.arg(selection_likelihood)
  initialized <- .initialize_mv_object(
    matched_call_unevaluated            = match.call(expand.dots = FALSE),
    matched_call                        = match.call(),
    envir                               = parent.frame(),
    caller                              = "RoBMA.mv()",
    object_class                        = c(
      "RoBMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma"
    ),
    dots                                = list(...),
    missing_measure                     = missing(measure),
    measure                             = measure,
    known_v_residual_fraction_specified = !missing(
      known_v_residual_fraction
    ),
    R                                   = R,
    Rscale                              = Rscale,
    standardize_continuous_predictors   = standardize_continuous_predictors,
    set_contrast_factor_predictors      = set_contrast_factor_predictors,
    known_v_parameterization            = known_v_parameterization,
    known_v_residual_fraction           = known_v_residual_fraction,
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
    effect_direction                    = effect_direction
  )
  object <- initialized[["object"]]
  dots   <- initialized[["dots"]]
  if (isTRUE(dots[["only_data"]])) {
    return(object)
  }

  object[["priors"]] <- .check_and_list_priors.RoBMA(
    prior_effect               = prior_effect,
    prior_heterogeneity        = prior_heterogeneity,
    prior_mods                 = prior_mods,
    prior_scale                = prior_scale,
    prior_bias                 = prior_bias,
    prior_effect_null          = prior_effect_null,
    prior_heterogeneity_null   = prior_heterogeneity_null,
    prior_mods_null            = prior_mods_null,
    prior_scale_null           = prior_scale_null,
    prior_bias_null            = prior_bias_null,
    rescale_priors             = rescale_priors,
    prior_unit_information_sd  = prior_unit_information_sd,
    prior_informed_field       = prior_informed_field,
    prior_informed_subfield    = prior_informed_subfield,
    data                       = object[["data"]],
    model_type                 = model_type,
    random_component_averaging = TRUE
  )

  .finalize_mv_object(
    object                     = object,
    selection_likelihood       = selection_likelihood,
    selection_control          = selection_control,
    marginalize_estimate_level = marginalize_estimate_level,
    only_priors                = isTRUE(dots[["only_priors"]])
  )
}
