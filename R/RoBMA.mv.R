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
#' Omitted or NULL `random` specifies no heterogeneity or random-effect
#' inclusion mixture. Publication-bias and fixed-effect alternatives remain
#' available. Heterogeneity priors and scale formulas require an explicit
#' random structure; the known sampling covariance remains unchanged.
#'
#' `RoBMA.mv()` combines the publication-bias product space of [RoBMA()] with
#' the known-covariance, formula, random-effect, and prediction machinery of
#' [brma.mv()]. Selection models, PET, PEESE, and the unadjusted branch are
#' averaged together according to `model_type` or custom `prior_bias` and
#' `prior_bias_null` specifications. Selection is defined over the estimates
#' in each publication event. PET uses `sqrt(diag(V))` and PEESE uses `diag(V)` as their
#' bias predictors while the likelihood retains the full known sampling
#' covariance.
#'
#' Random structures supplied through `random` use the same independent
#' component gates as [BMA.mv()]. A Dirichlet prior allocates the slab total
#' variance across multiple top-level components and each allocation is
#' multiplied by its own inclusion indicator without renormalizing the
#' remaining allocations. The default prior inclusion probability is 0.5 for
#' every component. Public `tau_total` and `tau2_total` draws are the realized
#' gated aggregate, including the all-off zero branch. Public
#' `tau2_prop(...)` draws are realized shares conditional on positive total
#' heterogeneity; excluded components have zero share, and the all-off branch
#' is undefined. The positive slab scale and raw Dirichlet weights remain
#' internal coordinates.
#'
#' The constructor's `selection` specification applies to every generated
#' weightfunction. It integrates estimate-level random effects and the complete
#' sampling error by default, and conditions on other random effects.
#' Explicit priors carry their own [selection_model()]
#' settings, including publication grouping and product or best-p-value weights.
#' Active selection branches must share one conditioning cell and publication
#' partition. Their bins, weight priors, and weighting rules may differ.
#' Non-selection branches retain the corresponding Gaussian, PET, or PEESE
#' contribution with the same contextual source representation.
#'
#' Whole sampling-error selection, grouping, and the Gaussian source
#' likelihood follow [bselmodel.mv()]. Numerical backends do not change these
#' model choices. All three source choices are explicit in `selection`.
#' Without an active selection branch, supported estimate-level Gaussian
#' random intercepts are automatically integrated into the likelihood variance.
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
#'   selection = selection_model(group = study),
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
    set_contrast_factor_predictors = "meandif",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,
    model_type = "PSMA",

    # selection likelihood
    selection = BayesTools::selection_model(),
    selection_control = set_selection_likelihood_control(),

    # MCMC fitting settings
    known_v_parameterization = "auto",
    sample = 5000, burnin = 2000, adapt = 500,
    chains = 3, thin = 1, parallel = FALSE,
    autofit = FALSE, autofit_control = set_autofit_control(),
    convergence_checks = set_convergence_checks(),

    # additional settings
    seed = NULL, silent, ...,
    vi = NULL, sei = NULL) {

  BayesTools::check_selection_model(selection, name = "selection")
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
    R                                   = R,
    Rscale                              = Rscale,
    standardize_continuous_predictors   = standardize_continuous_predictors,
    set_contrast_factor_predictors      = set_contrast_factor_predictors,
    known_v_parameterization            = known_v_parameterization,
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
    random_component_averaging = TRUE,
    weightfunction_model       = selection
  )

  .finalize_mv_object(
    object                     = object,
    selection_control          = selection_control,
    only_priors                = isTRUE(dots[["only_priors"]])
  )
}
