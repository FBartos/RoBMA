#' @title Bayesian Model-Averaged Multivariate Meta-Analysis
#'
#' @description Fits Bayesian model-averaged multivariate and multilevel
#' meta-analytic models with known sampling covariance and formula random
#' effects. Random structures supplied through `random` are model-averaged with
#' independent component gates.
#'
#' @inheritParams brma.mv
#' @inheritParams RoBMA_prior_specification
#'
#' @details
#' Omitted or NULL `random` specifies no heterogeneity and no random-effect
#' inclusion mixture. Heterogeneity priors and scale formulas require an
#' explicit random structure; the known sampling covariance remains unchanged.
#'
#' `BMA.mv()` combines the product-space model-averaging workflow of [BMA()]
#' with the likelihood, formula, random-effect, and known-covariance machinery
#' of [brma.mv()]. For multiple top-level random components, a Dirichlet prior
#' allocates the slab total variance and an independent Bernoulli gate multiplies
#' each allocated contribution. The allocation is not renormalized when a gate
#' is off. Thus component \eqn{j} contributes
#' \deqn{I_j w_j \tau^2,}
#' where \eqn{I_j} is its inclusion indicator, \eqn{w_j} is its slab allocation,
#' and \eqn{\tau} is the slab total SD. With one component, no artificial
#' Dirichlet weight is created and the contribution is \eqn{I_1 \tau^2}.
#'
#' By default, every top-level component has an independent 0.5 inclusion
#' probability. `prior_heterogeneity` defines the shared positive slab and
#' `prior_heterogeneity_null` defines exclusion at exactly zero. Their prior
#' weights determine the common gate probability; multiple alternative priors
#' retain their relative weights within the slab. Setting
#' `prior_heterogeneity_null = NULL` or `FALSE` fixes all component gates on,
#' while `prior_heterogeneity = NULL` or `FALSE` fixes them off. Nested random
#' terms are gated as one top-level component before their internal Dirichlet
#' split.
#'
#' A partial [BayesTools::prior_random()] supplied through
#' `prior_heterogeneity` may override contrasts, correlations, monitoring, and
#' parameterization. `BMA.mv()` owns the gated scale architecture, so a custom
#' `prior_random()` must not supply an SD, SD source, term-specific SD,
#' covariance-owned SD, or variance allocation.
#'
#' Random-effect inclusion probabilities are reported separately from fixed
#' effect and scale-regression inclusion. Averaged random-component SDs retain
#' the excluded zero branch; conditional random summaries condition each
#' component on its own gate. `tau_total` and `tau2_total` are the realized gated
#' aggregate, including zero when every component is excluded. `tau2_prop(j)` is
#' the realized share \eqn{I_j w_j / \sum_k I_k w_k}, conditional on positive
#' total heterogeneity. Excluded components therefore have zero share, and the
#' all-off branch is omitted only from variance-proportion summaries. The
#' positive slab scale and raw Dirichlet weights remain internal coordinates.
#'
#' Product-space marginal likelihood and bridge-sampling methods are not
#' available. Predictive comparison through [loo.brma()] and [waic.brma()]
#' remains available. If `V` is singular, every allowed product-space branch
#' must contain structural variance sufficient to regularize its null space;
#' an allowed all-components-off branch is therefore rejected.
#'
#' The function does not add publication-bias models. Use [RoBMA.mv()] to add
#' model-averaged selection, PET, and PEESE publication-bias adjustments.
#'
#' @return A fitted object of class
#' `c("BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma")`.
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
#' fit <- BMA.mv(
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
#' @seealso [BMA()], [RoBMA.mv()], [brma.mv()], [summary.brma()],
#'   [summary_models()]
#'
#' @export
BMA.mv <- function(
    # input specification
    yi, V, ni,
    mods, scale, random,
    R = NULL, Rscale = "cor",
    data, slab, subset,
    measure,

    # prior specification
    prior_effect, prior_heterogeneity, prior_mods, prior_scale,
    prior_effect_null, prior_heterogeneity_null,
    prior_mods_null, prior_scale_null,
    standardize_continuous_predictors = TRUE,
    set_contrast_factor_predictors = "meandif",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,

    # MCMC fitting settings
    known_v_parameterization = "auto",
    sample = 5000, burnin = 2000, adapt = 500,
    chains = 3, thin = 1, parallel = FALSE,
    autofit = FALSE, autofit_control = set_autofit_control(),
    convergence_checks = set_convergence_checks(),

    # additional settings
    seed = NULL, silent, ...,
    vi = NULL, sei = NULL) {

  initialized <- .initialize_mv_object(
    matched_call_unevaluated            = match.call(expand.dots = FALSE),
    matched_call                        = match.call(),
    envir                               = parent.frame(),
    caller                              = "BMA.mv()",
    object_class                        = c(
      "BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma"
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
    silent                              = silent
  )
  object <- initialized[["object"]]
  dots   <- initialized[["dots"]]
  if (isTRUE(dots[["only_data"]])) {
    return(object)
  }

  object$priors <- .check_and_list_priors.RoBMA(
    prior_effect               = prior_effect,
    prior_heterogeneity        = prior_heterogeneity,
    prior_mods                 = prior_mods,
    prior_scale                = prior_scale,
    prior_effect_null          = prior_effect_null,
    prior_heterogeneity_null   = prior_heterogeneity_null,
    prior_mods_null            = prior_mods_null,
    prior_scale_null           = prior_scale_null,
    rescale_priors             = rescale_priors,
    prior_unit_information_sd  = prior_unit_information_sd,
    prior_informed_field       = prior_informed_field,
    prior_informed_subfield    = prior_informed_subfield,
    data                        = object[["data"]],
    random_component_averaging = TRUE
  )

  .finalize_mv_object(
    object                     = object,
    only_priors                = isTRUE(dots[["only_priors"]])
  )
}
