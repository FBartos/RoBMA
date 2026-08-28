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
#' component on its own gate. `sd_total` and `var_prop(...)` describe the
#' pre-gate slab allocation.
#'
#' Product-space marginal likelihood and bridge-sampling methods are not
#' available. Predictive comparison through [loo.brma()] and [waic.brma()]
#' remains available. If `V` is singular, every allowed product-space branch
#' must contain structural variance sufficient to regularize its null space;
#' an allowed all-components-off branch is therefore rejected.
#'
#' The function does not add publication-bias models. A multivariate extension
#' of [RoBMA()] is not implemented.
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
#' @seealso [BMA()], [brma.mv()], [summary.brma()], [summary_models()]
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
    set_contrast_factor_predictors = "treatment",
    prior_unit_information_sd, rescale_priors = 1,
    prior_informed_field, prior_informed_subfield,

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
    caller = "BMA.mv()"
  )

  dots                                <- list(...)
  missing_measure                     <- missing(measure)
  known_v_residual_fraction_specified <- !missing(known_v_residual_fraction)
  matched_call                        <- match.call()

  if (missing_measure && !isTRUE(dots[["only_data"]])) {
    .stop_missing_measure("BMA.mv()")
  }
  if (missing_measure) {
    measure <- "GEN"
  }
  random_group_covariance <- .brma_mv_make_group_covariance(
    R      = R,
    Rscale = Rscale
  )
  known_v_parameterization <- match.arg(
    known_v_parameterization,
    c("auto", "latent", "whitened", "block_mvn")
  )
  dots <- .validate_constructor_dots(
    dots   = dots,
    caller = "BMA.mv()"
  )

  object <- .createObject(
    dots = dots,
    class = c("BMA.mv", "RoBMA", "brma.mv", "brma.norm", "brma"),
    chains = chains, adapt = adapt, burnin = burnin, sample = sample,
    thin = thin, autofit = autofit, parallel = parallel, silent = silent,
    seed = seed, autofit_control = autofit_control,
    convergence_checks = convergence_checks
  )

  object$data <- .check_and_list_data(
    .call                              = matched_call,
    .envir                             = parent.frame(),
    class                              = "mv",
    set_contrast_factor_predictors    = set_contrast_factor_predictors,
    standardize_continuous_predictors = standardize_continuous_predictors,
    measure                            = measure,
    random_group_covariance            = random_group_covariance,
    known_v_parameterization            = known_v_parameterization,
    known_v_residual_fraction           = known_v_residual_fraction,
    known_v_residual_fraction_specified = known_v_residual_fraction_specified
  )
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

  .finalize_brma_mv_object(
    object                     = object,
    marginalize_estimate_level = marginalize_estimate_level,
    only_priors                = isTRUE(dots[["only_priors"]])
  )
}
