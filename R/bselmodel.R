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
#' @param selection specification created by [selection_model()] for
#' automatically constructed weightfunction priors. The default integrates
#' estimate-level random effects and the complete sampling error, and conditions
#' on other random effects. Each source can instead be conditioned upon or
#' integrated through the corresponding selection-model setting. The default
#' `weight_rule = "product"` needs no publication groups and ignores `group`.
#' Publication grouping is active only for `weight_rule = "best"`. Explicit
#' `prior_bias` objects retain their own selection specification; this argument
#' does not overwrite it. Conditioning on sampling variation retains the entire
#' sampling-error realization, including in ordinary univariate models.
#' This sampling setting requires omitted or unit observation `weights`.
#' @param selection_control numerical integration settings created by
#' [set_selection_likelihood_control()]. These affect numerical evaluation,
#' without changing the conditioning model or the vector weighting rule.
#'
#' @details
#' `bselmodel()` is a normal/effect-size selection-model constructor. Custom
#' `prior_bias` can be a weightfunction prior or a supported BayesTools
#' selection-kernel prior; p-hacking kernels are not supported in active RoBMA.
#' By default, the model conditions on contextual cluster effects and uses
#' the product of estimate weights. For independent rows conditional on these
#' effects, this gives the usual selected-normal likelihood with within-cluster
#' heterogeneity integrated. Use `selection = selection_model(...)` to choose
#' `weight_rule = "best"`, which uses the weight at the smallest p-value in
#' each publication group. For `"best"`, `group` identifies a data column;
#' when `group = NULL`, the specialized `cluster` argument supplies publication
#' groups, or an unclustered model uses one group per estimate. Product models
#' use neither input to define selection groups.
#'
#' All conditioning choices use the same Gaussian source model. Conditioned
#' sources remain latent; their population distributions stay outside selection
#' normalization. Integrated sources may be reweighted by selection. With all
#' sources conditioned and positive weights, the weights cancel and the
#' observed law is the ordinary Gaussian model. See [bselmodel.mv()] for the
#' complete source and covariance contract.
#'
#' Product normalizers factorize only for conditionally independent
#' rows. Dependent events use supported covariance-factor quadrature or fixed
#' randomized quasi-Monte Carlo integration with explicit error diagnostics.
#' Non-unit observation `weights` require
#' `known_sampling_variance = "integrate"` and independent product factors.
#' Unit weights are equivalent to omitting `weights`.
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

  # selection likelihood
  selection = BayesTools::selection_model(),
  selection_control = set_selection_likelihood_control(),

  # MCMC fitting settings
  sample = 5000, burnin = 2000, adapt = 500,
  chains = 3, thin = 1, parallel = FALSE,
  autofit = FALSE, autofit_control = set_autofit_control(),
  convergence_checks = set_convergence_checks(),

  # additional settings
  seed = NULL, silent, ...
) {

  BayesTools::check_selection_model(selection, name = "selection")

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
    effect_direction = effect_direction, measure = measure,
    selection_binding = !isTRUE(dots[["only_data"]]))
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
    data = object[["data"]], bias_type = "selmodel", steps = steps,
    weightfunction_model = selection)
  object <- .prepare_selection_model_object(object)
  object <- .prepare_selection_likelihood_object(
    object            = object,
    selection_control = selection_control
  )
  .fit_and_finalize_object(
    object,
    only_priors = isTRUE(dots[["only_priors"]])
  )
}
