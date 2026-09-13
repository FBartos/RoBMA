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
#' Omitted or NULL `random` specifies a fixed-effect selection model with
#' no implicit heterogeneity. Both random-effect selection axes are then
#' inapplicable; the sampling-error selection axis still applies.
#' Heterogeneity priors and scale formulas require an explicit random structure.
#'
#' `bselmodel.mv()` combines the selection-model interface of [bselmodel()] with
#' the known-covariance and random-formula model of [brma.mv()]. Each
#' weightfunction carries a [selection_model()] specifying whether estimate-level
#' random effects, other random effects, and the complete sampling error are
#' conditioned upon or integrated before selection normalization. The default
#' integrates estimate-level random effects and sampling error, conditions on
#' other random effects, and uses the product of estimate weights.
#' `weight_rule = "best"` instead uses
#' the weight at the smallest actual p-value in the publication group.
#'
#' With one cutoff and all three sources integrated, `"best"`
#' corresponds to the relaxed report-all rule of
#' \insertCite{vanaertinpressmultivariate;textual}{RoBMA}: a favorable result
#' can select the publication vector, after which all measured outcomes are
#' reported. The other conditioning cells and multiple-cutoff best rules are
#' extensions of that specification.
#'
#' A publication group is resolved from an explicit `group` column in the
#' prior or a supported constructor cluster. Multivariate random structures
#' require an explicit group column. Covariance matrices, including standard
#' [metafor::vcalc()] results, and random-effect grouping factors do not
#' establish publication identities.
#'
#' Sampling error is \eqn{e \sim N(0,V)}. With positive sampling standard
#' errors it can equivalently be written as \eqn{e = Sz}, with
#' \eqn{z \sim N(0,R)} and \eqn{V = SRS}. Here \eqn{S} contains the original
#' sampling standard errors and \eqn{R} preserves their complete correlation
#' structure. Conditioning retains the entire realized vector \eqn{e} during
#' selection normalization. No residual sampling variation remains inside that
#' normalizer. Integration instead contributes the full \eqn{V} to candidate
#' variation. These choices apply equally to univariate, diagonal, and
#' correlated sampling covariance, including `vi` and `sei` inputs.
#'
#' The covariance `V` remains authoritative. A diagonal-plus-factor input from
#' [known_v_factor()] can support more efficient calculations,
#' but equivalent representations define the same selection model when their
#' publication groups agree. Selection thresholds always use the original
#' `sqrt(diag(V))`.
#'
#' The likelihood and post-fit methods use the same resolved conditioning
#' model. `known_v_parameterization` and numerical integration settings cannot
#' change it. The `estimate_random_effects` setting controls estimate-level
#' true-effect variation. Shared covariance-factor and dense numerical paths use
#' explicit integration diagnostics.
#'
#' The source settings define the reporting model. An optional comparison with
#' integrating retained contexts can be requested with
#' [selection_sensitivity_diagnostics()] or stored with
#' [add_selection_sensitivity_diagnostics()]. Fitting does not compute this
#' separate comparison automatically.
#'
#' `estimate_random_effects`, `other_random_effects`, and
#' `known_sampling_variance` each accept `"condition"` or `"integrate"`.
#' Both `"product"` and `"best"` support the resulting choices for supported
#' Gaussian sources contained within publication groups.
#'
#' Retained contexts remain unknown and are estimated. Their population mixing
#' laws stay outside the selection normalizer; integrating a context instead
#' allows selection to reweight its mixing law. An absent context makes its
#' two settings coincide. BayesTools identifies an estimate-level random term
#' by a distinct grouping level for every retained estimate. At most one such
#' term is allowed. Every other declared term is an other-level source,
#' including terms with independent coefficient supports inside repeated groups.
#' These roles do not split coefficient families or discard known group
#' covariance; each choice applies to the complete declared random term.
#'
#' With every outcome-generating source conditioned upon, candidate outcomes
#' are deterministic within a retained context. If weights are positive almost
#' surely, they cancel under conditional normalization: the observed law is
#' the ordinary Gaussian model and carries no information about the selection
#' weights. A zero acceptance probability on a positive-probability set of
#' retained contexts makes this selection process invalid. It must not be
#' handled by silently discarding those contexts. Ordinary population-level
#' publication selection instead integrates all outcome-generating sources.
#'
#' Product weights also support the existing Gaussian dependency paths spanning
#' publication groups. Best weights require integrated dependencies to remain
#' within each publication group; an unsupported source is identified before
#' fitting. Conditioned sources may connect groups. Integrated candidate
#' covariance can be singular while the observed law remains nondegenerate
#' after averaging over retained sources. No residual noise is added to change
#' this boundary. Ensembles require a common conditioning cell and publication
#' partition across active selection branches, while allowing different bins,
#' weight priors, and supported weight rules.
#' Non-unit observation `weights` require
#' `known_sampling_variance = "integrate"` and independent product factors.
#' Unit weights are equivalent to omitting `weights`.
#'
#' Joint likelihoods, single-model bridge evidence, conditional-density
#' estimation, selected prediction, latent summaries, and z-plots use this
#' resolved model. Estimate-unit LOO, WAIC, and LOO-PIT integrate deleted
#' outcomes within the original publication event. Their conditional row scores
#' must not be summed and interpreted as a joint group likelihood. Cluster or
#' block deletion remains unavailable for arbitrary `brma.mv()` random formulas.
#' Product-space ensembles do not provide bridge marginal likelihoods.
#' With conditioned sampling, estimate deletion integrates the deleted sampling
#' error conditional on the retained sampling errors of the remaining estimates.
#'
#' Marginal selected predictions draw new retained contexts from their original
#' mixing laws before drawing selected outcomes and reconstructing latent true
#' effects. At estimate depth, latent predictions draw fitted true effects;
#' [blup()] and [fitted.brma()] instead return conditional means. Response
#' prediction at estimate depth describes a new reporting event around fitted
#' truth, with new sampling error. Prediction's `conditioning_depth` is a
#' separate choice from the selection-model source settings.
#' Cluster depth is available only for the
#' specialized `cluster` interface. Known-covariance latent and response
#' predictions with `newdata` require `V_new` and the resolved publication-group
#' reference; cross-covariance with fitted outcomes and non-marginal `newdata`
#' are unavailable. See [predict.brma()] for prediction targets and labels.
#'
#' When the complete sampling error is conditioned upon, total fitted truth
#' equals the observed estimate minus its fitted sampling error. Conditional
#' means and estimate-depth latent draws then coincide within each posterior
#' draw; uncertainty remains across posterior draws. Blockwise [ranef()] retains
#' its conditional-mean target for each random-effect block.
#'
#' Partial-vector calculations retain the full selection event. Densities need
#' a nondegenerate observed Gaussian law, except for supported one-coordinate
#' projections of explicitly declared rank-one sources. Multicoordinate singular
#' Lebesgue densities and general singular conditional latent posteriors are
#' unavailable; zero integrated true-effect covariance gives deterministic
#' latent effects. No general public multivariate CDF is provided. LOO-PIT and
#' z-plots use their respective deletion and marginal projection targets.
#'
#' @return A fitted object of class
#' `c("bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma")`.
#'
#' @examples \dontrun{
#' dat <- data.frame(yi = c(0.10, 0.20), paper = c("A", "A"))
#' V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
#' fit <- bselmodel.mv(
#'   yi      = yi,
#'   V       = V,
#'   data    = dat,
#'   measure = "GEN",
#'   prior_bias = prior_weightfunction(
#'     "one-sided", 0.025, model = selection_model(group = paper)
#'   ),
#'   seed    = 1,
#'   silent  = TRUE
#' )
#' summary(fit)
#' }
#'
#' @seealso [bselmodel()], [brma.mv()], [set_selection_likelihood_control()],
#'   [summary.brma()]
#'
#' @references \insertAllCited{}
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
    caller                              = "bselmodel.mv()",
    object_class                        = c(
      "bselmodel.mv", "bselmodel", "brma.mv", "brma.norm", "brma"
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
    steps                     = steps,
    weightfunction_model      = selection
  )

  .finalize_mv_object(
    object                     = object,
    selection_control          = selection_control,
    only_priors                = isTRUE(dots[["only_priors"]])
  )
}
