#' @title Bayesian Meta-Analysis With Known Sampling Covariance
#'
#' @description Fits normal-likelihood Bayesian meta-analytic models with a
#' known working variance-covariance matrix `V`. The interface supports
#' location moderators, BayesTools random-effect formulas using
#' [random-effect formula structure tags][random_effect_formula_tags], and
#' `brma.mv()`-specific known-`V` backends.
#'
#' @inheritParams brma
#' @inheritParams data_input
#' @inheritParams prior_specification
#' @inheritParams fitting_specification
#' @param R optional known random-effect group covariance/correlation matrix,
#' or a named list of matrices. Matrix row and column names identify grouping
#' levels. This is distinct from the known sampling covariance `V`.
#' @param Rscale scaling applied to each supplied `R` matrix before fitting.
#' Supported values are `"cor"`, `"none"`, `"cor0"`, `"cov0"`, with
#' metafor-compatible aliases `TRUE`, `FALSE`, `0`, `1`, `2`, and `3`.
#' @param marginalize_estimate_level logical; should a one-to-one estimate-level
#' random-intercept block be integrated into the likelihood variance instead of
#' sampled as latent random coefficients? Defaults to `TRUE`.
#'
#' @details
#' `brma.mv()` represents known sampling dependence through a latent `D + BB'`
#' decomposition of `V` (`known_v_parameterization = "latent"`), through an
#' eigen-rotated likelihood (`known_v_parameterization = "whitened"`), or
#' through an exact native block multivariate-normal likelihood
#' (`known_v_parameterization = "block_mvn"`). The default
#' `known_v_parameterization = "auto"` chooses the whitened backend whenever the
#' extra diagonal variance is scalar, otherwise chooses block-MVN when the
#' largest dense covariance block is at or below the configured feasibility
#' threshold, and falls back to the latent backend for large blocks. The automatic
#' block-MVN threshold defaults to 128 rows per covariance block and can be
#' changed with `options(RoBMA.known_v_block_mvn_max_block_size = value)`, where
#' `value` is a positive integer or `Inf`.
#'
#' The `random` argument uses BayesTools
#' [random-effect formula structure tags][random_effect_formula_tags].
#' Plain shorthand is accepted only for random intercepts such as
#' `~ 1 | study` or nested intercepts such as `~ 1 | study/effect`; random
#' slopes require an explicit structure tag or the `||` diagonal shorthand.
#' The optional `R` argument supplies known covariance or correlation matrices
#' across random-effect grouping levels, following `metafor::rma.mv()` naming.
#' `R` is separate from the known sampling covariance `V`: `V` describes
#' sampling-error dependence between estimates, while `R` describes the
#' latent group-axis covariance of random intercepts. Current `R` support is
#' limited to random-intercept blocks. `Rscale` controls how each supplied
#' `R` matrix is scaled before fitting and accepts
#' `"cor"`, `"none"`, `"cor0"`, `"cov0"`, plus metafor-compatible aliases
#' `TRUE`/`1` for `"cor"`, `FALSE`/`0` for `"none"`, `2` for `"cor0"`,
#' and `3` for `"cov0"`.
#' When `random` is used, shared sampled
#' random-effect blocks enter the conditional mean through BayesTools. With
#' `marginalize_estimate_level = TRUE`, a single one-to-one random-intercept
#' block is integrated into the likelihood variance; otherwise the likelihood
#' remains conditional on sampled random effects and uses only the known
#' sampling covariance `V`. With a single top-level random component, `scale`
#' models the row-wise total random-effect SD consumed by that component. With
#' multiple top-level random components, use a named list of `scale` formulas
#' whose names uniquely match the random components. Component names can be
#' ordinary text; RoBMA normalizes them internally for JAGS. Each formula
#' models the row-wise total SD for its component. A single `scale` formula is
#' rejected for multiple top-level random components because the target is
#' ambiguous.
#' When `scale` is a named list, `prior_scale` must either be omitted or use the
#' same named-list shape with one prior list per scale component, using the
#' same user-facing names. Post-fit scale predictions for random-formula
#' models keep these component scales separate rather than collapsing them into
#' one total \eqn{\tau}.
#' If `marginalize_estimate_level = TRUE`, a single random-intercept block that
#' maps one-to-one to estimates is compiled as a marginalized variance component
#' instead of sampled latent coefficients. This preserves estimate-level
#' heterogeneity semantics while keeping shared higher-level random effects
#' conditional. For known `R`, this automatic marginalization is available only
#' when BayesTools validates a one-column random-intercept block whose
#' `Z K Z'` row-space contribution is diagonal and one-to-one with the fitted
#' rows; RoBMA then adds the BayesTools-prepared row multiplier to the known-`V`
#' extra variance. Other known-`R` blocks remain sampled or are rejected by
#' BayesTools validation. Ambiguous cases with multiple one-to-one blocks are
#' rejected.
#'
#' The implementation intentionally omits the `weights` and `cluster` arguments
#' from `metafor::rma.mv()`. Use `random` for multilevel structures.
#' Estimate-unit `logLik()`, `add_loo()`, `add_waic()`, `residuals()`, and
#' `rstudent()` are available for known-`V` models. Bridge-sampling marginal
#' likelihoods via `add_marglik()` are available for fitted single-model
#' known-`V` objects. `rstandard()` is available for known-`V` models at
#' supported estimate and marginal conditioning depths.
#'
#' @section Diagnostic target semantics:
#' The main `brma.mv()` post-fit methods intentionally use different targets:
#' \itemize{
#'   \item `logLik(unit = "estimate")`, `add_loo()`/`loo()`,
#'     `add_waic()`/`waic()`, and `rstudent()` use estimate-unit log-score
#'     targets. When known `V` induces dependence, these contributions are
#'     Schur-complement conditionals \eqn{p(y_i \mid y_{-i}, \theta)} within
#'     dependency blocks. They are existing-estimate diagnostics, not
#'     independent new-estimate prediction targets.
#'     When known random-effect covariance `R` is supplied, sampled random
#'     effects remain conditioned in these estimate-depth targets. The known
#'     `R` matrix shapes the posterior and prior for those random effects, but
#'     it is not added again as a marginal \eqn{ZGZ'} covariance term. Supported
#'     known-`R` blocks compiled as marginalized instead enter through the
#'     diagonal extra variance term \eqn{\tau^2 r_i} prepared by BayesTools.
#'   \item `dfbetas()` and `covratio()` use estimate-unit PSIS weights. With
#'     correlated known `V`, their interpretation is parameter influence under
#'     conditional estimate deletion, not independent-study deletion.
#'   \item `add_marglik()`/`bridge_sampler()` uses the full joint fitted
#'     likelihood corresponding to the selected known-`V` backend. It does not
#'     use the estimate-wise conditional target from LOO/WAIC. Known `R`
#'     enters this full joint target either through the latent random-effect
#'     prior for sampled blocks or through the diagonal marginalized
#'     known-`R` variance for supported marginalized blocks.
#'   \item `hatvalues()`, marginal `rstandard()`, and `vif()` use a marginal
#'     GLS covariance target based on `V + ZGZ'`, where formula random effects
#'     are marginalized through the BayesTools covariance metadata.
#'   \item `dffits()` and `cooks.distance()` use fixed-location fitted-value
#'     influence targets. Cluster/block LOO is deferred for `brma.mv()` because
#'     its deletion target is a joint dependency block, not an independent row.
#' }
#' Default scalar heterogeneity summaries are not reported for `brma.mv()`
#' because known-`V` and random-formula structures need component-specific
#' summaries rather than a single \eqn{\tau}/\eqn{I^2} table.
#'
#' @return A fitted object of class `c("brma.mv", "brma.norm", "brma")`.
#'
#' @examples \dontrun{
#' V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
#' fit <- brma.mv(
#'   yi                        = c(0.10, 0.20),
#'   V                         = V,
#'   measure                   = "GEN",
#'   prior_unit_information_sd = 1,
#'   seed                      = 1,
#'   silent                    = TRUE
#' )
#' summary(fit)
#' }
#'
#' @seealso \code{\link{random_effect_formula_tags}},
#'   \code{\link{random_effect_prior_specification}}, [brma()],
#'   [summary.brma()], [predict.brma()]
#'
#' @export
brma.mv <- function(
    # input specification
    yi, V, ni,
    mods, scale, random,
    R = NULL, Rscale = "cor",
    data, slab, subset,
    measure,

    # prior specification
    prior_effect, prior_heterogeneity, prior_mods, prior_scale,
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
    seed = NULL, silent, ...
    ) {

  matched_call_unevaluated            <- match.call(expand.dots = FALSE)
  .brma_mv_reject_unsupported_dots(matched_call_unevaluated)

  dots                                <- list(...)
  missing_measure                     <- missing(measure)
  known_v_residual_fraction_specified <- !missing(known_v_residual_fraction)
  matched_call                        <- match.call()

  if (missing_measure && !isTRUE(dots[["only_data"]])) {
    .stop_missing_measure("brma.mv()")
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
    dots          = dots,
    caller        = "brma.mv()",
    allowed_extra = c("vi", "sei")
  )

  object <- .createObject(
    dots = dots, class = c("brma.mv", "brma.norm", "brma"),
    chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin,
    autofit = autofit, parallel = parallel, silent = silent, seed = seed,
    autofit_control = autofit_control, convergence_checks = convergence_checks
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
  if (isTRUE(dots[["only_data"]]))
    return(object)

  object$priors <- .check_and_list_priors.brma(
    prior_effect                   = prior_effect,
    prior_heterogeneity            = prior_heterogeneity,
    prior_mods                     = prior_mods,
    prior_scale                    = prior_scale,
    rescale_priors                 = rescale_priors,
    prior_unit_information_sd      = prior_unit_information_sd,
    prior_informed_field           = prior_informed_field,
    prior_informed_subfield        = prior_informed_subfield,
    data                           = object[["data"]])

  object <- .prepare_brma_mv_random_effects_compile(
    object                     = object,
    marginalize_estimate_level = marginalize_estimate_level
  )
  if (isTRUE(dots[["only_priors"]]))
    return(.set_only_priors_class(object))

  object$fit <- .fit(object)

  object$summary      <- .object_summary(object)
  object$coefficients <- .object_coefficients(object)

  object <- .autocompute_brma(
    object,
    loo     = TRUE,
    waic    = TRUE,
    marglik = FALSE
  )

  return(object)
}

.brma_mv_reject_unsupported_dots <- function(matched_call) {

  dots <- matched_call[["..."]]
  if (is.null(dots) || length(dots) == 0L) {
    return(invisible(TRUE))
  }

  unsupported <- intersect(
    names(dots),
    c(
      "cluster", "weights", "prior_bias", "prior_PET", "prior_PEESE",
      "model_type"
    )
  )
  if (length(unsupported) == 0L) {
    return(invisible(TRUE))
  }

  messages <- c(
    cluster = paste0(
      "'cluster' is not supported in brma.mv(); use the dedicated ",
      "'random' argument for multilevel structures."
    ),
    weights = "'weights' are not supported in brma.mv().",
    prior_bias = paste0(
      "Selection/publication-bias priors are not supported in brma.mv() yet."
    ),
    prior_PET = paste0(
      "PET publication-bias priors are not supported in brma.mv() yet."
    ),
    prior_PEESE = paste0(
      "PEESE publication-bias priors are not supported in brma.mv() yet."
    ),
    model_type = paste0(
      "'model_type' is not supported in brma.mv(); publication-bias model ",
      "types are not available for brma.mv() yet."
    )
  )
  stop(
    paste(unname(messages[unsupported]), collapse = " "),
    call. = FALSE
  )
}
