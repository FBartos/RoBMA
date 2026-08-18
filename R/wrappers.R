# ============================================================================ #
# brma.wrappers.R
# ============================================================================ #
#
# This file contains user-friendly wrapper functions for common prediction
# tasks with brma model fits. These wrappers simplify the interface to
# predict.brma() for typical use cases:
#
# - fitted.brma(): In-sample fitted values
# - pooled_effect.brma(): Aggregated pooled effect size (mu)
# - pooled_heterogeneity.brma(): Aggregated pooled heterogeneity (tau)
# - blup.brma(): Best Linear Unbiased Predictions (true effects)
#
# Design principles:
# - Simple interface: minimal required arguments
# - Consistent with predict.brma: use same underlying machinery
# - Quiet operation: suppress aggregation messages for cleaner output
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# nobs and coef methods
# ---------------------------------------------------------------------------- #

#' @title Number of Observations for brma Objects
#'
#' @description Extract the number of observations (studies/estimates) from
#' a fitted brma object.
#'
#' @param object a fitted brma object
#' @param ... additional arguments (currently ignored)
#'
#' @return An integer giving the number of observations.
#'
#' @seealso [brma()], [brma.glmm()]
#' @export
nobs.brma <- function(object, ...) {
  return(length(.outcome_data_yi(object)))
}


#' @title Extract Model Coefficients for brma Objects
#'
#' @description Extract model coefficients (posterior means) from a
#' fitted brma object.
#'
#' @param object a fitted brma object
#' @param ... additional arguments (currently ignored)
#'
#' @return A named numeric vector of posterior mean coefficients.
#'
#' @seealso [summary.brma()]
#' @export
coef.brma <- function(object, ...) {
  return(object[["coefficients"]])
}


#' @title Fitted Values for brma Objects
#'
#' @description Extract in-sample fitted values from a fitted \code{brma}
#' object.
#'
#' @param object a fitted brma object.
#' @param unit output unit. Only \code{"estimate"} is implemented currently.
#' @param conditioning_depth conditioning depth for location fitted values.
#' \code{"marginal"} uses fixed effects only, \code{"cluster"} conditions on
#' cluster-level random effects, and \code{"estimate"} conditions on the full
#' estimate-level fitted value.
#' @param component fitted component to return. \code{"location"} (alias
#' \code{"mods"}) returns location fitted values, \code{"scale"} returns
#' fitted heterogeneity \eqn{\tau_i}, and \code{"all"} returns both as a
#' named list.
#' @param bias_adjusted whether location fitted values should adjust for
#' publication bias. Defaults to \code{FALSE}.
#' @param conditional whether to return fitted values from conditional posterior
#' predictions for RoBMA product-space objects.
#' @inheritParams predict.brma
#' @param ... additional arguments. Currently only \code{quiet} is honored.
#'
#' @details
#' This method is a compact adapter around \code{\link{predict.brma}}. It
#' summarizes posterior prediction draws by their column means and returns a
#' base numeric vector, matching the usual \code{\link[stats]{fitted}} contract.
#' Use \code{predict()} directly when posterior draws or intervals are needed.
#'
#' The default \code{conditioning_depth = "marginal"} corresponds to
#' \code{predict(object, type = "terms")} and matches the usual fitted-value
#' convention for meta-regression. For normal models,
#' \code{conditioning_depth = "estimate"} corresponds to BLUP means for the
#' observed estimates.
#'
#' For \code{component = "all"}, \code{conditioning_depth},
#' \code{output_measure}, and \code{transform} apply only to the
#' \code{location} component. The \code{scale} component always returns
#' fitted \eqn{\tau_i} values.
#'
#' @return A numeric vector of fitted values, or a named list with
#' \code{location} and \code{scale} components when \code{component = "all"}.
#'
#' @seealso [predict.brma()], [residuals.brma()], [blup.brma()]
#' @export
fitted.brma <- function(object, unit = "estimate",
                        conditioning_depth = "marginal",
                        component = "location",
                        bias_adjusted = FALSE,
                        output_measure = NULL,
                        transform = NULL,
                        conditional = FALSE, ...) {

  dots                         <- list(...)
  quiet                        <- .fitted_brma_quiet(dots)
  conditioning_depth_specified <- !missing(conditioning_depth)

  unit                         <- .normalize_unit(unit)
  conditioning_depth           <- .normalize_conditioning_depth(conditioning_depth)
  component                    <- .fitted_component_normalize(component)

  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(conditional, "conditional")

  .check_unit_conditioning_depth(
    object             = object,
    unit               = unit,
    conditioning_depth = conditioning_depth,
    caller             = "fitted()"
  )

  if (unit == "cluster") {
    .check_cluster_unit_deferred("fitted()")
  }

  if (component == "scale" && conditioning_depth_specified &&
      conditioning_depth != "marginal") {
    stop(
      "'conditioning_depth' is not available for component = 'scale'.",
      call. = FALSE
    )
  }

  if (component == "scale" &&
      (!is.null(output_measure) || !is.null(transform))) {
    stop(
      "'output_measure' and 'transform' are only available for location ",
      "fitted values.",
      call. = FALSE
    )
  }

  if (component == "location") {
    return(.fitted_brma_location(
      object             = object,
      conditioning_depth = conditioning_depth,
      bias_adjusted      = bias_adjusted,
      output_measure     = output_measure,
      transform          = transform,
      conditional        = conditional,
      quiet              = quiet
    ))
  }

  if (component == "scale") {
    return(.fitted_brma_scale(
      object      = object,
      conditional = conditional,
      quiet       = quiet
    ))
  }

  return(list(
    location = .fitted_brma_location(
      object             = object,
      conditioning_depth = conditioning_depth,
      bias_adjusted      = bias_adjusted,
      output_measure     = output_measure,
      transform          = transform,
      conditional        = conditional,
      quiet              = quiet
    ),
    scale = .fitted_brma_scale(
      object      = object,
      conditional = conditional,
      quiet       = quiet
    )
  ))
}


.fitted_brma_quiet <- function(dots) {

  if (is.null(dots[["quiet"]])) {
    return(TRUE)
  }

  BayesTools::check_bool(dots[["quiet"]], "quiet")
  return(dots[["quiet"]])
}


.fitted_brma_location <- function(object, conditioning_depth,
                                  bias_adjusted, output_measure,
                                  transform, conditional, quiet) {

  pred_type <- switch(conditioning_depth,
    "marginal" = "terms",
    "cluster"  = "cluster",
    "estimate" = "blup"
  )

  samples <- predict.brma(
    object         = object,
    newdata        = NULL,
    type           = pred_type,
    output_measure = output_measure,
    transform      = transform,
    bias_adjusted  = bias_adjusted,
    conditional    = conditional,
    quiet          = quiet
  )

  return(.fitted_brma_col_means(samples, object))
}


.fitted_brma_scale <- function(object, conditional, quiet) {

  samples <- predict.brma(
    object      = object,
    newdata     = NULL,
    type        = "terms.scale",
    conditional = conditional,
    quiet       = quiet
  )

  if (is.list(samples) && !is.data.frame(samples)) {
    out <- lapply(samples, .fitted_brma_col_means, object = object)
    names(out) <- names(samples)
    return(out)
  }

  return(.fitted_brma_col_means(samples, object))
}


.fitted_brma_col_means <- function(samples, object) {

  out <- colMeans(as.matrix(samples))
  out <- .diagnostic_set_names(out, object)

  return(out)
}


# ---------------------------------------------------------------------------- #
# pooled_effect generic and brma method
# ---------------------------------------------------------------------------- #

#' @title Pooled Effect Size
#'
#' @description Computes the pooled (aggregated) effect size estimate
#' from a fitted model.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of the pooled effect size.
#'
#' @seealso [predict.brma()]
#' @export
pooled_effect <- function(object, ...) {
  UseMethod("pooled_effect")
}


#' @title Pooled Effect Size for brma Objects
#'
#' @description Computes the pooled effect and prediction interval at the
#' average expanded moderator design of a fitted brma object.
#'
#' @param object a fitted brma object
#' @param bias_adjusted whether to adjust for publication bias. Defaults to
#' \code{TRUE}. For PET/PEESE models this removes the regression bias term from
#' the pooled location effect. Selection-model weighting affects response
#' predictions, not this fixed-location posterior summary.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param conditional whether to return the pooled effect conditional on the
#' effect component for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @inheritParams predict.brma
#' @param ... reserved for internal posterior-sample reuse.
#'
#' @details
#' This function evaluates the fixed-location model at the average expanded
#' moderator design. Because the location predictor is linear in its design,
#' this equals averaging its fitted-design posterior draws directly.
#'
#' For meta-regression models, the pooled effect averages the effect size
#' estimate across moderator levels proportionately to the levels observed
#' in the data. This provides an estimate representative of the sample of
#' studies.
#'
#' For models without moderators, this returns the single mu parameter. The
#' prediction interval is marginal over one new true effect at the same average
#' location, scale, and random-effect design. It therefore uses
#' \code{pooled_heterogeneity()}, not the observed-design RMS heterogeneity
#' reported by \code{summary_heterogeneity()}.
#'
#' @return A \code{brma_samples} object containing posterior samples. When printed,
#' displays a one-row summary table whose \code{PI} columns contain posterior
#' prediction quantiles. Use \code{summary()} or \code{as.data.frame()} to
#' obtain the table directly.
#' The samples can be converted to \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   pooled_effect(fit)
#' }
#' }
#'
#' @seealso [predict.brma()], [pooled_heterogeneity()], [blup()]
#' @export
pooled_effect.brma <- function(object, bias_adjusted = TRUE,
                               output_measure = NULL, transform = NULL,
                               probs = c(.025, .975),
                               conditional = FALSE, ...) {
  dots <- list(...)
  .check_unused_dots(
    dots    = dots,
    allowed = ".posterior_samples",
    caller  = "pooled_effect.brma()"
  )
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' pooled effects are available only for RoBMA objects.",
         call. = FALSE)
  }

  data              <- object[["data"]]
  priors            <- object[["priors"]]
  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )
  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (.is_mods(object)) {
      .create_fit_formula_list(data = data, "mods")
    } else {
      NULL
    },
    mods_priors       = if (.is_random(object)) {
      priors[["location"]]
    } else {
      priors[["mods"]]
    },
    is_mods           = .is_mods(object),
    is_PET            = .is_PET(object),
    is_PEESE          = .is_PEESE(object),
    effect_direction  = .effect_direction(object),
    bias_adjusted     = bias_adjusted,
    K                 = nrow(data[["outcome"]]),
    posterior_samples = posterior_samples,
    priors            = priors
  )
  samples <- matrix(rowMeans(mu_samples), ncol = 1L)
  colnames(samples) <- "mu"

  tau_samples <- .pooled_heterogeneity_total_samples(
    object            = object,
    posterior_samples = posterior_samples
  )
  prediction_samples <- matrix(
    stats::rnorm(
      n    = nrow(samples),
      mean = samples[, 1L],
      sd   = tau_samples[, 1L]
    ),
    ncol     = 1L,
    dimnames = list(NULL, "mu")
  )

  chain_info <- .brma_samples_chain_info(
    fit       = object[["fit"]],
    n_samples = nrow(posterior_samples)
  )
  effect_transform <- .effect_output_setup(
    object         = object,
    output_measure = output_measure,
    transform      = transform
  )
  out <- .new_effect_brma_samples(
    samples            = samples,
    n_chains           = chain_info[["n_chains"]],
    n_iter             = chain_info[["n_iter"]],
    title              = "Pooled Effect Size",
    probs              = probs,
    data               = NULL,
    effect_transform   = effect_transform,
    prediction_samples = prediction_samples
  )
  return(.condition_prediction_samples(
    object            = object,
    samples           = out,
    conditional       = conditional,
    parameters        = .conditional_effect_parameters(object),
    posterior_samples = posterior_samples,
    quiet             = TRUE
  ))
}


# ---------------------------------------------------------------------------- #
# pooled_heterogeneity generic and brma method
# ---------------------------------------------------------------------------- #

#' @title Pooled Heterogeneity
#'
#' @description Computes the pooled (aggregated) heterogeneity estimate (tau)
#' from a fitted model.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of the pooled heterogeneity.
#'
#' @seealso [predict.brma()]
#' @export
pooled_heterogeneity <- function(object, ...) {
  UseMethod("pooled_heterogeneity")
}


#' @title Pooled Heterogeneity for brma Objects
#'
#' @description Computes heterogeneity (tau) at the average expanded scale and
#' random-effect design of a fitted brma object.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param conditional whether to return the pooled heterogeneity conditional on
#' the heterogeneity component for RoBMA product-space objects. Defaults to
#' \code{FALSE}.
#' @param component heterogeneity component to return for \code{brma.mv()}
#' models. Defaults to \code{"all"}. Use \code{"total"} for the
#' variance-additive total heterogeneity.
#' @param ... reserved for internal posterior-sample forwarding.
#'
#' @details
#' For location-scale models, the function averages the fitted log-scale linear
#' predictor and then applies the inverse link. This is heterogeneity at the
#' average expanded scale design, rather than RMS heterogeneity across observed
#' rows. Use \code{summary_heterogeneity()} for the latter.
#'
#' For models without scale regression, this returns the single tau parameter.
#'
#' For multilevel (3-level) models, the returned tau is the total heterogeneity:
#' \code{tau = sqrt(tau_within^2 + tau_between^2)}.
#'
#' For \code{brma.mv()} random-formula models, \code{component = "all"}
#' returns one \code{brma_samples} object when there is a single heterogeneity
#' component and a named list when there are multiple components.
#' \code{component = "total"} computes the variance-additive total at the
#' average random-effect design. Random slopes are evaluated as
#' \eqn{\bar q^T G \bar q}; known-\eqn{R} group multipliers are averaged over
#' the fitted rows for the one-study marginal variance target.
#'
#' @return A \code{brma_samples} object containing posterior samples. When printed,
#' displays a summary table. For decomposed \code{brma.mv()} models, a named
#' list of \code{brma_samples} objects is returned. Use \code{summary()} on an
#' individual object or \code{as.data.frame()} on the complete result to obtain
#' a summary data frame with component and parameter identifiers. Use
#' \code{as.data.frame(format = "list")} to retain separate component tables.
#' The samples can be converted to \pkg{posterior} draws formats using
#' \code{as_draws()}.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   pooled_heterogeneity(fit)
#' }
#' }
#'
#' @seealso [predict.brma()], [pooled_effect()], [blup()]
#' @export
pooled_heterogeneity.brma <- function(object, probs = c(.025, .975),
                                      conditional = FALSE, component = "all",
                                      ...) {
  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' pooled heterogeneity is available only for RoBMA objects.",
         call. = FALSE)
  }
  if (inherits(object, "brma.mv")) {
    return(
      .pooled_heterogeneity_brma_mv(
        object      = object,
        probs       = probs,
        conditional = conditional,
        component   = component,
        ...
      )
    )
  }

  .check_univariate_heterogeneity_component(component)

  dots <- list(...)
  .check_unused_dots(
    dots    = dots,
    allowed = ".posterior_samples",
    caller  = "pooled_heterogeneity.brma()"
  )
  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )
  samples <- .pooled_heterogeneity_total_samples(
    object            = object,
    posterior_samples = posterior_samples
  )
  colnames(samples) <- "tau"

  chain_info <- .brma_samples_chain_info(
    fit       = object[["fit"]],
    n_samples = nrow(posterior_samples)
  )
  out <- .new_brma_samples(
    samples  = samples,
    n_chains = chain_info[["n_chains"]],
    n_iter   = chain_info[["n_iter"]],
    title    = "Pooled Heterogeneity",
    probs    = probs,
    data     = NULL
  )
  return(.condition_prediction_samples(
    object            = object,
    samples           = out,
    conditional       = conditional,
    parameters        = .conditional_heterogeneity_parameters(object),
    posterior_samples = posterior_samples,
    quiet             = TRUE
  ))
}


.pooled_heterogeneity_total_samples <- function(object, posterior_samples) {

  if (inherits(object, "brma.mv")) {
    components <- .pooled_brma_mv_heterogeneity_components(
      object            = object,
      posterior_samples = posterior_samples
    )
    return(.total_brma_mv_heterogeneity_samples(components))
  }

  tau_result <- .evaluate.brma.pooled_tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) {
      .create_fit_formula_list(data = object[["data"]], "scale")
    } else {
      NULL
    },
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = .is_multilevel(object),
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )

  return(tau_result[["tau_total"]])
}


# ---------------------------------------------------------------------------- #
# blup generic and brma method
# ---------------------------------------------------------------------------- #

#' @title Best Linear Unbiased Predictions (BLUPs)
#'
#' @description Computes the estimated true effects (theta) from a
#' fitted model. These correspond to Best Linear Unbiased Predictions (BLUPs)
#' or empirical Bayes estimates.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of BLUP or empirical-Bayes true-effect summaries.
#'
#' @seealso [predict.brma()]
#' @export
blup <- function(object, ...) {
  UseMethod("blup")
}


#' @title Best Linear Unbiased Predictions for brma Objects
#'
#' @description Computes the estimated true effects (theta) for a
#' fitted brma object. These correspond to Best Linear Unbiased Predictions
#' (BLUPs) or empirical Bayes estimates.
#'
#' @param object a fitted brma object
#' @param bias_adjusted whether to adjust for publication bias. Defaults to
#' \code{FALSE}, which returns estimates including publication bias effects
#' (i.e., what we expect the true effects to be given the biased observations).
#' Set to \code{TRUE} to obtain bias-corrected estimates.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @inheritParams predict.brma
#' @param ... additional arguments passed to \code{\link{predict.brma}}; wrapper
#' arguments such as \code{newdata}, \code{type}, \code{quiet},
#' \code{output_measure}, and \code{transform} are fixed by this method.
#'
#' @details
#' This function is a convenience wrapper around \code{predict.brma(...,
#' type = "blup", newdata = NULL)}. Unlike \code{predict(..., type =
#' "estimate", conditioning_depth = "estimate")}, it returns conditional
#' location/BLUP means rather than adding conditional latent-effect uncertainty.
#'
#' For unweighted two-level normal models, true effects are computed using
#' empirical Bayes shrinkage:
#' \deqn{\theta_i = \lambda_i \cdot y_i + (1 - \lambda_i) \cdot \mu_i}
#' where \eqn{\lambda_i = \tau^2 / (\tau^2 + se_i^2)}.
#' With likelihood weights, \eqn{se_i^2} is replaced by the weighted sampling
#' variance \eqn{se_i^2 / w_i}.
#'
#' For GLMM models (binomial, Poisson), the estimate-level random effects
#' are extracted directly from the posterior samples.
#'
#' For multilevel (3-level) normal models, cluster-level effects are estimated
#' jointly within cluster blocks and estimate-level effects are then shrunk
#' conditional on those cluster effects.
#'
#' @return A \code{brma_samples} object containing posterior draws of BLUP or
#' empirical-Bayes true-effect summaries with one column per estimate. For
#' existing normal data, these are conditional BLUP means, not simulated
#' latent-effect draws. When printed, displays a summary table. Use
#' \code{summary()} or \code{as.data.frame()} to obtain the summary table
#' directly. The samples can be
#' converted to \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   blup(fit)
#' }
#' }
#'
#' @seealso [predict.brma()], [pooled_effect()], [pooled_heterogeneity()],
#' [true_effects()]
#' @export
blup.brma <- function(object, bias_adjusted = FALSE,
                      output_measure = NULL, transform = NULL,
                      probs = c(.025, .975), ...) {
  out <- predict.brma(
    object         = object,
    newdata        = NULL,
    type           = "blup",
    output_measure = output_measure,
    transform      = transform,
    probs          = probs,
    bias_adjusted  = bias_adjusted,
    quiet          = TRUE,
    ...
  )
  attr(out, "title") <- .effect_output_title(
    title            = "True Effects (BLUP Means)",
    effect_transform = attr(out, "effect_transform")
  )
  return(out)
}


# ---------------------------------------------------------------------------- #
# true_effects generic and brma method (alias for blup)
# ---------------------------------------------------------------------------- #

#' @title True Effects
#'
#' @description Computes the estimated true effects (theta) from a
#' fitted model. This is a separate S3 generic whose \code{brma} method
#' delegates to \code{\link{blup.brma}}.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of BLUP or empirical-Bayes true-effect summaries.
#'
#' @seealso [blup()], [predict.brma()]
#' @export
true_effects <- function(object, ...) {
  UseMethod("true_effects")
}


#' @title True Effects for brma Objects
#'
#' @description Computes the estimated true effects (theta) for a
#' fitted brma object. This is an alias for \code{\link{blup.brma}}.
#'
#' @inheritParams blup.brma
#'
#' @details
#' This function is identical to \code{\link{blup.brma}}. See that function
#' for full details on how true effects are computed.
#'
#' @return A \code{brma_samples} object containing posterior draws of BLUP or
#' empirical-Bayes true-effect summaries with one column per estimate. For
#' existing normal data, these are conditional BLUP means, not simulated
#' latent-effect draws. When printed, displays a summary table. Use
#' \code{summary()} or \code{as.data.frame()} to obtain the summary table
#' directly. The samples can be
#' converted to \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   true_effects(fit)
#' }
#' }
#'
#' @seealso [blup.brma()], [predict.brma()], [pooled_effect()],
#' [pooled_heterogeneity()]
#' @export
true_effects.brma <- function(object, bias_adjusted = FALSE,
                              output_measure = NULL, transform = NULL,
                              probs = c(.025, .975), ...) {
  blup.brma(
    object         = object,
    bias_adjusted  = bias_adjusted,
    output_measure = output_measure,
    transform      = transform,
    probs          = probs,
    ...
  )
}


# ---------------------------------------------------------------------------- #
# ranef generic and brma method
# ---------------------------------------------------------------------------- #

#' @title Random Effects
#'
#' @description Extracts the estimated random effects (deviations from the
#' fixed-effect predictions) from a fitted model.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically posterior samples of the
#' random effect deviations.
#'
#' @seealso [blup()], [predict.brma()]
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}


#' @title Random Effects for brma Objects
#'
#' @description Extracts random effect deviations from a fitted brma object.
#' These are posterior-sample offsets from the fixed-effect predictions,
#' analogous to random-effect deviations returned by \code{metafor::ranef()}.
#'
#' @param object a fitted brma object
#' @param bias_adjusted whether to adjust for publication bias. Defaults to
#' \code{FALSE}. See \code{\link{blup.brma}} for details.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param component random-effect component to return. Defaults to
#' \code{"all"}. For decomposed models, use a component name such as
#' \code{"cluster"}, \code{"estimate"}, or a \code{brma.mv()} random-effect
#' block name. Use \code{"total"} for the summed random-effect deviation.
#' @param simplify whether \code{component = "all"} should return a single
#' \code{brma_samples} object instead of a one-element list. Defaults to
#' \code{TRUE}, matching standard 2-level \code{brma()} behavior.
#' @param expand whether to repeat random-effect contributions for every fitted
#' observation. Defaults to \code{FALSE}, returning one column per unique
#' grouping-level contribution in first fitted-observation order, consistently
#' with \code{metafor::ranef()}. Indicator-coded random coefficients retain one
#' column per observed grouping-level and coefficient combination. Set to
#' \code{TRUE} for observation-aligned output.
#' @param ... additional arguments forwarded to \code{\link{predict.brma}} for
#' supported options such as \code{conditional}. \code{newdata}, \code{type},
#' \code{quiet}, \code{output_measure}, and \code{transform} are controlled by
#' this method.
#'
#' @details
#' Random effects are computed as the difference between the true
#' effects (BLUPs) and the fixed-effect predictions:
#' \deqn{u_i = \hat{\theta}_i - \hat{\mu}_i}
#'
#' For standard (2-level) models, returns a single \code{brma_samples}
#' object with the estimate-level random effects.
#'
#' For multilevel (3-level) models, returns a list with two
#' \code{brma_samples} matrices:
#' \describe{
#'   \item{\code{cluster}}{Cluster-level random effects
#'     (\eqn{\gamma_j \cdot \tau_{between}}), representing between-cluster
#'     deviations from the fixed effects.}
#'   \item{\code{estimate}}{Estimate-level random effects
#'     (\eqn{\theta_i - \mu_i - \gamma_j \cdot \tau_{between}}),
#'     representing within-cluster deviations from the cluster means.}
#' }
#'
#' For \code{brma.mv()} random-formula models, decomposes the Gaussian
#' conditional means by random-effect block. The result is invariant to whether
#' each block was sampled or marginalized during fitting. If there is only one
#' block and \code{simplify = TRUE}, returns a single \code{brma_samples} object.
#' Multiple blocks are returned in one flat list under their canonical fitted
#' block names; no additional location layer or component-name prefix is added.
#' Unique-level output is available when a block's contribution is constant
#' within each grouping level or observed grouping-level and indicator-
#' coefficient combination. Other random slopes and row-varying random-effect
#' scales generally require \code{expand = TRUE}. Summed \code{component =
#' "total"} output across multiple blocks also requires \code{expand = TRUE}
#' because different blocks need not share grouping levels.
#'
#' @return A \code{brma_samples} object for a single selected/simplified
#' component, or a named list of \code{brma_samples} objects for decomposed
#' \code{component = "all"} output. Use \code{as.data.frame()} to obtain the
#' displayed summary table in long form, or
#' \code{as.data.frame(format = "list")} to retain separate component tables.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   ranef(fit)
#' }
#' }
#'
#' @seealso [blup.brma()], [predict.brma()], [pooled_effect()]
#' @export
ranef.brma <- function(object, bias_adjusted = FALSE,
                       probs = c(.025, .975), component = "all",
                       simplify = TRUE, expand = FALSE, ...) {

  BayesTools::check_bool(expand, "expand")
  dots <- list(...)
  if (!is.null(dots[["output_measure"]]) || !is.null(dots[["transform"]])) {
    stop(
      "'output_measure' and 'transform' are not available for random-effect ",
      "deviations.",
      call. = FALSE
    )
  }

  is_multilevel <- .is_multilevel(object)

  if (inherits(object, "brma.mv") && .is_random(object)) {
    return(.ranef_brma_mv_random(
      object        = object,
      bias_adjusted = bias_adjusted,
      probs         = probs,
      component     = component,
      simplify      = simplify,
      expand        = expand,
      ...
    ))
  }

  if (is_multilevel && .outcome_type(object) == "norm" &&
      !.is_weightfunction(object)) {
    return(.ranef_brma_multilevel_normal(
      object        = object,
      bias_adjusted = bias_adjusted,
      probs         = probs,
      component     = component,
      simplify      = simplify,
      expand        = expand,
      dots          = dots
    ))
  }

  # get BLUPs (fixed + all random effects)
  blup_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "blup",
    probs         = probs,
    bias_adjusted = bias_adjusted,
    quiet         = TRUE,
    ...
  )

  # extract MCMC chain info from the returned samples
  n_chains <- attr(blup_samples, "nchains")
  n_iter   <- attr(blup_samples, "niter")
  data     <- object[["data"]]
  labels   <- .get_estimate_labels(object)

  # get fixed-effect predictions only
  terms_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "terms",
    probs         = probs,
    bias_adjusted = bias_adjusted,
    quiet         = TRUE,
    ...
  )

  if (!is_multilevel) {

    # 2-level: random effects = BLUP - fixed effects
    ranef_mat <- unclass(blup_samples) - unclass(terms_samples)
    K         <- ncol(ranef_mat)
    colnames(ranef_mat) <- paste0("u[", labels[seq_len(K)], "]")

    estimate_ranef <- .new_brma_samples(
      samples  = ranef_mat,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = "Random Effects:",
      probs    = probs,
      data     = data
    )

    return(.select_ranef_components(
      components = list(estimate = estimate_ranef),
      component  = component,
      simplify   = simplify,
      expand     = expand,
      labels     = labels,
      n_chains   = n_chains,
      n_iter     = n_iter,
      probs      = probs,
      data       = data
    ))

  } else {

    # 3-level: decompose into cluster-level and estimate-level
    cluster_samples <- predict.brma(
      object        = object,
      newdata       = NULL,
      type          = "cluster",
      probs         = probs,
      bias_adjusted = bias_adjusted,
      quiet         = TRUE,
      ...
    )

    # cluster-level random effects: cluster predictions - fixed effects
    cluster_ranef_mat <- unclass(cluster_samples) - unclass(terms_samples)
    K                 <- ncol(cluster_ranef_mat)
    cluster          <- data[["outcome"]][["cluster"]]
    cluster_labels   <- .get_cluster_labels(object)[as.character(cluster)]
    cluster_labels   <- unname(cluster_labels)
    cluster_missing  <- is.na(cluster_labels)
    cluster_labels[cluster_missing] <- as.character(cluster[cluster_missing])
    if (any(duplicated(cluster_labels))) {
      cluster_names <- paste0(cluster_labels, "|", labels[seq_len(K)])
    } else {
      cluster_names <- cluster_labels
    }
    if (expand) {
      colnames(cluster_ranef_mat) <- paste0("u_cluster[", cluster_names, "]")
    } else {
      cluster_levels         <- unique(cluster)
      cluster_level_map      <- match(cluster, cluster_levels)
      cluster_level_labels   <- .get_cluster_labels(object)[
        as.character(cluster_levels)
      ]
      cluster_level_labels   <- unname(cluster_level_labels)
      cluster_level_missing <- is.na(cluster_level_labels)
      cluster_level_labels[cluster_level_missing] <- as.character(
        cluster_levels[cluster_level_missing]
      )
      cluster_ranef_mat <- .ranef_unique_level_samples(
        samples      = cluster_ranef_mat,
        group_map    = cluster_level_map,
        group_levels = cluster_level_labels,
        block        = "cluster"
      )
    }

    cluster_ranef <- .new_brma_samples(
      samples  = cluster_ranef_mat,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = "Cluster-Level Random Effects:",
      probs    = probs,
      data     = data
    )

    # estimate-level random effects: BLUP - cluster predictions
    estimate_ranef_mat <- unclass(blup_samples) - unclass(cluster_samples)
    colnames(estimate_ranef_mat) <- paste0("u_estimate[", labels[seq_len(K)], "]")

    estimate_ranef <- .new_brma_samples(
      samples  = estimate_ranef_mat,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = "Estimate-Level Random Effects:",
      probs    = probs,
      data     = data
    )

    out <- list(
      cluster  = cluster_ranef,
      estimate = estimate_ranef
    )
    return(.select_ranef_components(
      components = out,
      component  = component,
      simplify   = simplify,
      expand     = expand,
      labels     = labels,
      n_chains   = n_chains,
      n_iter     = n_iter,
      probs      = probs,
      data       = data
    ))
  }
}


.ranef_brma_multilevel_normal <- function(
    object, bias_adjusted, probs, component, simplify, expand, dots) {

  .check_unused_dots(
    dots    = dots,
    allowed = c("conditional", ".posterior_samples"),
    caller  = "ranef.brma()"
  )
  conditional <- dots[["conditional"]]
  if (is.null(conditional)) {
    conditional <- FALSE
  }

  context <- .predict_brma_context(
    object                       = object,
    newdata                      = NULL,
    V_new                        = NULL,
    type                         = "blup",
    conditioning_depth           = "estimate",
    conditioning_depth_specified = TRUE,
    as_measure                   = TRUE,
    output_measure               = NULL,
    transform                    = NULL,
    probs                        = probs,
    bias_adjusted                = bias_adjusted,
    quiet                        = TRUE,
    conditional                  = conditional,
    dots                         = list(.posterior_samples = dots[[".posterior_samples"]])
  )
  scale_state    <- .predict_brma_scale_state(context)
  location_state <- .predict_brma_location_state(context, scale_state)
  components     <- location_state[["multilevel_blup"]]

  if (is.null(components)) {
    stop("The exact multilevel normal BLUP components are unavailable.",
         call. = FALSE)
  }

  labels          <- .get_estimate_labels(object)
  cluster         <- context[["new_data"]][["outcome"]][["cluster"]]
  cluster_labels  <- .get_cluster_labels(object)[as.character(cluster)]
  cluster_labels  <- unname(cluster_labels)
  cluster_missing <- is.na(cluster_labels)
  cluster_labels[cluster_missing] <- as.character(cluster[cluster_missing])
  if (any(duplicated(cluster_labels))) {
    cluster_names <- paste0(cluster_labels, "|", labels)
  } else {
    cluster_names <- cluster_labels
  }

  if (expand) {
    colnames(components[["cluster"]]) <- paste0(
      "u_cluster[", cluster_names, "]"
    )
  } else {
    cluster_levels         <- unique(cluster)
    cluster_level_map      <- match(cluster, cluster_levels)
    cluster_level_labels   <- .get_cluster_labels(object)[
      as.character(cluster_levels)
    ]
    cluster_level_labels   <- unname(cluster_level_labels)
    cluster_level_missing <- is.na(cluster_level_labels)
    cluster_level_labels[cluster_level_missing] <- as.character(
      cluster_levels[cluster_level_missing]
    )
    components[["cluster"]] <- .ranef_unique_level_samples(
      samples      = components[["cluster"]],
      group_map    = cluster_level_map,
      group_levels = cluster_level_labels,
      block        = "cluster"
    )
  }
  colnames(components[["estimate"]]) <- paste0(
    "u_estimate[", labels, "]"
  )

  out <- list(
    cluster = .new_brma_samples(
      samples  = components[["cluster"]],
      n_chains = context[["n_chains"]],
      n_iter   = context[["n_iter"]],
      title    = "Cluster-Level Random Effects:",
      probs    = probs,
      data     = context[["new_data"]]
    ),
    estimate = .new_brma_samples(
      samples  = components[["estimate"]],
      n_chains = context[["n_chains"]],
      n_iter   = context[["n_iter"]],
      title    = "Estimate-Level Random Effects:",
      probs    = probs,
      data     = context[["new_data"]]
    )
  )

  parameters <- if (conditional) {
    .conditional_effect_parameters(object)
  } else {
    NULL
  }
  out <- .condition_prediction_samples(
    object            = object,
    samples           = out,
    conditional       = conditional,
    parameters        = parameters,
    posterior_samples = context[["posterior_samples"]],
    quiet             = TRUE
  )

  n_chains <- attr(out[[1L]], "nchains")
  n_iter   <- attr(out[[1L]], "niter")
  return(.select_ranef_components(
    components = out,
    component  = component,
    simplify   = simplify,
    expand     = expand,
    labels     = labels,
    n_chains   = n_chains,
    n_iter     = n_iter,
    probs      = probs,
    data       = context[["new_data"]]
  ))
}


.ranef_brma_mv_random <- function(object, bias_adjusted = FALSE,
                                  probs = c(.025, .975), component = "all",
                                  simplify = TRUE, expand = FALSE, ...) {

  dots              <- list(...)
  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )

  terms_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "terms",
    probs         = probs,
    bias_adjusted = bias_adjusted,
    quiet         = TRUE,
    ...
  )

  n_chains <- attr(terms_samples, "nchains")
  n_iter   <- attr(terms_samples, "niter")
  data     <- object[["data"]]
  labels   <- .get_estimate_labels(object)
  K        <- ncol(terms_samples)

  bias_offset <- NULL
  if (bias_adjusted && (.is_PET(object) || .is_PEESE(object))) {
    bias_offset <- .evaluate.brma.bias_offset(
      fit               = object[["fit"]],
      outcome_data      = data[["outcome"]],
      is_PET            = .is_PET(object),
      is_PEESE          = .is_PEESE(object),
      effect_direction  = .effect_direction(object),
      K                 = K,
      posterior_samples = posterior_samples,
      priors            = object[["priors"]]
    )
  }

  components <- .evaluate.brma.mv_random_blup.norm(
    object            = object,
    mu_samples        = unclass(terms_samples),
    posterior_samples = posterior_samples,
    bias_offset       = bias_offset,
    by_block          = TRUE
  )

  design       <- .fitted_formula_design(object, "mu", required = TRUE)
  random_terms <- design[["random_effects"]]
  block_order  <- vapply(
    random_terms,
    `[[`,
    character(1),
    "block_name"
  )
  names(random_terms) <- block_order
  block_order <- block_order[block_order %in% names(components)]
  components  <- components[block_order]

  out <- lapply(names(components), function(block) {
    mat <- components[[block]]
    if (expand) {
      colnames(mat) <- paste0("u_", block, "[", labels[seq_len(K)], "]")
    } else {
      term     <- random_terms[[block]]
      metadata <- .ranef_unique_level_term(term, block)
      mat <- .ranef_unique_level_samples(
        samples      = mat,
        group_map    = metadata[["group_map"]],
        group_levels = metadata[["group_levels"]],
        block        = block
      )
    }
    .new_brma_samples(
      samples  = mat,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = paste0("Random Effects: ", block),
      probs    = probs,
      data     = data
    )
  })
  names(out) <- names(components)

  return(.select_ranef_components(
    components = out,
    component  = component,
    simplify   = simplify,
    expand     = expand,
    labels     = labels,
    n_chains   = n_chains,
    n_iter     = n_iter,
    probs      = probs,
    data       = data
  ))
}


.select_ranef_components <- function(components, component = "all",
                                     simplify = TRUE, expand = FALSE,
                                     labels, n_chains,
                                     n_iter, probs, data) {

  BayesTools::check_bool(simplify, "simplify")
  BayesTools::check_bool(expand, "expand")
  if (is.null(component)) {
    component <- "all"
  }
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    stop("'component' must be a single component name.", call. = FALSE)
  }

  component_names <- names(components)
  if (component == "all") {
    if (isTRUE(simplify) && length(components) == 1L) {
      return(components[[1L]])
    }
    return(.new_brma_samples_list(components))
  }

  if (component == "total") {
    .check_ranef_total_expansion(components, expand)
    total <- Reduce(`+`, lapply(components, as.matrix))
    colnames(total) <- paste0("u[", labels[seq_len(ncol(total))], "]")
    return(.new_brma_samples(
      samples  = total,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = "Random Effects:",
      probs    = probs,
      data     = data
    ))
  }

  if (!component %in% component_names) {
    stop(
      "Unknown random-effect component '", component, "'. Available components: ",
      paste(c(component_names, "all", "total"), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  return(components[[component]])
}


# Collapse observation-aligned random-effect samples to grouping levels.
.ranef_unique_level_samples <- function(samples, group_map, group_levels,
                                        block) {

  group_map    <- as.integer(group_map)
  group_levels <- as.character(group_levels)
  n_groups     <- length(group_levels)
  if (length(group_map) != ncol(samples) || n_groups < 1L ||
      anyNA(group_map) || any(group_map < 1L | group_map > n_groups)) {
    stop("Internal error: invalid random-effect grouping metadata.",
         call. = FALSE)
  }

  group_order <- unique(group_map)
  if (length(group_order) != n_groups) {
    stop("Internal error: random-effect grouping levels are missing fitted rows.",
         call. = FALSE)
  }
  level_rows <- match(group_order, group_map)

  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(samples))
  for (group in seq_len(n_groups)) {
    rows <- which(group_map == group)
    if (length(rows) < 2L) {
      next
    }
    reference  <- samples[, rows[[1L]]]
    difference <- abs(sweep(samples[, rows, drop = FALSE], 1L, reference))
    if (any(difference > tolerance)) {
      stop(
        "Unique-level random effects are unavailable for block '", block,
        "' because it has row-specific contributions. Use 'expand = TRUE'.",
        call. = FALSE
      )
    }
  }

  out <- samples[, level_rows, drop = FALSE]
  colnames(out) <- paste0(
    "u_", block, "[", group_levels[group_order], "]"
  )
  return(out)
}


# Resolve the unique random-effect cells represented by a fitted block.
.ranef_unique_level_term <- function(term, block) {

  model_matrix <- term[["model_matrix"]]
  is_intercept <- is.matrix(model_matrix) && ncol(model_matrix) == 1L &&
    all(abs(model_matrix[, 1L] - 1) <= sqrt(.Machine$double.eps))
  if (is_intercept) {
    return(list(
      group_map    = term[["group_map"]],
      group_levels = term[["group_levels"]]
    ))
  }

  is_indicator <- is.matrix(model_matrix) && is.numeric(model_matrix) &&
    ncol(model_matrix) > 1L && all(is.finite(model_matrix)) &&
    all(model_matrix == 0 | model_matrix == 1) &&
    all(rowSums(model_matrix) == 1)
  if (!is_indicator) {
    stop(
      "Unique-level random effects are unavailable for random-slope block '",
      block, "'. Use 'expand = TRUE'.",
      call. = FALSE
    )
  }

  coefficient_labels <- term[["sd_leaves"]][["leaf_terms_by_column"]]
  if (!is.character(coefficient_labels) ||
      length(coefficient_labels) != ncol(model_matrix) ||
      anyNA(coefficient_labels) || any(!nzchar(coefficient_labels)) ||
      anyDuplicated(coefficient_labels)) {
    stop(
      "Unique-level random effects are unavailable for block '", block,
      "' because its semantic coefficient labels are missing. Use ",
      "'expand = TRUE'.",
      call. = FALSE
    )
  }

  group_map         <- term[["group_map"]]
  group_levels      <- term[["group_levels"]]
  coefficient_index <- max.col(model_matrix, ties.method = "first")
  cell_key          <- paste(group_map, coefficient_index, sep = ":")
  cell_rows         <- !duplicated(cell_key)
  cell_map          <- match(cell_key, cell_key[cell_rows])
  cell_levels       <- paste0(
    coefficient_labels[coefficient_index[cell_rows]],
    " | ",
    group_levels[group_map[cell_rows]]
  )

  return(list(
    group_map    = cell_map,
    group_levels = cell_levels
  ))
}


# Require common observation rows before summing multiple random-effect blocks.
.check_ranef_total_expansion <- function(components, expand) {

  if (!expand && length(components) > 1L) {
    stop(
      "'component = \"total\"' across multiple random-effect blocks requires ",
      "'expand = TRUE'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
