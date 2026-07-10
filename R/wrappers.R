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
    "estimate" = "estimate"
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
#' @description Computes the pooled (aggregated) effect size estimate
#' from a fitted brma object by averaging across the moderator model matrix.
#'
#' @param object a fitted brma object
#' @param bias_adjusted whether to adjust for publication bias. Defaults to
#' \code{TRUE}. For PET/PEESE models this removes the regression bias term from
#' the pooled location effect. Selection-model weighting affects response
#' predictions, not this \code{type = "terms"} wrapper.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param conditional whether to return the pooled effect conditional on the
#' effect component for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @inheritParams predict.brma
#' @param ... additional arguments passed to \code{\link{predict.brma}}; wrapper
#' arguments such as \code{newdata}, \code{type}, and \code{quiet} are fixed.
#'
#' @details
#' This function is a convenience wrapper around \code{predict.brma(...,
#' type = "terms", newdata = TRUE, bias_adjusted = TRUE, quiet = TRUE)}.
#'
#' For meta-regression models, the pooled effect averages the effect size
#' estimate across moderator levels proportionately to the levels observed
#' in the data. This provides an estimate representative of the sample of
#' studies.
#'
#' For models without moderators, this returns the single mu parameter.
#'
#' @return A \code{brma_samples} object containing posterior samples. When printed,
#' displays a summary table. Use \code{summary()} to obtain the summary table directly.
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
  out <- predict.brma(
    object         = object,
    newdata        = TRUE,
    type           = "terms",
    output_measure = output_measure,
    transform      = transform,
    probs          = probs,
    bias_adjusted  = bias_adjusted,
    conditional    = conditional,
    quiet          = TRUE,
    ...
  )
  attr(out, "title") <- .effect_output_title(
    title            = if (conditional) {
      "Conditional Pooled Effect Size"
    } else {
      "Pooled Effect Size"
    },
    effect_transform = attr(out, "effect_transform")
  )
  return(out)
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
#' @description Computes the pooled (aggregated) heterogeneity estimate (tau)
#' from a fitted brma object by pooling heterogeneity variances across the
#' scale model matrix.
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
#' @param ... additional arguments passed to \code{\link{predict.brma}}; wrapper
#' arguments such as \code{newdata}, \code{type}, and \code{quiet} are fixed.
#'
#' @details
#' For location-scale models (with scale regression), the pooled heterogeneity
#' is the square root of the average observation-specific heterogeneity
#' variance, \eqn{\sqrt{mean(\tau_i^2)}}. This differs from
#' \code{predict(..., type = "terms.scale", newdata = TRUE)}, which reports the
#' average predicted standard deviation.
#'
#' For models without scale regression, this returns the single tau parameter.
#'
#' For multilevel (3-level) models, the returned tau is the total heterogeneity:
#' \code{tau = sqrt(tau_within^2 + tau_between^2)}.
#'
#' For \code{brma.mv()} random-formula models, \code{component = "all"}
#' returns one \code{brma_samples} object when there is a single heterogeneity
#' component and a named list when there are multiple components.
#' \code{component = "total"} computes the variance-additive row-level total
#' \eqn{\sqrt{\sum_j \tau_{ij}^2}} and returns its root-mean-square value
#' across rows. Shared allocation nodes can be selected with the same
#' component aliases as \code{summary_heterogeneity()}.
#'
#' @return A \code{brma_samples} object containing posterior samples. When printed,
#' displays a summary table. For decomposed \code{brma.mv()} models, a named
#' list of \code{brma_samples} objects is returned. Use \code{summary()} to
#' obtain the summary table directly. The samples can be converted to
#' \pkg{posterior} draws formats using \code{as_draws()}.
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

  scale_samples <- predict.brma(
    object      = object,
    newdata     = NULL,
    type        = "terms.scale",
    probs       = probs,
    conditional = conditional,
    quiet       = TRUE,
    ...
  )
  samples <- matrix(sqrt(rowMeans(as.matrix(scale_samples)^2)), ncol = 1L)
  colnames(samples) <- "tau"

  out <- .new_brma_samples(
    samples  = samples,
    n_chains = attr(scale_samples, "nchains"),
    n_iter   = attr(scale_samples, "niter"),
    title    = if (conditional) {
      "Conditional Pooled Heterogeneity"
    } else {
      "Pooled Heterogeneity"
    },
    probs    = probs,
    data     = NULL
  )
  return(out)
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
#' type = "effect", newdata = NULL)}.
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
#' \code{summary()} to obtain the summary table directly. The samples can be
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
    type           = "effect",
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
#' \code{summary()} to obtain the summary table directly. The samples can be
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
#' For multilevel (3-level) models, returns a list with two observation-aligned
#' \code{brma_samples} matrices, one column per estimate row:
#' \describe{
#'   \item{\code{cluster}}{Cluster-level random effects
#'     (\eqn{\gamma_j \cdot \tau_{between}}), representing between-cluster
#'     deviations from the fixed effects.}
#'   \item{\code{estimate}}{Estimate-level random effects
#'     (\eqn{\theta_i - \mu_i - \gamma_j \cdot \tau_{between}}),
#'     representing within-cluster deviations from the cluster means.}
#' }
#'
#' For \code{brma.mv()} random-formula models, decomposes by random-effect
#' block. Sampled blocks use fitted latent random effects; blocks compiled as
#' marginalized use Gaussian BLUP means. If there is only one block and
#' \code{simplify = TRUE}, returns a single \code{brma_samples} object.
#'
#' @return A \code{brma_samples} object for a single selected/simplified
#' component, or a named list of \code{brma_samples} objects for decomposed
#' \code{component = "all"} output.
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
                       simplify = TRUE, ...) {

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
      ...
    ))
  }

  # get BLUPs (fixed + all random effects)
  blup_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "estimate",
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
    colnames(cluster_ranef_mat) <- paste0("u_cluster[", cluster_names, "]")

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
      labels     = labels,
      n_chains   = n_chains,
      n_iter     = n_iter,
      probs      = probs,
      data       = data
    ))
  }
}


.ranef_brma_mv_random <- function(object, bias_adjusted = FALSE,
                                  probs = c(.025, .975), component = "all",
                                  simplify = TRUE, ...) {

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

  sampled_components <- .evaluate.brma.random_effects_components(
    object            = object,
    data              = data,
    posterior_samples = posterior_samples,
    same_data         = TRUE
  )
  sampled_total <- .sum_random_effect_components(
    components = sampled_components,
    S          = nrow(posterior_samples),
    K          = K
  )

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

  marginalized_components <- .evaluate.brma.mv_marginalized_random_blup.norm(
    object                 = object,
    mu_samples             = unclass(terms_samples),
    sampled_random_samples = sampled_total,
    posterior_samples      = posterior_samples,
    bias_offset            = bias_offset
  )

  components <- c(sampled_components, marginalized_components)
  design     <- .fitted_formula_design(object, "mu", required = TRUE)
  block_order <- vapply(
    design[["random_effects"]],
    `[[`,
    character(1),
    "block_name"
  )
  block_order <- block_order[block_order %in% names(components)]
  components  <- components[block_order]

  out <- lapply(names(components), function(block) {
    mat <- components[[block]]
    colnames(mat) <- paste0("u_", block, "[", labels[seq_len(K)], "]")
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

  return(.select_ranef_mv_components(
    components = list(location = out),
    component  = component,
    simplify   = simplify,
    labels     = labels,
    n_chains   = n_chains,
    n_iter     = n_iter,
    probs      = probs,
    data       = data
  ))
}


.select_ranef_mv_components <- function(components, component = "all",
                                        simplify = TRUE, labels, n_chains,
                                        n_iter, probs, data) {

  BayesTools::check_bool(simplify, "simplify")
  if (is.null(component)) {
    component <- "all"
  }
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    stop("'component' must be a single component name.", call. = FALSE)
  }

  component <- .normalize_ranef_mv_component(component)
  flat      <- .flatten_ranef_mv_components(components)

  if (component == "all") {
    if (isTRUE(simplify) && length(flat) == 1L) {
      return(flat[[1L]])
    }
    return(components)
  }

  if (component == "total") {
    total <- Reduce(`+`, lapply(flat, as.matrix))
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

  if (component %in% names(components)) {
    selected <- components[[component]]
    if (isTRUE(simplify) && length(selected) == 1L) {
      return(selected[[1L]])
    }
    return(selected)
  }

  block_matches <- which(names(flat) == component)
  if (length(block_matches) == 1L) {
    return(flat[[block_matches]])
  }
  if (length(block_matches) > 1L) {
    stop(
      "Random-effect component '", component, "' is ambiguous across ",
      "model components. Select the parent component first.",
      call. = FALSE
    )
  }

  stop(
    "Unknown random-effect component '", component, "'. Available components: ",
    paste(unique(c(names(components), names(flat), "all", "total")),
          collapse = ", "),
    ".",
    call. = FALSE
  )
}


.normalize_ranef_mv_component <- function(component) {

  if (component %in% c("mods", "mu")) {
    return("location")
  }

  component
}


.flatten_ranef_mv_components <- function(components) {

  flat <- unlist(components, recursive = FALSE, use.names = FALSE)
  names(flat) <- unlist(
    lapply(components, names),
    recursive = FALSE,
    use.names = FALSE
  )
  flat
}


.sum_random_effect_components <- function(components, S, K) {

  if (length(components) == 0L) {
    return(matrix(0, nrow = S, ncol = K))
  }

  Reduce(`+`, components)
}


.select_ranef_components <- function(components, component = "all",
                                     simplify = TRUE, labels, n_chains,
                                     n_iter, probs, data) {

  BayesTools::check_bool(simplify, "simplify")
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
    return(components)
  }

  if (component == "total") {
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
