# ============================================================================ #
# brma.wrappers.R
# ============================================================================ #
#
# This file contains user-friendly wrapper functions for common prediction
# tasks with brma model fits. These wrappers simplify the interface to
# predict.brma() for typical use cases:
#
# - pooled_effect.brma(): Aggregated pooled effect size (mu)
# - pooled_heterogeneity.brma(): Aggregated pooled heterogeneity (tau)
# - blup.brma(): Best Linear Unbiased Predictions (true study effects)
#
# Design principles:
# - Simple interface: minimal required arguments
# - Consistent with predict.brma: use same underlying machinery
# - Quiet operation: suppress aggregation messages for cleaner output
#
# ============================================================================ #


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
#' \code{TRUE}, which returns bias-corrected estimates. Set to \code{FALSE}
#' to obtain estimates that include publication bias effects.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' This function is a convenience wrapper around \code{predict.brma(...,
#' type = "terms", newdata = TRUE)}.
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
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get pooled effect
#' pooled_effect(fit)
#' }
#'
#' @seealso [predict.brma()], [pooled_heterogeneity()], [blup()]
#' @export
pooled_effect.brma <- function(object, bias_adjusted = TRUE,
                               probs = c(.025, .975), ...) {
  out <- predict.brma(
    object        = object,
    newdata       = TRUE,
    type          = "terms",
    probs         = probs,
    bias_adjusted = bias_adjusted,
    quiet         = TRUE,
    ...
  )
  attr(out, "title") <- "Pooled Effect Size"
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
#' from a fitted brma object by averaging across the scale model matrix.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' This function is a convenience wrapper around \code{predict.brma(...,
#' type = "terms.scale", newdata = TRUE)}.
#'
#' For location-scale models (with scale regression), the pooled heterogeneity
#' averages tau across the scale model matrix proportionately to the levels
#' observed in the data.
#'
#' For models without scale regression, this returns the single tau parameter.
#'
#' For multilevel (3-level) models, the returned tau is the total heterogeneity:
#' \code{tau = sqrt(tau_within^2 + tau_between^2)}.
#'
#' @return A \code{brma_samples} object containing posterior samples. When printed,
#' displays a summary table. Use \code{summary()} to obtain the summary table directly.
#' The samples can be converted to \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get pooled heterogeneity
#' pooled_heterogeneity(fit)
#' }
#'
#' @seealso [predict.brma()], [pooled_effect()], [blup()]
#' @export
pooled_heterogeneity.brma <- function(object, probs = c(.025, .975), ...) {
  out <- predict.brma(
    object     = object,
    newdata    = TRUE,
    type       = "terms.scale",
    probs      = probs,
    quiet      = TRUE,
    ...
  )
  attr(out, "title") <- "Pooled Heterogeneity"
  return(out)
}


# ---------------------------------------------------------------------------- #
# blup generic and brma method
# ---------------------------------------------------------------------------- #

#' @title Best Linear Unbiased Predictions (BLUPs)
#'
#' @description Computes the estimated true study effects (theta) from a
#' fitted model. These correspond to Best Linear Unbiased Predictions (BLUPs)
#' or empirical Bayes estimates.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of the true study effects.
#'
#' @seealso [predict.brma()]
#' @export
blup <- function(object, ...) {
  UseMethod("blup")
}


#' @title Best Linear Unbiased Predictions for brma Objects
#'
#' @description Computes the estimated true study effects (theta) for a
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
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' This function is a convenience wrapper around \code{predict.brma(...,
#' type = "effect", newdata = NULL)}.
#'
#' For normal models, true effects are computed using empirical Bayes shrinkage:
#' \deqn{\theta_i = \lambda_i \cdot y_i + (1 - \lambda_i) \cdot \mu_i}
#' where \eqn{\lambda_i = \tau^2 / (\tau^2 + se_i^2)}.
#'
#' For GLMM models (binomial, Poisson), the estimate-level random effects
#' are extracted directly from the posterior samples.
#'
#' For multilevel (3-level) models, the true effects incorporate both
#' study-level (\eqn{\gamma}) and estimate-level random effects.
#'
#' @return A \code{brma_samples} object containing posterior samples with one
#' column per study. When printed, displays a summary table. Use \code{summary()}
#' to obtain the summary table directly. The samples can be converted to
#' \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get BLUPs (true study effects)
#' blup(fit)
#' }
#'
#' @seealso [predict.brma()], [pooled_effect()], [pooled_heterogeneity()],
#' [true_effects()]
#' @export
blup.brma <- function(object, bias_adjusted = FALSE,
                      probs = c(.025, .975), ...) {
  out <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "effect",
    probs         = probs,
    bias_adjusted = bias_adjusted,
    quiet         = TRUE,
    ...
  )
  attr(out, "title") <- "True Study Effects (BLUPs)"
  return(out)
}


# ---------------------------------------------------------------------------- #
# true_effects generic and brma method (alias for blup)
# ---------------------------------------------------------------------------- #

#' @title True Study Effects
#'
#' @description Computes the estimated true study effects (theta) from a
#' fitted model. This is an alias for \code{\link{blup}}.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of the true study effects.
#'
#' @seealso [blup()], [predict.brma()]
#' @export
true_effects <- function(object, ...) {
  UseMethod("true_effects")
}


#' @title True Study Effects for brma Objects
#'
#' @description Computes the estimated true study effects (theta) for a
#' fitted brma object. This is an alias for \code{\link{blup.brma}}.
#'
#' @inheritParams blup.brma
#'
#' @details
#' This function is identical to \code{\link{blup.brma}}. See that function
#' for full details on how true effects are computed.
#'
#' @return A \code{brma_samples} object containing posterior samples with one
#' column per study. When printed, displays a summary table. Use \code{summary()}
#' to obtain the summary table directly. The samples can be converted to
#' \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get true effects (equivalent to blup())
#' true_effects(fit)
#' }
#'
#' @seealso [blup.brma()], [predict.brma()], [pooled_effect()],
#' [pooled_heterogeneity()]
#' @export
true_effects.brma <- function(object, bias_adjusted = FALSE,
                              probs = c(.025, .975), ...) {
  blup.brma(
    object        = object,
    bias_adjusted = bias_adjusted,
    probs         = probs,
    ...
  )
}
