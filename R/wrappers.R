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
#' @inheritParams predict.brma
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
                               output_measure = NULL, transform = NULL,
                               probs = c(.025, .975), ...) {
  out <- predict.brma(
    object         = object,
    newdata        = TRUE,
    type           = "terms",
    output_measure = output_measure,
    transform      = transform,
    probs          = probs,
    bias_adjusted  = bias_adjusted,
    quiet          = TRUE,
    ...
  )
  attr(out, "title") <- .effect_output_title(
    title            = "Pooled Effect Size",
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
#' @description Computes the estimated true effects (theta) from a
#' fitted model. These correspond to Best Linear Unbiased Predictions (BLUPs)
#' or empirical Bayes estimates.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of the true effects.
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
#' cluster-level (\eqn{\gamma}) and estimate-level random effects.
#'
#' @return A \code{brma_samples} object containing posterior samples with one
#' column per estimate. When printed, displays a summary table. Use \code{summary()}
#' to obtain the summary table directly. The samples can be converted to
#' \pkg{posterior} draws formats using \code{as_draws()}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get BLUPs (true effects)
#' blup(fit)
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
    title            = "True Effects (BLUPs)",
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
#' fitted model. This is an alias for \code{\link{blup}}.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value, typically a summary table or
#' posterior samples of the true effects.
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
#' @return A \code{brma_samples} object containing posterior samples with one
#' column per estimate. When printed, displays a summary table. Use \code{summary()}
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
#' These are the offsets from the fixed-effect predictions, corresponding
#' to what \code{metafor::ranef()} returns.
#'
#' @param object a fitted brma object
#' @param bias_adjusted whether to adjust for publication bias. Defaults to
#' \code{FALSE}. See \code{\link{blup.brma}} for details.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' Random effects are computed as the difference between the true
#' effects (BLUPs) and the fixed-effect predictions:
#' \deqn{u_i = \hat{\theta}_i - \hat{\mu}_i}
#'
#' For standard (2-level) models, returns a single \code{brma_samples}
#' object with the estimate-level random effects.
#'
#' For multilevel (3-level) models, returns a list with two components:
#' \describe{
#'   \item{\code{cluster}}{Cluster-level random effects
#'     (\eqn{\gamma_j \cdot \tau_{between}}), representing between-cluster
#'     deviations from the fixed effects.}
#'   \item{\code{estimate}}{Estimate-level random effects
#'     (\eqn{\theta_i - \mu_i - \gamma_j \cdot \tau_{between}}),
#'     representing within-cluster deviations from the cluster means.}
#' }
#'
#' @return For 2-level models, a \code{brma_samples} object. For 3-level
#' models, a named list of \code{brma_samples} objects (one per variance
#' component).
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # extract random effects (deviations from mu)
#' ranef(fit)
#' }
#'
#' @seealso [blup.brma()], [predict.brma()], [pooled_effect()]
#' @export
ranef.brma <- function(object, bias_adjusted = FALSE,
                       probs = c(.025, .975), ...) {

  dots <- list(...)
  if (!is.null(dots[["output_measure"]]) || !is.null(dots[["transform"]])) {
    stop(
      "'output_measure' and 'transform' are not available for random-effect ",
      "deviations.",
      call. = FALSE
    )
  }

  is_multilevel <- .is_multilevel(object)

  # extract MCMC chain info
  n_chains <- length(object[["fit"]][["mcmc"]])
  n_iter   <- object[["fit"]][["sample"]]
  data     <- object[["data"]]
  labels   <- .get_estimate_labels(object)

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

    return(.new_brma_samples(
      samples  = ranef_mat,
      n_chains = n_chains,
      n_iter   = n_iter,
      title    = "Random Effects:",
      probs    = probs,
      data     = data
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
    return(out)
  }
}
