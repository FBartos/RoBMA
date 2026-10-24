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

  coefficients <- object[["coefficients"]]
  if (is.null(coefficients)) {
    return(NULL)
  }
  names(coefficients) <- .location_repair_intercept_labels(
    labels = names(coefficients),
    object = object
  )

  return(coefficients)
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
  .warn_unused_dots(dots, allowed = "quiet", caller = "fitted()")
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
#' the pooled location effect. This reference declares no new-publication
#' standard error or selection event; explicit unadjusted \code{predict()} rows
#' instead target a selected new publication. Selection-model weighting affects response
#' predictions, not this fixed-location posterior summary.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param conditional whether to return the pooled effect conditional on the
#' effect component for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @inheritParams predict.brma
#' @param ... additional arguments passed to \code{\link{predict.brma}}
#' (\code{as_measure}, \code{V_new}, and \code{conditioning_depth}); wrapper
#' arguments such as \code{newdata}, \code{type}, and \code{quiet} are fixed,
#' and other arguments are unused and trigger a warning.
#'
#' @details
#' This function is a convenience wrapper around \code{predict.brma(...,
#' type = "terms", newdata = TRUE, bias_adjusted = TRUE, quiet = TRUE)} that
#' adds the prediction interval columns.
#'
#' The fixed-location model is evaluated at the average expanded
#' moderator design. Because the location predictor is linear in its design,
#' this equals averaging its fitted-design posterior draws directly.
#'
#' For meta-regression models, the pooled effect averages the effect size
#' estimate across moderator levels proportionately to the levels observed
#' in the data. This provides an estimate representative of the sample of
#' studies.
#'
#' For \code{brma.mv()} models and models with factor-specific or random-slope
#' terms, the average expanded design can be synthetic and need not correspond
#' to a single realizable outcome or study. The pooled effect and its prediction
#' interval then describe that linear composite and should be interpreted only
#' when the composite has a meaningful scientific interpretation. Use
#' \code{predict(..., newdata = ...)} for a specified realizable design.
#'
#' For models without moderators, this returns the single mu parameter. The
#' prediction interval is marginal over one new true effect at the same average
#' location, scale, and random-effect design. It therefore uses
#' \code{pooled_heterogeneity(component = "total")}, rather than the observed-design RMS heterogeneity
#' reported by \code{summary_heterogeneity()}. For known-\eqn{R} random
#' effects this is a unit-kernel reference interval; use \code{predict()} for
#' an interval tied to a fitted group and its known covariance multiplier.
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

  forwarded <- .split_aggregate_dots(list(...))

  # the aggregated prediction at the average fitted design, with the draws of
  # one new true effect at the pooled reference (the PI columns)
  return(do.call(.predict_brma_aggregate, c(
    list(
      object              = object,
      type                = "terms",
      probs               = probs,
      bias_adjusted       = bias_adjusted,
      output_measure      = output_measure,
      transform           = transform,
      conditional         = conditional,
      quiet               = TRUE,
      caller              = "pooled_effect.brma()",
      dots                = forwarded[["unused"]],
      title               = "Pooled Effect Size",
      prediction_interval = TRUE
    ),
    forwarded[["predict"]]
  )))
}


# Split the dots of a pooled_*() wrapper: the arguments of predict() that the
# wrapper does not fix are forwarded to it (as in 4.0.0), the posterior-sample
# reuse argument is internal, and everything else is unused.
.split_aggregate_dots <- function(dots) {

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }
  predict_arguments <- c("as_measure", "V_new", "conditioning_depth")

  return(list(
    predict = dots[dot_names %in% predict_arguments],
    unused  = dots[!dot_names %in% predict_arguments]
  ))
}


# Fixed location at the average fitted design: the row mean of the fitted-design
# location draws (S x 1, column "mu").
.pooled_location_samples <- function(object, bias_adjusted, posterior_samples) {

  data   <- object[["data"]]
  priors <- object[["priors"]]
  mu_samples <- .evaluate.brma.mu(
    fit               = .object_formula_fit(object),
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

  return(samples)
}


# ---------------------------------------------------------------------------- #
# pooled_heterogeneity generic and brma method
# ---------------------------------------------------------------------------- #

#' @title Pooled Heterogeneity
#'
#' @description Computes method-specific pooled (aggregated) heterogeneity from
#' a fitted model.
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
#' @param ... additional arguments passed to \code{\link{predict.brma}}
#' (\code{as_measure}, \code{V_new}, and \code{conditioning_depth}); wrapper
#' arguments such as \code{newdata}, \code{type}, and \code{quiet} are fixed,
#' and other arguments are unused and trigger a warning.
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
#' For \code{brma.mv()} models, one ordinary selected component is named
#' \code{tau}; a genuine variance-additive total is named \code{tau_total}; and a
#' mean-variance allocation anchor is named \code{tau_common}. These quantities
#' are available when a random formula is declared. For random-formula models,
#' \code{component = "all"}
#' returns one \code{brma_samples} object when there is a single heterogeneity
#' component and a named list when there are multiple components.
#' \code{component = "total"} computes the variance-additive total at the
#' average random-effect design. Random slopes are evaluated as
#' \eqn{\bar q^T G \bar q}, including within-block correlations. Default and
#' named allocation components retain their declared scale anchors; explicit
#' totals and pooled-effect prediction intervals use the full covariance of
#' the average design, including allocation-owned blocks.
#' With a known group covariance \eqn{R}, the
#' reported \code{sd} is the fitted multiplier of that kernel. It need not equal
#' every row's marginal standard deviation when \eqn{\mathrm{diag}(R)} is not
#' one. Row-specific marginal standard deviations are available from
#' \code{predict(type = "terms.scale")}.
#'
#' For \code{brma.mv()} models and models with factor-specific or random-slope
#' terms, the average expanded random-effect design can be synthetic and need
#' not correspond to a single realizable outcome or study. The result is the
#' heterogeneity of that linear composite, not typical heterogeneity across the
#' observed designs. Interpret it only when the composite has a meaningful
#' scientific interpretation; use \code{summary_heterogeneity()} for
#' observed-design heterogeneity or \code{predict(..., newdata = ...)} for a
#' specified realizable design.
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

  forwarded <- .split_aggregate_dots(list(...))
  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' predictions are available only for RoBMA objects.",
         call. = FALSE)
  }
  if (inherits(object, "brma.mv")) {
    .warn_unused_dots(
      dots    = forwarded[["unused"]],
      allowed = ".posterior_samples",
      caller  = "pooled_heterogeneity.brma()"
    )
    return(
      .pooled_heterogeneity_brma_mv(
        object             = object,
        probs              = probs,
        conditional        = conditional,
        component          = component,
        .posterior_samples = forwarded[["unused"]][[".posterior_samples"]]
      )
    )
  }

  .check_univariate_heterogeneity_component(component)

  return(do.call(.predict_brma_aggregate, c(
    list(
      object                  = object,
      type                    = "terms.scale",
      probs                   = probs,
      bias_adjusted           = FALSE,
      output_measure          = NULL,
      transform               = NULL,
      conditional             = conditional,
      quiet                   = TRUE,
      caller                  = "pooled_heterogeneity.brma()",
      dots                    = forwarded[["unused"]],
      title                   = "Pooled Heterogeneity",
      heterogeneity_component = "heterogeneity"
    ),
    forwarded[["predict"]]
  )))
}


.pooled_heterogeneity_total_samples <- function(object, posterior_samples) {

  if (inherits(object, "brma.mv")) {
    components <- .pooled_brma_mv_heterogeneity_components(
      object            = object,
      posterior_samples = posterior_samples,
      covariance        = TRUE
    )
    return(.total_brma_mv_heterogeneity_samples(components))
  }

  return(.pooled_tau_components(object, posterior_samples)[["tau_total"]])
}

# Total, within-, and between-cluster heterogeneity at the average fitted scale
# design of a univariate model.
.pooled_tau_components <- function(object, posterior_samples) {

  tau_result <- .evaluate.brma.pooled_tau(
    fit               = .object_formula_fit(object),
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) {
      .create_fit_formula_list(data = object[["data"]], "scale")
    } else {
      NULL
    },
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = .is_multilevel(object),
    is_random         = .is_random(object),
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )

  return(tau_result)
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
#' For known-\code{V} models, the conditional means require applying the
#' inverse of the marginal covariance. A covariance that is too ill-conditioned
#' for that (a reciprocal condition number of about \code{2e-10} or less in any
#' posterior draw used) is unavailable: the call stops with an error of class
#' \code{RoBMA_known_v_precision_unavailable}, whose parent class is
#' \code{RoBMA_numerical_unavailable}, instead of returning a rounded number.
#' The condition carries the draw positions (\code{rows}) and their
#' \code{reciprocal_condition}.
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
#' observation. Defaults to \code{TRUE} (one column per estimate; the cluster
#' block of a multilevel model labels its columns
#' \code{u_cluster[<cluster>|<estimate>]}), except for \code{brma.mv()}
#' random-formula models, where it defaults to \code{FALSE}. \code{FALSE}
#' returns one column per unique
#' grouping-level contribution in the fitted grouping order used by
#' \code{metafor::ranef()} (\code{u_cluster[<cluster>]} for a multilevel
#' model). Nested groups follow their original outer-to-inner
#' order. Indicator-coded random coefficients retain one column per observed
#' grouping-level and coefficient combination. Set to \code{TRUE} for
#' observation-aligned output.
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
#' coefficient combination. Other random slopes generally require
#' \code{expand = TRUE}. When the random-effect standard deviation of a block
#' varies across the rows of one returned cell (a grouping level, or a grouping
#' level and indicator coefficient) through a \code{scale} formula
#' (also for the cluster component of a clustered model with a \code{scale}
#' formula, whether normal, GLMM, or with a weight function), the unique-level
#' output is refused with an error of class
#' \code{RoBMA_ranef_expand_required} (also of class
#' \code{RoBMA_unsupported_request}); use \code{expand = TRUE}. Only the
#' components that are returned are judged: \code{component = "estimate"} of a
#' clustered model and a block name whose own scale does not vary are returned,
#' whereas \code{"all"} and \code{"total"} include every block. Summed
#' \code{component = "total"} output across multiple blocks also requires
#' \code{expand = TRUE} because different blocks need not share grouping levels.
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
                       simplify = TRUE, expand = TRUE, ...) {

  expand_default <- missing(expand)
  BayesTools::check_bool(expand, "expand")
  dots <- list(...)
  controlled <- intersect(names(dots), c("newdata", "type", "quiet"))
  if (length(controlled) > 0L) {
    # the arguments the method fixes are unused: warn and continue without them
    .warn_unused_dots(dots[controlled], allowed = character(), caller = "ranef()")
    return(do.call(ranef.brma, c(
      list(
        object = object, bias_adjusted = bias_adjusted, probs = probs,
        component = component, simplify = simplify
      ),
      if (!expand_default) list(expand = expand),
      dots[setdiff(names(dots), controlled)]
    )))
  }
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
      expand        = if (expand_default) FALSE else expand,
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
      samples   = ranef_mat,
      n_chains  = n_chains,
      n_iter    = n_iter,
      title     = "Random Effects:",
      component = "random",
      probs     = probs,
      data      = data
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
    cluster_labels <- .ranef_cluster_labels(object, cluster)
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
      cluster_level_labels <- .ranef_cluster_labels(object, cluster_levels)
      .check_ranef_cluster_unique_level(object, component, cluster_level_map)
      cluster_ranef_mat <- .ranef_unique_level_samples(
        samples      = cluster_ranef_mat,
        group_map    = cluster_level_map,
        group_levels = cluster_level_labels,
        block        = "cluster"
      )
    }

    cluster_ranef <- .new_brma_samples(
      samples   = cluster_ranef_mat,
      n_chains  = n_chains,
      n_iter    = n_iter,
      title     = "Cluster-Level Random Effects:",
      component = "random/cluster",
      probs     = probs,
      data      = data
    )

    # estimate-level random effects: BLUP - cluster predictions
    estimate_ranef_mat <- unclass(blup_samples) - unclass(cluster_samples)
    colnames(estimate_ranef_mat) <- paste0("u_estimate[", labels[seq_len(K)], "]")

    estimate_ranef <- .new_brma_samples(
      samples   = estimate_ranef_mat,
      n_chains  = n_chains,
      n_iter    = n_iter,
      title     = "Estimate-Level Random Effects:",
      component = "random/estimate",
      probs     = probs,
      data      = data
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

  .warn_unused_dots(
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
  cluster_labels <- .ranef_cluster_labels(object, cluster)
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
    cluster_level_labels <- .ranef_cluster_labels(object, cluster_levels)
    .check_ranef_cluster_unique_level(object, component, cluster_level_map)
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
      samples   = components[["cluster"]],
      n_chains  = context[["n_chains"]],
      n_iter    = context[["n_iter"]],
      title     = "Cluster-Level Random Effects:",
      component = "random/cluster",
      probs     = probs,
      data      = context[["new_data"]]
    ),
    estimate = .new_brma_samples(
      samples   = components[["estimate"]],
      n_chains  = context[["n_chains"]],
      n_iter    = context[["n_iter"]],
      title     = "Estimate-Level Random Effects:",
      component = "random/estimate",
      probs     = probs,
      data      = context[["new_data"]]
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
  dots[[".posterior_samples"]] <- posterior_samples

  # Terms, bias offsets and the random posterior are computed on every
  # posterior row; the conditional row selection is applied to the
  # decomposition last, as the prediction methods do.
  conditional <- dots[["conditional"]]
  if (is.null(conditional)) {
    conditional <- FALSE
  }
  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' predictions are available only for RoBMA objects.",
         call. = FALSE)
  }
  dots[["conditional"]] <- NULL

  terms_samples <- do.call(
    predict.brma,
    c(
      list(
        object        = object,
        newdata       = NULL,
        type          = "terms",
        probs         = probs,
        bias_adjusted = bias_adjusted,
        quiet         = TRUE
      ),
      dots
    )
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

  components <- .predict_brma_mv_random_posterior(
    object            = object,
    mu_samples        = unclass(terms_samples),
    posterior_samples = posterior_samples,
    bias_offset       = bias_offset,
    type              = "mean",
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

  # Only the blocks that the requested component returns are converted to
  # unique levels and judged. A component that returns no block is unknown (or
  # malformed) and is refused by the selection before any block is judged.
  returned <- vapply(
    block_order, .ranef_block_returned, logical(1L), component = component
  )
  if (!any(returned)) {
    return(.select_ranef_components(
      components = components[block_order],
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
  block_order <- block_order[returned]
  components  <- components[block_order]

  out <- lapply(names(components), function(block) {
    mat <- components[[block]]
    if (expand) {
      colnames(mat) <- paste0("u_", block, "[", labels[seq_len(K)], "]")
    } else {
      term     <- random_terms[[block]]
      metadata <- .ranef_unique_level_term(term, block)
      if (.ranef_block_sd_varies_within_cells(
          object, term, metadata[["group_map"]])) {
        .stop_ranef_expand_required(
          "Unique-level random effects are unavailable for random-effect ",
          "block '", block, "' because its standard deviation varies across ",
          "the rows of a grouping level. Use 'expand = TRUE'."
        )
      }
      mat <- .ranef_unique_level_samples(
        samples      = mat,
        group_map    = metadata[["group_map"]],
        group_levels = metadata[["group_levels"]],
        block        = block
      )
    }
    .new_brma_samples(
      samples   = mat,
      n_chains  = n_chains,
      n_iter    = n_iter,
      title     = paste0("Random Effects: ", block),
      component = paste("random", block, sep = "/"),
      probs     = probs,
      data      = data
    )
  })
  names(out) <- names(components)

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
    posterior_samples = posterior_samples,
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
    if (expand) {
      colnames(total) <- paste0("u[", labels[seq_len(ncol(total))], "]")
    } else {
      # Only names carrying the block prefix (u_<block>[label]) are stripped;
      # names that are already u[label] stay as they are.
      prefix   <- paste0("u_", names(components)[[1L]])
      prefixed <- startsWith(colnames(total), prefix)
      colnames(total)[prefixed] <- paste0(
        "u", substring(colnames(total)[prefixed], nchar(prefix) + 1L)
      )
    }
    return(.new_brma_samples(
      samples   = total,
      n_chains  = n_chains,
      n_iter    = n_iter,
      title     = "Random Effects:",
      component = "random",
      probs     = probs,
      data      = data
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


.ranef_cluster_labels <- function(object, cluster) {

  labels  <- unname(.get_cluster_labels(object)[as.character(cluster)])
  missing <- is.na(labels)
  labels[missing] <- as.character(cluster[missing])

  return(labels)
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

  group_order <- seq_len(n_groups)
  level_rows  <- match(group_order, group_map)
  if (anyNA(level_rows)) {
    stop("Internal error: random-effect grouping levels are missing fitted rows.",
         call. = FALSE)
  }

  out <- samples[, level_rows, drop = FALSE]
  colnames(out) <- paste0(
    "u_", block, "[", group_levels[group_order], "]"
  )
  return(out)
}


# Restrict a declared grouping factor to the levels that own fitted rows, in
# declared order. A declared level without rows has no random effect to report.
.ranef_observed_group_levels <- function(term) {

  group_map    <- as.integer(term[["group_map"]])
  group_levels <- as.character(term[["group_levels"]])
  observed     <- term[["group_observed_levels"]]
  keep <- if (is.null(observed)) {
    seq_along(group_levels) %in% group_map
  } else {
    group_levels %in% as.character(observed)
  }

  return(list(
    group_map    = match(group_map, which(keep)),
    group_levels = group_levels[keep]
  ))
}


# Unique-level output keeps one column per returned cell (a grouping level, or
# a grouping level and indicator coefficient). It is refused with an error of
# class "RoBMA_ranef_expand_required" whenever the random-effect SD varies
# across the rows mapped to one returned cell, because the first row would
# silently stand for all of them.
.stop_ranef_expand_required <- function(...) {

  message <- paste(unlist(lapply(list(...), as.character)), collapse = "")
  stop(structure(
    class = c("RoBMA_ranef_expand_required", "RoBMA_unsupported_request",
              "error", "condition"),
    list(message = message, call = NULL)
  ))
}


# Whether 'component' returns the output of the block 'block': every block for
# "all" and "total", the block itself for its name. The unique-level conversion
# is judged only for returned blocks; an unknown component returns none.
.ranef_block_returned <- function(block, component) {

  if (is.null(component)) {
    return(TRUE)
  }

  return(is.character(component) && length(component) == 1L &&
           !is.na(component) && component %in% c("all", "total", block))
}


# Refuse the unique-level cluster output of a clustered model whose scale design
# varies across the rows of a cluster, when the cluster component is returned.
.check_ranef_cluster_unique_level <- function(object, component, cluster_map) {

  if (.ranef_block_returned("cluster", component) && .is_scale(object) &&
      .ranef_scale_design_varies_within_cells(
        object, .data_scale_formula_parameters(object[["data"]]), cluster_map)) {
    .stop_ranef_expand_required(
      "Unique-level random effects are unavailable for the cluster ",
      "component because the 'scale' design varies across the rows of a ",
      "cluster. Use 'expand = TRUE'."
    )
  }

  return(invisible(TRUE))
}


# Whether the design rows of the scale parameters that feed a block's
# row-shaped SD sources differ among the rows mapped to one returned output cell
# ('cell_map', as returned by .ranef_unique_level_term()). Structural metadata
# only; posterior draws are never consulted.
.ranef_block_sd_varies_within_cells <- function(object, term, cell_map) {

  binding <- term[["sd_binding"]]
  sources <- c(
    list(binding[["source"]]),
    if (is.list(binding[["sources_by_column"]])) binding[["sources_by_column"]],
    if (is.list(binding[["allocations"]])) lapply(binding[["allocations"]], `[[`, "source")
  )
  is_row <- vapply(sources, function(source) {
    !is.null(source) && identical(source[["shape"]], "row")
  }, logical(1L))
  if (!any(is_row)) {
    return(FALSE)
  }

  # Row-shaped sources are the scale components; an unmapped source cannot be
  # inspected and is refused.
  specs       <- .data_scale_component_specs(object[["data"]])
  parameters  <- stats::setNames(
    vapply(specs, `[[`, character(1L), "parameter"),
    vapply(specs, `[[`, character(1L), "source")
  )
  source_names <- unique(vapply(sources[is_row], `[[`, character(1L), "name"))
  mapped       <- parameters[source_names]
  if (anyNA(mapped)) {
    return(anyDuplicated(as.integer(cell_map)) > 0L)
  }

  return(.ranef_scale_design_varies_within_cells(
    object, unname(mapped), cell_map
  ))
}


# TRUE when the design of any listed scale parameter differs (exactly) between
# the rows of one cell, or cannot be inspected for a cell with several rows.
.ranef_scale_design_varies_within_cells <- function(object, parameters,
                                                    cell_map) {

  cell_map <- as.integer(cell_map)
  first    <- match(cell_map, cell_map)
  if (all(first == seq_along(cell_map))) {
    return(FALSE)
  }

  for (parameter in parameters) {
    design <- .fitted_formula_design(object, parameter, required = FALSE)
    rows   <- design[["model_matrix"]]
    if (!is.matrix(rows) || nrow(rows) != length(cell_map)) {
      return(TRUE)
    }
    if (any(rows != rows[first, , drop = FALSE])) {
      return(TRUE)
    }
  }

  return(FALSE)
}


# Resolve the unique random-effect cells represented by a fitted block.
.ranef_unique_level_term <- function(term, block) {

  model_matrix <- term[["model_matrix"]]
  is_intercept <- is.matrix(model_matrix) && is.numeric(model_matrix) &&
    ncol(model_matrix) == 1L && all(is.finite(model_matrix)) &&
    all(model_matrix[, 1L] == 1)
  if (is_intercept) {
    return(.ranef_observed_group_levels(term))
  }

  is_indicator <- is.matrix(model_matrix) && is.numeric(model_matrix) &&
    ncol(model_matrix) > 1L && all(is.finite(model_matrix)) &&
    all(model_matrix == 0 | model_matrix == 1) &&
    all(rowSums(model_matrix) == 1)
  if (!is_indicator) {
    .stop_ranef_expand_required(
      "Unique-level random effects are unavailable for random-slope block '",
      block, "'. Use 'expand = TRUE'."
    )
  }

  coefficient_labels <- if (identical(term[["contrast_owner"]], "structure")) {
    term[["column_names"]]
  } else {
    term[["sd_leaves"]][["leaf_terms_by_column"]]
  }
  if (!is.character(coefficient_labels) ||
      length(coefficient_labels) != ncol(model_matrix) ||
      anyNA(coefficient_labels) || any(!nzchar(coefficient_labels)) ||
      anyDuplicated(coefficient_labels)) {
    .stop_ranef_expand_required(
      "Unique-level random effects are unavailable for block '", block,
      "' because its semantic coefficient labels are missing. Use ",
      "'expand = TRUE'."
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
    .stop_ranef_expand_required(
      "'component = \"total\"' across multiple random-effect blocks requires ",
      "'expand = TRUE'."
    )
  }

  invisible(TRUE)
}
