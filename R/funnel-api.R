# ============================================================================ #
# funnel-api.R
# ============================================================================ #

### funnel plot functions ----
#' @export
funnel <- function(x, ...) UseMethod("funnel")


#' @name funnel
#' @title Frequentist-Style and Bayesian Funnel Plots for brma Objects
#'
#' @description \code{funnel.brma} creates fast plug-in funnel contours.
#' Continuous parameters are fixed at their within-model posterior means. For
#' model-averaged fits, complete joint models retain their posterior model
#' probabilities and their plug-in sampling CDFs are averaged before the
#' contours are inverted. This makes the ordinary normal-model construction
#' directly comparable to a conventional plug-in funnel such as
#' \code{metafor::funnel(..., addtau2 = TRUE)}.
#'
#' \code{bfunnel.brma} instead averages the sampling CDF over posterior draws
#' before inversion. It therefore incorporates continuous parameter uncertainty
#' as well as model uncertainty.
#'
#' @param x a fitted brma object
#' @param residual whether to use residual mode. Defaults to not specified,
#' which means the function automatically determines the mode:
#' \itemize{
#'   \item not specified: Displays outcomes only for intercept-only models with
#'     one common marginal heterogeneity distribution. Otherwise it
#'     automatically uses LOO-PIT residual mode.
#'   \item \code{FALSE}: explicitly requests outcome mode. This is available
#'     only for intercept-only models with common marginal heterogeneity.
#'   \item \code{TRUE}: explicitly requests residual mode, displaying residuals
#'     on the x-axis and using \code{type} to determine how those residuals are
#'     computed.
#' }
#' @param type the type of residuals to use when in residual mode.
#' Options are:
#' \itemize{
#'   \item \code{"LOO-PIT"} (alias: \code{"rstudent"}; default): Leave-one-out
#'     probability integral transform residuals returned by
#'     \code{\link{rstudent.brma}}.
#'   \item \code{"rstandard"}: Internally standardized residuals using
#'     \code{\link{rstandard.brma}}. Only available for normal outcome models.
#'   \item \code{"outcome"}: Raw outcome residuals from
#'     \code{\link{residuals.brma}} with \code{type = "outcome"}. Its funnel
#'     bands are descriptive sampling-error-only reference bands.
#' }
#' Only used when funnel is in residual mode.
#' @param unit output unit for residual mode. Only \code{"estimate"} is
#' implemented in this pass.
#' @param conditioning_depth residual conditioning depth for residual mode.
#' Options are \code{"marginal"}, \code{"cluster"}, and \code{"estimate"}.
#' The default LOO-PIT residual path is available only with marginal
#' conditioning depth.
#' @param sampling_heterogeneity whether heterogeneity should be incorporated
#' into the sampling distribution funnel. Defaults to \code{TRUE}. Only used
#' in outcome mode and ignored in residual mode.
#' @param sampling_bias whether publication bias should be incorporated into the
#' sampling distribution funnel. Defaults to \code{TRUE}. Only used when
#' \code{residual = FALSE} or when automatic mode selects outcome mode. Ignored
#' in residual mode. When \code{TRUE} and the model
#' includes selection models (weightfunction), uses selected-normal quantiles.
#' When \code{TRUE} and the model includes PET/PEESE, incorporates the expected
#' skew from these regression adjustments.
#' @param max_samples maximum number of posterior draws used by
#' \code{bfunnel()}. Defaults to \code{10000}; use \code{Inf} for all draws.
#' \code{funnel()} uses all draws to estimate conditional plug-in values and
#' posterior model probabilities, then evaluates only one CDF per joint model.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional graphical arguments to customize the plot appearance:
#' \describe{
#'   \item{xlim, ylim}{numeric vectors of length 2 specifying axis limits}
#'   \item{xlab, ylab}{character strings for axis labels}
#'   \item{main}{character string for plot title (default: no title)}
#'   \item{pch}{point symbol (default: 21, filled circle). Use standard R pch values.}
#'   \item{col}{point border color (default: "black")}
#'   \item{bg}{point fill/background color (default: "#A6A6A6")}
#'   \item{cex}{point size multiplier for base graphics (default: 1)}
#'   \item{size}{point size for ggplot2 (default: 2)}
#'   \item{las}{axis-label style for base graphics (default: 1)}
#'   \item{back}{background region color (default: "grey"). Set to \code{NA} to suppress.}
#'   \item{shade}{funnel region fill color (default: "white"). Set to \code{NA} to suppress.}
#'   \item{lty}{line type for funnel edges and center line (default: "dotted")}
#'   \item{col.line}{color for funnel edge lines (default: "black")}
#'   \item{refline}{numeric override for the reference line. By default,
#'   residual mode uses 0, while outcome mode uses the center of the fitted
#'   sampling distribution, which may be curved when PET/PEESE bias adjustment
#'   is incorporated.}
#'   \item{col.refline}{color of vertical reference line (default: "black")}
#'   \item{as_data}{if \code{TRUE}, returns plot data instead of creating a plot}
#' }
#'
#' @details
#' The ordinary \code{funnel()} has outcome and residual modes. If
#' \code{residual} is not specified, outcomes are used only when the fitted
#' model is intercept-only and its marginal heterogeneity is identical across
#' studies. The latter check uses fitted row-marginal random-effect variances,
#' including complete \eqn{ZGZ'} contributions for \code{brma.mv()} random-
#' formula models. All other models use residual mode.
#'
#' \strong{Outcome mode} (intercept-only models with common heterogeneity):
#' Displays observed effect sizes on the x-axis and standard errors on the
#' y-axis. The reference line follows the center of the fitted sampling
#' distribution. When \code{sampling_bias = FALSE}, this center is the plug-in
#' mixture median (the pooled effect in a single model); when PET/PEESE bias
#' adjustment is incorporated, the center line can vary with the standard
#' error. The funnel region represents the central 95\%
#' region of the plug-in sampling distribution, optionally incorporating
#' heterogeneity and publication bias. Continuous parameters are replaced by
#' posterior means separately within every complete joint model. The resulting
#' model-specific CDFs are averaged with posterior model probabilities and only
#' then inverted; neither parameters nor contour endpoints are averaged across
#' models.
#' For normal models, these contours use the marginal response target:
#' heterogeneity represents newly realized latent effects and fitted study
#' effects are not retained.
#' For correlated known-\code{V} \code{brma.mv()} models, outcome-mode funnels
#' are descriptive scalar-SE displays based on the diagonal of \code{V}; residual
#' mode follows the fitted estimate-unit residual target.
#'
#' \code{bfunnel()} is available only for intercept-only normal models with
#' common marginal heterogeneity. It averages the same sampling CDF over
#' posterior draws and then inverts the averaged CDF. Thus it incorporates
#' posterior uncertainty in the effect, heterogeneity, PET/PEESE coefficients,
#' and selection parameters in addition to model uncertainty. For models
#' without a common outcome-scale target, use the default LOO-PIT residual
#' funnel from \code{funnel()}.
#'
#'
#' \strong{Residual mode} (models without a common outcome-scale target):
#' Displays residuals on the x-axis and standard errors on the
#' y-axis. The funnel region is bounded by the 0.025 and 0.975 normal quantiles
#' times the displayed standard error. It is a descriptive reference band, not
#' a posterior-predictive coverage interval. With \code{type = "LOO-PIT"}, the
#' plotted x-coordinate is the PIT-normalized \code{z} value returned by
#' \code{\link{rstudent.brma}} multiplied by its LOO predictive standard
#' deviation. Consequently, each plotted \code{x / se} is exactly its LOO-PIT
#' \code{z}. This reduces to the ordinary raw deleted residual when the LOO
#' predictive distribution is normal and extends the same funnel geometry to
#' non-normal continuous predictive mixtures. With
#' \code{type = "rstandard"}, the plotted values are internally standardized
#' residual companions from \code{\link{rstandard.brma}}. With
#' \code{type = "outcome"}, these are raw outcome residuals and the band uses
#' sampling error only; it does not incorporate fitted heterogeneity.
#'
#' The \code{type} argument controls how residuals are computed in residual
#' mode. See \code{\link{residuals.brma}} for details on each type.
#' The \code{sampling_heterogeneity} and \code{sampling_bias} arguments are
#' ignored in residual mode.
#'
#' For GLMM models, outcome-mode observed effect sizes are computed from the raw
#' frequency data using formulas equivalent to \code{metafor::escalc}. Their
#' contours use the corresponding continuity-corrected effect-size estimates
#' and approximate standard errors in a normal reference calculation. The
#' default x-axis label and a warning identify this as a descriptive
#' approximation; it is not a coverage interval from the fitted binomial or
#' Poisson likelihood and does not use the prior-nuisance marginal GLMM response
#' target from \code{predict()}. Exact discrete-likelihood coverage cannot be
#' indexed by a scalar standard error alone because it also depends on arm
#' totals or exposures and nuisance rates.
#' LOO-PIT residual mode is unavailable because the fitted predictive
#' distribution is discrete and no discrete PIT convention has been defined.
#' Use \code{type = "outcome"} for a descriptive effect-size-scale residual
#' funnel when needed. Its sampling-error-only bands are not exact diagnostics
#' for the fitted count likelihood.
#'
#' @return If \code{as_data = TRUE}, \code{funnel.brma} and
#' \code{bfunnel.brma} return a list with the
#' data used for plotting, including the plotted points, funnel polygons,
#' plotting limits, labels, reference line, and outcome-funnel estimand. GLMM
#' outcome-mode data also
#' contain an \code{approximation} record identifying the descriptive normal
#' effect-size approximation and stating that it is not fitted discrete-
#' likelihood coverage. Otherwise, it returns
#' \code{NULL} invisibly if \code{plot_type = "base"} or a ggplot object if
#' \code{plot_type = "ggplot"}.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE) &&
#'     requireNamespace("metafor", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'   dat <- metafor::escalc(
#'     measure = "RR",
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
#'     data    = dat.bcg
#'   )
#'
#'   fit <- brma(yi = yi, vi = vi, data = dat, measure = "RR")
#'   funnel(fit)
#'   bfunnel(fit)
#'   funnel(fit, pch = 19, col = "blue", bg = "lightblue")
#'
#'   fit_reg <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     mods    = ~ ablat + year,
#'     data    = dat,
#'     measure = "RR"
#'   )
#'   fit_reg <- add_loo(fit_reg)
#'   funnel(fit_reg)
#'   funnel(fit_reg, type = "outcome")
#'
#'   funnel_data <- funnel(fit, as_data = TRUE)
#'   funnel(fit, plot_type = "ggplot")
#' }
#' }
#'
#' @seealso [residuals.brma()], [rstandard.brma()], [rstudent.brma()],
#'   [predict.brma()], [zplot.brma()]
#' @aliases funnel
#' @export
#' @exportS3Method metafor::funnel
#' @rdname funnel
funnel.brma <- function(x, residual, type = "LOO-PIT",
                        unit = "estimate", conditioning_depth = "marginal",
                        sampling_heterogeneity = TRUE, sampling_bias = TRUE,
                        max_samples = 10000,
                        plot_type = "base", ...) {
  # input validation
  conditioning_depth_specified <- !missing(conditioning_depth)
  dots                         <- list(...)
  .check_legacy_level_arg(dots, "funnel()")

  BayesTools::check_bool(sampling_heterogeneity, "sampling_heterogeneity")
  BayesTools::check_bool(sampling_bias, "sampling_bias")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  max_samples <- .normalize_funnel_max_samples(max_samples)

  # set up graphical arguments with defaults
  dots <- .set_dots_funnel(dots)

  # get model characteristics
  is_mods              <- .is_mods(x)
  is_scale             <- .is_scale(x)
  common_heterogeneity <- NULL

  # determine mode: residual mode if mods/scale present, or if explicitly requested

  if (missing(residual)) {
    is_residual <- is_mods || is_scale
    if (!is_residual) {
      common_heterogeneity <- .funnel_common_heterogeneity(x)
      is_residual          <- !common_heterogeneity[["common"]]
    }
  } else {
    BayesTools::check_bool(residual, "residual")
    is_residual <- residual
  }

  if (!is_residual && (is_mods || is_scale)) {
    stop(
      "Outcome-mode funnel plots are not supported for models with location ",
      "or scale predictors because their fitted location or heterogeneity is ",
      "study-specific. Use 'residual = TRUE' (or omit 'residual') to plot a ",
      "residual funnel.",
      call. = FALSE
    )
  }

  if (!is_residual && is.null(common_heterogeneity)) {
    common_heterogeneity <- .funnel_common_heterogeneity(x)
  }
  if (!is_residual && !common_heterogeneity[["common"]]) {
    .funnel_stop_no_common_heterogeneity("funnel()")
  }

  # generate funnel data based on mode
  if (is_residual) {
    BayesTools::check_char(type, "type", allow_values = c("outcome", "rstandard", "LOO-PIT", "rstudent"))
    unit               <- .normalize_unit(unit)
    conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
    .check_unit_conditioning_depth(
      object             = x,
      unit               = unit,
      conditioning_depth = conditioning_depth,
      caller             = "funnel()"
    )

    if (unit == "cluster") {
      .check_cluster_unit_deferred("funnel()")
    }

    if ((type == "LOO-PIT" || type == "rstudent") && conditioning_depth_specified) {
      stop(
        "LOO-PIT residuals use the estimate-unit LOO target; ",
        "do not set 'conditioning_depth'.",
        call. = FALSE
      )
    }

    # sampling heterogeneity/bias ignored in residual mode - explicitly set to FALSE
    sampling_heterogeneity <- FALSE
    sampling_bias          <- FALSE
    funnel_data <- .funnel_data_residual(
      x                  = x,
      type               = type,
      unit               = unit,
      conditioning_depth = conditioning_depth,
      dots               = dots
    )
  } else {
    funnel_data <- .funnel_data_outcome(
      x                      = x,
      sampling_heterogeneity = sampling_heterogeneity,
      sampling_bias          = sampling_bias,
      max_samples            = max_samples,
      estimand               = "plugin",
      common_heterogeneity   = common_heterogeneity,
      dots                   = dots
    )
  }

  # allow data return for programmatic access
  if (isTRUE(dots[["as_data"]])) {
    return(funnel_data)
  }

  # create plot
  if (plot_type == "ggplot") {
    return(.funnel_plot_ggplot(funnel_data, dots))
  } else {
    .funnel_plot_base(funnel_data, dots)
    return(invisible(NULL))
  }
}


#' @rdname funnel
#' @export
bfunnel <- function(x, ...) UseMethod("bfunnel")


#' @rdname funnel
#' @export
bfunnel.default <- function(x, ...) {

  stop("'bfunnel' is available only for fitted brma objects.", call. = FALSE)
}


#' @rdname funnel
#' @export
bfunnel.brma <- function(x, sampling_heterogeneity = TRUE,
                         sampling_bias = TRUE, max_samples = 10000,
                         plot_type = "base", ...) {

  dots <- list(...)
  .check_legacy_level_arg(dots, "bfunnel()")

  BayesTools::check_bool(sampling_heterogeneity, "sampling_heterogeneity")
  BayesTools::check_bool(sampling_bias, "sampling_bias")
  BayesTools::check_char(
    plot_type,
    "plot_type",
    allow_values = c("base", "ggplot")
  )
  max_samples <- .normalize_funnel_max_samples(max_samples)
  dots        <- .set_dots_funnel(dots)

  if (.outcome_type(x) != "norm") {
    stop(
      "'bfunnel()' is available only for normal outcome models because a ",
      "scalar standard error does not index the fitted predictive ",
      "distribution of a discrete outcome model.",
      call. = FALSE
    )
  }
  if (.is_mods(x) || .is_scale(x)) {
    stop(
      "'bfunnel()' requires an intercept-only model with common marginal ",
      "heterogeneity. Use 'funnel()' for a LOO-PIT residual funnel.",
      call. = FALSE
    )
  }

  common_heterogeneity <- .funnel_common_heterogeneity(x)
  if (!common_heterogeneity[["common"]]) {
    .funnel_stop_no_common_heterogeneity("bfunnel()")
  }

  funnel_data <- .funnel_data_outcome(
    x                      = x,
    sampling_heterogeneity = sampling_heterogeneity,
    sampling_bias          = sampling_bias,
    max_samples            = max_samples,
    estimand               = "posterior_predictive",
    common_heterogeneity   = common_heterogeneity,
    dots                   = dots
  )

  if (isTRUE(dots[["as_data"]])) {
    return(funnel_data)
  }

  if (plot_type == "ggplot") {
    return(.funnel_plot_ggplot(funnel_data, dots))
  }

  .funnel_plot_base(funnel_data, dots)
  return(invisible(NULL))
}


.funnel_stop_no_common_heterogeneity <- function(caller) {

  stop(
    "'", caller, "' outcome contours require one common marginal ",
    "heterogeneity distribution across studies. This model has study-specific ",
    "marginal heterogeneity; use a LOO-PIT residual funnel with 'funnel()'.",
    call. = FALSE
  )
}
