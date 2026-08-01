# ============================================================================ #
# funnel-api.R
# ============================================================================ #

### funnel plot functions ----
#' @export
funnel <- function(x, ...) UseMethod("funnel")


#' @title Funnel Plot for brma Object
#'
#' @description \code{funnel.brma} creates a funnel plot for a fitted brma object.
#' For intercept-only models without scale regression, the default outcome mode
#' displays observed effect sizes against the fitted sampling distribution. For
#' models with location or scale moderators, the default residual mode displays
#' residuals against a standard-error funnel.
#'
#' @param x a fitted brma object
#' @param residual whether to use residual mode. Defaults to not specified,
#' which means the function automatically determines the mode:
#' \itemize{
#'   \item not specified: For intercept-only models without scale
#'     regression, displays observed effect sizes against the fitted sampling
#'     distribution funnel; for models with moderators or scale regression,
#'     automatically uses residual mode.
#'   \item \code{FALSE}: explicitly requests outcome mode.
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
#'     \code{\link{residuals.brma}} with \code{type = "outcome"}.
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
#' @param max_samples maximum number of posterior samples used for
#' model-averaged publication-bias funnel contours. Defaults to \code{10000}.
#' Use \code{Inf} to use all posterior samples.
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
#' The funnel plot has two modes. If \code{residual} is not specified, the mode
#' is chosen automatically from the fitted model: intercept-only models without
#' scale regression use outcome mode, whereas models with location or scale
#' moderators use residual mode.
#'
#' \strong{Outcome mode} (intercept-only models without scale regression):
#' Displays observed effect sizes on the x-axis and standard errors on the
#' y-axis. The reference line follows the center of the fitted sampling
#' distribution. When \code{sampling_bias = FALSE}, this center is the pooled
#' effect; when PET/PEESE bias adjustment is incorporated, the center line can
#' vary with the standard error. The funnel region represents the central 95\%
#' region of the sampling distribution, optionally incorporating heterogeneity
#' and publication bias.
#' For correlated known-\code{V} \code{brma.mv()} models, outcome-mode funnels
#' are descriptive scalar-SE displays based on the diagonal of \code{V}; residual
#' mode follows the fitted estimate-unit residual target.
#'
#'
#' \strong{Residual mode} (models with moderators or scale regression):
#' Displays residuals on the x-axis and standard errors on the
#' y-axis. The funnel region represents the central 95\% region of
#' \eqn{N(0, \mathrm{SE}^2)}. With \code{type = "LOO-PIT"}, the plotted
#' residuals and standard errors are the raw-scale LOO predictive companions
#' returned by \code{\link{rstudent.brma}}; the PIT-normalized \code{z} values
#' are used by \code{\link{qqnorm.brma}} and influence diagnostics. With
#' \code{type = "rstandard"}, the plotted values are internally standardized
#' residual companions from \code{\link{rstandard.brma}}. With
#' \code{type = "outcome"}, these are raw outcome residuals. Under a correctly
#' specified model, most points should fall within this region.
#'
#' The \code{type} argument controls how residuals are computed in residual
#' mode. See \code{\link{residuals.brma}} for details on each type.
#' The \code{sampling_heterogeneity} and \code{sampling_bias} arguments are
#' ignored in residual mode.
#'
#' For GLMM models, outcome-mode observed effect sizes are computed from the raw
#' frequency data using formulas equivalent to \code{metafor::escalc}.
#' LOO-PIT residual mode is unavailable because the fitted predictive
#' distribution is discrete and no discrete PIT convention has been defined.
#' Use \code{type = "outcome"} for a descriptive effect-size-scale residual
#' funnel when needed.
#'
#' @return If \code{as_data = TRUE}, \code{funnel.brma} returns a list with the
#' data used for plotting, including the plotted points, funnel polygons,
#' plotting limits, labels, and reference line. Otherwise, it returns
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
#'   [predict.brma()]
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
  is_mods  <- .is_mods(x)
  is_scale <- .is_scale(x)

  # determine mode: residual mode if mods/scale present, or if explicitly requested

  if (missing(residual)) {
    is_residual <- is_mods || is_scale
  } else {
    BayesTools::check_bool(residual, "residual")
    is_residual <- residual
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
