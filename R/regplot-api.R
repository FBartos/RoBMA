### regplot functions ----
#' @export
regplot <- function(x, ...) UseMethod("regplot")


#' @title Regression Plot (Bubble Plot) for brma Object
#'
#' @description \code{regplot.brma} creates a regression plot (also known as
#' bubble plot) for a fitted brma object with location or scale moderators. The
#' plot displays observed effect sizes against a moderator variable, with point
#' sizes proportional to study precision.
#'
#' @param x a fitted brma object with location or scale moderators
#' @param mod index or name of the moderator variable to plot on the x-axis.
#' If not specified and only one moderator is present, that moderator is used.
#' If multiple moderators are present, this argument is required.
#' @param pred logical; whether to show the prediction line. Defaults to \code{TRUE}.
#' @param ci logical; whether to show credible interval bands. Defaults to \code{TRUE}.
#' @param pi logical; whether to show prediction interval bands. Defaults to \code{FALSE}.
#' @param si logical; whether to show sampling interval bands. Defaults to \code{FALSE}.
#' The sampling interval shows the expected range of observed effect sizes,
#' incorporating both heterogeneity (tau) and a representative level of
#' sampling error (median SE across studies by default). When the model includes
#' publication bias adjustments and \code{sampling_bias = TRUE}, the
#' sampling interval incorporates the expected distortion from the
#' selection process.
#' @param level numeric; credible/prediction interval level in percent.
#' Defaults to \code{95}.
#' @param at numeric vector; for continuous moderators, values at which to
#' evaluate the prediction. If not specified, uses a sequence across the
#' observed range.
#' @param digits integer; number of decimal places for labels. Defaults to \code{2}.
#' @param transf function; transformation to apply to the y-axis (effect sizes).
#' Defaults to \code{NULL} (no transformation).
#' @param atransf reserved for axis-label transformations. Currently not
#' implemented and must be `NULL`.
#' @param targs reserved for additional transformation arguments. Currently not
#' implemented and must be `NULL`.
#' @param refline numeric; position of horizontal reference line.
#' Defaults to \code{NULL} (no reference line).
#' @param psize numeric vector or \code{NULL}; point sizes for each study.
#' If scalar, it is recycled to all studies. If \code{NULL} (default), sizes
#' are computed based on inverse sampling variance.
#' @param plim numeric vector of length 2; range for point sizes.
#' Defaults to \code{c(0.5, 3)}.
#' @param by character; name of a moderator variable to use for separate
#' lines/colors. Defaults to \code{NULL}. If omitted and the selected
#' moderator enters exactly one two-way interaction, the other variable in the
#' interaction is used automatically. Continuous \code{by} moderators are shown
#' at mean - SD, mean, and mean + SD.
#' @param legend logical; whether to show legend when \code{by} is specified.
#' Defaults to \code{TRUE}.
#' @param xlab character; x-axis label. Defaults to the moderator name.
#' @param ylab character; y-axis label. Defaults to "Observed Effect Size".
#' @param xlim numeric vector of length 2; x-axis limits. Defaults to data range.
#' @param ylim numeric vector of length 2; y-axis limits. Defaults to data range.
#' @param sampling_bias whether publication bias should be incorporated into
#' plotted predictions and sampling intervals. Defaults to \code{TRUE}. For
#' PET/PEESE models, this includes the expected regression bias in predictions.
#' For selection models, sampling-bias adjustment applies to sampling intervals
#' when \code{si = TRUE}; the mean prediction is unchanged.
#' @param sei single positive numeric value used as the reference standard
#' error for sampling-bias and sampling-interval calculations. Defaults to the
#' median observed standard error. When \code{mod = "sei"} or \code{mod =
#' "vi"}, the plotted moderator supplies the pointwise standard error and an
#' explicit \code{sei} override is not allowed.
#' @param max_samples maximum number of posterior samples used for prediction
#' summaries and interval bands. Defaults to \code{10000}. Use \code{Inf} to
#' use all posterior samples.
#' @param plot_type character; whether to use base R graphics (\code{"base"})
#' or ggplot2 (\code{"ggplot"}). Defaults to \code{"base"}.
#' @param as_data logical; if \code{TRUE}, returns plot data instead of
#' creating plot. Defaults to \code{FALSE}.
#' @param ... additional graphical arguments:
#' \describe{
#'   \item{main}{character string for plot title}
#'   \item{pch}{point symbol (default: 21)}
#'   \item{col}{point border color (default: "black")}
#'   \item{bg}{point fill color (default: "#A6A6A6")}
#'   \item{lcol}{line color (default: "black")}
#'   \item{lwd}{line width (default: 2)}
#'   \item{shade}{CI/PI/SI band shading (default: TRUE)}
#'   \item{col.ci}{CI band color (default: "gray70")}
#'   \item{col.pi}{PI band color (default: "gray85")}
#'   \item{col.si}{SI band color (default: "gray92")}
#'   \item{alpha.ci}{CI band transparency (default: 0.4)}
#'   \item{alpha.pi}{PI band transparency (default: 0.2)}
#'   \item{alpha.si}{SI band transparency (default: 0.15)}
#'   \item{jitter}{jitter amount for categorical moderators (default: 0.2)}
#'   \item{box.width}{box width for categorical interval summaries (default: 0.5)}
#' }
#'
#' @details
#' The regression plot (bubble plot) is a standard visualization for
#' meta-regression results. It displays:
#' \itemize{
#'   \item Observed effect sizes (y-axis) against moderator values (x-axis)
#'   \item Point sizes proportional to study precision (inverse variance)
#'   \item Prediction line showing the estimated regression relationship
#'   \item Confidence bands showing uncertainty in the mean prediction
#'   \item Optional prediction bands showing expected range of true effects
#'   \item Optional sampling interval bands showing expected range of observed outcomes
#' }
#'
#' For continuous moderators, predictions are computed across the observed
#' range of the moderator. For categorical moderators (factors), predictions
#' are computed at each factor level with optional jittering of points.
#' A moderator that appears only in the scale formula produces a constant mean
#' prediction while prediction and sampling intervals reflect its effect on
#' heterogeneity.
#'
#' The \code{by} argument allows displaying separate regression lines for
#' different levels of a second moderator, useful for visualizing interactions.
#' For \code{brma.mv()} known-\code{V} models, prediction and sampling bands are
#' pointwise descriptive bands. They summarize row-wise marginal random-effect
#' variation and, for sampling intervals, a representative sampling standard
#' error. They are not full joint predictive regions for correlated
#' known-\code{V} outcomes.
#'
#' @return \code{regplot.brma} returns \code{NULL} invisibly if
#' \code{plot_type = "base"} or a ggplot object if \code{plot_type = "ggplot"}.
#' If \code{as_data = TRUE}, returns a list with plot data components.
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
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     mods    = ~ ablat + year,
#'     data    = dat,
#'     measure = "RR"
#'   )
#'   regplot(fit, mod = "ablat")
#'   regplot(fit, mod = "year", pi = TRUE, si = TRUE)
#'   regplot(fit, mod = "ablat", plot_type = "ggplot")
#'
#'   fit_cat <- brma(yi = yi, vi = vi, mods = ~ alloc, data = dat, measure = "RR")
#'   regplot(fit_cat)
#' }
#' }
#'
#' @seealso [funnel.brma()], [predict.brma()]
#' @aliases regplot
#' @export
#' @exportS3Method metafor::regplot
#' @rdname regplot
regplot.brma <- function(x, mod = NULL, pred = TRUE, ci = TRUE, pi = FALSE, si = FALSE,
                         level = 95, at = NULL, digits = 2,
                         transf = NULL, atransf = NULL, targs = NULL,
                         refline = NULL, psize = NULL, plim = c(0.5, 3),
                         by = NULL, legend = TRUE,
                         xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                         sampling_bias = TRUE, sei = NULL, max_samples = 10000,
                         plot_type = "base", as_data = FALSE, ...) {

  # input validation
  BayesTools::check_bool(pred, "pred")
  BayesTools::check_bool(ci, "ci")
  BayesTools::check_bool(pi, "pi")
  BayesTools::check_bool(si, "si")
  BayesTools::check_real(level, "level", lower = 50, upper = 99.9)
  .check_plot_numeric(at, "at", allow_null = TRUE)
  BayesTools::check_int(digits, "digits", lower = 0)
  .check_plot_function(transf, "transf")
  .check_plot_numeric(refline, "refline", check_length = 1, allow_null = TRUE)
  .check_plot_numeric(plim, "plim", check_length = 2, lower = 0, allow_null = FALSE)
  if (plim[1] >= plim[2]) {
    stop("'plim' must be an increasing numeric vector of length 2.", call. = FALSE)
  }
  BayesTools::check_bool(legend, "legend")
  .check_plot_label(xlab, "xlab")
  .check_plot_label(ylab, "ylab")
  .check_plot_limits(xlim, "xlim")
  .check_plot_limits(ylim, "ylim")
  BayesTools::check_bool(sampling_bias, "sampling_bias")
  max_samples <- .normalize_max_samples(max_samples, "max_samples")
  if (!is.null(sei)) {
    BayesTools::check_real(sei, "sei", lower = 0)
    if (sei <= 0) {
      stop("'sei' must be positive.", call. = FALSE)
    }
  }
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(as_data, "as_data")

  # not yet implemented arguments
  if (!is.null(atransf))
    stop("The 'atransf' argument is not yet implemented.", call. = FALSE)
  if (!is.null(targs))
    stop("The 'targs' argument is not yet implemented.", call. = FALSE)

  # check that model has location or scale moderators
  predictor_info <- .regplot_predictor_info(x)
  if (length(predictor_info) == 0L) {
    stop("regplot requires a model with moderators. ",
         "Use funnel() for intercept-only models.", call. = FALSE)
  }

  # set up graphical arguments with defaults
  dots <- .set_dots_regplot(...)

  # identify the moderator to plot
  mod_info <- .regplot_get_moderator(mod, predictor_info)
  mod_name <- mod_info$name
  mod_type <- mod_info$type
  mod_data <- mod_info$data
  if (mod_name %in% c("sei", "vi") && !is.null(sei)) {
    stop(
      "'sei' cannot be supplied when the plotted moderator is 'sei' or 'vi'.",
      call. = FALSE
    )
  }
  if (mod_type == "categorical" && !is.null(at)) {
    stop("'at' is only available for continuous moderators.", call. = FALSE)
  }

  K <- length(.outcome_data_yi(x))
  if (!is.null(psize)) {
    .check_plot_numeric(psize, "psize", lower = 0, allow_null = FALSE)
    if (!length(psize) %in% c(1L, K)) {
      stop("'psize' must be either a scalar or have one value per study.", call. = FALSE)
    }
  }

  # identify grouping variable (explicit or auto-detected interaction)
  by_info <- .regplot_get_by(x, by, mod_name, predictor_info)

  # generate plot data
  regplot_data <- .regplot_data(
    x              = x,
    mod_name       = mod_name,
    mod_type       = mod_type,
    mod_data       = mod_data,
    by_info        = by_info,
    predictor_info = predictor_info,
    pred           = pred,
    ci             = ci,
    pi             = pi,
    si             = si,
    level          = level,
    at             = at,
    digits         = digits,
    psize          = psize,
    plim           = plim,
    transf         = transf,
    xlim           = xlim,
    ylim           = ylim,
    xlab           = xlab,
    ylab           = ylab,
    refline        = refline,
    sampling_bias  = sampling_bias,
    max_samples    = max_samples,
    reference_sei  = sei,
    dots           = dots
  )

  # allow data return for programmatic access
  if (isTRUE(as_data)) {
    return(regplot_data)
  }

  # create plot
  if (plot_type == "ggplot") {
    return(.regplot_plot_ggplot(regplot_data, dots, legend = legend))
  } else {
    .regplot_plot_base(regplot_data, dots, legend = legend)
    return(invisible(NULL))
  }
}
