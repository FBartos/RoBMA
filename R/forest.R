# ============================================================================ #
# brma.forest.R
# ============================================================================ #
#
# Forest plot adapters for brma objects using metafor's public plotting API.
#
# ============================================================================ #


#' @title Prepare brma Data for metafor Forest Plots
#'
#' @description \code{as_metafor_forest} extracts observed effect sizes,
#' standard errors, study labels, and RoBMA posterior summary intervals in a
#' structure that can be passed to \pkg{metafor}'s forest plotting functions.
#'
#' @param x a fitted brma object.
#' @param level credible/confidence level in percent. Defaults to \code{95}.
#' @param addpred logical; whether to compute a posterior prediction interval
#'   for the pooled effect. Defaults to \code{FALSE}.
#' @param bias_adjusted logical; whether pooled and prediction summaries should
#'   adjust for publication bias. Defaults to \code{TRUE}.
#' @param conditional logical; whether to return conditional posterior
#'   summaries for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @param ... additional named arguments stored in \code{forest_args} and later
#'   passed to \code{\link[metafor]{forest.default}}. Shared graphical
#'   arguments are also stored in \code{addpoly_args}. This includes
#'   \pkg{metafor}'s auxiliary-column arguments \code{ilab},
#'   \code{ilab.lab}, \code{ilab.xpos}, and \code{ilab.pos}.
#'
#' @return A list of class \code{metafor_forest.brma} with study-level data,
#' pooled posterior summaries, and ready-to-use argument lists for
#' \code{\link[metafor]{forest.default}} and
#' \code{\link[metafor]{addpoly.default}}. The \code{forest_args} include
#' enough vertical space for the default summary row. The returned object can
#' also be plotted directly with \code{\link[metafor]{forest}}.
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
#'   count_columns <- as.matrix(dat[, c("tpos", "tneg", "cpos", "cneg")])
#'   forest_data <- as_metafor_forest(
#'     fit,
#'     addpred   = TRUE,
#'     slab      = paste(dat$author, dat$year),
#'     ilab      = count_columns,
#'     ilab.lab  = c("T+", "T-", "C+", "C-"),
#'     ilab.xpos = c(-7, -6, -5, -4),
#'     xlim      = c(-10, 4),
#'     atransf   = exp
#'   )
#'   metafor::forest(forest_data, predstyle = "bar")
#'
#'   # Or draw the metafor calls directly:
#'   do.call(metafor::forest.default, forest_data$forest_args)
#'   do.call(metafor::addpoly.default, forest_data$addpoly_args)
#' }
#' }
#'
#' @export
as_metafor_forest <- function(x, ...) {
  UseMethod("as_metafor_forest")
}


#' @rdname as_metafor_forest
#' @export
as_metafor_forest.brma <- function(x, level = 95, addpred = FALSE,
                                   bias_adjusted = TRUE,
                                   conditional = FALSE, ...) {

  level <- .forest_normalize_level(level)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  dots  <- list(...)

  BayesTools::check_bool(addpred, "addpred")
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(conditional, "conditional")
  .forest_check_reserved_dots(dots)
  if (conditional && !.is_RoBMA(x)) {
    stop("'conditional' summaries are available only for RoBMA objects.",
         call. = FALSE)
  }

  yi       <- .outcome_data_yi(x)
  sei      <- .outcome_data_sei(x)
  slab     <- .get_estimate_labels(x)
  measure  <- .measure(x)
  z_value  <- stats::qnorm(1 - (1 - level) / 2)
  yi_plot  <- yi
  attr(yi_plot, "measure") <- measure
  attr(yi_plot, "slab")    <- slab

  studies <- data.frame(
    yi       = yi,
    sei      = sei,
    vi       = sei^2,
    ci.lb    = yi - z_value * sei,
    ci.ub    = yi + z_value * sei,
    slab     = slab,
    stringsAsFactors = FALSE
  )
  if (.is_weights(x)) {
    studies[["likelihood_weight"]] <- .outcome_data_weights(x)
  }

  pooled <- .forest_summary_row(
    pooled_effect(
      object        = x,
      probs         = probs,
      bias_adjusted = bias_adjusted,
      conditional   = conditional
    )
  )

  pooled[["pi.lb"]] <- NA_real_
  pooled[["pi.ub"]] <- NA_real_

  if (addpred) {
    prediction <- .forest_summary_row(
      predict.brma(
        object        = x,
        newdata       = TRUE,
        type          = "estimate",
        probs         = probs,
        bias_adjusted = bias_adjusted,
        conditional   = conditional,
        quiet         = TRUE
      )
    )
    pooled[["pi.lb"]] <- prediction[["ci.lb"]]
    pooled[["pi.ub"]] <- prediction[["ci.ub"]]
  }

  forest_args <- list(
    x     = yi_plot,
    sei   = studies[["sei"]],
    slab  = studies[["slab"]],
    level = 100 * level
  )
  forest_args       <- .forest_merge_args(forest_args, dots)
  default_predstyle <- if (addpred) "bar" else "line"
  default_ylim      <- is.null(forest_args[["ylim"]])
  if (default_ylim) {
    forest_args[["ylim"]] <- .forest_default_ylim(
      k         = length(yi),
      row       = -1,
      dots      = forest_args,
      predstyle = default_predstyle
    )
  }

  addpoly_args <- list(
    x         = pooled[["estimate"]],
    ci.lb     = pooled[["ci.lb"]],
    ci.ub     = pooled[["ci.ub"]],
    rows      = -1,
    level     = 100 * level,
    predstyle = default_predstyle,
    mlab      = "Pooled Effect"
  )
  if (addpred) {
    prediction_se <- .forest_prediction_se(
      pi.lb = pooled[["pi.lb"]],
      pi.ub = pooled[["pi.ub"]],
      level = level
    )
    pi_lb <- pooled[["pi.lb"]]
    attr(pi_lb, "level") <- 1 - level
    attr(pi_lb, "dist")  <- "norm"
    attr(pi_lb, "se")    <- prediction_se

    addpoly_args[["pi.lb"]] <- pi_lb
    addpoly_args[["pi.ub"]] <- pooled[["pi.ub"]]
  }
  addpoly_args <- .forest_addpoly_shared_args(addpoly_args, forest_args)

  out <- list(
    studies       = studies,
    pooled        = pooled,
    forest_args   = forest_args,
    addpoly_args  = addpoly_args,
    measure       = measure,
    level         = 100 * level,
    addpred       = addpred,
    bias_adjusted = bias_adjusted,
    conditional   = conditional,
    default_ylim  = default_ylim
  )
  class(out) <- c("metafor_forest.brma", "list")

  return(out)
}


#' @title Forest Plot for brma Objects Using metafor
#'
#' @description \code{metafor::forest()} creates a forest plot for a fitted
#' brma object by adapting RoBMA results to \pkg{metafor}'s public forest
#' plotting API. Study rows are drawn by
#' \code{\link[metafor]{forest.default}} and the RoBMA pooled posterior
#' summary is added with
#' \code{\link[metafor]{addpoly.default}}.
#'
#' @param x a fitted brma object or an object returned by
#'   \code{as_metafor_forest()}.
#' @param addfit logical; whether to add the RoBMA pooled-effect summary
#'   polygon. Defaults to \code{TRUE}.
#' @param addpred logical; whether to add a posterior prediction interval to
#'   the pooled-effect summary polygon. Defaults to \code{FALSE}.
#' @param predstyle character; prediction interval style passed to
#'   \code{\link[metafor]{addpoly.default}}. Defaults to \code{"line"}.
#' @param mlab label for the pooled-effect summary row. Defaults to
#'   \code{"Pooled Effect"}.
#' @param row numeric row position for the pooled-effect summary polygon.
#'   Defaults to \code{-1}.
#' @param level credible/confidence level in percent. Defaults to \code{95}.
#' @param bias_adjusted logical; whether pooled and prediction summaries should
#'   adjust for publication bias. Defaults to \code{TRUE}.
#' @param conditional logical; whether to return conditional posterior
#'   summaries for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @param as_data logical; if \code{TRUE}, returns the prepared data instead of
#'   creating a plot.
#' @param predlim argument passed to
#'   \code{\link[metafor]{addpoly.default}} for prediction interval rendering.
#' @param border,constarea arguments passed to
#'   \code{\link[metafor]{addpoly.default}} for summary polygon rendering.
#' @param ... additional arguments passed to
#'   \code{\link[metafor]{forest.default}}. This includes auxiliary-column
#'   arguments such as \code{ilab}, \code{ilab.lab}, \code{ilab.xpos}, and
#'   \code{ilab.pos}. Shared graphical arguments such as \code{annotate},
#'   \code{digits}, \code{transf}, \code{atransf}, \code{efac}, \code{col},
#'   \code{lty}, \code{fonts}, and \code{cex} are also forwarded to
#'   \code{\link[metafor]{addpoly.default}}.
#'
#' @details
#' This method deliberately uses metafor's public vector-based forest plotting
#' API instead of constructing a fake metafor model object. The observed study
#' rows are shown as usual Wald intervals based on the supplied effect sizes
#' and standard errors. The model summary row uses the posterior mean and
#' credible interval from \code{\link{pooled_effect}}. When
#' \code{addpred = TRUE}, the prediction interval is computed from
#' \code{\link{predict.brma}} with \code{type = "estimate"} and
#' \code{newdata = TRUE}. For \pkg{metafor}'s \code{predstyle = "shade"} and
#' \code{predstyle = "dist"}, the displayed predictive distribution uses
#' metafor's normal approximation implied by the posterior prediction interval
#' endpoints.
#'
#' Additional study-level columns can be added with metafor's \code{ilab}
#' interface. When supplied to \code{as_metafor_forest()} or to
#' \code{metafor::forest(..., as_data = TRUE)}, those arguments are stored in
#' the adapter object's \code{forest_args}; when supplied to
#' \code{metafor::forest()} at plot time, they are forwarded directly.
#'
#' Use \code{as_metafor_forest()} to prepare the adapter object without
#' requiring a \pkg{metafor} installation.
#'
#' @return Invisibly returns the object returned by
#' \code{\link[metafor]{forest.default}}. If \code{as_data = TRUE}, returns a
#' \code{metafor_forest.brma} object.
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
#'   metafor::forest(fit)
#'   metafor::forest(fit, addpred = TRUE, atransf = exp)
#'
#'   count_columns <- as.matrix(dat[, c("tpos", "tneg", "cpos", "cneg")])
#'   metafor::forest(
#'     fit,
#'     addpred   = TRUE,
#'     slab      = paste(dat$author, dat$year),
#'     ilab      = count_columns,
#'     ilab.lab  = c("T+", "T-", "C+", "C-"),
#'     ilab.xpos = c(-8.5, -7.5, -6.5, -5.5),
#'     xlim      = c(-10, 4),
#'     atransf   = exp
#'   )
#' }
#' }
#'
#' @seealso [as_metafor_forest()], [pooled_effect()], [predict.brma()]
#' @exportS3Method metafor::forest
#' @rdname forest
forest.brma <- function(x, addfit = TRUE, addpred = FALSE,
                        predstyle = "line", mlab = NULL, row = -1,
                        level = 95, bias_adjusted = TRUE,
                        conditional = FALSE, as_data = FALSE,
                        predlim, border, constarea, ...) {

  BayesTools::check_bool(addfit, "addfit")
  BayesTools::check_bool(addpred, "addpred")
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(as_data, "as_data")
  BayesTools::check_char(predstyle, "predstyle",
                         allow_values = c("line", "bar", "shade", "dist"))
  BayesTools::check_real(row, "row", check_length = 1, allow_NA = FALSE)

  if (addpred && !addfit) {
    stop("'addpred = TRUE' requires 'addfit = TRUE'.", call. = FALSE)
  }
  if (predstyle != "line" && !addpred) {
    stop("'predstyle' values other than 'line' require 'addpred = TRUE'.",
         call. = FALSE)
  }

  if (as_data) {
    return(as_metafor_forest.brma(
      x             = x,
      level         = level,
      addpred       = addpred,
      bias_adjusted = bias_adjusted,
      conditional   = conditional,
      ...
    ))
  }

  forest_data <- as_metafor_forest.brma(
    x             = x,
    level         = level,
    addpred       = addpred,
    bias_adjusted = bias_adjusted,
    conditional   = conditional
  )

  return(.forest_plot_metafor(
    forest_data = forest_data,
    addfit      = addfit,
    dots        = list(...),
    mlab        = mlab,
    row         = row,
    predstyle   = predstyle,
    predlim     = if (missing(predlim)) NULL else predlim,
    border      = if (missing(border)) NULL else border,
    constarea   = if (missing(constarea)) NULL else constarea
  ))
}


#' @rdname forest
#' @exportS3Method metafor::forest
forest.metafor_forest.brma <- function(x, addfit = TRUE,
                                        predstyle = "line", mlab = NULL,
                                        row = -1, predlim, border,
                                        constarea, ...) {

  return(.forest_plot_metafor(
    forest_data = x,
    addfit      = addfit,
    dots        = list(...),
    mlab        = mlab,
    row         = row,
    predstyle   = predstyle,
    predlim     = if (missing(predlim)) NULL else predlim,
    border      = if (missing(border)) NULL else border,
    constarea   = if (missing(constarea)) NULL else constarea
  ))
}


.forest_require_metafor <- function() {

  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop(
      "Package 'metafor' is required for forest.brma(). ",
      "Install it with install.packages(\"metafor\").",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


.forest_normalize_level <- function(level) {

  BayesTools::check_real(level, "level", check_length = 1, allow_NA = FALSE)
  if (level > 1) {
    level <- level / 100
  }
  if (level <= 0 || level >= 1) {
    stop("'level' must be greater than 0 and less than 100.", call. = FALSE)
  }

  return(level)
}


.forest_summary_row <- function(samples) {

  summary_table <- as.data.frame(summary(samples))
  if (nrow(summary_table) != 1L) {
    stop("Expected a single-row posterior summary.", call. = FALSE)
  }
  if (ncol(summary_table) < 4L ||
      !"Mean" %in% colnames(summary_table)) {
    stop("Posterior summary does not contain the expected columns.",
         call. = FALSE)
  }

  out <- data.frame(
    estimate = summary_table[1L, "Mean"],
    ci.lb    = summary_table[1L, 3L],
    ci.ub    = summary_table[1L, 4L],
    stringsAsFactors = FALSE
  )

  return(out)
}


.forest_prediction_se <- function(pi.lb, pi.ub, level) {

  crit <- stats::qnorm(1 - (1 - level) / 2)
  se   <- (pi.ub - pi.lb) / (2 * crit)

  if (!is.finite(se) || se <= 0) {
    return(NA_real_)
  }

  return(se)
}


.forest_plot_metafor <- function(forest_data, addfit, dots, mlab, row,
                                 predstyle, predlim, border, constarea) {

  BayesTools::check_bool(addfit, "addfit")
  BayesTools::check_char(predstyle, "predstyle",
                         allow_values = c("line", "bar", "shade", "dist"))
  BayesTools::check_real(row, "row", check_length = 1, allow_NA = FALSE)

  if (addfit && predstyle != "line" && !isTRUE(forest_data[["addpred"]])) {
    stop("'predstyle' values other than 'line' require a forest adapter ",
         "created with 'addpred = TRUE'.", call. = FALSE)
  }

  .forest_require_metafor()
  .forest_check_reserved_dots(dots)

  custom_ylim <- "ylim" %in% names(dots)
  forest_args <- .forest_merge_args(forest_data[["forest_args"]], dots)

  if ((is.null(forest_args[["ylim"]]) ||
       (isTRUE(forest_data[["default_ylim"]]) && !custom_ylim)) &&
      addfit) {
    forest_args[["ylim"]] <- .forest_default_ylim(
      k         = nrow(forest_data[["studies"]]),
      row       = row,
      dots      = forest_args,
      predstyle = predstyle
    )
  }

  plot <- do.call(metafor::forest.default, forest_args)

  if (addfit) {
    addpoly_args <- .forest_addpoly_args(
      forest_data = forest_data,
      dots        = forest_args,
      mlab        = mlab,
      row         = row,
      predstyle   = predstyle,
      predlim     = predlim,
      border      = border,
      constarea   = constarea
    )
    do.call(metafor::addpoly.default, addpoly_args)
  }

  return(invisible(plot))
}


.forest_merge_args <- function(defaults, dots) {

  if (length(dots) == 0L) {
    return(defaults)
  }

  replace_names <- intersect(names(defaults), names(dots))
  for (name in replace_names) {
    defaults[[name]] <- dots[[name]]
  }

  return(c(defaults, dots[setdiff(names(dots), names(defaults))]))
}


.forest_check_reserved_dots <- function(dots) {

  if (length(dots) > 0L &&
      (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All arguments passed through '...' must be named.", call. = FALSE)
  }

  reserved <- intersect(
    names(dots),
    c(
      "x", "vi", "sei", "ci.lb", "ci.ub", "level",
      "addfit", "addpred", "predstyle", "mlab", "row", "predlim",
      "border", "constarea", "bias_adjusted", "conditional", "as_data"
    )
  )
  if (length(reserved) > 0L) {
    stop(
      "The following arguments are controlled by the RoBMA forest adapter ",
      "and cannot be passed through '...': ",
      paste(reserved, collapse = ", "), ".",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


.forest_default_ylim <- function(k, row, dots, predstyle) {

  top <- dots[["top"]]
  if (is.null(top)) {
    top <- 3
  }
  BayesTools::check_real(top, "top", check_length = 1, allow_NA = FALSE)

  study_rows <- dots[["rows"]]
  if (is.null(study_rows)) {
    study_rows <- k:1
  } else if (length(study_rows) == 1L) {
    study_rows <- study_rows:(study_rows - k + 1L)
  }

  lower <- min(row - ifelse(predstyle == "line", 1, 2), 0)
  upper <- max(study_rows, na.rm = TRUE) + top

  return(c(lower, upper))
}


.forest_addpoly_args <- function(forest_data, dots, mlab, row, predstyle,
                                 predlim, border, constarea) {

  addpoly_args <- forest_data[["addpoly_args"]]
  addpoly_args[["rows"]]      <- row
  addpoly_args[["predstyle"]] <- predstyle
  addpoly_args[["mlab"]]      <- if (is.null(mlab)) "Pooled Effect" else mlab

  if (!is.null(predlim)) {
    addpoly_args[["predlim"]] <- predlim
  }
  if (!is.null(border)) {
    addpoly_args[["border"]] <- border
  }
  if (!is.null(constarea)) {
    addpoly_args[["constarea"]] <- constarea
  }

  addpoly_args <- .forest_addpoly_shared_args(addpoly_args, dots)

  return(addpoly_args)
}


.forest_addpoly_shared_args <- function(addpoly_args, dots) {

  shared <- c(
    "annotate", "digits", "width", "transf", "atransf", "targs",
    "efac", "col", "lty", "fonts", "cex", "lwd"
  )
  for (name in intersect(shared, names(dots))) {
    addpoly_args[[name]] <- dots[[name]]
  }

  return(addpoly_args)
}
