# ============================================================================ #
# brma.regplot.R
# ============================================================================ #
#
# This file contains the regression plot (bubble plot) functions for brma
# objects with moderators.
#
# Bubble plots display observed effect sizes against a moderator variable,
# with point sizes proportional to study precision (inverse variance).
# For continuous moderators, a regression line with CI/PI bands is shown.
# For categorical moderators, grouped points with optional jitter are shown.
#
# ============================================================================ #


### regplot functions ----
#' @export
regplot <- function(x, ...) UseMethod("regplot")


#' @title Regression Plot (Bubble Plot) for brma Object
#'
#' @description \code{regplot.brma} creates a regression plot (also known as
#' bubble plot) for a fitted brma object with moderators. The plot displays
#' observed effect sizes against a moderator variable, with point sizes
#' proportional to study precision.
#'
#' @param x a fitted brma object with moderators
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
#' median observed standard error.
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
#'
#' The \code{by} argument allows displaying separate regression lines for
#' different levels of a second moderator, useful for visualizing interactions.
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
                         sampling_bias = TRUE, sei = NULL,
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

  # check that model has moderators
  if (!.is_mods(x)) {
    stop("regplot requires a model with moderators. ",
         "Use funnel() for intercept-only models.", call. = FALSE)
  }

  # set up graphical arguments with defaults
  dots <- .set_dots_regplot(...)

  # identify the moderator to plot
  mod_info <- .regplot_get_moderator(x, mod)
  mod_name <- mod_info$name
  mod_type <- mod_info$type
  mod_data <- mod_info$data
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
  by_info <- .regplot_get_by(x, by, mod_name)

  # generate plot data
  regplot_data <- .regplot_data(
    x             = x,
    mod_name      = mod_name,
    mod_type      = mod_type,
    mod_data      = mod_data,
    by_info       = by_info,
    pred          = pred,
    ci            = ci,
    pi            = pi,
    si            = si,
    level         = level,
    at            = at,
    digits        = digits,
    psize         = psize,
    plim          = plim,
    transf        = transf,
    xlim          = xlim,
    ylim          = ylim,
    xlab          = xlab,
    ylab          = ylab,
    refline       = refline,
    sampling_bias = sampling_bias,
    reference_sei = sei,
    dots          = dots
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


# ---------------------------------------------------------------------------- #
# .regplot_get_moderator
# ---------------------------------------------------------------------------- #
#
# Identify and validate the moderator variable for plotting.
#
# @param x     brma object
# @param mod   index or name of moderator (NULL for auto-detect)
#
# @return list with name, type, and data for the moderator
#
# ---------------------------------------------------------------------------- #
.regplot_get_moderator <- function(x, mod) {

  mods_data <- x[["data"]][["mods"]]
  design    <- .fitted_formula_design(x, "mu", required = TRUE)

  if (is.null(mods_data)) {
    stop("No moderator data found in the model.", call. = FALSE)
  }

  mod_names <- design[["predictors"]]

  if (is.null(mod)) {
    # auto-detect: use single moderator or error if multiple
    if (length(mod_names) == 1) {
      mod_name <- mod_names[1]
    } else {
      stop("Multiple moderators present. Please specify 'mod' argument. ",
           "Available moderators: ", paste(mod_names, collapse = ", "),
           call. = FALSE)
    }
  } else if (is.numeric(mod)) {
    # index-based selection
    if (mod < 1 || mod > length(mod_names)) {
      stop("Invalid moderator index. Must be between 1 and ",
           length(mod_names), ".", call. = FALSE)
    }
    mod_name <- mod_names[mod]
  } else if (is.character(mod)) {
    # name-based selection
    if (!mod %in% mod_names) {
      stop("Moderator '", mod, "' not found. ",
           "Available moderators: ", paste(mod_names, collapse = ", "),
           call. = FALSE)
    }
    mod_name <- mod
  } else {
    stop("'mod' must be NULL, a character name, or a numeric index.",
         call. = FALSE)
  }

  # extract moderator data and determine type
  mod_values <- mods_data[[mod_name]]
  predictor_type <- design[["predictor_types"]][[mod_name]]

  if (identical(predictor_type, "factor") ||
      is.factor(mod_values) || is.character(mod_values)) {
    mod_type <- "categorical"
    mod_values <- factor(
      mod_values,
      levels = if (!is.null(design[["xlevels"]][[mod_name]])) {
        design[["xlevels"]][[mod_name]]
      } else {
        levels(factor(mod_values))
      }
    )
  } else {
    mod_type <- "continuous"
  }

  return(list(
    name = mod_name,
    type = mod_type,
    data = mod_values
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_get_by
# ---------------------------------------------------------------------------- #
#
# Identify and validate a grouping moderator for interaction displays.
#
# ---------------------------------------------------------------------------- #
.regplot_get_by <- function(x, by, mod_name) {

  mods_data <- x[["data"]][["mods"]]
  design    <- .fitted_formula_design(x, "mu", required = TRUE)

  if (is.null(by)) {
    by <- .regplot_detect_interaction_by(x, mod_name)
  }

  if (is.null(by)) {
    return(NULL)
  }

  mod_names <- design[["predictors"]]

  if (!is.character(by) || length(by) != 1L) {
    stop("'by' must be NULL or a single moderator name.", call. = FALSE)
  }
  if (by == mod_name) {
    stop("'by' must be different from 'mod'.", call. = FALSE)
  }
  if (!by %in% mod_names) {
    stop("Moderator '", by, "' not found. ",
         "Available moderators: ", paste(mod_names, collapse = ", "),
         call. = FALSE)
  }

  by_values <- mods_data[[by]]
  predictor_type <- design[["predictor_types"]][[by]]

  if (identical(predictor_type, "factor") ||
      is.factor(by_values) || is.character(by_values)) {
    by_type <- "categorical"
    by_values <- factor(
      by_values,
      levels = if (!is.null(design[["xlevels"]][[by]])) {
        design[["xlevels"]][[by]]
      } else {
        levels(factor(by_values))
      }
    )
  } else {
    by_type <- "continuous"
  }

  return(list(
    name = by,
    type = by_type,
    data = by_values
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_detect_interaction_by
# ---------------------------------------------------------------------------- #
#
# Auto-select the second variable from a unique two-way interaction containing
# the plotted moderator.
#
# ---------------------------------------------------------------------------- #
.regplot_detect_interaction_by <- function(x, mod_name) {

  mods_data <- x[["data"]][["mods"]]

  if (is.null(mods_data)) {
    return(NULL)
  }

  term_labels <- .fitted_formula_terms(
    object            = x,
    parameter         = "mu",
    include_intercept = FALSE,
    display           = TRUE
  )
  term_labels <- term_labels[grepl(":", term_labels, fixed = TRUE)]

  if (length(term_labels) == 0L) {
    return(NULL)
  }

  candidates <- character()

  for (term in term_labels) {
    variables <- strsplit(term, ":", fixed = TRUE)[[1]]
    variables <- trimws(variables)

    if (length(variables) == 2L && mod_name %in% variables) {
      candidates <- c(candidates, setdiff(variables, mod_name))
    }
  }

  candidates <- unique(candidates[candidates %in% names(mods_data)])

  if (length(candidates) == 0L) {
    return(NULL)
  }
  if (length(candidates) > 1L) {
    stop("The selected moderator enters multiple interactions. ",
         "Please specify 'by'. Candidate moderators: ",
         paste(candidates, collapse = ", "), ".", call. = FALSE)
  }

  return(candidates)
}


# ---------------------------------------------------------------------------- #
# .regplot_band_data_continuous
# ---------------------------------------------------------------------------- #
#
# Build continuous interval-band data with one row per polygon vertex while
# preserving the underlying interval bounds for programmatic access.
#
# @param xpred numeric vector of prediction x values
# @param lower numeric vector of lower interval bounds
# @param upper numeric vector of upper interval bounds
#
# @return data.frame with consistent polygon and interval metadata
#
# ---------------------------------------------------------------------------- #
.regplot_band_data_continuous <- function(xpred, lower, upper, group = "All",
                                          group_id = 1L) {

  band_x     <- c(xpred, rev(xpred))
  band_lower <- c(lower, rev(lower))
  band_upper <- c(upper, rev(upper))

  data.frame(
    x     = unname(band_x),
    y     = unname(c(lower, rev(upper))),
    lower = unname(band_lower),
    upper = unname(band_upper),
    xpred = unname(band_x),
    group = rep(group, length(band_x)),
    group_id = rep(group_id, length(band_x))
  )
}


# ---------------------------------------------------------------------------- #
# .regplot_prediction_grid
# ---------------------------------------------------------------------------- #
#
# Build moderator values for predictions. The selected moderator forms the
# x-axis; an optional `by` moderator expands the grid into interaction lines.
#
# ---------------------------------------------------------------------------- #
.regplot_prediction_grid <- function(x, mod_name, mod_type, mod_data,
                                     by_info, at, digits) {

  mods_data <- x[["data"]][["mods"]]
  design    <- .fitted_formula_design(x, "mu", required = TRUE)

  if (mod_type == "continuous") {
    if (is.null(at)) {
      x_range <- range(mod_data)
      at_pred <- seq(x_range[1], x_range[2], length.out = 101)
    } else {
      at_pred <- at
    }
    x_values  <- at_pred
    x_numeric <- at_pred
    x_levels  <- NULL
  } else {
    x_levels  <- levels(mod_data)
    at_pred   <- x_levels
    x_values  <- factor(at_pred, levels = x_levels)
    x_numeric <- seq_along(at_pred)
  }

  if (is.null(by_info)) {
    by_values <- list(
      values = NULL,
      labels = "All",
      levels = NULL
    )
  } else {
    by_values <- .regplot_by_values(by_info, digits)
  }

  grid <- expand.grid(
    x_id     = seq_along(at_pred),
    group_id = seq_along(by_values[["labels"]]),
    KEEP.OUT.ATTRS = FALSE
  )
  grid[["group"]] <- by_values[["labels"]][grid[["group_id"]]]
  grid[["x"]]     <- x_numeric[grid[["x_id"]]]
  if (mod_type == "categorical") {
    grid[["level"]] <- at_pred[grid[["x_id"]]]
  }

  n_pred       <- nrow(grid)
  newdata_list <- list()

  for (nm in design[["predictors"]]) {
    if (nm == mod_name) {
      if (mod_type == "continuous") {
        newdata_list[[nm]] <- x_values[grid[["x_id"]]]
      } else {
        newdata_list[[nm]] <- factor(at_pred[grid[["x_id"]]], levels = x_levels)
      }
    } else if (!is.null(by_info) && nm == by_info[["name"]]) {
      if (by_info[["type"]] == "continuous") {
        newdata_list[[nm]] <- by_values[["values"]][grid[["group_id"]]]
      } else {
        newdata_list[[nm]] <- factor(
          by_values[["values"]][grid[["group_id"]]],
          levels = by_values[["levels"]]
        )
      }
    } else {
      other_vals <- mods_data[[nm]]
      if (is.factor(other_vals) || is.character(other_vals)) {
        other_levels <- if (!is.null(design[["xlevels"]][[nm]])) {
          design[["xlevels"]][[nm]]
        } else {
          levels(factor(other_vals))
        }
        newdata_list[[nm]] <- factor(
          rep(other_levels[1], n_pred),
          levels = other_levels
        )
      } else {
        newdata_list[[nm]] <- rep(mean(other_vals), n_pred)
      }
    }
  }

  return(list(
    newdata  = as.data.frame(newdata_list),
    grid     = grid,
    at_pred  = at_pred,
    groups   = by_values[["labels"]],
    x_levels = x_levels
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_by_values
# ---------------------------------------------------------------------------- #
#
# Values used for the grouping moderator.
#
# ---------------------------------------------------------------------------- #
.regplot_by_values <- function(by_info, digits) {

  by_data <- by_info[["data"]]

  if (by_info[["type"]] == "continuous") {
    by_mean <- mean(by_data)
    by_sd   <- stats::sd(by_data)
    values  <- c(by_mean - by_sd, by_mean, by_mean + by_sd)
    labels  <- c("Mean - 1 SD", "Mean", "Mean + 1 SD")

    if (!is.finite(by_sd) || by_sd == 0) {
      values <- by_mean
      labels <- format(round(by_mean, digits = digits), trim = TRUE)
    }

    labels <- paste0(by_info[["name"]], ": ", labels)
    return(list(values = values, labels = labels, levels = NULL))
  }

  by_factor <- factor(by_data)

  return(list(
    values = levels(by_factor),
    labels = levels(by_factor),
    levels = levels(by_factor)
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_observed_groups
# ---------------------------------------------------------------------------- #
#
# Assign observed studies to plotting groups.
#
# ---------------------------------------------------------------------------- #
.regplot_observed_groups <- function(by_info, groups, n) {

  if (is.null(by_info)) {
    return(list(
      group    = rep("All", n),
      group_id = rep(1L, n)
    ))
  }

  by_data <- by_info[["data"]]

  if (by_info[["type"]] == "categorical") {
    group <- as.character(factor(by_data, levels = groups))
  } else {
    by_mean <- mean(by_data)
    by_sd   <- stats::sd(by_data)

    if (!is.finite(by_sd) || by_sd == 0 || length(groups) == 1L) {
      group <- rep(groups[1], length(by_data))
    } else {
      by_values <- c(by_mean - by_sd, by_mean, by_mean + by_sd)
      breaks    <- c(-Inf, mean(by_values[1:2]), mean(by_values[2:3]), Inf)
      group     <- as.character(cut(by_data, breaks = breaks, labels = groups))
    }
  }

  return(list(
    group    = group,
    group_id = match(group, groups)
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_add_dummy_outcome
# ---------------------------------------------------------------------------- #
#
# Add outcome columns for prediction paths that use sampling bias inputs.
#
# ---------------------------------------------------------------------------- #
.regplot_add_dummy_outcome <- function(x, newdata, n_pred, sei) {

  outcome_type <- .outcome_type(x)

  if (outcome_type == "norm") {
    newdata[["yi"]]  <- rep(0, n_pred)
    newdata[["sei"]] <- rep(sei, n_pred)
  } else if (outcome_type == "bin") {
    newdata[["ai"]]  <- rep(0, n_pred)
    newdata[["ci"]]  <- rep(0, n_pred)
    newdata[["n1i"]] <- rep(0, n_pred)
    newdata[["n2i"]] <- rep(0, n_pred)
  } else if (outcome_type == "pois") {
    newdata[["x1i"]] <- rep(0, n_pred)
    newdata[["x2i"]] <- rep(0, n_pred)
    newdata[["t1i"]] <- rep(0, n_pred)
    newdata[["t2i"]] <- rep(0, n_pred)
  }

  return(newdata)
}


# ---------------------------------------------------------------------------- #
# .regplot_band_data_categorical
# ---------------------------------------------------------------------------- #
#
# Categorical intervals are rendered as box-style summaries.
#
# ---------------------------------------------------------------------------- #
.regplot_band_data_categorical <- function(grid, middle, lower, upper) {

  out <- data.frame(
    x        = grid[["x"]],
    y        = middle,
    middle   = middle,
    lower    = lower,
    upper    = upper,
    group    = grid[["group"]],
    group_id = grid[["group_id"]],
    level    = grid[["level"]],
    stringsAsFactors = FALSE
  )

  return(out)
}


# ---------------------------------------------------------------------------- #
# .regplot_data
# ---------------------------------------------------------------------------- #
#
# Generate data for regression plot.
#
# @param x        brma object
# @param mod_name name of moderator variable
# @param mod_type "continuous" or "categorical"
# @param mod_data moderator data values
# @param by_info  list with grouping variable info (or NULL)
# @param pred     logical; show prediction line
# @param ci       logical; show CI bands
# @param pi       logical; show PI bands
# @param si       logical; show SI bands
# @param level    confidence level (0-100)
# @param at       values at which to evaluate predictions
# @param psize    point sizes (or NULL for auto)
# @param plim     point size range
# @param transf   transformation function
# @param xlim     x-axis limits
# @param ylim     y-axis limits
# @param xlab     x-axis label
# @param ylab     y-axis label
# @param refline       reference line position
# @param sampling_bias logical; incorporate bias into predictions
# @param reference_sei numeric; reference standard error for sampling paths
# @param dots          graphical parameters
#
# @return list with plot data components
#
# ---------------------------------------------------------------------------- #
.regplot_data <- function(x, mod_name, mod_type, mod_data, by_info,
                          pred, ci, pi, si, level, at, digits, psize, plim,
                          transf, xlim, ylim, xlab, ylab, refline,
                          sampling_bias, reference_sei, dots) {

  yi      <- .outcome_data_yi(x)
  sei_obs <- .outcome_data_sei(x)
  vi      <- sei_obs^2
  K       <- length(yi)
  se_rep  <- if (is.null(reference_sei)) stats::median(sei_obs) else reference_sei

  alpha <- (100 - level) / 100
  probs <- c(alpha / 2, 1 - alpha / 2)

  if (is.null(psize)) {
    weights   <- 1 / vi
    weights_n <- (weights - min(weights)) / (max(weights) - min(weights) + 1e-10)
    psize     <- plim[1] + weights_n * (plim[2] - plim[1])
  } else if (length(psize) == 1) {
    psize <- rep(psize, K)
  }

  yi_plot <- yi
  if (!is.null(transf)) {
    yi_plot <- transf(yi)
  }

  prediction_grid <- .regplot_prediction_grid(
    x         = x,
    mod_name  = mod_name,
    mod_type  = mod_type,
    mod_data  = mod_data,
    by_info   = by_info,
    at        = at,
    digits    = digits
  )
  grid   <- prediction_grid[["grid"]]
  groups <- prediction_grid[["groups"]]
  n_pred <- nrow(grid)

  prediction_se <- if (sampling_bias) se_rep else 0
  newdata <- .regplot_add_dummy_outcome(
    x       = x,
    newdata = prediction_grid[["newdata"]],
    n_pred  = n_pred,
    sei     = prediction_se
  )

  pred_samples <- predict.brma(
    object        = x,
    newdata       = newdata,
    type          = "terms",
    bias_adjusted = !sampling_bias,
    quiet         = TRUE
  )
  pred_samples <- as.matrix(pred_samples)

  pred_mean  <- colMeans(pred_samples)
  pred_lower <- apply(pred_samples, 2, stats::quantile, probs = probs[1], names = FALSE)
  pred_upper <- apply(pred_samples, 2, stats::quantile, probs = probs[2], names = FALSE)

  if (pi || si) {
    is_scale      <- .is_scale(x)
    is_multilevel <- .is_multilevel(x)

    scale_data    <- NULL
    scale_formula <- NULL
    if (is_scale) {
      new_data_prepared <- .prepare_newdata(object = x, newdata = newdata, type = "estimate")
      scale_data        <- new_data_prepared[["scale"]]
      scale_formula     <- .create_fit_formula_list(data = new_data_prepared, "scale")
    }

    posterior_samples <- .get_posterior_samples(x[["fit"]])
    tau_result <- .evaluate.brma.tau(
      fit               = x[["fit"]],
      scale_data        = scale_data,
      scale_formula     = scale_formula,
      scale_priors      = x[["priors"]][["scale"]],
      is_scale          = is_scale,
      is_multilevel     = is_multilevel,
      K                 = n_pred,
      posterior_samples = posterior_samples
    )
    tau_total <- sqrt(tau_result[["tau_within"]]^2 + tau_result[["tau_between"]]^2)
  }

  if (pi) {
    pi_bounds <- .regplot_mixture_interval_quantiles(
      mean_samples = pred_samples,
      sd_samples   = tau_total,
      probs        = probs
    )
    pi_lower <- pi_bounds[["lower"]]
    pi_upper <- pi_bounds[["upper"]]
  } else {
    pi_lower <- NULL
    pi_upper <- NULL
  }

  if (si) {
    sd_si <- sqrt(tau_total^2 + se_rep^2)

    if (sampling_bias && .is_weightfunction(x)) {
      si_bounds <- .regplot_selection_mixture_interval_quantiles(
        x            = x,
        mean_samples = pred_samples,
        sd_samples   = sd_si,
        se           = se_rep,
        probs        = probs
      )
    } else {
      si_bounds <- .regplot_mixture_interval_quantiles(
        mean_samples = pred_samples,
        sd_samples   = sd_si,
        probs        = probs
      )
    }

    si_lower <- si_bounds[["lower"]]
    si_upper <- si_bounds[["upper"]]
  } else {
    si_lower <- NULL
    si_upper <- NULL
  }

  if (!is.null(transf)) {
    pred_mean  <- transf(pred_mean)
    ci_lo      <- transf(pred_lower)
    ci_hi      <- transf(pred_upper)
    pred_lower <- pmin(ci_lo, ci_hi)
    pred_upper <- pmax(ci_lo, ci_hi)

    if (pi) {
      pi_lo    <- transf(pi_lower)
      pi_hi    <- transf(pi_upper)
      pi_lower <- pmin(pi_lo, pi_hi)
      pi_upper <- pmax(pi_lo, pi_hi)
    }
    if (si) {
      si_lo    <- transf(si_lower)
      si_hi    <- transf(si_upper)
      si_lower <- pmin(si_lo, si_hi)
      si_upper <- pmax(si_lo, si_hi)
    }
  }

  if (is.null(xlim)) {
    if (mod_type == "continuous") {
      xlim <- range(pretty(range(mod_data)))
    } else {
      xlim <- c(0.5, length(levels(mod_data)) + 0.5)
    }
  }

  if (is.null(ylim)) {
    all_y <- c(yi_plot, pred_mean)
    if (ci) all_y <- c(all_y, pred_lower, pred_upper)
    if (pi) all_y <- c(all_y, pi_lower, pi_upper)
    if (si) all_y <- c(all_y, si_lower, si_upper)
    ylim <- range(pretty(range(all_y, na.rm = TRUE)))
  }

  if (is.null(xlab)) xlab <- mod_name
  if (is.null(ylab)) ylab <- "Observed Effect Size"

  observed_groups <- .regplot_observed_groups(
    by_info = by_info,
    groups  = groups,
    n       = K
  )

  df_points <- data.frame(
    x        = if (mod_type == "continuous") mod_data else as.numeric(mod_data),
    y        = yi_plot,
    size     = psize,
    group    = observed_groups[["group"]],
    group_id = observed_groups[["group_id"]],
    stringsAsFactors = FALSE
  )

  if (mod_type == "categorical") {
    df_points$level <- mod_data
  }

  df_pred <- NULL
  if (pred) {
    df_pred <- data.frame(
      x        = grid[["x"]],
      y        = pred_mean,
      group    = grid[["group"]],
      group_id = grid[["group_id"]],
      stringsAsFactors = FALSE
    )
    if (mod_type == "categorical") {
      df_pred[["level"]] <- grid[["level"]]
    }
  }

  df_ci <- NULL
  if (ci) {
    if (mod_type == "continuous") {
      df_ci <- do.call(rbind, lapply(seq_along(groups), function(i) {
        rows <- grid[["group_id"]] == i
        .regplot_band_data_continuous(
          xpred    = grid[["x"]][rows],
          lower    = pred_lower[rows],
          upper    = pred_upper[rows],
          group    = groups[i],
          group_id = i
        )
      }))
      rownames(df_ci) <- NULL
    } else {
      df_ci <- .regplot_band_data_categorical(grid, pred_mean, pred_lower, pred_upper)
    }
  }

  df_pi <- NULL
  if (pi) {
    if (mod_type == "continuous") {
      df_pi <- do.call(rbind, lapply(seq_along(groups), function(i) {
        rows <- grid[["group_id"]] == i
        .regplot_band_data_continuous(
          xpred    = grid[["x"]][rows],
          lower    = pi_lower[rows],
          upper    = pi_upper[rows],
          group    = groups[i],
          group_id = i
        )
      }))
      rownames(df_pi) <- NULL
    } else {
      df_pi <- .regplot_band_data_categorical(grid, pred_mean, pi_lower, pi_upper)
    }
  }

  df_si <- NULL
  if (si) {
    if (mod_type == "continuous") {
      df_si <- do.call(rbind, lapply(seq_along(groups), function(i) {
        rows <- grid[["group_id"]] == i
        .regplot_band_data_continuous(
          xpred    = grid[["x"]][rows],
          lower    = si_lower[rows],
          upper    = si_upper[rows],
          group    = groups[i],
          group_id = i
        )
      }))
      rownames(df_si) <- NULL
    } else {
      df_si <- .regplot_band_data_categorical(grid, pred_mean, si_lower, si_upper)
    }
  }

  df_refline <- NULL
  if (!is.null(refline)) {
    df_refline <- data.frame(y = refline)
  }

  return(list(
    points   = df_points,
    pred     = df_pred,
    ci       = df_ci,
    pi       = df_pi,
    si       = df_si,
    refline  = df_refline,
    xlim     = xlim,
    ylim     = ylim,
    xlab     = xlab,
    ylab     = ylab,
    mod_type = mod_type,
    mod_name = mod_name,
    by_name  = if (!is.null(by_info)) by_info[["name"]] else NULL,
    by_type  = if (!is.null(by_info)) by_info[["type"]] else NULL,
    groups   = groups,
    levels   = if (mod_type == "categorical") levels(mod_data) else NULL,
    sei      = se_rep
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_mixture_interval_quantiles
# ---------------------------------------------------------------------------- #
#
# Deterministic quantiles of posterior mixtures of normal distributions.
#
# ---------------------------------------------------------------------------- #
.has_native_regplot_mixture <- function() {

  return(is.loaded("RoBMA_regplot_normal_mixture_interval", PACKAGE = "RoBMA"))
}

.has_native_regplot_selection_mixture <- function() {

  return(is.loaded("RoBMA_regplot_selnorm_mixture_interval", PACKAGE = "RoBMA"))
}


.regplot_mixture_interval_quantiles <- function(mean_samples, sd_samples, probs) {

  mean_samples <- as.matrix(mean_samples)
  sd_samples   <- as.matrix(sd_samples)

  if (.has_native_regplot_mixture()) {
    return(.Call(
      "RoBMA_regplot_normal_mixture_interval",
      .native_numeric_matrix(mean_samples),
      .native_numeric_matrix(sd_samples),
      .native_numeric_vector(probs),
      PACKAGE = "RoBMA"
    ))
  }

  return(.regplot_mixture_interval_quantiles_r(
    mean_samples = mean_samples,
    sd_samples   = sd_samples,
    probs        = probs
  ))
}


.regplot_mixture_interval_quantiles_r <- function(mean_samples, sd_samples, probs) {

  lower <- numeric(ncol(mean_samples))
  upper <- numeric(ncol(mean_samples))

  for (i in seq_len(ncol(mean_samples))) {
    cdf_fun <- function(q) {
      .regplot_normal_mixture_cdf(
        q    = q,
        mean = mean_samples[, i],
        sd   = sd_samples[, i]
      )
    }

    lower[i] <- .regplot_mixture_quantile(probs[1], mean_samples[, i], sd_samples[, i], cdf_fun)
    upper[i] <- .regplot_mixture_quantile(probs[2], mean_samples[, i], sd_samples[, i], cdf_fun)
  }

  return(list(lower = lower, upper = upper))
}


# ---------------------------------------------------------------------------- #
# .regplot_selection_mixture_interval_quantiles
# ---------------------------------------------------------------------------- #
#
# Deterministic quantiles of observed-effect mixtures with selection branches.
#
# ---------------------------------------------------------------------------- #
.regplot_selection_mixture_interval_quantiles <- function(x, mean_samples,
                                                          sd_samples, se,
                                                          probs) {

  mean_samples     <- as.matrix(mean_samples)
  sd_samples       <- as.matrix(sd_samples)
  setup            <- .regplot_selection_setup(x)
  effect_direction <- .effect_direction(x)

  if (!is.null(setup[["selection"]]) &&
      se > 0 &&
      .has_native_regplot_selection_mixture()) {
    return(.regplot_selnorm_mixture_interval_quantiles(
      mean_samples      = mean_samples,
      sd_samples        = sd_samples,
      se                = se,
      probs             = probs,
      selection_context = setup[["selection"]]
    ))
  }

  return(.regplot_selection_mixture_interval_quantiles_r(
    mean_samples     = mean_samples,
    sd_samples       = sd_samples,
    se               = se,
    probs            = probs,
    setup            = setup,
    effect_direction = effect_direction
  ))
}

.regplot_selnorm_mixture_interval_quantiles <- function(mean_samples,
                                                        sd_samples, se,
                                                        probs,
                                                        selection_context) {

  .selection_require_step_evaluable(
    selection_context,
    ".regplot_selection_mixture_interval_quantiles()"
  )

  return(.Call(
    "RoBMA_regplot_selnorm_mixture_interval",
    .native_numeric_matrix(mean_samples),
    .native_numeric_matrix(sd_samples),
    .native_numeric_vector(se),
    .native_numeric_vector(probs),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    .native_numeric_vector(selection_context[["z_lower"]]),
    .native_numeric_vector(selection_context[["z_upper"]]),
    .native_integer_vector(selection_context[["sign"]]),
    .native_integer_vector(selection_context[["phack_q"]]),
    .native_numeric_vector(selection_context[["phack_z_source"]]),
    .native_numeric_vector(selection_context[["phack_z_dest"]]),
    .native_numeric_vector(selection_context[["segments"]][["bounds"]]),
    .native_integer_vector(selection_context[["segments"]][["step_bin"]]),
    .native_integer_vector(selection_context[["segments"]][["phack_region"]]),
    PACKAGE = "RoBMA"
  ))
}


.regplot_selection_mixture_interval_quantiles_r <- function(mean_samples,
                                                            sd_samples, se,
                                                            probs, setup,
                                                            effect_direction) {

  lower <- numeric(ncol(mean_samples))
  upper <- numeric(ncol(mean_samples))

  for (i in seq_len(ncol(mean_samples))) {
    cdf_fun <- function(q) {
      .regplot_selection_mixture_cdf(
        q                = q,
        mean             = mean_samples[, i],
        sd               = sd_samples[, i],
        se               = se,
        setup            = setup,
        effect_direction = effect_direction
      )
    }

    lower[i] <- .regplot_mixture_quantile(probs[1], mean_samples[, i], sd_samples[, i], cdf_fun)
    upper[i] <- .regplot_mixture_quantile(probs[2], mean_samples[, i], sd_samples[, i], cdf_fun)
  }

  return(list(lower = lower, upper = upper))
}


# ---------------------------------------------------------------------------- #
# .regplot_mixture_quantile
# ---------------------------------------------------------------------------- #
#
# Invert a posterior-mixture CDF.
#
# ---------------------------------------------------------------------------- #
.regplot_mixture_quantile <- function(p, mean, sd, cdf_fun) {

  eps_sd <- sqrt(.Machine$double.eps)

  if (all(sd < eps_sd)) {
    return(unname(stats::quantile(mean, probs = p, names = FALSE, type = 8)))
  }

  spread <- pmax(sd, eps_sd)
  lower  <- min(mean - 10 * spread, na.rm = TRUE)
  upper  <- max(mean + 10 * spread, na.rm = TRUE)

  if (!is.finite(lower) || !is.finite(upper)) {
    return(NA_real_)
  }
  if (lower >= upper) {
    lower <- lower - 1
    upper <- upper + 1
  }

  obj_fun     <- function(q) cdf_fun(q) - p
  lower_value <- obj_fun(lower)
  upper_value <- obj_fun(upper)
  step        <- max(spread, na.rm = TRUE)

  if (!is.finite(step) || step <= 0) {
    step <- max(1, abs(mean), na.rm = TRUE)
  }

  for (i in seq_len(25)) {
    if (lower_value <= 0 && upper_value >= 0) {
      break
    }
    if (lower_value > 0) {
      lower       <- lower - step
      lower_value <- obj_fun(lower)
    }
    if (upper_value < 0) {
      upper       <- upper + step
      upper_value <- obj_fun(upper)
    }
    step <- step * 2
  }

  if (lower_value > 0 || upper_value < 0) {
    return(.regplot_grid_quantile(p, lower, upper, cdf_fun))
  }

  out <- tryCatch(
    stats::uniroot(obj_fun, interval = c(lower, upper), tol = 1e-6)[["root"]],
    error = function(e) NA_real_
  )

  if (is.na(out)) {
    out <- .regplot_grid_quantile(p, lower, upper, cdf_fun)
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .regplot_normal_mixture_cdf
# ---------------------------------------------------------------------------- #
.regplot_normal_mixture_cdf <- function(q, mean, sd) {

  eps_sd     <- sqrt(.Machine$double.eps)
  cdf_values <- rep(NA_real_, length(mean))
  zero_sd    <- sd < eps_sd

  if (any(zero_sd)) {
    cdf_values[zero_sd] <- as.numeric(q >= mean[zero_sd])
  }
  if (any(!zero_sd)) {
    cdf_values[!zero_sd] <- stats::pnorm(
      q,
      mean = mean[!zero_sd],
      sd   = sd[!zero_sd]
    )
  }

  cdf_values <- pmin(pmax(cdf_values, 0), 1)
  return(base::mean(cdf_values))
}


# ---------------------------------------------------------------------------- #
# .regplot_selection_mixture_cdf
# ---------------------------------------------------------------------------- #
.regplot_selection_mixture_cdf <- function(q, mean, sd, se, setup,
                                           effect_direction) {

  eps_sd     <- sqrt(.Machine$double.eps)
  cdf_values <- rep(NA_real_, length(mean))
  zero_sd    <- sd < eps_sd

  if (any(zero_sd)) {
    cdf_values[zero_sd] <- as.numeric(q >= mean[zero_sd])
  }

  normal_rows <- !setup[["is_weightfunction"]] & !zero_sd
  if (any(normal_rows)) {
    cdf_values[normal_rows] <- stats::pnorm(
      q,
      mean = mean[normal_rows],
      sd   = sd[normal_rows]
    )
  }

  selected_rows <- setup[["is_weightfunction"]] & !zero_sd
  if (any(selected_rows)) {
    setup[["mu"]] <- mean
    rows <- which(selected_rows)
    cdf_values[rows] <- .funnel_selected_cdf(
      q                = q,
      rows             = rows,
      se               = se,
      total_sd         = sd,
      setup            = setup,
      effect_direction = effect_direction
    )
  }

  cdf_values <- pmin(pmax(cdf_values, 0), 1)
  return(base::mean(cdf_values))
}


# ---------------------------------------------------------------------------- #
# .regplot_selection_setup
# ---------------------------------------------------------------------------- #
.regplot_selection_setup <- function(x) {

  posterior_samples <- .get_posterior_samples(x[["fit"]])
  S                 <- nrow(posterior_samples)
  bias_indicator           <- .extract_bias_indicator(x, posterior_samples = posterior_samples)
  selection                <- .selection_context(
    object            = x,
    posterior_samples = posterior_samples
  )
  use_normal <- if (is.null(selection)) {
    rep(TRUE, S)
  } else {
    selection[["use_normal"]]
  }

  return(list(
    bias_indicator    = bias_indicator,
    is_weightfunction = !use_normal,
    selection         = selection
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_grid_quantile
# ---------------------------------------------------------------------------- #
.regplot_grid_quantile <- function(p, lower, upper, cdf_fun) {

  grid <- seq(lower, upper, length.out = 1000)
  cdf  <- vapply(grid, cdf_fun, numeric(1))
  index <- which(cdf >= p)[1]

  if (is.na(index)) {
    return(grid[length(grid)])
  }

  return(grid[index])
}


# ---------------------------------------------------------------------------- #
# Plotting helpers
# ---------------------------------------------------------------------------- #
.regplot_has_shade <- function(shade) {
  return(isTRUE(shade))
}

.regplot_palette <- function(groups, default_col) {

  n <- length(groups)
  if (n == 1L) {
    return(stats::setNames(default_col[1], groups))
  }

  if (length(default_col) >= n && length(unique(default_col)) > 1L) {
    cols <- default_col[seq_len(n)]
  } else {
    cols <- grDevices::hcl.colors(n, palette = "Dark 3")
  }

  return(stats::setNames(cols, groups))
}

.regplot_dodge_x <- function(x, group_id, n_groups, dodge_width) {

  if (n_groups <= 1L) {
    return(x)
  }

  offset <- (group_id - (n_groups + 1) / 2) * dodge_width / n_groups
  return(x + offset)
}

.regplot_jitter_values <- function(n, amount) {

  if (amount <= 0) {
    return(rep(0, n))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  set.seed(42)
  return(stats::runif(n, -amount, amount))
}

.regplot_band_edges <- function(df_band) {

  out <- do.call(rbind, lapply(split(df_band, df_band[["group"]]), function(df_group) {
    n <- nrow(df_group) / 2
    data.frame(
      x        = rep(df_group[["x"]][seq_len(n)], 2),
      y        = c(df_group[["lower"]][seq_len(n)], df_group[["upper"]][seq_len(n)]),
      edge     = rep(c("lower", "upper"), each = n),
      group    = rep(df_group[["group"]][1], 2 * n),
      group_id = rep(df_group[["group_id"]][1], 2 * n),
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  return(out)
}

.regplot_draw_continuous_band_base <- function(df_band, fill_col, alpha,
                                               border_col, shade, lwd) {

  if (is.null(df_band)) {
    return(invisible(NULL))
  }

  for (group in unique(df_band[["group"]])) {
    df_group <- df_band[df_band[["group"]] == group, , drop = FALSE]
    n        <- nrow(df_group) / 2

    if (.regplot_has_shade(shade)) {
      graphics::polygon(
        df_group[["x"]],
        df_group[["y"]],
        col    = grDevices::adjustcolor(fill_col[group], alpha.f = alpha),
        border = NA
      )
    }

    graphics::lines(
      df_group[["x"]][seq_len(n)],
      df_group[["lower"]][seq_len(n)],
      col = border_col[group],
      lwd = max(1, lwd / 2)
    )
    graphics::lines(
      df_group[["x"]][seq_len(n)],
      df_group[["upper"]][seq_len(n)],
      col = border_col[group],
      lwd = max(1, lwd / 2)
    )
  }

  return(invisible(NULL))
}

.regplot_draw_categorical_band_base <- function(df_band, fill_col, alpha,
                                                border_col, shade, width,
                                                n_groups, dodge_width, lwd) {

  if (is.null(df_band)) {
    return(invisible(NULL))
  }

  x <- .regplot_dodge_x(df_band[["x"]], df_band[["group_id"]], n_groups, dodge_width)
  w <- width / max(1, n_groups)

  for (i in seq_len(nrow(df_band))) {
    group <- df_band[["group"]][i]
    graphics::rect(
      xleft   = x[i] - w / 2,
      ybottom = df_band[["lower"]][i],
      xright  = x[i] + w / 2,
      ytop    = df_band[["upper"]][i],
      col     = if (.regplot_has_shade(shade)) grDevices::adjustcolor(fill_col[group], alpha.f = alpha) else NA,
      border  = border_col[group],
      lwd     = max(1, lwd / 2)
    )
    graphics::segments(
      x0  = x[i] - w / 2,
      y0  = df_band[["middle"]][i],
      x1  = x[i] + w / 2,
      y1  = df_band[["middle"]][i],
      col = border_col[group],
      lwd = max(1, lwd / 2)
    )
  }

  return(invisible(NULL))
}


.regplot_add_continuous_band_ggplot <- function(out, df_band, fill_col, alpha,
                                                border_col, shade, has_by) {

  if (is.null(df_band)) {
    return(out)
  }

  if (.regplot_has_shade(shade)) {
    if (has_by) {
      out <- out +
        ggplot2::geom_polygon(
          data    = df_band,
          mapping = ggplot2::aes(
            x     = .data[["x"]],
            y     = .data[["y"]],
            group = .data[["group"]],
            fill  = .data[["group"]]
          ),
          alpha   = alpha
        )
    } else {
      out <- out +
        ggplot2::geom_polygon(
          data    = df_band,
          mapping = ggplot2::aes(
            x     = .data[["x"]],
            y     = .data[["y"]],
            group = .data[["group"]]
          ),
          fill    = fill_col[1],
          alpha   = alpha
        )
    }
  }

  df_edges <- .regplot_band_edges(df_band)
  if (has_by) {
    out <- out +
      ggplot2::geom_line(
        data    = df_edges,
        mapping = ggplot2::aes(
          x      = .data[["x"]],
          y      = .data[["y"]],
          group  = interaction(.data[["group"]], .data[["edge"]]),
          colour = .data[["group"]]
        ),
        linewidth = 0.35
      )
  } else {
    out <- out +
      ggplot2::geom_line(
        data    = df_edges,
        mapping = ggplot2::aes(
          x     = .data[["x"]],
          y     = .data[["y"]],
          group = .data[["edge"]]
        ),
        colour  = border_col[1],
        linewidth = 0.35
      )
  }

  return(out)
}

.regplot_add_categorical_band_ggplot <- function(out, df_band, fill_col, alpha,
                                                 border_col, shade, has_by,
                                                 width, n_groups, dodge_width) {

  if (is.null(df_band)) {
    return(out)
  }

  df_band[["x_plot"]] <- .regplot_dodge_x(
    df_band[["x"]],
    df_band[["group_id"]],
    n_groups,
    dodge_width
  )
  width <- width / max(1, n_groups)

  if (has_by && .regplot_has_shade(shade)) {
    out <- out +
      ggplot2::geom_crossbar(
        data    = df_band,
        mapping = ggplot2::aes(
          x      = .data[["x_plot"]],
          y      = .data[["middle"]],
          ymin   = .data[["lower"]],
          ymax   = .data[["upper"]],
          colour = .data[["group"]],
          fill   = .data[["group"]]
        ),
        width   = width,
        alpha   = alpha,
        linewidth = 0.35
      )
  } else if (has_by) {
    out <- out +
      ggplot2::geom_crossbar(
        data    = df_band,
        mapping = ggplot2::aes(
          x      = .data[["x_plot"]],
          y      = .data[["middle"]],
          ymin   = .data[["lower"]],
          ymax   = .data[["upper"]],
          colour = .data[["group"]]
        ),
        width   = width,
        fill    = NA,
        alpha   = 1,
        linewidth = 0.35
      )
  } else {
    out <- out +
      ggplot2::geom_crossbar(
        data    = df_band,
        mapping = ggplot2::aes(
          x    = .data[["x_plot"]],
          y    = .data[["middle"]],
          ymin = .data[["lower"]],
          ymax = .data[["upper"]]
        ),
        width   = width,
        colour  = border_col[1],
        fill    = if (.regplot_has_shade(shade)) fill_col[1] else NA,
        alpha   = if (.regplot_has_shade(shade)) alpha else 1,
        linewidth = 0.35
      )
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .regplot_plot_base
# ---------------------------------------------------------------------------- #
#
# Create regression plot using base R graphics.
#
# @param data   list of plot data from .regplot_data
# @param dots   list of graphical parameters
# @param legend logical; show legend
#
# @return NULL invisibly
#
# ---------------------------------------------------------------------------- #
.regplot_plot_base <- function(data, dots, legend = TRUE) {

  pch       <- dots[["pch"]]
  col       <- dots[["col"]]
  bg        <- dots[["bg"]]
  lcol      <- dots[["lcol"]]
  lwd       <- dots[["lwd"]]
  col_ci    <- dots[["col.ci"]]
  col_pi    <- dots[["col.pi"]]
  col_si    <- dots[["col.si"]]
  alpha_ci  <- dots[["alpha.ci"]]
  alpha_pi  <- dots[["alpha.pi"]]
  alpha_si  <- dots[["alpha.si"]]
  shade     <- dots[["shade"]]
  main      <- dots[["main"]]
  jitter_am <- dots[["jitter"]]
  box_width <- dots[["box.width"]]
  dodge_w   <- dots[["dodge.width"]]
  las       <- dots[["las"]]

  df_points  <- data$points
  df_pred    <- data$pred
  df_ci      <- data$ci
  df_pi      <- data$pi
  df_si      <- data$si
  df_refline <- data$refline
  mod_type   <- data$mod_type
  groups     <- data$groups
  n_groups   <- length(groups)
  has_by     <- n_groups > 1L

  line_cols <- .regplot_palette(groups, lcol)
  point_cols <- if (has_by) line_cols else stats::setNames(rep(col, n_groups), groups)
  point_bgs  <- if (has_by) {
    stats::setNames(grDevices::adjustcolor(line_cols, alpha.f = 0.45), groups)
  } else {
    stats::setNames(rep(bg, n_groups), groups)
  }

  if (has_by) {
    ci_cols <- pi_cols <- si_cols <- line_cols
  } else {
    ci_cols <- stats::setNames(rep(col_ci, n_groups), groups)
    pi_cols <- stats::setNames(rep(col_pi, n_groups), groups)
    si_cols <- stats::setNames(rep(col_si, n_groups), groups)
  }

  if (mod_type == "categorical") {
    graphics::plot(
      NA, NA,
      xlim = data$xlim,
      ylim = data$ylim,
      xlab = data$xlab,
      ylab = data$ylab,
      main = main,
      type = "n",
      xaxt = "n",
      las  = las
    )
    graphics::axis(
      1,
      at     = seq_along(data$levels),
      labels = data$levels,
      las    = las
    )
  } else {
    graphics::plot(
      NA, NA,
      xlim = data$xlim,
      ylim = data$ylim,
      xlab = data$xlab,
      ylab = data$ylab,
      main = main,
      type = "n",
      las  = las
    )
  }

  if (!is.null(df_refline)) {
    graphics::abline(h = df_refline$y, lty = "dashed", col = "gray50")
  }

  if (!is.null(df_si)) {
    if (mod_type == "continuous") {
      .regplot_draw_continuous_band_base(df_si, si_cols, alpha_si, si_cols, shade, lwd)
    } else {
      .regplot_draw_categorical_band_base(
        df_si, si_cols, alpha_si, si_cols, shade,
        box_width * 1.4, n_groups, dodge_w, lwd
      )
    }
  }

  if (!is.null(df_pi)) {
    if (mod_type == "continuous") {
      .regplot_draw_continuous_band_base(df_pi, pi_cols, alpha_pi, pi_cols, shade, lwd)
    } else {
      .regplot_draw_categorical_band_base(
        df_pi, pi_cols, alpha_pi, pi_cols, shade,
        box_width * 1.15, n_groups, dodge_w, lwd
      )
    }
  }

  if (!is.null(df_ci)) {
    if (mod_type == "continuous") {
      .regplot_draw_continuous_band_base(df_ci, ci_cols, alpha_ci, ci_cols, shade, lwd)
    } else {
      .regplot_draw_categorical_band_base(
        df_ci, ci_cols, alpha_ci, ci_cols, shade,
        box_width, n_groups, dodge_w, lwd
      )
    }
  }

  if (!is.null(df_pred)) {
    if (mod_type == "continuous") {
      for (group in groups) {
        df_group <- df_pred[df_pred[["group"]] == group, , drop = FALSE]
        graphics::lines(df_group$x, df_group$y, col = line_cols[group], lwd = lwd)
      }
    } else {
      x_pred <- .regplot_dodge_x(df_pred$x, df_pred$group_id, n_groups, dodge_w)
      graphics::points(
        x_pred,
        df_pred$y,
        pch = 18,
        col = line_cols[df_pred$group],
        cex = 1.5
      )
    }
  }

  if (mod_type == "categorical") {
    x_points <- .regplot_dodge_x(df_points$x, df_points$group_id, n_groups, dodge_w)
    x_points <- x_points + .regplot_jitter_values(
      n      = nrow(df_points),
      amount = jitter_am / max(1, n_groups)
    )
    graphics::points(
      x_points,
      df_points$y,
      pch = pch,
      col = point_cols[df_points$group],
      bg  = point_bgs[df_points$group],
      cex = df_points$size
    )
  } else {
    graphics::points(
      df_points$x,
      df_points$y,
      pch = pch,
      col = point_cols[df_points$group],
      bg  = point_bgs[df_points$group],
      cex = df_points$size
    )
  }

  if (has_by && legend) {
    graphics::legend(
      "topright",
      legend = groups,
      col    = line_cols,
      pt.bg  = point_bgs,
      pch    = 21,
      lty    = 1,
      lwd    = lwd,
      bty    = "n"
    )
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .regplot_plot_ggplot
# ---------------------------------------------------------------------------- #
#
# Create regression plot using ggplot2.
#
# @param data   list of plot data from .regplot_data
# @param dots   list of graphical parameters
# @param legend logical; show legend
#
# @return ggplot object
#
# ---------------------------------------------------------------------------- #
.regplot_plot_ggplot <- function(data, dots, legend = TRUE) {

  pch       <- dots[["pch"]]
  col       <- dots[["col"]]
  bg        <- dots[["bg"]]
  lcol      <- dots[["lcol"]]
  lwd       <- dots[["lwd"]]
  col_ci    <- dots[["col.ci"]]
  col_pi    <- dots[["col.pi"]]
  col_si    <- dots[["col.si"]]
  alpha_ci  <- dots[["alpha.ci"]]
  alpha_pi  <- dots[["alpha.pi"]]
  alpha_si  <- dots[["alpha.si"]]
  shade     <- dots[["shade"]]
  main      <- dots[["main"]]
  jitter_am <- dots[["jitter"]]
  box_width <- dots[["box.width"]]
  dodge_w   <- dots[["dodge.width"]]

  df_points  <- data$points
  df_pred    <- data$pred
  df_ci      <- data$ci
  df_pi      <- data$pi
  df_si      <- data$si
  df_refline <- data$refline
  mod_type   <- data$mod_type
  groups     <- data$groups
  n_groups   <- length(groups)
  has_by     <- n_groups > 1L

  line_cols <- .regplot_palette(groups, lcol)
  if (has_by) {
    ci_cols <- pi_cols <- si_cols <- line_cols
  } else {
    ci_cols <- stats::setNames(rep(col_ci, n_groups), groups)
    pi_cols <- stats::setNames(rep(col_pi, n_groups), groups)
    si_cols <- stats::setNames(rep(col_si, n_groups), groups)
  }

  out <- ggplot2::ggplot()

  if (!is.null(df_refline)) {
    out <- out +
      ggplot2::geom_hline(
        yintercept = df_refline$y,
        linetype   = "dashed",
        colour     = "gray50"
      )
  }

  if (!is.null(df_si)) {
    if (mod_type == "continuous") {
      out <- .regplot_add_continuous_band_ggplot(out, df_si, si_cols, alpha_si, si_cols, shade, has_by)
    } else {
      out <- .regplot_add_categorical_band_ggplot(
        out, df_si, si_cols, alpha_si, si_cols, shade, has_by,
        box_width * 1.4, n_groups, dodge_w
      )
    }
  }

  if (!is.null(df_pi)) {
    if (mod_type == "continuous") {
      out <- .regplot_add_continuous_band_ggplot(out, df_pi, pi_cols, alpha_pi, pi_cols, shade, has_by)
    } else {
      out <- .regplot_add_categorical_band_ggplot(
        out, df_pi, pi_cols, alpha_pi, pi_cols, shade, has_by,
        box_width * 1.15, n_groups, dodge_w
      )
    }
  }

  if (!is.null(df_ci)) {
    if (mod_type == "continuous") {
      out <- .regplot_add_continuous_band_ggplot(out, df_ci, ci_cols, alpha_ci, ci_cols, shade, has_by)
    } else {
      out <- .regplot_add_categorical_band_ggplot(
        out, df_ci, ci_cols, alpha_ci, ci_cols, shade, has_by,
        box_width, n_groups, dodge_w
      )
    }
  }

  if (!is.null(df_pred)) {
    if (mod_type == "continuous") {
      if (has_by) {
        out <- out +
          ggplot2::geom_line(
            data    = df_pred,
            mapping = ggplot2::aes(
              x      = .data[["x"]],
              y      = .data[["y"]],
              colour = .data[["group"]],
              group  = .data[["group"]]
            ),
            linewidth = lwd / 2
          )
      } else {
        out <- out +
          ggplot2::geom_line(
            data    = df_pred,
            mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
            colour  = lcol[1],
            linewidth = lwd / 2
          )
      }
    } else {
      df_pred[["x_plot"]] <- .regplot_dodge_x(df_pred$x, df_pred$group_id, n_groups, dodge_w)
      if (has_by) {
        out <- out +
          ggplot2::geom_point(
            data    = df_pred,
            mapping = ggplot2::aes(
              x      = .data[["x_plot"]],
              y      = .data[["y"]],
              colour = .data[["group"]]
            ),
            shape   = 18,
            size    = 4
          )
      } else {
        out <- out +
          ggplot2::geom_point(
            data    = df_pred,
            mapping = ggplot2::aes(x = .data[["x_plot"]], y = .data[["y"]]),
            shape   = 18,
            colour  = lcol[1],
            size    = 4
          )
      }
    }
  }

  if (mod_type == "categorical") {
    df_points[["x_plot"]] <- .regplot_dodge_x(df_points$x, df_points$group_id, n_groups, dodge_w)
    df_points[["x_plot"]] <- df_points[["x_plot"]] + .regplot_jitter_values(
      n      = nrow(df_points),
      amount = jitter_am / max(1, n_groups)
    )
  } else {
    df_points[["x_plot"]] <- df_points[["x"]]
  }

  if (has_by) {
    out <- out +
      ggplot2::geom_point(
        data    = df_points,
        mapping = ggplot2::aes(
          x      = .data[["x_plot"]],
          y      = .data[["y"]],
          size   = .data[["size"]],
          colour = .data[["group"]],
          fill   = .data[["group"]]
        ),
        shape   = pch
      )
  } else {
    out <- out +
      ggplot2::geom_point(
        data    = df_points,
        mapping = ggplot2::aes(
          x    = .data[["x_plot"]],
          y    = .data[["y"]],
          size = .data[["size"]]
        ),
        shape   = pch,
        colour  = col,
        fill    = bg
      )
  }

  out <- out +
    ggplot2::scale_size_identity()

  if (mod_type == "categorical") {
    out <- out +
      ggplot2::scale_x_continuous(
        breaks = seq_along(data$levels),
        labels = data$levels,
        limits = data$xlim,
        name   = data$xlab
      )
  } else {
    out <- out +
      ggplot2::scale_x_continuous(
        limits = data$xlim,
        name   = data$xlab
      )
  }

  out <- out +
    ggplot2::scale_y_continuous(
      limits = data$ylim,
      name   = data$ylab
    )

  if (has_by) {
    out <- out +
      ggplot2::scale_colour_manual(values = line_cols, name = data$by_name) +
      ggplot2::scale_fill_manual(values = line_cols, name = data$by_name)
  }

  if (!is.null(main)) {
    out <- out + ggplot2::ggtitle(main)
  }

  out <- out + ggplot2::guides(size = "none")
  if (has_by && !legend) {
    out <- out + ggplot2::guides(colour = "none", fill = "none")
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .set_dots_regplot
# ---------------------------------------------------------------------------- #
#
# Process dots arguments for regression plot with sensible defaults.
#
# @param ... additional graphical arguments
#
# @return list of graphical parameters with defaults applied
#
# ---------------------------------------------------------------------------- #
.set_dots_regplot <- function(...) {

  dots <- list(...)

  # point styling
  dots <- .plot_point_style_defaults(dots)

  # line styling
  if (is.null(dots[["lcol"]]))     dots[["lcol"]]     <- "black"
  if (is.null(dots[["lwd"]]))      dots[["lwd"]]      <- 2
  .check_plot_positive_scalar(dots[["lwd"]], "lwd")

  # band styling
  if (is.null(dots[["shade"]]))    dots[["shade"]]    <- TRUE
  if (is.null(dots[["col.ci"]]))   dots[["col.ci"]]   <- "gray70"
  if (is.null(dots[["col.pi"]]))   dots[["col.pi"]]   <- "gray85"
  if (is.null(dots[["col.si"]]))   dots[["col.si"]]   <- "gray92"
  if (is.null(dots[["alpha.ci"]])) dots[["alpha.ci"]] <- 0.4
  if (is.null(dots[["alpha.pi"]])) dots[["alpha.pi"]] <- 0.2
  if (is.null(dots[["alpha.si"]])) dots[["alpha.si"]] <- 0.15
  BayesTools::check_bool(dots[["shade"]], "shade")
  .check_plot_numeric(dots[["alpha.ci"]], "alpha.ci", check_length = 1, lower = 0, upper = 1, allow_null = FALSE)
  .check_plot_numeric(dots[["alpha.pi"]], "alpha.pi", check_length = 1, lower = 0, upper = 1, allow_null = FALSE)
  .check_plot_numeric(dots[["alpha.si"]], "alpha.si", check_length = 1, lower = 0, upper = 1, allow_null = FALSE)

  # categorical moderator jitter
  if (is.null(dots[["jitter"]]))      dots[["jitter"]]      <- 0.2
  if (is.null(dots[["box.width"]]))   dots[["box.width"]]   <- 0.5
  if (is.null(dots[["dodge.width"]])) dots[["dodge.width"]] <- 0.75
  .check_plot_numeric(dots[["jitter"]], "jitter", check_length = 1, lower = 0, allow_null = FALSE)
  .check_plot_positive_scalar(dots[["box.width"]], "box.width")
  .check_plot_numeric(dots[["dodge.width"]], "dodge.width", check_length = 1, lower = 0, allow_null = FALSE)

  # title (NULL = no title by default)
  if (is.null(dots[["main"]]))     dots[["main"]]     <- NULL
  if (is.null(dots[["las"]]))      dots[["las"]]      <- 1
  .check_plot_label(dots[["main"]], "main")
  BayesTools::check_int(dots[["las"]], "las", lower = 0)
  if (!dots[["las"]] %in% 0:3) {
    stop("'las' must be one of 0, 1, 2, or 3.", call. = FALSE)
  }
  .check_plot_positive_scalar(dots[["cex"]], "cex")
  .check_plot_positive_scalar(dots[["size"]], "size")

  return(dots)
}
