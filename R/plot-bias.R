#' @title Plots Weight Function of brma Object
#'
#' @description \code{plot_weightfunction.brma} visualizes the posterior
#' (and optionally prior) publication-bias weight function of a brma object.
#'
#' @param x a fitted \code{brma} object with a weightfunction/selection
#' component, such as \code{bselmodel()} or a \code{RoBMA()} object with
#' weightfunction priors.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param rescale_p_values whether to rescale p-values to the interval shown
#' by the weightfunction plot. Defaults to \code{TRUE}.
#' @param show_data whether observed one-sided p-values should be shown as rug
#' marks on the weightfunction axis. Defaults to \code{TRUE}.
#' @param dots_data list of additional graphical arguments for observed
#' p-value rug marks. Supported arguments include \code{col}/\code{color},
#' \code{alpha}, \code{lwd}/\code{linewidth}/\code{size},
#' \code{side}/\code{rug_side}, and
#' \code{height}/\code{rug_height}/\code{ticksize}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bselmodel(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   plot_weightfunction(fit)
#'   plot_weightfunction(fit, prior = TRUE)
#' }
#' }
#'
#'
#' @return \code{plot_weightfunction.brma} returns either \code{NULL} if
#' \code{plot_type = "base"} or a \code{ggplot2} object if
#' \code{plot_type = "ggplot"}.
#' The method errors for fitted objects without a weightfunction component.
#'
#' @seealso [RoBMA()]
#' @export
plot_weightfunction <- function(x, ...)  UseMethod("plot_weightfunction")

#' @export
#' @rdname plot_weightfunction
plot_weightfunction.brma  <- function(
    x, rescale_p_values = TRUE,
    prior = FALSE, plot_type = "base", show_data = TRUE,
    dots_data = NULL, dots_prior = NULL, ...) {

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(rescale_p_values, "rescale_p_values")
  BayesTools::check_bool(show_data, "show_data")
  BayesTools::check_list(dots_data, "dots_data", allow_NULL = TRUE)
  dots_raw <- list(...)
  if (is.null(dots_data)) {
    dots_data <- list()
  }

  if (!.is_weightfunction(x)) {
    stop("'plot_weightfunction' is available only for models with a weightfunction component.", call. = FALSE)
  }

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = .as_mixed_posteriors_parameters(x, "omega")
  )

  ### set up plotting arguments
  dots       <- do.call(.set_dots_plot, dots_raw)
  dots_prior <- .set_dots_prior(dots_prior)
  data_p     <- if (show_data) .weightfunction_observed_p_values(x) else NULL

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- "weightfunction"
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- FALSE
  args$rescale_x                <- rescale_p_values
  args$show_data                <- show_data
  args$data                     <- data_p
  args$dots_data                <- dots_data

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

.weightfunction_observed_p_values <- function(x) {

  outcome_data <- x[["data"]][["outcome"]]
  selection    <- .selection_spec(
    priors           = x[["priors"]],
    yi               = outcome_data[["yi"]],
    sei              = outcome_data[["sei"]],
    effect_direction = .effect_direction(x),
    signed_data      = FALSE
  )

  if (is.null(selection)) {
    return(NULL)
  }

  return(stats::pnorm(
    selection[["sign"]] * outcome_data[["yi"]] / outcome_data[["sei"]],
    lower.tail = FALSE
  ))
}


#' @title Plot PET-PEESE Fit of brma Object
#'
#' @description \code{plot_pet_peese} visualizes posterior
#' (and prior) PET-PEESE fit of a brma object.
#'
#' @param x a fitted brma object with a PET or PEESE component.
#' @param show_data whether the observed effect sizes should be added to
#' figure. Defaults to \code{TRUE}.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim},
#' \code{pch}, \code{bg}, \code{cex}, and \code{size}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, y-axis range, and observed data
#' point style respectively.
#'
#' @details The plot shows observed \code{yi} values against \code{sei}. PET
#' regression uses \eqn{\mu + PET \cdot se_i}; PEESE regression uses
#' \eqn{\mu + PEESE \cdot se_i^2}, with the fitted effect direction.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- bPET(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   plot_pet_peese(fit)
#' }
#' }
#'
#'
#' @return \code{plot_pet_peese} returns either \code{NULL} invisibly if
#' \code{plot_type = "base"} or a ggplot2 object if
#' \code{plot_type = "ggplot"}.
#'
#' @seealso [RoBMA()]
#' @export
plot_pet_peese <- function(x, ...)  UseMethod("plot_pet_peese")

#' @export
#' @rdname plot_pet_peese
plot_pet_peese.brma  <- function(
    x, show_data = TRUE,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(show_data, "show_data")

  if (!(.is_PET(x) || .is_PEESE(x))) {
    stop("'plot_pet_peese' is available only for models with a PET or PEESE component.", call. = FALSE)
  }

  ### obtain posterior samples in the plotting format
  location_parameter <- .check_and_select_plot_parameter(
    parameter       = "mu",
    parameter_mods  = NULL,
    parameter_scale = NULL,
    object          = x
  )
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = .as_mixed_posteriors_parameters(
      x,
      c(location_parameter, if (.is_PET(x)) "PET", if (.is_PEESE(x)) "PEESE")
    )
  )

  ### set up plotting arguments
  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  # set plotting range (make sure all effect sizes are within the plotting range)
  data_outcome <- x[["data"]][["outcome"]]
  if (is.null(dots[["ylim"]])) {
    dots[["ylim"]] <- range(pretty(c(0, range(data_outcome$yi))))
  }

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- "PETPEESE"
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- FALSE
  args$effect_direction         <- .effect_direction(x)

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # include data
  if (show_data) {
    data_dots <- .plot_point_style_defaults(dots)

    if (plot_type == "ggplot") {
      plot <- plot + ggplot2::geom_point(
        ggplot2::aes(
          x  = data_outcome$sei,
          y  = data_outcome$yi),
        size   = data_dots[["size"]],
        shape  = data_dots[["pch"]],
        colour = data_dots[["col"]],
        fill   = data_dots[["bg"]]
      )
    } else if (plot_type == "base") {
      graphics::points(
        data_outcome$sei, data_outcome$yi,
        cex = data_dots[["cex"]],
        pch = data_dots[["pch"]],
        col = data_dots[["col"]],
        bg  = data_dots[["bg"]]
      )
    }
  }

  # return the plots
  if (plot_type == "base") {
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}
