#' @title Plots a fitted RoBMA object
#'
#' @description \code{plot.RoBMA} allows to visualize
#' different \code{"RoBMA"} object parameters in various
#' ways. See \code{type} for the different model types.
#'
#' @param x a fitted RoBMA object
#' @param parameter a parameter to be plotted. Defaults to
#' \code{"mu"} (for the effect size). The additional options
#' are \code{"tau"} (for the heterogeneity),
#' \code{"weightfunction"} (for the estimated weightfunction),
#' or \code{"PET-PEESE"} (for the PET-PEESE regression).
#' @param conditional whether conditional estimates should be
#' plotted. Defaults to \code{FALSE} which plots the model-averaged
#' estimates. Note that both \code{"weightfunction"} and
#' \code{"PET-PEESE"} are always ignoring the other type of
#' publication bias adjustment.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param output_scale transform the effect sizes and the meta-analytic
#' effect size estimate to a different scale. Defaults to \code{NULL}
#' which returns the same scale as the model was estimated on.
#' @param rescale_x whether the x-axis of the \code{"weightfunction"}
#' should be re-scaled to make the x-ticks equally spaced.
#' Defaults to \code{FALSE}.
#' @param show_data whether the study estimates and standard
#' errors should be show in the \code{"PET-PEESE"} plot.
#' Defaults to \code{TRUE}.
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
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # the 'plot' function allows to visualize the results of a fitted RoBMA object, for example;
#' # the model-averaged effect size estimate
#' plot(fit, parameter = "mu")
#'
#' # and show both the prior and posterior distribution
#' plot(fit, parameter = "mu", prior = TRUE)
#'
#' # conditional plots can by obtained by specifying
#' plot(fit, parameter = "mu", conditional = TRUE)
#'
#' # plotting function also allows to visualize the weight function
#' plot(fit, parameter = "weightfunction")
#'
#' # re-scale the x-axis
#' plot(fit, parameter = "weightfunction", rescale_x = TRUE)
#'
#' # or visualize the PET-PEESE regression line
#' plot(fit, parameter = "PET-PEESE")
#' }
#'
#'
#' @return \code{plot.RoBMA} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot.RoBMA  <- function(x, parameter = "mu",
                        conditional = FALSE, plot_type = "base", prior = FALSE, output_scale = NULL,
                        rescale_x = FALSE, show_data = TRUE, dots_prior = NULL, ...){

  # forwarded options checked within the .extract_posterior.RoBMA function
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(rescale_x, "rescale_x")
  BayesTools::check_bool(show_data, "show_data")

  # TODO: add implementation of those in BayesTools
  # don't allow prior and posterior PET-PEESE for ss algorithm (the samples specification is not available)
  if(prior && parameter == "PETPEESE" && x$add_info[["algorithm"]] == "ss"){
    stop("The prior and posterior distribution for the PET-PEESE regression cannot be plotted for the ss algorithm.")
  }
  # don't allow conditional PET-PEESE for ss algorithm (the samples specification is not available)
  if(conditional && parameter == "PETPEESE" && x$add_info[["algorithm"]] == "ss"){
    stop("The conditional distribution for the PET-PEESE regression cannot be plotted for the ss algorithm.")
  }

  # get the name of the parameter
  parameter         <- .check_and_transform_parameter_name(parameter, x$add_info[["predictors"]])
  parameter_samples <- .check_and_transform_parameter_samples_name(parameter, x$add_info[["predictors"]])

  ### obtain posterior samples
  samples <- .extract_posterior.RoBMA(x, parameter, conditional)

  ### manage transformations
  # (transformation passed to the plotting function so prior density can be also transformed)
  # get the settings
  results_scale <- x$add_info[["output_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale, estimation = FALSE)
  }
  # set the transformations
  if(parameter != "omega" && results_scale != output_scale){
    if(parameter == "PEESE"){
      # the transformation is inverse for PEESE
      transformation <- .get_transform_mu(output_scale, results_scale, fun = FALSE)
    }else if(parameter == "PET"){
      # PET is scale invariant
      transformation <- NULL
    }else if(parameter == "rho"){
      # rho is scale invariant
      transformation <- NULL
    }else if(parameter %in% c("mu", "PETPEESE") || parameter %in% x$add_info[["predictors"]]){
      transformation <- .get_transform_mu(results_scale, output_scale, fun = FALSE)
    }else if(parameter == "tau"){
      transformation <- .get_scale(results_scale, output_scale, fun = FALSE)
    }
  }else{
    transformation <- NULL
  }

  # remove PET/PEESE prior class (otherwise PET-PEESE density is produced in BayesTools)
  if(parameter %in% c("PET", "PEESE")){
    for(i in seq_along(attr(samples[["PET"]], "prior_list"))){
      if(is.prior.PET(attr(samples[["PET"]], "prior_list")[[i]])){
        class(attr(samples[["PET"]], "prior_list")[[i]])   <- class(attr(samples[["PET"]], "prior_list")[[i]])[!class(attr(samples[["PET"]], "prior_list")[[i]]) %in% "prior.PET"]
      }else if(!is.prior.point(attr(samples[["PET"]], "prior_list")[[i]])){
        attr(samples[["PET"]], "prior_list")[[i]] <- prior("point", parameter = list("location" = 0), prior_weights = attr(samples[["PET"]], "prior_list")[[i]][["prior_weights"]])
      }
    }
    for(i in seq_along(attr(samples[["PEESE"]], "prior_list"))){
      if(is.prior.PEESE(attr(samples[["PEESE"]], "prior_list")[[i]])){
        class(attr(samples[["PEESE"]], "prior_list")[[i]])   <- class(attr(samples[["PEESE"]], "prior_list")[[i]])[!class(attr(samples[["PEESE"]], "prior_list")[[i]]) %in% "prior.PEESE"]
      }else if(!is.prior.point(attr(samples[["PEESE"]], "prior_list")[[i]])){
        attr(samples[["PEESE"]], "prior_list")[[i]] <- prior("point", parameter = list("location" = 0), prior_weights = attr(samples[["PEESE"]], "prior_list")[[i]][["prior_weights"]])
      }
    }
  }


  if(parameter %in% x$add_info[["predictors"]]){
    if(inherits(samples[[parameter_samples]], "mixed_posteriors.factor")){
      if(attr(samples[[parameter_samples]],"orthonormal") || attr(samples[[parameter_samples]],"meandif")){
        n_levels <- length(attr(samples[[parameter_samples]],"level_names"))
      }else if(attr(samples[[parameter_samples]],"treatment")){
        n_levels <- length(attr(x$add_info[["predictors"]],"level_names")) - 1
      }
    }else{
      n_levels <- 1
    }
  }else{
    n_levels <- 1
  }


  dots       <- .set_dots_plot(..., n_levels = n_levels)
  dots_prior <- .set_dots_prior(dots_prior)

  if(parameter == "PETPEESE" & show_data){
    data <- x[["data"]]
    data <- combine_data(data = data, transformation = .transformation_invar(output_scale))
    if(is.null(dots[["xlim"]])){
      dots[["xlim"]] <- range(pretty(c(0, data$se)))
    }
  }

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- parameter_samples
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$transformation           <- transformation
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$rescale_x                <- rescale_x
  args$par_name                 <- if(parameter %in% c("mu", "tau", "PET", "PEESE", x$add_info[["predictors"]])) .plot.RoBMA_par_names(parameter, x, output_scale)[[1]]
  args$dots_prior               <- dots_prior

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))


  if(parameter == "PETPEESE" & show_data){

    if(plot_type == "ggplot"){
      plot <- plot + ggplot2::geom_point(
          ggplot2::aes(
            x  = data$se,
            y  = data$y),
          size  = if(!is.null(dots[["cex"]])) dots[["cex"]] else 2,
          shape = if(!is.null(dots[["pch"]])) dots[["pch"]] else 18
        )
      # make sure all effect sizes are within the plotting range
      if(is.null(dots[["ylim"]]) & any(data$y < plot$plot_env$ylim[1] | data$y > plot$plot_env$ylim[2])){
        new_y_range <- range(c(data$y, plot$plot_env$ylim))
        plot <- suppressMessages(plot + ggplot2::scale_y_continuous(breaks = pretty(new_y_range), limits = range(pretty(new_y_range)), oob = scales::oob_keep))
      }

    }else if(plot_type == "base"){
      graphics::points(data$se, data$y, cex = if(!is.null(dots[["cex"]])) dots[["cex"]] else 2, pch = if(!is.null(dots[["pch"]])) dots[["pch"]] else 18)
    }
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

#' @title Forest plot for a RoBMA object
#'
#' @description \code{forest} creates a forest plot for
#' a \code{"RoBMA"} object.
#'
#' @param order order of the studies. Defaults to \code{NULL} -
#' ordering as supplied to the fitting function. Studies
#' can be ordered either \code{"increasing"} or \code{"decreasing"} by
#' effect size, or by labels \code{"alphabetical"}.
#' @inheritParams plot.RoBMA
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # the forest function creates a forest plot for a fitted RoBMA object, for example,
#' # the forest plot for the individual studies and the model-averaged effect size estimate
#' forest(fit)
#'
#' # the conditional effect size estimate
#' forest(fit, conditional = TRUE)
#'
#' # or transforming the effect size estimates to Fisher's z
#' forest(fit, output_scale = "fishers_z")
#' }
#'
#' @return \code{forest} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @export
forest <- function(x, conditional = FALSE, plot_type = "base", output_scale = NULL, order = NULL, ...){

  # apply version changes to RoBMA object
  x <- .update_object(x)

  # check whether plotting is possible
  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  #check settings
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing", "alphabetical"))


  ### get the posterior samples & data
  if(conditional){
    samples_mu <- x[["RoBMA"]][["posteriors_conditional"]][["mu"]]
    if(is.null(samples_mu))
      stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect.")
  }else{
    samples_mu <- x[["RoBMA"]][["posteriors"]][["mu"]]
  }

  if(is.RoBMA.reg(x)){
    data <- x[["data"]][["outcome"]]
  }else{
    data <- x[["data"]]
  }


  ### manage transformations
  # get the settings
  results_scale <- x$add_info[["output_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale, estimation = FALSE)
  }

  # set the transformations
  if(results_scale != output_scale){
    samples_mu <- .transform_mu(samples_mu, results_scale, output_scale)
  }

  # obtain the posterior estimates
  est_mu <- mean(samples_mu)
  lCI_mu <- unname(stats::quantile(samples_mu, .025))
  uCI_mu <- unname(stats::quantile(samples_mu, .975))

  # get the CIs (+add transformation if necessary)
  if(is.BiBMA(x)){
    data <- .combine_data_bi(data = data, transformation = .transformation_invar(output_scale, estimation = FALSE), return_all = TRUE, estimation = FALSE)
  }else{
    data <- combine_data(data = data, transformation = .transformation_invar(output_scale, estimation = FALSE), return_all = TRUE, estimation = FALSE)
  }

  # simplify the data structure
  data$y <- data[,output_scale]
  data   <- data[,c("y", "lCI", "uCI", "study_names")]


  # add ordering
  if(!is.null(order)){
    if(order == "increasing"){
      data <- data[order(data$y, decreasing = FALSE),]
    }else if(order == "decreasing"){
      data <- data[order(data$y, decreasing = TRUE),]
    }else if(order == "alphabetical"){
      data <- data[order(data$study_names),]
    }
  }

  data$x <- (nrow(data)+2):3

  # additional plot settings
  y_at      <- c(1, data$x)
  y_labels  <- c(ifelse(conditional, "Conditional", "Model-Averaged"), data$study_names)
  y_labels2 <- paste0(
    format(round(c(est_mu, data$y), 2), nsmall = 2),
    " [", format(round(c(lCI_mu, data$lCI), 2), nsmall = 2), ", ",
    format(round(c(uCI_mu, data$uCI), 2), nsmall = 2), "]")
  ylim      <- c(0, max(data$x) + 1)
  xlab      <- .plot.RoBMA_par_names("mu", x, output_scale)[[1]]
  x_labels  <- pretty(range(c(data$lCI, data$uCI, lCI_mu, uCI_mu)))
  xlim      <- range(pretty(range(c(data$lCI, data$uCI, lCI_mu, uCI_mu))))


  ### do the plotting
  dots <- .set_dots_plot(...)

  if(plot_type == "base"){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4, max(nchar(y_labels)) * 1/2, 0, max(nchar(y_labels2)) * 1/2))
    }else{
      graphics::par(...)
    }

    graphics::plot(NA, bty = "n", las = 1, xlab = xlab, ylab = "", main = "", yaxt = "n", ylim = ylim, xlim = xlim)
    graphics::axis(2, at = y_at, labels = y_labels,  las = 1, col = NA)
    graphics::axis(4, at = y_at, labels = y_labels2, las = 1, col = NA, hadj = 0)
    if(output_scale == "y"){
      graphics::abline(v = 0, lty = 3)
    }else{
      graphics::abline(v = .transform_mu(0, "d", output_scale), lty = 3)
    }

    graphics::arrows(x0 = data$lCI, x1 = data$uCI, y0 = data$x, code = 3, angle = 90, length = 0.1)
    graphics::points(x = data$y, y = data$x, pch = 15)

    graphics::polygon(
      x = c(lCI_mu, est_mu , uCI_mu, est_mu),
      y = c(1, 1.25, 1, 0.75),
      col = "black"
    )

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()

    # add the studies
    plot <- plot +ggplot2::geom_errorbarh(
      mapping = ggplot2::aes(
        xmin   = data$lCI,
        xmax   = data$uCI,
        y      = data$x),
      color   = "black",
      height  = .25)
    plot <- plot +ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = data$y,
        y = data$x),
      shape = 15)

    # add the overall estimate
    plot <- plot + ggplot2::geom_polygon(
      mapping = ggplot2::aes(
        x = c(lCI_mu, est_mu , uCI_mu, est_mu),
        y = c(1, 1.25, 1, 0.75)),
      fill = "black")

    # add the vertical line
    if(output_scale == "y"){
      plot <- plot + ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = c(0,0),
          y = ylim),
        linetype = "dotted")
    }else{
      plot <- plot + ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = .transform_mu(c(0,0), "d", output_scale),
          y = ylim),
        linetype = "dotted")
    }

    # add all the other stuff
    plot <- plot + ggplot2::scale_y_continuous(
      name = "", breaks = y_at, labels = y_labels, limits = ylim,
      sec.axis = ggplot2::sec_axis( ~ ., breaks = y_at, labels = y_labels2))
    plot <- plot + ggplot2::scale_x_continuous(
      name = xlab, breaks = x_labels, labels = x_labels, limits = xlim)
    plot <- plot + ggplot2::theme(
      axis.title.y      = ggplot2::element_blank(),
      axis.line.y       = ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_text(hjust = 0, color = "black"),
      axis.text.y.right = ggplot2::element_text(hjust = 1, color = "black"))

    attr(plot, "sec_axis") <- TRUE
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


#' @title Models plot for a RoBMA object
#'
#' @description \code{plot_models} plots individual models'
#' estimates for a \code{"RoBMA"} object.
#'
#' @param parameter a parameter to be plotted. Defaults to
#' \code{"mu"} (for the effect size). The additional option
#' is \code{"tau"} (for the heterogeneity).
#' @param order how the models should be ordered.
#' Defaults to \code{"decreasing"} which orders them in decreasing
#' order in accordance to \code{order_by} argument. The alternative is
#' \code{"increasing"}.
#' @param order_by what feature should be use to order the models.
#' Defaults to \code{"model"} which orders the models according to
#' their number. The alternatives are \code{"estimate"} (for the effect
#' size estimates), \code{"probability"} (for the posterior model probability),
#' and \code{"BF"} (for the inclusion Bayes factor).
#' @inheritParams plot.RoBMA
#' @inheritParams forest
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # the plot_models function creates a plot for of the individual models' estimates, for example,
#' # the effect size estimates from the individual models can be obtained with
#' plot_models(fit)
#'
#' # and effect size estimates from only the conditional models
#' plot_models(fit, conditional = TRUE)
#' }
#'
#'
#' @return \code{plot_models} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @export
plot_models <- function(x, parameter = "mu", conditional = FALSE, output_scale = NULL, plot_type = "base", order = "decreasing", order_by = "model", ...){

  # apply version changes to RoBMA object
  x <- .update_object(x)

  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  if(x$add_info[["algorithm"]] == "ss")
    stop("The estimated model using the spike and slab style model-averaging (`algorithm = 'ss'`). As such, no models are directly iterated over and the individual model estimates cannot be displayed.")

  #check settings
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing"))
  if(!parameter %in% c("mu", "tau")){
    stop("The passed parameter is not supported for plotting. See '?plot.RoBMA' for more details.")
  }
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing"))
  BayesTools::check_char(order_by, "order_by", allow_NULL = TRUE, allow_values = c("model", "estimate", "probability", "BF"))


  ### manage transformations
  # get the settings
  results_scale <- x$add_info[["output_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale, estimation = FALSE)
  }
  # set the transformations
  if(results_scale != output_scale){
    if(parameter == "PETPEESE"){
      # the transformation is inverse for PEESE
      transformation <- .get_scale(output_scale, results_scale, fun = FALSE)
    }else if(parameter == "PET"){
      # PET is scale invariant
      transformation <- NULL
    }else if(parameter == "mu"){
      transformation <- .get_transform_mu(results_scale, output_scale, fun = FALSE)
    }else if(parameter == "tau"){
      transformation <- .get_scale(results_scale, output_scale, fun = FALSE)
    }
  }else{
    transformation <- NULL
  }


  ### prepare input
  if(is.RoBMA.reg(x) && parameter == "mu"){

    if(conditional){

      model_list <- x[["models"]]
      samples    <- x[["RoBMA"]][["posteriors_predictors_conditional"]]
      inference  <- x[["RoBMA"]][["inference_conditional"]]

      # check whether the input exists
      if(parameter == "mu_intercept" && is.null(samples[["mu_intercept"]]))
        stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect.")

    }else{

      model_list <- x[["models"]]
      samples    <- x[["RoBMA"]][["posteriors_predictors"]]
      inference  <- x[["RoBMA"]][["inference"]]

    }

    parameter <- "mu_intercept"
    names(samples)[names(samples) == "mu_intercept"]     <- "mu_intercept"
    names(inference)[names(inference) == "Effect"]       <- "mu_intercept"

  }else{

    if(conditional){

      model_list <- x[["models"]]
      samples    <- x[["RoBMA"]][["posteriors_conditional"]]
      inference  <- x[["RoBMA"]][["inference_conditional"]]

      # check whether the input exists
      if(parameter == "mu" && is.null(samples[["mu"]]))
        stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect.")
      if(parameter == "tau" && is.null(samples[["tau"]]))
        stop("The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the heterogeneity. Please, verify that you specified at least one model assuming the presence of the heterogeneity.")

    }else{

      model_list <- x[["models"]]
      samples    <- x[["RoBMA"]][["posteriors"]]
      inference  <- x[["RoBMA"]][["inference"]]

    }

    # deal with the non-matching names
    names(inference)[names(inference) == "Effect"]        <- "mu"
    names(inference)[names(inference) == "Heterogeneity"] <- "tau"
  }


  dots <- list(...)

  # prepare the argument call
  args                          <- dots
  args$model_list               <- model_list
  args$samples                  <- samples
  args$inference                <- inference
  args$parameter                <- parameter
  args$plot_type                <- plot_type
  args$prior                    <- FALSE
  args$conditional              <- conditional
  args$order                    <- list(order, order_by)
  args$transformation           <- transformation
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$par_name                 <- .plot.RoBMA_par_names(parameter, x, output_scale)[[1]]

  plot <- do.call(BayesTools::plot_models, args)


  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}



.set_dots_plot        <- function(..., n_levels = 1){

  dots <- list(...)
  if(is.null(dots[["col"]]) & n_levels == 1){
    dots[["col"]]      <- "black"
  }else if(is.null(dots[["col"]]) & n_levels > 1){
    dots[["col"]]      <- grDevices::palette.colors(n = n_levels + 1, palette = "Okabe-Ito")[-1]
  }
  if(is.null(dots[["col.fill"]])){
    dots[["col.fill"]] <- "#4D4D4D4C" # scales::alpha("grey30", .30)
  }

  return(dots)
}
.set_dots_prior       <- function(dots_prior){

  if(is.null(dots_prior)){
    dots_prior <- list()
  }

  if(is.null(dots_prior[["col"]])){
    dots_prior[["col"]]      <- "grey60"
  }
  if(is.null(dots_prior[["lty"]])){
    dots_prior[["lty"]]      <- 1
  }
  if(is.null(dots_prior[["col.fill"]])){
    dots_prior[["col.fill"]] <- "#B3B3B34C" # scales::alpha("grey70", .30)
  }

  return(dots_prior)
}
.plot.RoBMA_par_names <- function(par, fit, output_scale){

  if(par %in% c("mu", "mu_intercept")){

    par_names <- list(switch(
      output_scale,
      "r"     = expression(rho),
      "d"     = expression("Cohen's"~italic(d)),
      "z"     = expression("Fisher's"~italic(z)),
      "logOR" = expression("logOR"),
      "OR"    = expression("OR"),
      "y"     = expression(mu)
    ))

  }else if(par == "tau"){

    par_names <- list(switch(
      output_scale,
      "r"     = expression(tau~(rho)),
      "d"     = expression(tau~("Cohen's"~italic(d))),
      "z"     = expression(tau~("Fisher's"~italic(z))),
      "logOR" = expression(tau~("logOR")),
      "OR"    = expression(tau~("logOR")),
      "y"     = expression(tau)
    ))

  }else if(par == "omega"){

    stop("should not be used")
    # summary_info <- summary(fit)
    # sum_all      <- summary_info[["estimates"]]
    # omega_names  <- rownames(sum_all)[grepl(par, rownames(sum_all))]
    # par_names    <- vector("list", length = length(omega_names))
    # for(i in 1:length(par_names)){
    #   par_names[[i]] <- bquote(~omega[~.(substr(omega_names[i],6,nchar(omega_names[i])))])
    # }

  }else if(par == "theta"){

    stop("should not be used")
    # add_type <- if(!is.null(type)) paste0("(",type,")") else NULL
    # par_names <- vector("list", length = length(study_names))
    # for(i in 1:length(par_names)){
    #   par_names[i] <- list(switch(
    #     output_scale,
    #     "r"     = bquote(~rho[.(paste0("[",study_names[i],"]"))]),
    #     "d"     = bquote("Cohen's"~italic(d)[.(paste0("[",study_names[i],"]"))]),
    #     "z"     = bquote("Fisher's"~italic(z)[.(paste0("[",study_names[i],"]"))]),
    #     "logOR" = bquote("log"(italic("OR"))[.(paste0("[",study_names[i],"]"))]),
    #     "OR"    = bquote(~italic("OR")[.(paste0("[",study_names[i],"]"))]),
    #     "y"     = bquote(~mu[.(paste0("[",study_names[i],"]"))])
    #   ))
    # }

  }else if(par == "PET"){

    par_names <- list(switch(
      output_scale,
      "r"     = expression("PET"~(rho)),
      "d"     = expression("PET"~("Cohen's"~italic(d))),
      "z"     = expression("PET"~("Fisher's"~italic(z))),
      "logOR" = expression("PET"~("logOR")),
      "OR"    = expression("PET"~("OR")),
      "y"     = expression("PET")
    ))

    par_names <- list(bquote("PET"))

  }else if(par == "PEESE"){

    par_names <- list(switch(
      output_scale,
      "r"     = expression("PEESE"~(rho)),
      "d"     = expression("PEESE"~("Cohen's"~italic(d))),
      "z"     = expression("PEESE"~("Fisher's"~italic(z))),
      "logOR" = expression("PEESE"~("logOR")),
      "OR"    = expression("PEESE"~("OR")),
      "y"     = expression("PEESE")
    ))

  }else{

    par_names <- list(switch(
      output_scale,
      "r"     = bquote(.(par)~(rho)),
      "d"     = bquote(.(par)~("Cohen's"~italic(d))),
      "z"     = bquote(.(par)~("Fisher's"~italic(z))),
      "logOR" = bquote(.(par)~("logOR")),
      "OR"    = bquote(.(par)~("OR")),
      "y"     = bquote(.(par))
    ))
  }

  return(par_names)
}
