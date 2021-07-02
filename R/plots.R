#' @title Plots a fitted RoBMA object
#'
#' @description \code{plot.RoBMA} allows to visualize
#' different \code{"RoBMA"} object parameters in various
#' ways. See \code{type} for the different model types.
#'
#' @param x a fitted RoBMA object
#' @param parameter a parameter to be plotted. Either
#' \code{"mu"}, \code{"tau"}, \code{"theta"}, or
#' \code{"omega"}. A bivariate plot for model-averaged
#' estimates of "mu" and "tau" can be obtained by
#' \code{c("mu","tau")} if \code{type = "averaged"}. In
#' addition, a forest plot with the original estimates can
#' be obtained by \code{"forest"} or added to the theta
#' estimates by \code{c("theta", "forest")}.
#' @param type what type of estimates should be plotted.
#' Options are \code{"averaged"} for the model-averaged
#' estimates, \code{"conditional"} for the conditional
#' estimates, or \code{"individual"} for the
#' individual models estimates. The options
#' \code{c("individual", "conditional")} can be supplied
#' together to show only coditional individual models.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot2"} for plotting. The later
#' requires \pkg{ggplot2} package to be installed.
#' @param mean whether the mean should be plotted.
#' @param median whether the median should be plotted.
#' @param CI width of the confidence intervals.
#' @param prior add prior density to the plot. Only available
#' for \code{type = "averaged"} or \code{type = "conditional"}.
#' Defaults to \code{FALSE}.
#' @param digits_estimates number of decimals to be displayed for
#' \code{parameter = "theta"}, \code{parameter = "forest"}, and
#' \code{type = "individual"} plot.
#' @param show_figures which figures should be returned in the case
#' when multiple plots are generated. Useful when
#' \code{parameter = "omega", type = "individual"} which generates
#' a figure for each weights cut-off. Defaults to \code{-1} which
#' omits the first weight. Set to \code{NULL} to show all figures
#' or to \code{c(1,3)} to show only the first and third one.
#' @param weights whether the weights or weight function should
#' be returned. Only applicable when \code{parameter = "omega"}.
#' Defaults to \code{FALSE} - the weight function is plotted.
#' @param rescale_x whether the x-axis should be rescaled in order
#' to make the x-ticks equally spaced. Available only for the
#' weightfunction plot. Defaults to \code{FALSE}.
#' @param ... additional arguments to be passed to
#' \link[graphics]{par} if \code{plot_type = "base"}. Especially
#' useful for \code{parameter == "theta"},
#' \code{parameter == "forest"} or \code{type = "individual"}
#' where automatic margins might cut out parts of the labels.
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # plot function allows to visualize the results of a fitted RoBMA object, for example,
#' # the model-averaged mean parameter estimate
#' plot(fit, parameter = "mu")
#'
#' # or show both the prior and posterior distribution
#' plot(fit, parameter = "mu", prior = TRUE)
#'
#' # condtional plots might by obtained by specifying
#' plot(fit, parameter = "mu", type = "conditional")
#'
#' # plotting function also allows to visualize the weight function
#' # (or individual weights by adding 'weights = TRUE')
#' plot(fit, parameter = "omega")
#'
#' # or the forest plot (the estimated study effects can be shown by setting 'parameter = "theta"')
#' plot(fit, parameter = "forest")
#'
#' # it is also possible to compare the individual model estimates
#' # and order them by the posterior probability
#' plot(fit, parameter = "mu", type = "individual", order = "prob")
#'
#' }
#'
#' @seealso [RoBMA()]
#' @export
plot.RoBMA  <- function(x, parameter,
                        conditional = FALSE, plot_type = "base", output_scale = NULL, prior = FALSE,
                        rescale_x = FALSE, show_data = TRUE, dots_prior = NULL, ...){

  # check whether plotting is possible
  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  #check settings
  BayesTools::check_char(parameter, "parameter")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  if(plot_type == "ggplot" && !try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(rescale_x, "rescale_x")
  BayesTools::check_bool(show_data, "show_data")

  # deal with bad parameter names for PET-PEESE, weightfunction
  if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
    parameter <- "omega"
  }else if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) == "petpeese"){
    parameter <- "PETPEESE"
  }else if(!parameter %in% c("mu", "tau")){
    stop("The passed parameter is not supported for plotting. See '?plot.RoBMA' for more details.")
  }


  ### manage transformations
  # get the settings
  prior_scale <- x$add_info[["prior_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale)
  }
  # set the transformations
  if(prior_scale != output_scale){
    # TODO: figure out transformations, including jacobinas for effect sizes transformations - needed for the priors
    stop("Plotting output on a different than prior scale is not possible yet.")
    if(parameter == "PETPEESE"){

      # the transformation is inverse for PEESE
      transformation           = "lin"
      transformation_arguments = list(a = 0, b = .get_scale_b(output_scale, priors_scale))

    }else if(parameter == "mu"){

      # this transformation needs jacobian!

    }else if(parameter == "tau"){

      transformation           = "lin"
      transformation_arguments = list(a = 0, b = .get_scale_b(priors_scale, output_scale))

    }
  }else{
    transformation           <- NULL
    transformation_arguments <- NULL
  }


  # choose the samples
  if(conditional && parameter == "PETPEESE"){
    # get model-averaged posterior across PET and PEESE parameters
    models <- x[["models"]]

    effect <- sapply(x[["models"]], function(model)!.is_component_null(model[["priors"]], "effect"))
    PET    <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
    PEESE  <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

    if(any(PET) & any(PEESE)){
      is_conditional  <- effect & (PET | PEESE)
      if(sum(is_conditional) == 0)
        stop("There is no model that assumes presence of the effect and PET-PEESE bias correction")
      parameters      <- c("mu", "PET", "PEESE")
      parameters_null <- c("mu" = list(!is_conditional), "PET" = list(!is_conditional), "PEESE" = list(!is_conditional))
    }else if(any(PET)){
      is_conditional  <- effect & PET
      if(sum(is_conditional) == 0)
        stop("There is no model that assumes presence of the effect and PET-PEESE bias correction")
      parameters      <- c("mu", "PET")
      parameters_null <- c("mu" = list(!is_conditional), "PET"   = list(!is_conditional))
    }else if(any(PEESE)){
      is_conditional  <- effect & PEESE
      if(sum(is_conditional) == 0)
        stop("There is no model that assumes presence of the effect and PET-PEESE bias correction")
      parameters      <- c("mu", "PEESE")
      parameters_null <- c("mu" = list(!is_conditional), "PEESE" = list(!is_conditional))
    }

    samples <- BayesTools::mix_posteriors(
      model_list   = models,
      parameters   = parameters,
      is_null_list = parameters_null,
      seed         = x$add_info[["seed"]],
      conditional  = TRUE
    )

  }else{
    if(conditional){
      samples <- x[["RoBMA"]][["posteriors_conditional"]]
    }else{
      samples <- x[["RoBMA"]][["posteriors"]]
    }
  }

  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  plot <- do.call(BayesTools::plot_posterior, c(
    samples                  = list(samples),
    parameter                = parameter,
    plot_type                = plot_type,
    prior                    = prior,
    n_points                 = 1000,
    n_samples                = 10000,
    force_samples            = FALSE,
    transformation           = transformation,
    transformation_arguments = transformation_arguments,
    transformation_settings  = FALSE,
    rescale_x                = rescale_x,
    par_name                 = NULL,
    dots_prior               = list(dots_prior),
    dots))


  if(parameter == "PETPEESE" & show_data){

    data <- x[["data"]]
    # TODO: change from 'prior_scale' to 'output_scale' when plotting scale can be changed
    data <- combine_data(data = data, transformation = .transformation_invar(prior_scale))

    if(plot_type == "ggplot"){
      plot <- plot + ggplot2::geom_point(
          ggplot2::aes(
            x  = data$se,
            y  = data$y),
          size  = if(!is.null(dots[["cex"]])) dots[["cex"]] else 2,
          shape = if(!is.null(dots[["pch"]])) dots[["pch"]] else 18
        )
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
#' @param order order of the studies. Defaults to \code{NULL} -
#' ordering as supplied to the fitting function. Studies
#' can be ordered either \code{"ascending"} or \code{"descending"} by
#' effect size, or by labels \code{"alphabetical"}.
#' @inheritParams plot.RoBMA
#'
#' @export
forest <- function(x, conditional = FALSE, plot_type = "base", output_scale = NULL, prior = FALSE, order = NULL, ...){

  # check whether plotting is possible
  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  #check settings
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  if(plot_type == "ggplot" && !try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing", "alphabetical"))


  ### get the posterior samples & data
  if(conditional){
    samples_mu <- x[["RoBMA"]][["posteriors_conditional"]][["mu"]]
  }else{
    samples_mu <- x[["RoBMA"]][["posteriors"]][["mu"]]
  }
  data <- x[["data"]]


  ### manage transformations
  # get the settings
  prior_scale <- x$add_info[["prior_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale)
  }

  # set the transformations
  if(prior_scale != output_scale){
    samples_mu <- .transform_mu(samples_mu, prior_scale, output_scale)
  }

  # obtain the posterior estimates
  est_mu <- mean(samples_mu)
  lCI_mu <- unname(quantile(samples_mu, .025))
  uCI_mu <- unname(quantile(samples_mu, .975))

  # get the CIs (+add transformation if necessary)
  data <- combine_data(data = data, transformation = .transformation_invar(output_scale), return_all = TRUE)

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
      graphics::par(mar = c(4, max(nchar(y_labels)) * 2/3, 0, max(nchar(y_labels2)) * 1/2))
    }else{
      graphics::par(...)
    }

    graphics::plot(NA, bty = "n", las = 1, xlab = xlab, ylab = "", main = "", yaxt = "n", ylim = ylim, xlim = xlim)
    graphics::axis(2, at = y_at, labels = y_labels,  las = 1, col = NA)
    graphics::axis(4, at = y_at, labels = y_labels2, las = 1, col = NA, hadj = 0)
    graphics::abline(v = .transform_mu(0, "d", output_scale), lty = 3)

    arrows(x0 = data$lCI, x1 = data$uCI, y0 = data$x, code = 3, angle = 90, length = 0.1)
    graphics::points(x = data$y, y = data$x, pch = 15)

    graphics::polygon(
      x = c(lCI_mu, est_mu , uCI_mu, est_mu),
      y = c(1, 1.25, 1, 0.75),
      col = "black"
    )

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot(data = data)

    # add the studies
    plot <- plot +ggplot2::geom_errorbarh(
      mapping = ggplot2::aes(
        xmin   = lCI,
        xmax   = uCI,
        y      = x),
      color   = "black",
      height  = .25)
    plot <- plot +ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = y,
        y = x),
      shape = 15)

    # add the overall estimate
    plot <- plot + ggplot2::geom_polygon(
      mapping = ggplot2::aes(
        x = x,
        y = y),
      data = data.frame(
        x = c(lCI_mu, est_mu , uCI_mu, est_mu),
        y = c(1, 1.25, 1, 0.75)),
      fill = "black")

    # add the vertical line
    plot <- plot + ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = x,
        y = y),
      data     = data.frame(
        x = .transform_mu(c(0,0), "d", output_scale),
        y = ylim
      ),
      linetype = "dotted")

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
#' @param order either of the studies. Defaults to \code{NULL} -
#' ordering as supplied to the fitting function. Studies
#' can be ordered either \code{"ascending"} or \code{"descending"} by
#' effect size, or by labels \code{"alphabetical"}.
#' @inheritParams plot.RoBMA
#'
#' @export
plot_models <- function(x, parameter, conditional = FALSE, plot_type = "base", output_scale = NULL, prior = FALSE, order = "decreasing", order_by = "model", ...){

  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  #check settings
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  if(plot_type == "ggplot" && !try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing", "alphabetical"))
  if(!parameter %in% c("mu", "tau")){
    stop("The passed parameter is not supported for plotting. See '?plot.RoBMA' for more details.")
  }
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing", "alphabetical"))
  BayesTools::check_char(order_by, "order_by", allow_NULL = TRUE, allow_values = c("model", "estimate", "probability", "BF"))


  ### manage transformations
  # get the settings
  prior_scale <- x$add_info[["prior_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale)
  }
  # set the transformations
  if(prior_scale != output_scale){
    # TODO: figure out transformations, including jacobinas for effect sizes transformations - needed for the priors
    stop("Plotting output on a different than prior scale is not possible yet.")
    if(parameter == "PEESE"){

      # the transformation is inverse for PEESE
      transformation           = "lin"
      transformation_arguments = list(a = 0, b = .get_scale_b(output_scale, priors_scale))

    }else if(parameter == "mu"){

      # this transformation needs jacobian!

    }else if(parameter == "tau"){

      transformation           = "lin"
      transformation_arguments = list(a = 0, b = .get_scale_b(priors_scale, output_scale))

    }
  }else{
    transformation           <- NULL
    transformation_arguments <- NULL
  }


  ### prepare input
  if(conditional){

    model_list <- x[["models"]]
    samples    <- x[["RoBMA"]][["posteriors_conditional"]]
    inference  <- x[["RoBMA"]][["inference_conditional"]]

  }else{

    model_list <- x[["models"]]
    samples    <- x[["RoBMA"]][["posteriors"]]
    inference  <- x[["RoBMA"]][["inference"]]

  }

  # deal with the non-matching names
  names(inference)[names(inference) == "Effect"]        <- "mu"
  names(inference)[names(inference) == "Heterogeneity"] <- "tau"

  dots <- list(...)

  plot <- do.call(BayesTools::plot_models, c(
    model_list               = list(model_list),
    samples                  = list(samples),
    inference                = list(inference),
    parameter                = parameter,
    plot_type                = plot_type,
    prior                    = prior,
    condtional               = conditional,
    order                    = list(list(order, order_by)),
    transformation           = transformation,
    transformation_arguments = transformation_arguments,
    transformation_settings  = FALSE,
    par_name                 = .plot.RoBMA_par_names(parameter, x, output_scale)[[1]],
    dots))


  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}



.set_dots_plot        <- function(...){

  dots <- list(...)
  if(is.null(dots[["col"]])){
    dots[["col"]]      <- "black"
  }
  if(is.null(dots[["col.fill"]])){
    dots[["col.fill"]] <- "#CCCCCC33" # scales::alpha("grey80", .20)
  }

  return(dots)
}
.set_dots_prior       <- function(dots_prior){

  if(is.null(dots_prior)){
    dots_prior <- list()
  }

  if(is.null(dots_prior[["col"]])){
    dots_prior[["col"]]      <- "grey30"
  }
  if(is.null(dots_prior[["lty"]])){
    dots_prior[["lty"]]      <- 2
  }
  if(is.null(dots_prior[["col.fill"]])){
    dots_prior[["col.fill"]] <- "#4D4D4D33" # scales::alpha("grey30", .20)
  }

  return(dots_prior)
}
.plot.RoBMA_par_names <- function(par, fit, output_scale, type = NULL){

  add_type  <- if(!is.null(type)) paste0("(",type,")") else NULL

  if(par == "mu"){

    par_names <- list(switch(
      output_scale,
      "r"     = bquote(~rho~.(add_type)),
      "d"     = bquote("Cohen's"~italic(d)~.(add_type)),
      "z"     = bquote("Fisher's"~italic(z)~.(add_type)),
      "logOR" = bquote("log"(italic("OR"))~.(add_type)),
      "OR"    = bquote(~italic("OR")~.(add_type)),
      "y"     = bquote(~mu~.(add_type))
    ))

  }else if(par == "tau"){

    par_names <- list(switch(
      output_scale,
      "r"     = bquote(~tau~("Fisher's"~italic(z)~"scale;"~.(type))),
      "OR"    = bquote(~tau~("log"(italic("OR"))~"scale;"~.(type))),
      bquote(~tau~.(add_type))
    ))

  }else if(par == "omega"){

    summary_info <- summary(fit)
    sum_all      <- summary_info[["averaged"]]
    omega_names  <- rownames(sum_all)[grepl(par, rownames(sum_all))]
    par_names    <- vector("list", length = length(omega_names))
    for(i in 1:length(par_names)){
      par_names[[i]] <- bquote(~omega[~.(substr(omega_names[i],6,nchar(omega_names[i])))])
    }

  }else if(par == "theta"){

    add_type <- if(!is.null(type)) paste0("(",type,")") else NULL
    par_names <- vector("list", length = length(study_names))
    for(i in 1:length(par_names)){
      par_names[i] <- list(switch(
        output_scale,
        "r"     = bquote(~rho[.(paste0("[",study_names[i],"]"))]~.(add_type)),
        "d"     = bquote("Cohen's"~italic(d)[.(paste0("[",study_names[i],"]"))]~.(add_type)),
        "z"     = bquote("Fisher's"~italic(z)[.(paste0("[",study_names[i],"]"))]~.(add_type)),
        "logOR" = bquote("log"(italic("OR"))[.(paste0("[",study_names[i],"]"))]~.(add_type)),
        "OR"    = bquote(~italic("OR")[.(paste0("[",study_names[i],"]"))]~.(add_type)),
        "y"     = bquote(~mu[.(paste0("[",study_names[i],"]"))]~.(add_type))
      ))
    }

  }

  return(par_names)
}
