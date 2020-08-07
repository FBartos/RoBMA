#' @title Checks a fitted RoBMA object
#'
#' @description \code{diagnostics} creates visual
#' checks of individual models convergence. Numerical
#' overview of individual models can be obtained by
#' \code{summary(object, type = "models", diagnostics = TRUE)},
#' or even more detailed information by
#' \code{summary(object, type = "individual")}.
#'
#' @param fit a fitted RoBMA object
#' @param parameter a parameter to be plotted. Either
#' \code{"mu"}, \code{"tau"}, \code{"theta"}, or
#' \code{"omega"}.
#' @param type what type of model check should be plotted.
#' Options are \code{"chains"} for the chains trace plots,
#' \code{"autocorrelation"} for autocorrelation of the
#' chains, and \code{"densities"} for the overlaying
#' densities of the individual chains.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot2"} for plotting. The later
#' requires \pkg{ggplot2} package to be installed.
#' @param show_figures which figures should be returned in case of
#' multiple plots are generated. Useful when
#' \code{parameter = "omega"} when a plot for each parameter would be
#' generated. Can be also used for \code{parameter = "theta"} to
#' obtain only a specific subset of thetas. Set to \code{NULL}
#' to show all parameters (default for \code{parameter = "theta"}).
#' @param show_models diagnostics for which models should be produced.
#' Defaults to \code{NULL} that shows diagnostics to all models.
#' @param par_transform whether the figures should be produced for the
#' par_transform parameters. Defaults to \code{TRUE}.
#' @param title whether the model number should be displayed in title.
#' Defaults to \code{TRUE} when more than one model is selected.
#' @param lags number of lags to be shown for
#' \code{type = "autocorrelation"}. Defaults to \code{30}.
#' @param ... additional arguments to be passed to
#' \link[graphics]{par} if \code{plot_type = "base"}.
#' @details The visualization functions are based on
#' \link[rstan]{stan_plot} function and its color schemes.
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # diagnostics function allows to visualize diagnostics of a fitted RoBMA object, for example,
#' # the trace plot for the mean parameter in each model model
#' diagnostics(fit, parameter = "mu", type = "chain")
#'
#' # in order to show the trace plot only for the 11th model, add show_models parameter
#' diagnostics(fit, parameter = "mu", type = "chain", show_models = 11)
#'
#' # furthermore, the autocorrelations
#' diagnostics(fit, parameter = "mu", type = "autocorrelation")
#'
#' # and overlying densities for each plot can also be visualize
#' diagnostics(fit, parameter = "mu", type = "densities")
#' }
#' @export diagnostics
#' @seealso [RoBMA()], [summary.RoBMA()]
diagnostics <- function(fit, parameter, type, plot_type = "base",
                  show_figures = if(parameter == "omega") -1 else NULL, show_models = NULL, par_transform = TRUE,
                  lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){


  # check & correct parameters
  if(class(fit) != "RoBMA")stop("Diagnostics are available only for RoBMA models.")
  if(fit$add_info$save == "min")stop("Diagnostics cannot be produced because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBMA for more details).")

  if(length(parameter) != 1)stop("Only one parameter can be selected.")
  if(!parameter %in% c("mu", "tau", "omega", "theta"))stop("Selected parameter is not supported. Choose one of 'mu', 'tau', 'omega', or 'theta'.")
  if(length(type) != 1)stop("Only one type can be selected.")
  if(substr(type, 1, 1) == "c")type <- "chains"
  if(substr(type, 1, 1) == "t")type <- "chains" # for trace
  if(substr(type, 1, 1) == "d")type <- "densities"
  if(substr(type, 1, 1) == "a")type <- "autocorrelation"
  if(!type %in% c("chains", "autocorrelation", "densities"))stop("Selected diagnostic type is not supported. Choose one of 'chains', 'autocorrelation', or 'densities'.")
  if(plot_type == "ggplot2")plot_type <- "ggplot"
  if(!plot_type %in% c("base", "ggplot"))stop("The passed plot_type is not supported for plotting. See '?diagnostics' for more details.")

  # check availability of ggplot
  if(plot_type == "ggplot"){
    if(!try(requireNamespace("ggplot2")))stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  }


  # do the plotting
  out <- list()

  models_ind <- 1:length(fit$models)
  if(!is.null(show_models)){
    models_ind <- models_ind[show_models]
  }

  # a message with info about muliple plots
  if(plot_type == "base" & (length(models_ind) > 1 | parameter == "omega"))message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")

  for(m in models_ind){

    temp_out  <- NULL
    temp_data <- .diagnostics_plot_data(fit, m, parameter, par_transform)

    # deal with no parameter in model
    if(is.null(temp_data)){
      if(length(models_ind) == 1){
        message("Selected model does not containt the parameter of interest.")
        return(invisible(NULL))
      }else{
        out[m] <- temp_out
        next
      }
    }

    # make the plots
    par_ind <- 1:length(temp_data)
    if(!is.null(show_figures)){
      par_ind <- par_ind[show_figures]
    }

    for(i in par_ind){

      if(type == "chains")temp_out <- c(temp_out, list(.diagnostics_plot_trace(temp_data[[i]], plot_type, if(title) m else NULL, ...)))
      if(type == "densities")temp_out <- c(temp_out, list(.diagnostics_plot_density(temp_data[[i]], plot_type, if(title) m else NULL, parameter, ...)))
      if(type == "autocorrelation")temp_out <- c(temp_out, list(.diagnostics_plot_ac(temp_data[[i]], plot_type, if(title) m else NULL, lags, ...)))

    }

    if(length(temp_out) == 1)temp_out <- temp_out[[1]]

    if(length(models_ind) == 1){
      out    <- temp_out
    }else{
      out[m] <- list(temp_out)
    }

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible(NULL))
  }else if(plot_type == "ggplot"){
    return(out)
  }
}


.diagnostics_plot_trace   <- function(plot_data, plot_type, title, ...){

  if(plot_type == "base"){

    # save plotting settings
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4, 4, if(!is.null(title)) 3 else 1, 1))
    }else{
      graphics::par(list(...))
    }

    graphics::plot(NA, type = "n", xlim = range(plot_data$samp$iteration), ylim = range(plot_data$samp$value),
                   xlab = "", ylab = "", bty = "n", las = 1)
    for(i in as.numeric(unique(plot_data$samp$chain))){
      graphics::lines(plot_data$samp$iteration[plot_data$samp$chain == i], plot_data$samp$value[plot_data$samp$chain == i],
                      col = .diagnostics_color(plot_data$nchains)[i])
    }
    if(!is.null(title)){
      graphics::mtext(paste0("Model ",title), side = 3, line = 1, cex = 1.25)
    }
    graphics::mtext(plot_data$parameter, side = 2, line = 2.5, cex = 1.25)

    graph <- NULL

  }else if(plot_type == "ggplot"){

    graph <- ggplot2::ggplot(plot_data$samp, ggplot2::aes_string(x = "iteration", y = "value", color = "chain")) +
      ggplot2::geom_path() +
      ggplot2::scale_color_manual(values = .diagnostics_color(plot_data$nchains))
    temp_x_range <- range(plot_data$samp$iteration)
    temp_y_range <- range(plot_data$samp$value)
    graph <- graph + ggplot2::scale_x_continuous(
        name   = "Iterations",
        limits = temp_x_range,
        breaks = pretty(temp_x_range, n = 3),
        labels = pretty(temp_x_range, n = 3)
      ) +
      ggplot2::scale_y_continuous(
        name   = plot_data$parameter,
        limits = temp_y_range,
        breaks = pretty(temp_y_range),
        labels = pretty(temp_y_range)
      )
    if(!is.null(title))graph <- graph + ggplot2::ggtitle(paste0("Model ",title))
  }

  return(graph)
}
.diagnostics_plot_density <- function(plot_data, plot_type, title, par, ...){

  if(plot_type == "base"){

    # save plotting settings
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4, 4, if(!is.null(title)) 3 else 1, 1))
    }else{
      graphics::par(list(...))
    }
    with_trunc <- list()
    if(!is.infinite(plot_data$lower))with_trunc$from <- plot_data$lower
    if(!is.infinite(plot_data$upper))with_trunc$to   <- plot_data$upper

    temp_den <- vector(mode = "list", length = length(unique(plot_data$samp$chain)))
    for(i in as.numeric(unique(plot_data$samp$chain))){
      # deal with first weights if requested
      if(all(plot_data$samp$value[plot_data$samp$chain == i] == 1) & par == "omega"){
        temp_den[[i]] <- NULL
      }else{
        temp_den[[i]] <- do.call(stats::density, c(list(x = plot_data$samp$value[plot_data$samp$chain == i]), with_trunc))
      }
    }

    graphics::plot(NA, type = "n",
                   xlim = if(all(sapply(temp_den, is.null))) c(0, 1) else range(sapply(1:length(temp_den), function(i)temp_den[[i]]$x)),
                   ylim = if(all(sapply(temp_den, is.null))) c(0, 1) else c(0, max(sapply(1:length(temp_den), function(i)temp_den[[i]]$y))),
                   xlab = "", ylab = "", bty = "n", las = 1)
    for(i in 1:length(temp_den)){
      if(is.null(temp_den[[i]]) & par == "omega"){
        graphics::arrows(x0 = 1, y0 = 0, y1 = 1, lwd = 2, lty = 1, col = .diagnostics_color(plot_data$nchains)[i])
      }else{
        graphics::lines(temp_den[[i]], col = .diagnostics_color(plot_data$nchains)[i])
        graphics::polygon(x = c(if(!is.infinite(plot_data$lower)) plot_data$lower, temp_den[[i]]$x, if(!is.infinite(plot_data$upper)) plot_data$upper),
                          y = c(if(!is.infinite(plot_data$lower)) 0,               temp_den[[i]]$y, if(!is.infinite(plot_data$upper)) 0),
                          border = .diagnostics_color(plot_data$nchains)[i],
                          col    = scales::alpha(.diagnostics_color(plot_data$nchains)[i], alpha = .5))
      }
    }
    if(!is.null(title)){
      graphics::mtext(paste0("Model ",title), side = 3, line = 1, cex = 1.25)
    }
    graphics::mtext(if(all(sapply(temp_den, is.null))) "Probability" else "Density", side = 2, line = 2.5, cex = 1.25)
    graphics::mtext(plot_data$parameter,                                             side = 1, line = 2.5, cex = 1.25)

    graph <- NULL

  }else if(plot_type == "ggplot"){

    graph <-  ggplot2::ggplot(plot_data$samp, ggplot2::aes_string(x = "value")) +
      ggplot2::geom_density(mapping = ggplot2::aes_string(fill = "chain"), color = "black", alpha = .5) +
      ggplot2::scale_fill_manual(values = .diagnostics_color(plot_data$nchains))
    temp_y_max   <- max(ggplot2::ggplot_build(graph)$data[[1]]$density)
    temp_x_range <- if(par == "omega") c(0, 1) else range(plot_data$samp$value)
    graph <- graph +  ggplot2::scale_y_continuous(
        name   = "Density",
        limits = range(pretty(c(0, temp_y_max))),
        breaks = pretty(c(0, temp_y_max)),
        labels = pretty(c(0, temp_y_max))
      ) +
      ggplot2::scale_x_continuous(
        name   = plot_data$parameter,
        limits = range(pretty(temp_x_range)),
        breaks = pretty(temp_x_range),
        labels = pretty(temp_x_range)
      )
    if(!is.null(title))graph <- graph + ggplot2::ggtitle(paste0("Model ",title))

  }

  return(graph)
}
.diagnostics_plot_ac      <- function(plot_data, plot_type, title, lags = 30, ...){

  ac_dat   <- .diagnostics_ac_data(dat = plot_data$samp, lags = lags)

  if(plot_type == "base"){

    # save plotting settings
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    # set up margins
    if(length(list(...)) == 0){
      graphics::par(mar = c(4,4,3,1))
    }else{
      graphics::par(list(...))
    }

    temp_dat <- as.numeric(by(ac_dat$ac, ac_dat$lag, mean))
    temp_dat[is.nan(temp_dat)] <- 1
    graphics::barplot(temp_dat, names.arg = unique(ac_dat$lag), col = "#B2001D", las = 1)
    graphics::mtext("Lag",                  side = 1, line = 2.5, cex = 1.25)
    graphics::mtext("Avg. autocorrelation", side = 2, line = 2.5, cex = 1.25)
    if(!is.null(title)){
      graphics::mtext(bquote(.(paste0("Model ",title))~": "~.(plot_data$parameter)), side = 3, line = 1, cex = 1.25)
    }else{
      graphics::mtext(plot_data$parameter, side = 3, line = 1, cex = 1.25)
    }

    graph <- NULL

  }else if(plot_type == "ggplot"){
    graph     <- ggplot2::ggplot(ac_dat, ggplot2::aes_string(x = "lag", y = "ac")) +
      ggplot2::geom_bar(size = .5, color = "black", fill = "#B2001D", position = "dodge", stat = "summary", fun.y = "mean") +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
      ggplot2::labs(x = "Lag", y = "Avg. autocorrelation")
    if(!is.null(title)){
      graph <- graph + ggplot2::ggtitle( bquote(.(paste0("Model ",title))~": "~.(plot_data$parameter)))
    }else{
      graph <- graph + ggplot2::ggtitle(plot_data$parameter)
    }
  }

  return(graph)
}
.diagnostics_ac_data      <- function(dat, lags){
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))

  ac_list <- tapply(dat$value, INDEX = ch, FUN = function(x)stats::acf(x, lag.max = lags, plot = FALSE)$acf[, , 1L], simplify = FALSE)

  nl <- lags + 1
  ch <- factor(rep(1:nc, each = nl), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(0, lags), nc)

  return(data.frame(chains = ch, ac = do.call(c, ac_list), lag = ll))
}
.diagnostics_color        <- function(n){
  return(rep_len(c("#E66101", "#998EC3", "#542788", "#F1A340", "#D8DAEB", "#FEE0B6"), n))
}
.diagnostics_plot_data    <- function(fit, model, par, par_transform){

  if(length(fit$models[[model]]$fit) == 0){

    return(NULL)

  }else{

    samples <- coda::as.array.mcmc.list(fit$models[[model]]$fit$mcmc, drop = FALSE)

    if(!any(grepl(par, dimnames(samples)$var)))return(NULL)

    if(par %in% c("mu", "tau")){
      ind <- c(1:length(dimnames(samples)$var))[par == dimnames(samples)$var]
    }else{
      ind <- c(1:length(dimnames(samples)$var))[grepl(par, dimnames(samples)$var)]
    }

    if(par == "omega")ind <- rev(ind)


    # create parameter names
    par_names <- .plot.RoBMA_par_names(par, fit)


    plot_data <- list()


    for(i in 1:length(ind)){
      plot_data[[dimnames(samples)$var[ind[i]]]] <- list(
        samp = data.frame(
          value     = as.vector(samples[,ind[i],]),
          parameter = dimnames(samples)$var[ind[i]],
          chain     = as.factor(c(unlist(sapply(1:dim(samples)[3], function(x)rep(x,dim(samples)[1]))))),
          iteration = rep(1:dim(samples)[1], dim(samples)[3])
        ),
        nchains   = dim(samples)[3],
        nparams   = 1,
        warmup    = 0,
        parameter = par_names[[i]],
        lower     = if(par == "omega") 0 else fit$models[[model]]$priors[[par]]$truncation$lower,
        upper     = if(par == "omega") 1 else fit$models[[model]]$priors[[par]]$truncation$upper
      )
    }

    # transform the values if requested
    if(par_transform){
      if(par %in% c("mu", "theta") & fit$add_info$effect_size %in% c("r", "OR")){
        plot_data[[1]]$samp$value <- .transform(plot_data[[1]]$samp$value, fit$add_info$effect_size, fit$add_info$mu_transform)
      }
    }

  }

  return(plot_data)
}
