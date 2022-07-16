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
#' \code{"mu"}, \code{"tau"}, \code{"omega"}, \code{"PET"},
#' or \code{"PEESE"}.
#' @param type type of MCMC diagnostic to be plotted.
#' Options are \code{"chains"} for the chains' trace plots,
#' \code{"autocorrelation"} for autocorrelation of the
#' chains, and \code{"densities"} for the overlaying
#' densities of the individual chains. Can be abbreviated to
#' first letters.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param show_models MCMC diagnostics of which models should be
#' plotted. Defaults to \code{NULL} which plots MCMC diagnostics
#' for a specified parameter for every model that is part of the
#' ensemble.
#' @param title whether the model number should be displayed in title.
#' Defaults to \code{TRUE} when more than one model is selected.
#' @param lags number of lags to be shown for
#' \code{type = "autocorrelation"}. Defaults to \code{30}.
#' @param ... additional arguments to be passed to
#' \link[graphics]{par} if \code{plot_type = "base"}.
#'
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
#'
#'
#' @return \code{diagnostics} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object/list of objects (depending on the number of parameters to be plotted)
#' of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()], [summary.RoBMA()]
#' @name diagnostics
#' @aliases diagnostics_autocorrelation diagnostics_trace diagnostics_density
#' @export diagnostics
#' @export diagnostics_density
#' @export diagnostics_autocorrelation
#' @export diagnostics_trace

#' @rdname diagnostics
diagnostics <- function(fit, parameter, type, plot_type = "base", show_models = NULL,
                  lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){

  # check settings
  if(!is.RoBMA(fit))
    stop("Diagnostics are available only for RoBMA models.")
  if(fit$add_info$save == "min")
    stop("Diagnostics cannot be produced because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' while fitting the model (see ?RoBMA for more details).")
  BayesTools::check_char(parameter, "parameter")
  BayesTools::check_char(type, "type")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  # deal with bad type names
  if(substr(type, 1, 1) == "c"){
    type <- "trace"
  }else if(substr(type, 1, 1) == "t"){
    type <- "trace" # for trace
  }else if(substr(type, 1, 1) == "d"){
    type <- "density"
  }else if(substr(type, 1, 1) == "a"){
    type <- "autocorrelation"
  }else{
    stop("Unsupported diagnostics type. See '?diagnostics' for more details.")
  }

  # deal with bad parameter names for PET-PEESE, weightfunction
  if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
    parameter         <- "omega"
    parameter_samples <- "omega"
  }else if(parameter %in% c("tau", "rho", "PET", "PEESE")){
    parameter         <- parameter
    parameter_samples <- parameter
  }else if(parameter == "mu"){
    parameter         <- parameter
    parameter_samples <- if(is.RoBMA.reg(fit)) "mu_intercept" else "mu"
  }else if(is.RoBMA.reg(fit) && parameter %in% fit$add_info[["predictors"]]){
    parameter         <- parameter
    parameter_samples <- .BayesTools_parameter_name(parameter)
  }else{
    if(is.RoBMA.reg(fit)){
      stop(paste0("The passed parameter does not correspond to any of main model parameter ('mu', 'tau', 'omega', 'PET', 'PEESE') or any of the specified predictors: ", paste0("'", fit$add_info[["predictors"]], "'", collapse = ", "), ". See '?plot.RoBMA' for more details."))
    }else{
      stop(paste0("The passed parameter does not correspond to any of main model parameter ('mu', 'tau', 'omega', 'PET', 'PEESE'). See '?plot.RoBMA' for more details."))
    }
  }

  # omit the first figure for publication bias weights (it is constants for all interesting weightfunctions)
  if(parameter == "omega"){
    show_figures <-  -1
  }else{
    show_figures <-  NULL
  }


  # do the plotting
  out <- list()

  models_ind <- 1:length(fit$models)
  if(!is.null(show_models)){
    models_ind <- models_ind[show_models]
  }

  # a message with info about multiple plots
  if(plot_type == "base" & (length(models_ind) > 1 | parameter == "omega"))
    message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")

  dots  <- .set_dots_diagnostics(..., type = type, chains = fit[["fit_control"]][["chains"]])
  plots <- list()

  for(i in models_ind){

    model_parameters <- c(names(attr(fit$models[[i]][["fit"]], "prior_list")))

    if(!parameter_samples %in% model_parameters){

      plots[[i]] <- NULL

    }else if(inherits(fit$models[[i]][["fit"]], "null_model")){

      plots[[i]] <- NULL

    }else{

      # get the parameter name
      args                          <- dots
      args$fit                      <- fit$models[[i]][["fit"]]
      args$parameter                <- parameter_samples
      args$parameter_names          <- if(parameter %in% c("mu", "tau")) .plot.RoBMA_par_names(parameter, fit, fit$add_info[["prior_scale"]])[[1]]
      args$type                     <- type
      args$plot_type                <- plot_type
      args$lags                     <- lags
      args$transformations          <- NULL
      args$transform_orthonormal    <- TRUE
      args$short_name               <- FALSE
      args$parameter_names          <- FALSE
      args$formula_prefix           <- FALSE

      if(!is.null(title) && title){
        args$main <- paste0("Model ", i)
      }

      plots[[i]] <- do.call(BayesTools::JAGS_diagnostics, args)
    }
  }



  # return the plots
  if(plot_type == "base"){
    return(invisible(plots))
  }else if(plot_type == "ggplot"){
    if(length(plots) == 1){
      plots <- plots[[1]]
    }
    return(plots)
  }
}


#' @rdname diagnostics
diagnostics_autocorrelation <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "autocorrelation", plot_type = plot_type, show_models = show_models, lags = lags, title = title, ...)
}

#' @rdname diagnostics
diagnostics_trace           <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "trace", plot_type = plot_type, show_models = show_models, title = title, ...)
}

#' @rdname diagnostics
diagnostics_density         <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "density", plot_type = plot_type, show_models = show_models, title = title, ...)
}

.set_dots_diagnostics  <- function(..., type, chains){

  dots <- list(...)
  if(is.null(dots[["col"]])){
    dots[["col"]]      <- if(type == "autocorrelation") "black" else rev(scales::viridis_pal()(chains))
  }

  return(dots)
}
