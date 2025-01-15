#' @title Summarize marginal estimates of a fitted RoBMA regression object
#'
#' @description \code{marginal_summary} creates summary tables for
#' marginal estimates of a RoBMA regression model.
#'
#' @param object a fitted RoBMA regression object
#' @param conditional show the conditional estimates (assuming that the
#' alternative is true).
#' @inheritParams summary.RoBMA
#'
#'
#' @return \code{marginal_summary} returns a list of tables of class 'BayesTools_table'.
#'
#' @seealso [RoBMA()], [summary.RoBMA()], [diagnostics()], [check_RoBMA()]
#' @export
marginal_summary <- function(object, conditional = FALSE,
                             output_scale = NULL, probs = c(.025, .975), logBF = FALSE, BF01 = FALSE){

  # apply version changes to RoBMA object
  object <- .update_object(object)

  if(!is.RoBMA.reg(object))
    stop("'marginal_summary' function is available only for RoBMA regression models")
  if(sum(.get_model_convergence(object)) == 0)
    stop("There is no converged model in the ensemble.")
  if(!is.null(.check_predictors_scaling(object)))
    stop("'marginal_summary' function requires standardized predictors")

  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(BF01,  "BF01")
  BayesTools::check_bool(logBF, "logBF")

  # check the scales
  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale) != "y"){
    stop("Models estimated using the generall effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale)
  }


  # transform the estimates if needed
  if(object$add_info[["output_scale"]] != output_scale){
    object <- .transform_posterior(object, object$add_info[["output_scale"]], output_scale)
  }

  # obtain table
  estimates <- BayesTools::marginal_estimates_table(
    samples        = object$RoBMA[["inference_marginal"]][["averaged"]],
    inference      = object$RoBMA[["inference_marginal"]][["inference"]],
    parameters     = names(object$RoBMA[["inference_marginal"]][["inference"]]),
    logBF          = logBF,
    BF01           = BF01,
    formula_prefix = FALSE,
    probs          = probs,
    title          = "Model-averaged marginal estimates:",
    footnotes      = .scale_note(object$add_info[["prior_scale"]], output_scale, marginal = TRUE),
    warnings       = .collect_errors_and_warnings(object)
  )


  # create the output object
  output <- list(
    call       = object[["call"]],
    title      = "Robust Bayesian meta-analysis",
    estimates  = estimates
  )

  if(conditional){

    estimates_conditional <- BayesTools::marginal_estimates_table(
      samples        = object$RoBMA[["inference_marginal"]][["conditional"]],
      inference      = object$RoBMA[["inference_marginal"]][["inference"]],
      parameters     = names(object$RoBMA[["inference_marginal"]][["inference"]]),
      logBF          = logBF,
      BF01           = BF01,
      formula_prefix = FALSE,
      probs          = probs,
      title          = "Conditional marginal estimates:",
      footnotes      = .scale_note(object$add_info[["prior_scale"]], output_scale, marginal = TRUE),
      warnings       = .collect_errors_and_warnings(object)
    )

    output$estimates_conditional <- estimates_conditional
  }


  class(output) <- "marginal_summary.RoBMA"

  return(output)
}


#' @title Prints marginal_summary object for RoBMA method
#'
#' @param x a summary of a RoBMA object
#' @param ... additional arguments
#'
#'
#' @return \code{print.marginal_summary.RoBMA} invisibly returns the print statement.
#'
#' @seealso [RoBMA()]
#' @export
print.marginal_summary.RoBMA <- function(x, ...){

  cat("Call:\n")
  print(x[["call"]])

  cat("\n")
  cat(x[["title"]])

  cat("\n")
  print(x[["estimates"]])

  if(!is.null(x[["estimates_conditional"]])){
    cat("\n")
    print(x[["estimates_conditional"]])
  }


  return(invisible())
}


#' @title Plots marginal estimates of a fitted RoBMA regression object
#'
#' @description \code{marginal_plot} allows to visualize prior and
#' posterior distributions of marginal estimates of a RoBMA regression model.
#'
#' @param x a fitted RoBMA regression object
#' @param parameter regression parameter to be plotted
#' @param conditional whether conditional marginal estimates should be
#' plotted. Defaults to \code{FALSE} which plots the model-averaged
#' estimates.
#' @inheritParams plot.RoBMA
#'
#'
#' @return \code{plot.RoBMA} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
marginal_plot  <- function(x, parameter, conditional = FALSE, plot_type = "base", prior = FALSE, output_scale = NULL, dots_prior = NULL, ...){

  # apply version changes to RoBMA object
  x <- .update_object(x)

  # check whether plotting is possible
  if(!is.RoBMA.reg(x))
    stop("'marginal_plot' function is available only for RoBMA regression models")
  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")
  if(!is.null(.check_predictors_scaling(x)))
    stop("'marginal_plot' function requires standardized predictors")

  # check settings
  BayesTools::check_char(parameter, "parameter", allow_values = x$add_info[["predictors"]])
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_bool(prior, "prior")

  ### manage transformations
  # get the settings
  results_scale <- x$add_info[["output_scale"]]
  if(is.null(output_scale)){
    output_scale <- x$add_info[["output_scale"]]
  }else{
    output_scale <- .transformation_var(output_scale)
  }
  # transform the estimates if needed
  if(x$add_info[["output_scale"]] != output_scale){
    x <- .transform_posterior(x, x$add_info[["output_scale"]], output_scale)
  }


  # choose the samples
  if(conditional){
    samples <- x$RoBMA[["inference_marginal"]][["conditional"]]
  }else{
    samples <- x$RoBMA[["inference_marginal"]][["averaged"]]
  }

  dots       <- .set_dots_plot(..., n_levels = length(samples[[.BayesTools_parameter_name(parameter)]]))
  dots_prior <- .set_dots_prior_marginal(dots_prior, n_levels = length(samples[[.BayesTools_parameter_name(parameter)]]))

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- .BayesTools_parameter_name(parameter)
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$transformation           <- NULL
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$par_name                 <- .plot.RoBMA_par_names(parameter, x, output_scale)[[1]]
  args$dots_prior               <- dots_prior

  plot <- suppressMessages(do.call(BayesTools::plot_marginal, args))


  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


.set_dots_prior_marginal <- function(dots_prior, n_levels){

  if(is.null(dots_prior)){
    dots_prior <- list()
  }

  if(is.null(dots_prior[["col"]]) & n_levels == 1){
    dots_prior[["col"]]      <- "black"
  }else if(is.null(dots_prior[["col"]]) & n_levels > 1){
    dots_prior[["col"]]      <- grDevices::palette.colors(n = n_levels + 1, palette = "Okabe-Ito")[-1]
  }
  if(is.null(dots_prior[["lty"]])){
    dots_prior[["lty"]]      <- 2
  }

  return(dots_prior)
}
