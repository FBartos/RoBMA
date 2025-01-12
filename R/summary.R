#' @title Prints a fitted RoBMA object
#'
#' @param x a fitted RoBMA object.
#' @param ... additional arguments.
#'
#'
#' @return \code{print.RoBMA} invisibly returns the print statement.
#'
#' @seealso [RoBMA()]
#' @export
print.RoBMA <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


#' @title Summarize fitted RoBMA object
#'
#' @description \code{summary.RoBMA} creates summary tables for a
#' RoBMA object.
#'
#' @param object a fitted RoBMA object
#' @param type whether to show the overall RoBMA results (\code{"ensemble"}),
#' an overview of the individual models (\code{"models"}), an overview of
#' the individual models MCMC diagnostics (\code{"diagnostics"}), or a detailed summary
#' of the individual models (\code{"individual"}). Can be abbreviated to first letters.
#' @param conditional show the conditional estimates (assuming that the
#' alternative is true). Defaults to \code{FALSE}. Only available for
#' \code{type == "ensemble"}.
#' @param output_scale transform the  meta-analytic estimates to a different
#' scale. Defaults to \code{NULL} which returns the same scale as the model was estimated on.
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .975)}
#' @param logBF show log of Bayes factors. Defaults to \code{FALSE}.
#' @param BF01 show Bayes factors in support of the null hypotheses. Defaults to
#' \code{FALSE}.
#' @param short_name whether priors names should be shortened to the first
#' (couple) of letters. Defaults to \code{FALSE}.
#' @param remove_spike_0 whether spike prior distributions with location at zero should
#' be omitted from the summary. Defaults to \code{FALSE}.
#' @param ... additional arguments
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' # summary can provide many details about the model
#' summary(fit)
#'
#' # estimates from the conditional models can be obtained with
#' summary(fit, conditional = TRUE)
#'
#' # overview of the models and their prior and posterior probability, marginal likelihood,
#' # and inclusion Bayes factor can be obtained with
#' summary(fit, type = "models")
#'
#' # diagnostics overview, containing the maximum R-hat, minimum ESS, maximum MCMC error, and
#' # maximum MCMC error / sd across parameters for each individual model can be obtained with
#' summary(fit, type = "diagnostics")
#'
#' # summary of individual models and their parameters can be further obtained by
#' summary(fit, type = "individual")
#' }
#'
#' @note See [diagnostics()] for visual convergence checks of the individual models.
#'
#'
#' @return \code{summary.RoBMA} returns a list of tables of class 'BayesTools_table'.
#'
#' @seealso [RoBMA()], [diagnostics()], [check_RoBMA()]
#' @export
summary.RoBMA       <- function(object, type = "ensemble", conditional = FALSE,
                                output_scale = NULL, probs = c(.025, .975), logBF = FALSE, BF01 = FALSE,
                                short_name = FALSE, remove_spike_0 = FALSE, ...){

  # apply version changes to RoBMA object
  object <- .update_object(object)

  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(type, "type")
  BayesTools::check_char(output_scale, "output_scale", allow_NULL = TRUE)
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(BF01,  "BF01")
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(short_name, "short_name")
  BayesTools::check_bool(remove_spike_0, "remove_spike_0")

  # check the scales
  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale, estimation = FALSE) != "y"){
    stop("Models estimated using the generall effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale, estimation = FALSE)
  }

  # print diagnostics if all models fail to converge
  if(!any(.get_model_convergence(object))){
    if(substr(type,1,1) != "d")
      warning("All models failed to converge. Model diagnostics were printed instead.")
    type        <- "diagnostics"
  }

  # dispatch the spike and slab summary function if needed
  if(object[["add_info"]][["algorithm"]] == "ss"){
    return(.summary.RoBMA.ss(object, type, conditional, output_scale, probs, logBF, BF01, remove_spike_0, ...))
  }

  if(substr(type,1,1) == "e"){

    # transform the estimates if needed
    if(object$add_info[["output_scale"]] != output_scale){
      object <- .transform_posterior(object, object$add_info[["output_scale"]], output_scale)
    }

    # obtain components overview
    components <- BayesTools::ensemble_inference_table(
      inference  = object$RoBMA[["inference"]],
      parameters = names(object$RoBMA[["inference"]])[names(object$RoBMA[["inference"]]) %in% c("Effect", "Heterogeneity", "Bias", "Hierarchical", "Baseline")],
      logBF      = logBF,
      BF01       = BF01,
      title      = "Components summary:"
    )

    # obtain estimates tables
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBMA[["posteriors"]],
      parameters = names(object$RoBMA[["posteriors"]]),
      probs      = probs,
      title      = "Model-averaged estimates:",
      footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object)),
      warnings   = .collect_errors_and_warnings(object)
    )

    # deal with possibly empty table in case of no alternative models
    if(is.null(object$RoBMA[["posteriors_conditional"]])){
      estimates_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates_conditional)          <- c("Mean", "Median", probs)
      class(estimates_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_conditional))
      attr(estimates_conditional, "type")      <- rep("estimate", ncol(estimates_conditional))
      attr(estimates_conditional, "rownames")  <- TRUE
      attr(estimates_conditional, "title")     <- "Conditional estimates:"
      attr(estimates_conditional, "footnotes") <- c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object))
      attr(estimates_conditional, "warnings")  <- .collect_errors_and_warnings(object)
    }else{
      estimates_conditional <- BayesTools::ensemble_estimates_table(
        samples    = object$RoBMA[["posteriors_conditional"]],
        parameters = names(object$RoBMA[["posteriors_conditional"]]),
        probs      = probs,
        title      = "Conditional estimates:",
        footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object)),
        warnings   = .collect_errors_and_warnings(object)
      )
    }

    # create the output object
    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      components = components,
      estimates  = estimates
    )

    if(conditional){
      output$estimates_conditional <- estimates_conditional
    }

    # add meta-regression summaries
    if(is.RoBMA.reg(object)){

      # rename the inference components
      for(i in seq_along(object$RoBMA[["inference_predictors"]])){
        attr(object$RoBMA[["inference_predictors"]][[i]], "parameter_name") <- gsub("(mu) ", "", attr(object$RoBMA[["inference_predictors"]][[i]], "parameter_name"), fixed = TRUE)
      }

      if(!is.null(object$RoBMA[["inference_predictors"]])){
        output$components_predictors <- BayesTools::ensemble_inference_table(
          inference  = object$RoBMA[["inference_predictors"]],
          parameters = names(object$RoBMA[["inference_predictors"]]),
          logBF      = logBF,
          BF01       = BF01,
          title      = "Meta-regression components summary:"
        )
      }

      if(!is.null(object$RoBMA[["posteriors_predictors"]])){
        # obtain estimates tables
        output$estimates_predictors <- BayesTools::ensemble_estimates_table(
          samples        = object$RoBMA[["posteriors_predictors"]],
          parameters     = names(object$RoBMA[["posteriors_predictors"]]),
          probs          = probs,
          title          = "Model-averaged meta-regression estimates:",
          formula_prefix = FALSE,
          footnotes      = .scale_note(object$add_info[["prior_scale"]], output_scale),
          warnings       = .collect_errors_and_warnings(object)
        )
      }


      # deal with possibly empty table in case of no alternative models
      if(conditional){
        if(is.null(object$RoBMA[["posteriors_predictors_conditional"]])){
          estimates_predictors_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
          colnames(estimates_predictors_conditional)          <- c("Mean", "Median", probs)
          class(estimates_predictors_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_predictors_conditional))
          attr(estimates_predictors_conditional, "type")      <- rep("estimate", ncol(estimates_predictors_conditional))
          attr(estimates_predictors_conditional, "rownames")  <- TRUE
          attr(estimates_predictors_conditional, "title")     <- "Conditional meta-regression estimates:"
          attr(estimates_predictors_conditional, "footnotes") <- .scale_note(object$add_info[["prior_scale"]], output_scale)
          attr(estimates_predictors_conditional, "warnings")  <- .collect_errors_and_warnings(object)
        }else{
          estimates_predictors_conditional <- BayesTools::ensemble_estimates_table(
            samples    = object$RoBMA[["posteriors_predictors_conditional"]],
            parameters = names(object$RoBMA[["posteriors_predictors_conditional"]]),
            probs      = probs,
            title      = "Conditional meta-regression estimates:",
            formula_prefix = FALSE,
            footnotes  = .scale_note(object$add_info[["prior_scale"]], output_scale),
            warnings   = .collect_errors_and_warnings(object)
          )
        }
        output$estimates_predictors_conditional <- estimates_predictors_conditional
      }

    }

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "ensemble"

    return(output)

  }else if(substr(type,1,1) == "m"){

    components <- names(object$RoBMA[["inference"]])[names(object$RoBMA[["inference"]]) %in% c("Effect", "Heterogeneity", "Bias", "Hierarchical", "Baseline")]

    parameters <- list()
    if(any(components == "Effect")){
      if(.is_regression(object)){
        parameters[["intercept"]] <- "mu_intercept"
        for(i in seq_along(object$add_info[["predictors"]])){
          parameters[[object$add_info[["predictors"]][i]]] <- .BayesTools_parameter_name(object$add_info[["predictors"]][i])
        }
      }else{
        parameters[["Effect"]] <- "mu"
      }
    }
    if(any(components == "Heterogeneity")){
      parameters[["Heterogeneity"]] <- "tau"
      if(.is_multivariate(object)){
        parameters[["Hierarchical"]] <- "rho"
      }
    }
    if(any(components == "Bias")){
      parameters[["Bias"]] <- c("PET", "PEESE", "omega")
    }
    if(is.BiBMA(object)){
      parameters[["Baseline"]] <- "pi"
    }

    summary <- BayesTools::ensemble_summary_table(
      models         = object[["models"]],
      parameters     = parameters,
      title          = "Models overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )

    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      summary    = summary
    )

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "models"

    return(output)

  }else if(substr(type,1,1) == "d"){

    components <- names(object$RoBMA[["inference"]])[names(object$RoBMA[["inference"]]) %in% c("Effect", "Heterogeneity", "Bias", "Hierarchical", "Baseline")]

    parameters <- list()
    if(any(components == "Effect")){
      if(.is_regression(object)){
        parameters[["intercept"]] <- "mu_intercept"
        for(i in seq_along(object$add_info[["predictors"]])){
          parameters[[object$add_info[["predictors"]][i]]] <- .BayesTools_parameter_name(object$add_info[["predictors"]][i])
        }
      }else{
        parameters[["Effect"]] <- "mu"
      }
    }
    if(any(components == "Heterogeneity")){
      parameters[["Heterogeneity"]] <- "tau"
      if(.is_multivariate(object)){
        parameters[["Hierarchical"]] <- "rho"
      }
    }
    if(any(components == "Bias")){
      parameters[["Bias"]] <- c("PET", "PEESE", "omega")
    }
    if(is.BiBMA(object)){
      parameters[["Baseline"]] <- "pi"
    }

    diagnostics <- BayesTools::ensemble_diagnostics_table(
      models         = object[["models"]],
      parameters     = parameters,
      title          = "Diagnostics overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )

    output <- list(
      call        = object[["call"]],
      title       = .object_title(object),
      diagnostics = diagnostics
    )

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "diagnostics"

    return(output)

  }else if(substr(type, 1, 1) == "i"){

    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      models     = list()
    )

    for(i in seq_along(object[["models"]])){

      summary  <- BayesTools::model_summary_table(
        model             = object[["models"]][[i]],
        remove_parameters = "gamma",
        short_name        = short_name,
        remove_spike_0    = remove_spike_0
      )
      if(output_scale == "y"){
        estimates <- object[["models"]][[i]][["fit_summary"]]
        attr(estimates, "warnings")  <- object[["models"]][[i]][["warnings"]]
        attr(estimates, "title")     <- "Parameter estimates:"
      }else{
        estimates <- object[["models"]][[i]][["fit_summaries"]][[output_scale]]
        attr(estimates, "footnotes") <- .scale_note(object[["models"]][[i]][["prior_scale"]], output_scale)
        attr(estimates, "warnings")  <- object[["models"]][[i]][["warnings"]]
        attr(estimates, "title")     <- "Parameter estimates:"
      }

      output[["models"]][[i]] <- list(
        summary   = summary,
        estimates = estimates
      )
    }

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "individual"

    return(output)

  }else{
    stop(paste0("Unknown summary type: '", type, "'."))
  }
}

.summary.RoBMA.ss    <- function(object, type, conditional, output_scale, probs, logBF, BF01, remove_spike_0, ...){

  if(substr(type,1,1) == "e"){

    # transform the estimates if needed
    if(object$add_info[["output_scale"]] != output_scale){
      object <- .transform_posterior(object, object$add_info[["output_scale"]], output_scale)
    }

    # obtain components overview
    components <- update(
      object[["RoBMA"]][["inference"]],
      title      = "Components summary:",
      BF01       = BF01,
      logBF      = logBF
    )

    # obtain estimates tables
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBMA[["posteriors"]],
      parameters = names(object$RoBMA[["posteriors"]]),
      probs      = probs,
      title      = "Model-averaged estimates:",
      footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object)),
      warnings   = .collect_errors_and_warnings(object)
    )

    # deal with possibly empty table in case of no alternative models
    if(is.null(object$RoBMA[["posteriors_conditional"]])){
      estimates_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates_conditional)          <- c("Mean", "Median", probs)
      class(estimates_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_conditional))
      attr(estimates_conditional, "type")      <- rep("estimate", ncol(estimates_conditional))
      attr(estimates_conditional, "rownames")  <- TRUE
      attr(estimates_conditional, "title")     <- "Conditional estimates:"
      attr(estimates_conditional, "footnotes") <- c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object))
      attr(estimates_conditional, "warnings")  <- .collect_errors_and_warnings(object)
    }else{
      estimates_conditional <- BayesTools::ensemble_estimates_table(
        samples    = object$RoBMA[["posteriors_conditional"]],
        parameters = names(object$RoBMA[["posteriors_conditional"]]),
        probs      = probs,
        title      = "Conditional estimates:",
        footnotes  = c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object)),
        warnings   = .collect_errors_and_warnings(object)
      )
    }

    # create the output object
    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      components = components,
      estimates  = estimates
    )

    if(conditional){
      output$estimates_conditional <- estimates_conditional
    }

    # add meta-regression summaries
    if(is.RoBMA.reg(object)){

      if(!is.null(object$RoBMA[["inference_predictors"]])){
        output$components_predictors <- update(
          object[["RoBMA"]][["inference_predictors"]],
          title      = "Meta-regression components summary:",
          BF01       = BF01,
          logBF      = logBF
        )
      }

      if(!is.null(object$RoBMA[["posteriors_predictors"]])){
        # obtain estimates tables
        output$estimates_predictors <- BayesTools::ensemble_estimates_table(
          samples        = object$RoBMA[["posteriors_predictors"]],
          parameters     = names(object$RoBMA[["posteriors_predictors"]]),
          probs          = probs,
          title          = "Model-averaged meta-regression estimates:",
          formula_prefix = FALSE,
          footnotes      = .scale_note(object$add_info[["prior_scale"]], output_scale),
          warnings       = .collect_errors_and_warnings(object)
        )
      }

      # deal with possibly empty table in case of no alternative models
      if(conditional){
        if(is.null(object$RoBMA[["posteriors_predictors_conditional"]])){
          estimates_predictors_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
          colnames(estimates_predictors_conditional)          <- c("Mean", "Median", probs)
          class(estimates_predictors_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_predictors_conditional))
          attr(estimates_predictors_conditional, "type")      <- rep("estimate", ncol(estimates_predictors_conditional))
          attr(estimates_predictors_conditional, "rownames")  <- TRUE
          attr(estimates_predictors_conditional, "title")     <- "Conditional meta-regression estimates:"
          attr(estimates_predictors_conditional, "footnotes") <- .scale_note(object$add_info[["prior_scale"]], output_scale)
          attr(estimates_predictors_conditional, "warnings")  <- .collect_errors_and_warnings(object)
        }else{
          estimates_predictors_conditional <- BayesTools::ensemble_estimates_table(
            samples    = object$RoBMA[["posteriors_predictors_conditional"]],
            parameters = names(object$RoBMA[["posteriors_predictors_conditional"]]),
            probs      = probs,
            title      = "Conditional meta-regression estimates:",
            formula_prefix = FALSE,
            footnotes  = .scale_note(object$add_info[["prior_scale"]], output_scale),
            warnings   = .collect_errors_and_warnings(object)
          )
        }
        output$estimates_predictors_conditional <- estimates_predictors_conditional
      }
    }

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "ensemble"

    return(output)

  }else{

    # transform the estimates if needed
    if(object$add_info[["output_scale"]] != output_scale){
      transformations <- .runjags_summary_list(NULL, priors = object$model[["priors"]], prior_scale = object$add_info[["prior_scale"]], measures = output_scale, transformations_only = TRUE)
    }

    estimates <- suppressMessages(BayesTools::runjags_estimates_table(
      object[["model"]][["fit"]],
      transformations   = if(object$add_info[["output_scale"]] != output_scale) transformations else NULL,
      transform_factors = TRUE,
      remove_spike_0    = remove_spike_0,
      formula_prefix    = FALSE,
      remove_inclusion  = TRUE,
      conditional       = FALSE,
      title             = "Model-averaged estimates:",
      warnings          = .collect_errors_and_warnings(object),
      footnotes         = c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object)),
    ))

    if(conditional){
      estimates_conditional <- suppressMessages(BayesTools::runjags_estimates_table(
        object[["model"]][["fit"]],
        transformations   = if(object$add_info[["output_scale"]] != output_scale) transformations else NULL,
        transform_factors = TRUE,
        remove_spike_0    = remove_spike_0,
        formula_prefix    = FALSE,
        remove_inclusion  = TRUE,
        conditional       = TRUE,
        title             = "Conditional estimates:",
        warnings          = .collect_errors_and_warnings(object),
        footnotes         = c(.scale_note(object$add_info[["prior_scale"]], output_scale), .note_omega(object)),
      ))
    }

    # create the output object
    output <- list(
      call       = object[["call"]],
      title      = .object_title(object),
      estimates  = estimates
    )

    if(conditional){
      output$estimates_conditional <- estimates_conditional
    }

    class(output) <- "summary.RoBMA"
    attr(output, "type") <- "ensemble"

    return(output)

  }
}

#' @title Prints summary object for RoBMA method
#'
#' @param x a summary of a RoBMA object
#' @param ... additional arguments
#'
#'
#' @return \code{print.summary.RoBMA} invisibly returns the print statement.
#'
#' @seealso [RoBMA()]
#' @export
print.summary.RoBMA <- function(x, ...){

  cat("Call:\n")
  print(x[["call"]])

  cat("\n")
  cat(x[["title"]])


  if(attr(x, "type") == "ensemble"){

    for(type in c("components", "components_predictors", "estimates", "estimates_predictors", "estimates_conditional", "estimates_predictors_conditional")){
      if(!is.null(x[[type]])){
        cat("\n")
        print(x[[type]])
      }
    }

    return(invisible())

  }else if(attr(x, "type") == "models"){

    cat("\n")
    print(x[["summary"]])

    return(invisible())

  }else if(attr(x, "type") == "diagnostics"){

    cat("\n")
    print(x[["diagnostics"]])

    return(invisible())

  }else if(attr(x, "type") == "individual"){

    for(i in seq_along(x[["models"]])){

      if(i > 1){
        cat("\n")
      }
      print(x[["models"]][[i]][["summary"]])

      cat("\n")
      print(x[["models"]][[i]][["estimates"]])
    }

    return(invisible())
  }
}


#' @title Reports whether x is a RoBMA object
#'
#' @param x an object to test
#'
#'
#' @return returns a boolean.
#'
#' @name is.RoBMA
#' @aliases is.RoBMA.reg is.NoBMA is.NoBMA.reg is.BiBMA
#' @export is.RoBMA
#' @export is.RoBMA.reg
#' @export is.NoBMA
#' @export is.NoBMA.reg
#' @export is.BiBMA

#' @rdname is.RoBMA
is.RoBMA            <- function(x){
  inherits(x, "RoBMA")
}
#' @rdname is.RoBMA
is.RoBMA.reg        <- function(x){
  inherits(x, "RoBMA.reg")
}
#' @rdname is.RoBMA
is.NoBMA            <- function(x){
  inherits(x, "NoBMA")
}
#' @rdname is.RoBMA
is.NoBMA.reg        <- function(x){
  inherits(x, "NoBMA.reg")
}
#' @rdname is.RoBMA
is.BiBMA            <- function(x){
  inherits(x, "BiBMA")
}


.object_title <- function(object){
  if(is.NoBMA(object)){
    return("Bayesian model-averaged meta-analysis (normal-normal model)")
  }else if(is.NoBMA.reg(object)){
    return("Bayesian model-averaged meta-regression (normal-normal model)")
  }else if(is.BiBMA(object)){
    return("Bayesian model-averaged meta-analysis (binomial-normal model)")
  }else if(is.RoBMA.reg(object)){
    return("Robust Bayesian meta-regression")
  }else if(is.RoBMA(object)){
    return("Robust Bayesian meta-analysis")
  }
}



#' @title Interprets results of a RoBMA model.
#'
#' @description \code{interpret} creates a brief textual summary
#' of a fitted RoBMA object.
#'
#' @inheritParams summary.RoBMA
#'
#'
#' @return \code{interpret} returns a character.
#'
#' @export
interpret           <- function(object, output_scale = NULL){

  .check_is_any_RoBMA_object(object)

  if(is.null(output_scale)){
    output_scale <- object$add_info[["output_scale"]]
  }else if(object$add_info[["output_scale"]] == "y" & .transformation_var(output_scale) != "y"){
    stop("Models estimated using the generall effect size scale 'y' / 'none' cannot be transformed to a different effect size scale.")
  }else{
    output_scale <- .transformation_var(output_scale)
  }

  if(object$add_info[["output_scale"]] != output_scale){
    object <- .transform_posterior(object, object$add_info[["output_scale"]], output_scale)
  }

  if(object[["add_info"]][["algorithm"]] == "bridge"){

    specification <- list()

    if(any(names(object$RoBMA[["inference"]]) == "Effect")){
      specification[["Effect"]] <- list(
        inference           = "Effect",
        samples             = "mu",
        inference_name      = "effect",
        inference_BF_name   = "BF_10",
        samples_name        = .transformation_names(object$add_info[["output_scale"]])
      )
    }
    if(any(names(object$RoBMA[["inference"]]) == "Heterogeneity")){
      specification[["Heterogeneity"]] <- list(
        inference           = "Heterogeneity",
        samples             = "tau",
        inference_name      = "heterogeneity",
        inference_BF_name   = "BF^rf",
        samples_name        = "tau"
      )
    }
    if(any(names(object$RoBMA[["inference"]]) == "Bias")){
      specification[["Bias"]] <- list(
        inference           = "Bias",
        inference_name      = "publication bias",
        inference_BF_name   = "BF_pb"
      )
    }

    text <- BayesTools::interpret(
      inference     = object$RoBMA[["inference"]],
      samples       = object$RoBMA[["posteriors"]],
      specification = specification,
      method        = .object_title(object)
    )

  }else{

    specification <- list()
    components    <- summary(object)[["components"]]

    if(any(rownames(components) == "Effect")){
      specification[["Effect"]] <- list(
        inference_name        = "effect",
        inference_BF_name     = "BF_10",
        inference_BF          = components[rownames(components) == "Effect", "inclusion_BF"],
        estimate_name         = .transformation_names(object$add_info[["output_scale"]]),
        estimate_samples      = as.numeric(object$RoBMA$posteriors[["mu"]]),
        estimate_conditional  = FALSE
      )
    }
    if(any(rownames(components) == "Heterogeneity")){
      specification[["Heterogeneity"]] <- list(
        inference_name        = "heterogeneity",
        inference_BF_name     = "BF^rf",
        inference_BF          = components[rownames(components) == "Heterogeneity", "inclusion_BF"],
        estimate_name         = "tau",
        estimate_samples      = as.numeric(object$RoBMA$posteriors[["tau"]]),
        estimate_conditional  = FALSE
      )
    }
    if(any(rownames(components) == "Bias")){
      specification[["Bias"]] <- list(
        inference_name      = "publication bias",
        inference_BF_name   = "BF_pb",
        inference_BF        = components[rownames(components) == "Bias", "inclusion_BF"]
      )
    }

    text <- BayesTools::interpret2(specification, method = .object_title(object))

  }


  return(text)
}

