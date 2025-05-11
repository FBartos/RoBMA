# the main functions
.fit_RoBMA_model       <- function(object, i, extend = FALSE){

  model              <- object[["models"]][[i]]
  priors             <- model[["priors"]]
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  add_info           <- object[["add_info"]]

  errors   <- NULL
  warnings <- NULL

  if(!fit_control[["silent"]]){
    cat(paste0("\nFitting model [", i, "]\n"))
  }

  # don't sample the complete null model
  if(!.is_model_constant(priors)){

    if(.is_model_multivariate(model)){
      object[["data"]] <- .order_data.mv(object[["data"]], .is_model_regression(model))
    }


    # deal with regression vs basic models
    if(.is_model_regression(model)){
      data_outcome       <- object[["data"]][["outcome"]]
      fit_priors         <- priors[names(priors) != "terms"]
      formula_list       <- .generate_model_formula_list(object[["formula"]])
      formula_data_list  <- .generate_model_formula_data_list(object[["data"]])
      formula_prior_list <- .generate_model_formula_prior_list(priors)
    }else if(inherits(model, "RoBMA.model")){
      data_outcome       <- object[["data"]]
      fit_priors         <- priors
      formula_list       <- NULL
      formula_data_list  <- NULL
      formula_prior_list <- NULL
    }

    # deal with multivariate vs univariate models
    if(.is_model_multivariate(model)){
      # generate the model syntax
      model_syntax <- .generate_model_syntax.mv(
        priors           = fit_priors,
        effect_direction = add_info[["effect_direction"]],
        prior_scale      = add_info[["prior_scale"]],
        effect_measure   = add_info[["effect_measure"]],
        data             = data_outcome,
        regression       = .is_model_regression(model)
        )

      # remove unnecessary objects from data to mitigate warnings
      fit_data     <- .fit_data.mv(
        data             = data_outcome,
        priors           = fit_priors,
        effect_direction = add_info[["effect_direction"]],
        prior_scale      = add_info[["prior_scale"]]
      )

    }else{
      # generate the model syntax
      model_syntax <- .generate_model_syntax(
        priors           = fit_priors,
        effect_direction = add_info[["effect_direction"]],
        prior_scale      = add_info[["prior_scale"]],
        effect_measure   = add_info[["effect_measure"]],
        weighted         = attr(model, "weighted"),
        regression       = .is_model_regression(model)
      )

      # remove unnecessary objects from data to mitigate warnings
      fit_data     <- .fit_data(
        data             = data_outcome,
        priors           = fit_priors,
        effect_direction = add_info[["effect_direction"]],
        prior_scale      = add_info[["prior_scale"]],
        weighted         = attr(model, "weighted"),
        weighted_type    = attr(model, "weighted_type")
      )
    }

    # fit the model
    if(!extend || length(model[["fit"]]) == 0){

      fit <- BayesTools::JAGS_fit(
        model_syntax       = model_syntax,
        data               = fit_data,
        prior_list         = fit_priors,
        formula_list       = formula_list,
        formula_data_list  = formula_data_list,
        formula_prior_list = formula_prior_list,
        chains             = fit_control[["chains"]],
        adapt              = fit_control[["adapt"]],
        burnin             = fit_control[["burnin"]],
        sample             = fit_control[["sample"]],
        thin               = fit_control[["thin"]],
        autofit            = fit_control[["autofit"]],
        autofit_control    = autofit_control,
        parallel           = fit_control[["parallel"]],
        cores              = fit_control[["cores"]],
        silent             = fit_control[["silent"]],
        seed               = fit_control[["seed"]],
        required_packages  = "RoBMA"
      )

    }else{

      fit <- BayesTools::JAGS_extend(
        fit                = model[["fit"]],
        autofit_control    = autofit_control,
        parallel           = fit_control[["parallel"]],
        cores              = fit_control[["cores"]],
        silent             = fit_control[["silent"]],
        seed               = fit_control[["seed"]]
      )

    }


    # assess the model fit and deal with errors
    if(inherits(fit, "error")){

      if(grepl("Unknown function", fit$message))
        stop("The RoBMA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.")

      fit                     <- list()
      attr(fit, "prior_list") <- fit_priors

      converged      <- FALSE
      has_posterior  <- FALSE
      errors         <- c(errors, fit$message)

      # deal with failed models
      marglik        <- list()
      marglik$logml  <- NA
      class(marglik) <- "bridge"

    }else{

      has_posterior <- TRUE
      check_fit     <- BayesTools::JAGS_check_convergence(
        fit          = fit,
        prior_list   = attr(fit, "prior_list"),
        max_Rhat     = convergence_checks[["max_Rhat"]],
        min_ESS      = convergence_checks[["min_ESS"]],
        max_error    = convergence_checks[["max_error"]],
        max_SD_error = convergence_checks[["max_SD_error"]]
      )
      warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
      if(convergence_checks[["remove_failed"]] && !check_fit){
        converged <- FALSE
      }else{
        converged <- TRUE
      }

    }

    # compute marginal likelihood
    if(length(fit) != 0){

      marglik <- BayesTools::JAGS_bridgesampling(
        fit                = fit,
        data               = fit_data,
        prior_list         = fit_priors,
        formula_list       = formula_list,
        formula_data_list  = formula_data_list,
        formula_prior_list = formula_prior_list,
        log_posterior      = if(.is_model_multivariate(model)) .marglik_function.mv  else .marglik_function,
        maxiter            = 50000,
        silent             = fit_control[["silent"]],
        priors             = priors,
        effect_direction   = add_info[["effect_direction"]],
        prior_scale        = add_info[["prior_scale"]],
        effect_measure     = add_info[["effect_measure"]]
      )

      # deal with failed marginal likelihoods
      if(inherits(marglik, "error")){

        error_message  <- marglik$message
        errors         <- c(errors, error_message)
        converged      <- FALSE
        marglik        <- list()
        marglik$logml  <- NA
        class(marglik) <- "bridge"
        attr(marglik, "errors") <- error_message

      }else{

        # forward warnings if present
        warnings <- c(warnings, attr(marglik, "warnings"))

      }
    }


  }else{

    # deal with regression vs basic models
    if(.is_model_regression(model)){

      # check that all terms but intercept are spikes at zero
      if(any(sapply(priors[["terms"]][names(priors[["terms"]]) != "intercept"], function(prior) prior$parameters[["location"]] != 0)))
        stop("All constant model can include only non-zero intercept parameter.")

      data_outcome    <- object[["data"]][["outcome"]]
      const_location  <- priors[["terms"]][["intercept"]]$parameters[["location"]]
      fit_priors      <- c(priors[names(priors) != "terms"], priors[["terms"]])
      names(fit_priors)[names(fit_priors) %in% names(priors[["terms"]])] <- paste0("mu_", names(fit_priors)[names(fit_priors) %in% names(priors[["terms"]])])
      fit_priors      <- .add_priors_levels(fit_priors, object[["data"]][["predictors"]])
    }else if(inherits(model, "RoBMA.model")){
      data_outcome    <- object[["data"]]
      fit_priors      <- priors
      const_location  <- priors$mu$parameters[["location"]]
    }

    if(fit_priors[["tau"]]$parameters[["location"]] != 0)
      stop("All constant model cannot include non zero heterogeneity parameter.")

    fit_data                <- .fit_data(
      data             = data_outcome,
      priors           = priors,
      effect_direction = add_info[["effect_direction"]],
      prior_scale      = add_info[["prior_scale"]],
      weighted         = attr(model, "weighted"),
      weighted_type    = attr(model, "weighted_type")
    )
    converged               <- TRUE
    has_posterior           <- FALSE
    fit                     <- list()
    attr(fit, "prior_list") <- fit_priors
    class(fit)              <- "null_model"
    marglik                 <- list()

    # weighted vs unweighted models
    if(attr(model, "weighted")){
      marglik$logml <- sum(stats::dnorm(fit_data[["y"]], const_location, fit_data[["se"]], log = TRUE) * fit_data[["weight"]])
    }else{
      marglik$logml <- sum(stats::dnorm(fit_data[["y"]], const_location, fit_data[["se"]], log = TRUE))
    }

    class(marglik) <- "bridge"

  }

  # add model summaries
  if(has_posterior){
    fit_summary   <- suppressMessages(BayesTools::runjags_estimates_table(fit = fit, warnings = warnings, transform_factors = TRUE, formula_prefix = FALSE))
    if(add_info[["prior_scale"]] != "y"){
      fit_summaries <- .runjags_summary_list(fit, attr(fit, "prior_list"), add_info[["prior_scale"]], warnings)
    }else{
      fit_summaries <- NULL
    }
  }else{
    fit_summary    <- BayesTools::runjags_estimates_empty_table()
    if(add_info[["prior_scale"]] != "y"){
      fit_summaries  <- list(
        "d"     = BayesTools::runjags_estimates_empty_table(),
        "r"     = BayesTools::runjags_estimates_empty_table(),
        "logOR" = BayesTools::runjags_estimates_empty_table(),
        "z"     = BayesTools::runjags_estimates_empty_table()
      )
    }else{
      fit_summaries <- NULL
    }
  }


  # add results
  model$fit           <- fit
  model$fit_summary   <- fit_summary
  model$fit_summaries <- fit_summaries
  model$marglik       <- marglik
  model$errors        <- errors
  model$warnings      <- warnings
  model$converged     <- converged
  model$has_posterior <- has_posterior
  model$output_scale  <- add_info[["prior_scale"]]
  model$prior_scale   <- add_info[["prior_scale"]]

  return(model)
}
.fit_RoBMA_model_ss    <- function(object,    extend = FALSE){

  model              <- object[["model"]]
  priors             <- object[["model"]][["priors"]]
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  add_info           <- object[["add_info"]]

  errors   <- NULL
  warnings <- NULL

  # deal with regression vs basic models
  if(.is_model_regression(model)){
    data_outcome       <- object[["data"]][["outcome"]]
    fit_priors         <- priors[names(priors) != "terms"]
    formula_list       <- .generate_model_formula_list(object[["formula"]])
    formula_data_list  <- .generate_model_formula_data_list(object[["data"]])
    formula_prior_list <- .generate_model_formula_prior_list(priors)
  }else{
    data_outcome       <- object[["data"]]
    fit_priors         <- priors
    formula_list       <- NULL
    formula_data_list  <- NULL
    formula_prior_list <- NULL
  }

  # generate the model syntax
  model_syntax <- .generate_model_syntax_ss(
    priors           = priors,
    effect_direction = add_info[["effect_direction"]],
    prior_scale      = add_info[["prior_scale"]],
    effect_measure   = add_info[["effect_measure"]],
    weighted         = attr(model, "weighted"),
    regression       = .is_model_regression(model),
    multivariate     = .is_model_multivariate(model)
  )

  # remove unnecessary objects from data to mitigate warnings
  fit_data     <- .fit_data_ss(
    data             = data_outcome,
    priors           = priors,
    effect_direction = add_info[["effect_direction"]],
    prior_scale      = add_info[["prior_scale"]],
    weighted         = attr(model, "weighted"),
    weighted_type    = attr(model, "weighted_type"),
    multivariate     = .is_model_multivariate(model)
  )

  # fit the model
  if(!extend || length(model[["fit"]]) == 0){

    fit <- BayesTools::JAGS_fit(
      model_syntax          = model_syntax,
      data                  = fit_data,
      prior_list            = fit_priors,
      formula_list          = formula_list,
      formula_data_list     = formula_data_list,
      formula_prior_list    = formula_prior_list,
      chains                = fit_control[["chains"]],
      adapt                 = fit_control[["adapt"]],
      burnin                = fit_control[["burnin"]],
      sample                = fit_control[["sample"]],
      thin                  = fit_control[["thin"]],
      autofit               = fit_control[["autofit"]],
      autofit_control       = autofit_control,
      parallel              = fit_control[["parallel"]],
      cores                 = fit_control[["cores"]],
      silent                = fit_control[["silent"]],
      seed                  = fit_control[["seed"]],
      required_packages     = "RoBMA"
    )

  }else{

    fit <- BayesTools::JAGS_extend(
      fit                = model[["fit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]]
    )

  }


  # assess the model fit and deal with errors
  if(inherits(fit, "error")){

    if(grepl("Unknown function", fit$message))
      stop("The RoBMA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.")

    fit                     <- list()
    attr(fit, "prior_list") <- fit_priors

    converged      <- FALSE
    has_posterior  <- FALSE
    errors         <- c(errors, fit$message)

    # deal with failed models
    marglik        <- list()
    marglik$logml  <- NA
    class(marglik) <- "bridge"

  }else{

    has_posterior <- TRUE
    check_fit     <- BayesTools::JAGS_check_convergence(
      fit          = fit,
      prior_list   = attr(fit, "prior_list"),
      max_Rhat     = convergence_checks[["max_Rhat"]],
      min_ESS      = convergence_checks[["min_ESS"]],
      max_error    = convergence_checks[["max_error"]],
      max_SD_error = convergence_checks[["max_SD_error"]]
    )
    warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
    if(convergence_checks[["remove_failed"]] && !check_fit){
      converged <- FALSE
    }else{
      converged <- TRUE
    }

  }


  # add results
  model$fit           <- fit
  model$errors        <- errors
  model$warnings      <- warnings
  model$converged     <- converged
  model$has_posterior <- has_posterior
  model$output_scale  <- add_info[["prior_scale"]]
  model$prior_scale   <- add_info[["prior_scale"]]

  return(model)
}
.fit_BiBMA_model       <- function(object, i, extend = FALSE){

  model              <- object[["models"]][[i]]
  priors             <- model[["priors"]]
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  add_info           <- object[["add_info"]]

  errors   <- NULL
  warnings <- NULL

  if(!fit_control[["silent"]]){
    cat(paste0("\nFitting model [", i, "]\n"))
  }

  ### the model is never constant as there needs to be a prior for pi
  # deal with regression vs basic models
  if(.is_model_regression(model)){
    data_outcome       <- object[["data"]][["outcome"]]
    fit_priors         <- priors[names(priors) != "terms"]
    formula_list       <- .generate_model_formula_list(object[["formula"]])
    formula_data_list  <- .generate_model_formula_data_list(object[["data"]])
    formula_prior_list <- .generate_model_formula_prior_list(priors)
  }else if(inherits(model, "BiBMA.model")){
    data_outcome       <- object[["data"]]
    fit_priors         <- priors
    formula_list       <- NULL
    formula_data_list  <- NULL
    formula_prior_list <- NULL
  }

  model_syntax <- .generate_model_syntax.bi(
    priors           = fit_priors,
    random           = attr(model, "random"),
    weighted         = attr(model, "weighted"),
    regression       = .is_model_regression(model),
    multivariate     = .is_model_multivariate(model)
  )

  # remove unnecessary objects from data to mitigate warnings
  fit_data     <- .fit_data.bi(
    data             = data_outcome,
    weighted         = attr(model, "weighted"),
    weighted_type    = attr(model, "weighted_type"),
    multivariate     = .is_model_multivariate(model)
  )

  # fit the model
  if(!extend || length(model[["fit"]]) == 0){

    fit <- BayesTools::JAGS_fit(
      model_syntax       = model_syntax,
      data               = fit_data,
      prior_list         = fit_priors,
      formula_list       = formula_list,
      formula_data_list  = formula_data_list,
      formula_prior_list = formula_prior_list,
      chains             = fit_control[["chains"]],
      adapt              = fit_control[["adapt"]],
      burnin             = fit_control[["burnin"]],
      sample             = fit_control[["sample"]],
      thin               = fit_control[["thin"]],
      autofit            = fit_control[["autofit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]],
      required_packages  = "RoBMA"
    )

  }else{

    fit <- BayesTools::JAGS_extend(
      fit             = model[["fit"]],
      autofit_control = autofit_control,
      parallel        = fit_control[["parallel"]],
      cores           = fit_control[["cores"]],
      silent          = fit_control[["silent"]],
      seed            = fit_control[["seed"]]
    )

  }

  # assess the model fit and deal with errors
  if(inherits(fit, "error")){

    if(grepl("Unknown function", fit$message))
      stop("The RoBMA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.")

    fit                     <- list()
    attr(fit, "prior_list") <- fit_priors

    converged      <- FALSE
    has_posterior  <- FALSE
    errors         <- c(errors, fit$message)

    # deal with failed models
    marglik        <- list()
    marglik$logml  <- NA
    class(marglik) <- "bridge"

  }else{

    has_posterior <- TRUE
    check_fit     <- BayesTools::JAGS_check_convergence(
      fit          = fit,
      prior_list   = attr(fit, "prior_list"),
      max_Rhat     = convergence_checks[["max_Rhat"]],
      min_ESS      = convergence_checks[["min_ESS"]],
      max_error    = convergence_checks[["max_error"]],
      max_SD_error = convergence_checks[["max_SD_error"]]
    )
    warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
    if(convergence_checks[["remove_failed"]] && !check_fit){
      converged <- FALSE
    }else{
      converged <- TRUE
    }

  }

  # compute marginal likelihood
  if(length(fit) != 0){

    marglik <- BayesTools::JAGS_bridgesampling(
      fit                = fit,
      data               = fit_data,
      prior_list         = fit_priors,
      formula_list       = formula_list,
      formula_data_list  = formula_data_list,
      formula_prior_list = formula_prior_list,
      log_posterior      = .marglik_function.bi,
      maxiter            = 50000,
      silent             = fit_control[["silent"]],
      priors             = priors,
      random             = attr(model, "random")
    )

    # deal with failed marginal likelihoods
    if(inherits(marglik, "error")){

      errors         <- c(errors, marglik$message)
      converged      <- FALSE
      marglik        <- list()
      marglik$logml  <- NA
      class(marglik) <- "bridge"

    }else{

      # forward warnings if present
      warnings <- c(warnings, attr(marglik, "warnings"))

    }
  }

  # add model summaries
  if(has_posterior){
    fit_summary   <- BayesTools::runjags_estimates_table(
      fit                = fit,
      warnings           = warnings,
      transform_factors  = TRUE,
      formula_prefix     = FALSE,
      remove_parameters  = c("pi", if(attr(model, "random")) "gamma")
    )
    fit_summaries <- .runjags_summary_list(
      fit               = fit,
      priors            = attr(fit, "prior_list"),
      prior_scale       = add_info[["prior_scale"]],
      warnings          = warnings,
      remove_parameters = c("pi", if(attr(model, "random")) "gamma")
    )
  }else{
    fit_summary    <- BayesTools::runjags_estimates_empty_table()
    fit_summaries  <- list(
      "d"     = BayesTools::runjags_estimates_empty_table(),
      "r"     = BayesTools::runjags_estimates_empty_table(),
      "logOR" = BayesTools::runjags_estimates_empty_table(),
      "z"     = BayesTools::runjags_estimates_empty_table()
    )
  }


  # add results
  model$fit           <- fit
  model$fit_summary   <- fit_summary
  model$fit_summaries <- fit_summaries
  model$marglik       <- marglik
  model$errors        <- errors
  model$warnings      <- warnings
  model$converged     <- converged
  model$has_posterior <- has_posterior
  model$output_scale  <- add_info[["prior_scale"]]
  model$prior_scale   <- add_info[["prior_scale"]]

  return(model)
}
.fit_BiBMA_model_ss    <- function(object,    extend = FALSE){

  model              <- object[["model"]]
  priors             <- object[["model"]][["priors"]]
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  add_info           <- object[["add_info"]]

  errors   <- NULL
  warnings <- NULL

  # deal with regression vs basic models
  if(.is_model_regression(model)){
    data_outcome       <- object[["data"]][["outcome"]]
    fit_priors         <- priors[names(priors) != "terms"]
    formula_list       <- .generate_model_formula_list(object[["formula"]])
    formula_data_list  <- .generate_model_formula_data_list(object[["data"]])
    formula_prior_list <- .generate_model_formula_prior_list(priors)
  }else if(inherits(model, "BiBMA.model_ss")){
    data_outcome       <- object[["data"]]
    fit_priors         <- priors
    formula_list       <- NULL
    formula_data_list  <- NULL
    formula_prior_list <- NULL
  }

  model_syntax <- .generate_model_syntax.bi(
    priors           = fit_priors,
    random           = attr(model, "random"),
    weighted         = attr(model, "weighted"),
    regression       = .is_model_regression(model),
    multivariate     = .is_model_multivariate(model)
  )

  # remove unnecessary objects from data to mitigate warnings
  fit_data     <- .fit_data.bi(
    data             = data_outcome,
    weighted         = attr(model, "weighted"),
    weighted_type    = attr(model, "weighted_type"),
    multivariate     = .is_model_multivariate(model)
  )


  # fit the model
  if(!extend || length(model[["fit"]]) == 0){

    fit <- BayesTools::JAGS_fit(
      model_syntax          = model_syntax,
      data                  = fit_data,
      prior_list            = fit_priors,
      formula_list          = formula_list,
      formula_data_list     = formula_data_list,
      formula_prior_list    = formula_prior_list,
      chains                = fit_control[["chains"]],
      adapt                 = fit_control[["adapt"]],
      burnin                = fit_control[["burnin"]],
      sample                = fit_control[["sample"]],
      thin                  = fit_control[["thin"]],
      autofit               = fit_control[["autofit"]],
      autofit_control       = autofit_control,
      parallel              = fit_control[["parallel"]],
      cores                 = fit_control[["cores"]],
      silent                = fit_control[["silent"]],
      seed                  = fit_control[["seed"]],
      required_packages     = "RoBMA"
    )

  }else{

    fit <- BayesTools::JAGS_extend(
      fit                = model[["fit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]]
    )

  }


  # assess the model fit and deal with errors
  if(inherits(fit, "error")){

    if(grepl("Unknown function", fit$message))
      stop("The RoBMA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.")

    fit                     <- list()
    attr(fit, "prior_list") <- fit_priors

    converged      <- FALSE
    has_posterior  <- FALSE
    errors         <- c(errors, fit$message)

    # deal with failed models
    marglik        <- list()
    marglik$logml  <- NA
    class(marglik) <- "bridge"

  }else{

    has_posterior <- TRUE
    check_fit     <- BayesTools::JAGS_check_convergence(
      fit          = fit,
      prior_list   = attr(fit, "prior_list"),
      max_Rhat     = convergence_checks[["max_Rhat"]],
      min_ESS      = convergence_checks[["min_ESS"]],
      max_error    = convergence_checks[["max_error"]],
      max_SD_error = convergence_checks[["max_SD_error"]]
    )
    warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
    if(convergence_checks[["remove_failed"]] && !check_fit){
      converged <- FALSE
    }else{
      converged <- TRUE
    }

  }

  # add results
  model$fit           <- fit
  model$errors        <- errors
  model$warnings      <- warnings
  model$converged     <- converged
  model$has_posterior <- has_posterior
  model$output_scale  <- add_info[["prior_scale"]]
  model$prior_scale   <- add_info[["prior_scale"]]

  return(model)
}

# tools
.fit_data                 <- function(data, priors, effect_direction, prior_scale, weighted, weighted_type){

  # unlist the data.frame
  original_measure <- attr(data, "original_measure")
  effect_measure   <- attr(data, "effect_measure")

  fit_data <- list()
  # change the effect size direction (important for one-sided selection and PET/PEESE)
  if(effect_direction == "negative"){
    fit_data$y <- - data[,"y"]
  }else{
    fit_data$y <- data[,"y"]
  }
  fit_data$se <- data[,"se"]
  fit_data$K  <- length(data[["y"]])

  # add critical y-values
  if(!is.null(priors[["omega"]])){
    fit_data$crit_y  <- t(.get_cutoffs(fit_data[["y"]], fit_data[["se"]], priors[["omega"]], original_measure, effect_measure))
  }

  # add weights proportional to the number of estimates from a study
  if(weighted){
    fit_data$weight <- .get_id_weights(data, weighted_type)
  }

  return(fit_data)
}
.fit_data.mv              <- function(data, priors, effect_direction, prior_scale){

  # unlist the data.frame
  original_measure <- attr(data, "original_measure")
  effect_measure   <- attr(data, "effect_measure")

  fit_data <- list()
  data     <- data[order(data[,"study_ids"]),]

  ### deal with the univariate data
  # change the effect size direction (important for one-sided selection and PET/PEESE)
  if(any(is.na(data[,"study_ids"]))){
    if(effect_direction == "negative"){
      fit_data$y <- - data[is.na(data[,"study_ids"]),"y"]
    }else{
      fit_data$y <- data[is.na(data[,"study_ids"]),"y"]
    }
    fit_data$se <- data[is.na(data[,"study_ids"]),"se"]
    fit_data$K  <- length(fit_data[["y"]])

    # add critical y-values
    if(!is.null(priors[["omega"]])){
      fit_data$crit_y  <- t(.get_cutoffs(fit_data[["y"]], fit_data[["se"]], priors[["omega"]], original_measure[is.na(data[,"study_ids"])], effect_measure))
    }
  }


  ### add the multivariate part
  if(effect_direction == "negative"){
    fit_data$y_v <- - data[!is.na(data[,"study_ids"]),"y"]
  }else{
    fit_data$y_v <- data[!is.na(data[,"study_ids"]),"y"]
  }
  fit_data$se2_v <- data[!is.na(data[,"study_ids"]),"se"]^2
  fit_data$K_v   <- length(fit_data[["y_v"]])

  # add critical y-values
  if(!is.null(priors[["omega"]])){
    fit_data$crit_y_v  <- t(.get_cutoffs(fit_data[["y_v"]], data[!is.na(data[,"study_ids"]),"se"], priors[["omega"]], original_measure[!is.na(data[,"study_ids"])], effect_measure))
  }else if(!is.null(priors[["PET"]])){
    fit_data$se_v  <- data[!is.na(data[,"study_ids"]),"se"]
  }

  fit_data$indx_v <- c((1:fit_data[["K_v"]])[!duplicated(data[!is.na(data[,"study_ids"]),"study_ids"])][-1] - 1, fit_data[["K_v"]])

  return(fit_data)
}
.fit_data.bi              <- function(data, weighted, weighted_type, multivariate){

  # unlist the data frame
  fit_data <- list()
  fit_data$x1 <- data[["x1"]]
  fit_data$x2 <- data[["x2"]]
  fit_data$n1 <- data[["n1"]]
  fit_data$n2 <- data[["n2"]]
  fit_data$K  <- nrow(data)

  # add weights proportional to the number of estimates from a study
  if(weighted){
    fit_data$weight <- .get_id_weights(data, weighted_type)
  }

  if(multivariate){
    # rewritting ids as previous formating used NAs for unique ones
    data[,"study_ids"][is.na(data[,"study_ids"])] <- (max(data[,"study_ids"], na.rm = TRUE) + 1):(max(data[,"study_ids"], na.rm = TRUE) + sum(is.na(data[,"study_ids"])))
    fit_data$study_ids <- data[,"study_ids"]
  }

  return(fit_data)
}
.fit_data_ss              <- function(data, priors, effect_direction, prior_scale, weighted, weighted_type, multivariate){

  # unlist the data.frame
  original_measure <- attr(data, "original_measure")
  effect_measure   <- attr(data, "effect_measure")

  fit_data <- list()
  # change the effect size direction (important for one-sided selection and PET/PEESE)
  if(effect_direction == "negative"){
    fit_data$y <- - data[,"y"]
  }else{
    fit_data$y <- data[,"y"]
  }
  fit_data$se <- data[,"se"]
  fit_data$K  <- length(data[["y"]])

  # add weights proportional to the number of estimates from a study
  if(weighted){
    fit_data$weight <- .get_id_weights(data, weighted_type)
  }

  if(any(sapply(priors[["bias"]], is.prior.weightfunction))){

    # create the weightfunction mapping for effect size thresholds
    steps  <- BayesTools::weightfunctions_mapping(priors[["bias"]][sapply(priors[["bias"]], is.prior.weightfunction)], cuts_only = TRUE, one_sided = TRUE)
    steps  <- rev(steps)[c(-1, -length(steps))]
    crit_y <- .get_cutoffs(fit_data[["y"]], fit_data[["se"]], list(distribution = "one.sided", parameters = list(steps = steps)), original_measure, effect_measure)

    # create the weightfunction mapping to weights (transform all weight functions to one-sided)
    crit_y_mapping <- matrix(0, nrow = length(priors[["bias"]]), ncol = ncol(crit_y))
    crit_y_mapping_max <- rep(0, length(priors[["bias"]]))
    for(i in seq_along(priors[["bias"]])){
      if(is.prior.weightfunction(priors[["bias"]][[i]])){
        ### the following subsetting allows us "merge" steps with equal weights due to the construction
        # specify indexes of the relevant steps
        this_steps <- .get_one_sided_cuts(priors[["bias"]][[i]])
        crit_y_mapping[i,1:length(this_steps)] <- which(steps %in% this_steps)
        crit_y_mapping_max[i] <- length(this_steps)
      }
    }

    fit_data$crit_y             <- t(crit_y)
    fit_data$crit_y_mapping     <- t(crit_y_mapping)
    fit_data$crit_y_mapping_max <- crit_y_mapping_max
  }

  if(multivariate){
    # rewritting ids as previous formating used NAs for unique ones
    data[,"study_ids"][is.na(data[,"study_ids"])] <- (max(data[,"study_ids"], na.rm = TRUE) + 1):(max(data[,"study_ids"], na.rm = TRUE) + sum(is.na(data[,"study_ids"])))
    fit_data$study_ids <- data[,"study_ids"]
  }

  return(fit_data)
}
.order_data.mv            <- function(data, regression){
  # prepares data in a better order for the subsequent vectorization of multivariate distributions

  if(regression){
    ids <- data[["outcome"]]$study_ids
  }else{
    ids <- data$study_ids
  }

  # first independent and then dependent estimates
  ordering <- order(ifelse(is.na(ids), -1, ids))

  # re-order the data set and predictors
  if(regression){
    data[["outcome"]] <- data[["outcome"]][ordering,]
    for(i in seq_along(data[["predictors"]])){
      data[["predictors"]][[i]] <- data[["predictors"]][[i]][ordering]
    }
  }else{
    data <- data[ordering,]
  }

  return(data)
}
.generate_model_syntax    <- function(priors, effect_direction, prior_scale, effect_measure, weighted, regression){

  model_syntax <- "model{\n"

  ### prior transformations
  # the precise transformation for heterogeneity is not used due the inability to re-scale large variances
  # instead, approximate linear scaling is employed in the same way as in metaBMA package
  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "mu[i]",  "mu_transformed[i]"))
    model_syntax <- paste0(model_syntax, "}\n")
  }else{
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "mu",  "mu_transformed"))
  }

  model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "tau", "tau_transformed"))
  if(!is.null(priors[["PET"]])){
    model_syntax <- paste0(model_syntax, paste0("PET_transformed = PET\n"))
  }else if(!is.null(priors[["PEESE"]])){
    # don't forget that the transformation is inverse for PEESE
    model_syntax <- paste0(model_syntax, .JAGS_scale(effect_measure, prior_scale, "PEESE", "PEESE_transformed"))
  }


  ### model
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")

  # marginalized random effects and the effect size
  prec <- "1 / ( pow(se[i],2) + pow(tau_transformed,2) )"

  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed[i]", "mu_transformed[i]")
  }else{
    eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed", "mu_transformed")
  }

  # add PET/PEESE
  if(!is.null(priors[["PET"]])){
    eff <- paste0("(", eff, " + PET_transformed * se[i])")
  }else if(!is.null(priors[["PEESE"]])){
    eff <- paste0("(", eff, " + PEESE_transformed * pow(se[i],2))")
  }

  # the observed data
  if(weighted){
    if(is.null(priors[["omega"]])){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm(",     eff, ",", prec, ", weight[i])\n")
    }else if(grepl("one.sided", priors[["omega"]]$distribution)){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwwnorm_1s(", eff, ",", prec, ", crit_y[,i], omega, weight[i]) \n")
    }else if(grepl("two.sided", priors[["omega"]]$distribution)){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwwnorm_2s(", eff, ",", prec, ", crit_y[,i], omega, weight[i]) \n")
    }
  }else{
    if(is.null(priors[["omega"]])){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dnorm(",     eff, ",", prec, " )\n")
    }else if(grepl("one.sided", priors[["omega"]]$distribution)){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_1s(", eff, ",", prec, ", crit_y[,i], omega) \n")
    }else if(grepl("two.sided", priors[["omega"]]$distribution)){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_2s(", eff, ",", prec, ", crit_y[,i], omega) \n")
    }
  }

  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.generate_model_syntax.mv <- function(priors, effect_direction, prior_scale, effect_measure, data, regression){

  model_syntax <- "model{\n"

  ### prior transformations
  # the precise transformation for heterogeneity is not used due the inability to re-scale large variances
  # instead, approximate linear scaling is employed in the same way as in metaBMA package
  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    model_syntax <- paste0(model_syntax, "for(i in 1:(K+K_v)){\n")
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "mu[i]",  "mu_transformed[i]"))
    model_syntax <- paste0(model_syntax, "}\n")
  }else{
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "mu",  "mu_transformed"))
  }
  model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "tau", "tau_transformed"))

  if(!is.null(priors[["PET"]])){
    model_syntax <- paste0(model_syntax, paste0("PET_transformed = PET\n"))
  }else if(!is.null(priors[["PEESE"]])){
    # don't forget that the transformation is inverse for PEESE
    model_syntax <- paste0(model_syntax, .JAGS_scale(effect_measure, prior_scale, "PEESE", "PEESE_transformed"))
  }


  ### deal with the univariate data
  if(any(is.na(data[,"study_ids"]))){

    # marginalized random effects and the effect size
    prec     <- "1 / ( pow(se[i],2) + pow(tau_transformed,2) )"

    # deal with mu as a vector or scalar based on whether it is regression or not
    if(regression){
      eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed[i]", "mu_transformed[i]")
    }else{
      eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed", "mu_transformed")
    }

    # add PET/PEESE
    if(!is.null(priors[["PET"]])){
      eff <- paste0("(", eff, " + PET_transformed * se[i])")
    }else if(!is.null(priors[["PEESE"]])){
      eff <- paste0("(", eff, " + PEESE_transformed * pow(se[i],2))")
    }

    # the observed data
    model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")
    if(is.null(priors[["omega"]])){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dnorm(",     eff, ",", prec, " )\n")
    }else if(grepl("one.sided", priors[["omega"]]$distribution)){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_1s(", eff, ",", prec, ", crit_y[,i], omega) \n")
    }else if(grepl("two.sided", priors[["omega"]]$distribution)){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_2s(", eff, ",", prec, ", crit_y[,i], omega) \n")
    }
    model_syntax <- paste0(model_syntax, "}\n")

  }


  ### deal with the multivariate data
  model_syntax <- paste0(model_syntax, paste0("tau_transformed2 = pow(tau_transformed, 2)\n"))

  # create the mean vector
  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    eff_v <- ifelse(effect_direction == "negative", "-1 * mu_transformed[K+i]", "mu_transformed[K+i]")
  }else{
    eff_v <- ifelse(effect_direction == "negative", "-1 * mu_transformed", "mu_transformed")
  }
  model_syntax <- paste0(model_syntax, "for(i in 1:K_v){\n")
  if(!is.null(priors[["PET"]])){
    eff_v <- paste0(eff_v, " + PET_transformed * se_v[i]")
  }else if(!is.null(priors[["PEESE"]])){
    eff_v <- paste0(eff_v, " + PEESE_transformed * se2_v[i]")
  }
  model_syntax <- paste0(model_syntax, paste0("  eff_v[i] = ", eff_v, "\n"))
  model_syntax <- paste0(model_syntax, "}\n")

  # the observed data
  if(is.null(priors[["omega"]])){
    model_syntax <- paste0(model_syntax, "y_v ~ dmnorm_v(eff_v, se2_v, tau_transformed2, rho, indx_v)\n")
  }else if(grepl("one.sided", priors[["omega"]]$distribution)){
    model_syntax <- paste0(model_syntax, "y_v ~ dwmnorm_1s_v(eff_v, se2_v, tau_transformed2, rho, crit_y_v, omega, indx_v) \n")
  }else if(grepl("two.sided", priors[["omega"]]$distribution)){
    model_syntax <- paste0(model_syntax, "y_v ~ dwmnorm_2s_v(eff_v, se2_v, tau_transformed2, rho, crit_y_v, omega, indx_v) \n")
  }

  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.generate_model_syntax.bi <- function(priors, random, weighted, regression, multivariate){

  model_syntax <- "model{\n"

  ### priors and model are always specified on log(OR) scale
  # no need to apply transformations here

  if(random && multivariate){
    model_syntax <- paste0(model_syntax, "tau_within  = tau * sqrt(rho)\n")
    model_syntax <- paste0(model_syntax, "tau_between = tau * sqrt(1-rho)\n")
  }

  ### model
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")
  # using non-central parameterization
  if(regression){
    eff <- "mu[i]"
  }else{
    eff <- "mu"
  }
  if(random && multivariate){
    eff <- paste0(eff, " + epsilon[study_ids[i]] * tau_within + gamma[i] * tau_between")
  }else if(random && !multivariate){
    eff <- paste0(eff, " + gamma[i] * tau")
  }else if(!random && multivariate){
    eff <- paste0(eff, " + epsilon[study_ids[i]] * tau")
  }

  # transform the parameters to the probability scale
  model_syntax <- paste0(model_syntax, "  logit(p1[i]) = logit(pi[i]) + 0.5 * (", eff, ")\n")
  model_syntax <- paste0(model_syntax, "  logit(p2[i]) = logit(pi[i]) - 0.5 * (", eff, ")\n")

  # the observed data
  if(weighted){
    stop("This is not enough to take fraction of the likelihood properly -- also the non-marginalized random effects need to be repeated among the same observations.")
    model_syntax <- paste0(model_syntax, "  x1[i] ~ dwbinom(p1[i], n1[i], weight[i])\n")
    model_syntax <- paste0(model_syntax, "  x2[i] ~ dwbinom(p2[i], n2[i], weight[i])\n")
  }else{
    model_syntax <- paste0(model_syntax, "  x1[i] ~ dbinom(p1[i], n1[i])\n")
    model_syntax <- paste0(model_syntax, "  x2[i] ~ dbinom(p2[i], n2[i])\n")
  }

  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.generate_model_syntax_ss <- function(priors, effect_direction, prior_scale, effect_measure, weighted, regression, multivariate){

  ### extract prior information
  if(is.prior.mixture(priors[["bias"]])){
    is_PET            <- sapply(priors[["bias"]], is.prior.PET)
    is_PEESE          <- sapply(priors[["bias"]], is.prior.PEESE)
    is_weightfunction <- sapply(priors[["bias"]], is.prior.weightfunction)
  }else{
    is_PET            <- is.prior.PET(priors[["bias"]])
    is_PEESE          <- is.prior.PEESE(priors[["bias"]])
    is_weightfunction <- is.prior.weightfunction(priors[["bias"]])
  }


  # create the model syntax
  model_syntax <- "model{\n"

  ### prior transformations
  # the precise transformation for heterogeneity is not used due the inability to re-scale large variances
  # instead, approximate linear scaling is employed in the same way as in metaBMA package
  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "mu[i]",  "mu_transformed[i]"))
    model_syntax <- paste0(model_syntax, "}\n")
  }else{
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "mu",  "mu_transformed"))
  }

  if(multivariate){
    model_syntax <- paste0(model_syntax, "tau_within  = tau * sqrt(rho)\n")
    model_syntax <- paste0(model_syntax, "tau_between = tau * sqrt(1-rho)\n")
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "tau", "tau_transformed"))
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "tau_within",  "tau_within_transformed"))
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "tau_between", "tau_between_transformed"))
  }else{
    model_syntax <- paste0(model_syntax, .JAGS_scale(prior_scale, effect_measure, "tau", "tau_transformed"))
  }

  if(any(is_PET)){
    model_syntax <- paste0(model_syntax, paste0("PET_transformed = PET\n"))
  }
  if(any(is_PEESE)){
    # don't forget that the transformation is inverse for PEESE
    model_syntax <- paste0(model_syntax, .JAGS_scale(effect_measure, prior_scale, "PEESE", "PEESE_transformed"))
  }

  # marginalized random effects and the effect size
  prec <- "1 / ( pow(se[i],2) + pow(tau_transformed,2) )"

  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed[i]", "mu_transformed[i]")
  }else{
    eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed", "mu_transformed")
  }

  ### model
  model_syntax <- paste0(model_syntax, "\nfor(i in 1:K){\n")

  # marginalized random effects and the effect size
  if(multivariate){
    tau2 <- "( pow(se[i],2) + pow(tau_between_transformed,2) )"
  }else{
    tau2 <- "( pow(se[i],2) + pow(tau_transformed,2) )"
  }

  # deal with mu as a vector or scalar based on whether it is regression or not
  if(regression){
    eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed[i]", "mu_transformed[i]")
  }else{
    eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed", "mu_transformed")
  }

  # add PET/PEESE
  if(any(is_PET)){
    eff <- paste0(eff, " + PET_transformed * se[i]")
  }
  if(any(is_PEESE)){
    eff <- paste0(eff, " + PEESE_transformed * pow(se[i],2)")
  }

  # add hierarchical
  if(multivariate){
    eff <- paste0(eff, " + gamma[study_ids[i]] * tau_within_transformed")
  }

  if(any(is_weightfunction)){
    if(weighted){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwwnorm_mix(", eff, ",", "sqrt", tau2, ", crit_y[,i], omega, crit_y_mapping[,bias_indicator], crit_y_mapping_max[bias_indicator], weight[i])\n")
    }else{
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_mix(", eff, ",", "sqrt", tau2, ", crit_y[,i], omega, crit_y_mapping[,bias_indicator], crit_y_mapping_max[bias_indicator])\n")
    }
  }else{
    if(weighted){
      model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm(", eff, ",", "1/", tau2, ", weight[i])\n")
    }else{
      model_syntax <- paste0(model_syntax, "  y[i] ~ dnorm(",  eff, ",", "1/", tau2, ")\n")
    }
  }

  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.marglik_function         <- function(parameters, data, priors, effect_direction, prior_scale, effect_measure){

  # extract parameters
  mu  <- parameters[["mu"]]
  tau <- parameters[["tau"]]
  if(!is.null(priors[["PET"]])){
    PET    <- parameters[["PET"]]
  }else if(!is.null(priors[["PEESE"]])){
    PEESE  <- parameters[["PEESE"]]
  }else if(!is.null(priors[["omega"]])){
    omega  <- parameters[["omega"]]
  }

  ### re-scale parameters
  if(prior_scale != effect_measure){
    mu_transformed  <- do.call(
      .get_scale(prior_scale, effect_measure),
      args = list(mu)
    )
    tau_transformed <- do.call(
      .get_scale(prior_scale, effect_measure),
      args = list(tau)
    )
    if(!is.null(priors[["PET"]])){
      PET_transformed   <- PET
    }else if(!is.null(priors[["PEESE"]])){
      # don't forget that the transformation is inverse for PEESE
      PEESE_transformed <- do.call(
        .get_scale(effect_measure, prior_scale),
        args = list(PEESE)
      )
    }
  }else{
    mu_transformed  <- mu
    tau_transformed <- tau
    if(!is.null(priors[["PET"]])){
      PET_transformed   <- PET
    }else if(!is.null(priors[["PEESE"]])){
      PEESE_transformed <- PEESE
    }
  }

  ### model
  # marginalized random effects and the effect size
  pop_sd  <- sqrt(data[["se"]]^2 + tau_transformed^2)

  # add PET/PEESE
  eff <- ifelse(effect_direction == "negative", -1, 1) * mu_transformed
  if(!is.null(priors[["PET"]])){
    eff <- eff + PET_transformed   * data[["se"]]
  }else if(!is.null(priors[["PEESE"]])){
    eff <- eff + PEESE_transformed * data[["se"]]^2
  }

  ### compute the marginal log_likelihood
  log_lik <- 0

  # the individual studies
  if(!is.null(data[["weight"]])){
    if(is.null(priors[["omega"]])){
      log_lik <- log_lik + sum(stats::dnorm(data[["y"]], mean = eff, sd = pop_sd, log = TRUE) * data[["weight"]])
    }else if(priors[["omega"]]$distribution == "one.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = t(data[["crit_y"]]), type = "one.sided", log = TRUE) * data[["weight"]])
    }else if(priors[["omega"]]$distribution == "two.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = t(data[["crit_y"]]), type = "two.sided", log = TRUE) * data[["weight"]])
    }
  }else{
    if(is.null(priors[["omega"]])){
      log_lik <- log_lik + sum(stats::dnorm(data[["y"]], mean = eff, sd = pop_sd, log = TRUE))
    }else if(priors[["omega"]]$distribution == "one.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = t(data[["crit_y"]]), type = "one.sided", log = TRUE))
    }else if(priors[["omega"]]$distribution == "two.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = t(data[["crit_y"]]), type = "two.sided", log = TRUE))
    }
  }


  return(log_lik)
}
.marglik_function.mv      <- function(parameters, data, priors, effect_direction, prior_scale, effect_measure){

  # extract parameters
  mu  <- parameters[["mu"]]
  tau <- parameters[["tau"]]
  if(!is.null(priors[["PET"]])){
    PET    <- parameters[["PET"]]
  }else if(!is.null(priors[["PEESE"]])){
    PEESE  <- parameters[["PEESE"]]
  }else if(!is.null(priors[["omega"]])){
    omega  <- parameters[["omega"]]
  }
  if(!is.null(priors[["rho"]])){
    rho <- parameters[["rho"]]
  }

  ### re-scale parameters
  if(prior_scale != effect_measure){
    mu_transformed  <- do.call(
      .get_scale(prior_scale, effect_measure),
      args = list(mu)
    )
    tau_transformed <- do.call(
      .get_scale(prior_scale, effect_measure),
      args = list(tau)
    )
    if(!is.null(priors[["PET"]])){
      PET_transformed   <- PET
    }else if(!is.null(priors[["PEESE"]])){
      # don't forget that the transformation is inverse for PEESE
      PEESE_transformed <- do.call(
        .get_scale(effect_measure, prior_scale),
        args = list(PEESE)
      )
    }
  }else{
    mu_transformed  <- mu
    tau_transformed <- tau
    if(!is.null(priors[["PET"]])){
      PET_transformed   <- PET
    }else if(!is.null(priors[["PEESE"]])){
      PEESE_transformed <- PEESE
    }
  }

  ### model


  ### compute the marginal log_likelihood
  log_lik <- 0

  ### the independent studies
  if(!is.null(data[["y"]])){

    # marginalized random effects and the effect size
    pop_sd  <- sqrt(data[["se"]]^2 + tau_transformed^2)

    # add PET/PEESE
    eff <- ifelse(effect_direction == "negative", -1, 1) * mu_transformed
    if(!is.null(priors[["PET"]])){
      eff <- eff + PET_transformed   * data[["se"]]
    }else if(!is.null(priors[["PEESE"]])){
      eff <- eff + PEESE_transformed * data[["se"]]^2
    }

    if(is.null(priors[["omega"]])){
      log_lik <- log_lik + sum(stats::dnorm(data[["y"]], mean = eff, sd = pop_sd, log = TRUE))
    }else if(priors[["omega"]]$distribution == "one.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = t(data[["crit_y"]]), type = "one.sided", log = TRUE))
    }else if(priors[["omega"]]$distribution == "two.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = t(data[["crit_y"]]), type = "two.sided", log = TRUE))
    }
  }


  ### the dependent studies
  # construct the mean vector & add PET/PEESE
  temp_eff <- ifelse(effect_direction == "negative", -1, 1) * mu_transformed
  if(!is.null(priors[["PET"]])){
    temp_eff <- temp_eff + PET_transformed   * data[["se_v"]]
  }else if(!is.null(priors[["PEESE"]])){
    temp_eff <- temp_eff + PEESE_transformed * data[["se2_v"]]
  }else{
    temp_eff <- rep(temp_eff, data[["K_v"]])
  }

  # the observed data
  if(is.null(priors[["omega"]])){
    log_lik <- log_lik + .dwmnorm_v_fast(data[["y_v"]], mean_v = temp_eff, se2_v = data[["se2_v"]], tau2 = tau_transformed^2, rho = rho, type = "none", indx_v = data[["indx_v"]], log = TRUE)
  }else if(grepl("one.sided", priors[["omega"]]$distribution)){
    log_lik <- log_lik + .dwmnorm_v_fast(data[["y_v"]], mean_v = temp_eff, se2_v = data[["se2_v"]], tau2 = tau_transformed^2, rho = rho, crit_x_v = data[["crit_y_v"]], omega = omega, indx_v = data[["indx_v"]], type = "one.sided", log = TRUE)
  }else if(grepl("two.sided", priors[["omega"]]$distribution)){
    log_lik <- log_lik + .dwmnorm_v_fast(data[["y_v"]], mean_v = temp_eff, se2_v = data[["se2_v"]], tau2 = tau_transformed^2, rho = rho, crit_x_v = data[["crit_y_v"]], omega = omega, indx_v = data[["indx_v"]], type = "two.sided", log = TRUE)
  }

  return(log_lik)
}
.marglik_function.bi      <- function(parameters, data, priors, random){

  # extract parameters
  mu    <- parameters[["mu"]]
  pi    <- parameters[["pi"]]
  if(random){
    delta <- mu + parameters[["gamma"]] * parameters[["tau"]]
  }else{
    delta <- mu
  }

  ### compute the marginal log_likelihood
  log_lik <- 0

  # transform to probabilities
  p1 = .inv_logit( .logit(pi) + 0.5 * delta)
  p2 = .inv_logit( .logit(pi) - 0.5 * delta)

  if(!is.null(data[["weight"]])){
    log_lik <- log_lik + sum(stats::dbinom(x = data[["x1"]], prob = p1, size = data[["n1"]], log = TRUE) * data[["weight"]])
    log_lik <- log_lik + sum(stats::dbinom(x = data[["x2"]], prob = p2, size = data[["n2"]], log = TRUE) * data[["weight"]])
  }else{
    log_lik <- log_lik + sum(stats::dbinom(x = data[["x1"]], prob = p1, size = data[["n1"]], log = TRUE))
    log_lik <- log_lik + sum(stats::dbinom(x = data[["x2"]], prob = p2, size = data[["n2"]], log = TRUE))
  }

  return(log_lik)
}

# additional tools
.fitting_priority       <- function(models){

  # model fitting difficulty using the following heuristic:
  # selection models > random effects | PET/PEESE > non-null models
  fitting_difficulty <- sapply(models, function(model){

    difficulty <- 0

    if(inherits(model, "RoBMA.model") | inherits(model, "BiBMA.model")){
      if(is.prior.simple(model$priors[["mu"]])){
        difficulty <- difficulty + 1
      }
    }else if(.is_model_regression(model) | inherits(model, "BiBMA.reg.model")){
      difficulty <-  difficulty + sum(sapply(model$priors[["terms"]], function(prior){
        if(is.prior.point(prior)){
          return(0)
        }else if(is.prior.factor(prior)){
          return(1.5)
        }else if(is.prior.simple(prior)){
          return(1)
        }
      }))
    }
    if(is.prior.simple(model$priors[["tau"]])){
      difficulty <- difficulty + 3
    }
    if(is.null(model$priors[["PET"]])){
      difficulty <- difficulty + 1
    }else if(is.null(model$priors[["PEESE"]])){
      difficulty <- difficulty + 1
    }else if(is.null(model$priors[["omega"]])){
      difficulty <- difficulty + 5
    }
  })

  return(order(fitting_difficulty, decreasing = TRUE))
}
.runjags_summary_list   <- function(fit, priors, prior_scale, warnings, remove_parameters = NULL, measures = c("d", "r", "z", "logOR", "OR"), transformations_only = FALSE){

  summary_list <- list()

  for(measure in measures){

    # prepare transformations if necessary
    if(measure != prior_scale){
      transformations <- list()
      # effect size transformation
      if("mu" %in% names(priors) && ((is.prior.point(priors[["mu"]]) && priors[["mu"]][["parameters"]][["location"]] != 0) || !is.prior.point(priors[["mu"]]))){
        transformations[["mu"]] <- .get_transform_mu(prior_scale, measure, fun = FALSE)
      }
      # meta-regression transformations
      for(i in seq_along(priors[!names(priors) %in% c("mu", "tau", "omega", "PET", "PEESE", "bias")])){
        transformations[[names(priors[!names(priors) %in% c("mu", "tau", "omega", "PET", "PEESE")])[i]]] <- .get_transform_mu(prior_scale, measure, fun = FALSE)
      }
      # TODO: delete
      # for(i in seq_along(priors[!names(priors) %in% c("mu", "tau", "omega", "PET", "PEESE")])){
      #   transformations[[names(priors[!names(priors) %in% c("mu", "tau", "omega", "PET", "PEESE")])[i]]] <- list(
      #     "fun" = .transform_mu,
      #     "arg" = list(
      #       "from" = prior_scale,
      #       "to"   = measure
      #     )
      #   )
      # }

      # heterogeneity transformation
      if("tau" %in% names(priors) && ((is.prior.point(priors[["tau"]]) && priors[["tau"]][["parameters"]][["location"]] != 0) || !is.prior.point(priors[["tau"]]))){
        transformations[["tau"]] <- .get_scale(prior_scale, measure, fun = FALSE)
      }
      # PEESE transformation (PET is invariant)
      if("PEESE" %in% names(priors) || (!is.null(priors[["bias"]]) && any(sapply(priors[["bias"]], is.prior.PEESE)))){
        # the transformation is inverse for PEESE
        transformations[["PEESE"]] <- .get_scale(measure, prior_scale, fun = FALSE)
      }
      if(length(transformations) == 0){
        transformations <- NULL
      }
    }else{
      transformations <- NULL
    }

    if(transformations_only){
      if(length(measures) > 1)
        stop("Only one measure can be transformed at a time.")
      return(transformations)
    }

    summary_list[[measure]] <- suppressMessages(BayesTools::runjags_estimates_table(
      fit               = fit,
      transformations   = transformations,
      transform_factors = TRUE,
      formula_prefix    = FALSE,
      remove_inclusion  = TRUE,
      warnings          = warnings,
      footnotes         = .scale_note(prior_scale, measure),
      remove_parameters = remove_parameters
    ))
  }

  return(summary_list)
}
.get_id_weights         <- function(data, type){

  weight <- rep(NA, nrow(data))

  # create table of number of estimates per study
  if(!is.null(data[,"weight"]) && !anyNA(data[,"weight"])){
    weight <- data[,"weight"]
  }else{
    ids_weight <- data.frame(
      id     = names(table(data[,"study_ids"])),
      weight = switch(
        type,
        "inverse"      = 1/as.vector(table(data[,"study_ids"])),
        "inverse_sqrt" = 1/sqrt(as.vector(table(data[,"study_ids"])))
      )
    )

    # fill their weights
    for(i in seq_along(ids_weight$id)){
      weight[!is.na(data[,"study_ids"]) & data[,"study_ids"] == ids_weight$id[i]] <- ids_weight$weight[ids_weight$id == ids_weight$id[i]]
    }

    # assign all remaining studies weight 1
    weight[is.na(weight)] <- 1
  }


  return(weight)
}
.add_priors_levels      <- function(priors, data){

  for(v in names(data)[sapply(data, function(d) is.factor(d))]){
    attr(priors[[.BayesTools_parameter_name(v)]], "levels")      <- length(levels(data[[v]]))
    attr(priors[[.BayesTools_parameter_name(v)]], "level_names") <- levels(data[[v]])
    attr(priors[[.BayesTools_parameter_name(v)]], "parameter")   <- "mu"
  }

  return(priors)
}
.is_model_regression    <- function(model){
  return(inherits(model, "RoBMA.reg.model") || inherits(model, "RoBMA.reg.model_ss") || inherits(model, "BiBMA.reg.model") || inherits(model, "BiBMA.reg.model_ss"))
}
.is_model_multivariate  <- function(model){
  return(attr(model, "multivariate"))
}

.generate_model_formula_list       <- function(formula){

  # remove the left hand side
  if(attr(stats::terms(formula), "response") == 1){
    formula[2] <- NULL
  }
  formula <- list("mu" = formula)

  return(formula)
}
.generate_model_formula_data_list  <- function(data){

  if(length(data[["predictors"]]) == 0){
    data <- list("mu" = data.frame(matrix(ncol = 0, nrow = nrow(data[["outcome"]]))))
  }else{
    data <- list("mu" = data.frame(data[["predictors"]]))
  }

  return(data)
}
.generate_model_formula_prior_list <- function(priors){

  priors <- list("mu" = priors[["terms"]])

  return(priors)
}

# marginal likelihood computation function
.logit     <- function(x) log(x/(1-x))
.inv_logit <- function(x) 1/(1+exp(-x))

# JAGS tools for model building and marginal likelihood
.JAGS_transformation    <- function(from, to, from_par, to_par_name){

  return(paste0(to_par_name, " = ", .JAGS_transform(from, to, from_par), "\n"))

}
.JAGS_transform         <- function(from, to, par){

  if(from == to){
    return(par)
  }else{
    return(paste0(from, "2", to, "(", par,")"))
  }

}
.JAGS_transformation_se <- function(from, to, from_par, from_anc, to_par_name){

  return(paste0(to_par_name, " = ",.JAGS_transform_se(from, to, from_par, from_anc), "\n"))

}
.JAGS_transform_se      <- function(from, to, par, par_anc){

  if(from == to){
    return(par)
  }else{
    return(paste0("se_", from, "2", "se_", to, "(", par, if(!all(c(from, to) %in% c("logOR", "d"))) paste0(", ", par_anc),")"))
  }

}
.JAGS_scale             <- function(from, to, from_par, to_par_name){

  if(from == to){
    return(paste0(to_par_name, " = ", from_par, "\n"))
  }else{
    return(paste0(to_par_name, " = ", "scale_", from, "2", to, "(", from_par,")", "\n"))
  }

}
.JAGS_means_vector      <- function(study_id, mean){

  return(paste0(
    "for (i in 1:K_", study_id,"){\n",
    "  eff_", study_id, "[i] = ", mean, "\n",
    "}\n"
  ))

}
.JAGS_covariance_matrix <- function(study_id){

  return(paste0(
    "for (i in 1:K_", study_id,"){\n",
    "  sigma_", study_id, "[i,i] = pow(se_", study_id, "[i],2) + tau_transformed2\n",
    "  for(j in 1:(i-1)){\n",
    "    sigma_", study_id, "[i,j] = cov\n",
    "  }\n",
    "  for(j in (i+1):K_", study_id,"){\n",
    "    sigma_", study_id, "[i,j] = cov\n",
    "  }\n",
    "}\n"
  ))

}
