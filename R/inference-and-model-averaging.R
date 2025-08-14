# Internal function to balance component probabilities when models fail to converge:
# Purpose: Redistributes prior weights from failed models to successful models with same structure
#
# Process:
# 1. Identifies converged vs failed models  
# 2. Characterizes each model by its component structure (effect, heterogeneity, bias, etc.)
# 3. For each failed model, finds converged models with identical component structure
# 4. Redistributes failed model's prior weight proportionally among matched converged models
# 5. Ensures balanced representation of component hypotheses despite convergence failures
#
# Key concepts:
# - Component structure: which components (effect, bias, etc.) are null vs alternative
# - Prior weight balancing: maintains intended relative hypothesis weights
# - Regression terms: handled separately if present (predictors_test)
# - Fallback: if no matching converged model exists, weight is lost (affects inference)
#
# This prevents convergence failures from biasing results toward particular hypotheses
# This function is available only for bridge sampling algorithm (no ballancing is possible for spike-and-slab)
.balance_component_probability   <- function(object){

  converged <- .get_model_convergence(object)
  if(all(!converged))
    stop("All models included in the ensemble failed to converge.")

  # Characterize each model by its component hypothesis structure
  # TRUE = alternative hypothesis, FALSE = null hypothesis for each component
  component_types <- cbind.data.frame(
    effect         = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "effect")),
    heterogeneity  = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "heterogeneity")),
    bias           = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "bias")),
    baseline       = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "baseline")),
    hierarchical   = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "hierarchical"))
  )
  
  # Add regression predictor tests if present
  if(!is.null(object$add_info[["predictors_test"]])){
    for(predictor_test in object$add_info[["predictors_test"]]){
      component_types[[predictor_test]] <- sapply(object[["models"]], function(model) predictor_test %in% model[["terms_test"]])
    }
  }
  
  # Extract original prior weights set by user
  prior_weights  <- sapply(object[["models"]], function(model) model[["prior_weights_set"]])

  # Redistribute weights from each failed model to converged models with same structure
  for(i in seq_along(object[["models"]])[!converged]){

    # Find converged models with identical component structure
    temp_same_structure <- sapply(seq_along(object$models)[converged], function(j) all(component_types[j,] == component_types[i,]))
    temp_same_structure <- seq_along(object$models)[converged][temp_same_structure]

    # Redistribute the failed model's prior weight proportionally
    if(length(temp_same_structure) >= 1){

      prior_weights[temp_same_structure] <- prior_weights[temp_same_structure] + prior_weights[i] / length(temp_same_structure)
      prior_weights[i] <- 0
      object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], "Some of the models failed to converge. However, there were other models with the same combination of model components and their prior probability was increased to account for the failed models.")

    }else{

      prior_weights[i] <- 0
      object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], "Some of the models failed to converge and their prior probability couldn't be balanced over models with the same combination of model components since they don't exist.")

    }
  }

  for(i in seq_along(object[["models"]])){
    object[["models"]][[i]][["prior_weights"]] <- prior_weights[i]
  }

  return(object)
}
.restore_component_probability   <- function(object){

  # extract the prior odds set by user
  prior_weights  <- sapply(object[["models"]], function(model) model[["prior_weights_set"]])

  for(i in seq_along(object[["models"]])){
    object[["models"]][[i]][["prior_weights"]] <- prior_weights[i]
  }

  return(object)

}

# Internal function to perform Bayesian model averaging across ensemble:
# Purpose: Combines posterior distributions from all converged models using marginal likelihoods
#
# Process:
# 1. Filters to converged models with non-zero prior weights
# 2. Characterizes models by component and parameter types  
# 3. Computes posterior model probabilities using marginal likelihoods
# 4. Performs model averaging for each parameter/component
# 5. Computes inclusion Bayes factors and posterior probabilities
# 6. Handles special cases (PET-PEESE combinations, regression terms)
#
# Key concepts:
# - Model averaging: weighted combination of posteriors using posterior model probabilities
# - Inclusion BF: evidence for component being included vs excluded
# - Component classification: effect, heterogeneity, bias, hierarchical, baseline
# - Parameter types: weight functions, PET, PEESE for publication bias
#
# Returns: List with averaged posteriors, inclusion BFs, and model probabilities
# .as_ensemble_inference is used for spike-and-slab algorithm (mimics output of .ensemble_inference for bridge sampling)
.ensemble_inference     <- function(object){

  # Filter to converged models with positive prior weights for inference
  prior_weights <- sapply(object[["models"]], function(model) model[["prior_weights"]])
  models        <- object[["models"]][.get_model_convergence(object) & prior_weights > 0]

  # Characterize models by component hypothesis types (null vs alternative)
  effect         <- sapply(models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(models, function(model)!.is_component_null(model[["priors"]], "bias"))
  hierarchical   <- sapply(models, function(model)!.is_component_null(model[["priors"]], "hierarchical"))
  baseline       <- sapply(models, function(model)!.is_component_null(model[["priors"]], "baseline"))

  # Characterize models by publication bias parameter types
  weightfunctions <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # define inference options: always effect and heterogeneity
  components      <- c("Effect", "Heterogeneity")
  parameters      <- c("mu", "tau")
  components_null <- list("Effect" = !effect, "Heterogeneity" = !heterogeneity)
  parameters_null <- list("mu"     = !effect, "tau"           = !heterogeneity)

  if(any(bias)){
    components      <- c(components,      "Bias")
    components_null <- c(components_null, "Bias" = list(!bias))
  }
  if(any(hierarchical)){
    components      <- c(components,      "Hierarchical")
    parameters      <- c(parameters,      "rho")
    components_null <- c(components_null, "Hierarchical" = list(!hierarchical))
    parameters_null <- c(parameters_null, "rho" = list(!hierarchical))
  }
  if(any(baseline)){
    components      <- c(components,      "Baseline")
    components_null <- c(components_null, "Baseline" = list(!baseline))
  }
  if(any(weightfunctions)){
    components      <- c(components,      "bias.selection-models")
    parameters      <- c(parameters,      "omega")
    components_null <- c(components_null, "bias.selection-models" = list(!(weightfunctions & bias)))
    parameters_null <- c(parameters_null, "omega" = list(!weightfunctions))
  }
  if(any(PET) & any(PEESE)){
    components      <- c(components,      "bias.PET-PEESE")
    parameters      <- c(parameters,      "PET", "PEESE")
    components_null <- c(components_null, "bias.PET-PEESE" = list(!((PET | PEESE) & bias)))
    parameters_null <- c(parameters_null, "PET" = list(!PET), "PEESE" = list(!PEESE))
  }else if(any(PET)){
    components      <- c(components,      "bias.PET")
    parameters      <- c(parameters,      "PET")
    components_null <- c(components_null, "bias.PET" = list(!(PET & bias)))
    parameters_null <- c(parameters_null, "PET" = list(!PET))
  }else if(any(PEESE)){
    components      <- c(components,      "bias.PEESE")
    parameters      <- c(parameters,      "PEESE")
    components_null <- c(components_null, "bias.PEESE" = list(!(PEESE & bias)))
    parameters_null <- c(parameters_null, "PEESE" = list(!PEESE))
  }

  # deal with meta-regression
  if(!is.null(object[["formula"]])){
    # use the intercept for the effect
    parameters[parameters == "mu"]                         <- "mu_intercept"
    names(parameters_null)[names(parameters_null) == "mu"] <- "mu_intercept"

    # add the terms
    model_predictors      <- lapply(models, function(model) model[["terms"]])
    model_predictors_test <- lapply(models, function(model) model[["terms_test"]])

    predictors      <- object$add_info[["predictors"]]
    predictors_test <- object$add_info[["predictors_test"]]

    # define inference options
    components_predictors      <- NULL
    parameters_predictors      <- "mu_intercept"
    parameters_marginal        <- "mu_intercept"
    components_predictors_null <- list()
    parameters_predictors_null <- list("mu_intercept" = !effect)
    parameters_marginal_null   <- list("mu_intercept" = !effect)

    components_predictors_distributions      <- NULL
    components_predictors_distributions_null <- list()


    # predictors
    for(i in seq_along(predictors_test)){
      components_predictors <- c(components_predictors, .BayesTools_parameter_name(predictors_test[i]))
      components_predictors_null[[.BayesTools_parameter_name(predictors_test[i])]] <-
        sapply(model_predictors_test, function(x) if(length(x) == 0) TRUE else !(predictors_test[i] %in% x))
    }

    for(i in seq_along(predictors)){
      parameters_predictors <- c(parameters_predictors, .BayesTools_parameter_name(predictors[i]))
      parameters_marginal   <- c(parameters_marginal,   .BayesTools_parameter_name(predictors[i]))
      parameters_predictors_null[[.BayesTools_parameter_name(predictors[i])]] <-
        sapply(model_predictors_test, function(x) !(predictors[i] %in% x))
      parameters_marginal_null[[.BayesTools_parameter_name(predictors[i])]] <-
        sapply(seq_along(models), function(j) !((predictors[i] %in% model_predictors_test[[j]]) || effect[j]))
    }


    # get models inference
    if(is.null(components_predictors)){
      inference_predictors <- NULL
    }else{
      inference_predictors <- BayesTools::ensemble_inference(
        model_list   = models,
        parameters   = components_predictors,
        is_null_list = components_predictors_null,
        conditional  = FALSE
      )
    }
    # deal with the possibility of only null models models
    if(all(sapply(parameters_predictors_null, all))){
      inference_predictors_conditional <- NULL
    }else{
      inference_predictors_conditional <- BayesTools::ensemble_inference(
        model_list   = models,
        parameters   = parameters_predictors[!sapply(parameters_predictors_null, all)],
        is_null_list = parameters_predictors_null[!sapply(parameters_predictors_null, all)],
        conditional  = TRUE
      )
    }

    # get model-averaged posteriors
    if(is.null(parameters_predictors)){
      posteriors_predictors <- NULL
    }else{
      posteriors_predictors <- BayesTools::mix_posteriors(
        model_list   = models,
        parameters   = parameters_predictors,
        is_null_list = parameters_predictors_null,
        seed         = object$add_info[["seed"]],
        conditional  = FALSE
      )
      posteriors_predictors <- BayesTools::transform_factor_samples(posteriors_predictors)
    }

    # deal with the possibility of only null models models
    if(all(sapply(parameters_predictors_null, all))){
      posteriors_predictors_conditional <- NULL
    }else{
      posteriors_predictors_conditional <- BayesTools::mix_posteriors(
        model_list   = models,
        parameters   = parameters_predictors[!sapply(parameters_predictors_null, all)],
        is_null_list = parameters_predictors_null[!sapply(parameters_predictors_null, all)],
        seed         = object$add_info[["seed"]],
        conditional  = TRUE
      )
      posteriors_predictors_conditional <- BayesTools::transform_factor_samples(posteriors_predictors_conditional)
    }

    # create marginal estimates and summary
    if(all(sapply(parameters_marginal_null, all))){
      inference_marginal <- NULL
    }else{
      inference_marginal <- BayesTools::marginal_inference(
        model_list          = models,
        marginal_parameters = parameters_marginal[!sapply(parameters_marginal_null, all)],
        parameters          = parameters_marginal,
        is_null_list        = parameters_marginal_null,
        formula             = object[["formula"]],
        seed                = object$add_info[["seed"]],
        silent              = TRUE
      )
    }
  }else{
    # create empty objects in case of no predictors
    inference_predictors              <- NULL
    inference_predictors_conditional  <- NULL
    posteriors_predictors             <- NULL
    posteriors_predictors_conditional <- NULL
    inference_marginal                <- NULL
  }

  ### get models inference
  inference <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components,
    is_null_list = components_null,
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    inference_conditional <- NULL
  }else{
    inference_conditional <- BayesTools::ensemble_inference(
      model_list   = models,
      parameters   = components[!sapply(components_null, all)],
      is_null_list = components_null[!sapply(components_null, all)],
      conditional  = TRUE
    )
  }


  ### get model-averaged posteriors
  posteriors <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters,
    is_null_list = parameters_null,
    seed         = object$add_info[["seed"]],
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    posteriors_conditional <- NULL
  }else{
    posteriors_conditional <- BayesTools::mix_posteriors(
      model_list   = models,
      parameters   = parameters[!sapply(parameters_null, all)],
      is_null_list = parameters_null[!sapply(parameters_null, all)],
      seed         = object$add_info[["seed"]],
      conditional  = TRUE
    )
  }

  # rename mu_intercept back to mu
  if(any(names(posteriors) == "mu_intercept")){
    attr(posteriors[["mu_intercept"]], "parameter")        <- "mu"
    names(posteriors)[names(posteriors) == "mu_intercept"] <- "mu"
  }
  if(any(names(posteriors_conditional) == "mu_intercept")){
    attr(posteriors_conditional[["mu_intercept"]], "parameter")                    <- "mu"
    names(posteriors_conditional)[names(posteriors_conditional) == "mu_intercept"] <- "mu"
  }

  # return the results
  output <- list(
    inference                         = inference,
    inference_conditional             = inference_conditional,
    inference_predictors              = inference_predictors,
    inference_predictors_conditional  = inference_predictors_conditional,
    inference_marginal                = inference_marginal,
    posteriors                        = posteriors,
    posteriors_conditional            = posteriors_conditional,
    posteriors_predictors             = posteriors_predictors,
    posteriors_predictors_conditional = posteriors_predictors_conditional
  )
  return(output)
}
.as_ensemble_inference  <- function(object){

  # use only converged models with prior weights > 0 for inference about parameters
  model  <- object[["model"]]
  priors <- object[["model"]][["priors"]]

  # obtain the component type
  effect         <- if(.is_model_regression(model)) !.is.prior_null(priors[["terms"]][["intercept"]]) else !.is.prior_null(priors[["mu"]])
  heterogeneity  <- if(!is.null(priors[["tau"]]))          !.is.prior_null(priors[["tau"]])
  bias           <- if(!is.null(priors[["bias"]]))         !.is.prior_null(priors[["bias"]])
  hierarchical   <- if(!is.null(priors[["rho"]]))          !.is.prior_null(priors[["rho"]])
  baseline       <- if(!is.null(priors[["baseline"]]))     !.is.prior_null(priors[["baseline"]])

  # obtain the parameter types
  if(length(priors[["bias"]]) > 0){
    weightfunctions <- sapply(priors[["bias"]], is.prior.weightfunction)
    PET             <- sapply(priors[["bias"]], is.prior.PET)
    PEESE           <- sapply(priors[["bias"]], is.prior.PEESE)
  }

  ### prepare ensemble inference objects
  # further separation between predictors and main components etc is done directly in the summary function
  inference <- BayesTools::runjags_inference_table(model[["fit"]], formula_prefix = FALSE)

  # rename the main component types
  rownames(inference)[rownames(inference) == "tau"]  <- "Heterogeneity"
  rownames(inference)[rownames(inference) == "bias"] <- "Bias"
  rownames(inference)[rownames(inference) == "rho"]  <- "Hierarchical"
  rownames(inference)[rownames(inference) == "pi"]   <- "Baseline"

  if(.is_model_regression(model)){

    # intercept = average effect for meta-regression
    rownames(inference)[rownames(inference) == "intercept"] <- "Effect"

    # add marginal inference
    parameters          <- NULL
    marginal_parameters <- NULL
    conditional_list    <- list()

    for(i in seq_along(priors[["terms"]])){
      parameters <- c(parameters, .BayesTools_parameter_name(names(priors[["terms"]])[i]))
      marginal_parameters <- c(marginal_parameters, .BayesTools_parameter_name(names(priors[["terms"]])[i]))
      conditional_list[[.BayesTools_parameter_name(names(priors[["terms"]])[i])]] <- c(
        if(names(priors[["terms"]])[i] != "intercept") .BayesTools_parameter_name("intercept"),
        .BayesTools_parameter_name(names(priors[["terms"]])[i])
      )
    }

    inference_marginal <- suppressWarnings(BayesTools::as_marginal_inference(
      model               = model[["fit"]],
      marginal_parameters = marginal_parameters,
      parameters          = parameters,
      conditional_list    = conditional_list,
      conditional_rule    = "OR",
      formula             = object[["formula"]],
      silent              = TRUE,
      force_plots         = TRUE
    ))
  }else{

    # rename effect
    rownames(inference)[rownames(inference) == "mu"]        <- "Effect"

    # add empty marginal inference object
    inference_marginal <- NULL
  }

  # split the component and predictors inference
  inference_predictors <- inference[!rownames(inference) %in% c("Effect", "Heterogeneity", "Bias", "Hierarchical", "Baseline"), ]
  inference            <- inference[ rownames(inference) %in% c("Effect", "Heterogeneity", "Bias", "Hierarchical", "Baseline"), ]
  if(nrow(inference_predictors) == 0){
    inference_predictors <- NULL
  }


  ### prepare model-averaged posteriors
  # define inference options: always effect and heterogeneity
  parameters   <- c("mu", "tau")
  conditional  <- list()

  if(any(effect)){
    conditional <- c(conditional, "mu" = "mu")
  }
  if(any(heterogeneity)){
    conditional <- c(conditional, "tau" = "tau")
  }
  if(!is.null(priors[["bias"]]) && !is.prior.none(priors[["bias"]])){
    parameters  <- c(parameters, "bias")
  }
  if(any(bias) && any(weightfunctions)){
    conditional <- c(conditional, "bias" = "omega")
  }
  if(any(bias) && any(PET)){
    conditional <- c(conditional, "bias" = "PET")
  }
  if(any(bias) && any(PEESE)){
    conditional <- c(conditional, "bias" = "PEESE")
  }
  if(any(hierarchical)){
    parameters  <- c(parameters,  "rho")
    conditional <- c(conditional, "rho" = "rho")
  }
  if(any(baseline)){
    parameters  <- c(parameters,  "baseline")
    conditional <- c(conditional, "baseline" = "baseline")
  }

  # deal with meta-regression
  if(.is_model_regression(model)){
    # use the intercept for the effect
    parameters[parameters == "mu"]                 <- "mu_intercept"
    names(conditional)[names(conditional) == "mu"] <- "mu_intercept"
    conditional[["mu_intercept"]]                  <- "mu_intercept"
  }

  # get model-averaged posteriors
  posteriors <- BayesTools::as_mixed_posteriors(
    model        = model[["fit"]],
    parameters   = parameters
  )

  # simplify the model-averaged bias distribution
  if("bias" %in% parameters){
    if(any(weightfunctions)){
      posteriors[["omega"]] <- posteriors[["bias"]][,grepl("omega", colnames(posteriors[["bias"]])),drop=FALSE]
      attributes(posteriors[["omega"]]) <- c(
        attributes(posteriors[["omega"]])[names(attributes(posteriors[["omega"]])) %in% c("dim", "dimnames")],
        attributes(posteriors[["bias"]])[!names(attributes(posteriors[["bias"]])) %in% c("dim", "dimnames")]
      )
    }
    if(any(PET)){
      posteriors[["PET"]] <- posteriors[["bias"]][,grepl("PET", colnames(posteriors[["bias"]])),drop=FALSE]
      attributes(posteriors[["PET"]]) <- c(
        attributes(posteriors[["PET"]])[names(attributes(posteriors[["PET"]])) %in% c("dim", "dimnames")],
        attributes(posteriors[["bias"]])[!names(attributes(posteriors[["bias"]])) %in% c("dim", "dimnames")]
      )
    }
    if(any(PEESE)){
      posteriors[["PEESE"]] <- posteriors[["bias"]][,grepl("PEESE", colnames(posteriors[["bias"]])),drop=FALSE]
      attributes(posteriors[["PEESE"]]) <- c(
        attributes(posteriors[["PEESE"]])[names(attributes(posteriors[["PEESE"]])) %in% c("dim", "dimnames")],
        attributes(posteriors[["bias"]])[!names(attributes(posteriors[["bias"]])) %in% c("dim", "dimnames")]
      )
    }
    posteriors[["bias"]] <- NULL
  }

  # conditional posteriors
  if(length(conditional) > 0){
    posteriors_conditional <- list()
    for(i in seq_along(conditional)){
      posteriors_conditional[[conditional[[i]]]] <- BayesTools::as_mixed_posteriors(
        model        = model[["fit"]],
        parameters   = names(conditional)[i],
        conditional  = conditional[[i]],
        force_plots  = TRUE
      )[[names(conditional)[i]]]
    }
  }else{
    posteriors_conditional <- NULL
  }

  # rename mu_intercept back to mu
  if(any(names(posteriors) == "mu_intercept")){
    attr(posteriors[["mu_intercept"]], "parameter")        <- "mu"
    names(posteriors)[names(posteriors) == "mu_intercept"] <- "mu"
  }
  if(any(names(posteriors_conditional) == "mu_intercept")){
    attr(posteriors_conditional[["mu_intercept"]], "parameter")                    <- "mu"
    names(posteriors_conditional)[names(posteriors_conditional) == "mu_intercept"] <- "mu"
  }

  ### model-averaged posteriors for predictors
  if(.is_model_regression(model)){

    parameters_predictors  <- NULL
    conditional_predictors <- list()

    # add the terms
    for(i in seq_along(priors[["terms"]])){
      parameters_predictors <- c(parameters_predictors, .BayesTools_parameter_name(names(priors[["terms"]])[i]))
      conditional_predictors[[.BayesTools_parameter_name(names(priors[["terms"]])[i])]] <- .BayesTools_parameter_name(names(priors[["terms"]])[i])
    }

    posteriors_predictors <- BayesTools::as_mixed_posteriors(
      model        = model[["fit"]],
      parameters   = parameters_predictors
    )
    posteriors_predictors <- BayesTools::transform_factor_samples(posteriors_predictors)

    # conditional posteriors
    if(length(conditional_predictors) > 0){
      posteriors_predictors_conditional <- list()
      for(i in seq_along(conditional_predictors)){
        # suppress the not a conditional parameter warning in order to get "conditional estimates"
        posteriors_predictors_conditional[[conditional_predictors[[i]]]] <- suppressWarnings(BayesTools::as_mixed_posteriors(
          model        = model[["fit"]],
          parameters   = names(conditional_predictors)[i],
          conditional  = conditional_predictors[[i]]
        ))[[names(conditional_predictors)[i]]]
      }
      posteriors_predictors_conditional <- BayesTools::transform_factor_samples(posteriors_predictors_conditional)
    }else{
      posteriors_predictors_conditional <- NULL
    }

  }else{
    posteriors_predictors             <- NULL
    posteriors_predictors_conditional <- NULL
  }


  # return the results
  output <- list(
    inference                         = inference,
    inference_predictors              = inference_predictors,
    inference_marginal                = inference_marginal,
    posteriors                        = posteriors,
    posteriors_conditional            = posteriors_conditional,
    posteriors_predictors             = posteriors_predictors,
    posteriors_predictors_conditional = posteriors_predictors_conditional
  )
  return(output)
}
.compute_coeficients    <- function(RoBMA){
  if(!is.null(RoBMA[["posteriors_predictors"]])){
    coefs <- do.call(c, unname(lapply(RoBMA[["posteriors_predictors"]], function(posterior){
      if(inherits(posterior, "mixed_posteriors.factor")){
        out        <- apply(posterior, 2, mean)
        names(out) <- .output_parameter_names(names(out))
      }else{
        out        <- mean(posterior)
        names(out) <- .output_parameter_names(attr(posterior,"parameter"))
      }
      return(out)
    })))
  }else{
    coefs        <- c("mu" = mean(RoBMA$posteriors[["mu"]]))
  }
  return(c(
    coefs,
    "tau"    = mean(RoBMA$posteriors[["tau"]]),
    "rho"    = if(!is.null(RoBMA$posteriors[["rho"]]))   mean(RoBMA$posteriors[["rho"]]),
    if(!is.null(RoBMA$posteriors[["omega"]])) apply(RoBMA$posteriors[["omega"]], 2, mean),
    "PET"    = if(!is.null(RoBMA$posteriors[["PET"]]))   mean(RoBMA$posteriors[["PET"]]),
    "PEESE"  = if(!is.null(RoBMA$posteriors[["PEESE"]])) mean(RoBMA$posteriors[["PEESE"]])
  ))
}
