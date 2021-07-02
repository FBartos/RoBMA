.balance_probability   <- function(object){

  converged  <- object$add_info[["converged"]]
  # assess the component type
  effect         <- sapply(object[["models"]], function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(object[["models"]], function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(object[["models"]], function(model)!.is_component_null(model[["priors"]], "bias"))
  # extract the prior odds set by user
  prior_weights     <- sapply(object[["models"]], function(model)model[["prior_weights_set"]])


  # check whether there is a comparable model for each non-converged models
  for(i in seq_along(object[["models"]])[!converged]){

    temp_ind  <- seq_along(object[["models"]])[-i]
    temp_same <- temp_ind[effect[-i] == effect[i] & heterogeneity[-i] == heterogeneity[i] & bias[-i] == bias[i] & converged[-i]]

    # if yes, transfer the prior odds
    if(length(temp_same) >= 1){

      prior_weights[temp_same] <- prior_weights[temp_same] + prior_weights[i] / length(temp_same)
      prior_weights[i] <- 0
      object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], "Some of the models failed to converge. However, there were other models with the same combination of presence/absence of effect/heterogeneity/publication bias and their prior probability was increased to account for the failed models.")

    }else{

      prior_weights[i] <- 0
      object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], "Some of the models failed to converge and their prior probability couldn't be balanced over models with the same combination of presence/absence of effect/heterogeneity/publication bias since they don't exist.")

    }
  }

  for(i in seq_along(object[["models"]])[!converged]){
    object[["models"]][[i]][["prior_weights"]] <- prior_weights[i]
  }

  return(object)
}
.ensemble_inference    <- function(object){

  # use only converged models with prior weights > 0 for inference about parameters
  prior_weights <- sapply(object[["models"]], function(model) model[["prior_weights"]])
  models        <- object[["models"]][.get_model_convergence(object) & prior_weights > 0]

  # obtain the component type
  effect         <- sapply(models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(models, function(model)!.is_component_null(model[["priors"]], "bias"))

  # obtain the parameter types
  weightfunctions <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # define inference options
  components      <- c("Effect", "Heterogeneity", "Bias")
  parameters      <- c("mu", "tau")
  components_null <- list("Effect" = !effect, "Heterogeneity" = !heterogeneity, "Bias" = !bias)
  parameters_null <- list("mu"     = !effect, "tau"           = !heterogeneity)

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

  # get models inference
  inference <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components,
    is_null_list = components_null,
    conditional  = FALSE
  )
  inference_conditional <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components,
    is_null_list = components_null,
    conditional  = TRUE
  )

  # get model-averaged posteriors
  posteriors <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters,
    is_null_list = parameters_null,
    seed         = object$add_info[["seed"]],
    conditional  = FALSE
  )
  posteriors_conditional <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters,
    is_null_list = parameters_null,
    seed         = object$add_info[["seed"]],
    conditional  = TRUE
  )


  # return the results
  output <- list(
    inference              = inference,
    inference_conditional  = inference_conditional,
    posteriors             = posteriors,
    posteriors_conditional = posteriors_conditional
  )
  return(output)
}

.compute_coeficients   <- function(RoBMA){
  return(c(
    "mu"     = mean(RoBMA$posteriors[["mu"]]),
    "tau"    = mean(RoBMA$posteriors[["tau"]]),
    if(!is.null(RoBMA$posteriors[["omega"]])) apply(RoBMA$posteriors[["omega"]], 2, mean),
    "PET"    = if(!is.null(RoBMA$posteriors[["PET"]]))   mean(RoBMA$posteriors[["PET"]]),
    "PEESE"  = if(!is.null(RoBMA$posteriors[["PEESE"]])) mean(RoBMA$posteriors[["PEESE"]])
  ))
}
