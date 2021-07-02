### functions for creating model objects
.check_and_list_priors            <- function(model_type, priors_effect_null, priors_effect, priors_heterogeneity_null, priors_heterogeneity, priors_bias_null, priors_bias, prior_scale){

  if(!is.null(model_type)){
    # precanned models
    if(model_type == "psma"){
      priors_effect         <- prior(distribution = "normal",    parameters = list(mean = 0,  sd = 1))
      priors_heterogeneity  <- prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15))
      priors_bias           <- list(
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
        prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
        prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),             prior_weights = 1/4),
        prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),             prior_weights = 1/4)
      )
      priors_effect_null        <- prior(distribution = "point", parameters = list(location = 0))
      priors_heterogeneity_null <- prior(distribution = "point", parameters = list(location = 0))
      priors_bias_null          <- prior_none()
    }else if(model_type == "pp"){
      priors_effect          <- prior(distribution = "normal",    parameters = list(mean = 0, sd = 1))
      priors_heterogeneity   <- prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15))
      priors_bias            <- list(
        prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),  prior_weights = 1/2),
        prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),  prior_weights = 1/2)
      )
      priors_effect_null         <- prior(distribution = "point", parameters = list(location = 0))
      priors_heterogeneity_null  <- prior(distribution = "point", parameters = list(location = 0))
      priors_bias_null           <- prior_none()
    }else if(model_type == "2w"){
      priors_effect              <- prior(distribution = "normal",    parameters = list(mean = 0, sd = 1))
      priors_heterogeneity       <- prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15))
      priors_bias                <- list(
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),        prior_weights = 1/2),
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),  prior_weights = 1/2)
      )
      priors_effect_null         <- prior(distribution = "point", parameters = list(location = 0))
      priors_heterogeneity_null  <- prior(distribution = "point", parameters = list(location = 0))
      priors_bias_null           <- prior_none()
    }else{
      stop("Unknown 'model_type'.")
    }
  }

  priors     <- list()
  priors$effect         <- .check_and_list_component_priors(priors_effect_null,          priors_effect,         "effect")
  priors$heterogeneity  <- .check_and_list_component_priors(priors_heterogeneity_null,   priors_heterogeneity,  "heterogeneity")
  priors$bias           <- .check_and_list_component_priors(priors_bias_null,            priors_bias,           "bias")

  return(priors)
}
.check_and_list_component_priors  <- function(priors_null, priors_alt, component){

  # check that at least one prior is specified (either null or alternative)
  if(is.null(priors_null) & is.null(priors_alt))
    stop(paste0("At least one prior needs to be specified for the ", component," parameter (either null or alternative)."))

  # create an empty list if user didn't specified priors
  if(is.null(priors_null)){
    priors_null <- list()
  }else{
    # make sure that priors are passed as a list
    if(is.prior(priors_null)){
      priors_null <- list(priors_null)
    }
    # mark the priors as null
    for(p in seq_along(priors_null)){
      priors_null[[p]]$is_null <- TRUE
    }
  }
  if(is.null(priors_alt)){
    priors_alt <- list()
  }else{
    # make sure that priors are passed as a list
    if(is.prior(priors_alt)){
      priors_alt <- list(priors_alt)
    }
    # mark the priors as alternative
    for(p in seq_along(priors_alt)){
      priors_alt[[p]]$is_null <- FALSE
    }
  }

  # join null and alternative priors
  priors <- c(priors_null, priors_alt)


  ### check that the specified priors are valid
  # check that all objects of the priors list are a 'RoBMA.prior'
  if(!all(sapply(priors, is.prior)))
    stop(paste0("Argument priors_", component, " does not contain valid prior distributions. The prior distributions need to be passed as a list of objects specified using 'prior()' function. See '?prior' for more information." ))


  # check that the passed priors are supported for the component (and replace none placeholders)
  if(component %in% c("effect", "heterogeneity")){

    for(p in seq_along(priors)){

      # replace empty priors by spike at 0
      if(is.prior.none(priors[[p]])){
        priors[[p]] <- prior("point", list("location" = 0))
      }

      # check for allowed priors
      if(!(is.prior.simple(priors[[p]]) | is.prior.point(priors[[p]])))
        stop(paste0("'", print(priors[[p]], silent = TRUE),"' prior distribution is not supported for the ", component," component."))

      # check for support in the case of heterogeneity
      if(component == "heterogeneity"){
        if(priors[[p]][["distribution"]] == "point" && priors[[p]]$parameters[["location"]] < 0){
          stop("The location of a point prior distribution for the heterogeneity component cannot be negative.")
        }else if(priors[[p]][["distribution"]] == "uniform" && (priors[[i]]$parameters[["a"]] < 0 | priors[[p]]$parameters[["b"]] < 0)){
          stop("The uniform prior distribution for the heterogeneity component cannot be defined on negative numbers.")
        }else if(priors[[p]]$truncation[["lower"]] < 0){
          priors[[p]]$truncation[["lower"]] <- 0
          warning(paste0("The range of a prior distribution for ", component, " component cannot be negative. The lower truncation point was set to zero."), immediate. = TRUE)
        }
      }

    }


  }else if(component == "bias"){

    for(p in seq_along(priors)){

      # check for allowed priors
      if(!(is.prior.PET(priors[[p]]) | is.prior.PEESE(priors[[p]]) | is.prior.weightfunction(priors[[p]]) | is.prior.none(priors[[p]])))
        stop(paste0("'", print(priors[[p]], silent = TRUE),"' prior distribution is not supported for the bias component."))
    }
  }

  return(priors)
}
.make_models  <- function(priors){

  # create models according to the set priors
  models <- NULL
  for(effect in priors[["effect"]]){
    for(heterogeneity in priors[["heterogeneity"]]){
      for(bias in priors[["bias"]]){
        models <- c(
          models,
          list(.make_model(effect, heterogeneity, bias, effect[["prior_weights"]] * heterogeneity[["prior_weights"]] * bias[["prior_weights"]]))
        )
      }
    }
  }

  return(models)
}
.make_model   <- function(prior_effect, prior_heterogeneity, prior_bias, prior_weights){

  priors <- list()

  priors$mu    <- prior_effect
  priors$tau   <- prior_heterogeneity
  if(is.prior.PET(prior_bias)){
    priors$PET    <- prior_bias
  }else if(is.prior.PEESE(prior_bias)){
    priors$PEESE  <- prior_bias
  }else if(is.prior.weightfunction(prior_bias)){
    priors$omega  <- prior_bias
  }

  model <- list(
    priors         = priors,
    prior_weights     = prior_weights,
    prior_weights_set = prior_weights
  )
  class(model) <- "RoBMA.model"

  return(model)
}
