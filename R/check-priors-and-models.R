### functions for creating model objects
.check_and_list_priors      <- function(model_type, priors_effect_null, priors_effect, priors_heterogeneity_null, priors_heterogeneity, priors_bias_null, priors_bias, priors_hierarchical_null, priors_hierarchical, prior_scale){

  # format the model-type (in RoBMA.reg, the add_info check is run only after .check_and_list_priors)
  model_type <- .check_and_set_model_type(model_type, prior_scale)

  if(!is.null(model_type) & length(model_type == 1)){
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
        prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),   prior_weights = 1/4),
        prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf), prior_weights = 1/4)
      )
      priors_hierarchical       <- NULL
      priors_effect_null        <- prior(distribution = "point", parameters = list(location = 0))
      priors_heterogeneity_null <- prior(distribution = "point", parameters = list(location = 0))
      priors_bias_null          <- prior_none()
      priors_hierarchical_null  <- NULL
    }else if(model_type == "pp"){
      priors_effect          <- prior(distribution = "normal",    parameters = list(mean = 0, sd = 1))
      priors_heterogeneity   <- prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15))
      priors_bias            <- list(
        prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),    prior_weights = 1/2),
        prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),  prior_weights = 1/2)
      )
      priors_hierarchical        <- NULL
      priors_effect_null         <- prior(distribution = "point", parameters = list(location = 0))
      priors_heterogeneity_null  <- prior(distribution = "point", parameters = list(location = 0))
      priors_bias_null           <- prior_none()
      priors_hierarchical_null   <- NULL
    }else if(model_type == "2w"){
      priors_effect              <- prior(distribution = "normal",    parameters = list(mean = 0, sd = 1))
      priors_heterogeneity       <- prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15))
      priors_bias                <- list(
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),        prior_weights = 1/2),
        prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),  prior_weights = 1/2)
      )
      priors_hierarchical        <- NULL
      priors_effect_null         <- prior(distribution = "point", parameters = list(location = 0))
      priors_heterogeneity_null  <- prior(distribution = "point", parameters = list(location = 0))
      priors_bias_null           <- prior_none()
      priors_hierarchical_null   <- NULL
    }else{
      stop("Unknown 'model_type'.")
    }
  }

  priors     <- list()
  priors$effect         <- .check_and_list_component_priors(priors_effect_null,          priors_effect,         "effect")
  priors$heterogeneity  <- .check_and_list_component_priors(priors_heterogeneity_null,   priors_heterogeneity,  "heterogeneity")
  priors$bias           <- .check_and_list_component_priors(priors_bias_null,            priors_bias,           "bias")
  priors$hierarchical   <- .check_and_list_component_priors(priors_hierarchical_null,    priors_hierarchical,   "hierarchical")

  return(priors)
}
.check_and_list_priors.reg  <- function(priors, data, model_type, test_predictors, prior_scale,
                                        priors_effect_null, priors_effect,
                                        priors_heterogeneity_null, priors_heterogeneity,
                                        priors_bias_null, priors_bias,
                                        priors_hierarchical_null, priors_hierarchical,
                                        prior_covariates_null, prior_covariates,
                                        prior_factors_null, prior_factors){

  priors_main <- .check_and_list_priors(
    model_type                = model_type,
    priors_effect_null        = priors_effect_null,
    priors_effect             = priors_effect,
    priors_heterogeneity_null = priors_heterogeneity_null,
    priors_heterogeneity      = priors_heterogeneity,
    priors_bias_null          = priors_bias_null,
    priors_bias               = priors_bias,
    priors_hierarchical_null  = priors_hierarchical_null,
    priors_hierarchical       = priors_hierarchical,
    prior_scale               = prior_scale
  )

  predictors      <- attr(data[["predictors"]],"terms")
  predictors_type <- attr(data[["predictors"]],"terms_type")

  # check the input
  if(!is.prior.simple(prior_covariates_null) || is.prior.factor(prior_covariates_null))
    stop("The default prior for covariates (null) is not a valid prior distribution.", call. = FALSE)
  if(!is.prior.simple(prior_covariates) || is.prior.factor(prior_covariates))
    stop("The default prior for covariates is not a valid prior distribution.", call. = FALSE)
  if(!is.prior.factor(prior_factors_null) & !is.prior.point(prior_factors_null))
    stop("The default prior for factors (null) is not a valid prior distribution.", call. = FALSE)
  if(!is.prior.factor(prior_factors) & !is.prior.point(prior_factors))
    stop("The default prior for factors is not a valid prior distribution.", call. = FALSE)

  # check for reserved words
  if(any(names(priors) %in% .reserved_words()))
    stop(paste0("The following prior names are internally reserved keywords and cannot be used: ",
                paste0(" '", names(priors)[names(priors) %in% .reserved_words()], "' ", collapse = ", ")), call. = FALSE)

  # completely the prior distribution specification
  if(is.null(priors) && (!is.null(test_predictors) && length(test_predictors) == 1 && isFALSE(test_predictors))){
    # default estimation if no priors and test_predictors is false
    test_predictors <- NULL

  }else if(is.null(priors) && is.null(test_predictors)){
    # complete default - tests all predictors with default priors
    test_predictors <- predictors

  }else if(!is.null(priors)){
    # find whether user specified some parameter priors, if not - tests all predictors with default priors
    predictors_prior_info <- unlist(sapply(predictors, function(p){
      if(is.null(priors[[p]])){
        return("no-priors-are-specified")
      }else if(is.prior(priors[[p]])){
        return("one-prior-is-specified")
      }else if(length(priors[[p]]) == 2 && all(names(priors[[p]]) %in% c("null", "alt"))){
        if(all(sapply(priors[[p]], is.prior))){
          return(p)
        }else{
          stop(paste0("The prior distribution for '",p,"' is specified incorrectly."))
        }
      }else{
        stop(paste0("The prior distribution for '",p,"' is specified incorrectly."))
      }
    }))
    if(isFALSE(test_predictors)){
      test_predictors <- NULL
    }else if(is.null(test_predictors)){
      test_predictors <- predictors[!predictors %in% c("one-prior-is-specified", "no-priors-are-specified")]
    }
  }


  if(is.null(priors)){

    priors <- list()

    to_test <- predictors[predictors %in% test_predictors]
    no_test <- predictors[!predictors %in% test_predictors]

    for(i in seq_along(to_test)){
      priors[[to_test[i]]] <- list(
        null = if(predictors_type[to_test[i]] == "factor") prior_factors_null else prior_covariates_null,
        alt  = if(predictors_type[to_test[i]] == "factor") prior_factors      else prior_covariates
      )
    }
    for(i in seq_along(no_test)){
      priors[[no_test[i]]] <- list(
        alt  = if(predictors_type[no_test[i]] == "factor") prior_factors else prior_covariates
      )
    }

  }else{

    if(any(!names(priors) %in% predictors))
      stop(paste0("The following priors do not corresponds to any predictor or additional parameter: '", paste(names(priors)[!names(priors) %in% predictors], collapse = "', '", sep = ""), "'"))


    to_test <- predictors[predictors %in% test_predictors]
    no_test <- predictors[!predictors %in% test_predictors]

    for(i in seq_along(to_test)){
      if(is.null(priors[[to_test[i]]])){
        priors[[to_test[i]]] <- list(
          null = if(predictors_type[to_test[i]] == "factor") prior_factors_null else prior_covariates_null,
          alt  = if(predictors_type[to_test[i]] == "factor") prior_factors      else prior_covariates
        )
      }else if(is.prior(priors[[to_test[i]]])){
        priors[[to_test[i]]] <- list(
          null  = if(predictors_type[to_test[i]] == "factor") prior_factors_null else prior_covariates_null,
          alt   = priors[[to_test[i]]]
        )
      }else if(length(priors[[to_test[i]]]) == 2 && all(names(priors[[to_test[i]]]) %in% c("null", "alt"))){
        priors[[to_test[i]]] <- list(
          null  = priors[[to_test[i]]][["null"]],
          alt   = priors[[to_test[i]]][["alt"]]
        )
      }else{
        stop(paste0("The predictor '", to_test[i], "' is supposed to be used for testing and the prior distributions are not specified properly"))
      }
    }

    for(i in seq_along(no_test)){
      if(is.null(priors[[no_test[i]]])){
        priors[[no_test[i]]] <- list(
          alt  = if(predictors_type[no_test[i]] == "factor") prior_factors  else prior_covariates
        )
      }else{
        if(is.prior(priors[[no_test[i]]])){
          priors[[no_test[i]]] <- list(
            alt  = priors[[no_test[i]]]
          )
        }else{
          # should be stopped before
          stop(paste0("The predictor '", no_test[i], "' is supposed to be used for testing and the prior distributions are not specified properly"))
        }
      }
    }
  }

  priors_main$terms <- priors[predictors]

  attr(priors_main, "terms")          <- predictors
  attr(priors_main, "terms_test")     <- if(length(test_predictors) == 1 && test_predictors == "") NULL else test_predictors

  return(priors_main)
}
.check_and_list_component_priors  <- function(priors_null, priors_alt, component){

  # check that at least one prior is specified (either null or alternative)
  if(component != "hierarchical" && (is.null(priors_null) & is.null(priors_alt)))
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
  if(component %in% c("effect", "heterogeneity", "hierarchical")){

    for(p in seq_along(priors)){

      # replace empty priors by spike at 0
      if(is.prior.none(priors[[p]])){
        temp_is_null         <- priors[[p]]$is_null
        priors[[p]]          <- prior("point", list("location" = 0))
        priors[[p]]$is_null  <- temp_is_null
      }

      # check for allowed priors
      if(!(is.prior.simple(priors[[p]]) | is.prior.point(priors[[p]])))
        stop(paste0("'", print(priors[[p]], silent = TRUE),"' prior distribution is not supported for the ", component," component."))

      # check for support in the case of heterogeneity
      if(component == "heterogeneity"){
        if(priors[[p]][["distribution"]] == "point" && priors[[p]]$parameters[["location"]] < 0){
          stop("The location of a point prior distribution for the heterogeneity component cannot be negative.")
        }else if(priors[[p]][["distribution"]] == "uniform" && (priors[[p]]$parameters[["a"]] < 0 | priors[[p]]$parameters[["b"]] < 0)){
          stop("The uniform prior distribution for the heterogeneity component cannot be defined on negative numbers.")
        }else if(priors[[p]]$truncation[["lower"]] < 0){
          priors[[p]]$truncation[["lower"]] <- 0
          warning(paste0("The range of a prior distribution for ", component, " component cannot be negative. The lower truncation point was set to zero."), immediate. = TRUE)
        }
      }else if(component == "hierarchical"){
        if(priors[[p]][["distribution"]] == "point" && abs(priors[[p]]$parameters[["location"]]) > 1){
          stop("The location of a point prior distribution for the hierarchical correlation must be within [-1, 1] interval.")
        }else if(priors[[p]][["distribution"]] == "uniform" && (priors[[p]]$parameters[["a"]] < -1 | priors[[p]]$parameters[["b"]] > 1)){
          stop("The uniform prior distribution for the hierarchical correlation cannot be defined outsied of the [-1, 1] interval.")
        }

        if(priors[[p]]$truncation[["lower"]] < -1){
          priors[[p]]$truncation[["lower"]] <- -1
          warning("The range of a prior distribution for the hierarchical correlation cannot be lower than -1. The lower truncation point was set to -1.", immediate. = TRUE)
        }
        if(priors[[p]]$truncation[["upper"]] > 1){
          priors[[p]]$truncation[["lower"]] <- 1
          warning("The range of a prior distribution for the hierarchical correlation cannot be higher than 1. The upper truncation point was set to 1.", immediate. = TRUE)
        }
        if(priors[[p]]$truncation[["lower"]] > priors[[p]]$truncation[["upper"]]){
          stop("Invalid lower and upper truncation points for the hierarchical correlation.", immediate. = TRUE)
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

.make_models     <- function(priors, multivariate, weighted){

  # create models according to the set priors
  models <- NULL
  for(effect in priors[["effect"]]){
    for(heterogeneity in priors[["heterogeneity"]]){
      for(bias in priors[["bias"]]){
        if(!is.null(priors[["hierarchical"]]) && multivariate){
          for(hierarchical in priors[["hierarchical"]]){
            models <- c(
              models,
              list(.make_model(effect, heterogeneity, bias, hierarchical, effect[["prior_weights"]] * heterogeneity[["prior_weights"]] * bias[["prior_weights"]] * hierarchical[["prior_weights"]], multivariate, weighted))
            )
          }
        }else{
          models <- c(
            models,
            list(.make_model(effect, heterogeneity, bias, NULL, effect[["prior_weights"]] * heterogeneity[["prior_weights"]] * bias[["prior_weights"]], multivariate, weighted))
          )
        }
      }
    }
  }

  return(models)
}
.make_model      <- function(prior_effect, prior_heterogeneity, prior_bias, prior_hierarchical, prior_weights, multivariate, weighted){

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
  # add 3 level structure only if there is heterogeneity
  if(!(prior_heterogeneity[["distribution"]] == "point" && prior_heterogeneity$parameters[["location"]] == 0) && !is.null(prior_hierarchical)){
    priors$rho <- prior_hierarchical
  }

  model <- list(
    priors            = priors,
    prior_weights     = prior_weights,
    prior_weights_set = prior_weights
  )
  class(model) <- "RoBMA.model"

  attr(model, "multivariate")  <- multivariate && !is.null(priors$rho)
  attr(model, "weighted")      <- weighted
  attr(model, "weighted_type") <- attr(weighted, "type")

  return(model)
}

.make_models.reg <- function(priors, multivariate, weighted){

  models_base <- .make_models(priors = priors, multivariate = multivariate, weighted = weighted)

  ### create grid of the models
  grid <- list(
    model = seq_along(models_base)
  )

  no_test <- attr(priors, "terms")[!attr(priors, "terms") %in% attr(priors, "terms_test")]
  to_test <- attr(priors, "terms")[ attr(priors, "terms") %in% attr(priors, "terms_test")]

  for(i in seq_along(no_test)){
    grid[[no_test[i]]] <- "alt"
  }
  for(i in seq_along(to_test)){
    grid[[to_test[i]]] <- c("null", "alt")
  }

  grid <- do.call(expand.grid, grid)

  if(nrow(grid) > 50)
    warning(sprintf("You are about to estimate %i models based on the model formula and prior specification.", nrow(grid)), immediate. = TRUE, call. = FALSE)

  ### create empty models objects for fitting
  models <- lapply(1:nrow(grid), function(i) .make_model.reg(models_base[[grid[i,1]]], grid[i,-1,drop=FALSE], priors))

  return(models)
}
.make_model.reg  <- function(model_base, grid_row, priors, multivariate, weighted){

  model_priors  <- model_base[["priors"]]
  prior_weights <- model_base[["prior_weights"]]
  terms         <- attr(priors, "terms")

  ### add priors for the terms
  model_priors[["terms"]] <- list()

  # rename mu to the intercept
  model_priors[["terms"]][["intercept"]] <- model_priors[["mu"]]
  model_priors[["mu"]] <- NULL

  terms_test <- NULL
  for(i in seq_along(terms)){
    model_priors[["terms"]][[terms[i]]]  <- priors[["terms"]][[terms[i]]][[grid_row[,terms[i]]]]
    prior_weights                        <- prior_weights * priors[["terms"]][[terms[i]]][[grid_row[,terms[i]]]]$prior_weights
    model_priors[["terms"]][[terms[i]]][["is_null"]] <- grid_row[,terms[i]] == "null"
    if(grid_row[,terms[i]] == "alt"){
      terms_test <- c(terms_test, terms[i])
    }
  }

  model <- list(
    priors       = model_priors,
    terms        = terms,
    terms_test   = terms_test,
    prior_weights     = prior_weights,
    prior_weights_set = prior_weights
  )

  class(model) <- "RoBMA.reg.model"
  attr(model, "multivariate")  <- attr(model_base, "multivariate")
  attr(model, "weighted")      <- attr(model_base, "weighted")
  attr(model, "weighted_type") <- attr(model_base, "weighted_type")

  return(model)
}
