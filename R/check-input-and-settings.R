#' @title Prints summary of \code{"RoBMA"} ensemble implied by the specified priors
#'
#' @description \code{check_setup} prints summary of \code{"RoBMA"} ensemble
#' implied by the specified prior distributions. It is useful for checking
#' the ensemble configuration prior to fitting all of the models.
#'
#' @inheritParams RoBMA
#' @param models should the models' details be printed.
#' @param silent do not print the results.
#'
#' @return \code{check_setup} invisibly returns list of summary tables.
#'
#' @seealso [check_setup.reg()] [RoBMA()]
#' @aliases check_setup.RoBMA
#' @export
check_setup <- function(
  model_type   = NULL,
  priors_effect         = prior(distribution = "normal",    parameters = list(mean  = 0, sd = 1)),
  priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15)),
  priors_bias           = list(
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
    prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
    prior_PET(distribution   = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),  prior_weights = 1/4),
    prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),  prior_weights = 1/4)
  ),
  priors_effect_null         = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0)),
  priors_bias_null           = prior_none(),
  priors_hierarchical        = prior("beta", parameters = list(alpha = 1, beta = 1)),
  priors_hierarchical_null   = NULL,
  models = FALSE, silent = FALSE){


  object <- RoBMA(
    # data specification
    d = 0, se = 1,

    # prior specification
    model_type   = model_type,
    priors_effect         = priors_effect,
    priors_heterogeneity  = priors_heterogeneity,
    priors_bias           = priors_bias,
    priors_effect_null         = priors_effect_null,
    priors_heterogeneity_null  = priors_heterogeneity_null,
    priors_bias_null           = priors_bias_null,
    priors_hierarchical        = priors_hierarchical,
    priors_hierarchical_null   = priors_hierarchical_null,

    # MCMC fitting settings
    do_not_fit = TRUE)

  object <- list()
  object$priors   <- .check_and_list_priors(model_type, priors_effect_null, priors_effect, priors_heterogeneity_null, priors_heterogeneity, priors_bias_null, priors_bias, priors_hierarchical, priors_hierarchical_null, "d")
  object$models   <- .make_models(object[["priors"]], multivariate = FALSE, weighted = FALSE)

  ### model types overview
  effect         <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "bias"))

  # obtain the parameter types
  weightfunctions <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # number of model types
  n_models    <- c(
    mu    = sum(effect),
    tau   = sum(heterogeneity),
    omega = sum(bias)
  )

  # extract model weights
  prior_weights   <- sapply(object$models, function(m)m$prior_weights)
  # standardize model weights
  prior_weights   <- prior_weights / sum(prior_weights)
  # conditional model weights
  models_prior <- c(
    mu    <- sum(prior_weights[effect]),
    tau   <- sum(prior_weights[heterogeneity]),
    omega <- sum(prior_weights[bias])
  )

  # create overview table
  components <- data.frame(
    "models"     = n_models,
    "prior_prob" = models_prior
  )
  rownames(components) <- c("Effect", "Heterogeneity", "Bias")

  class(components)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(components))
  attr(components, "type")      <- c("n_models", "prior_prob")
  attr(components, "rownames")  <- TRUE
  attr(components, "n_models")  <- length(object$models)
  attr(components, "title")     <- "Components summary:"
  attr(components, "footnotes") <- NULL
  attr(components, "warnings")  <- NULL

  object$components <- components

  ### model details
  if(models){
    priors_effect        <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$mu, silent = TRUE))
    priors_heterogeneity <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_bias          <- sapply(1:length(object$models), function(i){
      if(weightfunctions[i]){
        print(object$models[[i]]$priors$omega, silent = TRUE)
      }else if(PET[i]){
        print(object$models[[i]]$priors$PET, silent = TRUE)
      }else if(PEESE[i]){
        print(object$models[[i]]$priors$PEESE, silent = TRUE)
      }else{
        ""
      }
    })
    prior_weights  <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_weights)
    prior_prob     <- prior_weights / sum(prior_weights)

    summary <- cbind.data.frame(
      "Model"         = 1:length(object$models),
      "Effect"        = priors_effect,
      "Heterogeneity" = priors_heterogeneity,
      "Bias"          = priors_bias,
      "prior_prob"    = prior_prob
    )
    class(summary)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(summary))
    attr(summary, "type")      <- c("integer", rep("prior", 3), "prior_prob")
    attr(summary, "rownames")  <- FALSE
    attr(summary, "title")     <- "Models overview:"
    attr(summary, "footnotes") <- NULL
    attr(summary, "warnings")  <- NULL

    object$summary <- summary
  }


  if(!silent){
    cat("Robust Bayesian meta-analysis (set-up)\n")
    print(components, quote = FALSE, right = TRUE)

    if(models){
      cat("\n")
      print(summary, quote = FALSE, right = TRUE)
    }
  }

  return(invisible(object))
}

#' @rdname check_setup
check_setup.RoBMA <- check_setup

#' @title Prints summary of \code{"RoBMA.reg"} ensemble implied by the specified priors
#' and formula
#'
#' @description \code{check_setup} prints summary of \code{"RoBMA.reg"} ensemble
#' implied by the specified prior distributions. It is useful for checking
#' the ensemble configuration prior to fitting all of the models.
#'
#' @inheritParams check_setup
#' @inheritParams RoBMA.reg
#'
#' @return \code{check_setup.reg} invisibly returns list of summary tables.
#'
#' @seealso [check_setup()] [RoBMA.reg()]
#' @aliases check_setup.RoBMA.reg
#' @export
check_setup.reg <- function(
    formula, data, test_predictors = TRUE, study_names = NULL, study_ids = NULL,
    transformation     = if(any(colnames(data) != "y")) "fishers_z" else "none",
    prior_scale        = if(any(colnames(data) != "y")) "cohens_d"  else "none",
    standardize_predictors = TRUE,
    effect_direction       = "positive",

    # prior specification
    priors       = NULL,
    model_type   = NULL,

    priors_effect         = prior(distribution = "normal",    parameters = list(mean  = 0, sd = 1)),
    priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15)),
    priors_bias           = list(
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
      prior_PET(distribution   = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),  prior_weights = 1/4),
      prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),  prior_weights = 1/4)
    ),
    priors_effect_null         = prior(distribution = "point", parameters = list(location = 0)),
    priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0)),
    priors_bias_null           = prior_none(),
    priors_hierarchical        = prior("beta", parameters = list(alpha = 1, beta = 1)),
    priors_hierarchical_null   = NULL,

    prior_covariates       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
    prior_covariates_null  = prior("spike",  parameters = list(location = 0)),
    prior_factors          = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25), contrast = "meandif"),
    prior_factors_null     = prior("spike",  parameters = list(location = 0)),
    models = FALSE, silent = FALSE, ...){


  # use 'do_not_fit' with RoBMA.reg
  object <- RoBMA.reg(
    ### passed arguments
    formula = formula, data = data, test_predictors = test_predictors, study_names = study_names, study_ids = study_ids,
    transformation     = transformation,
    prior_scale        = prior_scale,
    standardize_predictors = standardize_predictors,
    effect_direction       = effect_direction,

    priors       = priors,
    model_type   = model_type,

    priors_effect         = priors_effect,
    priors_heterogeneity  = priors_heterogeneity,
    priors_bias           = priors_bias,
    priors_effect_null         = priors_effect_null,
    priors_heterogeneity_null  = priors_heterogeneity_null,
    priors_bias_null           = priors_bias_null,
    priors_hierarchical        = priors_hierarchical,
    priors_hierarchical_null   = priors_hierarchical_null,

    prior_covariates       = prior_covariates,
    prior_covariates_null  = prior_covariates_null,
    prior_factors          = prior_factors,
    prior_factors_null     = prior_factors_null,

    do_not_fit = TRUE)


  ### Components summary
  effect         <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "bias"))

  # obtain the parameter types
  weightfunctions <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # number of model types
  n_models    <- c(
    mu    = sum(effect),
    tau   = sum(heterogeneity),
    omega = sum(bias)
  )

  # extract model weights
  prior_weights   <- sapply(object$models, function(m) m$prior_weights)
  # standardize model weights
  prior_weights   <- prior_weights / sum(prior_weights)
  # conditional model weights
  models_prior <- c(
    mu    = sum(prior_weights[effect]),
    tau   = sum(prior_weights[heterogeneity]),
    omega = sum(prior_weights[bias])
  )

  # create overview table
  components.tab <- data.frame(
    "models"     = n_models,
    "prior_prob" = models_prior
  )
  rownames(components.tab) <- c("Effect", "Heterogeneity", "Bias")

  class(components.tab)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(components.tab))
  attr(components.tab, "type")      <- c("n_models", "prior_prob")
  attr(components.tab, "rownames")  <- TRUE
  attr(components.tab, "n_models")  <- length(object$models)
  attr(components.tab, "title")     <- "Components summary:"
  attr(components.tab, "footnotes") <- NULL
  attr(components.tab, "warnings")  <- NULL


  ### Meta-regression components summary
  model_predictors      <- lapply(object$models, function(model) model[["terms"]])
  model_predictors_test <- lapply(object$models, function(model) model[["terms_test"]])

  predictors      <- object$add_info[["predictors"]]
  predictors_test <- object$add_info[["predictors_test"]]

  # define inference options
  components_predictors      <- NULL
  parameters_predictors      <- "mu_intercept"
  components_predictors_null <- list()
  parameters_predictors_null <- list("mu_intercept" = !effect)

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
    parameters_predictors_null[[.BayesTools_parameter_name(predictors[i])]] <-
      sapply(model_predictors_test, function(x) if(length(x) == 0) TRUE else !(predictors[i] %in% x))
  }

  # create overview table
  if(length(components_predictors_null) > 0){

    components_predictors.tab <- data.frame(
      "models"     = sapply(components_predictors_null, sum),
      "prior_prob" = sapply(seq_along(components_predictors_null), function(i) sum(prior_weights[components_predictors_null[[i]]]))
    )
    rownames(components_predictors.tab) <- .output_parameter_names(components_predictors)

    class(components_predictors.tab)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(components_predictors.tab))
    attr(components_predictors.tab, "type")      <- c("n_models", "prior_prob")
    attr(components_predictors.tab, "rownames")  <- TRUE
    attr(components_predictors.tab, "n_models")  <- length(object$models)
    attr(components_predictors.tab, "title")     <- "Meta-regression components summary:"
    attr(components_predictors.tab, "footnotes") <- NULL
    attr(components_predictors.tab, "warnings")  <- NULL

  }else{

    components_predictors.tab <- BayesTools::ensemble_inference_empty_table(title = "Meta-regression components summary:")
    components_predictors.tab <- BayesTools::remove_column(components_predictors.tab, 4)
    components_predictors.tab <- BayesTools::remove_column(components_predictors.tab, 3)

  }


  ### store summaries in the object
  object$components            <- components.tab
  object$components_predictors <- components_predictors.tab


  ### model details
  if(models){
    priors_heterogeneity <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_bias          <- sapply(1:length(object$models), function(i){
      if(weightfunctions[i]){
        print(object$models[[i]]$priors$omega, silent = TRUE)
      }else if(PET[i]){
        print(object$models[[i]]$priors$PET, silent = TRUE)
      }else if(PEESE[i]){
        print(object$models[[i]]$priors$PEESE, silent = TRUE)
      }else{
        ""
      }
    })
    prior_weights  <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_weights)
    prior_prob     <- prior_weights / sum(prior_weights)

    summary <- data.frame("Model" = 1:length(object$models))
    for(p in seq_along(parameters_predictors)){
      summary <- cbind(summary, sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$terms[[.output_parameter_names(parameters_predictors[p])]], silent = TRUE)))
      colnames(summary)[p+1] <- .output_parameter_names(parameters_predictors[i])
    }
    summary <- cbind(
      summary,
      "Heterogeneity" = priors_heterogeneity,
      "Bias"          = priors_bias,
      "prior_prob"    = prior_prob
    )
    class(summary)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(summary))
    attr(summary, "type")      <- c("integer", rep("prior", 2 + length(parameters_predictors)), "prior_prob")
    attr(summary, "rownames")  <- FALSE
    attr(summary, "title")     <- "Models overview:"
    attr(summary, "footnotes") <- NULL
    attr(summary, "warnings")  <- NULL

    object$summary <- summary
  }


  if(!silent){
    cat("Robust Bayesian meta-regression (set-up)\n")
    print(components.tab, quote = FALSE, right = TRUE)

    cat("\n")
    print(components_predictors.tab, quote = FALSE, right = TRUE)

    if(models){
      cat("\n")
      print(summary, quote = FALSE, right = TRUE)
    }
  }

  return(invisible(object))
}

#' @rdname check_setup.reg
check_setup.RoBMA.reg <- check_setup.reg


#' @title Prints summary of \code{"BiBMA.reg"} ensemble implied by the specified priors
#' and formula
#'
#' @description \code{check_setup} prints summary of \code{"RoBMA.reg"} ensemble
#' implied by the specified prior distributions. It is useful for checking
#' the ensemble configuration prior to fitting all of the models.
#'
#' @inheritParams BiBMA
#' @inheritParams check_setup
#'
#' @return \code{check_setup.reg} invisibly returns list of summary tables.
#'
#' @seealso [check_setup()] [BiBMA()]
#' @export
check_setup.BiBMA <- function(
    # prior specification
    priors_effect         = prior(distribution = "student",   parameters = list(location = 0, scale = 0.58, df = 4)),
    priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1.77, scale = 0.55)),

    priors_effect_null         = prior(distribution = "point", parameters = list(location = 0)),
    priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0)),

    priors_baseline        = NULL,
    priors_baseline_null   = prior_factor("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent"),

    models = FALSE, silent = FALSE, ...){

  object <- BiBMA(
    # add empty data (they do not determine model structure)
    x1 = 0, x2 = 0, n1 = 1, n2 = 1,

    # prior specification -- passed from the user
    priors_effect         = priors_effect,
    priors_heterogeneity  = priors_heterogeneity,

    priors_effect_null         = priors_effect_null,
    priors_heterogeneity_null  = priors_heterogeneity_null,

    priors_baseline        = priors_baseline,
    priors_baseline_null   = priors_baseline_null,

    # return model structure only
    do_not_fit = TRUE
  )

  ### model types overview
  effect         <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  baseline       <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "baseline"))


  # number of model types
  n_models    <- c(
    mu  = sum(effect),
    tau = sum(heterogeneity)
  )
  if(any(baseline)){
    n_models <- c(n_models, baseline = sum(baseline))
  }

  # extract model weights
  prior_weights   <- sapply(object$models, function(m)m$prior_weights)
  # standardize model weights
  prior_weights   <- prior_weights / sum(prior_weights)
  # conditional model weights
  models_prior <- c(
    mu    = sum(prior_weights[effect]),
    tau   = sum(prior_weights[heterogeneity])
  )
  if(any(baseline)){
    models_prior <- c(models_prior, baseline = sum(prior_weights[baseline]))
  }

  # create overview table
  components <- data.frame(
    "models"     = n_models,
    "prior_prob" = models_prior
  )
  rownames(components) <- c("Effect", "Heterogeneity", if(any(baseline)) "Baseline")

  class(components)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(components))
  attr(components, "type")      <- c("n_models", "prior_prob")
  attr(components, "rownames")  <- TRUE
  attr(components, "n_models")  <- length(object$models)
  attr(components, "title")     <- "Components summary:"
  attr(components, "footnotes") <- NULL
  attr(components, "warnings")  <- NULL

  object$components <- components

  ### model details
  if(models){
    priors_effect        <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$mu, silent = TRUE))
    priors_heterogeneity <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_baseline      <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$pi, silent = TRUE))

    prior_weights <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_weights)
    prior_prob    <- prior_weights / sum(prior_weights)

    summary <- cbind.data.frame(
      "Model"         = 1:length(object$models),
      "Effect"        = priors_effect,
      "Heterogeneity" = priors_heterogeneity,
      "Baseline"      = priors_baseline,
      "prior_prob"    = prior_prob
    )
    class(summary)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(summary))
    attr(summary, "type")      <- c("integer", rep("prior", 3), "prior_prob")
    attr(summary, "rownames")  <- FALSE
    attr(summary, "title")     <- "Models overview:"
    attr(summary, "footnotes") <- NULL
    attr(summary, "warnings")  <- NULL

    object$summary <- summary
  }


  if(!silent){
    cat("Bayesian model-averaged meta-analysis (binomial model) (set-up)\n")
    print(components, quote = FALSE, right = TRUE)

    if(models){
      cat("\n")
      print(summary, quote = FALSE, right = TRUE)
    }
  }

  return(invisible(object))
}

#' @title Prints summary of \code{"RoBMA.reg"} ensemble implied by the specified priors
#' and formula
#'
#' @description \code{check_setup} prints summary of \code{"RoBMA.reg"} ensemble
#' implied by the specified prior distributions. It is useful for checking
#' the ensemble configuration prior to fitting all of the models.
#'
#' @inheritParams check_setup
#' @inheritParams RoBMA.reg
#'
#' @return \code{check_setup.reg} invisibly returns list of summary tables.
#'
#' @seealso [check_setup()] [RoBMA.reg()]
#' @export
check_setup.reg <- function(
    formula, data, test_predictors = TRUE, study_names = NULL, study_ids = NULL,
    transformation     = if(any(colnames(data) != "y")) "fishers_z" else "none",
    prior_scale        = if(any(colnames(data) != "y")) "cohens_d"  else "none",
    standardize_predictors = TRUE,
    effect_direction       = "positive",

    # prior specification
    priors       = NULL,
    model_type   = NULL,

    priors_effect         = prior(distribution = "normal",    parameters = list(mean  = 0, sd = 1)),
    priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15)),
    priors_bias           = list(
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),       prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),             prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.025, 0.05)),      prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.5)),        prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
      prior_PET(distribution   = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),  prior_weights = 1/4),
      prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),  prior_weights = 1/4)
    ),
    priors_effect_null         = prior(distribution = "point", parameters = list(location = 0)),
    priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0)),
    priors_bias_null           = prior_none(),
    priors_hierarchical        = prior("beta", parameters = list(alpha = 1, beta = 1)),
    priors_hierarchical_null   = NULL,

    prior_covariates       = prior("normal", parameters = list(mean = 0, sd = 0.25)),
    prior_covariates_null  = prior("spike",  parameters = list(location = 0)),
    prior_factors          = prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25), contrast = "meandif"),
    prior_factors_null     = prior("spike",  parameters = list(location = 0)),
    models = FALSE, silent = FALSE, ...){


  # use 'do_not_fit' with RoBMA.reg
  object <- RoBMA.reg(
    ### passed arguments
    formula = formula, data = data, test_predictors = test_predictors, study_names = study_names, study_ids = study_ids,
    transformation     = transformation,
    prior_scale        = prior_scale,
    standardize_predictors = standardize_predictors,
    effect_direction       = effect_direction,

    priors       = priors,
    model_type   = model_type,

    priors_effect         = priors_effect,
    priors_heterogeneity  = priors_heterogeneity,
    priors_bias           = priors_bias,
    priors_effect_null         = priors_effect_null,
    priors_heterogeneity_null  = priors_heterogeneity_null,
    priors_bias_null           = priors_bias_null,
    priors_hierarchical        = priors_hierarchical,
    priors_hierarchical_null   = priors_hierarchical_null,

    prior_covariates       = prior_covariates,
    prior_covariates_null  = prior_covariates_null,
    prior_factors          = prior_factors,
    prior_factors_null     = prior_factors_null,

    do_not_fit = TRUE,

    ### fitting settings defaults
    # MCMC fitting settings
    chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
    autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

    # additional settings
    save = "all", seed = NULL, silent = TRUE)


  ### Components summary
  effect         <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "bias"))

  # obtain the parameter types
  weightfunctions <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # number of model types
  n_models    <- c(
    mu    = sum(effect),
    tau   = sum(heterogeneity),
    omega = sum(bias)
  )

  # extract model weights
  prior_weights   <- sapply(object$models, function(m) m$prior_weights)
  # standardize model weights
  prior_weights   <- prior_weights / sum(prior_weights)
  # conditional model weights
  models_prior <- c(
    mu    <- sum(prior_weights[effect]),
    tau   <- sum(prior_weights[heterogeneity]),
    omega <- sum(prior_weights[bias])
  )

  # create overview table
  components.tab <- data.frame(
    "models"     = n_models,
    "prior_prob" = models_prior
  )
  rownames(components.tab) <- c("Effect", "Heterogeneity", "Bias")

  class(components.tab)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(components.tab))
  attr(components.tab, "type")      <- c("n_models", "prior_prob")
  attr(components.tab, "rownames")  <- TRUE
  attr(components.tab, "n_models")  <- length(object$models)
  attr(components.tab, "title")     <- "Components summary:"
  attr(components.tab, "footnotes") <- NULL
  attr(components.tab, "warnings")  <- NULL


  ### Meta-regression components summary
  model_predictors      <- lapply(object$models, function(model) model[["terms"]])
  model_predictors_test <- lapply(object$models, function(model) model[["terms_test"]])

  predictors      <- object$add_info[["predictors"]]
  predictors_test <- object$add_info[["predictors_test"]]

  # define inference options
  components_predictors      <- NULL
  parameters_predictors      <- "mu_intercept"
  components_predictors_null <- list()
  parameters_predictors_null <- list("mu_intercept" = !effect)

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
    parameters_predictors_null[[.BayesTools_parameter_name(predictors[i])]] <-
      sapply(model_predictors_test, function(x) if(length(x) == 0) TRUE else !(predictors[i] %in% x))
  }

  # create overview table
  if(length(components_predictors_null) > 0){

    components_predictors.tab <- data.frame(
      "models"     = sapply(components_predictors_null, sum),
      "prior_prob" = sapply(seq_along(components_predictors_null), function(i) sum(prior_weights[components_predictors_null[[i]]]))
    )
    rownames(components_predictors.tab) <- .output_parameter_names(components_predictors)

    class(components_predictors.tab)             <- c("BayesTools_table", "BayesTools_ensemble_inference", class(components_predictors.tab))
    attr(components_predictors.tab, "type")      <- c("n_models", "prior_prob")
    attr(components_predictors.tab, "rownames")  <- TRUE
    attr(components_predictors.tab, "n_models")  <- length(object$models)
    attr(components_predictors.tab, "title")     <- "Meta-regression components summary:"
    attr(components_predictors.tab, "footnotes") <- NULL
    attr(components_predictors.tab, "warnings")  <- NULL

  }else{

    components_predictors.tab <- BayesTools::ensemble_inference_empty_table(title = "Meta-regression components summary:")
    components_predictors.tab <- BayesTools::remove_column(components_predictors.tab, 4)
    components_predictors.tab <- BayesTools::remove_column(components_predictors.tab, 3)

  }


  ### store summaries in the object
  object$components            <- components.tab
  object$components_predictors <- components_predictors.tab


  ### model details
  if(models){
    priors_heterogeneity <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_bias          <- sapply(1:length(object$models), function(i){
      if(weightfunctions[i]){
        print(object$models[[i]]$priors$omega, silent = TRUE)
      }else if(PET[i]){
        print(object$models[[i]]$priors$PET, silent = TRUE)
      }else if(PEESE[i]){
        print(object$models[[i]]$priors$PEESE, silent = TRUE)
      }else{
        ""
      }
    })
    prior_weights  <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_weights)
    prior_prob     <- prior_weights / sum(prior_weights)

    summary <- data.frame("Model" = 1:length(object$models))
    for(p in seq_along(parameters_predictors)){
      summary <- cbind(summary, sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$terms[[.output_parameter_names(parameters_predictors[p])]], silent = TRUE)))
      colnames(summary)[p+1] <- .output_parameter_names(parameters_predictors[i])
    }
    summary <- cbind(
      summary,
      "Heterogeneity" = priors_heterogeneity,
      "Bias"          = priors_bias,
      "prior_prob"    = prior_prob
    )
    class(summary)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(summary))
    attr(summary, "type")      <- c("integer", rep("prior", 2 + length(parameters_predictors)), "prior_prob")
    attr(summary, "rownames")  <- FALSE
    attr(summary, "title")     <- "Models overview:"
    attr(summary, "footnotes") <- NULL
    attr(summary, "warnings")  <- NULL

    object$summary <- summary
  }


  if(!silent){
    cat("Robust Bayesian meta-regression (set-up)\n")
    print(components.tab, quote = FALSE, right = TRUE)

    cat("\n")
    print(components_predictors.tab, quote = FALSE, right = TRUE)

    if(models){
      cat("\n")
      print(summary, quote = FALSE, right = TRUE)
    }
  }

  return(invisible(object))
}

#' @title Control MCMC fitting process
#'
#' @description Controls settings for the autofit
#' process of the MCMC JAGS sampler (specifies termination
#' criteria), and values for the convergence checks.
#'
#' @param max_Rhat maximum value of the R-hat diagnostic.
#' Defaults to \code{1.05}.
#' @param min_ESS minimum estimated sample size.
#' Defaults to \code{500}.
#' @param max_error maximum value of the MCMC error.
#' Defaults to \code{NULL}. Be aware that PEESE publication bias
#' adjustment can have estimates on different scale than the rest of
#' the output, resulting in relatively large max MCMC error.
#' @param max_SD_error maximum value of the proportion of MCMC error
#' of the estimated SD of the parameter.
#' Defaults to \code{NULL}.
#' @param max_time list with the time and unit specifying the maximum
#' autofitting process per model. Passed to \link[base]{difftime} function
#' (possible units are \code{"secs", "mins", "hours", "days", "weeks", "years"}).
#' Defaults to \code{list(time = 60, unit = "mins")}.
#' @param sample_extend number of samples to extend the fitting process if
#' the criteria are not satisfied.
#' Defaults to \code{1000}.
#' @param remove_failed whether models not satisfying the convergence checks should
#' be removed from the inference. Defaults to \code{FALSE} - only a warning is raised.
#' @param balance_probability whether prior model probability should be balanced
#' across the combinations of models with the same H0/H1 for effect / heterogeneity / bias
#' in the case of non-convergence. Defaults to \code{TRUE}.
#' @param restarts number of times new initial values should be generated in case a
#' model fails to initialize. Defaults to \code{10}.
#'
#'
#' @return \code{set_autofit_control} returns a list of autofit control settings
#' and \code{set_convergence_checks} returns a list of convergence checks settings.
#'
#' @export set_autofit_control
#' @export set_convergence_checks
#' @name RoBMA_control
#' @aliases set_autofit_control, set_convergence_checks
#'
#' @seealso [RoBMA], [update.RoBMA]
NULL

#' @rdname RoBMA_control
set_autofit_control     <- function(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10){

  autofit_settings <- list(
    max_Rhat      = max_Rhat,
    min_ESS       = min_ESS,
    max_error     = max_error,
    max_SD_error  = max_SD_error,
    max_time      = max_time,
    sample_extend = sample_extend,
    restarts      = restarts
  )
  autofit_settings <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_settings, call = "Checking 'autofit_control':\n\t")

  return(autofit_settings)
}
#' @rdname RoBMA_control
set_convergence_checks  <- function(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL, remove_failed = FALSE, balance_probability = TRUE){

  convergence_checks <- list(
    max_Rhat            = max_Rhat,
    min_ESS             = min_ESS,
    max_error           = max_error,
    max_SD_error        = max_SD_error,
    remove_failed       = remove_failed,
    balance_probability = balance_probability
  )
  # allows NULL arguments so it can be used in this way too
  convergence_checks <- .check_and_list_convergence_checks(convergence_checks)

  return(convergence_checks)
}



.update_fit_control     <- function(old_fit_control, chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed){

  if(is.null(chains)){
    chains <- old_fit_control[["chains"]]
  }
  if(is.null(adapt)){
    adapt  <- old_fit_control[["adapt"]]
  }
  if(is.null(burnin)){
    burnin <- old_fit_control[["burnin"]]
  }
  if(is.null(sample)){
    sample <- old_fit_control[["sample"]]
  }
  if(is.null(thin)){
    thin  <- old_fit_control[["thin"]]
  }
  if(is.null(autofit)){
    autofit  <- old_fit_control[["autofit"]]
  }
  if(is.null(parallel)){
    parallel <- old_fit_control[["parallel"]]
  }
  if(is.null(silent)){
    silent <- old_fit_control[["silent"]]
  }
  if(is.null(seed)){
    seed   <- old_fit_control[["seed"]]
  }

  new_fit_control <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)

  return(new_fit_control)
}
.update_autofit_control <- function(old_autofit_control, autofit_control){

  if(!is.null(autofit_control[["max_Rhat"]])){
    max_Rhat <- autofit_control[["max_Rhat"]]
  }else{
    max_Rhat <- old_autofit_control[["max_Rhat"]]
  }
  if(!is.null(autofit_control[["min_ESS"]])){
    min_ESS <- autofit_control[["min_ESS"]]
  }else{
    min_ESS <- old_autofit_control[["min_ESS"]]
  }
  if(!is.null(autofit_control[["max_error"]])){
    max_error <- autofit_control[["max_error"]]
  }else{
    max_error <- old_autofit_control[["max_error"]]
  }
  if(!is.null(autofit_control[["max_SD_error"]])){
    max_SD_error <- autofit_control[["max_SD_error"]]
  }else{
    max_SD_error <- old_autofit_control[["max_SD_error"]]
  }
  if(!is.null(autofit_control[["max_time"]])){
    max_time <- autofit_control[["max_time"]]
  }else{
    max_time <- old_autofit_control[["max_time"]]
  }
  if(!is.null(autofit_control[["sample_extend"]])){
    sample_extend <- autofit_control[["sample_extend"]]
  }else{
    sample_extend <- old_autofit_control[["sample_extend"]]
  }
  if(!is.null(autofit_control[["restarts"]])){
    restarts <- autofit_control[["restarts"]]
  }else{
    restarts <- old_autofit_control[["restarts"]]
  }

  new_autofit_control <- set_autofit_control(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, max_time = max_time, sample_extend = sample_extend, restarts = restarts)
  new_autofit_control <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = new_autofit_control)

  return(new_autofit_control)
}
.update_convergence_checks <- function(old_convergence_checks, convergence_checks){

  if(!is.null(convergence_checks[["max_Rhat"]])){
    max_Rhat <- convergence_checks[["max_Rhat"]]
  }else{
    max_Rhat <- old_convergence_checks[["max_Rhat"]]
  }
  if(!is.null(convergence_checks[["min_ESS"]])){
    min_ESS <- convergence_checks[["min_ESS"]]
  }else{
    min_ESS <- old_convergence_checks[["min_ESS"]]
  }
  if(!is.null(convergence_checks[["max_error"]])){
    max_error <- convergence_checks[["max_error"]]
  }else{
    max_error <- old_convergence_checks[["max_error"]]
  }
  if(!is.null(convergence_checks[["max_SD_error"]])){
    max_SD_error <- convergence_checks[["max_SD_error"]]
  }else{
    max_SD_error <- old_convergence_checks[["max_SD_error"]]
  }
  if(!is.null(convergence_checks[["remove_failed"]])){
    remove_failed <- convergence_checks[["remove_failed"]]
  }else{
    remove_failed <- old_convergence_checks[["remove_failed"]]
  }
  if(!is.null(convergence_checks[["balance_probability"]])){
    balance_probability <- convergence_checks[["balance_probability"]]
  }else{
    balance_probability <- old_convergence_checks[["balance_probability"]]
  }

  new_convergence_checks <- set_convergence_checks(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, remove_failed = remove_failed, balance_probability = balance_probability)
  new_convergence_checks <- .check_and_list_convergence_checks(new_convergence_checks)
}


.check_and_list_convergence_checks <- function(convergence_checks){

  remove_failed       <- convergence_checks[["remove_failed"]]
  balance_probability <- convergence_checks[["balance_probability"]]
  convergence_checks["remove_failed"]       <- NULL
  convergence_checks["balance_probability"] <- NULL
  convergence_checks <- BayesTools::JAGS_check_and_list_autofit_settings(convergence_checks, skip_sample_extend = TRUE, call = "Checking 'convergence_checks':\n\t")

  BayesTools::check_bool(remove_failed,       "remove_failed",       call = "Checking 'convergence_checks':\n\t")
  BayesTools::check_bool(balance_probability, "balance_probability", call = "Checking 'convergence_checks':\n\t")
  convergence_checks[["remove_failed"]]       <- remove_failed
  convergence_checks[["balance_probability"]] <- balance_probability
  return(convergence_checks)
}
.check_and_list_add_info  <- function(model_type, predictors = NULL, predictors_test = NULL, prior_scale, output_scale, effect_measure, effect_direction, algorithm, standardize_predictors = NULL, seed, save, warnings, errors){

  BayesTools::check_char(effect_direction, "effect_direction", allow_values = c("positive", "negative"))
  BayesTools::check_real(seed, "seed", allow_NULL = TRUE)
  BayesTools::check_char(save, "save", allow_values = c("min", "all"))
  model_type <- .check_and_set_model_type(model_type, prior_scale)
  BayesTools::check_char(predictors, "predictors", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_char(predictors_test, "predictors_test", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(standardize_predictors, "standardize_predictors", allow_NULL = TRUE)
  BayesTools::check_char(algorithm, "algorithm", allow_values = c("bridge", "ss"))

  if((prior_scale == "y" & effect_measure != "y") | (prior_scale != "y" & effect_measure == "y"))
    stop("Prior / effect size transformations are not available for unstandardized effect sizes.", call. = FALSE)

  add_info <- list(
    model_type             = model_type,
    predictors             = predictors,
    predictors_test        = predictors_test,
    prior_scale            = prior_scale,
    output_scale           = output_scale,
    effect_measure         = effect_measure,
    effect_direction       = effect_direction,
    standardize_predictors = standardize_predictors,
    algorithm              = algorithm,
    seed                   = seed,
    save                   = save,
    warnings               = warnings,
    errors                 = errors,
    version                = utils::packageVersion("RoBMA")
  )

  return(add_info)
}
.check_and_set_model_type <- function(model_type, prior_scale){

  if(length(model_type) != 0){
    model_type <- tolower(model_type)
    BayesTools::check_char(model_type, "model_type", allow_values = c("psma", "2w", "pp"))
    if(prior_scale != "d")
      stop("The pre-specified models can be used only with prior distributions specified on Cohen's d.", call. = FALSE)
  }

  return(model_type)
}
.check_effect_direction   <- function(object){

  warnings <- NULL

  if(!is.null(object$data[["outcome"]])){
    data <- object$data[["outcome"]]
  }else{
    data <- object[["data"]]
  }

  # check whether majority of effect sizes are in expected direction. throw warning if not.
  if(any(sapply(object$priors[["bias"]], function(p) p[["distribution"]]) == "one.sided") |
     any(grepl("PET",   sapply(object$priors[["bias"]], function(p) p[["distribution"]])))  |
     any(grepl("PEESE", sapply(object$priors[["bias"]], function(p) p[["distribution"]])))){
    if(stats::median(data[["y"]]) > 0 & object$add_info[["effect_direction"]] == "negative" |
       stats::median(data[["y"]]) < 0 & object$add_info[["effect_direction"]] == "positive"){
      warnings <- "The majority of effect sizes is in the oposite direction than expected. The direction of effect sizes is important for the one-sided weight functions. Please, check the 'effect_direction' argument in 'RoBMA' fitting function."
    }
  }

  # the actual effect size direction changes are done prior and after fitting using the '.fit_data' and '.change_direction' functions

  return(warnings)
}
.check_predictors_scaling <- function(object){

  warnings   <- NULL
  predictors <- object[["data"]][["predictors"]]

  # check that all continuous predictors are centered and scale
  for(i in seq_along(predictors)){

    if(attr(predictors, "variables_info")[[i]][["type"]] == "continuous"){
      if(!(isTRUE(all.equal(mean(predictors[[i]]), 0)) && isTRUE(all.equal(stats::sd(predictors[[i]]), 1)))){
        warnings <- c(warnings, paste0("The continuous predictor '",names(predictors[i]),"' is not standardized. Be cafefull about the specified prior distribution and hypothesis test."))
      }
    }

  }

  return(warnings)
}

# some functions for the JASP implementation
.RoBMA_collect_dots      <- function(...){

  dots <- list(...)

  known_dots <- c("is_JASP", "is_JASP_prefix", "weighted", "do_not_fit", "weighted_type")
  if(any(!names(dots) %in% known_dots))
    stop(paste0("Uknown arguments to 'RoBMA': ", paste("'", names(dots)[!names(dots) %in% known_dots], "'", collapse = ", "), "."), call. = FALSE)

  if(is.null(dots[["is_JASP"]])){
    dots[["is_JASP"]] <- FALSE
  }else{
    dots[["is_JASP"]] <- TRUE
  }

  if(is.null(dots[["weighted"]])){
    dots[["weighted"]] <- FALSE
  }else{
    BayesTools::check_bool(dots[["weighted"]], "weighted")

    # select weight type
    if(is.null(dots[["weighted_type"]])){
      attr(dots[["weighted"]], "type") <- "inverse"
    }else{
      BayesTools::check_char(dots[["weighted_type"]], "weighted", allow_values = c("inverse", "inverse_sqrt", "custom"))
      attr(dots[["weighted"]], "type") <- dots[["weighted_type"]]
      dots[["weighted_type"]] <- NULL
    }
  }

  if(is.null(dots[["do_not_fit"]])){
    dots[["do_not_fit"]] <- FALSE
  }else{
    BayesTools::check_bool(dots[["do_not_fit"]], "do_not_fit")
  }

  return(dots)
}
.JASP_progress_bar_start <- function(n){
  eval(expr = parse(text = 'startProgressbar(n)'))
}
.JASP_progress_bar_tick  <- function(){
  eval(expr = parse(text = 'progressbarTick()'))
}
