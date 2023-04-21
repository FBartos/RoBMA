#' @title Estimate a Robust Bayesian Meta-Analysis Meta-Regression
#'
#' @description \code{RoBMA} is used to estimate a Robust Bayesian
#' Meta-Analysis. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @param formula a formula for the meta-regression model
#' @param test_predictors vector of predictor names that will be test
#' (i.e., assigned both the null and alternative prior distributions).
#' Defaults to \code{NULL}, no parameters are tested and only used for
#' adjustment.
#' @param priors named list of prior distributions for each predictor
#' (with names corresponding to the predictors). It allows users to
#' specify both the null and alternative hypothesis prior distributions
#' for each predictor by assigning the corresponding element of the named
#' list with another named list (with \code{"null"} and
#' \code{"alt"}).
#' @param prior_covariates a prior distributions for the regression parameter
#' of continuous covariates on the effect size under the alternative hypothesis
#' (unless set explicitly in \code{priors}). Defaults to a relatively wide normal
#' distribution \code{prior(distribution = "normal", parameters = list(mean = 0, sd = 0.25))}.
#' @param prior_covariates_null a prior distributions for the regression parameter
#' of continuous covariates on the effect size under the null hypothesis
#' (unless set explicitly in \code{priors}). Defaults to a no effect
#' \code{prior("spike",  parameters = list(location = 0))}.
#' @param prior_factors a prior distributions for the regression parameter
#' of categorical covariates on the effect size under the alternative hypothesis
#' (unless set explicitly in \code{priors}). Defaults to a relatively wide
#' multivariate normal distribution specifying differences from the mean contrasts
#' \code{prior_factor("mnormal", parameters = list(mean = 0, sd = 0.25), contrast = "meandif")}.
#' @param prior_factors_null a prior distributions for the regression parameter
#' of categorical covariates on the effect size under the null hypothesis
#' (unless set explicitly in \code{priors}). Defaults to a no effect
#' \code{prior("spike",  parameters = list(location = 0))}.
#' @param standardize_predictors whether continuous predictors should be standardized prior to
#' estimating the model. Defaults to \code{TRUE}.
#' @inheritParams RoBMA
#' @inheritParams combine_data
#'
#' @details See [RoBMA()] for more details.
#'
#' Note that these default prior distributions are relatively wide and more informed
#' prior distributions for testing for the presence of moderation should be considered.
#'
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return \code{RoBMA.reg} returns an object of class 'RoBMA.reg'.
#'
#' @seealso [RoBMA()] [summary.RoBMA()], [update.RoBMA()], [check_setup.reg()]
#' @export
RoBMA.reg <- function(
    formula, data, test_predictors = NULL, study_names = NULL, study_ids = NULL,
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

    # MCMC fitting settings
    chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
    autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

    # additional settings
    save = "all", seed = NULL, silent = TRUE, ...){


  dots         <- .RoBMA_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  object$data    <- .combine_data.reg(formula, data, standardize_predictors, transformation, study_names, study_ids)
  object$formula <- formula

  # switch between multivariate and weighted models
  if(attr(object$data[["outcome"]], "weighted"))
    .weighted_warning()

  if(.is_multivariate(object))
    .multivariate_warning()


  ### check MCMC settings
  object$fit_control        <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)


  ### prepare and check the settings
  object$priors     <- .check_and_list_priors.reg(
    priors = priors, data = object[["data"]], model_type = model_type, test_predictors = test_predictors, prior_scale = .transformation_var(prior_scale),
    priors_effect_null = priors_effect_null, priors_effect = priors_effect,
    priors_heterogeneity_null = priors_heterogeneity_null, priors_heterogeneity = priors_heterogeneity,
    priors_bias_null = priors_bias_null, priors_bias = priors_bias,
    priors_hierarchical_null = priors_hierarchical_null, priors_hierarchical = priors_hierarchical,
    prior_covariates_null = prior_covariates_null, prior_covariates = prior_covariates,
    prior_factors_null = prior_factors_null, prior_factors = prior_factors)
  object$models     <- .make_models.reg(object[["priors"]], .is_multivariate(object), .is_weighted(object), dots[["do_not_fit"]])


  ### additional information
  object$add_info <- .check_and_list_add_info(
    model_type             = model_type,
    predictors             = attr(object[["priors"]], "terms"),
    predictors_test        = attr(object[["priors"]], "terms_test"),
    prior_scale            = .transformation_var(prior_scale),
    output_scale           = .transformation_var(prior_scale),
    effect_measure         = attr(object$data[["outcome"]], "effect_measure"),
    effect_direction       = effect_direction,
    standardize_predictors = standardize_predictors,
    seed                   = seed,
    save                   = save,
    warnings               = NULL,
    errors                 = NULL
  )

  # the check requires the 'add_info' object already created
  object$add_info[["warnings"]] <- c(.check_effect_direction(object), .check_predictors_scaling(object))


  if(dots[["do_not_fit"]]){
    return(object)
  }


  ### fit the models and compute marginal likelihoods
  if(!object$fit_control[["parallel"]]){

    if(dots[["is_JASP"]]){
      .JASP_progress_bar_start(length(object[["models"]]))
    }

    for(i in seq_along(object[["models"]])){
      object$models[[i]] <- .fit_RoBMA_model(object, i)
      if(dots[["is_JASP"]]){
        .JASP_progress_bar_tick()
      }
    }

  }else{

    fitting_order <- .fitting_priority(object[["models"]])

    cl <- parallel::makePSOCKcluster(floor(RoBMA.get_option("max_cores") / object$fit_control[["chains"]]))
    parallel::clusterEvalQ(cl, {library("RoBMA")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::parLapplyLB(cl, fitting_order, .fit_RoBMA_model, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)

  }


  # create ensemble only if at least one model converged
  if(any(.get_model_convergence(object))){

    # balance probability of non-converged models
    if(object$convergence_checks[["balance_probability"]] && !all(.get_model_convergence(object))){
      object <- .balance_probability(object)
    }

    ### compute the model-space results
    object$models        <- BayesTools::models_inference(object[["models"]])
    object$RoBMA         <- .ensemble_inference(object)
    object$coefficients  <- .compute_coeficients(object[["RoBMA"]])
  }


  ### collect and print errors and warnings
  object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   .get_model_errors(object))
  object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], .get_model_warnings(object))
  .print_errors_and_warnings(object)


  ### remove model posteriors if asked to
  if(save == "min"){
    object <- .remove_model_posteriors(object)
    object <- .remove_model_margliks(object)
  }


  class(object) <- c("RoBMA", "RoBMA.reg")
  return(object)
}


.combine_data.reg <- function(formula, data, standardize_predictors, transformation, study_names, study_ids){

  if(!is.language(formula))
    stop("The 'formula' is not specidied as a formula.")
  if(!is.data.frame(data))
    stop("'data' must be an object of type data.frame.")
  BayesTools::check_bool(standardize_predictors, "standardize_predictors")


  ### deal with the effect sizes
  data_outcome <- combine_data(
    d      = if("d"      %in%  colnames(data)) data[,"d"],
    r      = if("r"      %in%  colnames(data)) data[,"r"],
    z      = if("z"      %in%  colnames(data)) data[,"z"],
    logOR  = if("logOR"  %in%  colnames(data)) data[,"logOR"],
    t      = if("t"      %in%  colnames(data)) data[,"t"],
    y      = if("y"      %in%  colnames(data)) data[,"y"],
    se     = if("se"     %in%  colnames(data)) data[,"se"],
    v      = if("v"      %in%  colnames(data)) data[,"v"],
    n      = if("n"      %in%  colnames(data)) data[,"n"],
    lCI    = if("lCI"    %in%  colnames(data)) data[,"lCI"],
    uCI    = if("uCI"    %in%  colnames(data)) data[,"uCI"],
    weight = if("weight" %in%  colnames(data)) data[,"weight"],
    study_names    = study_names,
    study_ids      = study_ids,
    transformation = transformation,
    return_all     = FALSE)

  ### obtain the predictors part
  data_predictors <- data[,!colnames(data) %in% c("d", "r", "z", "logOR", "t", "y", "se", "v", "n", "lCI", "uCI", "weight"), drop = FALSE]

  if(attr(stats::terms(formula), "response") == 1){
    formula[2] <- NULL
  }
  rhs             <- formula[c(1,2)]
  model_frame     <- stats::model.frame(rhs, data = data_predictors)

  # check that intercept is specified
  if(attr(attr(model_frame, "terms"),"intercept") == 0)
    stop("Intercept cannot be ommited from the model (you can set the coefficient to zero via 'priors_effect'). ")

  # change characters into factors
  for(i in seq_along(attr(attr(model_frame, "terms"), "dataClasses"))){
    if(attr(attr(model_frame, "terms"), "dataClasses")[[i]] == "character"){
      model_frame[,names(attr(attr(model_frame, "terms"), "dataClasses"))[i]] <-
        as.factor(model_frame[,names(attr(attr(model_frame, "terms"), "dataClasses"))[i]])
      attr(attr(model_frame, "terms"), "dataClasses")[[i]] <- "factor"
    }
  }


  model_frame     <- as.list(model_frame)
  if(length(model_frame) == 0){
    data_predictors <- list()
  }else{
    data_predictors <- model_frame[1:length(model_frame)]
  }
  attr(data_predictors, "variables")  <- attr(attr(model_frame, "terms"), "term.labels")[attr(attr(model_frame, "terms"), "order") == 1]
  attr(data_predictors, "terms")      <- attr(attr(model_frame, "terms"), "term.labels")
  attr(data_predictors, "terms_type") <- attr(attr(model_frame, "terms"), "dataClasses")


  # add additional information about the predictors
  data_predictors_info <- list()
  to_warn              <- NULL
  for(i in seq_along(data_predictors)){
    if(attr(data_predictors, "terms_type")[i] == "numeric"){

      data_predictors_info[[names(data_predictors)[i]]] <- list(
        type = "continuous",
        mean = mean(data_predictors[[names(data_predictors)[i]]]),
        sd   = stats::sd(data_predictors[[names(data_predictors)[i]]])
      )

      if(standardize_predictors){
        data_predictors[[names(data_predictors)[i]]] <- .pred_scale(data_predictors[[names(data_predictors)[i]]], data_predictors_info[[names(data_predictors)[i]]])
      }else if(RoBMA.get_option("check_scaling") && (abs(mean(data_predictors[[names(data_predictors)[i]]])) > 0.01 | abs(1 - stats::sd(data_predictors[[names(data_predictors)[i]]])) > 0.01)){
        to_warn <- c(to_warn, names(data_predictors)[i])
      }

    }else if(attr(data_predictors, "terms_type")[i] == "factor"){

      data_predictors_info[[names(data_predictors)[i]]] <- list(
        type    = "factor",
        default = levels(data_predictors[[names(data_predictors)[i]]])[1],
        levels  = levels(data_predictors[[names(data_predictors)[i]]])
      )

    }
  }
  attr(data_predictors, "variables_info") <- data_predictors_info


  output <- list(
    outcome    = data_outcome,
    predictors = data_predictors
  )

  # throw warnings and errors
  if(length(to_warn) > 0){
    scaling_warning <- paste0("The continuous predictors ", paste0("'", to_warn, "'", collapse = ", "), " are not scaled. Note that extra care need to be taken when specifying prior distributions for unscaled predictors.")
    warning(scaling_warning, immediate. = TRUE, call. = FALSE)
    warning("You can suppress this and following warnings via 'RoBMA.options(check_scaling = FALSE)'. To automatically rescale predictors set 'standardize_predictors = TRUE'.", immediate. = TRUE, call. = FALSE)
    attr(output, "warnings") <- scaling_warning
  }

  # check for reserved words
  if(any(attr(data_predictors, "terms") %in% .reserved_words()))
    stop(paste0("The following variable names are internally reserved keywords and cannot be used: ",
                paste0(" '", attr(data_predictors, "terms")[attr(data_predictors, "terms") %in% .reserved_words()], "' ", collapse = ", ")))

  return(output)
}
.pred_scale   <- function(x, predictor_info){
  (x - predictor_info[["mean"]]) / predictor_info[["sd"]]
}
.pred_unscale <- function(x, predictor_info){
  x * predictor_info[["sd"]] + predictor_info[["mean"]]
}
