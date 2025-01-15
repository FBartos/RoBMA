#' @title Estimate a Bayesian Model-Averaged Meta-Analysis of Binomial Data
#'
#' @description \code{BiBMA} estimate a binomial-normal Bayesian
#' model-averaged meta-analysis. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @param x1 a vector with the number of successes in the first group
#' @param x2 a vector with the number of successes in the second group
#' @param n1 a vector with the number of observations in the first group
#' @param n2 a vector with the number of observations in the second group
#' @param priors_baseline prior distributions for the alternative hypothesis about
#' intercepts (\code{pi}) of each study. Defaults to \code{NULL}.
#' @param priors_baseline_null prior distributions for the null hypothesis about
#' intercepts (\code{pi}) for each study. Defaults to an independent uniform prior distribution
#' for each intercept \code{prior("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent")}.
#' @param priors_effect list of prior distributions for the effect size (\code{mu})
#' parameter that will be treated as belonging to the alternative hypothesis. Defaults to
#' \code{prior(distribution = "student",   parameters = list(location = 0, scale = 0.58, df = 4))},
#' based on logOR meta-analytic estimates from the Cochrane Database of Systematic Reviews
#' \insertCite{bartos2023empirical}{RoBMA}.
#' @param priors_heterogeneity list of prior distributions for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the alternative hypothesis. Defaults to
#' \code{prior(distribution = "invgamma",  parameters = list(shape = 1.77, scale = 0.55))} that
#' is based on heterogeneities of logOR estimates from the Cochrane Database of Systematic Reviews
#' \insertCite{bartos2023empirical}{RoBMA}.
#' @inheritParams RoBMA
#' @inheritParams combine_data
#'
#' @details The \code{BiBMA()} function estimates the binomial-normal Bayesian model-averaged
#' meta-analysis described in \insertCite{bartos2023empirical;textual}{RoBMA}. See
#' \href{../doc/MedicineBiBMA.html}{\code{vignette("MedicineBiBMA", package = "RoBMA")}}
#' vignette for a reproduction of the \insertCite{oduwole2018honey;textual}{RoBMA} example.
#' Also [RoBMA()] for additional details.
#'
#' Generic [summary.RoBMA()], [print.RoBMA()], and [plot.RoBMA()] functions are
#' provided to facilitate manipulation with the ensemble. A visual check of the
#' individual model diagnostics can be obtained using the [diagnostics()] function.
#' The fitted model can be further updated or modified by [update.RoBMA()] function.
#'
#' @examples \dontrun{
#' # using the example data from Oduwole (2018) and reproducing the example from
#' # Bartos et al. (2023) with domain specific informed prior distributions
#'
#' fit <- BiBMA(
#'   x1          = c(5, 2),
#'   x2          = c(0, 0),
#'   n1          = c(35, 40),
#'   n2          = c(39, 40),
#'   priors_effect        = prior_informed(
#'       "Acute Respiratory Infections",
#'       type = "logOR", parameter = "effect"),
#'   priors_heterogeneity = prior_informed(
#'       "Acute Respiratory Infections",
#'       type = "logOR", parameter = "heterogeneity")
#'  )
#'
#'  summary(fit)
#'
#'  # produce summary on OR scale
#'  summary(fit, output_scale = "OR")
#'
#' }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return \code{NoBMA} returns an object of class 'RoBMA'.
#'
#' @seealso [RoBMA()], [summary.RoBMA()], [update.RoBMA()], [check_setup()]
#' @export
BiBMA <- function(
  # data specification
  x1, x2, n1, n2, study_names = NULL, study_ids = NULL,

  # prior specification
  priors_effect         = prior(distribution = "student",   parameters = list(location = 0, scale = 0.58, df = 4)),
  priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1.77, scale = 0.55)),

  priors_effect_null         = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0)),

  priors_baseline        = NULL,
  priors_baseline_null   = prior_factor("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent"),

  # MCMC fitting settings
  chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
  autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

  # additional settings
  save = "all", seed = NULL, silent = TRUE, ...){

  dots         <- .RoBMA_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  object$data <- .combine_data.bi(x1 = x1, x2 = x2, n1 = n1, n2 = n2, study_names = study_names, study_ids = study_ids, weight = NULL)

  # switch between multivariate and weighted models
  if(attr(object$data, "weighted"))
    .weighted_warning()

  if(.is_multivariate(object))
    stop("Multivariate outcomes are not implemented for binomial outcomes.")


  ### check MCMC settings
  object$fit_control        <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)


  ### prepare and check the settings
  object$priors     <- .check_and_list_priors.bi(
    priors_effect_null = priors_effect_null, priors_effect = priors_effect,
    priors_heterogeneity_null = priors_heterogeneity_null, priors_heterogeneity = priors_heterogeneity,
    priors_baseline_null = priors_baseline_null, priors_baseline = priors_baseline)
  object$models     <- .make_models.bi(object[["priors"]], nrow(object$data), .is_weighted(object))


  ### additional information
  object$add_info <- .check_and_list_add_info(
    model_type             = NULL,
    predictors             = attr(object[["priors"]], "terms"),
    predictors_test        = attr(object[["priors"]], "terms_test"),
    prior_scale            = .transformation_var("logOR"),
    output_scale           = .transformation_var("logOR"),
    effect_measure         = attr(object[["data"]], "effect_measure"),
    effect_direction       = "positive",
    algorithm              = "bridge",
    seed                   = seed,
    save                   = save,
    warnings               = NULL,
    errors                 = NULL
  )

  if(dots[["do_not_fit"]]){
    return(object)
  }


  ### fit the models and compute marginal likelihoods
  if(!object$fit_control[["parallel"]]){

    if(dots[["is_JASP"]]){
      .JASP_progress_bar_start(length(object[["models"]]))
    }

    for(i in seq_along(object[["models"]])){
      object$models[[i]] <- .fit_BiBMA_model(object, i)
      if(dots[["is_JASP"]]){
        .JASP_progress_bar_tick()
      }
    }

  }else{

    fitting_order <- .fitting_priority(object[["models"]])

    cl <- parallel::makePSOCKcluster(floor(RoBMA.get_option("max_cores") / object$fit_control[["chains"]]))
    parallel::clusterEvalQ(cl, {library("RoBMA")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::parLapplyLB(cl, fitting_order, .fit_BiBMA_model, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)

  }

  # create ensemble only if at least one model converged
  if(any(.get_model_convergence(object))){

    # balance probability of non-converged models
    if(object$convergence_checks[["balance_probability"]] && !all(.get_model_convergence(object))){
      object <- .balance_component_probability(object)
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


  class(object) <- c("BiBMA", "RoBMA")
  return(object)
}



#' @title Updates a fitted BiBMA object
#'
#' @description \code{update.BiBMA} can be used to
#' \enumerate{
#'   \item{add an additional model to an existing \code{"BiBMA"} object by
#'    specifying either a null or alternative prior for each parameter
#'    and the prior odds of the model (\code{prior_weights}), see the
#'    \code{vignette("CustomEnsembles")} vignette,}
#'   \item{change the prior odds of fitted models by specifying a vector
#'   \code{prior_weights} of the same length as the fitted models,}
#'   \item{refitting models that failed to converge with updated settings
#'   of control parameters,}
#'   \item{or changing the convergence criteria and recalculating the ensemble
#'   results by specifying new \code{control} argument and setting
#'   \code{refit_failed == FALSE}.}
#' }
#'
#' @param object a fitted BiBMA object
#' @param prior_effect prior distribution for the effect size (\code{mu})
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_heterogeneity prior distribution for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_baseline prior distribution for the intercepts (\code{pi}) of each study
#' that will be treated as belonging to the alternative hypothesis. Defaults to \code{NULL}.
#' @param prior_effect_null prior distribution for the effect size (\code{mu})
#' parameter that will be treated as belonging to the null hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_heterogeneity_null prior distribution for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the null hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_baseline_null prior distribution for the intercepts (\code{pi}) of each study
#' that will be treated as belonging to the null hypothesis. Defaults to \code{NULL}.
#' @inheritParams BiBMA
#' @inheritParams update.RoBMA
#' @param ... additional arguments.
#'
#' @details See [BiBMA()] for more details.
#'
#' @return \code{BiBMA} returns an object of class 'BiBMA'.
#'
#' @seealso [BiBMA()], [summary.RoBMA()], [prior()], [check_setup()]
#' @export
update.BiBMA <- function(object, refit_failed = TRUE, extend_all = FALSE,
                         prior_effect = NULL,      prior_heterogeneity = NULL,      prior_baseline = NULL,      prior_weights = NULL,
                         prior_effect_null = NULL, prior_heterogeneity_null = NULL, prior_baseline_null = NULL,
                         study_names = NULL,
                         chains = NULL, adapt = NULL, burnin = NULL, sample = NULL, thin = NULL, autofit = NULL, parallel = NULL,
                         autofit_control = NULL, convergence_checks = NULL,
                         save = "all", seed = NULL, silent = TRUE, ...){

  BayesTools::check_bool(refit_failed, "refit_failed")
  BayesTools::check_bool(extend_all, "extend_all")

  dots         <- .RoBMA_collect_dots(...)

  if(object$add_info$save == "min")
    stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?BiBMA for more details).")

  # add study names if supplied
  if(!is.null(study_names)){
    if(length(study_names) != nrow(object[["data"]]))
      stop("The study names do not match the length of supplied data.")
    object[["data"]][,"study_names"] <-  as.character(study_names)
  }


  ### choose proper action based on the supplied input
  if((!is.null(prior_effect)         | !is.null(prior_effect_null))  &
     (!is.null(prior_heterogeneity)  | !is.null(prior_heterogeneity_null)) &
     (!is.null(prior_baseline)       | !is.null(prior_baseline_null))){

    what_to_do <- "fit_new_model"
    message("Fitting a new model with specified priors.")
    new_priors <- .check_and_list_priors.bi(
      priors_effect_null        = prior_effect_null,        priors_effect        = prior_effect,
      priors_heterogeneity_null = prior_heterogeneity_null, priors_heterogeneity = prior_heterogeneity,
      priors_baseline_null      = prior_baseline_null,      priors_baseline      = prior_baseline)

    object$models[length(object$models) + 1]  <- list(.make_models.bi(new_priors, nrow(object$data), .is_weighted(object))[[1]])

    if(!is.null(prior_weights)){
      object$models[[length(object$models)]]$prior_weights     <- prior_weights
      object$models[[length(object$models)]]$prior_weights_set <- prior_weights
    }


  }else if(!is.null(prior_weights)){

    what_to_do <- "update_prior_weights"
    message("Updating prior odds for the fitted models.")
    if(length(prior_weights) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.BiBMA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_weights     <- prior_weights[i]
      object$models[[i]]$prior_weights_set <- prior_weights[i]
    }

  }else if(extend_all){

    what_to_do <- "extend_all"
    message("Extending all models with additional samples.")

  }else if(refit_failed & any(!.get_model_convergence(object))){

    what_to_do <- "refit_failed_models"
    message("Refitting models that failed to converge.")

  }else if(!is.null(convergence_checks)){

    # dispatches separately from the rest of the settings
    # recomputing convergence can take a bit of time
    what_to_do <- "update_convergence_checks"
    message("Updating convergence checks.")

  }else{

    what_to_do <- "update_settings"
    message("Updating fitting settings.")

  }


  ### update control settings if any change is specified
  object[["fit_control"]]        <- .update_fit_control(object[["fit_control"]], chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object[["autofit_control"]]    <- .update_autofit_control(object[["autofit_control"]], autofit_control)
  object[["convergence_checks"]] <- .update_convergence_checks(object[["convergence_checks"]], convergence_checks)


  ### clean errors and warnings
  object$add_info[["errors"]]   <- NULL
  object$add_info[["warnings"]] <- .check_effect_direction(object)


  ### do the stuff
  if(what_to_do == "fit_new_model"){

    object[["models"]][[length(object$models)]] <- .fit_BiBMA_model(object, length(object$models))

  }else if(what_to_do %in% c("refit_failed_models", "extend_all")){

    models_to_update <- switch(
      what_to_do,
      "refit_failed_models" = seq_along(object$models)[!.get_model_convergence(object)],
      "extend_all"          = seq_along(object$models)
    )

    if(!object$fit_control[["parallel"]]){

      if(dots[["is_JASP"]]){
        .JASP_progress_bar_start(length(models_to_update))
      }

      for(i in models_to_update){
        object$models[[i]] <- .fit_BiBMA_model(object, i, extend = TRUE)
        if(dots[["is_JASP"]]){
          .JASP_progress_bar_tick()
        }
      }

    }else{

      cl <- parallel::makePSOCKcluster(floor(RoBMA.get_option("max_cores") / object$fit_control[["chains"]]))
      parallel::clusterEvalQ(cl, {library("RoBMA")})
      parallel::clusterExport(cl, "object", envir = environment())
      object$models[models_to_update] <- parallel::parLapplyLB(cl, models_to_update, .fit_BiBMA_model, object = object, extend = TRUE)
      parallel::stopCluster(cl)

    }

  }else if(what_to_do == "update_convergence"){

    # propagate settings changes to the individual models
    for(i in seq_along(object$models)){
      object$models[[i]] <- .update_model_checks(object$models[[i]], object[["convergence_checks"]])
    }

  }

  # restore original prior model probabilities (possibly changed by previous balancing)
  object <- .restore_component_probability(object)

  # create ensemble only if at least one model converged
  if(any(.get_model_convergence(object))){

    # balance probability of non-converged models
    if(object$convergence_checks[["balance_probability"]] && !all(.get_model_convergence(object))){
      object <- .balance_component_probability(object)
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

  return(object)
}
