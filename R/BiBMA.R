#' @title Estimate a Bayesian Model-Averaged Meta-Analysis of Binomial Data
#'
#' @description \code{BiBMA} estimate a Binomial Bayesian
#' model-averaged meta-analysis. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @param x1 a vector with the number of successes in the first group
#' @param x2 a vector with the number of successes in the second group
#' @param n1 a vector with the number of observations in the first group
#' @param n2 a vector with the number of observations in the second group
#' @param priors_baseline prior distributions for independent intercepts (\code{pi})
#' for each study. Defaults to a uniform prior distribution
#' \code{prior("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent")}.
#' @param priors_baseline_null prior distributions for the null hypothesis about
#' independent intercepts (\code{pi}) for each study. Defaults to \code{NULL}.
#' @inheritParams RoBMA
#' @inheritParams combine_data
#'
#' @details See [RoBMA()] for more details.
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
  model_type   = NULL,
  priors_effect         = prior(distribution = "normal",    parameters = list(mean  = 0, sd = 1)),
  priors_heterogeneity  = prior(distribution = "invgamma",  parameters = list(shape = 1, scale = .15)),

  priors_effect_null         = prior(distribution = "point", parameters = list(location = 0)),
  priors_heterogeneity_null  = prior(distribution = "point", parameters = list(location = 0)),

  priors_baseline        = prior_factor("beta", parameters = list(alpha = 1, beta = 1), contrast = "independent"),
  priors_baseline_null   = NULL,

  # MCMC fitting settings
  chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
  autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

  # additional settings
  save = "all", seed = NULL, silent = TRUE, ...){

  dots         <- .RoBMA_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  object$data <- .combine_data.bi(x1 = x1, x2 = x2, n1 = n1, n2 = n2, study_names = study_names, study_ids = study_ids)

  # switch between multivariate and weighted models
  if(.is_multivariate(object)){
    if(dots[["weighted"]]){
      .weighted_warning()
      attr(object$data, "all_independent") <- TRUE
      attr(object$data, "weighted")        <- TRUE
    }else{
      stop("Multivariate outcomes are not implemented for binomial outcomes.")
      .multivariate_warning()
    }
  }

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

  object <<- object
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


  class(object) <- c("BiBMA", "RoBMA")
  return(object)
}

