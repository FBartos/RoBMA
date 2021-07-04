#' @title Estimate a Robust Bayesian Meta-Analysis
#'
#' @description \code{RoBMA} is used to estimate a Robust Bayesian
#' Meta-Analysis. Either t-statistics (\code{t}) and sample sizes of
#' the original studies (\code{n} or \code{n1} and \code{n2}), or
#' effect sizes (\code{d}) and standard errors (\code{se}) can be
#' used to estimate the model.
#'
#' @param effect_direction the expected direction of the effect. The one-sided
#' selection sets the weights omega to 1 to significant results in the expected
#' direction. Defaults to \code{"positive"} (another option is \code{"negative"}).
#' @param prior_scale a scale used to define priors. Defaults to \code{"cohens_d"}.
#' Other options are \code{"fishers_z"}, correlation coefficient \code{"r"},
#' and \code{"logOR"}. The prior scale does not need to match the effect sizes measure -
#' the samples from prior distributions are internally transformed to match the
#' \code{transformation} of the data. The \code{prior_scale} are corresponds to
#' the scale of default output, but can be changed within the summary function.
#' @param priors_effect list of prior distributions for the \code{mu} parameter that
#' will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "normal",   parameters = list(mean = 0, sd = 1))}.
#' @param priors_heterogeneity list of prior distributions for the \code{tau} parameter that
#' will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))}.
#' @param priors_bias list of prior weight functions for the \code{omega}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{list(
#' prior(distribution = "two.sided", parameters = list(alpha = c(1, 1),     steps = c(.05)),      prior_weights = 1/2),
#' prior(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),  steps = c(.05, .10)), prior_weights = 1/2)
#' )}.
#' @param priors_effect_null list of prior distributions for the \code{mu} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param priors_heterogeneity_null list of prior distributions for the \code{tau} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param priors_bias_null list of prior weight functions for the \code{omega} parameter
#' that will be treated as belonging to the null hypothesis. Defaults to point
#' distribution with location at 1 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param chains a number of chains of the MCMC algorithm.
#' @param sample a number of sampling iterations of the MCMC algorithm.
#' Defaults to \code{10000}, with a minimum of \code{4000}.
#' @param burnin a number of burnin iterations of the MCMC algorithm.
#' Defaults to \code{5000}.
#' @param thin a thinning of the chains of the MCMC algorithm. Defaults to
#' \code{1}.
#' @param control a list of additional arguments for the MCMC algorithm.
#' Possible options are:
#' \describe{
#'   \item{autofit}{Whether the models should be refitted until convergence.
#'   Defaults to \code{FALSE}}
#'   \item{max_error}{The target MCMC error for the autofit function. The
#'   argument is passed to \link[coda]{raftery.diag} as 'r'. Defaults to
#'   \code{.01}.}
#'   \item{max_rhat}{The target Rhat error for the autofit function. The
#'   argument is passed to \link[runjags]{add.summary} as 'psrf.target'.
#'   Defaults to \code{1.05}.}
#'   \item{max_time}{A string specifying the maximum fitting time in case
#'   of autofit. Defaults to \code{Inf}. Can be specified as a number and
#'   a unit (Acceptable units include ’seconds’, ’minutes’, ’hours’, ’days’,
#'   ’weeks’, or the first letter(s) of each), i.e. \code{"1hr"}.}
#'   \item{adapt}{A number of iterations used for MCMC adaptation. Defaults
#'   to \code{1000}.}
#'   \item{bridge_max_iter}{Maximum number of iterations for the
#'   \link[bridgesampling]{bridge_sampler} function. Defaults to \code{10000}}
#'   \item{allow_max_error}{Maximum allowed MCMC error for a model to be taken
#'   into consideration. The model will be removed from the ensemble if it fails to
#'   achieve the set MCMC error. Defaults to \code{NULL} - no model will be
#'   removed based on MCMC error.}
#'   \item{allow_max_rhat}{Maximum allowed Rhat for a model to be taken into
#'   consideration. Model will be removed from the ensemble if it fails to
#'   achieve the set Rhat. Defaults to \code{NULL} - no model will be removed
#'   based on Rhat.}
#'   \item{allow_min_ESS}{Minimum allowed ESS for a model to be taken into
#'   consideration. Model will be removed from the ensemble if it fails to
#'   achieve the set ESS. Defaults to \code{NULL} - no model will be removed
#'   based on ESS.}
#'   \item{balance_prob}{Whether the prior probability of the removed model
#'   should be redistributed to other models with the same type if possible
#'    (crossing of effect / heterogeneity / publication bias). Defaults to
#'    \code{TRUE}.}
#'   \item{silent}{Whether all fitting messages should be suppressed. Defaults
#'   to \code{FALSE}. Ideal for getting rid of the "full precision may not have
#'   been achieved in pnt{final}'" warning that cannot be suppressed in any
#'   other way.}
#'   \item{boost}{Whether the likelihood functions implemented using the boost
#'   C++ library should be used as the first option. The higher precision of
#'   boost allows to estimate models in difficult cases. Defaults to \code{FALSE}.
#'   The R distributions are used as default and boost is used only if they fail.
#'   Warning: the estimation using boost takes about twice as long.}
#'   \item{cores}{Maximum number of cores to be used for parallel computation. If
#'   \code{parallel == TRUE}, the default number is equal to number of cores - 1,
#'   and 1 (no parallel processing otherwise).}
#' }
#' @param save whether all models posterior distributions should be kept
#' after obtaining a model-averaged result. Defaults to \code{"all"} which
#' does not remove anything. Set to \code{"min"} to significantly reduce
#' the size of final object, however, some model diagnostics [check()] will
#' not be available.
#' @param seed a seed to be set before model fitting, marginal likelihood
#' computation, and posterior mixing for exact results reproducibility. Defaults
#' to \code{NULL} - no seed is set.
#' @param parallel whether the individual models should be fitted in parallel.
#' Defaults to \code{FALSE}. The \code{cores} argument within the \code{control}
#' list will overwrite the setting if specified to a number higher than 1.
#'
#' @details The default settings with either t-statistics / Cohen's d effect
#' sizes and sample sizes / standard errors correspond to the ensemble proposed by
#' \insertCite{maier2020}{RoBMA}. The \code{vignette("CustomEnsembles")} and
#' \code{vignette("ReproducingBMA")} vignettes describe how to use [RoBMA()] to fit
#' custom meta-analytic ensembles (see [prior()] for more information about prior
#' distributions). To get help with the error and warning messages,
#' see \code{vignette("WarningsAndErrors")}.
#'
#' The RoBMA function first generates models from a combination of the
#' provided priors for each of the model parameters. Then, the individual models
#' are fitted using \link[runjags]{autorun.jags} function. A marginal likelihood
#' is computed using \link[bridgesampling]{bridge_sampler} function. The individual
#' models are then combined into an ensemble using the posterior model probabilities.
#'
#' Generic [summary.RoBMA()], [print.RoBMA()], and [plot.RoBMA()] functions are
#' provided to facilitate manipulation with the ensemble. A visual check of the
#' individual model diagnostics can be obtained using the [diagnostics()] function.
#' The fitted model can be further updated or modified by [update.RoBMA()] function.
#'
#' @return \code{RoBMA} returns an object of \link[base]{class} \code{"RoBMA"}.
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' # in order to speed up the process, we can reduce the default number of chains, iteration,
#' # and disable the autofit functionality (see ?RoBMA for all possible settings)
#' fit_faster <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels,
#' chains = 2, iter = 5000, control = list(autofit = FALSE))
#'
#' # RoBMA function allows to use different prior specifications
#' # for example, change the prior for tau to be half normal and specify one-sided selection only
#' # on significant p-values (see '?.prior' for all options regarding prior distributions)
#' fit1 <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels,
#'               priors_heterogeneity = prior("normal",
#'                                  parameters = list(mean = 0, sd = 1),
#'                                  truncation = list(lower = 0, upper = Inf)),
#'               priors_bias = prior("one-sided",
#'                                    parameters = list(cuts = c(.05), alpha = c(1, 1))))
#'
#' # the priors for the null models can be modified or even omitted in a similar manner,
#' # allowing to test different (non-nill-null) hypotheses
#' fit2 <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels,
#'               priors_effect_null  = prior("normal",
#'                                  parameters = list(mean = 0, sd = .1),
#'                                  truncation = list(lower = -0.1, upper = 0.1)))
#'
#' # an already fitted RoBMA model can be further updated or modified by using the update function
#' # for example, the prior model probabilities can be changed after the fitting by
#' # (but see '?update.RoBMA' for other possibilities including refitting or adding more models)
#' fit3 <- update(fit2, prior_weights = c(10,1,1,1,1,1,1,1,1,1,1,1))
#'
#' # we can get a quick overview of the model coefficients just by printing the model
#' fit
#'
#' # a more detailed overview using the summary function (see '?summary.RoBMA' for all options)
#' summary(fit)
#'
#' # results of the models can be visualized using the plot function (see ?plot.RoBMA for all options)
#' # for example, the model-averaged mean estimate
#' plot(fit, parameter = "mu")
#'
#' # diagnostics for the individual parameters in individual models can be obtained using diagnostics
#' # function (see 'diagnostics' for all options)
#' diagnostics(fit, parameter = "mu", type = "chains")
#' }
#'
#' @references
#' \insertAllCited{}
#' @inheritParams combine_data
#' @export RoBMA
#' @seealso [summary.RoBMA()], [update.RoBMA()], [prior()], [check_setup()]
RoBMA <- function(
  # data specification
  d = NULL, r = NULL, logOR = NULL, z = NULL, y = NULL,
  se = NULL, v = NULL, n = NULL, lCI = NULL, uCI = NULL, t = NULL, study_names = NULL,
  data = NULL,
  transformation   = if(is.null(y)) "fishers_z" else "none",
  prior_scale      = if(is.null(y)) "cohens_d"  else "none",
  effect_direction = "positive",

  # prior specification
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

  # MCMC fitting settings
  chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
  autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

  # additional settings
  save = "all", seed = NULL, silent = TRUE){

  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  if("data.RoBMA" %in% class(data)){
    object$data <- data
  }else{
    object$data <- combine_data(d = d, r = r, z = z, logOR = logOR, t = t, y = y, se = se, v = v, n = n, lCI = lCI, uCI = uCI, study_names = study_names, data = data, transformation = transformation)
  }


  ### check MCMC settings
  object$fit_control        <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)


  ### additional information
  object$add_info <- .check_and_list_add_info(
    model_type       = model_type,
    prior_scale      = .transformation_var(prior_scale),
    output_scale     = .transformation_var(prior_scale),
    effect_measure   = attr(object$data, "effect_measure"),
    effect_direction = effect_direction,
    seed             = seed,
    save             = save,
    warnings         = NULL,
    errors           = NULL
  )


  ### prepare and check the settings
  object$priors   <- .check_and_list_priors(object$add_info[["model_type"]], priors_effect_null, priors_effect, priors_heterogeneity_null, priors_heterogeneity, priors_bias_null, priors_bias, object$add_info[["prior_scale"]])
  object$models   <- .make_models(object[["priors"]])
  object$add_info$warnings <- c(object$add_info[["warnings"]], .check_effect_direction(object))


  ### fit the models and compute marginal likelihoods
  if(!object$fit_control[["parallel"]]){

    for(i in seq_along(object[["models"]])){
      object$models[[i]] <- .fit_RoBMA_model(object, i)
    }

  }else{

    fitting_order <- .fitting_priority(object[["models"]])

    cl <- parallel::makePSOCKcluster(floor(parallel::detectCores() - 1 / object$fit_control[["chains"]]))
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
  }


  class(object) <- "RoBMA"
  return(object)
}




#' @title Updates a fitted RoBMA object
#'
#' @description \code{update.RoBMA} can be used to
#' \enumerate{
#'   \item{add an additional model to an existing \code{"RoBMA"} object by
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
#' @param object a fitted RoBMA object.
#' @param prior_mu a prior distribution for the \code{mu} parameter that
#' will be treated as belonging to the alternative hypothesis.
#' @param prior_tau a prior distribution for the \code{tau} parameter that
#' will be treated as belonging to the alternative hypothesis.
#' @param prior_omega a prior weight function for the \code{omega}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' @param prior_mu_null list of prior distribution for the \code{mu} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param prior_tau_null a prior distribution for the \code{tau} parameter that
#' will be treated as belonging to the null hypothesis.
#' @param prior_omega_null a prior weight function for the \code{omega} parameter
#' that will be treated as belonging to the null hypothesis.
#' @param prior_weights either a single value specifying prior model odds
#' of a newly specified model using priors argument, or a vector of the
#' same length as already fitted models to update their prior odds.
#' @param refit_failed whether failed models should be refitted. Relevant only
#' if new priors or \code{prior_weights} are not supplied. Defaults to \code{TRUE}.
#' @inheritParams RoBMA
#' @param ... additional arguments.
#'
#' @details See [RoBMA()] for more details.
#'
#' @return \code{RoBMA} returns an object of \link[base]{class} \code{"RoBMA"}.
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' # the update function allows us to change the prior model probability of each model
#' fit1 <- update(fit, prior_weights = c(10,1,1,1,1,1,1,1,1,1,1,1))
#'
#' # add an additional model with different priors specification (see '?prior' for more information)
#' fit2 <- update(fit,
#'                priors_effect_null = prior("point", parameters = list(location = 0)),
#'                priors_heterogeneity = prior("normal",
#'                                   parameters = list(mean = 0, sd = 1),
#'                                   truncation = list(lower = 0, upper = Inf)),
#'                priors_bias = prior("one-sided",
#'                                     parameters = list(cuts = c(.05), alpha = c(1, 1))))
#'
#' # change the model convergence criteria to mark models with ESS lower than 2000 as non-covergent
#' fit3 <- update(fit, control = list(allow_min_ESS = 2000))
#'
#' # and refit them failed models with increased number of burnin iterations
#' fit4 <- update(fit3, burnin = 10000)
#'
#' }
#'
#' @method update RoBMA
#' @export update.RoBMA
#' @rawNamespace S3method(update, RoBMA)
#'
#' @seealso [RoBMA()], [summary.RoBMA()], [prior()], [check_setup()]
update.RoBMA <- function(object, refit_failed = TRUE, output_scale = NULL,
                         prior_effect = NULL,      prior_heterogeneity = NULL,      prior_bias = NULL, prior_weights = NULL,
                         prior_effect_null = NULL, prior_heterogeneity_null = NULL, prior_bias_null = NULL,
                         study_names = NULL,
                         chains = NULL, adapt = NULL, burnin = NULL, sample = NULL, thin = NULL, autofit = NULL, parallel = NULL, cores = NULL, autofit_control = NULL, convergence_checks = NULL,
                         save = "all", seed = NULL, silent = TRUE, ...){
  prior_sigma <- NULL

  if(object$add_info$save == "min")
    stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBMA for more details).")

  # add study names if supplied
  if(!is.null(study_names)){
    if(length(study_names) != nrow(object[["data"]]))
      stop("The study names do not match the length of supplied data.")
    object[["data"]][,"study_names"] <-  as.character(study_names)
  }


  ### choose proper action based on the supplied input
  if((!is.null(prior_effect)         | !is.null(prior_effect_null))  &
     (!is.null(prior_heterogeneity)  | !is.null(prior_heterogeneity_null)) &
     (!is.null(prior_bias)           | !is.null(prior_bias_null))){

    what_to_do <- "fit_new_model"
    new_priors <- .check_and_list_priors(NULL, prior_effect_null, prior_effect, prior_heterogeneity_null, prior_heterogeneity, prior_bias_null, prior_bias, object$add_info[["prior_scale"]])

    object$models[length(object$models) + 1]  <- list(.make_models(new_priors)[[1]])

    if(!is.null(prior_weights)){
      object$models[[length(object$models)]]$prior_weights     <- prior_weights
      object$models[[length(object$models)]]$prior_weights_set <- prior_weights
    }


  }else if(!is.null(prior_weights)){

    what_to_do <- "update_prior_weights"
    if(length(prior_weights) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.RoBMA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_weights     <- prior_weights[i]
      object$models[[i]]$prior_weights_set <- prior_weights[i]
    }

  }else if(!is.null(output_scale)){

    stop("this functionality is not currently implemented")
    what_to_do <- "transform_estimates"

  }else if(refit_failed & any(!.get_model_convergence(object))){

    what_to_do <- "refit_failed_models"

  }else{

    what_to_do <- "update_settings"

  }


  ### update control settings if any change is specified
  object[["fit_control"]]        <- .update_fit_control(object[["fit_control"]], chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object[["autofit_control"]]    <- .update_autofit_control(object[["autofit_control"]], autofit_control)
  object[["convergence_checks"]] <- .update_convergence_checks(object[["convergence_checks"]], convergence_checks)


  ### do the stuff
  if(what_to_do == "fit_new_model"){

    object[["models"]][[length(object$models)]] <- .fit_RoBMA_model(object, length(object$models))

  }else if(what_to_do == "refit_failed_models"){

    for(i in c(1:length(object$models))[!.get_model_convergence(object)]){
      object[["models"]][[i]] <- .fit_RoBMA_model(object, i)
    }

  }else if(what_to_do == "transform_estimates"){

    # TODO: implement
    stop("Not implemented.")
    for(i in c(1:length(object$models))){
      object$models[[i]] <- .transform_posterior(object$models[[i]], object$add_info$output_scale, .transformation_var(output_scale))
    }
    object <- .transform_posterior(object, object$add_info$output_scale, .transformation_var(output_scale))
    object$add_info$output_scale <- .transformation_var(output_scale)

    return(object)
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
  }

  return(object)
}
