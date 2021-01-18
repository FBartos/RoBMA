#' @title Estimate a Robust Bayesian Meta-Analysis
#'
#' @description \code{RoBMA} is used to estimate a Robust Bayesian
#' Meta-Analysis. Either t-statistics (\code{t}) and sample sizes of
#' the original studies (\code{n} or \code{n1} and \code{n2}), or
#' effect sizes (\code{d}) and standard errors (\code{se}) can be
#' used to estimate the model.
#'
#' @param t a vector of t-statistics.
#' @param d a vector of effect sizes measured as Cohen's d.
#' @param r a vector of effect sizes measured as correlations.
#' @param y a vector of unspecified effect sizes.
#' @param se a vector of standard errors of the effect sizes.
#' @param OR a vector of odds ratios.
#' @param n a vector of overall sample sizes.
#' @param n1 a vector of sample sizes for the first group.
#' @param n2 a vector of sample sizes for the second group.
#' @param lCI a vector of lower bounds of confidence intervals.
#' @param uCI a vector of upper bounds of confidence intervals.
#' @param test_type a type of test used in the original studies. Options
#' are \code{"two.sample"} (default) and \code{"one.sample"}. Only available
#' if \code{d} is supplied.
#' @param mu_transform transformation to be applied to the supplied
#' effect sizes before fitting the individual models. Only available if
#' correlations or odds ratios are supplied as input. Defaults to
#' \code{"cohens_d"} for correlations (another options is \code{"fishers_z"})
#' and \code{"log_OR"} for odds ratios (another options is \code{"cohens_d"}).
#' Note that priors are specified on the transformed scale and
#' estimates are transformed back (apart from tau).
#' @param effect_direction the expected direction of the effect. The one-sided
#' selection sets the weights omega to 1 to significant results in the expted
#' direction. Defaults to \code{"positive"} (another oprion is \code{"negative"}).
#' @param study_names an optional argument with the names of the studies.
#' @param priors_mu list of prior distributions for the \code{mu} parameter that
#' will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "normal",   parameters = list(mean = 0, sd = 1))}.
#' @param priors_tau list of prior distributions for the \code{tau} parameter that
#' will be treated as belonging to the alternative hypothesis. Defaults to \code{
#' prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))}.
#' @param priors_omega list of prior weight functions for the \code{omega}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{list(
#' prior(distribution = "two.sided", parameters = list(alpha = c(1, 1),     steps = c(.05)),      prior_odds = 1/2),
#' prior(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),  steps = c(.05, .10)), prior_odds = 1/2)
#' )}.
#' @param priors_mu_null list of prior distributions for the \code{mu} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param priors_tau_null list of prior distributions for the \code{tau} parameter that
#' will be treated as belonging to the null hypothesis. Defaults to point distribution
#' with location at 0 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param priors_omega_null list of prior weight functions for the \code{omega} parameter
#' that will be treated as belonging to the null hypothesis. Defaults to point
#' distribution with location at 1 (
#' \code{prior(distribution = "point", parameters = list(location = 0))}).
#' @param chains a number of chains of the MCMC algorithm.
#' @param iter a number of sampling iterations of the MCMC algorithm.
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
#'   \item{cores}{Number of cores to bu used fo parallel computation. If
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
#'               priors_tau = prior("normal",
#'                                  parameters = list(mean = 0, sd = 1),
#'                                  truncation = list(lower = 0, upper = Inf)),
#'               priors_omega = prior("one-sided",
#'                                    parameters = list(cuts = c(.05), alpha = c(1, 1))))
#'
#' # the priors for the null models can be modified or even omited in a similar manner,
#' # allowing to test different (non-nill-null) hypotheses
#' fit2 <- RoBMA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels,
#'               priors_mu_null  = prior("normal",
#'                                  parameters = list(mean = 0, sd = .1),
#'                                  truncation = list(lower = -0.1, upper = 0.1)))
#'
#' # an already fitted RoBMA model can be further updated or modified by using the update function
#' # for example, the prior model probabilities can be changed after the fitting by
#' # (but see '?update.RoBMA' for other posibilities including refitting or adding more models)
#' fit3 <- update(fit2, prior_odds = c(10,1,1,1,1,1,1,1,1,1,1,1))
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
#' @export RoBMA
#' @seealso [summary.RoBMA()], [update.RoBMA()], [prior()], [check_setup()]
RoBMA <- function(t = NULL, d = NULL, r = NULL, y = NULL, OR = NULL, se = NULL, n = NULL, n1 = NULL, n2 = NULL, lCI = NULL, uCI = NULL,
                  test_type = "two.sample", study_names = NULL,
                  mu_transform  = if(!is.null(r)) "cohens_d" else if (!is.null(OR)) "log_OR" else NULL,
                  effect_direction = "positive",
                  priors_mu    = prior(distribution = "normal",   parameters = list(mean = 0, sd = 1)),
                  priors_tau   = prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15)),
                  priors_omega = list(
                    prior(distribution = "two.sided", parameters = list(alpha = c(1, 1),     steps = c(.05)),      prior_odds = 1/2),
                    prior(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),  steps = c(.05, .10)), prior_odds = 1/2)
                  ),
                  priors_mu_null    = prior(distribution = "point", parameters = list(location = 0)),
                  priors_tau_null   = prior(distribution = "point", parameters = list(location = 0)),
                  priors_omega_null = prior(distribution = "point", parameters = list(location = 1)),
                  chains = 3, iter = 10000, burnin = 5000, thin = 1, parallel = FALSE,
                  control = NULL, save = "all", seed = NULL){

  prior_sigma <- NULL
  likelihood  <- "normal"

  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  object$data <- .prepare_data(t, d, r, y, OR, se, n, n1, n2, lCI, uCI, test_type, mu_transform)
  study_names <- .get_study_names(study_names, n_studies = length(object$data$t))


  ### add additional information
  object$add_info <- list(
    t                = object$data$t,
    d                = d,
    r                = r,
    y                = y,
    OR               = OR,
    n                = n,
    n1               = n1,
    n2               = n2,
    se               = se,
    lCI              = lCI,
    uCI              = uCI,
    effect_size      = object$data$effect_size,
    effect_direction = effect_direction,
    mu_transform     = if(object$data$effect_size %in% c("r","OR"))mu_transform,
    test_type        = test_type,
    study_names      = as.character(study_names),
    likelihood       = likelihood,
    seed             = seed,
    save             = save,
    warnings         = NULL
  )


  ### prepare and check the settings
  object$priors   <- .set_priors(priors_mu_null, priors_mu, priors_tau_null, priors_tau, priors_omega_null, priors_omega, prior_sigma, likelihood)
  object$models   <- .get_models(object$priors, likelihood)
  object$control  <- .set_control(control, chains, iter, burnin, thin, likelihood, seed, effect_direction, parallel)
  object$add_info$warnings <- c(object$add_info$warnings, .check_effect_direction(object))


  ### fit the models and compute marginal likelihoods
  if(object$control$cores < 2*object$control$chains){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in 1:length(object$models)){
      object$models[[i]] <- .fit_RoBMA_wrap(object, i)
    }

  }else{

    fitting_order <- .fitting_priority(object$models, object$control$likelihood)

    cl <- parallel::makePSOCKcluster(floor(object$control$cores / object$control$chains))
    parallel::clusterEvalQ(cl, {library("RoBMA")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::clusterApplyLB(cl, fitting_order, .fit_RoBMA_wrap, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)

  }


  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)


  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    # balance probability of non-converged models
    if(object$control$balance_prob & any(!object$add_info$converged))object <- .balance_prob(object, object$add_info$converged)


    ### compute the model-space results
    object$RoBMA         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBMA)
  }


  ### add warnings
  object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  if(save == "min"){
    for(i in 1:length(object$models)){
      if(length(object$models[[1]]$fit) != 1){
        object$models[[i]]$fit$mcmc <- NULL
      }
    }
  }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)warning(w)
  }
  if(sum(!object$add_info$converged) > 0)warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))
  if(!is.null(OR))warning("Analyzing odds ratios is an experimental feature. The performance of default prior distributions was not evaluated.")

  class(object) <- "RoBMA"
  return(object)
}

#' @title Updates a fitted RoBMA object
#'
#' @description \code{update.RoBMA} can be used to
#' \enumerate{
#'   \item{add an additional model to an existing \code{"RoBMA"} object by
#'    specifying either a null or alternative prior for each parameter
#'    and the prior odds of the model (\code{prior_odds}), see the
#'    \code{vignette("CustomEnsembles")} vignette,}
#'   \item{change the prior odds of fitted models by specifying a vector
#'   \code{prior_odds} of the same length as the fitted models,}
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
#' @param prior_odds either a single value specifying prior model odds
#' of a newly specified model using priors argument, or a vector of the
#' same length as already fitted models to update their prior odds.
#' @param refit_failed whether failed models should be refitted. Relevant only
#' if new priors or \code{prior_odds} are not supplied. Defaults to \code{TRUE}.
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
#' fit1 <- update(fit, prior_odds = c(10,1,1,1,1,1,1,1,1,1,1,1))
#'
#' # add an additional model with different priors specification (see '?prior' for more information)
#' fit2 <- update(fit,
#'                priors_mu_null = prior("point", parameters = list(location = 0)),
#'                priors_tau = prior("normal",
#'                                   parameters = list(mean = 0, sd = 1),
#'                                   truncation = list(lower = 0, upper = Inf)),
#'                priors_omega = prior("one-sided",
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
update.RoBMA <- function(object, refit_failed = TRUE,
                         prior_mu = NULL, prior_tau = NULL, prior_omega = NULL, prior_odds = NULL,
                         prior_mu_null = NULL, prior_tau_null = NULL, prior_omega_null = NULL,
                         prior_sigma = NULL,
                         study_names = NULL,
                         control = NULL, chains = NULL, iter = NULL, burnin = NULL, thin = NULL, parallel = NULL, seed = NULL, ...){

  if(object$add_info$save == "min")stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBMA for more details).")

  # add study names if supplied
  if(!is.null(study_names)){
    if(length(study_names) != length(object$data$t))stop("The study names do not match the length of supplied data.")
    object$add_info$study_names <-  as.character(study_names)
  }


  ### choose proper action based on the supplied input
  if((!is.null(prior_mu)    | !is.null(prior_mu_null))  &
     (!is.null(prior_tau)   | !is.null(prior_tau_null)) &
     (!is.null(prior_omega) | !is.null(prior_omega_null))){

    what_to_do <- "fit_new_model"
    object$models[length(object$models) + 1]  <- .get_models(list(
      mu    = .set_parameter_priors(prior_mu_null,    prior_mu,    "mu"),
      tau   = if(object$control$likelihood %in% c("t", "normal")) .set_parameter_priors(prior_tau_null,   prior_tau,   "tau"),
      omega = .set_parameter_priors(prior_omega_null, prior_omega, "omega"),
      sigma = if(object$control$likelihood == "wls") .set_parameter_priors(NULL, prior_sigma, "sigma"),
      likelihood = object$control$likelihood
    ), likelihood = object$control$likelihood)
    if(!is.null(prior_odds))object$models[[length(object$models)]]$prior_odds     <- prior_odds
    if(!is.null(prior_odds))object$models[[length(object$models)]]$prior_odds_set <- prior_odds


  }else if(!is.null(prior_odds)){

    what_to_do <- "update_prior_odds"
    if(length(prior_odds) != length(object$models))stop("The number of newly specified prior odds does not match the number of models. See '?update.RoBMA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_odds     <- prior_odds[i]
      object$models[[i]]$prior_odds_set <- prior_odds[i]
    }

  }else if(refit_failed & any(!object$add_info$converged)){

    what_to_do <- "refit_failed_models"

  }else{

    what_to_do <- "update_settings"

  }


  ### update control settings if any change is specified
  object$control  <- .update_control(object$control, control, chains, iter, burnin, thin, NULL, seed, NULL, parallel)
  object$add_info$warnings <- c(object$add_info$warnings, .check_effect_direction(object))

  ### do the stuff
  if(what_to_do == "fit_new_model"){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    object$models[[length(object$models)]] <- .fit_RoBMA_wrap(object, length(object$models))

  }else if(what_to_do == "refit_failed_models"){

    converged_models <- .get_converged_models(object)
    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in c(1:length(object$models))[!converged_models]){
      object$models[[i]] <- .fit_RoBMA_wrap(object, i)
    }

  }


  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)

  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    # balance probability
    object <- .balance_prob(object, object$add_info$converged)


    ### compute the model-space results
    object$RoBMA         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBMA)
  }


  ### add warnings
  object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  if(object$add_info$save == "min"){
    for(i in 1:length(object$models)){
      if(length(object$models[[1]]$fit) != 1){
        object$models[[i]]$fit$mcmc <- NULL
      }
    }
  }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)warning(w)
  }
  if(sum(!object$add_info$converged) > 0)warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))

  return(object)
}


### data preparation
.prepare_data          <- function(t, d, r, y, OR, se, n, n1, n2, lCI, uCI, test_type, mu_transform){

  data <- list()

  # some general input checks
  if(sum(c(!is.null(t), !is.null(d), !is.null(r), !is.null(y), !is.null(OR))) != 1)stop("One effect size measure needs to be specified.")
  if((is.null(r) & is.null(OR)) & !is.null(mu_transform))stop("'mu_transform' is available only if correlations or odds ratios are supplied as input.")
  if(!is.null(r))if(!mu_transform %in% c("cohens_d","fishers_z"))stop("'mu_transform' must be either 'cohens_d' or 'fishers_z' for correlations.")
  if(!is.null(OR))if(!mu_transform %in% c("cohens_d","log_OR"))stop("'mu_transform' must be either 'log_OR' or 'cohens_d' for odds ratios")
  if(!is.null(test_type))if(!test_type %in% c("one.sample","two.sample"))stop("'test_type' must be either 'one.sample' or 'two.sample'.")
  if((!is.null(lCI) & is.null(uCI)) | is.null(lCI) & !is.null(uCI))stop("Both bounds of the confidence intervals 'lCI' and 'uCI' needs to be supplied together.")


  # check for NA or Inf
  if(!is.null(t))if(any(is.na(t))     | any(is.infinite(t)))stop("NAs or Inf are not allowed in 't'.")
  if(!is.null(d))if(any(is.na(d))     | any(is.infinite(d)))stop("NAs or Inf are not allowed in 'd'.")
  if(!is.null(r))if(any(is.na(r))     | any(is.infinite(r)))stop("NAs or Inf are not allowed in 'r'.")
  if(!is.null(y))if(any(is.na(y))     | any(is.infinite(y)))stop("NAs or Inf are not allowed in 'y'.")
  if(!is.null(OR))if(any(is.na(OR))   | any(is.infinite(OR)))stop("NAs or Inf are not allowed in 'y'.")
  if(!is.null(se))if(any(is.na(se))   | any(is.infinite(se)))stop("NAs or Inf are not allowed in 'se'.")
  if(!is.null(n))if(any(is.na(n))     | any(is.infinite(n)))stop("NAs or Inf are not allowed in 'n'.")
  if(!is.null(n1))if(any(is.na(n1))   | any(is.infinite(n1)))stop("NAs or Inf are not allowed in 'n1'.")
  if(!is.null(n2))if(any(is.na(n2))   | any(is.infinite(n2)))stop("NAs or Inf are not allowed in 'n2'.")
  if(!is.null(lCI))if(any(is.na(lCI)) | any(is.infinite(lCI)))stop("NAs or Inf are not allowed in 'lCI'.")
  if(!is.null(uCI))if(any(is.na(uCI)) | any(is.infinite(uCI)))stop("NAs or Inf are not allowed in 'uCI'.")

  # check formating
  if(!is.null(t))if(!is.numeric(t)     | !is.vector(t))stop("'t' must be a numeric vector.")
  if(!is.null(d))if(!is.numeric(d)     | !is.vector(d))stop("'d' must be a numeric vector.")
  if(!is.null(r))if(!is.numeric(r)     | !is.vector(r))stop("'r' must be a numeric vector.")
  if(!is.null(y))if(!is.numeric(y)     | !is.vector(y))stop("'y' must be a numeric vector.")
  if(!is.null(OR))if(!is.numeric(OR)   | !is.vector(OR))stop("'OR' must be a numeric vector.")
  if(!is.null(se))if(!is.numeric(se)   | !is.vector(se))stop("'se' must be a numeric vector.")
  if(!is.null(n))if(!is.numeric(n)     | !is.vector(n))stop("'n' must be a numeric vector.")
  if(!is.null(n1))if(!is.numeric(n1)   | !is.vector(n1))stop("'n1' must be a numeric vector.")
  if(!is.null(n2))if(!is.numeric(n2)   | !is.vector(n2))stop("'n2' must be a numeric vector.")
  if(!is.null(lCI))if(!is.numeric(uCI) | !is.vector(uCI))stop("'uCI' must be a numeric vector.")
  if(!is.null(uCI))if(!is.numeric(lCI) | !is.vector(lCI))stop("'lCI' must be a numeric vector.")

  # some logical checks
  if(!is.null(r))if(!(all(r > -1) & all(r < 1)))stop("The correlation coefficients 'r' must be in range (-1, 1).")
  if(!is.null(n))if(!all(n > 1))stop("The sample sizes 'n' must be positive.")
  if(!is.null(n1))if(!all(n1 > 0))stop("The sample sizes 'n1' must be positive.")
  if(!is.null(n2))if(!all(n2 > 0))stop("The sample sizes 'n2' must be positive.")
  if(!is.null(se))if(!all(se > 0))stop("The standard errors 'se' must be positive.")
  if(!is.null(lCI) & !is.null(uCI))if(!all(uCI - lCI > 0))stop("The upper confidence interval bounds 'uCI' must be higher than the lower confidence interval bounds 'lCI'.")
  if(!is.null(OR))if(!all(OR > 0))stop("The odds ratios 'OR' must be positive.")
  if(!is.null(OR))if(!all(lCI > 0))stop("The lower confidence interval bounds for odds ratios 'lCI' must be positive.")
  if(!is.null(OR))if(!all(uCI > 0))stop("The upper confidence interval bounds for odds ratios 'uCI' must be positive.")


  if((!is.null(t) | !is.null(d)) & ( !is.null(n) | (!is.null(n1) & !is.null(n2)) | !is.null(se) |  (!is.null(lCI) & !is.null(uCI))) ){

    # obtain test statistics and the ncp multiplactors
    if(!is.null(n)){
      data$df <- n - ifelse(test_type == "two.sample", 2, 1)
      # multiplicator for converting effect sizes into ncp
      if(test_type == "two.sample"){
        if(is.null(t))  t  <- psych::d2t(d = d, n = n)
        if(is.null(d))  d  <- psych::t2d(t = t, n = n)
        if(is.null(se)) se <- (psych::d.ci(d, n = n)[,3] - psych::d.ci(d, n = n)[,1])/(2*stats::qnorm(.975))
        data$ncp_mlp    <- sqrt(n)/2
      }else if(test_type == "one.sample"){
        if(is.null(t))  t  <- psych::d2t(d = d, n1 = n)
        if(is.null(d))  d  <- psych::t2d(t = t, n1 = n)
        if(is.null(se)) se <- (psych::d.ci(d, n1 = n)[,3] - psych::d.ci(d, n1 = n)[,1])/(2*stats::qnorm(.975))
        data$ncp_mlp    <- sqrt(n)
      }
    }else if(!is.null(n1) & !is.null(n2)){
      if(is.null(t))  t  <- psych::d2t(d = d, n1 = n1, n2 = n2)
      if(is.null(d))  d  <- psych::t2d(t = t, n1 = n1, n2 = n2)
      if(is.null(se)) se <- (psych::d.ci(d, n1 = n1, n2 = n2)[,3] - psych::d.ci(d, n1 = n1, n2 = n2)[,1])/(2*stats::qnorm(.975))
      data$ncp_mlp    <- 1/sqrt(1/n1 + 1/n2)
      data$df         <- n1 + n2 - ifelse(test_type == "two.sample", 2, 1)
    }else if(!is.null(se) | (!is.null(lCI) & !is.null(uCI))){
      if(!is.null(lCI) & !is.null(uCI))se <- (uCI - lCI)/(2*stats::qnorm(.975))
      n               <- .get_n_for_d(d, se)
      if(is.null(t))t <- d/se
      data$df         <- n - 2
      data$ncp_mlp    <- sqrt(n)/2
    }

    data$y  <- d
    data$se <- se
    data$t  <- t
    data$K  <- length(t)

    data$effect_size <- "d"

  }else if(!is.null(r) & ((!is.null(n)) | (!is.null(lCI) & !is.null(uCI)))){

    if(!is.null(r) & !is.null(n)){
      if(mu_transform == "cohens_d"){

        # convert to cohen's d and compute it's test statistic
        # using n-2 leads to the actual t-statistic corresponding to the cor.test
        data$y       <- psych::r2d(r)
        data$se      <- (psych::d.ci(data$y, n = n-2)[,3] - psych::d.ci(data$y, n = n-2)[,1])/(2*stats::qnorm(.975))
        data$t       <- psych::d2t(data$y, n = n-2)
        data$ncp_mlp <- sqrt(n-2)/2
        data$df      <- n - 2
        data$K       <- length(n)

      }else if(mu_transform == "fishers_z"){

        # convert to fisher's z and compute it's test statistic
        data$y       <- psych::fisherz(r)
        data$se      <- (1/sqrt(n-3))
        data$t       <- data$y / data$se
        data$ncp_mlp <- sqrt(n-3)
        data$df      <- n - 2
        data$K       <- length(n)

      }
    }else if(!is.null(r) & !is.null(lCI) & !is.null(uCI)){
      if(mu_transform == "cohens_d"){

        # convert to cohen's d and se
        data$y          <- psych::r2d(r)
        data$se         <- (psych::r2d(uCI) - psych::r2d(lCI))/(2*stats::qnorm(.975))
        n               <- .get_n_for_d(data$y, data$se)
        if(is.null(t))t <- data$y/data$se
        data$df         <- n - 2
        data$ncp_mlp    <- sqrt(n)/2
        data$K          <- length(n)

      }else if(mu_transform == "fishers_z"){

        # convert to fisher's z and compute it's test statistic
        d               <- psych::r2d(r)
        d_se            <- (psych::r2d(uCI) - psych::r2d(lCI))/(2*stats::qnorm(.975))
        n               <- .get_n_for_d(d, d_se)
        data$se         <- (1/sqrt(n-3))
        data$y          <- psych::fisherz(r)
        data$t          <- data$y / data$se
        data$ncp_mlp    <- sqrt(n-3)
        data$df         <- n - 2
        data$K          <- length(n)

      }
    }

    data$effect_size <- "r"

  }else if(!is.null(y) & (!is.null(se) | (!is.null(lCI) & !is.null(uCI)))){

    if(!is.null(lCI) & !is.null(uCI))se <- (uCI - lCI)/(2*stats::qnorm(.975))
    # convert to z-statistics (and treat as t-statistics with df = Inf followingly)
    data$y       <- y
    data$K       <- length(y)
    data$t       <- y / se
    data$df      <- rep(999999, length(y)) # should be Inf, but JAGS have problem with that for some reason
    data$ncp_mlp <- 1 / se
    data$se      <- se

    data$effect_size <- "y"

  }else if(!is.null(OR) & !is.null(lCI) & !is.null(uCI)){

    if(mu_transform == "log_OR"){

      data$y           <- log(OR)
      data$se          <- (log(uCI) - log(lCI))/(2*stats::qnorm(.975))
      data$t           <- data$y/data$se
      data$df          <- rep(999999, length(OR))
      data$ncp_mlp     <- 1/data$se
      data$K           <- length(OR)

    }else if(mu_transform == "cohens_d"){
      # transform to Cohen's d
      # https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
      data$y           <- log(OR) * (sqrt(3)/pi)
      data$se          <- (log(uCI) * (sqrt(3)/pi) - log(lCI) * (sqrt(3)/pi))/(2*stats::qnorm(.975))
      n                <- .get_n_for_d(data$y, data$se)
      data$t           <- data$y/data$se
      data$df          <- n - 2
      data$ncp_mlp     <- sqrt(n)/2
      data$K           <- length(OR)
    }


    data$effect_size <- "OR"

  }else{
    stop("Insufficient input provided. Specify either the 't' / 'd' and 'n' / 'n1' & 'n2' / 'se' / 'lCI' & 'uCI', or, 'r' and 'n' / 'lCI' & 'uCI', or 'y' and 'se' / 'lCI' & 'uCU', or 'OR' and 'lCI' & 'uCI'.")
  }

  if(length(data$t) != data$K | length(data$df) != data$K | length(data$ncp_mlp) != data$K)stop("The length of specified data input does not match.")

  return(data)
}


### fitting function
.fit_RoBMA             <- function(object, i){

  model      <- object$models[[i]]
  priors     <- model$priors
  control    <- object$control
  refit_info <- NULL

  # don't sample the complete null model
  if(!.is_model_constant(priors)){

    # generate the model syntax
    model_syntax <- .generate_model_syntax(priors, model$likelihood, control$boost, object$control$effect_direction)


    # remove unnecessary objects from data to mitigate warnings
    fit_data          <- .fit_data(object$data, priors, model$likelihood, control$effect_direction)
    fit_inits         <- .fit_inits(priors, control$chains, control$seed)
    monitor_variables <- .to_monitor(priors, model$likelihood)


    # fit the model
    fit <- .fit_model_RoBMA_wrap(model_syntax, fit_data, fit_inits, monitor_variables, control)


    # deal with some fixable errors
    if(all(class(fit) %in% c("simpleError", "error", "condition"))){

      # problem with installing the RoBMA JAGS module
      if(grepl("Unknown distribution", fit$message)){
        stop("The RoBMA JAGS distributions could not be found. Please, check that the RoBMA package is properly installed.")
      }

      # create a new, data-tuned starting values if there is an outlier that fails the sampling
      if(any(names(unlist(fit_inits)) %in% c("mu", "inv_mu"))){

        fit_inits <- .fit_inits_update(fit_inits, priors, fit_data, model$likelihood)
        refit_info <- "empirical init"

        fit <- .fit_model_RoBMA_wrap(model_syntax, fit_data, fit_inits, monitor_variables, control)

      }

      # try boost library (if it wasn't the primary option)
      if(all(class(fit) %in% c("simpleError", "error", "condition")) & !control$boost){

        refit_info <- "refit with boost"

        model_syntax <- .generate_model_syntax(priors, model$likelihood, TRUE, control$effect_direction)
        fit <- .fit_model_RoBMA_wrap(model_syntax, fit_data, fit_inits, monitor_variables, control)

      }

    }

    # forward error if it's unfixable
    if(all(class(fit) %in% c("simpleError", "error", "condition")) & !control$boost){
      refit_info <- fit$message
    }


  }else{
    fit  <- NULL
  }

  # add the fit and summary to the main object
  model$fit      <- fit
  model$metadata <- list(
    i          = i,
    refit_info = refit_info)
  if(!is.null(fit) & !any(class(fit) %in% c("simpleError", "error", "condition"))){
    model$fit_summary <- .runjags.summary(fit)
  }

  return(model)
}
.fit_model_RoBMA_wrap  <- function(model_syntax, fit_data, fit_inits, monitor_variables, control){
  if(control$silent){
    fit <- callr::r(
      .fit_model_RoBMA,
      args = list(
        model_syntax      = model_syntax,
        fit_data          = fit_data,
        fit_inits         = fit_inits,
        monitor_variables = monitor_variables,
        control           = control
      )
    )
  }else{
    fit <- .fit_model_RoBMA(
      model_syntax      = model_syntax,
      fit_data          = fit_data,
      fit_inits         = fit_inits,
      monitor_variables = monitor_variables,
      control           = control
    )
  }
  return(fit)
}
.fit_model_RoBMA       <- function(model_syntax, fit_data, fit_inits, monitor_variables, control){

  model_call <- list(
    model           = model_syntax,
    data            = fit_data,
    inits           = fit_inits,
    monitor         = monitor_variables,
    n.chains        = control$chains,
    startburnin     = control$burnin,
    startsample     = control$iter,
    adapt           = control$adapt,
    thin            = control$thin,
    raftery.options = if(control$autofit)  list(r = control$max_error) else FALSE,
    psrf.target     = if(control$autofit)  control$max_rhat else Inf,
    max.time        = if(control$autofit)  control$max_time else Inf,
    method          = if(control$parallel) "rjparallel"     else "rjags",
    summarise       = FALSE
  )

  if(control$parallel){
    # the cluster needs to be created manually, because windows don't share the RoBMA JAGS module with the cluster by default
    cl <- parallel::makePSOCKcluster(if(control$chains > control$cores) control$cores else control$chains)
    parallel::clusterCall(cl, function(x) requireNamespace("RoBMA"))
    model_call$cl <- cl
  }else{
    # requires namespace in case that the fit is estimated in a separate R process (for the silent mode)
    requireNamespace("RoBMA")
  }

  if(!is.null(control$seed))set.seed(control$seed)
  fit <- tryCatch(do.call(runjags::autorun.jags, model_call), error = function(e)e)

  if(control$parallel){
    parallel::stopCluster(cl)
  }

  return(fit)
}
.marglik_RoBMA         <- function(object, i){

  model    <- object$models[[i]]
  fit      <- model$fit
  priors   <- model$priors
  control  <- object$control

  # don't sample the complete null model
  if(!.is_model_constant(priors)){


    # deal with failed model
    if(any(class(fit) %in% c("simpleError", "error"))){

      model$marg_lik <- .marglik_fail()

      return(model)
    }


    # compute marginal likelihood
    marglik_samples <- .marglik_prepare_data(fit, priors, model$likelihood, object$data)
    fit_data        <- .fit_data(object$data, priors, model$likelihood, control$effect_direction)

    if(!is.null(control$seed))set.seed(control$seed)
    marg_lik        <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
      samples          = marglik_samples$samples,
      data             = fit_data,
      log_posterior    = .marglik_function,
      likelihood       = model$likelihood,
      priors           = priors,
      effect_direction = control$effect_direction,
      lb               = marglik_samples$lb,
      ub               = marglik_samples$ub,
      maxiter          = control$bridge_max_iter,
      silent           = TRUE)),
      error = function(e)return(e))

  }else{
    # easy calculation of the marginal likelihood in case of null model
    marg_lik <- .marglik_null(object$data, model$likelihood, priors)
  }

  # handle errors
  if(any(class(marg_lik) %in% c("simpleError", "error"))){

    model$metadata$marg_lik <- marg_lik$message
    marg_lik <- .marglik_fail()

  }else if(is.na(marg_lik$logml)){

    model$metadata$marg_lik <- "not enough iterations"
    marg_lik <- .marglik_fail()

  }

  model$marg_lik <- marg_lik

  return(model)
}
.generate_model_syntax <- function(priors, likelihood, boost, effect_direction){

  # generate model syntax
  model_syntax <- "model{\n"

  ### mu priors
  model_syntax <- paste0(model_syntax, .JAGS_distribution("mu", priors$mu$distribution, priors$mu$truncation))

  # transformations
  if(effect_direction == "negative"){
    model_syntax <- paste0(model_syntax, "mu_neg = - mu\n")
  }


  ### tau priors
  if(likelihood %in% c("t", "normal")){
    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location > 0){
        model_syntax <- paste0(model_syntax, .JAGS_distribution("tau", priors$tau$distribution, priors$tau$truncation))
      }
    }else{
      model_syntax <- paste0(model_syntax, .JAGS_distribution("tau", priors$tau$distribution, priors$tau$truncation))
    }
  }


  ### omega priors
  model_syntax <- paste0(model_syntax, .JAGS_distribution_omega(priors$omega$distribution, priors$omega$truncation, priors$omega$parameters))


  ### sigma prior for wls
  if(likelihood == "wls"){
    model_syntax <- paste0(model_syntax, .JAGS_distribution("sigma", priors$sigma$distribution, priors$sigma$truncation))
  }


  ### model
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")

  # marginized random effects the effect size
  if(likelihood %in% c("t", "normal")){
    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location > 0){
        prec <- "1 / ( pow(se[i],2) + pow(tau,2) )"
      }else{
        prec <- "1 / pow(se[i],2)"
      }
    }else{
      prec <- "1 / ( pow(se[i],2) + pow(tau,2) )"
    }
  }else if(likelihood %in% "wls"){
    prec <- "1 / ( pow(se[i],2) * pow(sigma,2) )"
  }

  eff <- ifelse(effect_direction == "negative", "mu_neg", "mu")
  # add PET/PEESE
  if(grepl("PET", priors$omega$distribution)){
    eff <- paste0("(", eff, " + PET * se[i])")
  }else if(grepl("PEESE", priors$omega$distribution)){
    eff <- paste0("(", eff, " + PEESE * pow(se[i], 2))")
  }


  # the observed data
  if(likelihood == "t"){
    # convert to ncp
    ncp <- paste0(eff, "*ncp_mlp[i]")

    if(!priors$omega$distribution %in% c("one.sided", "two.sided")){
      model_syntax <- paste0(model_syntax, "t[i] ~ ",ifelse(boost, "dnt_boost", "dnt"),"(", ncp, ", 1, df[i])\n")
    }else if(priors$omega$distribution == "one.sided"){
      model_syntax <- paste0(model_syntax, "t[i] ~ ",ifelse(boost, "dwt_1s_boost", "dwt_1s"),"(df[i], ", ncp, ", crit_t[i,], omega) \n")
    }else if(priors$omega$distribution == "two.sided"){
      model_syntax <- paste0(model_syntax, "t[i] ~ ",ifelse(boost, "dwt_2s_boost", "dwt_2s"),"(df[i], ", ncp, ", crit_t[i,], omega) \n")
    }
  }else if(likelihood %in% c("normal", "wls")){
    if(!priors$omega$distribution %in% c("one.sided", "two.sided")){
      model_syntax <- paste0(model_syntax, "y[i] ~ ",ifelse(boost, "dnorm_boost", "dnorm"),"(", eff, ",", prec, " )\n")
    }else if(priors$omega$distribution == "one.sided"){
      model_syntax <- paste0(model_syntax, "y[i] ~ ",ifelse(boost, "dwnorm_1s_boost", "dwnorm_1s"),"(", eff, ",", prec, ", crit_x[i,], omega) \n")
    }else if(priors$omega$distribution == "two.sided"){
      model_syntax <- paste0(model_syntax, "y[i] ~ ",ifelse(boost, "dwnorm_2s_boost", "dwnorm_2s"),"(", eff, ",", prec, ", crit_x[i,], omega) \n")
    }
  }


  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.fit_data              <- function(data, priors, likelihood, effect_direction){

  # change the effect size direction (important for one-sided selection)
  if(likelihood == "t"){
    if(effect_direction == "negative"){
      data$t <- - data$t
    }
  }else if(likelihood %in% c("normal", "wls")){
    if(effect_direction == "negative"){
      data$y <- - data$y
    }
  }


  ### add settings for prior distribution
  for(var in names(priors)){

    # don't add parameters for null omega or null tau
    if(var == "tau"){
      if(priors[[var]]$distribution == "point"){
        if(priors[[var]]$parameters$location <= 0)next
      }
    }
    if(likelihood %in% c("t", "normal")){
      if(var == "omega"){
        if(priors[[var]]$distribution == "point")next
      }
    }


    for(par in names(priors[[var]]$parameters)){

      # add the crit_t in case of steps
      if(par == "steps"){

        crit_t <- matrix(ncol = 0, nrow = data$K)

        for(step in priors[[var]]$parameters$steps){
          if(priors[[var]]$distribution == "one.sided"){
            crit_t <- cbind(crit_t, stats::qt(step, data$df, 0, lower.tail = FALSE))
          }else if(priors[[var]]$distribution == "two.sided"){
            crit_t <- cbind(crit_t, stats::qt(step/2, data$df, 0, lower.tail = FALSE))
          }
        }

        data[["crit_t"]]  <- crit_t
        if(priors[[var]]$distribution == "one.sided"){
          if(all(names(priors[[var]]$parameters) %in% c("alpha1", "alpha2", "steps"))){
            data[["J1"]]  <- length(priors[[var]]$parameters$alpha1)
            data[["J2"]]  <- length(priors[[var]]$parameters$alpha2)
          }else if(all(names(priors[[var]]$parameters) %in% c("alpha", "steps"))){
            data[["J"]]   <- length(priors[[var]]$parameters$alpha)
          }
        }else if(priors[[var]]$distribution == "two.sided"){
          data[["J"]]   <- length(priors[[var]]$parameters$alpha)
        }

      }else{
        data[[paste0("prior_",var,"_",par)]] <- priors[[var]]$parameters[[par]]
      }

    }
  }

  # remove unnecessary stuff
  data$effect_size <- NULL
  if(priors$omega$distribution %in% c("one.sided", "two.sided", "point") & likelihood == "t"){
    data$se <- NULL
  }
  if(likelihood == "t"){
    data$y       <- NULL
  }
  if(likelihood %in% c("normal", "wls")){
    data$t       <- NULL
    data$df      <- NULL
    data$ncp_mlp <- NULL
    if(!is.null(data$crit_t)){
      data$crit_x  <- data[["crit_t"]] * data$se
      data$crit_t  <- NULL
    }
  }

  return(data)
}
.fit_inits             <- function(priors, chains, seed){

  if(is.null(seed)){
    seed <- sample(666666, 1)
  }
  inits <- vector(mode = "list", chains)
  set.seed(seed)

  for(i in 1:chains){

    temp_init <- list()

    temp_init <- c(temp_init, .fit_inits_mu_tau(priors$mu, "mu"))
    if(!is.null(priors$tau)){
      temp_init <- c(temp_init, .fit_inits_mu_tau(priors$tau, "tau"))
    }
    if(!is.null(priors$sigma)){
      temp_init <- c(temp_init, .fit_inits_mu_tau(priors$sigma, "sigma"))
    }
    temp_init <- c(temp_init, .fit_inits_omega(priors$omega))

    temp_init[[".RNG.seed"]] <- seed + i
    temp_init[[".RNG.name"]] <- "base::Super-Duper"

    inits[[i]] <- temp_init
  }

  return(inits)
}
.fit_inits_mu_tau      <- function(prior, par){

  temp_x  <- NULL

  # generate the value
  if(prior$distribution == "normal"){
    while(length(temp_x) != 1){
      temp_x <- stats::rnorm(1, mean = prior$parameters$mean, sd = prior$parameters$sd)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "t"){
    while(length(temp_x) != 1){
      temp_x <- extraDistr::rlst(1, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "gamma"){
    while(length(temp_x) != 1){
      temp_x <- stats::rgamma(1, shape = prior$parameters$shape, rate = prior$parameters$rate)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "invgamma"){
    while(length(temp_x) != 1){
      temp_x <- stats::rgamma(1, shape = prior$parameters$shape, rate = prior$parameters$scale)
      temp_x <- temp_x[temp_x >= prior$truncation$upper^-1 & temp_x <= prior$truncation$lower^-1]
    }
  }else if(prior$distribution == "uniform"){
    temp_x <- stats::runif(1, min = prior$parameters$a, max = prior$parameters$b)
  }

  # name the parameter
  if(prior$distribution == "invgamma"){
    names(temp_x) <- paste0("inv_", par)
  }else if(prior$distribution %in% c("normal", "t", "gamma", "uniform")){
    names(temp_x) <- par
  }

  return(temp_x)

}
.fit_inits_omega       <- function(prior){

  temp_x  <- list()
  # the rounding removes some random erros with init values - probably when stardandizing the omega


  if(grepl("PET", prior$distribution)){

    temp_x <- .fit_inits_mu_tau(prior, "PET")

  }else if(grepl("PEESE", prior$distribution)){

    temp_x <- .fit_inits_mu_tau(prior, "PEESE")

  }else if(all(names(prior$parameters) %in% c("alpha", "steps"))){

    temp_x$eta <- round(stats::rgamma(length(prior$parameters$alpha),   prior$parameters$alpha,  1),5)

  }else if(all(names(prior$parameters) %in% c("alpha1", "alpha2", "steps"))){

    temp_x$eta1 <- round(stats::rgamma(length(prior$parameters$alpha1), prior$parameters$alpha1, 1),5)
    temp_x$eta2 <- round(stats::rgamma(length(prior$parameters$alpha2), prior$parameters$alpha2, 1),5)

  }


  return(temp_x)

}
.fit_inits_update      <- function(fit_inits, priors, fit_data, likelihood){

  if(likelihood == "t"){
    new_mu <- mean(psych::t2d(fit_data$t, fit_data$df))
  }else if(likelihood %in% c("normal", "wls")){
    new_mu <- stats::weighted.mean(fit_data$y, 1/fit_data$se^2)
  }

  if(new_mu < priors$mu$truncation$lower){
    new_mu <- priors$mu$truncation$lower + .001
  }else if(new_mu > priors$mu$truncation$upper){
    new_mu <- priors$mu$truncation$upper - .001
  }

  if(any(names(unlist(fit_inits)) == "mu")){
    for(p in 1:length(fit_inits)){
      fit_inits[[p]]$mu <- new_mu
    }
  }else if(any(names(unlist(fit_inits)) == "inv_mu")){
    for(p in 1:length(fit_inits)){
      fit_inits[[p]]$mu <- 1/new_mu
    }
  }

  return(fit_inits)
}
.to_monitor            <- function(priors, likelihood){

  variables <- NULL

  # mu relevant
  if(priors$mu$distribution != "point"){
    if(priors$mu$distribution == "invgamma"){
      variables <- c(variables, "inv_mu")
    }
    variables <- c(variables, "mu")
  }

  if(likelihood %in% c("t", "normal")){
    # tau relevant
    if(priors$tau$distribution != "point"){
      if(priors$tau$distribution == "invgamma"){
        variables <- c(variables, "inv_tau")
      }
      variables <- c(variables, "tau")
    }
  }

  if(likelihood == "wls"){
    if(priors$sigma$distribution != "point"){
      if(priors$sigma$distribution == "invgamma"){
        variables <- c(variables, "inv_sigma")
      }
      variables <- c(variables, "sigma")
    }else if(priors$sigma$parameters$location > 0){
      variables <- c(variables, "sigma")
    }
  }

  # omega relevant
  if(priors$omega$distribution != "point"){
    if(grepl("PET", priors$omega$distribution)){
      variables <- c(variables, "PET")
    }else if(grepl("PEESE", priors$omega$distribution)){
      variables <- c(variables, "PEESE")
    }else if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){
      variables <- c(variables, "eta", "omega")
    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){
      variables <- c(variables, "eta1", "eta2", "omega")
    }
  }

  return(variables)
}
.marglik_function      <- function(samples.row, data, likelihood, priors, effect_direction){

  ### get parameters depending on the model type
  # mu
  if(priors$mu$distribution != "point"){
    if(priors$mu$distribution == "invgamma"){
      inv_mu <- samples.row[[ "inv_mu" ]]
      mu     <- 1/inv_mu
    }else{
      mu <- samples.row[[ "mu" ]]
    }
  }else{
    mu <- priors$mu$parameters$location
  }


  # tau
  if(likelihood %in% c("t", "normal")){
    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location != 0){
        tau      <- priors$tau$parameters$location
      }
    }else{
      if(priors$tau$distribution == "invgamma"){
        inv_tau <- samples.row[[ "inv_tau" ]]
        tau     <- 1/inv_tau
      }else{
        tau <- samples.row[[ "tau" ]]
      }
    }
  }


  # omega
  if(priors$omega$distribution != "point"){
    if(grepl("PET", priors$omega$distribution)){

      PET   <- samples.row[[ "PET" ]]

    }else if(grepl("PEESE", priors$omega$distribution)){

      PEESE <- samples.row[[ "PEESE" ]]

    }else if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){

      eta     <- samples.row[ paste0("eta[",1:data$J,"]") ]
      std_eta <- NULL
      omega   <- NULL
      for(j in 1:data$J){
        std_eta[j]  = eta[j] / sum(eta)
        omega[j]    = sum(std_eta[1:j])
      }

    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){

      eta1     <- samples.row[ paste0("eta1[",1:data$J1,"]") ]
      eta2     <- samples.row[ paste0("eta2[",1:data$J2,"]") ]
      std_eta1 <- NULL
      std_eta2 <- NULL
      omega    <- NULL

      for(j1 in 1:data$J1){
        std_eta1[j1]  <-  eta1[j1] / sum(eta1)
        omega[data$J2 - 1 + j1] = sum(std_eta1[1:j1])
      }
      for(j2 in 1:data$J2){
        std_eta2[j2]  = (eta2[j2] / sum(eta2)) * (1 - std_eta1[1])
      }
      for(j2 in 2:data$J2){
        omega[j2-1] = sum(std_eta2[j2:data$J2]) + std_eta1[1]
      }

    }
  }


  # sigma
  if(likelihood == "wls"){
    if(priors$sigma$distribution == "invgamma"){
      inv_sigma <- samples.row[[ "inv_sigma" ]]
      sigma     <- 1/inv_sigma
    }else{
      sigma <- samples.row[[ "sigma" ]]
    }
  }



  # sd for marginalized random effects
  if(likelihood %in% c("t", "normal")){
    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location != 0){
        pop_sd <- sqrt(data$se^2 + tau^2)
      }else{
        pop_sd <- data$se
      }
    }else{
      pop_sd <- sqrt(data$se^2 + tau^2)
    }
  }else if(likelihood == "wls"){
    pop_sd <- sigma*data$se
  }


  # add PET/PEESE
  eff <- ifelse(effect_direction == "negative", -1, 1)*mu
  if(grepl("PET", priors$omega$distribution)){
    eff <- eff + PET * data$se
  }else if(grepl("PEESE", priors$omega$distribution)){
    eff <- eff + PEESE * data$se^2
  }


  ### compute the marginal log_likelihood
  log_lik <- 0

  # mean
  log_lik <- log_lik + .marglik_distribution(mu, "mu", priors$mu$distribution, data, priors$mu$truncation)

  # tau
  if(likelihood %in% c("t", "normal")){
    log_lik <- log_lik + .marglik_distribution(tau, "tau", priors$tau$distribution, data, priors$tau$truncation)
  }

  # sigma
  if(likelihood == "wls"){
    log_lik <- log_lik + .marglik_distribution(sigma, "sigma", priors$sigma$distribution, data, priors$sigma$truncation)
  }


  # omega
  if(priors$omega$distribution != "point"){
    if(grepl("PET", priors$omega$distribution)){
      log_lik <- log_lik + .marglik_distribution(PET, "omega", substr(priors$omega$distribution, 5, nchar(priors$omega$distribution)), data, priors$omega$truncation)
    }else if(grepl("PEESE", priors$omega$distribution)){
      log_lik <- log_lik + .marglik_distribution(PEESE, "omega", substr(priors$omega$distribution, 7, nchar(priors$omega$distribution)), data, priors$omega$truncation)
    }else if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){
      log_lik <- log_lik + sum(stats::dgamma(eta, data$prior_omega_alpha, 1, log = TRUE))
    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){
      log_lik <- log_lik + sum(stats::dgamma(eta1, data$prior_omega_alpha1, 1, log = TRUE))
      log_lik <- log_lik + sum(stats::dgamma(eta2, data$prior_omega_alpha2, 1, log = TRUE))
    }
  }


  # the individual studies
  if(likelihood == "t"){
    # convert to ncp
    ncp <- eff * data$ncp_mlp

    if(!priors$omega$distribution %in% c("one.sided", "two.sided")){

      temp_t  <- stats::dt(data$t, df = data$df, ncp = ncp, log = TRUE)
      # shift to different t-distribution computation of the classical one returns -Inf
      if(any(is.infinite(temp_t))){
        temp_t[is.infinite(temp_t)] <- DPQ::dntJKBf(data$t[is.infinite(temp_t)], df = data$df[is.infinite(temp_t)], ncp = ncp[is.infinite(temp_t)], log = TRUE)
      }

      log_lik <- log_lik + sum(temp_t)

    }else if(priors$omega$distribution == "one.sided"){
      log_lik <- log_lik + sum(.dwt_fast(data$t, df = data$df, ncp = ncp, omega = omega, crit_t = data$crit_t, type = "one.sided", log = TRUE))
    }else if(priors$omega$distribution == "two.sided"){
      log_lik <- log_lik + sum(.dwt_fast(data$t, df = data$df, ncp = ncp, omega = omega, crit_t = data$crit_t, type = "two.sided", log = TRUE))
    }
  }else if(likelihood %in% c("normal","wls")){
    if(!priors$omega$distribution %in% c("one.sided", "two.sided")){
      log_lik <- log_lik + sum(stats::dnorm(data$y, mean = eff, sd = pop_sd, log = TRUE))
    }else if(priors$omega$distribution == "one.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data$y, mean = eff, sd = pop_sd, omega = omega, crit_x = data$crit_x, type = "one.sided", log = TRUE))
    }else if(priors$omega$distribution == "two.sided"){
      log_lik <- log_lik + sum(.dwnorm_fast(data$y, mean = eff, sd = pop_sd, omega = omega, crit_x = data$crit_x, type = "two.sided", log = TRUE))
    }
  }


  return(log_lik)
}
.marglik_prepare_data  <- function(fit, priors, likelihood, data){

  # posterior samples
  samples <- suppressWarnings(coda::as.mcmc(fit)) # gives a warning about merging chains

  ### get parameteres based on the model type
  pars <- NULL
  lb   <- NULL
  ub   <- NULL


  ### mu
  if(priors$mu$distribution != "point"){

    # add parameters and truncation
    if(priors$mu$distribution == "invgamma"){
      pars <- c(pars, "inv_mu")
      lb   <- c(lb, priors$mu$truncation$upper^-1)
      ub   <- c(ub, priors$mu$truncation$lower^-1)
    }else{
      pars <- c(pars, "mu")
      lb   <- c(lb, priors$mu$truncation$lower)
      ub   <- c(ub, priors$mu$truncation$upper)
    }
  }


  # tau
  if(likelihood %in% c("t", "normal")){
    if(priors$tau$distribution != "point"){
      if(priors$tau$distribution == "invgamma"){
        pars <- c(pars, "inv_tau")
        lb   <- c(lb, rep(-Inf, data$K), priors$tau$truncation$upper^-1)
        ub   <- c(ub, rep( Inf, data$K), priors$tau$truncation$lower^-1)
      }else{
        pars <- c(pars, "tau")
        lb   <- c(lb, rep(-Inf, data$K), priors$tau$truncation$lower)
        ub   <- c(ub, rep( Inf, data$K), priors$tau$truncation$upper)
      }
    }
  }else{
    if(priors$sigma$distribution == "invgamma"){
      pars <- c(pars, "inv_sigma")
      lb   <- c(lb,   priors$sigma$truncation$upper^-1)
      ub   <- c(ub,   priors$sigma$truncation$lower^-1)
    }else{
      pars <- c(pars, "sigma")
      lb   <- c(lb,   priors$sigma$truncation$lower)
      ub   <- c(ub,   priors$sigma$truncation$upper)
    }
  }



  # omega
  if(priors$omega$distribution != "point"){

    if(grepl("PET", priors$omega$distribution)){
      # parameters
      pars <- c(pars, "PET")
      # and truncation
      lb   <- c(lb, priors$omega$truncation$lower)
      ub   <- c(ub, priors$omega$truncation$upper)
    }else if(grepl("PEESE", priors$omega$distribution)){
      # parameters
      pars <- c(pars, "PEESE")
      # and truncation
      lb   <- c(lb, priors$omega$truncation$lower)
      ub   <- c(ub, priors$omega$truncation$upper)
    }else if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){
      # parameters
      pars <- c(pars, paste0("eta[",1:length(priors$omega$parameters$alpha),"]"))
      # and truncation
      lb   <- c(lb, rep(0,   length(priors$omega$parameters$alpha)))
      ub   <- c(ub, rep(Inf, length(priors$omega$parameters$alpha)))
    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){
      # parameters
      pars <- c(pars,
                paste0("eta1[",1:length(priors$omega$parameters$alpha1),"]"),
                paste0("eta2[",1:length(priors$omega$parameters$alpha2),"]"))
      # and truncation
      lb   <- c(lb,
                rep(0,   length(priors$omega$parameters$alpha1)),
                rep(0,   length(priors$omega$parameters$alpha2)))
      ub   <- c(ub,
                rep(Inf, length(priors$omega$parameters$alpha1)),
                rep(Inf, length(priors$omega$parameters$alpha2)))
    }

  }


  # keep only the samples needed for the bridge sampling function
  samples <- samples[,pars]

  # and make cure that they are in a matrix format
  if(!is.matrix(samples)){
    samples           <- matrix(samples, ncol = 1)
    colnames(samples) <- pars
  }

  # name the bounds
  names(ub) <- pars
  names(lb) <- pars

  out <- list(
    samples = samples,
    lb      = lb,
    ub      = ub
  )
  return(out)
}
.marglik_null          <- function(data, likelihood, priors){

  marg_lik        <- NULL
  if(likelihood == "t"){
    marg_lik$logml  <- sum(stats::dt(data$t, data$df, priors$mu$parameters$location, log = TRUE))
  }else if(likelihood =="normal"){
    marg_lik$logml  <- sum(stats::dnorm(data$y, priors$mu$parameters$location, data$se, log = TRUE))
  }else if(likelihood == "wls"){
    marg_lik$logml  <- sum(stats::dnorm(data$y, priors$mu$parameters$location, priors$sigma$parameters$location*data$se, log = TRUE))
  }

  class(marg_lik) <- "bridge"

  return(marg_lik)
}
.marglik_fail          <- function(){
  marg_lik        <- NULL
  marg_lik$logml  <- -Inf
  class(marg_lik) <- "bridge"
  return(marg_lik)
}
.fit_RoBMA_wrap        <- function(object, i){

  object$models[[i]] <- .fit_RoBMA(object, i)
  object$models[[i]] <- .marglik_RoBMA(object, i)
  if(!is.null(object$control$progress_tick))eval(parse(text = object$control$progress_tick))

  return(object$models[[i]])
}

# model building and marginal likelihood tools
.JAGS_distribution       <- function(par_name, distribution, truncation){

  if(par_name %in% c("PET", "PEESE")){
    par_data <- "omega"
  }else{
    par_data <- par_name
  }

  # distribution
  if(distribution == "point"){
    syntax <- paste0(par_name, " = prior_",par_data,"_location\n")
  }else if(distribution == "normal"){
    syntax <- paste0(par_name," ~ dnorm(prior_",par_data,"_mean, pow(prior_",par_data,"_sd, -2))")
  }else if(distribution == "t"){
    syntax <- paste0(par_name," ~ dt(prior_",par_data,"_location, pow(prior_",par_data,"_scale, -2), prior_",par_data,"_df)")
  }else if(distribution == "gamma"){
    syntax <- paste0(par_name," ~ dgamma(prior_",par_data,"_shape, prior_",par_data,"_rate)")
  }else if(distribution == "invgamma"){
    syntax <- paste0("inv_",par_name," ~ dgamma(prior_",par_data,"_shape, prior_",par_data,"_scale)")
  }else if(distribution == "uniform"){
    syntax <- paste0(par_name," ~ dunif(prior_",par_data,"_a, prior_",par_data,"_b)\n")
  }

  # truncation
  if(!distribution %in% c("point", "uniform")){
    if(!(is.infinite(truncation$lower)  & is.infinite(truncation$lower))){
      # the truncation for invgamma needs to be done the other way around since we sample from gamma
      if(distribution == "invgamma"){
        syntax <- paste0(syntax, "T(",
                         ifelse(is.infinite(truncation$upper^-1),"",truncation$upper^-1),
                         ",",
                         ifelse(is.infinite(truncation$lower^-1),"",truncation$lower^-1),
                         ")\n")
      }else{
        syntax <- paste0(syntax, "T(",
                         ifelse(is.infinite(truncation$lower),"",truncation$lower),
                         ",",
                         ifelse(is.infinite(truncation$upper),"",truncation$upper),
                         ")\n")
      }
    }else{
      syntax <- paste0(syntax, "\n")
    }
  }

  # transformations
  if(distribution == "invgamma"){
    syntax <- paste0(syntax, par_name," = pow(inv_",par_name,", -1)\n")
  }

  return(syntax)
}
.JAGS_distribution_omega <- function(distribution, truncation, parameters){

  if(distribution == "point"){
    return()
  }else if(grepl("PET", distribution)){
    syntax <- .JAGS_distribution("PET", substr(distribution, 5, nchar(distribution)), truncation)
  }else if(grepl("PEESE", distribution)){
    syntax <- .JAGS_distribution("PEESE", substr(distribution, 7, nchar(distribution)), truncation)
  }else if(all(names(parameters) %in% c("alpha", "steps"))){
    syntax <-
      "for(j in 1:J){
       eta[j] ~ dgamma(prior_omega_alpha[j], 1)
     }
       for(j in 1:J){
       std_eta[j]  = eta[j] / sum(eta)
       omega[j]    = sum(std_eta[1:j])
     }\n"
  }else if(all(names(parameters) %in% c("alpha1", "alpha2", "steps"))){
    syntax <-
      "for(j1 in 1:J1){
        eta1[j1] ~ dgamma(prior_omega_alpha1[j1], 1)
       }
       for(j2 in 1:J2){
         eta2[j2] ~ dgamma(prior_omega_alpha2[j2], 1)
       }
       for(j1 in 1:J1){
         std_eta1[j1]       = eta1[j1] / sum(eta1)
         omega[J2 - 1 + j1] = sum(std_eta1[1:j1])
       }
         for(j2 in 1:J2){
         std_eta2[j2]  = (eta2[j2] / sum(eta2)) * (1 - std_eta1[1])
       }
       for(j2 in 2:J2){
         omega[j2-1] = sum(std_eta2[j2:J2]) + std_eta1[1]
       }\n"
  }

  return(syntax)
}
.marglik_distribution    <- function(x, par_name, distribution, data, truncation){

  if(distribution == "normal"){
    log_lik <- stats::dnorm(x, mean = data[[paste0("prior_",par_name,"_mean")]], sd = data[[paste0("prior_",par_name,"_sd")]], log = TRUE) -
      log(
        stats::pnorm(truncation$upper, data[[paste0("prior_",par_name,"_mean")]], data[[paste0("prior_",par_name,"_sd")]], lower.tail = TRUE, log.p = FALSE) -
          stats::pnorm(truncation$lower, data[[paste0("prior_",par_name,"_mean")]], data[[paste0("prior_",par_name,"_sd")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "t"){
    log_lik <- extraDistr::dlst(x, df = data[[paste0("prior_",par_name,"_df")]], mu = data[[paste0("prior_",par_name,"_location")]], sigma = data[[paste0("prior_",par_name,"_scale")]], log = TRUE) -
      log(
        extraDistr::plst(truncation$upper, df = data[[paste0("prior_",par_name,"_df")]], mu = data[[paste0("prior_",par_name,"_location")]], sigma = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE) -
          extraDistr::plst(truncation$lower, df = data[[paste0("prior_",par_name,"_df")]], mu = data[[paste0("prior_",par_name,"_location")]], sigma = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "gamma"){
    log_lik <- stats::dgamma(x, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_rate")]], log = TRUE)  -
      log(
        stats::pgamma(truncation$upper, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_rate")]], lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(truncation$lower, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_rate")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "invgamma"){
    log_lik <- stats::dgamma(1/x, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_scale")]], log = TRUE) -
      log(
        stats::pgamma(truncation$lower^-1, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(truncation$upper^-1, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "uniform"){
    log_lik <- stats::dunif(x, min = data[[paste0("prior_",par_name,"_a")]], max = data[[paste0("prior_",par_name,"_b")]], log = TRUE)
  }else if(distribution == "point"){
    log_lik <- 0
  }

  return(log_lik)
}

### model inference functions
.model_inference            <- function(object, n_samples = 10000){

  models    <- object$models
  data      <- object$data
  add_info  <- object$add_info
  converged <- object$add_info$converged
  seed      <- object$control$seed

  # extract marginal likelihoods
  marg_liks <- sapply(models, function(x)x$marg_lik$logml)

  # determine the type of the models
  mm_mu     <- sapply(models, function(m)!.is_parameter_null(m$priors, "mu"))
  if(add_info$likelihood %in% c("t", "normal"))mm_tau <- sapply(models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega  <- sapply(models, function(m)!.is_parameter_null(m$priors, "omega"))

  # extract model weights
  prior_weights_all   <- sapply(models, function(m)m$prior_odds)
  prior_weights_mu    <- ifelse(mm_mu,    prior_weights_all, 0)
  if(add_info$likelihood %in% c("t", "normal"))prior_weights_tau   <- ifelse(mm_tau,   prior_weights_all, 0)
  prior_weights_omega <- ifelse(mm_omega, prior_weights_all, 0)

  # standardize model weights
  prior_weights_all   <- prior_weights_all   / sum(prior_weights_all)
  prior_weights_mu    <- prior_weights_mu    / sum(prior_weights_mu)
  if(add_info$likelihood %in% c("t", "normal"))prior_weights_tau   <- prior_weights_tau   / sum(prior_weights_tau)
  prior_weights_omega <- prior_weights_omega / sum(prior_weights_omega)


  ### compute model weights
  # overall
  weights_all   <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_all)
  if(any(mm_mu) & all(!is.nan(prior_weights_mu))){
    weights_mu  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_mu)
  }else{
    weights_mu <- NULL
  }
  if(add_info$likelihood %in% c("t", "normal")){
    if(any(mm_tau) & all(!is.nan(prior_weights_tau))){
      weights_tau <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_tau)
    }else{
      weights_tau <- NULL
    }
  }else{
    weights_tau <- NULL
  }
  if(any(mm_omega) & all(!is.nan(prior_weights_omega))){
    weights_omega <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_omega)
  }else{
    weights_omega <- NULL
  }


  ### compute inclusion BFs
  BF_effect        <- .inclusion_BF(prior_weights_all, weights_all, mm_mu)
  if(add_info$likelihood %in% c("t", "normal")){
    BF_heterogeneity <- .inclusion_BF(prior_weights_all, weights_all, mm_tau)
  }else if(add_info$likelihood == "wls"){
    BF_heterogeneity <- NULL
  }
  BF_bias          <- .inclusion_BF(prior_weights_all, weights_all, mm_omega)


  ### sample and mix the individual posteriors
  if(!is.null(seed))set.seed(seed)
  samples <- NULL
  samples$averaged    <- .mix_samples(models, weights_all, data, converged, add_info, n_samples)
  samples$conditional <- .mix_samples(models, list(mu = weights_mu, tau = weights_tau, omega = weights_omega), data, converged, add_info, n_samples)


  ### edit names
  names(weights_all)   <- names(models)
  names(weights_mu)    <- names(models)
  names(weights_tau)   <- names(models)
  names(weights_omega) <- names(models)


  # return the results
  output <- list(
    samples         = samples,
    BF = list(
      effect        = BF_effect,
      heterogeneity = if(object$control$likelihood %in% c("t", "normal"))BF_heterogeneity,
      bias          = BF_bias
    ),
    posterior_prob = list(
      all               = weights_all,
      conditional_mu    = weights_mu,
      conditional_tau   = weights_tau,
      conditional_omega = weights_omega
    )
  )
  return(output)
}
.mix_samples                <- function(models, weights, data, converged, add_info, n_samples){

  # metadata about model type
  mm_mu     <- sapply(models, function(m)!.is_parameter_null(m$priors, "mu"))
  if(add_info$likelihood %in% c("t", "normal"))mm_tau <- sapply(models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega  <- sapply(models, function(m)!.is_parameter_null(m$priors, "omega"))

  # the indicies for omega based on type of p-value cutpoints
  omega_ind  <- .get_omega_mapping(models)

  # prepare the object
  samples = list(
    mu    = NULL,
    tau   = NULL,
    omega = matrix(nrow = 0, ncol = if(is.null(omega_ind)) 0 else ncol(do.call(rbind, omega_ind))),
    PET   = NULL,
    PEESE = NULL,
    theta = if(add_info$likelihood %in% c("t", "normal"))matrix(nrow = 0, ncol = data$K),
    sigma = NULL
  )


  for(i in c(1:length(models))[converged]){

    # deal with the complete null model possibility
    if(models[[i]]$priors$mu$distribution == "point" &
       (if(add_info$likelihood == "wls"){FALSE}else if(models[[i]]$priors$tau$distribution == "point"){models[[i]]$priors$tau$parameters$location == 0}else{FALSE}) &
       models[[i]]$priors$omega$distribution == "point"){
      model_samples <- NULL
    }else{
      model_samples <- suppressWarnings(coda::as.mcmc(models[[i]]$fit))
      # and the posibility of only mu parameter
      if(!is.matrix(model_samples)){
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- "mu"
      }
    }

    # make either the conditional of model average weighting
    if(!is.list(weights)){

      # create sampling ind to keep the parameters correlated
      if(!is.null(model_samples))ind <- sample(nrow(model_samples), round(n_samples * weights[i]), replace = TRUE)

      # mu
      samples$mu <- c(samples$mu,
                      if(models[[i]]$priors$mu$distribution == "point"){
                        rep(models[[i]]$priors$mu$parameter$location, round(n_samples * weights[i]))
                      }else{
                        model_samples[ind, "mu"]
                      })

      if(add_info$likelihood %in% c("t", "normal")){
        # tau
        samples$tau <- c(samples$tau,
                         if(models[[i]]$priors$tau$distribution == "point"){
                           rep(models[[i]]$priors$tau$parameter$location, round(n_samples * weights[i]))
                         }else{
                           model_samples[ind, "tau"]
                         })
      }


      # omega
      samples$omega <- rbind(samples$omega,
                             if(models[[i]]$priors$omega$distribution %in% c("one.sided", "two.sided")){
                               model_samples[ind, paste0("omega[",omega_ind[[i]],"]")]
                             }else{
                               matrix(1, ncol = ncol(samples$omega), nrow = round(n_samples * weights[i]))
                             })

      # PET
      samples$PET   <- c(samples$PET,
                         if(grepl("PET", models[[i]]$priors$omega$distribution)){
                           model_samples[ind, "PET"]
                         }else{
                           rep(0, round(n_samples * weights[i]))
                         })

      # PEESE
      samples$PEESE <- c(samples$PEESE,
                         if(grepl( "PEESE", models[[i]]$priors$omega$distribution)){
                           model_samples[ind, "PEESE"]
                         }else{
                           rep(0, round(n_samples * weights[i]))
                         })

      # sigma
      if(add_info$likelihood == "wls"){
        samples$sigma <- c(samples$sigma,
                           if(models[[i]]$priors$sigma$distribution == "point"){
                             rep(models[[i]]$priors$sigma$parameter$location, round(n_samples * weights[i]))
                           }else{
                             model_samples[ind, "sigma"]
                           })
      }

    }else{

      ### notice that the conditional samples doesn't convey information across parameter types

      # mu
      if(!is.null(weights$mu)){
        samples$mu <- c(samples$mu,
                        if(mm_mu[[i]]){
                          if(models[[i]]$priors$mu$distribution == "point"){
                            rep(models[[i]]$priors$mu$parameter$location, round(n_samples * weights$mu[i]))
                          }else{
                            model_samples[sample(nrow(model_samples), round(n_samples * weights$mu[i]), replace = TRUE), "mu"]
                          }
                        })
      }


      if(add_info$likelihood %in% c("t", "normal")){
        # tau
        if(!is.null(weights$tau)){
          samples$tau <- c(samples$tau,
                           if(mm_tau[[i]]){
                             if(models[[i]]$priors$tau$distribution == "point"){
                               rep(models[[i]]$priors$tau$parameter$location, round(n_samples * weights$tau[i]))
                             }else{
                               model_samples[sample(nrow(model_samples), round(n_samples * weights$tau[i]), replace = TRUE), "tau"]
                             }
                           })
        }
      }


      # omega
      if(!is.null(weights$omega[i])){
        samples$omega <- rbind(samples$omega,
                               if(mm_omega[[i]]){
                                 if(models[[i]]$priors$omega$distribution %in% c("one.sided", "two.sided")){
                                   model_samples[sample(nrow(model_samples), round(n_samples * weights$omega[i]), replace = TRUE), paste0("omega[",omega_ind[[i]],"]")]
                                 }else{

                                   matrix(1, ncol = ncol(samples$omega), nrow = round(n_samples * weights$omega[i]))
                                 }
                               })
      }

      # PET
      if(!is.null(weights$omega[i])){
        samples$PET   <- c(samples$PET,
                           if(mm_omega[[i]]){
                             if(grepl("PET", models[[i]]$priors$omega$distribution)){
                               model_samples[sample(nrow(model_samples), round(n_samples * weights$omega[i]), replace = TRUE), "PET"]
                             }else{
                               rep(0, round(n_samples * weights$omega[i]))
                             }
                           })
      }

      # PEESE
      samples$PEESE   <- c(samples$PEESE,
                           if(mm_omega[[i]]){
                             if(grepl("PEESE", models[[i]]$priors$omega$distribution)){
                               model_samples[sample(nrow(model_samples), round(n_samples * weights$omega[i]), replace = TRUE), "PEESE"]
                             }else{
                               rep(0, round(n_samples * weights$omega[i]))
                             }
                           })


    }

  }


  # convert the transformed correlations back if needed
  if(add_info$effect_size == "r"){
    if(add_info$mu_transform == "cohens_d"){
      if(!is.null(samples$mu))    samples$mu    <- psych::d2r(samples$mu)
      if(!is.null(samples$theta)) samples$theta <- psych::d2r(samples$theta)
    }else if(add_info$mu_transform == "fishers_z"){
      if(!is.null(samples$mu))    samples$mu    <- psych::fisherz2r(samples$mu)
      if(!is.null(samples$theta)) samples$theta <- psych::fisherz2r(samples$theta)
    }
  }else if(add_info$effect_size == "OR"){
    if(add_info$mu_transform == "log_OR"){
      if(!is.null(samples$mu))    samples$mu    <- exp(samples$mu)
      if(!is.null(samples$theta)) samples$theta <- exp(samples$theta)
    }else if(add_info$mu_transform == "cohens_d"){
      if(!is.null(samples$mu))    samples$mu    <- .d2OR(samples$mu)
      if(!is.null(samples$theta)) samples$theta <- .d2OR(samples$theta)
    }

  }

  # fix omega names
  all_cuts                <- .get_omega_mapping(models, cuts_only = TRUE)
  if(!is.null(all_cuts))colnames(samples$omega) <- sapply(1:(length(all_cuts)-1), function(i)paste0("omega[",all_cuts[i],",",all_cuts[i+1],"]"))
  # fix theta names
  if(!is.null(samples$theta))colnames(samples$theta) <- paste0("theta[", 1:ncol(samples$theta), "]")

  # remove PET/PEESE if none of the models used it
  if(!any(grepl("PET", sapply(models, function(m)m$priors$omega$distribution)))){
    samples$PET   <- NULL
  }
  if(!any(grepl("PEESE", sapply(models, function(m)m$priors$omega$distribution)))){
    samples$PEESE <- NULL
  }
  if(!any(sapply(models, function(m)m$priors$omega$distribution) %in% c("one.sided", "two.sided"))){
    samples$omega <- NULL
  }
  if(add_info$likelihood != "wls"){
    samples$sigma <- NULL
  }

  return(samples)
}
.compute_coeficients        <- function(RoBMA){
  return(c(
    "mu"     = if(length(RoBMA$samples$averaged$mu) != 0)   mean(RoBMA$samples$averaged$mu),
    "tau"    = if(length(RoBMA$samples$averaged$tau) != 0)  mean(RoBMA$samples$averaged$tau),
    if(length(RoBMA$samples$averaged$omega)   != 0)         apply(RoBMA$samples$averaged$omega, 2, mean),
    "PET"   = if(length(RoBMA$samples$averaged$PET)   != 0) mean(RoBMA$samples$averaged$PET),
    "PEESE" = if(length(RoBMA$samples$averaged$PEESE) != 0) mean(RoBMA$samples$averaged$PEESE),
    "sigma" = if(length(RoBMA$samples$averaged$sigma) != 0) mean(RoBMA$samples$averaged$sigma)
  ))
}
.inclusion_BF               <- function(prior_weights, posterior_weights, conditional_models){
  (sum(posterior_weights[conditional_models])/sum(posterior_weights[!conditional_models]))  /
    (sum(prior_weights[conditional_models])/sum(prior_weights[!conditional_models]))
}
.get_converged_models       <- function(object){

  converged <- NULL

  # basic convergence checks
  for(i in 1:length(object$models)){
    if(!.is_model_constant(object$models[[i]]$priors)){

      if(any(class(object$models[[i]]$fit) %in% c("simpleError", "error")) | is.infinite(object$models[[i]]$marg_lik$logml) | is.na(object$models[[i]]$marg_lik$logml)){
        converged <- c(converged, FALSE)
      }else{
        converged <- c(converged, TRUE)
      }

    }else{
      converged <- c(converged, TRUE)
    }
  }

  object$models <- object$models[converged]

  # remove models with unsatisfactory performance
  if(!is.null(object$control$allow_max_error) |!is.null(object$control$allow_max_rhat) | !is.null(object$control$allow_min_ESS)){
    diagnostics_summary <- summary.RoBMA(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

    # deal with NAs for null models
    diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
    diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
    diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf


    if(!is.null(object$control$allow_max_error)){
      converged <- converged & (diagnostics_summary$"max(MCMC error)" < object$control$allow_max_error)
    }
    if(!is.null(object$control$allow_max_Rhat)){
      converged <- converged & diagnostics_summary$"max(Rhat)" < object$control$allow_max_rhat
    }
    if(!is.null(object$control$allow_min_ESS)){
      converged <- converged & diagnostics_summary$"min(ESS)"  > object$control$allow_min_ESS
    }
  }

  return(converged)
}
.balance_prob               <- function(object, converged_models){

  # extract data
  mm_mu      <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "mu"))
  mm_tau     <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega   <- sapply(object$models, function(m)!.is_parameter_null(m$priors, "omega"))
  prior_odds <- sapply(object$models, function(m)m$prior_odds_set)

  # check whether there is a comparable model for each non-converged models
  for(i in c(1:length(object$models))[!converged_models]){

    temp_ind  <- c(1:length(object$models))[-i]
    temp_same <- temp_ind[mm_mu[-i] == mm_mu[i] & mm_tau[-i] == mm_tau[i] & mm_omega[-i] == mm_omega[i] & converged_models[-i]]

    # if yes, transfer the prior odds
    if(length(temp_same) >= 1){
      prior_odds[temp_same] <- prior_odds[temp_same] + prior_odds[i] / length(temp_same)
      prior_odds[i] <- 0
      object$add_info$warnings <- c(object$add_info$warnings, "Some of the models failed to converge. However, there were other models with the same combination of presence/absence of effect/heterogeneity/publication bias and their prior probability was increased to account for the failed models.")
    }else{
      prior_odds[i] <- 0
      object$add_info$warnings <- c(object$add_info$warnings, "Some of the models failed to converge and their prior probability couldn't be balanced over models with the same combination of presence/absence of effect/heterogeneity/publication bias since they don't exist.")
    }
  }

  for(i in 1:length(object$models)){
    object$models[[i]]$prior_odds <- prior_odds[i]
  }

  return(object)
}
.model_refit_warnings       <- function(metadata){

  new_warn <- NULL

  # extract meta-data with fit-refit information
  refit_info <- t(sapply(metadata, function(x){
    if(is.null(x$refit_info)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$refit_info))
    }
  }))

  marglik_info <- t(sapply(metadata, function(x){
    if(is.null(x$marg_lik)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$marg_lik))
    }
  }))

  if(is.null(dim(refit_info)))  refit_info   <- matrix(refit_info,   ncol = 2)
  if(is.null(dim(marglik_info)))marglik_info <- matrix(marglik_info, ncol = 2)


  if(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data.",
      ifelse(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data and the model was refitted using boost likelihood function.",
      ifelse(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[!refit_info[, 2] %in% c("empirical init", "refit with boost") & !is.na(refit_info[,2]), 1]) > 0){
    refit_info_messages_i <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 1]
    refit_info_messages   <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 2]

    for(i in 1:length(refit_info_messages_i)){
      new_warn <- c(new_warn, paste0("Model ", refit_info_messages_i[i]," failed with the following error: ", refit_info_messages[i]))
    }
  }


  if(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Marginal likelihood computation of %1$s %2$s couldn't be completed within the specified number of iterations.",
      ifelse(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) == 1, "model", "models"),
      paste(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    marglik_info_messages_i <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]
    marglik_info_messages   <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 2]

    for(i in 1:length(marglik_info_messages_i)){
      new_warn <- c(new_warn, paste0("Marginal likelihood computation of model ", marglik_info_messages_i[i]," failed with the following error: ", marglik_info_messages[i]))
    }
  }

  return(new_warn)
}
.model_convergence_warnings <- function(object){

  new_warn <- NULL

  # used set values if specified by the user
  threshold_error <- ifelse(is.null(object$control$allow_max_error), Inf, object$control$allow_max_error)
  threshold_rhat  <- ifelse(is.null(object$control$allow_max_rhat), 1.05, object$control$allow_max_rhat)
  threshold_ESS   <- ifelse(is.null(object$control$allow_max_error), 100, object$control$allow_min_ESS)

  # get the diagnostics summary
  diagnostics_summary <- summary.RoBMA(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

  # deal with NAs for null models
  diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
  diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
  diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf

  # find the problematic models
  warning_error <- rownames(diagnostics_summary)[diagnostics_summary$"max(MCMC error)" > threshold_error]
  warning_rhat  <- rownames(diagnostics_summary)[diagnostics_summary$"max(Rhat)"       > threshold_rhat]
  warning_ESS   <- rownames(diagnostics_summary)[diagnostics_summary$"min(ESS)"        < threshold_ESS]

  # add warnings messages
  if(length(warning_error) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with MCMC error larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_error) == 1, "Model", "Models"),
      paste(warning_error, collapse = ", "),
      threshold_error
    ))
  }

  if(length(warning_rhat) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with R-hat larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_rhat) == 1, "Model", "Models"),
      paste(warning_rhat, collapse = ", "),
      threshold_rhat
    ))
  }

  if(length(warning_ESS) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with ESS lower than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_ESS) == 1, "Model", "Models"),
      paste(warning_ESS, collapse = ", "),
      threshold_ESS
    ))
  }

  return(new_warn)
}

### helper functions for settings
.set_priors             <- function(priors_mu_null, priors_mu, priors_tau_null, priors_tau, priors_omega_null, priors_omega, prior_sigma, likelihood){
  priors <- list()
  priors$mu <- .set_parameter_priors(priors_mu_null, priors_mu, "mu")
  if(likelihood %in% c("t", "normal")){
    priors$tau   <- .set_parameter_priors(priors_tau_null, priors_tau, "tau")
  }else if(likelihood %in% "wls"){
    priors$sigma <- .set_parameter_priors(NULL, prior_sigma, "sigma")
  }
  priors$omega <- .set_parameter_priors(priors_omega_null, priors_omega, "omega")

  return(priors)
}
.set_parameter_priors   <- function(priors_null, priors_alt, parameter){

  # check that at least one prior is specified (either null or alternative)
  if(is.null(priors_null) & is.null(priors_alt))stop(paste0("At least one prior needs to be specified for the ", parameter," parameter (either null or alternative)."))

  # create an empty list if user didn't specified priors
  if(is.null(priors_null)){
    priors_null <- list()
  }else{
    # make sure that priors are passed as a list
    if(class(priors_null) == "RoBMA.prior"){
      priors_null <- list(priors_null)
    }
    # mark the priors as null
    for(p in 1:length(priors_null)){
      priors_null[[p]]$is_null <- TRUE
    }
  }
  if(is.null(priors_alt)){
    priors_alt <- list()
  }else{
    # make sure that priors are passed as a list
    if(class(priors_alt) == "RoBMA.prior"){
      priors_alt <- list(priors_alt)
    }
    # mark the priors as alternative
    for(p in 1:length(priors_alt)){
      priors_alt[[p]]$is_null <- FALSE
    }
  }

  # join null and alternative priors
  priors <- c(priors_null, priors_alt)


  ### check that the speciefied ones are valid
  # check that all objects of the priors list are a 'RoBMA.prior'
  if(!all(sapply(priors, function(prior)class(prior) == "RoBMA.prior")))stop(paste0("Argument priors_", parameter, " does not contain valid prior distributions. The prior distributions need to be passed as a list of objects specified using 'prior()' function. See '?prior' for more information." ))


  if(parameter == "mu"){

    # check that the passed priors are supported for the parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "t", "gamma", "invgamma", "point", "uniform"))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the mu parameter. See '?prior' for further information."))
      }
    }


  }else if(parameter %in% c("tau", "sigma")){

    # check that the passed priors are supported for the tau/sigma parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "t", "gamma", "invgamma", "point", "uniform"))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the tau parameter. See '?prior' for further information."))
        if(priors[[i]]$distribution == "point"){
          if(priors[[i]]$parameters$location < 0){
            stop(paste0("The location of a point prior distribution for ", parameter, " parameter cannot be negative. See '?prior' for further information."))
          }
        }else if(priors[[i]]$distribution == "uniform"){
          if(priors[[i]]$parameters$a < 0 | priors[[i]]$parameters$b < 0 ){
            stop(paste0("The uniform prior distribution for ", parameter, " parameter cannot be defined on negative numbers. See '?prior' for further information."))
          }
        }else if(priors[[i]]$truncation$lower < 0){
          priors[[i]]$truncation$lower <- 0
          warning(paste0("The range of a prior distribution for ", parameter, " parameter cannot be negative. The lower truncation point was set to zero. See '?prior' for further information."))
        }
      }
    }


  }else if(parameter == "omega"){

    # check that the passed priors are supported for the parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!(priors[[i]]$distribution %in% c("two.sided", "one.sided", "point") | grepl("PET", priors[[i]]$distribution) | grepl("PEESE", priors[[i]]$distribution)))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the omega parameter. See '?prior' for further information."))
      }
    }
  }


  return(priors)

}
.get_models             <- function(priors, likelihood){

  # create models according to the set priors
  models <- NULL
  for(mu in priors$mu){
    if(likelihood %in% c("t", "normal")){
      for(tau in priors$tau){
        for(omega in priors$omega){
          models <- c(
            models,
            list(.create_model(mu, tau, omega, NULL, mu$prior_odds * tau$prior_odds * omega$prior_odds, likelihood))
          )
        }
      }
    }else if(likelihood == "wls"){
      for(omega in priors$omega){
        models <- c(
          models,
          list(.create_model(mu, NULL, omega, priors$sigma[[1]], mu$prior_odds * omega$prior_odds, likelihood))
        )
      }
    }

  }

  return(models)
}
.create_model           <- function(prior_mu, prior_tau, prior_omega, prior_sigma, prior_odds, likelihood){

  priors <- list()

  priors$mu    <- prior_mu
  if(!is.null(prior_tau)){
    priors$tau   <- prior_tau
  }
  if(!is.null(prior_sigma)){
    priors$sigma <- prior_sigma
  }
  priors$omega = prior_omega

  model <- list(
    priors         = priors,
    prior_odds     = prior_odds,
    prior_odds_set = prior_odds,
    likelihood     = likelihood
  )
  class(model) <- "RoBMA.model"

  return(model)

}
.set_control            <- function(control, chains, iter, burnin, thin, likelihood, seed, effect_direction, parallel){

  # set the control list
  if(is.null(control)){
    control$autofit         <- FALSE
    control$adapt           <- 1000
    control$bridge_max_iter <- 10000

    control$allow_max_error <- NULL
    control$allow_max_rhat  <- NULL
    control$allow_min_ESS   <- NULL
    control$allow_inc_theta <- FALSE
    control$balance_prob    <- TRUE

    control$silent          <- FALSE
    control$boost           <- FALSE

    if(parallel){
      control$cores         <- parallel::detectCores() - 1
    }else{
      control$cores         <- 1
    }

  }else{
    if(is.null(control[["max_error"]])){
      control$max_error       <- .01
    }
    if(is.null(control[["max_rhat"]])){
      control$max_rhat        <- 1.05
    }
    if(is.null(control[["max_time"]])){
      control$max_time        <- Inf
    }
    if(is.null(control[["autofit"]])){
      control$autofit         <- FALSE
    }
    if(is.null(control[["adapt"]])){
      control$adapt           <- 1000
    }
    if(is.null(control[["bridge_max_iter"]])){
      control$bridge_max_iter <- 10000
    }
    if(is.null(control[["allow_max_error"]])){
      control$allow_max_error <- NULL
    }
    if(is.null(control[["allow_max_rhat"]])){
      control$allow_max_rhat  <- NULL
    }
    if(is.null(control[["allow_min_ESS"]])){
      control$allow_min_ESS   <- NULL
    }
    if(is.null(control[["allow_inc_theta"]])){
      control$allow_inc_theta <- FALSE
    }
    if(is.null(control[["balance_prob"]])){
      control$balance_prob    <- TRUE
    }
    if(is.null(control[["silent"]])){
      control$silent          <- FALSE
    }
    if(is.null(control[["boost"]])){
      control$boost           <- FALSE
    }
    if(is.null(control[["cores"]])){
      if(parallel){
        control$cores         <- parallel::detectCores() - 1
      }else{
        control$cores         <- 1
      }
    }

  }

  if(control[["cores"]] > 1){
    parallel <- TRUE
  }

  # add the main MCMC settings
  control$chains    <- chains
  control$iter      <- iter
  control$burnin    <- burnin
  control$thin      <- thin
  control$seed      <- seed
  control$parallel  <- parallel
  control$effect_direction <- effect_direction
  control$likelihood       <- likelihood

  .check_control(control)
  return(control)
}
.update_control         <- function(control, control_new, chains, iter, burnin, thin, likelihood, seed, effect_direction, parallel){

  if(!is.null(control_new)){
    for(n in names(control_new)){
      control[[n]] <- control_new[[n]]
    }
  }

  if(!is.null(chains))   control[["chains"]]   <- chains
  if(!is.null(iter))     control[["iter"]]     <- iter
  if(!is.null(burnin))   control[["burnin"]]   <- burnin
  if(!is.null(thin))     control[["thin"]]     <- thin
  if(!is.null(parallel)) control[["parallel"]] <- parallel
  if(!is.null(seed))     control[["seed"]]     <- seed
  if(!is.null(effect_direction)) control[["effect_direction"]] <- effect_direction
  if(!is.null(likelihood))       control[["likelihood"]]       <- likelihood

  # stop if there is not enough samples planned for autojags package
  .check_control(control)

  return(control)
}
.check_control          <- function(control){
  # check whether only known controls were supplied
  known_controls <- c("chains", "iter", "burnin" , "adapt", "thin" ,"autofit", "max_error", "max_rhat", "max_time", "bridge_max_iter", "allow_max_error", "allow_max_rhat", "allow_min_ESS", "allow_inc_theta", "balance_prob", "silent", "progress_start", "progress_tick", "boost", "cores", "seed", "parallel", "effect_direction", "likelihood")
  if(any(!names(control) %in% known_controls))stop(paste0("The following control settings were not recognize: ", paste(names(control[!names(control) %in% known_controls]), collapse = ", ")))

  # check whether essential controls were supplied
  if(is.null(control[["chains"]])) stop("Number of chains must be defined.")
  if(is.null(control[["iter"]]))   stop("Number of iterations must be set.")
  if(is.null(control[["burnin"]])) stop("Number of burnin samples must be set.")
  if(is.null(control[["adapt"]]))  stop("Number of adaptation samples must be set.")
  if(is.null(control[["thin"]]))   stop("Thinning of the posterior samples must be set.")
  if(is.null(control[["effect_direction"]])) stop("The effect size direction must be set.")
  if(is.null(control[["likelihood"]]))       stop("The likelihood must be set.")

  if(!is.numeric(control[["chains"]]) | !control[["chains"]] >= 1)  stop("At least one chains must be set.")
  if(!is.numeric(control[["iter"]])   | !control[["iter"]] >= 1)    stop("Number of iterations must be a positive number.")
  if(!is.numeric(control[["burnin"]]) | !control[["burnin"]] >= 1)  stop("Number of burnin samples must be a positive number.")
  if(!is.numeric(control[["adapt"]])  | !control[["adapt"]] >= 1)   stop("Number of adaptation samples must be a positive number.")
  if(!is.numeric(control[["thin"]])   | !control[["thin"]] >= 1)    stop("Thinning of the posterior samples must be a positive number.")
  if(!is.logical(control[["parallel"]]))                            stop("The usage of parallel evaluation must be a logical argument.")
  if(!is.numeric(control[["cores"]])  | !control[["cores"]] >= 1)   stop("Number of cores must be a positive number.")
  if(!is.numeric(control[["seed"]])   & !is.null(control[["seed"]]))stop("Seed must be a numeric argument.")
  if(!control[["effect_direction"]] %in% c("positive", "negative")) stop("The effect size direction must be either positive or negative.")
  if(!control[["likelihood"]]       %in% c("t", "normal", "wls"))   stop("The likelihood must be either 't', 'normal', or 'wls'.")

  # stop if there is not enough samples planned for autojags package
  if(control[["iter"]]/control[["thin"]] < 4000)stop("At least 4000 iterations after thinning is required to compute the Raftery and Lewis's diagnostic.")
  if(!control[["chains"]] >= 2)stop("The number of chains must be at least 2 so that convergence can be assessed.")

  # check convergence criteria
  if(control$autofit)if(control[["max_error"]] >= 1 | control[["max_error"]] <= 0)stop("The target maximum MCMC error must be within 0 and 1.")
  if(control$autofit)if(control[["max_rhat"]] <= 1)stop("The target maximum R-hat must be higher than 1.")
  if(!is.null(control[["allow_max_error"]])) if(control[["allow_max_error"]] >= 1 | control[["allow_max_error"]] <= 0)stop("The maximum allowed MCMC error must be within 0 and 1.")
  if(!is.null(control[["allow_max_rhat"]]))  if(control[["allow_max_rhat"]] <= 1) stop("The maximum allowed R-hat must be higher than 1.")
  if(!is.null(control[["allow_min_ESS"]]))   if(control[["allow_min_ESS"]] <= 0)  stop("The minimum allowed ESS must be higher than 0.")

  if(control[["parallel"]]){
    if(!try(requireNamespace("parallel")))stop("parallel package needs to be installed for parallel processing. Run 'install.packages('parallel')'")
  }
  # now taken care of by the evaluation outside of R
  # runjags::runjags.options(silent.jags = control$silent, silent.runjags = control$silent)
}
.check_effect_direction <- function(object){

  if(!object$control$effect_direction %in% c("positive", "negative"))stop("'effect_direction' must be either 'positive' or 'negative'")
  warnings <- NULL

  # check whether majority of effect sizes are in expected direction. throw warning if not.
  if(any(sapply(object$priors$omega, function(p)p$distribution) == "one.sided")){
    if(stats::median(object$data$t) > 0 & object$control$effect_direction == "negative" |
       stats::median(object$data$t) < 0 & object$control$effect_direction == "positive"){
      warnings <- "The majority of effect sizes is in the oposite direction than expected. The direction of effect sizes is important for the one-sided weight functions. Please, check the 'effect_direction' argument in 'RoBMA' fitting function."
    }
  }

  # the actual effect size direction changes are done prior and after fitting using the '.fit_data' and '.change_direction' functions

  return(warnings)
}
.fitting_priority       <- function(models, likelihood){

  if(likelihood %in% c("t", "normal")){
    # model fitting difficulty using the following heuristic: random effects > weighted likelihood > non-null models
    fitting_difficulty <- sapply(models, function(model){
      ifelse(model$priors$mu$distribution    == "point", 0, 1) +
        ifelse(model$priors$tau$distribution   == "point", 0, 3) +
        ifelse(model$priors$omega$distribution == "point", 0, 5)
    })
  }else if(likelihood == "wls"){
    fitting_difficulty <- sapply(models, function(model){
      ifelse(model$priors$mu$distribution    == "point", 0, 1) +
        ifelse(model$priors$omega$distribution == "point", 0, 5)
    })
  }


  return(order(fitting_difficulty, decreasing = TRUE))
}

# general helper functions
.is_parameter_null <- function(priors, par){
  return(priors[[par]]$is_null)
}
.is_model_constant <- function(priors){

  constant <- NULL
  for(par in c("mu", "tau", "omega", "sigma")){
    if(!is.null(priors[[par]])){
      constant <- c(constant, priors[[par]]$distribution == "point")
    }
  }

  constant <- all(constant)

  return(constant)
}
.get_omega_mapping <- function(models, cuts_only = FALSE){

  # extract cuts and types
  p_cuts <- sapply(models, function(m)rev(m$priors$omega$parameters$steps), simplify = FALSE)
  p_type <- sapply(models, function(m)m$priors$omega$distribution)

  # remove point distributions, PET, and PEESE
  if(all( (p_type == "point" | grepl("PET", p_type) | grepl("PEESE", p_type) ) ))return(NULL)

  # get new cutpoint appropriate cut-points
  p_cuts_new <- p_cuts
  if(any(p_type == "one.sided")){

    # translate two.sided into one.sided
    for(p in 1:length(p_type)){
      if(p_type[p] == "two.sided")p_cuts_new[[p]] <- c(p_cuts[[p]]/2, 1 - rev(p_cuts[[p]])/2)
    }

  }

  # combine the steps
  all_cuts <- c(0, sort(unique(unlist(p_cuts_new))), 1)

  # return the naming for summary function if only asked for labels
  if(cuts_only){
    return(all_cuts)
  }


  # get lower and upper bounds + indicies
  omega_ind <- list()
  p_bound   <- list()
  for(p in 1:length(p_type)){
    if(!is.null(p_cuts_new[[p]])){

      p_bound[[p]] <- list(
        l = c(0, p_cuts_new[[p]]),
        u = c(p_cuts_new[[p]], 1))

      if(any(p_type == "one.sided")){

        if(p_type[p] == "two.sided"){
          omega_ind[[p]] <- rev(c( (length(p_cuts[[p]]) + 1):2, 1:(length(p_cuts[[p]]) + 1) ))
        }else if(p_type[p] == "one.sided"){
          omega_ind[[p]] <- rev(1:(length(p_cuts[[p]]) + 1))
        }

      }else{
        omega_ind[[p]] <- rev(1:(length(p_cuts[[p]]) + 1))
      }
    }
  }

  # create maping to weights
  omega_mapping <- list()
  for(p in 1:length(p_type)){
    if(!(p_type[p] == "point" | grepl("PET", p_type[p]) | grepl("PEESE", p_type[p]) )){
      omega_mapping[[p]] <- sapply(1:(length(all_cuts)-1), function(i)
        omega_ind[[p]][all_cuts[i] >= p_bound[[p]]$l & all_cuts[i+1] <= p_bound[[p]]$u]
      )
    }
  }


  return(omega_mapping)
}
.get_no_support    <- function(models, par){

  no_support  <- NULL

  all_support <- sapply(models, function(m)m$priors[[par]]$truncation, simplify = FALSE)
  all_support <- do.call(rbind.data.frame, all_support)

  if(!is.null(all_support)){

    # start
    if(!is.infinite(min(all_support$lower))){
      no_support <- c(no_support, list(list(lower = -Inf, upper = min(all_support$lower))))
      temp_end   <- min(all_support$lower)
    }else{
      temp_end   <- -Inf
    }

    # the middle
    all_support <- all_support[order(all_support$lower),]
    for(i in 1:nrow(all_support)){

      # prolong the current coverage
      if(all_support$lower[i] <= temp_end & all_support$upper[i] > temp_end){
        temp_end <- all_support$upper[i]
        next
      }

      # detect the gap
      if(all_support$lower[i] > temp_end){
        no_support <- c(no_support, list(list(lower = temp_end, upper = all_support$lower[i])))
        temp_end   <- all_support$lower[i]
      }

    }

    # the upper part
    if(!is.infinite(max(all_support$upper)))no_support <- c(no_support, list(list(lower = max(all_support$upper), upper = Inf)))
  }

  return(no_support)
}
.get_n_for_d       <- function(d, se){

  # according to https://stats.stackexchange.com/questions/144084/variance-of-cohens-d-statistic
  n <- (d^2 + 8) / (2 * se^2)

  if(any(n < 1))stop("The computed sample size of at least one of the original studies was lower than 1 (based on the Cohen's d and standard error). This does not seem to be correct. Please, check your input.")

  return(n)
}
.get_effect_and_ci <- function(add_info, CI){

  # extract the information
  t  <- add_info$t
  d  <- add_info$d
  r  <- add_info$r
  y  <- add_info$y
  OR <- add_info$OR
  n  <- add_info$n
  n1 <- add_info$n1
  n2 <- add_info$n2
  se <- add_info$se
  lCI<- add_info$lCI
  uCI<- add_info$uCI
  effect_size <- add_info$effect_size
  test_type   <- add_info$test_type
  study_names <- add_info$study_names

  # compute the mean and CI
  if(effect_size == "d"){

    if(!is.null(uCI) & !is.null(lCI))se <- (uCI - lCI)/(2*stats::qnorm(.975))
    if(is.null(n) & !is.null(d) & !is.null(se))n <- .get_n_for_d(d, se)

    if(test_type == "one.sample"){
      out <- psych::d.ci(psych::t2d(t = t, n1 = n), n1 = n, alpha = 1 - CI)
    }else if(test_type == "two.sample"){
      if(!is.null(n1) & !is.null(n2)){
        out <- psych::d.ci(psych::t2d(t = t, n1 = n1, n2 = n2), n1 = n1, n2 = n2, alpha = 1 - CI)
      }else{
        out <- psych::d.ci(psych::t2d(t = t, n = n), n = n, alpha = 1 - CI)
      }
    }
  }else if(effect_size == "r"){

    if(!is.null(uCI) & !is.null(lCI)){
      out <- cbind(lCI, r, uCI)
    }else{
      out <- matrix(suppressWarnings(psych::r.con(r = r, n = n, p = CI, twotailed = TRUE)), ncol = 2)
      out <- cbind(out[,1], r, out[,2])
    }

  }else if(effect_size == "y"){

    if(!is.null(uCI) & !is.null(lCI)){
      out <- cbind(lCI, y, uCI)
    }else{
      out <- cbind(y + stats::qnorm((1-CI)/2) * se, y, y - stats::qnorm((1-CI)/2) * se)
    }
  }else if(effect_size == "OR"){

    out <- cbind(lCI, OR, uCI)

  }

  if(is.null(study_names))study_names <- paste0("Study ", 1:nrow(out))
  out <- data.frame(out)
  colnames(out) <- c("lCI", "est", "uCI")
  out$name <- study_names

  return(out)
}
.get_study_names   <- function(study_names, n_studies){

  if(!is.null(study_names))if(length(study_names) != n_studies)stop("The study names do not match the length of supplied data.")
  if(is.null(study_names))study_names <- paste("Study", 1:n_studies)

  return(study_names)

}
.d2OR              <- function(d)exp(d*pi/sqrt(3))
.transform         <- function(x, effect_size, mu_transform){
  if(!is.null(effect_size)){
    if(effect_size == "r"){
      if(mu_transform == "cohens_d"){
        x   <- psych::d2r(x)
      }else if(mu_transform == "fishers_z"){
        x   <- psych::fisherz2r(x)
        x[is.nan(x)] <- ifelse(x[is.nan(x)] > 0, 1, -1)
      }
    }else if(effect_size == "OR"){
      if(mu_transform == "log_OR"){
        x   <- exp(x)
      }else if(mu_transform == "cohens_d"){
        x   <- .d2OR(x)
      }
    }
  }
  return(x)
}
.runjags.summary   <- function(fit){
  # the only reason for this function is that runjags summary function returns HPD instead of quantile intervals
  invisible(utils::capture.output(summary_fit <- summary(fit, silent.jags = T)))
  model_samples <- suppressWarnings(coda::as.mcmc(fit))

  for(i in 1:nrow(summary_fit)){
    summary_fit[i, "Lower95"] <- stats::quantile(model_samples[,i], .025)
    summary_fit[i, "Upper95"] <- stats::quantile(model_samples[,i], .975)
  }

  return(summary_fit)
}
