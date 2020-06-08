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
#' @param n a vector of overall sample sizes.
#' @param n1 a vector of sample sizes for first group.
#' @param n2 a vector of sample sizes for second group.
#' @param test_type a type of test used in the original studies. Options
#' are \code{"two.sample"} (default) and \code{"one.sample"}. Only available
#' if \code{d} is supplied.
#' @param mu_transform transformation to be applied to the supplied
#' effect sizes before fitting the individual models. Defaults to
#' \code{"cohens_d"} for correlations (another options \code{"fishers_z"}).
#' Note that priors are specified on the transformed scale and
#' estimates are transformed back (apart from tau).
#' @param study_names an optional argument with names of the studies.
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
#' Defaults to \code{10000}, with minimum of \code{4000}.
#' @param burnin a number of burnin iterations of the MCMC algorithm.
#' Defaults to \code{5000}.
#' @param thin a thinning of the chains of the MCMC algorithm. Defaults to
#' \code{1}.
#' @param control a list of additional arguments for the MCMC algorithm.
#' Possible options are:
#' \describe{
#'   \item{autofit}{Whether the models should be refitted until convergence.
#'   Defaults to \code{TRUE}}
#'   \item{max_error}{The target MCMC error for the autofit function. The
#'   argument is passed to \link[coda]{raftery.diag} as 'r'. Defaults to
#'   \code{.01}.}
#'   \item{max_time}{A string specifying the maximum fitting time in case
#'   of autofit. Defaults to \code{Inf}. Can be specified as a number and
#'   a unit (Acceptable units include ’seconds’, ’minutes’, ’hours’, ’days’,
#'   ’weeks’, or the first letter(s) of each), i.e. \code{"1hr"}.}
#'   \item{adapt}{A number of iterations used for MCMC adaptation. Defaults
#'   to \code{1000}.}
#'   \item{bridge_max_iter}{Maximum number of iterations for the
#'   \link[bridgesampling]{bridge_sampler} function. Defaults to \code{10000}}
#'   \item{allow_max_error}{Maximum allowed MCMC error for a model to be taken
#'   into consideration. Model will be removed from the ensemble if it fails to
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
#'   \item{allow_inc_theta}{Whether the diagnostics for theta should be
#'   included into model removal decision. Defaults to \code{NULL} - only
#'   'mu', 'tau', and 'omega' estimates will be taken into account.}
#'   \item{balance_prob}{Whether the prior probability of removed model
#'   should be redistributed to other models with the same type if possible
#'    (crossing of effect / heterogeneity / publication bias). Defaults to
#'    \code{TRUE}.}
#'   \item{silent}{Whether all fitting messages should be suppressed. Defaults
#'   to \code{FALSE}.}
#' }
#' @param save whether all models posterior distributions should be kept
#' after obtaining a model-averaged result. Defaults to \code{"all"} which
#' does not remove anything. Set to \code{"min"} to significantly reduce
#' the size of final object, however, some model diagnostics [check()] will
#' not be available.
#' @param seed a seed to be set before model fitting, marginal likelihood
#' computation, and posterior mixing for exact results reproducibility. Defaults
#' to \code{NULL} - no seed is set.
#'
#' @details The RoBMA function first generates models from combination of the
#' provided priors for each of the model parameter. Then, the individual models
#' are fitted using \link[runjags]{autorun.jags} function. A marginal likelihood
#' is computed using \link[bridgesampling]{bridge_sampler} function. The individual
#' models are then combined into an ensemble using posterior model probabilities.
#'
#' Generic [summary.RoBMA()], [print.RoBMA()], and [plot.RoBMA()] functions are
#' provided to facilitate manipulation with the ensemble. A visual check of the
#' individual model diagnostics can be obtained using the [diagnostics()] function.
#' The fitted model can be further updated or modified by [update.RoBMA()] function.
#'
#' See [prior()] for more information about possible prior specifications options.
#'
#' @return \code{RoBMA} returns an object of \link[base]{class} \code{"RoBMA"}.
#'
#' @examples \donttest{
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
#' @export RoBMA
#' @seealso [summary.RoBMA()], [update.RoBMA()], [prior()], [check_setup()]
RoBMA <- function(t = NULL, d = NULL, r = NULL, y = NULL, se = NULL, n = NULL, n1 = NULL, n2 = NULL,
                  test_type = "two.sample", study_names = NULL,
                  mu_transform  = if(!is.null(r)) "cohens_d" else NULL,
                  priors_mu    = prior(distribution = "normal",   parameters = list(mean = 0, sd = 1)),
                  priors_tau   = prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15)),
                  priors_omega = list(
                    prior(distribution = "two.sided", parameters = list(alpha = c(1, 1),     steps = c(.05)),      prior_odds = 1/2),
                    prior(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),  steps = c(.05, .10)), prior_odds = 1/2)
                  ),
                  priors_mu_null    = prior(distribution = "point", parameters = list(location = 0)),
                  priors_tau_null   = prior(distribution = "point", parameters = list(location = 0)),
                  priors_omega_null = prior(distribution = "point", parameters = list(location = 1)),
                  chains = 3, iter = 10000, burnin = 5000, thin = 1,
                  control = NULL, save = "all", seed = NULL){


  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  object$data <- .prepare_data(t, d, r, y, se, n, n1, n2, test_type, mu_transform)
  study_names <- .get_study_names(study_names, n_studies = length(object$data$t))


  ### prepare and check the settings
  object$priors  <- list(
    mu    = .set_parameter_priors(priors_mu_null,    priors_mu,    "mu"),
    tau   = .set_parameter_priors(priors_tau_null,   priors_tau,   "tau"),
    omega = .set_parameter_priors(priors_omega_null, priors_omega, "omega")
  )
  object$models  <- .get_models(object$priors)
  object$control <- .set_control(control, chains, iter, burnin, thin)


  ### add additional information
  object$add_info <- list(
    t            = t,
    d            = d,
    r            = r,
    y            = y,
    n            = n,
    n1           = n1,
    n2           = n2,
    se           = se,
    effect_size  = object$data$effect_size,
    mu_transform = if(object$data$effect_size == "r")mu_transform,
    test_type    = test_type,
    study_names  = as.character(study_names),
    seed         = seed,
    save         = save
  )


  ### fit the models and compute marginal likelihoods
  if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
  for(i in 1:length(object$models)){
    object <- .fit_RoBMA(object, i)
    object <- .marglik_RoBMA(object, i)
    if(!is.null(object$control$progress_tick))eval(parse(text = object$control$progress_tick))
  }


  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)


  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    # balance probability of non-converged models
    if(object$control$balance_prob & any(!object$add_info$converged))object <- .balance_prob(object, object$add_info$converged)


    ### compute the model-space results
    object$RoBMA         <- .model_inference(object$models, object$data, object$add_info$converged, ifelse(!is.null(object$add_info$r), object$add_info$mu_transform, FALSE), object$add_info$seed)
    object$coefficients  <- .compute_coeficients(object$RoBMA)
  }


  ### remove model posteriors if asked to
  if(save == "min"){
    for(i in 1:length(object$models)){
      if(length(object$models[[1]]$fit) != 1){
        object$models[[i]]$fit$mcmc <- NULL
      }
    }
  }


  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)warning(w)
  }
  if(sum(!object$add_info$converged) > 0)warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))

  class(object) <- "RoBMA"
  return(object)
}

#' @title Updates a fitted RoBMA object
#'
#' @description \code{update.RoBMA} can be used to
#' \enumerate{
#'   \item{add an additional model to an existing \code{"RoBMA"} object by
#'    specifying either a null or alternative prior for each parameter
#'    and the prior odds of the model (\code{prior_odds}),}
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
#' @examples \donttest{
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
                         study_names = NULL,
                         control = NULL, chains = NULL, iter = NULL, burnin = NULL, thin = NULL, ...){

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
      tau   = .set_parameter_priors(prior_tau_null,   prior_tau,   "tau"),
      omega = .set_parameter_priors(prior_omega_null, prior_omega, "omega")
    ))
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
  object$control <- .update_control(object$control, control, chains, iter, burnin, thin)


  ### do the stuff
  if(what_to_do == "fit_new_model"){

    object <- .fit_RoBMA(object, length(object$models))
    object <- .marglik_RoBMA(object, length(object$models))

  }else if(what_to_do == "refit_failed_models"){

    converged_models <- .get_converged_models(object)
    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in c(1:length(object$models))[!converged_models]){
      object <- .fit_RoBMA(object, i)
      object <- .marglik_RoBMA(object, i)
      if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_tick))
    }

  }


  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)

  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    # balance probability
    object <- .balance_prob(object, object$add_info$converged)


    ### compute the model-space results
    object$RoBMA         <- .model_inference(object$models, object$data, object$add_info$converged, ifelse(!is.null(object$add_info$r), object$add_info$mu_transform, FALSE), object$add_info$seed)
    object$coefficients  <- .compute_coeficients(object$RoBMA)
  }


  ### remove model posteriors if asked to
  if(object$add_info$save == "min"){
    for(i in 1:length(object$models)){
      if(length(object$models[[1]]$fit) != 1){
        object$models[[i]]$fit$mcmc <- NULL
      }
    }
  }


  if(sum(!object$add_info$converged) > 0)warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))

  return(object)
}


### data preparation
.prepare_data          <- function(t, d, r, y, se, n, n1, n2, test_type, mu_transform){

  data <- list()

  # some general input checks
  if(sum(c(!is.null(t), !is.null(d), !is.null(r), !is.null(y))) != 1)stop("One effect size measure needs to be specified.")
  if(is.null(r) & !is.null(mu_transform))stop("'mu_transform' is available only if correlations are supplied as input.")
  if(!is.null(r))if(!mu_transform %in% c("cohens_d","fishers_z"))stop("'mu_transform' must be either 'cohens_d' or 'fishers_z'")
  if(!is.null(test_type))if(!test_type %in% c("one.sample","two.sample"))stop("'test_type' must be either 'one.sample' or 'two.sample'.")

  # check for NA or Inf
  if(!is.null(t))if(any(is.na(t))   | any(is.infinite(t)))stop("NAs or Inf are not allowed in 't'.")
  if(!is.null(d))if(any(is.na(d))   | any(is.infinite(d)))stop("NAs or Inf are not allowed in 'd'.")
  if(!is.null(r))if(any(is.na(r))   | any(is.infinite(r)))stop("NAs or Inf are not allowed in 'r'.")
  if(!is.null(y))if(any(is.na(y))   | any(is.infinite(y)))stop("NAs or Inf are not allowed in 'y'.")
  if(!is.null(se))if(any(is.na(se)) | any(is.infinite(se)))stop("NAs or Inf are not allowed in 'se'.")
  if(!is.null(n))if(any(is.na(n))   | any(is.infinite(n)))stop("NAs or Inf are not allowed in 'n'.")
  if(!is.null(n1))if(any(is.na(n1)) | any(is.infinite(n1)))stop("NAs or Inf are not allowed in 'n1'.")
  if(!is.null(n2))if(any(is.na(n2)) | any(is.infinite(n2)))stop("NAs or Inf are not allowed in 'n2'.")


  if((!is.null(t) | !is.null(d)) & ( !is.null(n) | (!is.null(n1) & !is.null(n2)) | !is.null(se) ) ){
    # for either t-statistics or Cohen's d and a sample size
    if(!is.null(t))if(!is.numeric(t)  & !is.vector(t))stop("'t' must be a numeric vector.")
    if(!is.null(d))if(!is.numeric(d)  & !is.vector(d))stop("'d' must be a numeric vector.")
    if(!is.null(n))if(!is.numeric(n)  & !is.vector(n))stop("'n' must be a numeric vector.")
    if(!is.null(n1))if(!is.numeric(n1) & !is.vector(n1))stop("'n1' must be a numeric vector.")
    if(!is.null(n2))if(!is.numeric(n2) & !is.vector(n2))stop("'n2' must be a numeric vector.")
    if(!is.null(n))if(!all(n > 1))stop("The sample sizes 'n' must be positive.")
    if(!is.null(n1))if(!all(n1 > 0))stop("The sample sizes 'n1' must be positive.")
    if(!is.null(n2))if(!all(n2 > 0))stop("The sample sizes 'n2' must be positive.")
    if(!is.null(se))if(!all(se > 0))stop("The standard errors 'se' must be positive.")

    # obtain test statistics and the ncp multiplactors
    if(!is.null(n)){
      data$df <- n - ifelse(test_type == "two.sample", 2, 1)
      # multiplicator for converting effect sizes into ncp
      if(test_type == "two.sample"){
        if(is.null(t))t <- psych::d2t(d = d, n = n)
        data$ncp_mlp    <- sqrt(n)/2
      }else if(test_type == "one.sample"){
        if(is.null(t))t <- psych::d2t(d = d, n1 = n)
        data$ncp_mlp    <- sqrt(n)
      }
    }else if(!is.null(n1) & !is.null(n2)){
      if(is.null(t))t <- psych::d2t(d = d, n1 = n1, n2 = n2)
      data$ncp_mlp    <- 1/sqrt(1/n1 + 1/n2)
      data$df         <- n1 + n2 - ifelse(test_type == "two.sample", 2, 1)
    }else if(!is.null(se)){
      n               <- .get_n_for_d(d, se)
      if(is.null(t))t <- d/se
      data$df         <- n - 2
      data$ncp_mlp    <- sqrt(n)/2
    }

    data$t <- t
    data$K <- length(t)

    data$effect_size <- "d"

  }else if(!is.null(r) & !is.null(n)){

    # using correlation and sample sizes
    if(!is.numeric(r) & !is.vector(r))stop("'r' must be a numeric vector.")
    if(!is.numeric(n) & !is.vector(n))stop("'n' must be a numeric vector.")
    if(!all(n > 1))stop("The sample sizes 'n' must be positive.")
    if(!(all(r > -1) & all(r < 1)))stop("The correlation coefficients 'r' must be in range (-1, 1).")

    if(mu_transform == "cohens_d"){

      # convert to cohen's d and compute it's test statistic
      # using n-2 leads to the actual t-statistic corresponding to the cor.test
      data$t       <- psych::d2t(psych::r2d(r), n = n-2)
      data$ncp_mlp <- sqrt(n-2)/2
      data$df      <- n - 2
      data$K       <- length(n)

    }else if(mu_transform == "fishers_z"){

      # convert to fisher's z and compute it's test statistic
      data$t       <- psych::fisherz(r) / (1/sqrt(n-3))
      data$ncp_mlp <- sqrt(n-3)
      data$df      <- n - 2
      data$K       <- length(n)

    }
    data$effect_size <- "r"

  }else if(!is.null(y) & !is.null(se)){
    # for an effect sizes and standard errors
    if(is.numeric(y)  & !is.vector(y))stop("'y' must be a numeric vector.")
    if(is.numeric(se) & !is.vector(se))stop("'se' must be a numeric vector.")
    if(!all(se > 0))stop("The standard errors 'se' must be positive.")

    # convert to z-statistics (and treat as t-statistics with df = Inf followingly)
    data$K       <- length(y)
    data$t       <- y / se
    data$df      <- rep(999999, length(y)) # should be Inf, but JAGS have problem with that for some reason
    data$ncp_mlp <- 1 / se

    data$effect_size <- "y"

  }else{
    stop("Insufficient input provided. Specify either the 't' / 'sd' and 'n' / 'n1' & 'n2' / 'se', or, 'r' and 'n', or 'y' and 'se'.")
  }


  return(data)
}


### fitting function
.fit_RoBMA             <- function(object, i){

  priors   <- object$models[[i]]$priors
  control  <- object$control
  new_warn <- NULL

  # don't sample the complete null model
  if(!(priors$mu$distribution == "point" &
       (if(priors$tau$distribution == "point"){priors$tau$parameters$location == 0}else{FALSE}) &
       priors$omega$distribution == "point")){

    # genrate the model syntax
    model_syntax <- .generate_model_syntax(priors)


    # remove unneccessary objects from data to mittigate warnings
    fit_data          <- .fit_data(object$data, priors)
    fit_inits         <- .fit_inits(priors, control$chains, object$add_info$seed)
    monitor_variables <- .to_monitor(priors)


    # fit the model - it would be nice to ignore the impression warnings
    fit <- tryCatch(runjags::autorun.jags(
      model           = model_syntax,
      data            = fit_data,
      inits           = fit_inits,
      monitor         = monitor_variables,
      n.chains        = control$chains,
      startburnin     = control$burnin,
      startsample     = control$iter,
      adapt           = control$adapt,
      thin            = control$thin,
      raftery.options = if(control$autofit) list(r = control$max_error) else FALSE,
      max.time        = if(control$autofit) control$max_time else Inf,
      summarise       = TRUE
    ), error = function(e)e)

    # deal with some fixable errors
    if(all(class(fit) %in% c("simpleError", "error", "condition"))){

      # create a new, data-tuned starting values if there is an outlier that fails the sampling
      if(grepl("Node inconsistent with parents", fit$message, fixed = TRUE) & any(names(unlist(fit_inits)) %in% c("mu", "inv_mu"))){

        new_mu <- mean(psych::t2d(fit_data$t, fit_data$df))
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

        new_warn <- paste0("Model's ",i," initial fit failed due to incompatible starting values (most likely due to an outlier in the data and limited precission of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data.")

        fit <- tryCatch(runjags::autorun.jags(
          model           = model_syntax,
          data            = fit_data,
          inits           = fit_inits,
          monitor         = monitor_variables,
          n.chains        = control$chains,
          startburnin     = control$burnin,
          startsample     = control$iter,
          adapt           = control$adapt,
          thin            = control$thin,
          raftery.options = if(control$autofit) list(r = control$max_error) else FALSE,
          max.time        = if(control$autofit) control$max_time else Inf,
          summarise       = TRUE
        ), error = function(e)e)

      }

    }


  }else{
    fit  <- NULL
  }


  object$models[[i]]$fit <- fit
  if(!is.null(fit)){
    object$models[[i]]$fit_summary <- summary(fit)
  }
  object$add_info$warnings <- c(object$add_info$warnings, new_warn)

  return(object)

}
.marglik_RoBMA         <- function(object, i){

  fit      <- object$models[[i]]$fit
  priors   <- object$models[[i]]$priors
  new_warn <- NULL

  # don't sample the complete null model
  if(!(priors$mu$distribution == "point" &
       (if(priors$tau$distribution == "point"){priors$tau$parameters$location == 0}else{FALSE}) &
       priors$omega$distribution == "point")){


    # deal failed model
    if(any(class(fit) %in% c("simpleError", "error"))){
      new_warn <- paste0("Model ",i," failed with the following error: ", fit$message)

      object$models[[i]]$marg_lik <- .marglik_fail()
      object$add_info$warnings    <- c(object$add_info$warnings, new_warn)

      return(object)
    }


    # compute marginal likelihood
    marglik_samples <- .marglik_prepare_data(fit, priors, object$data)
    fit_data        <- .fit_data(object$data, priors)

    if(!is.null(object$control$seed))set.seed(object$control$seed)
    marg_lik        <- tryCatch(bridgesampling::bridge_sampler(
      samples       = marglik_samples$samples,
      data          = fit_data,
      log_posterior = .marglik_function,
      priors        = priors,
      lb            = marglik_samples$lb,
      ub            = marglik_samples$ub,
      maxiter       = object$control$bridge_max_iter,
      silent        = TRUE),
      error = function(e)return(e))

  }else{
    # easy calculation of the marginal likelihood in case of null model
    marg_lik <- .marglik_null(object$data)
  }

  # handle errors
  if(any(class(marg_lik) %in% c("simpleError", "error"))){

    new_warn <- paste0("Marginal likelihood computation of model ",i," failed with the following error: ", marg_lik$message)
    marg_lik <- .marglik_fail()

  }else if(is.na(marg_lik$logml)){

    new_warn <- paste0("Marginal likelihood computation of model ",i," couldn't be completed within the specified number of iterations.")
    marg_lik <- .marglik_fail()

  }

  object$models[[i]]$marg_lik <- marg_lik
  object$add_info$warnings    <- c(object$add_info$warnings, new_warn)

  return(object)
}
.generate_model_syntax <- function(priors){

  # generate model syntax
  model_syntax <- "model{"

  ### mu priors
  # distributions
  if(priors$mu$distribution == "point"){
    model_syntax <- paste0(model_syntax, "mu = prior_mu_location\n")
  }else if(priors$mu$distribution == "normal"){
    model_syntax <- paste0(model_syntax, "mu ~ dnorm(prior_mu_mean, pow(prior_mu_sd, -2))")
  }else if(priors$mu$distribution == "t"){
    model_syntax <- paste0(model_syntax, "mu ~ dt(prior_mu_location, pow(prior_mu_scale, -2), prior_mu_df)")
  }else if(priors$mu$distribution == "gamma"){
    model_syntax <- paste0(model_syntax, "mu ~ dgamma(prior_mu_shape, prior_mu_rate)")
  }else if(priors$mu$distribution == "invgamma"){
    model_syntax <- paste0(model_syntax, "inv_mu ~ dgamma(prior_mu_shape, prior_mu_scale)")
  }else if(priors$mu$distribution == "uniform"){
    model_syntax <- paste0(model_syntax, "mu ~ dunif(prior_mu_a, prior_mu_b)")
  }
  # truncation
  if(!priors$mu$distribution %in% c("point", "uniform")){
    if(!(is.infinite(priors$mu$truncation$lower)  & is.infinite(priors$mu$truncation$lower))){
      # the truncation for invgamma needs to be done the other way around since we sample from gamma
      if(priors$mu$distribution == "invgamma"){
        model_syntax <- paste0(model_syntax, "T(",
                               ifelse(is.infinite(priors$mu$truncation$upper^-1),"",priors$mu$truncation$upper^-1),
                               ",",
                               ifelse(is.infinite(priors$mu$truncation$lower^-1),"",priors$mu$truncation$lower^-1),
                               ")\n")
      }else{
        model_syntax <- paste0(model_syntax, "T(",
                               ifelse(is.infinite(priors$mu$truncation$lower),"",priors$mu$truncation$lower),
                               ",",
                               ifelse(is.infinite(priors$mu$truncation$upper),"",priors$mu$truncation$upper),
                               ")\n")
      }
    }else{
      model_syntax <- paste0(model_syntax, "\n")
    }
  }

  # transformations
  if(priors$mu$distribution == "invgamma"){
    model_syntax <- paste0(model_syntax, "mu = pow(inv_mu, -1)\n")
  }


  ### tau priors
  # distributions
  if(priors$tau$distribution == "point"){
    if(priors$tau$parameters$location > 0)model_syntax <- paste0(model_syntax, "tau = prior_tau_location\n")
  }else if(priors$tau$distribution == "normal"){
    model_syntax <- paste0(model_syntax, "tau ~ dnorm(prior_tau_mean, pow(prior_tau_sd, -2))")
  }else if(priors$tau$distribution == "t"){
    model_syntax <- paste0(model_syntax, "tau ~ dt(prior_tau_location, pow(prior_tau_scale, -2), prior_tau_df)")
  }else if(priors$tau$distribution == "gamma"){
    model_syntax <- paste0(model_syntax, "tau ~ dgamma(prior_tau_shape, prior_tau_rate)")
  }else if(priors$tau$distribution == "invgamma"){
    model_syntax <- paste0(model_syntax, "inv_tau ~ dgamma(prior_tau_shape, prior_tau_scale)")
  }else if(priors$tau$distribution == "uniform"){
    model_syntax <- paste0(model_syntax, "tau ~ dunif(prior_tau_a, prior_tau_b)")
  }
  # truncation
  if(!priors$tau$distribution %in% c("point", "uniform")){
    if(!(is.infinite(priors$tau$truncation$lower) & is.infinite(priors$tau$truncation$lower))){
      # the truncation for invgamma needs to be done the other way around since we sample from gamma
      if(priors$tau$distribution == "invgamma"){
        model_syntax <- paste0(model_syntax, "T(",
                               ifelse(is.infinite(priors$tau$truncation$upper^-1),"",priors$tau$truncation$upper^-1),
                               ",",
                               ifelse(is.infinite(priors$tau$truncation$lower^-1),"",priors$tau$truncation$lower^-1),
                               ")\n")
      }else{
        model_syntax <- paste0(model_syntax, "T(",
                               ifelse(is.infinite(priors$tau$truncation$lower),"",priors$tau$truncation$lower),
                               ",",
                               ifelse(is.infinite(priors$tau$truncation$upper),"",priors$tau$truncation$upper),
                               ")\n")
      }
    }
  }
  # transformations
  if(priors$tau$distribution == "invgamma"){
    model_syntax <- paste0(model_syntax, "tau = pow(inv_tau, -1)\n")
  }


  ### omega priors
  # distributions & transformations
  if(priors$omega$distribution != "point"){
    if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){
      model_syntax <- paste0(model_syntax,
                             "for(j in 1:J){
                                 eta[j] ~ dgamma(prior_omega_alpha[j], 1)
                              }
                              for(j in 1:J){
                                 std_eta[j]  = eta[j] / sum(eta)
                                 omega[j]    = sum(std_eta[1:j])
                              }\n")
    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){
      model_syntax <- paste0(model_syntax,
                             "for(j1 in 1:J1){
                                 eta1[j1] ~ dgamma(prior_omega_alpha1[j1], 1)
                              }
                              for(j2 in 1:J2){
                                 eta2[j2] ~ dgamma(prior_omega_alpha2[j2], 1)
                              }
                              for(j1 in 1:J1){
                                 std_eta1[j1]  <-  eta1[j1] / sum(eta1)
                                 omega[J2 - 1 + j1] = sum(std_eta1[1:j1])
                              }
                              for(j2 in 1:J2){
                                  std_eta2[j2]  = (eta2[j2] / sum(eta2)) * (1 - std_eta1[1])
                              }
                              for(j2 in 2:J2){
                                  omega[j2-1] = sum(std_eta2[j2:J2]) + std_eta1[1]
                              }\n")
    }
  }


  ### model
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")
  if(priors$tau$distribution != "point"){
    model_syntax <- paste0(model_syntax, "theta[i] ~ dnorm(mu, pow(tau, -2))\n")
  }else if(priors$tau$parameters$location > 0){
    model_syntax <- paste0(model_syntax, "theta[i] ~ dnorm(mu, pow(tau, -2))\n")
  }
  if(priors$omega$distribution == "point"){

    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location > 0){
        model_syntax <- paste0(model_syntax, "t[i] ~ dnt(theta[i]*ncp_mlp[i], 1, df[i])\n")
      }else{
        model_syntax <- paste0(model_syntax, "t[i] ~ dnt(mu*ncp_mlp[i], 1, df[i])\n")
      }
    }else{
      model_syntax <- paste0(model_syntax, "t[i] ~ dnt(theta[i]*ncp_mlp[i], 1, df[i])\n")
    }

  }else if(priors$omega$distribution == "one.sided"){

    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location > 0){
        model_syntax <- paste0(model_syntax, "t[i] ~ dwt_1s(df[i], theta[i]*ncp_mlp[i], crit_t[i,], omega) \n")
      }else{
        model_syntax <- paste0(model_syntax, "t[i] ~ dwt_1s(df[i], mu*ncp_mlp[i], crit_t[i,], omega) \n")
      }
    }else{
      model_syntax <- paste0(model_syntax, "t[i] ~ dwt_1s(df[i], theta[i]*ncp_mlp[i], crit_t[i,], omega) \n")
    }

  }else if(priors$omega$distribution == "two.sided"){

    if(priors$tau$distribution == "point"){
      if(priors$tau$parameters$location > 0){
        model_syntax <- paste0(model_syntax, "t[i] ~ dwt_2s(df[i], theta[i]*ncp_mlp[i], crit_t[i,], omega) \n")
      }else{
        model_syntax <- paste0(model_syntax, "t[i] ~ dwt_2s(df[i], mu*ncp_mlp[i], crit_t[i,], omega) \n")
      }
    }else{
      model_syntax <- paste0(model_syntax, "t[i] ~ dwt_2s(df[i], theta[i]*ncp_mlp[i], crit_t[i,], omega) \n")
    }

  }
  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.fit_data              <- function(data, priors){

  # remove unneccessary stuff
  data$effect_size <- NULL

  ### add settings for prior distribution
  for(var in names(priors)){

    # don't add parameters for null omega or null tau
    if(var == "omega"){
      if(priors[[var]]$distribution == "point")next
    }
    if(var == "tau"){
      if(priors[[var]]$distribution == "point"){
        if(priors[[var]]$parameters$location <= 0)next
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
    temp_init <- c(temp_init, .fit_inits_mu_tau(priors$tau, "tau"))
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

  if(all(names(prior$parameters) %in% c("alpha", "steps"))){

    temp_x$eta <- round(stats::rgamma(length(prior$parameters$alpha),   prior$parameters$alpha,  1),5)

  }else if(all(names(prior$parameters) %in% c("alpha1", "alpha2", "steps"))){

    temp_x$eta1 <- round(stats::rgamma(length(prior$parameters$alpha1), prior$parameters$alpha1, 1),5)
    temp_x$eta2 <- round(stats::rgamma(length(prior$parameters$alpha2), prior$parameters$alpha2, 1),5)

  }


  return(temp_x)

}
.to_monitor            <- function(priors){

  variables <- NULL

  # mu relevant
  if(priors$mu$distribution != "point"){
    if(priors$mu$distribution == "invgamma"){
      variables <- c(variables, "inv_mu")
    }
    variables <- c(variables, "mu")
  }

  # tau relevant
  if(priors$tau$distribution != "point"){
    if(priors$tau$distribution == "invgamma"){
      variables <- c(variables, "inv_tau")
    }
    variables <- c(variables, "tau", "theta")
  }else if(priors$tau$parameters$location > 0){
    variables <- c(variables, "theta")
  }

  # omega relevant
  if(priors$omega$distribution != "point"){
    if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){
      variables <- c(variables, "eta", "omega")
    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){
      variables <- c(variables, "eta1", "eta2", "omega")
    }
  }

  return(variables)
}
.marglik_function      <- function(samples.row, data, priors){

  ### get parameteres depending on the model type
  # mu
  if(priors$mu$distribution != "point"){
    if(priors$mu$distribution == "invgamma"){
      inv_mu <- samples.row[[ "inv_mu" ]]
      mu     <- 1/inv_mu
    }else{
      mu <- samples.row[[ "mu" ]]
    }
  }else{
    mu <- 0
  }

  # tau
  if(priors$tau$distribution == "point"){
    if(priors$tau$parameters$location != 0){
      tau      <- priors$tau$parameters$location
      theta    <- samples.row[ paste0("theta[", 1:data$K, "]") ]
    }
    ncp <- mu*data$ncp_mlp
  }else{
    theta    <- samples.row[ paste0("theta[", 1:data$K, "]") ]
    if(priors$tau$distribution == "invgamma"){
      inv_tau <- samples.row[[ "inv_tau" ]]
      tau     <- 1/inv_tau
    }else{
      tau <- samples.row[[ "tau" ]]
    }
    ncp <- theta*data$ncp_mlp
  }

  # omega
  if(priors$omega$distribution != "point"){
    if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){

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


  ### compute the marginal log_likelihood
  log_lik <- 0

  # mean
  if(priors$mu$distribution == "normal"){
    log_lik <- log_lik + stats::dnorm(mu, mean = data$prior_mu_mean, sd = data$prior_mu_sd, log = TRUE) -
      log(
        stats::pnorm(priors$mu$truncation$upper, data$prior_mu_mean, data$prior_mu_sd, lower.tail = TRUE, log.p = FALSE) -
          stats::pnorm(priors$mu$truncation$lower, data$prior_mu_mean, data$prior_mu_sd, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$mu$distribution == "t"){
    log_lik <- log_lik + extraDistr::dlst(mu, df = data$prior_mu_df, mu = data$prior_mu_location, sigma = data$prior_mu_scale, log = TRUE) -
      log(
        extraDistr::plst(priors$mu$truncation$upper, df = data$prior_mu_df, mu = data$prior_mu_location, sigma = data$prior_mu_scale, lower.tail = TRUE, log.p = FALSE) -
          extraDistr::plst(priors$mu$truncation$lower, df = data$prior_mu_df, mu = data$prior_mu_location, sigma = data$prior_mu_scale, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$mu$distribution == "gamma"){
    log_lik <- log_lik + stats::dgamma(mu, shape = data$prior_mu_shape, rate = data$prior_mu_rate, log = TRUE)  -
      log(
        stats::pgamma(priors$mu$truncation$upper, shape = data$prior_mu_shape, rate = data$prior_mu_rate, lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(priors$mu$truncation$lower, shape = data$prior_mu_shape, rate = data$prior_mu_rate, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$mu$distribution == "invgamma"){
    log_lik <- log_lik + stats::dgamma(inv_mu, shape = data$prior_mu_shape, rate = data$prior_mu_scale, log = TRUE) -
      log(
        stats::pgamma(priors$mu$truncation$lower^-1, shape = data$prior_mu_shape, rate = data$prior_mu_scale, lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(priors$mu$truncation$upper^-1, shape = data$prior_mu_shape, rate = data$prior_mu_scale, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$mu$distribution == "uniform"){
    log_lik <- log_lik + stats::dunif(mu, min = data$prior_mu_a, max = data$prior_mu_b, log = TRUE)
  }


  # tau
  if(priors$tau$distribution == "normal"){
    log_lik <- log_lik + stats::dnorm(tau, mean = data$prior_tau_mean, sd = data$prior_tau_sd, log = TRUE) -
      log(
        stats::pnorm(priors$tau$truncation$upper, data$prior_tau_mean, data$prior_tau_sd, lower.tail = TRUE, log.p = FALSE) -
          stats::pnorm(priors$tau$truncation$lower, data$prior_tau_mean, data$prior_tau_sd, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$tau$distribution == "t"){
    log_lik <- log_lik + extraDistr::dlst(tau, df = data$prior_tau_df, mu = data$prior_tau_location, sigma = data$prior_tau_scale, log = TRUE) -
      log(
        extraDistr::plst(priors$tau$truncation$upper, df = data$prior_tau_df, mu = data$prior_tau_location, sigma = data$prior_tau_scale, lower.tail = TRUE, log.p = FALSE) -
          extraDistr::plst(priors$tau$truncation$lower, df = data$prior_tau_df, mu = data$prior_tau_location, sigma = data$prior_tau_scale, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$tau$distribution == "gamma"){
    log_lik <- log_lik + stats::dgamma(tau, shape = data$prior_tau_shape, rate = data$prior_tau_rate, log = TRUE) -
      log(
        stats::pgamma(priors$tau$truncation$upper, shape = data$prior_tau_shape, rate = data$prior_tau_rate, lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(priors$tau$truncation$lower, shape = data$prior_tau_shape, rate = data$prior_tau_rate, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$tau$distribution == "invgamma"){
    log_lik <- log_lik + stats::dgamma(inv_tau, shape = data$prior_tau_shape, rate = data$prior_tau_scale, log = TRUE) -
      log(
        stats::pgamma(priors$tau$truncation$lower^-1, shape = data$prior_tau_shape, rate = data$prior_tau_scale, lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(priors$tau$truncation$upper^-1, shape = data$prior_tau_shape, rate = data$prior_tau_scale, lower.tail = TRUE, log.p = FALSE)
      )
  }else if(priors$tau$distribution == "uniform"){
    log_lik <- log_lik + stats::dunif(tau, min = data$prior_tau_a, max = data$prior_tau_b, log = TRUE)
  }


  # omega
  if(priors$omega$distribution != "point"){
    if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){

      log_lik <- log_lik + sum(stats::dgamma(eta, data$prior_omega_alpha, 1, log = TRUE))

    }else if(all(names(priors$omega$parameters) %in% c("alpha1", "alpha2", "steps"))){

      log_lik <- log_lik + sum(stats::dgamma(eta1, data$prior_omega_alpha1, 1, log = TRUE))
      log_lik <- log_lik + sum(stats::dgamma(eta2, data$prior_omega_alpha2, 1, log = TRUE))

    }
  }


  # the true effect sizes (in case of heterogeneity)
  if(priors$tau$distribution == "point"){
    if(priors$tau$parameters$location != 0){
      log_lik <- log_lik + sum(stats::dnorm(theta, mean = mu, sd = tau, log = TRUE))
    }
  }else{
    log_lik <- log_lik + sum(stats::dnorm(theta, mean = mu, sd = tau, log = TRUE))
  }


  # the observed t-statistics
  if(priors$omega$distribution == "point"){

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

  return(log_lik)
}
.marglik_prepare_data  <- function(fit, priors, data){

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
  if(priors$tau$distribution == "point"){
    if(priors$tau$parameters$location != 0){
      pars <- c(pars, paste0("theta[", 1:data$K, "]"))
      lb   <- c(lb, rep(-Inf, data$K))
      ub   <- c(ub, rep( Inf, data$K))
    }
  }else{
    pars <- c(pars, paste0("theta[", 1:data$K, "]"))
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


  # omega
  if(priors$omega$distribution != "point"){

    if(all(names(priors$omega$parameters) %in% c("alpha", "steps"))){
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
.marglik_null          <- function(data){

  marg_lik        <- NULL
  marg_lik$logml  <- sum(stats::dt(data$t, data$df, 0, log = TRUE))
  class(marg_lik) <- "bridge"

  return(marg_lik)
}
.marglik_fail          <- function(){
  marg_lik        <- NULL
  marg_lik$logml  <- -Inf
  class(marg_lik) <- "bridge"
  return(marg_lik)
}


### model inference functions
.model_inference       <- function(models, data, converged, mu_transform, seed, n_samples = 10000){

  # extract marginal likelihoods
  marg_liks <- sapply(models, function(x)x$marg_lik$logml)

  # determine the type of the models
  mm_mu     <- sapply(models, function(m)!.is_parameter_null(m$priors, "mu"))
  mm_tau    <- sapply(models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega  <- sapply(models, function(m)!.is_parameter_null(m$priors, "omega"))

  # extract model weights
  prior_weights_all   <- sapply(models, function(m)m$prior_odds)
  prior_weights_mu    <- ifelse(mm_mu,    prior_weights_all, 0)
  prior_weights_tau   <- ifelse(mm_tau,   prior_weights_all, 0)
  prior_weights_omega <- ifelse(mm_omega, prior_weights_all, 0)

  # standardize model weights
  prior_weights_all   <- prior_weights_all   / sum(prior_weights_all)
  prior_weights_mu    <- prior_weights_mu    / sum(prior_weights_mu)
  prior_weights_tau   <- prior_weights_tau   / sum(prior_weights_tau)
  prior_weights_omega <- prior_weights_omega / sum(prior_weights_omega)


  ### compute model weights
  # overall
  weights_all   <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_all)
  if(any(mm_mu) & all(!is.nan(prior_weights_mu))){
    weights_mu  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_mu)
  }else{
    weights_mu <- NULL
  }
  if(any(mm_tau) & all(!is.nan(prior_weights_tau))){
    weights_tau <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_tau)
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
  BF_heterogeneity <- .inclusion_BF(prior_weights_all, weights_all, mm_tau)
  BF_bias          <- .inclusion_BF(prior_weights_all, weights_all, mm_omega)


  ### sample and mix the individual posteriors
  if(!is.null(seed))set.seed(seed)
  samples <- NULL
  samples$averaged    <- .mix_samples(models, weights_all, data, converged, mu_transform, n_samples)
  samples$conditional <- .mix_samples(models, list(mu = weights_mu, tau = weights_tau, omega = weights_omega), data, converged, mu_transform, n_samples)


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
      heterogeneity = BF_heterogeneity,
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
.mix_samples           <- function(models, weights, data, converged, mu_transform, n_samples){

  # metadata about model type
  mm_mu     <- sapply(models, function(m)!.is_parameter_null(m$priors, "mu"))
  mm_tau    <- sapply(models, function(m)!.is_parameter_null(m$priors, "tau"))
  mm_omega  <- sapply(models, function(m)!.is_parameter_null(m$priors, "omega"))

  # the indicies for omega based on type of p-value cutpoints
  omega_ind  <- .get_omega_mapping(models)

  # prepare the object
  samples = list(
    mu    = NULL,
    tau   = NULL,
    omega = matrix(nrow = 0, ncol = if(is.null(omega_ind)) 0 else ncol(do.call(rbind, omega_ind))),
    theta = matrix(nrow = 0, ncol = data$K)
  )


  for(i in c(1:length(models))[converged]){

    # deal with the complete null model possibility
    if(models[[i]]$priors$mu$distribution == "point" &
         (if(models[[i]]$priors$tau$distribution == "point"){models[[i]]$priors$tau$parameters$location == 0}else{FALSE}) &
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

      # tau
      samples$tau <- c(samples$tau,
                       if(models[[i]]$priors$tau$distribution == "point"){
                         rep(models[[i]]$priors$tau$parameter$location, round(n_samples * weights[i]))
                       }else{
                         model_samples[ind, "tau"]
                       })

      # theta
      samples$theta <- rbind(samples$theta,
                             if(models[[i]]$priors$tau$distribution == "point"){
                               if(models[[i]]$priors$tau$parameter$location == 0){
                                 if(models[[i]]$priors$mu$distribution == "point"){
                                   matrix(models[[i]]$priors$mu$parameter$location, ncol = data$K, nrow = round(n_samples * weights[i]))
                                 }else{
                                   matrix(model_samples[ind,"mu"], ncol = data$K, nrow = round(n_samples * weights[i]))
                                 }
                               }else{
                                 model_samples[ind, paste0("theta[",1:data$K,"]")]
                               }
                             }else{
                               model_samples[ind, paste0("theta[",1:data$K,"]")]
                             })

      # omega
      samples$omega <- rbind(samples$omega,
                             if(models[[i]]$priors$omega$distribution == "point"){
                               matrix(1, ncol = ncol(samples$omega), nrow = round(n_samples * weights[i]))
                             }else{
                               model_samples[ind, paste0("omega[",omega_ind[[i]],"]")]
                             })


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


      # theta
      if(!is.null(weights$mu[i])){
        samples$theta <- rbind(samples$theta,
                               if(mm_tau[[i]]){
                                 if(models[[i]]$priors$tau$distribution == "point"){
                                   if(models[[i]]$priors$tau$parameter$location == 0){
                                     if(models[[i]]$priors$mu$distribution == "point"){
                                       matrix(models[[i]]$priors$mu$parameter$location, ncol = data$K, nrow = round(n_samples * weights$mu[i]))
                                     }else{
                                       matrix(model_samples[sample(nrow(model_samples), round(n_samples * weights$tau[i]), replace = TRUE),"mu"], ncol = data$K, nrow = round(n_samples * weights$mu[i]))
                                     }
                                   }else{
                                     model_samples[sample(nrow(model_samples), round(n_samples * weights$mu[i]), replace = TRUE), paste0("theta[",1:data$K,"]")]
                                   }
                                 }else{
                                   model_samples[sample(nrow(model_samples), round(n_samples * weights$mu[i]), replace = TRUE), paste0("theta[",1:data$K,"]")]
                                 }
                               })
      }


      # omega
      if(!is.null(weights$omega[i])){
        samples$omega <- rbind(samples$omega,
                               if(mm_omega[[i]]){
                                 if(models[[i]]$priors$omega$distribution == "point"){
                                   matrix(1, ncol = ncol(samples$omega), nrow = round(n_samples * weights$omega[i]))
                                 }else{
                                   model_samples[sample(nrow(model_samples), round(n_samples * weights$omega[i]), replace = TRUE), paste0("omega[",omega_ind[[i]],"]")]
                                 }
                               })
      }


    }

  }


  # convert the transformed correlations back if needed
  if(!isFALSE(mu_transform)){
    if(mu_transform == "cohens_d"){
      samples$mu    <- psych::d2r(samples$mu)
      samples$theta <- psych::d2r(samples$theta)
    }else if(mu_transform == "fishers_z"){
      samples$mu    <- psych::fisherz2r(samples$mu)
      samples$theta <- psych::fisherz2r(samples$theta)
    }
  }

  # fix omega names
  all_cuts                <- .get_omega_mapping(models, cuts_only = TRUE)
  if(!is.null(all_cuts))colnames(samples$omega) <- sapply(1:(length(all_cuts)-1), function(i)paste0("omega[",all_cuts[i],",",all_cuts[i+1],"]"))
  # fix theta names
  colnames(samples$theta) <- paste0("theta[", 1:ncol(samples$theta), "]")

  return(samples)
}
.compute_coeficients   <- function(RoBMA){
  return(c(
    "mu"     = if(length(RoBMA$samples$averaged$mu) != 0)mean(RoBMA$samples$averaged$mu),
    "tau"    = if(length(RoBMA$samples$averaged$tau) != 0)mean(RoBMA$samples$averaged$tau),
    if(ncol(RoBMA$samples$averaged$omega) != 0)apply(RoBMA$samples$averaged$omega, 2, mean)
  ))
}
.inclusion_BF          <- function(prior_weights, posterior_weights, conditional_models){
  (sum(posterior_weights[conditional_models])/sum(posterior_weights[!conditional_models]))  /
    (sum(prior_weights[conditional_models])/sum(prior_weights[!conditional_models]))
}
.get_converged_models  <- function(object){

  converged <- NULL

  # basic convergence checks
  for(i in 1:length(object$models)){
    if(!(object$models[[i]]$priors$mu$distribution == "point" &
       (if(object$models[[i]]$priors$tau$distribution == "point"){object$models[[i]]$priors$tau$parameters$location == 0}else{FALSE}) &
       object$models[[i]]$priors$omega$distribution == "point")){

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
.balance_prob          <- function(object, converged_models){

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
    }else{
      prior_odds[i] <- 0
      warning("Prior probability couldn't be balanced over models with same combination of presence/absence of effect/heterogeneity/publication bias since they don't exist.")
    }
  }

  for(i in 1:length(object$models)){
    object$models[[i]]$prior_odds <- prior_odds[i]
  }

  return(object)
}


### helper functions for settings
.set_parameter_priors  <- function(priors_null, priors_alt, parameter){

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


  }else if(parameter == "tau"){

    # check that the passed priors are supported for the tau parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("normal", "t", "gamma", "invgamma", "point", "uniform"))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the tau parameter. See '?prior' for further information."))
        if(priors[[i]]$distribution == "point"){
          if(priors[[i]]$parameters$location < 0){
            stop("The location of a point prior distribution for tau parameter cannot be negative. See '?prior' for further information.")
          }
        }else if(priors[[i]]$distribution == "uniform"){
          if(priors[[i]]$parameters$a < 0 | priors[[i]]$parameters$b < 0 ){
            stop("The uniform prior distribution for tau parameter cannot be defined on negative numbers. See '?prior' for further information.")
          }
        }else if(priors[[i]]$truncation$lower < 0){
          priors[[i]]$truncation$lower <- 0
          warning("The range of a prior distribution for tau parameter cannot be negative. The lower truncation point was set to zero. See '?prior' for further information.")
        }
      }
    }


  }else if(parameter == "omega"){

    # check that the passed priors are supported for the parameter
    if(length(priors) > 0){
      for(i in 1:length(priors)){
        if(!priors[[i]]$distribution %in% c("two.sided", "one.sided", "point"))stop(paste0(priors[[i]]$distribution," prior distribution is not supported for the omega parameter. See '?prior' for further information."))
      }
    }
  }


  return(priors)

}
.get_models            <- function(priors){

  # create models according to the set priors
  models <- NULL
  for(mu in priors$mu){
    for(tau in priors$tau){
      for(omega in priors$omega){
        models <- c(
          models,
          list(.create_model(mu, tau, omega, mu$prior_odds * tau$prior_odds * omega$prior_odds))
        )
      }
    }
  }

  return(models)
}
.create_model          <- function(prior_mu, prior_tau, priors_omega, prior_odds){

  model <- list(
    priors = list(
      mu    = prior_mu,
      tau   = prior_tau,
      omega = priors_omega
    ),
    prior_odds     = prior_odds,
    prior_odds_set = prior_odds
  )
  class(model) <- "RoBMA.model"

  return(model)

}
.set_control           <- function(control, chains, iter, burnin, thin){

  # set the control list
  if(is.null(control)){
    control$max_error       <- .01
    control$max_time        <- Inf
    control$autofit         <- TRUE
    control$adapt           <- 1000
    control$bridge_max_iter <- 10000

    control$allow_max_error <- NULL
    control$allow_max_rhat  <- NULL
    control$allow_min_ESS   <- NULL
    control$allow_inc_theta <- FALSE
    control$balance_prob    <- TRUE

    control$silent          <- FALSE

  }else{
    if(is.null(control[["max_error"]])){
      control$max_error       <- .01
    }
    if(is.null(control[["max_time"]])){
      control$max_time        <- Inf
    }
    if(is.null(control[["autofit"]])){
      control$autofit         <- TRUE
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
  }

  # add the main MCMC settings
  control$chains    <- chains
  control$iter      <- iter
  control$burnin    <- burnin
  control$thin      <- thin

  .check_control(control)
  return(control)
}
.update_control        <- function(control, control_new, chains, iter, burnin, thin){

  if(!is.null(control_new)){
    for(n in names(control_new)){
      control[[n]] <- control_new[[n]]
    }
  }

  if(!is.null(chains))  control[["chains"]] <- chains
  if(!is.null(iter))    control[["iter"]]   <- iter
  if(!is.null(burnin))  control[["burnin"]] <- burnin
  if(!is.null(thin))    control[["thin"]]   <- thin

  # stop if there is not enough samples planned for autojags package
  .check_control(control)

  return(control)
}
.check_control         <- function(control){
  # check whether only known controls were supplied
  known_controls <- c("chains", "iter", "burnin" , "adapt", "thin" ,"autofit", "max_error", "max_time", "bridge_max_iter", "allow_max_error", "allow_max_rhat", "allow_min_ESS", "allow_inc_theta", "balance_prob", "silent", "progress_start", "progress_tick")
  if(any(!names(control) %in% known_controls))stop(paste0("The following control settings were not recognize: ", paste(names(control[!names(control) %in% known_controls]), collapse = ", ")))

  # check whether essential controls were supplied
  if(is.null(control[["chains"]])) stop("Number of chains must be defined.")
  if(is.null(control[["iter"]]))   stop("Number of iterations must be set.")
  if(is.null(control[["burnin"]])) stop("Number of burnin samples must be set.")
  if(is.null(control[["adapt"]]))  stop("Number of adaptation samples must be set.")
  if(is.null(control[["thin"]]))   stop("Thinning of the posterior samples must be set.")

  if(!is.numeric(control[["chains"]]) | !control[["chains"]] >= 1) stop("At least one chains must be set.")
  if(!is.numeric(control[["iter"]])   | !control[["iter"]] >= 1)   stop("Number of iterations must be a positive number")
  if(!is.numeric(control[["burnin"]]) | !control[["burnin"]] >= 1) stop("Number of burnin samples must be a positive number")
  if(!is.numeric(control[["adapt"]])  | !control[["adapt"]] >= 1)  stop("Number of adaptation samples must be a positive number.")
  if(!is.numeric(control[["thin"]])   | !control[["thin"]] >= 1)   stop("Thinning of the posterior samples must be a positive number")

  # stop if there is not enough samples planned for autojags package
  if(control$iter/control$thin < 4000)stop("At least 4000 iterations after thinning is required to compute the Raftery and Lewis's diagnostic.")
  if(!control$chains >= 2)stop("The number of chains must be at least 2 so that convergence can be assessed.")

  # check convergence criteria
  if(control$autofit)if(control$max_error >= 1 | control$max_error <= 0)stop("The target maximum MCMC error must be within 0 and 1.")
  if(!is.null(control$allow_max_error))if(control$allow_max_error >= 1 | control$allow_max_error <= 0)stop("The maximum allowed MCMC error must be within 0 and 1.")
  if(!is.null(control$allow_max_rhat))if(control$allow_max_rhat <= 1)stop("The maximum allowed R-hat must be higher than 1.")
  if(!is.null(control$allow_min_ESS))if(control$allow_min_ESS <= 0)stop("The minimum allowed ESS must be higher than 0.")

  # why can't I change both from the runjags interface directly?
  runjags::runjags.options(silent.jags = control$silent, silent.runjags = control$silent)
}


# general helper functions
.is_parameter_null <- function(priors, par){
  return(priors[[par]]$is_null)
}
.get_omega_mapping <- function(models, cuts_only = FALSE){

  # extract cuts and types
  p_cuts <- sapply(models, function(m)rev(m$priors$omega$parameters$steps), simplify = FALSE)
  p_type <- sapply(models, function(m)m$priors$omega$distribution)

  # remove point distributions
  if(all(p_type == "point"))return(NULL)

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
    if(p_type[p] != "point"){
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

  if(any(n < 1))stop("The computed sample size of at least on of the original studies was lower than 1 (based on the Cohen's d and standard error). This does not seem to be correct. Please, check your input.")

  return(n)
}
.get_effect_and_ci <- function(add_info, CI){

  # extract the information
  t  <- add_info$t
  d  <- add_info$d
  r  <- add_info$r
  y  <- add_info$y
  n  <- add_info$n
  n1 <- add_info$n1
  n2 <- add_info$n2
  se <- add_info$se
  effect_size <- add_info$effect_size
  test_type   <- add_info$test_type
  study_names <- add_info$study_names

  # compute the mean and CI
  if(effect_size == "d"){
    if(test_type == "one.sample"){
      out <- psych::d.ci(psych::t2d(t = t, n1 = n), n1 = n, alpha = 1 - CI)
    }else if(test_type == "two.sample"){
      if(!is.null(n)){
        out <- psych::d.ci(psych::t2d(t = t, n = n), n = n, alpha = 1 - CI)
      }else if(!is.null(n1) & !is.null(n2)){
        out <- psych::d.ci(psych::t2d(t = t, n1 = n1, n2 = n2), n1 = n1, n2 = n2, alpha = 1 - CI)
      }
    }
  }else if(effect_size == "r"){
    out <- matrix(suppressWarnings(psych::r.con(r = r, n = n, p = CI, twotailed = TRUE)), ncol = 2)
    out <- cbind(out[,1], r, out[,2])
  }else if(effect_size == "y"){
    out <- cbind(y + stats::qnorm((1-CI)/2) * se, y, y - stats::qnorm((1-CI)/2) * se)
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
