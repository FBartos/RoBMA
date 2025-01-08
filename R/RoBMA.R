#' @title Estimate a Robust Bayesian Meta-Analysis
#'
#' @description \code{RoBMA} is used to estimate a robust Bayesian
#' meta-analysis. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @param data a data object created by the \code{combine_data} function. This is
#' an alternative input entry to specifying the \code{d}, \code{r}, \code{y}, etc...
#' directly. I.e., RoBMA function does not allow passing a data.frame and
#' referencing to the columns.
#' @param model_type string specifying the RoBMA ensemble. Defaults to \code{NULL}.
#' The other options are \code{"PSMA"}, \code{"PP"}, and \code{"2w"} which override
#' settings passed to the \code{priors_effect}, \code{priors_heterogeneity},
#' \code{priors_effect}, \code{priors_effect_null}, \code{priors_heterogeneity_null},
#' \code{priors_bias_null}, and \code{priors_effect}. See details for more information
#' about the different model types.
#' @param effect_direction the expected direction of the effect. Correctly specifying
#' the expected direction of the effect is crucial for one-sided selection models,
#' as they specify cut-offs using one-sided p-values. Defaults to \code{"positive"}
#' (another option is \code{"negative"}).
#' @param prior_scale an effect size scale used to define priors. Defaults to \code{"cohens_d"}.
#' Other options are \code{"fishers_z"}, correlation coefficient \code{"r"},
#' and \code{"logOR"}. The prior scale does not need to match the effect sizes measure -
#' the samples from prior distributions are internally transformed to match the
#' \code{transformation} of the data. The \code{prior_scale} corresponds to
#' the effect size scale of default output, but can be changed within the summary function.
#' @param rescale_priors a re-scaling factor for the prior distributions. The re-scaling
#' factor allows to adjust the width of all default priors simultaneously. Defaults to \code{1}.
#' @param priors_effect list of prior distributions for the effect size (\code{mu})
#' parameter that will be treated as belonging to the alternative hypothesis. Defaults to
#' a standard normal distribution
#' \code{prior(distribution = "normal", parameters = list(mean = 0, sd = 1))}.
#' @param priors_heterogeneity list of prior distributions for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the alternative hypothesis. Defaults to
#' \code{prior(distribution = "invgamma", parameters = list(shape = 1, scale = .15))} that
#' is based on heterogeneities estimates from psychology \insertCite{erp2017estimates}{RoBMA}.
#' @param priors_bias list of prior distributions for the publication bias adjustment
#' component that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{list(
#' prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),
#'     steps = c(0.05)),             prior_weights = 1/12),
#' prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),
#'     steps = c(0.05, 0.10)),       prior_weights = 1/12),
#' prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1),
#'      steps = c(0.05)),             prior_weights = 1/12),
#' prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),
#'      steps = c(0.025, 0.05)),      prior_weights = 1/12),
#' prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1),
#'      steps = c(0.05, 0.5)),        prior_weights = 1/12),
#' prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1),
#'      steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
#' prior_PET(distribution   = "Cauchy", parameters = list(0,1), truncation = list(0, Inf),
#'      prior_weights = 1/4),
#' prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf),
#'      prior_weights = 1/4)
#' )}, corresponding to the RoBMA-PSMA model introduce by \insertCite{bartos2021no;textual}{RoBMA}.
#' @param priors_effect_null list of prior distributions for the effect size (\code{mu})
#' parameter that will be treated as belonging to the null hypothesis. Defaults to
#' a point null hypotheses at zero,
#' \code{prior(distribution = "point", parameters = list(location = 0))}.
#' @param priors_heterogeneity_null list of prior distributions for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the null hypothesis. Defaults to
#' a point null hypotheses at zero (a fixed effect meta-analytic models),
#' \code{prior(distribution = "point", parameters = list(location = 0))}.
#' @param priors_bias_null list of prior weight functions for the \code{omega} parameter
#' that will be treated as belonging to the null hypothesis. Defaults no publication
#' bias adjustment, \code{prior_none()}.
#' @param priors_hierarchical list of prior distributions for the correlation of random effects
#' (\code{rho}) parameter that will be treated as belonging to the alternative hypothesis. This setting allows
#' users to fit a hierarchical (three-level) meta-analysis when \code{study_ids} are supplied.
#' Note that this is an experimental feature and see News for more details. Defaults to a beta distribution
#' \code{prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))}.
#' @param priors_hierarchical_null list of prior distributions for the correlation of random effects
#' (\code{rho}) parameter that will be treated as belonging to the null hypothesis. Defaults to \code{NULL}.
#' @param algorithm a string specifying the algorithm used for the model averaging. Defaults to \code{"bridge"}
#' which results in estimating individual models using JAGS and computing the marginal likelihood using bridge
#' sampling. An alternative is \code{"ss"} which uses spike and slab like parameterization to approximate the
#' Bayesian model averaging with a single model.
#' @param sample a number of sampling iterations of the MCMC algorithm.
#' Defaults to \code{5000}.
#' @param burnin a number of burnin iterations of the MCMC algorithm.
#' Defaults to \code{2000}.
#' @param adapt a number of adaptation iterations of the MCMC algorithm.
#' Defaults to \code{500}.
#' @param thin a thinning of the chains of the MCMC algorithm. Defaults to
#' \code{1}.
#' @param parallel whether the individual models should be fitted in parallel.
#' Defaults to \code{FALSE}. The implementation is not completely stable
#' and might cause a connection error.
#' @param autofit whether the model should be fitted until the convergence
#' criteria (specified in \code{autofit_control}) are satisfied. Defaults to
#' \code{TRUE}.
#' @param autofit_control allows to pass autofit control settings with the
#' [set_autofit_control()] function. See \code{?set_autofit_control} for
#' options and default settings.
#' @param convergence_checks automatic convergence checks to assess the fitted
#' models, passed with [set_convergence_checks()] function. See
#' \code{?set_convergence_checks} for options and default settings.
#' @param save whether all models posterior distributions should be kept
#' after obtaining a model-averaged result. Defaults to \code{"all"} which
#' does not remove anything. Set to \code{"min"} to significantly reduce
#' the size of final object, however, some model diagnostics and further
#' manipulation with the object will not be possible.
#' @param seed a seed to be set before model fitting, marginal likelihood
#' computation, and posterior mixing for reproducibility of results. Defaults
#' to \code{NULL} - no seed is set.
#' @param silent whether all print messages regarding the fitting process
#' should be suppressed. Defaults to \code{TRUE}. Note that \code{parallel = TRUE}
#' also suppresses all messages.
#' @param ... additional arguments.
#' @inheritParams combine_data
#'
#' @details The default settings of the RoBMA 2.0 package corresponds to the RoBMA-PSMA
#' ensemble proposed by \insertCite{bartos2021no;textual}{RoBMA}. The previous versions
#' of the package (i.e., RoBMA < 2.0) used specifications proposed by
#' \insertCite{maier2020robust;textual}{RoBMA} (this specification can be easily
#' obtained by setting \code{model_type = "2w"}. The RoBMA-PP specification from
#' \insertCite{bartos2021no;textual}{RoBMA} can be obtained by setting
#' \code{model_type = "PP"}. The complete list of default prior distributions is described at
#' [set_default_priors()].
#'
#' The \href{../doc/CustomEnsembles.html}{\code{vignette("CustomEnsembles", package = "RoBMA")}}
#' and \href{../doc/ReproducingBMA.html}{\code{vignette("ReproducingBMA", package = "RoBMA")}}
#' vignettes  describe how to use [RoBMA()] to fit custom meta-analytic ensembles (see [prior()],
#' [prior_weightfunction()], [prior_PET()], and [prior_PEESE()] for more information about prior
#' distributions).
#'
#' The RoBMA function first generates models from a combination of the
#' provided priors for each of the model parameters. Then, the individual models
#' are fitted using \link[runjags]{autorun.jags} function. A marginal likelihood
#' is computed using \link[bridgesampling]{bridge_sampler} function. The individual
#' models are then combined into an ensemble using the posterior model probabilities
#' using \link[BayesTools]{BayesTools} package.
#'
#' Generic [summary.RoBMA()], [print.RoBMA()], and [plot.RoBMA()] functions are
#' provided to facilitate manipulation with the ensemble. A visual check of the
#' individual model diagnostics can be obtained using the [diagnostics()] function.
#' The fitted model can be further updated or modified by [update.RoBMA()] function.
#'
#' @examples \dontrun{
#' # using the example data from Bem 2011 and fitting the default (RoBMA-PSMA) model
#' fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study)
#'
#' # in order to speed up the process, we can turn the parallelization on
#' fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study, parallel = TRUE)
#'
#' # we can get a quick overview of the model coefficients just by printing the model
#' fit
#'
#' # a more detailed overview using the summary function (see '?summary.RoBMA' for all options)
#' summary(fit)
#'
#' # the model-averaged effect size estimate can be visualized using the plot function
#' # (see ?plot.RoBMA for all options)
#' plot(fit, parameter = "mu")
#'
#' # forest plot can be obtained with the forest function (see ?forest for all options)
#' forest(fit)
#'
#' # plot of the individual model estimates can be obtained with the plot_models function
#' #  (see ?plot_models for all options)
#' plot_models(fit)
#'
#' # diagnostics for the individual parameters in individual models can be obtained using diagnostics
#' # function (see 'diagnostics' for all options)
#' diagnostics(fit, parameter = "mu", type = "chains")
#'
#' # the RoBMA-PP can be fitted with addition of the 'model_type' argument
#' fit_PP <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study, model_type = "PP")
#'
#' # as well as the original version of RoBMA (with two weightfunctions)
#' fit_original <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study,
#'                       model_type = "2w")
#'
#' # or different prior distribution for the effect size (e.g., a half-normal distribution)
#' # (see 'vignette("CustomEnsembles")' for a detailed guide on specifying a custom model ensemble)
#' fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study,
#'              priors_effect = prior("normal", parameters = list(0, 1),
#'                                    truncation = list(0, Inf)))
#' }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return \code{RoBMA} returns an object of class 'RoBMA'.
#'
#' @seealso [summary.RoBMA()], [update.RoBMA()], [check_setup()]
#' @export
RoBMA <- function(
  # data specification
  d = NULL, r = NULL, logOR = NULL, OR = NULL, z = NULL, y = NULL,
  se = NULL, v = NULL, n = NULL, lCI = NULL, uCI = NULL, t = NULL, study_names = NULL, study_ids = NULL,
  data = NULL, weight = NULL,
  transformation   = if(is.null(y)) "fishers_z" else "none",
  prior_scale      = if(is.null(y)) "cohens_d"  else "none",
  effect_direction = "positive",

  # prior specification
  model_type = NULL, rescale_priors = 1,

  priors_effect              = set_default_priors("effect",        rescale = rescale_priors),
  priors_heterogeneity       = set_default_priors("heterogeneity", rescale = rescale_priors),
  priors_bias                = set_default_priors("bias",          rescale = rescale_priors),
  priors_effect_null         = set_default_priors("effect",        null = TRUE),
  priors_heterogeneity_null  = set_default_priors("heterogeneity", null = TRUE),
  priors_bias_null           = set_default_priors("bias",          null = TRUE),
  priors_hierarchical        = set_default_priors("hierarchical"),
  priors_hierarchical_null   = set_default_priors("hierarchical", null = TRUE),

  # MCMC fitting settings
  algorithm = "bridge", chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
  autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

  # additional settings
  save = "all", seed = NULL, silent = TRUE, ...){

  dots         <- .RoBMA_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()

  ### prepare & check the data
  if("data.RoBMA" %in% class(data)){
    object$data <- data
  }else{
    object$data <- combine_data(d = d, r = r, z = z, logOR = logOR, OR = OR, t = t, y = y, se = se, v = v, n = n, lCI = lCI, uCI = uCI,
                                study_names = study_names, study_ids = study_ids, weight = weight, data = data, transformation = transformation)
  }

  # switch between multivariate and weighted models
  if(attr(object$data, "weighted"))
    .weighted_warning()

  if(.is_multivariate(object))
    .multivariate_warning()


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
    algorithm        = algorithm,
    seed             = seed,
    save             = save,
    warnings         = NULL,
    errors           = NULL
  )


  ### prepare and check the settings
  object$priors   <- .check_and_list_priors(
    model_type = object$add_info[["model_type"]],
    priors_effect_null = priors_effect_null, priors_effect = priors_effect,
    priors_heterogeneity_null = priors_heterogeneity_null, priors_heterogeneity = priors_heterogeneity,
    priors_bias_null = priors_bias_null, priors_bias = priors_bias,
    priors_hierarchical_null = priors_hierarchical_null, priors_hierarchical = priors_hierarchical,
    prior_scale = object$add_info[["prior_scale"]])
  if(algorithm == "bridge"){
    object$models <- .make_models(object[["priors"]], .is_multivariate(object), .is_weighted(object))
  }else if(algorithm == "ss"){
    object$model  <- .make_model_ss(object[["priors"]], .is_multivariate(object), .is_weighted(object))
  }
  object$add_info$warnings <- c(object$add_info[["warnings"]], .check_effect_direction(object))

  if(dots[["do_not_fit"]]){
    return(object)
  }


  ### fit the models and compute marginal likelihoods
  if(object$add_info[["algorithm"]] == "bridge"){

    if(!object$fit_control[["parallel"]]){

      # sequential model fitting using JAGS & bridge sampling
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

      # parallel model fitting using JAGS & bridge sampling
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
        object <- .balance_component_probability(object)
      }

      ### compute the model-space results
      object$models        <- BayesTools::models_inference(object[["models"]])
      object$RoBMA         <- .ensemble_inference(object)
      object$coefficients  <- .compute_coeficients(object[["RoBMA"]])
    }

  }else if(object$add_info[["algorithm"]] == "ss"){

    # model fitting using JAGS with spike and slab priors
    object$model         <- .fit_RoBMA_model.ss(object)
    object$RoBMA         <- .as_ensemble_inference(object)
    object$coefficients  <- .compute_coeficients(object[["RoBMA"]])

  }


  ### collect and print errors and warnings
  object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   .get_model_errors(object))
  object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], .get_model_warnings(object))
  .print_errors_and_warnings(object)


  ### remove model posteriors if asked to
  if(save == "min" && object$add_info[["algorithm"]] == "bridge"){
    object <- .remove_model_posteriors(object)
    object <- .remove_model_margliks(object)
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
#' @param object a fitted RoBMA object
#' @param prior_effect prior distribution for the effect size (\code{mu})
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_heterogeneity prior distribution for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_bias prior distribution for the publication bias adjustment
#' component that will be treated as belonging to the alternative hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_effect_null prior distribution for the effect size (\code{mu})
#' parameter that will be treated as belonging to the null hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_heterogeneity_null prior distribution for the heterogeneity \code{tau}
#' parameter that will be treated as belonging to the null hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_bias_null prior distribution for the publication bias adjustment
#' component that will be treated as belonging to the null hypothesis.
#' Defaults to \code{NULL}.
#' @param prior_hierarchical prior distribution for the correlation of random effects
#' (\code{rho}) parameter that will be treated as belonging to the alternative hypothesis. This setting allows
#' users to fit a hierarchical (three-level) meta-analysis when \code{study_ids} are supplied.
#' Note that this is an experimental feature and see News for more details. Defaults to a beta distribution
#' \code{prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))}.
#' @param prior_hierarchical_null prior distribution for the correlation of random effects
#' (\code{rho}) parameter that will be treated as belonging to the null hypothesis. Defaults to \code{NULL}.
#' @param prior_weights either a single value specifying prior model weight
#' of a newly specified model using priors argument, or a vector of the
#' same length as already fitted models to update their prior weights.
#' @param refit_failed whether failed models should be refitted. Relevant only
#' if new priors or \code{prior_weights} are not supplied. Defaults to \code{TRUE}.
#' @param extend_all extend sampling in all fitted models based on \code{"sample_extend"}
#' argument in [set_autofit_control()] function. Defaults to \code{FALSE}.
#' @inheritParams RoBMA
#' @param ... additional arguments.
#'
#' @details See [RoBMA()] for more details.
#'
#' @examples \dontrun{
#' # using the example data from Bem 2011 and fitting the default (RoBMA-PSMA) model
#' fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study)
#'
#' # the update function allows us to change the prior model weights of each model
#' fit1 <- update(fit, prior_weights = c(0, rep(1, 35)))
#'
#' # add an additional model with different priors specification
#' # (see '?prior' for more information)
#' fit2 <- update(fit,
#'                priors_effect_null = prior("point", parameters = list(location = 0)),
#'                priors_heterogeneity = prior("normal",
#'                                   parameters = list(mean = 0, sd = 1),
#'                                   truncation = list(lower = 0, upper = Inf)),
#'                priors_bias = prior_weightfunction("one-sided",
#'                                     parameters = list(cuts = c(.05, .10, .20),
#'                                                       alpha = c(1, 1, 1, 1))))
#'
#' # update the models with an increased number of sample iterations
#' fit3 <- update(fit, autofit_control = set_autofit_control(sample_extend = 1000), extend_all = TRUE)
#' }
#'
#'
#' @return \code{RoBMA} returns an object of class 'RoBMA'.
#'
#' @seealso [RoBMA()], [summary.RoBMA()], [prior()], [check_setup()]
#' @export
update.RoBMA <- function(object, refit_failed = TRUE, extend_all = FALSE,
                         prior_effect = NULL,      prior_heterogeneity = NULL,      prior_bias = NULL,      prior_hierarchical = NULL, prior_weights = NULL,
                         prior_effect_null = NULL, prior_heterogeneity_null = NULL, prior_bias_null = NULL, prior_hierarchical_null = NULL,
                         study_names = NULL,
                         chains = NULL, adapt = NULL, burnin = NULL, sample = NULL, thin = NULL, autofit = NULL, parallel = NULL,
                         autofit_control = NULL, convergence_checks = NULL,
                         save = "all", seed = NULL, silent = TRUE, ...){

  BayesTools::check_bool(refit_failed, "refit_failed")
  BayesTools::check_bool(extend_all, "extend_all")

  dots         <- .RoBMA_collect_dots(...)
  output_scale <- NULL   # TODO: add ability to change the output scale

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

    if(is.RoBMA.reg(object))
      stop("Adding a new model to the ensemble is not possible with RoBMA.reg models.")

    what_to_do <- "fit_new_model"
    message("Fitting a new model with specified priors.")
    new_priors <- .check_and_list_priors(
      model_type = NULL,
      priors_effect_null        = prior_effect_null,        priors_effect        = prior_effect,
      priors_heterogeneity_null = prior_heterogeneity_null, priors_heterogeneity = prior_heterogeneity,
      priors_bias_null          = prior_bias_null,          priors_bias          = prior_bias,
      priors_hierarchical_null  = prior_hierarchical_null,  priors_hierarchical  = prior_hierarchical,
      prior_scale = object$add_info[["prior_scale"]])

    object$models[length(object$models) + 1]  <- list(.make_models(new_priors, .is_multivariate(object), .is_weighted(object))[[1]])

    if(!is.null(prior_weights)){
      object$models[[length(object$models)]]$prior_weights     <- prior_weights
      object$models[[length(object$models)]]$prior_weights_set <- prior_weights
    }


  }else if(!is.null(prior_weights)){

    what_to_do <- "update_prior_weights"
    message("Updating prior odds for the fitted models.")
    if(length(prior_weights) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.RoBMA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_weights     <- prior_weights[i]
      object$models[[i]]$prior_weights_set <- prior_weights[i]
    }

  }else if(!is.null(output_scale)){

    stop("this functionality is not currently implemented")
    what_to_do <- "transform_estimates"

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

    object[["models"]][[length(object$models)]] <- .fit_RoBMA_model(object, length(object$models))

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
        object$models[[i]] <- .fit_RoBMA_model(object, i, extend = TRUE)
        if(dots[["is_JASP"]]){
          .JASP_progress_bar_tick()
        }
      }

    }else{

      cl <- parallel::makePSOCKcluster(floor(RoBMA.get_option("max_cores") / object$fit_control[["chains"]]))
      parallel::clusterEvalQ(cl, {library("RoBMA")})
      parallel::clusterExport(cl, "object", envir = environment())
      object$models[models_to_update] <- parallel::parLapplyLB(cl, models_to_update, .fit_RoBMA_model, object = object, extend = TRUE)
      parallel::stopCluster(cl)

    }

  }else if(what_to_do == "update_convergence_checks"){

    # propagate settings changes to the individual models
    for(i in seq_along(object$models)){
      object$models[[i]] <- .update_model_checks(object$models[[i]], object[["convergence_checks"]])
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
