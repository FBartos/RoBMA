#' @title Estimate a Robust Bayesian Meta-Analysis Meta-Regression
#'
#' @description \code{RoBMA} is used to estimate a robust Bayesian
#' meta-regression. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @param formula a formula for the meta-regression model
#' @param data a data.frame containing the data for the meta-regression. Note that the
#' column names have to correspond to the effect sizes (\code{d}, \code{logOR}, \code{OR},
#' \code{r}, \code{z}), a measure of sampling variability (\code{se}, \code{v}, \code{n},
#' \code{lCI}, \code{uCI}, \code{t}), and the predictors.
#' See [combine_data()] for a complete list of reserved names and additional information
#' about specifying input data.
#' @param test_predictors vector of predictor names to test for the presence
#' of moderation (i.e., assigned both the null and alternative prior distributions).
#' Defaults to \code{TRUE}, all predictors are tested using the default
#' prior distributions (i.e., \code{prior_covariates},
#' \code{prior_covariates_null}, \code{prior_factors}, and
#' \code{prior_factors_null}). To only estimate
#' and adjust for the effect of predictors use \code{FALSE}. If
#' \code{priors} is specified, any settings in \code{test_predictors}
#' is overridden.
#' @param priors named list of prior distributions for each predictor
#' (with names corresponding to the predictors). It allows users to
#' specify both the null and alternative hypothesis prior distributions
#' for each predictor by assigning the corresponding element of the named
#' list with another named list (with \code{"null"} and
#' \code{"alt"}).
#' If only one prior is specified for a given parameter, it is
#' assumed to correspond to the alternative hypotheses and the default null
#' hypothesis is specified (i.e., \code{prior_covariates_null} or
#' \code{prior_factors_null}).
#' If a named list with only one named prior distribution is provided (either
#' \code{"null"} or \code{"alt"}), only this prior distribution is used and no
#' default distribution is filled in.
#' Parameters without specified prior distributions are assumed to be only adjusted
#' for using the default alternative hypothesis prior distributions (i.e.,
#' \code{prior_covariates} or \code{prior_factors}).
#' If \code{priors} is specified, \code{test_predictors} is ignored.
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
#' estimating the model. Defaults to \code{TRUE}. Continuous predictor standardization is important
#' for applying the default prior distributions for continuous predictors. Note that the resulting
#' output corresponds to standardized meta-regression coefficients.
#' @inheritParams RoBMA
#' @inheritParams combine_data
#'
#' @details The \href{../doc/MetaRegression.html}{\code{vignette("/MetaRegression", package = "RoBMA")}}
#' vignette describes how to use [RoBMA.reg()] function to fit Bayesian meta-regression ensembles. See
#' \insertCite{bartos2023robust;textual}{RoBMA} for more details about the methodology and
#' [RoBMA()] for more details about the function options. By default, the function standardizes
#' continuous predictors. As such, the output should be interpreted as standardized meta-regression
#' coefficients.
#'
#' The RoBMA.reg function first generates models from a combination of the
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
#' Estimated marginal means can be computed by [marginal_summary()] function and
#' visualized by the [marginal_plot()] function.
#'
#' @examples \dontrun{
#' # using the example data from Andrews et al. (2021) and reproducing the example from
#' # Bartos et al. (2024) with measure and age covariate.
#'
#'  # note the the Andrews2021 data.frame columns identify the effect size "r" and
#'  # the standard error "se" of the effect size that are used to estimate the model
#'  fit_RoBMA <- RoBMA.reg(~ measure + age, data = Andrews2021, parallel = TRUE, seed = 1)
#'
#'  # summarize the results
#'  summary(fit_RoBMA, output_scale = "r")
#'
#'  # compute effect size estimates for each group
#'  marginal_summary(fit_RoBMA, output_scale = "r")
#'
#'  # visualize the effect size estimates for each group
#'  marginal_plot(fit_RoBMA, parameter = "measure", output_scale = "r", lwd = 2)
#' }
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
    formula, data, test_predictors = TRUE, study_names = NULL, study_ids = NULL,
    transformation     = if(any(colnames(data) != "y")) "fishers_z" else "none",
    prior_scale        = if(any(colnames(data) != "y")) "cohens_d"  else "none",
    standardize_predictors = TRUE,
    effect_direction       = "positive",

    # prior specification
    priors = NULL, model_type = NULL, rescale_priors = 1,

    priors_effect              = set_default_priors("effect",        rescale = rescale_priors),
    priors_heterogeneity       = set_default_priors("heterogeneity", rescale = rescale_priors),
    priors_bias                = set_default_priors("bias",          rescale = rescale_priors),
    priors_effect_null         = set_default_priors("effect",        null = TRUE),
    priors_heterogeneity_null  = set_default_priors("heterogeneity", null = TRUE),
    priors_bias_null           = set_default_priors("bias",          null = TRUE),
    priors_hierarchical        = set_default_priors("hierarchical"),
    priors_hierarchical_null   = set_default_priors("hierarchical", null = TRUE),

    prior_covariates       = set_default_priors("covariates", rescale = rescale_priors),
    prior_covariates_null  = set_default_priors("covariates", null = TRUE),
    prior_factors          = set_default_priors("factors", rescale = rescale_priors),
    prior_factors_null     = set_default_priors("factors", null = TRUE),

    # MCMC fitting settings
    algorithm = "bridge", chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
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

  ### additional information
  object$add_info <- .check_and_list_add_info(
    model_type             = model_type,
    predictors             = attr(object[["priors"]], "terms"),
    predictors_test        = attr(object[["priors"]], "terms_test"),
    prior_scale            = .transformation_var(prior_scale),
    output_scale           = .transformation_var(prior_scale),
    effect_measure         = attr(object$data[["outcome"]], "effect_measure"),
    effect_direction       = effect_direction,
    algorithm              = algorithm,
    standardize_predictors = standardize_predictors,
    seed                   = seed,
    save                   = save,
    warnings               = NULL,
    errors                 = NULL
  )

  ### make models
  if(algorithm == "bridge"){
    object$models <- .make_models.reg(object[["priors"]], .is_multivariate(object), .is_weighted(object))
  }else if(algorithm == "ss"){
    object$model  <- .make_model_ss.reg(object[["priors"]], .is_multivariate(object), .get_K(object), .is_weighted(object))
  }


  # the check requires the 'add_info' object already created
  object$add_info[["warnings"]] <- c(.check_effect_direction(object), .check_predictors_scaling(object))


  if(dots[["do_not_fit"]]){
    return(object)
  }


  ### fit the models and compute marginal likelihoods
  if(algorithm == "bridge"){

    # using the individual models & bridge sampling for model-averaging
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
    object$model         <- .fit_RoBMA_model_ss(object, dots)
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


  class(object) <- c("RoBMA", "RoBMA.reg")
  return(object)
}

