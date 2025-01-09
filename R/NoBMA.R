#' @title Estimate a Bayesian Model-Averaged Meta-Analysis
#'
#' @description \code{NoBMA} is a wrapper around [RoBMA()] that can
#' be used to estimate a publication bias unadjusted Bayesian
#' model-averaged meta-analysis. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @inheritParams RoBMA
#' @inheritParams combine_data
#'
#' @details See [RoBMA()] for more details.
#'
#' Note that these default prior distributions are relatively wide and more informed
#' prior distributions for testing for the presence of moderation should be considered.
#'
#' @return \code{NoBMA} returns an object of class 'RoBMA'.
#'
#' @seealso [RoBMA()], [summary.RoBMA()], [update.RoBMA()], [check_setup()]
#' @export
NoBMA <- function(
    # data specification
  d = NULL, r = NULL, logOR = NULL, OR = NULL, z = NULL, y = NULL,
  se = NULL, v = NULL, n = NULL, lCI = NULL, uCI = NULL, t = NULL, study_names = NULL, study_ids = NULL,
  data = NULL, weight = NULL,
  transformation   = if(is.null(y)) "fishers_z" else "none",
  prior_scale      = if(is.null(y)) "cohens_d"  else "none",

  # prior specification
  model_type   = NULL, rescale_priors = 1,
  priors_effect              = set_default_priors("effect",        rescale = rescale_priors),
  priors_heterogeneity       = set_default_priors("heterogeneity", rescale = rescale_priors),
  priors_effect_null         = set_default_priors("effect", null = TRUE),
  priors_heterogeneity_null  = set_default_priors("heterogeneity", null = TRUE),
  priors_hierarchical        = set_default_priors("hierarchical"),
  priors_hierarchical_null   = set_default_priors("hierarchical", null = TRUE),

  # MCMC fitting settings
  algorithm = "bridge", chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
  autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

  # additional settings
  save = "all", seed = NULL, silent = TRUE, ...){


  object <- RoBMA(
    # data specification
    d = d, r = r, logOR = logOR, OR = OR, z = z, y = y,
    se = se, v = v, n = n, lCI = lCI, uCI = uCI, t = t, study_names = study_names, study_ids = study_ids,
    data = data,
    transformation   = transformation,
    prior_scale      = prior_scale,
    effect_direction = "positive", # THIS IS PRESET

    # prior specification
    model_type            = model_type,
    priors_effect         = priors_effect,
    priors_heterogeneity  = priors_heterogeneity,
    priors_bias           = NULL, # THIS IS PRESET
    priors_effect_null         = priors_effect_null,
    priors_heterogeneity_null  = priors_heterogeneity_null,
    priors_bias_null           = prior_none(),  # THIS IS PRESET
    priors_hierarchical        = priors_hierarchical,
    priors_hierarchical_null   = priors_hierarchical_null,

    # MCMC fitting settings
    algorithm = algorithm, chains = chains, sample = sample, burnin = burnin, adapt = adapt, thin = thin, parallel = parallel,
    autofit = autofit, autofit_control = autofit_control, convergence_checks = convergence_checks,

    # additional settings
    save = save, seed = seed, silent = silent, ...)

  class(object) <- c("NoBMA", class(object))

  return(object)
}

#' @title Estimate a Bayesian Model-Averaged Meta-Regression
#'
#' @description \code{NoBMA.reg} is a wrapper around [RoBMA.reg()] that can
#' be used to estimate a publication bias unadjusted Bayesian
#' model-averaged meta-regression. The interface allows a complete customization of
#' the ensemble with different prior (or list of prior) distributions
#' for each component.
#'
#' @inheritParams RoBMA.reg
#' @inheritParams RoBMA
#' @inheritParams combine_data
#'
#' @details See [RoBMA.reg()] for more details.
#'
#' Note that these default prior distributions are relatively wide and more informed
#' prior distributions for testing for the presence of moderation should be considered.
#'
#' @return \code{NoBMA.reg} returns an object of class 'RoBMA'.
#'
#' @seealso [RoBMA()], [RoBMA.reg()], [summary.RoBMA()], [update.RoBMA()], [check_setup()]
#' @export
NoBMA.reg <- function(
    formula, data, test_predictors = TRUE, study_names = NULL, study_ids = NULL,
    transformation     = if(any(colnames(data) != "y")) "fishers_z" else "none",
    prior_scale        = if(any(colnames(data) != "y")) "cohens_d"  else "none",
    standardize_predictors = TRUE,

    # prior specification
    priors = NULL, model_type = NULL, rescale_priors = 1,
    priors_effect              = set_default_priors("effect",        rescale = rescale_priors),
    priors_heterogeneity       = set_default_priors("heterogeneity", rescale = rescale_priors),
    priors_effect_null         = set_default_priors("effect", null = TRUE),
    priors_heterogeneity_null  = set_default_priors("heterogeneity", null = TRUE),
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

  object <- RoBMA.reg(
    formula = formula, data = data, test_predictors = test_predictors, study_names = study_names, study_ids = study_ids,
    transformation     = transformation,
    prior_scale        = prior_scale,
    standardize_predictors = standardize_predictors,
    effect_direction       = "positive",  # THIS IS PRESET

    # prior specification
    priors       = priors,
    model_type   = model_type,

    priors_effect         = priors_effect,
    priors_heterogeneity  = priors_heterogeneity,
    priors_bias           = NULL,  # THIS IS PRESET
    priors_effect_null         = priors_effect_null,
    priors_heterogeneity_null  = priors_heterogeneity_null,
    priors_bias_null           = prior_none(),  # THIS IS PRESET
    priors_hierarchical        = priors_hierarchical,
    priors_hierarchical_null   = priors_hierarchical_null,

    prior_covariates       = prior_covariates,
    prior_covariates_null  = prior_covariates_null,
    prior_factors          = prior_factors,
    prior_factors_null     = prior_factors_null,

    # MCMC fitting settings
    algorithm = algorithm, chains = chains, sample = sample, burnin = burnin, adapt = adapt, thin = thin, parallel = parallel,
    autofit = autofit, autofit_control = autofit_control, convergence_checks = convergence_checks,

    # additional settings
    save = save, seed = seed, silent = silent, ...)

   class(object) <- c("NoBMA.reg", class(object))

   return(object)
}
