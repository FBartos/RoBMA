### RoBMA 4.0.0
# in progress -- important pieces missing
### compute marginal likelihood
marglik <- function(x) {

  # TODO: implement

}

.log_posterior <- function(
    parameters, data,
    is_mods, is_scale, is_multilevel, is_weights,
    is_PET, is_PEESE, is_weightfunction, effect_direction,
    outcome_type) {

  ### follow the same structure as the generated JAGS code
  # TODO: base on brma.evaluate.R and brma.pdf.R

  ### prepare heterogeneity parameter
  if (is_scale) {
    tau <- exp(parameters[["log_tau"]])
  } else {
    tau <- rep(parameters[["tau"]], data[["K"]])
  }
  if (is_multilevel) {
    tau_within   <- parameters[["tau"]] * sqrt(parameters[["rho"]])
    tau_between  <- parameters[["tau"]] * sqrt(1-parameters[["rho"]])
    tau_estimate <- tau_within^2 + data[["se"]]^2
  } else {
    tau_estimate <- tau^2 + data[["se"]]^2
  }

  ### prepare effect size parameter
  if (is_mods) {
    mu_estimate <- ifelse(effect_direction == "negative", -1, 1) * parameters[["mu"]]
  } else {
    mu_estimate <- ifelse(effect_direction == "negative", -1, 1) * rep(parameters[["mu"]], data[["K"]])
  }
  # add study-level effects
  if (is_multilevel) {
    mu_estimate <- mu_estimate + ifelse(effect_direction == "negative", -1, 1) * parameters[["gamma"]][data[["study_ids"]]] * tau_between
  }
  # add PET/PEESE
  if (is_PET) {
    mu_estimate <- mu_estimate + parameters[["PET"]] * data[["sei"]]
  }
  if (is_PEESE) {
    mu_estimate <- mu_estimate + parameters[["PEESE"]] * data[["sei"]]^2
  }

  ### specify model likelihood
  if (outcome_type == "norm") {

    # selection and normal models for norm data outcome
    if (is_weightfunction) {
      if (is_weights) {
        lik <- sum()
        model_syntax <- paste0(
          model_syntax, "  yi[i] ~ dwwnorm_mix(", mu_estimate, ",", "sqrt", tau2_estimate,
          ", crit_y[,i], omega, crit_y_mapping[,bias_indicator], crit_y_mapping_max[bias_indicator], weight[i])\n")
      } else {
        model_syntax <- paste0(
          model_syntax, "  yi[i] ~ dwnorm_mix(", mu_estimate, ",", "sqrt", tau2_estimate,
          ", crit_y[,i], omega, crit_y_mapping[,bias_indicator], crit_y_mapping_max[bias_indicator])\n")
      }
    } else {
      if (is_weights) {
        log_lik <- sum(stats::dnorm(data[["yi"]], mean = mu_estimate, sd = sqrt(tau_estimate), log = TRUE) * data[["weights"]])
      } else {
        log_lik <- sum(stats::dnorm(data[["yi"]], mean = mu_estimate, sd = sqrt(tau_estimate), log = TRUE))
      }
    }

  } else if (outcome_type %in% c("bin", "pois")) {

    # estimate-level variance is not marginalized in binomial-normal/poisson-normal models
    if (is_multilevel) {
      mu_estimate <- paste0(mu_estimate, " + theta[i] * tau_within")
    } else {
      mu_estimate <- paste0(mu_estimate, " + theta[i] * tau")
    }

    # specify appropriate link for each type of outcome
    if (outcome_type == "bin") {
      # transform the parameters to the probability scale
      model_syntax <- paste0(model_syntax, "  logit(p1[i]) = logit(pi[i]) + 0.5 * (", mu_estimate, ")\n")
      model_syntax <- paste0(model_syntax, "  logit(p2[i]) = logit(pi[i]) - 0.5 * (", mu_estimate, ")\n")
    } else if (outcome_type == "pois") {
      # transform the parameters to the probability scale
      model_syntax <- paste0(model_syntax, "  log(r1[i]) = phi[i] + 0.5 * (", mu_estimate, ") + log(t1i[i])\n")
      model_syntax <- paste0(model_syntax, "  log(r2[i]) = phi[i] - 0.5 * (", mu_estimate, ") + log(t2i[i])\n")
    }

    # specify appropriate likelihoods
    if (outcome_type == "bin") {
      # the observed data
      if (is_weights){
        stop("TODO: implement")
        model_syntax <- paste0(model_syntax, "  ai[i] ~ dwbinom(p1[i], n1i[i], weight[i])\n")
        model_syntax <- paste0(model_syntax, "  ci[i] ~ dwbinom(p2[i], n2i[i], weight[i])\n")
      } else {
        model_syntax <- paste0(model_syntax, "  ai[i] ~ dbinom(p1[i], n1i[i])\n")
        model_syntax <- paste0(model_syntax, "  ci[i] ~ dbinom(p2[i], n2i[i])\n")
      }
    } else if (outcome_type == "pois") {
      # the observed data
      if(is_weights){
        stop("TODO: implement")
        model_syntax <- paste0(model_syntax, "  x1i[i] ~ dwpois(r1[i], weight[i])\n")
        model_syntax <- paste0(model_syntax, "  x2i[i] ~ dwpois(r2[i], weight[i])\n")
      }else{
        model_syntax <- paste0(model_syntax, "  x1i[i] ~ dpois(r1[i])\n")
        model_syntax <- paste0(model_syntax, "  x2i[i] ~ dpois(r2[i])\n")
      }
    }

  }

  return(log_lik)
}

.marglik <- function(object) {

  fit                <- object$fit
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  data               <- object[["data"]]
  priors             <- object[["priors"]]

  errors   <- NULL
  warnings <- NULL

  ### create arguments to be passed to BayesTools::JAGS_fit
  fit_formula_list        <- list()
  fit_formula_data_list   <- list()
  fit_formula_prior_list  <- list()
  fit_formula_scale_list  <- list()

  ### create model base
  fit_priors <- .create_fit_priors(data = data, priors = priors)
  fit_data   <- .create_fit_data(data = data, priors = priors)

  ### add effect regressions
  if (.is_data_mods(data)) {
    fit_formula_list[["mu"]]       <- .create_fit_formula_list(data = data, parameter = "mods")
    fit_formula_data_list[["mu"]]  <- .create_fit_formula_data_list(data = data, parameter = "mods")
    fit_formula_prior_list[["mu"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "mods")
    fit_formula_scale_list[["mu"]] <- object[["standardize_continuous_predictors"]]
  }

  ### add heterogeneity regressions
  if (.is_data_scale(data)) {
    fit_formula_list[["log_tau"]]       <- .create_fit_formula_list(data = data, parameter = "scale")
    fit_formula_data_list[["log_tau"]]  <- .create_fit_formula_data_list(data = data, parameter = "scale")
    fit_formula_prior_list[["log_tau"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "scale")
    fit_formula_scale_list[["log_tau"]] <- object[["standardize_continuous_predictors"]]
  }

  ### generate the model syntax
  log_posterior <- .create_model_log_posterior(data = data, priors = priors)

  ### compute marginal likelihood
  fit <- BayesTools::JAGS_bridgesampling(
    fit                   = fit,
    log_posterior         = .log_posterior,
    data                  = fit_data,
    prior_list            = fit_priors,
    formula_list          = if (length(fit_formula_list)       > 0) fit_formula_list,
    formula_data_list     = if (length(fit_formula_data_list)  > 0) fit_formula_data_list,
    formula_prior_list    = if (length(fit_formula_prior_list) > 0) fit_formula_prior_list,
    formula_scale         = if (length(fit_formula_scale_list) > 0) fit_formula_scale_list,
    is_mods               = .is_data_mods(data),
    is_scale              = .is_data_scale(data),
    is_multilevel         = .is_data_multilevel(data),
    is_weights            = .is_data_weights(data),
    is_PET                = .is_priors_PET(priors),
    is_PEESE              = .is_priors_PEESE(priors),
    is_weightfunction     = .is_priors_weightfunction(priors),
    effect_direction      = "positive", # TODO: forward from models with publication bias
    outcome_type          = .data_outcome_type(data)
  )


  # assess the model fit and deal with errors
  if (inherits(fit, "error")) {

    if(grepl("Unknown function", fit$message))
      stop("The RoBMA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.", call. = FALSE)

    fit                     <- list()
    attr(fit, "prior_list") <- fit_priors

    converged      <- FALSE
    has_posterior  <- FALSE
    errors         <- c(errors, fit$message)

  } else {

    has_posterior <- TRUE
    check_fit     <- BayesTools::JAGS_check_convergence(
      fit          = fit,
      prior_list   = attr(fit, "prior_list"),
      max_Rhat     = convergence_checks[["max_Rhat"]],
      min_ESS      = convergence_checks[["min_ESS"]],
      max_error    = convergence_checks[["max_error"]],
      max_SD_error = convergence_checks[["max_SD_error"]]
    )
    warnings  <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
    converged <- check_fit

  }

  # add results
  fit$errors        <- errors
  fit$warnings      <- warnings
  fit$converged     <- converged
  fit$has_posterior <- has_posterior

  return(fit)
}
