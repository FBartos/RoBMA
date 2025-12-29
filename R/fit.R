### RoBMA 4.0.0

# Fitting functions -----
.create_fit_priors <- function(data, priors) {
  return(priors[["outcome"]])
}
.create_fit_data   <- function(data, priors) {

  # always include yi and sei
  fit_data <- list(
    yi  = data[["outcome"]][["yi"]],
    sei = data[["outcome"]][["sei"]]
  )

  # flip effect size direction (needed for selection models, PET, and PEESE models)
  # (done for everything for consistency)
  effect_direction  <- "positive" # TODO: forward from models with publication bias
  if (effect_direction == "negative") {
    fit_data[["yi"]] <- -1 * fit_data[["yi"]]
  }

  # add study_ids for 3lvl models
  if (.is_data_multilevel(data)) {
    fit_data[["study_ids"]] <- data[["outcome"]][["study_ids"]]
  }

  # add weights for weighted models
  if (.is_data_weights(data)) {
    fit_data[["weights"]] <- data[["outcome"]][["study_ids"]]
  }

  # add number of estimates
  fit_data[["K"]] <- length(data[["outcome"]][["yi"]])

  return(fit_data)
}
.create_fit_formula_list       <- function(data, parameter) {

  formula <- attr(data[[parameter]], "formula")

  # for scale regression - modify formula for exponential parameterization
  # (required for the exponential parameterization trick in BayesTools::JAGS_fit)
  if (parameter == "scale") {

    # check whether the intercept was removed from the formula
    # if so, warn the user and replace it with log(1) (intercepts cannot be omitted from scale models)
    if (attr(terms(formula), "intercept") == 0) {
      warning("Intercept cannot be omitted from scale models (the regression estimates a multiplicative constant for the intercept). The intercept removal term has been ignored.", call. = FALSE)
      formula <- .replace_intercept_with_log1(formula)
    } else {
      # add log(1) term (adds log of the intercept parameter prior)
      formula <- update(formula, . ~ . - 1 + log(1))
    }
  }

  return(formula)
}
.create_fit_formula_data_list  <- function(data, parameter) {
  return(data[[parameter]])
}
.create_fit_formula_prior_list <- function(priors, parameter) {
  return(priors[[parameter]])
}
.create_model_syntax           <- function(data, priors) {

  ### extract structural information about the model
  is_mods           <- .is_data_mods(data)
  is_scale          <- .is_data_scale(data)
  is_multilevel     <- .is_data_multilevel(data)
  is_weights        <- .is_data_weights(data)
  is_PET            <- .is_priors_PET(priors)
  is_PEESE          <- .is_priors_PEESE(priors)
  is_weightfunction <- .is_priors_weightfunction(priors)
  effect_direction  <- "positive" # TODO: forward from models with publication bias

  ### create the model syntax
  model_syntax <- "model{\n"

  ### the main model parameters are created automatically via BayesTools::JAGS_fit
  # - mu  (!is_mods)  / mu[i]  (is_mods)
  # - tau (!is_scale) / tau[i] (is_scale)
  # - rho (is_multilevel)
  # for publication bias
  # - PET            (is_PET)
  # - PEESE          (is_PEESE)
  # - weightfunction (is_weightfunction)

  ### prepare multilevel heterogeneity parameter (non-regression = defined outside of the for loop)
  # for multilevel: specify heterogeneity allocation
  # tau_within  - within-study variance (estimate-level variance)
  # tau_between - between-study variance (study-level variance)
  if (is_multilevel & !is_scale) {
    model_syntax <- paste0(model_syntax, "tau_within  = tau * sqrt(rho)\n")
    model_syntax <- paste0(model_syntax, "tau_between = tau * sqrt(1-rho)\n")
    tau_estimate <- "tau_within"
  }

  ### enter the main block
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")

  ### prepare effect size parameter
  # flip effect size direction (needed for selection models, PET, and PEESE models)
  # (done for everything for consistency, data are flipped within `.create_fit_data()`)
  if (is_mods) {
    mu_estimate <- ifelse(effect_direction == "negative", paste0("- mu[i]"), "mu[i]")
  } else {
    mu_estimate <- ifelse(effect_direction == "negative", "- mu", "mu")
  }
  # add study-level effects
  if (is_multilevel) {
    mu_estimate <- paste0(mu_estimate, ifelse(effect_direction == "negative", " - ", " + "),  "gamma[study_ids[i]] * tau_between")
  }
  # add PET/PEESE
  if (is_PET) {
    mu_estimate <- paste0(mu_estimate, " + PET * sei[i]")
  }
  if (is_PEESE) {
    mu_estimate <- paste0(mu_estimate, " + PEESE * pow(sei[i],2)")
  }


  ### prepare heterogeneity parameter (regression)
  # in the case of a regression: the scale model is specified on multiplicative scale:
  # tau = intercept * exp(sum(beta_i * x_i))
  # the model uses a trick to separately pass the tau_intercept and multiply it with
  # tau_exp scale regression with zero intercept
  if (is_scale) {
    model_syntax <- paste0(model_syntax, "  tau[i] = tau_intercept * exp(tau_exp[i])\n")
  }
  # for multilevel: specify heterogeneity allocation & dispatch estimate-specific/common parameter
  # tau_within  - within-study variance (estimate-level variance)
  # tau_between - between-study variance (study-level variance)
  # for rest: dispatch estimate-specific/common parameter
  if (is_multilevel && is_scale) {
    model_syntax <- paste0(model_syntax, "  tau_within[i]  = tau[i] * sqrt(rho)\n")
    model_syntax <- paste0(model_syntax, "  tau_between[i] = tau[i] * sqrt(1-rho)\n")
    tau_estimate <- "tau_within[i]"
  } else {
    if (is_scale) {
      tau_estimate <- "tau[i]"
    } else {
      tau_estimate <- "tau"
    }
  }
  # compute the total variance of the estimate
  tau2_estimate <- paste0("( pow(sei[i],2) + pow(",tau_estimate,",2) )")


  ### specify model likelihood
  # separate likelihoods for selection models ans simple models
  # subversion for weighted/unweighted likelihoods
  if (is_weightfunction) {
    if (is_weights) {
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
      model_syntax <- paste0(model_syntax, "  yi[i] ~ dwnorm(", mu_estimate, ",", "1/", tau2_estimate, ", weight[i])\n")
    } else {
      model_syntax <- paste0(model_syntax, "  yi[i] ~ dnorm(",  mu_estimate, ",", "1/", tau2_estimate, ")\n")
    }
  }

  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}

.fit <- function(object, extend = FALSE) {

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
  if (!is.null(priors[["mods"]])) {
    fit_formula_list[["mu"]]       <- .create_fit_formula_list(data = data, parameter = "mods")
    fit_formula_data_list[["mu"]]  <- .create_fit_formula_data_list(data = data, parameter = "mods")
    fit_formula_prior_list[["mu"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "mods")
    fit_formula_scale_list[["mu"]] <- object[["standardize_continuous_predictors"]]
  }

  ### add heterogeneity regressions
  if (!is.null(priors[["scale"]])) {
    fit_formula_list[["tau_exp"]]       <- .create_fit_formula_list(data = data, parameter = "scale")
    fit_formula_data_list[["tau_exp"]]  <- .create_fit_formula_data_list(data = data, parameter = "scale")
    fit_formula_prior_list[["tau_exp"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "scale")
    fit_formula_scale_list[["tau_exp"]] <- object[["standardize_continuous_predictors"]]
  }

  ### generate the model syntax
  model_syntax <- .create_model_syntax(data = data, priors = priors)

  ### fit the model
  if (!extend || length(model[["fit"]]) == 0) {

    fit <- BayesTools::JAGS_fit(
      model_syntax          = model_syntax,
      data                  = fit_data,
      prior_list            = fit_priors,
      formula_list          = if (length(fit_formula_list)       > 0) fit_formula_list,
      formula_data_list     = if (length(fit_formula_data_list)  > 0) fit_formula_data_list,
      formula_prior_list    = if (length(fit_formula_prior_list) > 0) fit_formula_prior_list,
      formula_scale         = if (length(fit_formula_scale_list) > 0) fit_formula_scale_list,
      chains                = fit_control[["chains"]],
      adapt                 = fit_control[["adapt"]],
      burnin                = fit_control[["burnin"]],
      sample                = fit_control[["sample"]],
      thin                  = fit_control[["thin"]],
      autofit               = fit_control[["autofit"]],
      autofit_control       = autofit_control,
      parallel              = fit_control[["parallel"]],
      cores                 = fit_control[["cores"]],
      silent                = fit_control[["silent"]],
      seed                  = fit_control[["seed"]],
      required_packages     = "RoBMA",
      is_JASP               = object[["is_JASP"]],
      is_JASP_prefix        = object[["is_JASP_prefix"]]
    )

  } else {

    fit <- BayesTools::JAGS_extend(
      fit                = object[["fit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]]
    )

  }


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


# Helper functions -----
.replace_intercept_with_log1 <- function(formula) {

  # converts formula with -1 or 0 (no intercept) back to formula with log(1)
  # by removing the -1 or 0 term and adding log(1)
  formula_str <- paste(deparse(formula), collapse = " ")

  # Remove various forms of -1, + -1, 0, or + 0
  formula_str <- gsub("\\s*\\+\\s*\\-\\s*1", "", formula_str)    # handles "+ - 1"
  formula_str <- gsub("\\s*\\-\\s*1\\s*", "", formula_str)        # handles "- 1"
  formula_str <- gsub("\\s*\\+\\s*0\\s*", "", formula_str)        # handles "+ 0"

  # Handle 0 at the start (e.g., "~ 0 + x")
  formula_str <- gsub("~\\s*0\\s*\\+\\s*", "~ ", formula_str)

  # Handle 0 alone (e.g., "~ 0")
  formula_str <- gsub("~\\s*0\\s*$", "~", formula_str)

  # Handle case where formula becomes empty (just "~")
  if (grepl("^\\s*~\\s*$", formula_str)) {
    formula_str <- "~ log(1)"
  } else {
    # Clean up any trailing + before adding log(1)
    formula_str <- gsub("\\+\\s*$", "", formula_str)
    # Add log(1) to the formula
    formula_str <- paste0(formula_str, " + log(1)")
  }

  # Clean up any double spaces
  formula_str <- gsub("\\s{2,}", " ", formula_str)
  formula_str <- trimws(formula_str)

  return(stats::as.formula(formula_str))
}

.is_priors_PET            <- function(priors) {

  if (is.null(priors[["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["bias"]]))
    return(any(sapply(priors[["bias"]], is.prior.PET)))

  return(is.prior.PET(priors[["bias"]]))
}
.is_priors_PEESE          <- function(priors) {

  if (is.null(priors[["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["bias"]]))
    return(any(sapply(priors[["bias"]], is.prior.PEESE)))

  return(is.prior.PEESE(priors[["bias"]]))
}
.is_priors_weightfunction <- function(priors) {

  if (is.null(priors[["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["bias"]]))
    return(any(sapply(priors[["bias"]], is.prior.weightfunction)))

  return(is.prior.weightfunction(priors[["bias"]]))
}
.is_data_multilevel       <- function(data) {
  return(attr(data, "study_ids"))
}
.is_data_mods             <- function(data) {
  return(attr(data, "mods"))
}
.is_data_scale            <- function(data) {
  return(attr(data, "scale"))
}
.is_data_weights          <- function(data) {
  return(attr(data, "weights"))
}

.is_PET            <- function(object) .is_priors_PET(object[["priors"]])
.is_PEESE          <- function(object) .is_priors_PEESE(object[["priors"]])
.is_weightfunction <- function(object) .is_priors_weightfunction(object[["priors"]])
.is_multilevel     <- function(object) .is_data_multilevel(object[["data"]])
.is_mods           <- function(object) .is_data_mods(object[["data"]])
.is_scale          <- function(object) .is_data_scale(object[["data"]])
.is_weights        <- function(object) .is_data_weights(object[["data"]])
