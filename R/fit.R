### RoBMA 4.0.0

# Fitting functions -----
.create_fit_priors <- function(data, priors) {

  # extract the common prior list
  prior_list <- priors[["outcome"]]

  # add levels to the baseline prior for outcomes
  if (.data_outcome_type(data) == "bin") {
    # encode number of levels for the baserate and random-effects prior
    attr(prior_list[["pi"]], "levels")    <- nrow(data[["outcome"]])
    attr(prior_list[["theta"]], "levels") <- nrow(data[["outcome"]])
  } else if (.data_outcome_type(data) == "pois") {
    # encode number of levels for the baserate and random-effects prior
    attr(prior_list[["phi"]], "levels")   <- nrow(data[["outcome"]])
    attr(prior_list[["theta"]], "levels") <- nrow(data[["outcome"]])
  }

  # add study-level indicators
  if (.is_data_multilevel(data)) {
    # encode number of levels for the random-effects prior
    attr(prior_list[["gamma"]], "levels") <- length(unique(data[["outcome"]][["study_ids"]]))
  }

  ### deal with non-prior mixture distributions (bPET, bPEESE, and bselmodel)
  # the prior name must match the output parameter for sample extraction in plots etc
  if (.is_priors_bias(priors)) {
    if (!BayesTools::is.prior.mixture(prior_list[["bias"]])) {
      if (BayesTools::is.prior.weightfunction(prior_list[["bias"]])) {
        names(prior_list)[names(prior_list) == "bias"] <- "omega"
      } else if (BayesTools::is.prior.PET(prior_list[["bias"]])) {
        names(prior_list)[names(prior_list) == "bias"] <- "PET"
      } else if (BayesTools::is.prior.PEESE(prior_list[["bias"]])) {
        names(prior_list)[names(prior_list) == "bias"] <- "PEESE"
      }
    }
  }


  return(prior_list)
}
.create_fit_data   <- function(data, priors) {

  ### add outcome specific data
  if (.data_outcome_type(data) == "norm") {
    # always include yi and sei
    fit_data <- list(
      yi  = data[["outcome"]][["yi"]],
      sei = data[["outcome"]][["sei"]]
    )

    # flip effect size direction (needed for selection models, PET, and PEESE models)
    # (done for everything for consistency)
    effect_direction  <- .data_effect_direction(data)
    if (effect_direction == "negative") {
      fit_data[["yi"]] <- -1 * fit_data[["yi"]]
    }

    # add publication bias weights mapping for selection models
    if (.is_priors_weightfunction(priors)) {

      # extract publication bias priors
      priors_bias <- priors[["outcome"]][["bias"]]

      # in the case of RoBMA, the prior is formated as prior_mixture (list of priors) while selection models, PET, and PEESE use a single prior
      # as such standardize them to always behave like a list
      if (!BayesTools::is.prior.mixture(priors_bias)) {
        priors_bias <- list(priors_bias)
        # all fitting is standardized to dwnorm_mix function
        # if only a simple prior is defined, we need to insert bias_indicator manually
        fit_data$bias_indicator <- 1
      }

      # create the weightfunction mapping for effect size thresholds
      steps   <- BayesTools::weightfunctions_mapping(priors_bias[sapply(priors_bias, BayesTools::is.prior.weightfunction)], cuts_only = TRUE, one_sided = TRUE)
      steps   <- rev(steps)[c(-1, -length(steps))]
      crit_yi <- .create_yi_cutoffs(fit_data[["yi"]], fit_data[["sei"]], list(distribution = "one.sided", parameters = list(steps = steps)))

      # create the weightfunction mapping to weights (transform all weight functions to one-sided)
      crit_yi_mapping     <- matrix(0, nrow = length(priors_bias), ncol = ncol(crit_yi))
      crit_yi_mapping_max <- rep(0, length(priors_bias))
      for (i in seq_along(priors_bias)) {
        if (BayesTools::is.prior.weightfunction(priors_bias[[i]])) {
          ### the following subsetting allows us "merge" steps with equal weights due to the construction
          # specify indexes of the relevant steps
          this_steps <- .get_one_sided_cuts(priors_bias[[i]])
          crit_yi_mapping[i,1:length(this_steps)] <- which(steps %in% this_steps)
          crit_yi_mapping_max[i] <- length(this_steps)
        }
      }

      fit_data$crit_yi             <- t(crit_yi)
      fit_data$crit_yi_mapping     <- t(crit_yi_mapping)
      fit_data$crit_yi_mapping_max <- crit_yi_mapping_max
    }

  } else if (.data_outcome_type(data) == "bin") {
    # always include ai, ci, n1i, and n2i
    fit_data <- list(
      ai  = data[["outcome"]][["ai"]],
      ci  = data[["outcome"]][["ci"]],
      n1i = data[["outcome"]][["n1i"]],
      n2i = data[["outcome"]][["n2i"]]
    )
  } else if (.data_outcome_type(data) == "pois") {
    # always include x1i, x2i, t1i, and t2i
    fit_data <- list(
      x1i = data[["outcome"]][["x1i"]],
      x2i = data[["outcome"]][["x2i"]],
      t1i = data[["outcome"]][["t1i"]],
      t2i = data[["outcome"]][["t2i"]]
    )
  }

  ### add common data
  # add number of estimates
  fit_data[["K"]] <- nrow(data[["outcome"]])

  # add study_ids for 3lvl models
  if (.is_data_multilevel(data)) {
    fit_data[["study_ids"]] <- data[["outcome"]][["study_ids"]]
  }

  # add weights for weighted models
  if (.is_data_weights(data)) {
    fit_data[["weights"]] <- data[["outcome"]][["study_ids"]]
  }



  return(fit_data)
}
.create_fit_formula_list       <- function(data, parameter) {

  formula <- attr(data[[parameter]], "formula")

  # for scale regression - modify formula for exponential parameterization
  # (required for the exponential parameterization trick in BayesTools::JAGS_fit)
  if (parameter == "scale") {

    # check whether the intercept was removed from the formula
    # if so, warn the user and return it back (intercepts cannot be omitted from scale models)
    if (attr(terms(formula), "intercept") == 0) {
      warning("Intercept cannot be omitted from scale models (the regression estimates a multiplicative constant for the intercept). The intercept removal term has been ignored.", call. = FALSE)
      formula <- BayesTools:::.add_intercept_to_formula(formula)
    }

    # add the corresponding attribute
    attr(formula, "log(intercept)") <- TRUE
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
  outcome_type      <- .data_outcome_type(data)
  effect_direction  <- .data_effect_direction(data)

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
    model_syntax <- paste0(model_syntax, "tau_within  = tau * sqrt(1-rho)\n")
    model_syntax <- paste0(model_syntax, "tau_between = tau * sqrt(rho)\n")
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
  # log_tau scale regression with zero intercept
  if (is_scale) {
    model_syntax <- paste0(model_syntax, "  tau[i] = exp(log_tau[i])\n")
  }
  # for multilevel: specify heterogeneity allocation & dispatch estimate-specific/common parameter
  # tau_within  - within-study variance (estimate-level variance)
  # tau_between - between-study variance (study-level variance)
  # for rest: dispatch estimate-specific/common parameter
  if (is_multilevel) {
    if (is_scale) {
      model_syntax <- paste0(model_syntax, "  tau_within[i]  = tau[i] * sqrt(1-rho)\n")
      model_syntax <- paste0(model_syntax, "  tau_between[i] = tau[i] * sqrt(rho)\n")
      tau_estimate <- "tau_within[i]"
    } else {
      # the variance allocation performed outside of the loop
      tau_estimate <- "tau_within"
    }
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
  # separate likelihoods for normal / selection / binomial models
  # subversion for weighted/unweighted likelihoods
  if (outcome_type == "norm") {

    # selection and normal models for norm data outcome
    if (is_weightfunction) {
      if (is_weights) {
        model_syntax <- paste0(
          model_syntax, "  yi[i] ~ dwwnorm_mix(", mu_estimate, ",", "sqrt", tau2_estimate,
          ", crit_yi[,i], omega, crit_yi_mapping[,bias_indicator], crit_yi_mapping_max[bias_indicator], weight[i])\n")
      } else {
        model_syntax <- paste0(
          model_syntax, "  yi[i] ~ dwnorm_mix(", mu_estimate, ",", "sqrt", tau2_estimate,
          ", crit_yi[,i], omega, crit_yi_mapping[,bias_indicator], crit_yi_mapping_max[bias_indicator])\n")
      }
    } else {
      if (is_weights) {
        model_syntax <- paste0(model_syntax, "  yi[i] ~ dwnorm(", mu_estimate, ",", "1/", tau2_estimate, ", weight[i])\n")
      } else {
        model_syntax <- paste0(model_syntax, "  yi[i] ~ dnorm(",  mu_estimate, ",", "1/", tau2_estimate, ")\n")
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
  if (.is_data_mods(data)) {
    fit_formula_list[["mu"]]       <- .create_fit_formula_list(data = data, parameter = "mods")
    fit_formula_data_list[["mu"]]  <- .create_fit_formula_data_list(data = data, parameter = "mods")
    fit_formula_prior_list[["mu"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "mods")
    fit_formula_scale_list[["mu"]] <- .data_standardize_continuous_predictors(data)
  }

  ### add heterogeneity regressions
  if (.is_data_scale(data)) {
    fit_formula_list[["log_tau"]]       <- .create_fit_formula_list(data = data, parameter = "scale")
    fit_formula_data_list[["log_tau"]]  <- .create_fit_formula_data_list(data = data, parameter = "scale")
    fit_formula_prior_list[["log_tau"]] <- .create_fit_formula_prior_list(priors = priors, parameter = "scale")
    fit_formula_scale_list[["log_tau"]] <- .data_standardize_continuous_predictors(data)
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
.create_yi_cutoffs <- function(yi, sei, prior) {

  # get a matrix of test-statistics for the critical values
  crit_zi <- matrix(ncol = 0, nrow = length(yi))

  # compute the thresholds
  for(step in prior$parameters$steps){
    if(prior$distribution == "one.sided"){
      crit_zi <- cbind(crit_zi, stats::qnorm(step,   lower.tail = FALSE))
    }else if(prior$distribution == "two.sided"){
      crit_zi <- cbind(crit_zi, stats::qnorm(step/2, lower.tail = FALSE))
    }
  }

  # transform the thresholds back to effect sizes
  crit_yi <- crit_zi * matrix(sei, ncol = ncol(crit_zi), nrow = nrow(crit_zi))

  return(crit_yi)
}

.is_priors_PET            <- function(priors) {

  if (is.null(priors[["outcome"]][["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["outcome"]][["bias"]]))
    return(any(sapply(priors[["outcome"]][["bias"]], is.prior.PET)))

  return(is.prior.PET(priors[["outcome"]][["bias"]]))
}
.is_priors_PEESE          <- function(priors) {

  if (is.null(priors[["outcome"]][["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["outcome"]][["bias"]]))
    return(any(sapply(priors[["outcome"]][["bias"]], is.prior.PEESE)))

  return(is.prior.PEESE(priors[["outcome"]][["bias"]]))
}
.is_priors_weightfunction <- function(priors) {

  if (is.null(priors[["outcome"]][["bias"]]))
    return(FALSE)

  if (is.prior.mixture(priors[["outcome"]][["bias"]]))
    return(any(sapply(priors[["outcome"]][["bias"]], is.prior.weightfunction)))

  return(is.prior.weightfunction(priors[["outcome"]][["bias"]]))
}
.is_priors_bias           <- function(priors) {
  return(.is_priors_PET(priors) || .is_priors_PEESE(priors) || .is_priors_weightfunction(priors))
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
.data_outcome_type        <- function(data) {
  return(attr(data, "outcome_type"))
}
.data_measure             <- function(data) {
  return(attr(data, "measure"))
}
.data_effect_direction    <- function(data) {
  return(attr(data, "effect_direction"))
}
.data_standardize_continuous_predictors <- function(data) {
  return(attr(data, "standardize_continuous_predictors"))
}

.is_PET            <- function(object) .is_priors_PET(object[["priors"]])
.is_PEESE          <- function(object) .is_priors_PEESE(object[["priors"]])
.is_weightfunction <- function(object) .is_priors_weightfunction(object[["priors"]])
.is_bias           <- function(object) .is_priors_bias(object[["priors"]])
.is_multilevel     <- function(object) .is_data_multilevel(object[["data"]])
.is_mods           <- function(object) .is_data_mods(object[["data"]])
.is_scale          <- function(object) .is_data_scale(object[["data"]])
.is_weights        <- function(object) .is_data_weights(object[["data"]])
.outcome_type      <- function(object) .data_outcome_type(object[["data"]])
.measure           <- function(object) .data_measure(object[["data"]])
.effect_direction  <- function(object) .data_effect_direction(object[["data"]])
.standardize_continuous_predictors  <- function(object) .data_standardize_continuous_predictors(object[["data"]])
