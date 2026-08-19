### check specified prior distributions whether they follow basic requirements
.check_prior.effect                   <- function(prior) {

  prior <- .check_prior.simple_or_point(prior, prior_name = "prior_effect")

  return(prior)
}
.check_prior.heterogeneity            <- function(prior) {

  if (.is_prior_random(prior)) {
    return(prior)
  }

  prior <- .check_prior.simple_or_point(prior, prior_name = "prior_heterogeneity")

  # check range restriction
  if (BayesTools::is.prior.point(prior) && mean(prior) < 0) {
    stop("Point prior distribution for 'prior_heterogeneity' must be at least 0.", call. = FALSE)
  }
  if (prior[["truncation"]][["lower"]] < 0) {
    warning("Lower truncation point for 'prior_heterogeneity' prior distribution must be at least 0. The prior distribution was modified.",
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["lower"]] <- 0
  }

  return(prior)
}
.is_prior_random                      <- function(prior) {

  inherits(prior, "prior_random")
}
.check_prior.simple_or_point          <- function(prior, prior_name) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop(sprintf("The '%s' argument must be a prior distribution.", prior_name), call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || BayesTools::is.prior.point(prior)))
    stop(sprintf("The '%s' argument must be a either a simple or point prior distribution.", prior_name), call. = FALSE)

  return(prior)
}
.check_prior.restricted_01            <- function(prior, prior_name) {

  prior <- .check_prior.simple_or_point(prior, prior_name)

  # check range restriction
  if (BayesTools::is.prior.point(prior) && (mean(prior) < 0 || mean(prior) > 1)) {
    stop(sprintf("Point prior distribution for '%s' must be within [0, 1].", prior_name), call. = FALSE)
  }
  if (prior[["truncation"]][["lower"]] < 0) {
    warning(sprintf("Lower truncation point for '%s' prior distribution must be at least 0. The prior distribution was modified.", prior_name),
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["lower"]] <- 0
  }

  if (prior[["truncation"]][["upper"]] > 1) {
    warning(sprintf("Upper truncation point for '%s' prior distribution must be at most 1. The prior distribution was modified.", prior_name),
            immediate. = TRUE, call. = FALSE)
    prior[["truncation"]][["upper"]] <- 1
  }

  return(prior)
}
.check_prior.term                     <- function(prior) {

  # check object type
  if (!BayesTools::is.prior(prior))
    stop("The 'prior_mods' / 'prior_scale' arguments must be a prior distribution.", call. = FALSE)
  if (!(BayesTools::is.prior.simple(prior) || BayesTools::is.prior.point(prior) || BayesTools::is.prior.factor(prior)))
    stop("The 'prior_mods' / 'prior_scale' arguments for predictors must be a either a simple, factor, or point prior distribution.", call. = FALSE)

  return()
}
.check_prior_specification_conflict   <- function(prior_unit_information_sd, prior_informed_field) {

  # check that both UISD and informed priors are not specified simultaneously
  if (!missing(prior_unit_information_sd) && !missing(prior_informed_field)) {
    stop(paste0(
      "Both 'prior_unit_information_sd' and 'prior_informed_field' were specified. ",
      "These arguments are mutually exclusive: 'prior_unit_information_sd' is used to scale default priors, ",
      "while 'prior_informed_field' uses pre-specified informed prior distributions that do not rely on UISD scaling."
    ), call. = FALSE)
  }

  return()
}
.check_prior_informed_subfield        <- function(prior_informed_field, prior_informed_subfield) {

  if (missing(prior_informed_field) && !missing(prior_informed_subfield)) {
    stop(
      "'prior_informed_subfield' can only be specified together with 'prior_informed_field'.",
      call. = FALSE
    )
  }

  return()
}
.check_glmm_nuisance_prior_match      <- function(data, prior_baserate, prior_lograte) {

  if (.data_outcome_type(data) == "bin" && !missing(prior_lograte)) {
    stop(
      "'prior_lograte' is only used for Poisson GLMM models with measure = 'IRR'. ",
      "Use 'prior_baserate' for binomial GLMM models with measure = 'OR'.",
      call. = FALSE
    )
  }
  if (.data_outcome_type(data) == "pois" && !missing(prior_baserate)) {
    stop(
      "'prior_baserate' is only used for binomial GLMM models with measure = 'OR'. ",
      "Use 'prior_lograte' for Poisson GLMM models with measure = 'IRR'.",
      call. = FALSE
    )
  }

  return()
}
.get_prior_bias_type                  <- function(prior) {
  if (BayesTools::is.prior.none(prior)) {
    bias_type <- "none"
  } else if (.prior_is_selection_kernel(prior)) {
    bias_type <- "selmodel"
  } else if (BayesTools::is.prior.PET(prior)) {
    bias_type <- "PET"
  } else if (BayesTools::is.prior.PEESE(prior)) {
    bias_type <- "PEESE"
  } else {
    stop("The publication bias prior was not recognized", call. = FALSE)
  }
  return(bias_type)
}

# helper function for defining the default lograte prior for IRR models
.get_default_lograte_prior            <- function(x1i, x2i, t1i, t2i) {

  if (anyNA(x1i) || anyNA(x2i) || anyNA(t1i) || anyNA(t2i)) {
    stop("Default 'prior_lograte' requires complete Poisson event counts and person-times.", call. = FALSE)
  }

  total_events <- sum(x1i + x2i)
  total_time   <- sum(t1i + t2i)
  if (!is.finite(total_events) || total_events <= 0) {
    stop("Default 'prior_lograte' requires at least one observed Poisson event. Specify 'prior_lograte' manually for all-zero event counts.", call. = FALSE)
  }
  if (!is.finite(total_time) || total_time <= 0) {
    stop("Default 'prior_lograte' requires positive finite total person-time.", call. = FALSE)
  }

  # compute the pooled crude log-rate for the mean
  mu_prior <- log(total_events / total_time)
  if (!is.finite(mu_prior)) {
    stop("Default 'prior_lograte' could not be computed from the supplied Poisson data.", call. = FALSE)
  }

  # define the prior object
  prior <- BayesTools::prior_factor("normal", parameters = list(mean = mu_prior, sd = RoBMA.get_option("default_lograte.sd")), contrast = "independent")

  return(prior)
}

# priors rescaling function
.rescale_prior_object <- function(prior, rescale_priors) {
  # re-scaling prior distributions

  # no need if re-scaling is not specified
  if (missing(rescale_priors) || rescale_priors == 1)
    return(prior)

  # no need for point / no prior distributions
  if (BayesTools::is.prior.point(prior) || BayesTools::is.prior.none(prior))
    return(prior)

  # only specific prior distributions can be rescaled
  # (UPGRADE: once BayesTools incorporates direct scaling, extend to all kinds of priors)
  can_be_rescaled <- c("normal", "mnormal", "cauchy", "mcauchy", "t", "mt", "invgamma")
  if (!is.element(prior[["distribution"]], can_be_rescaled))
    stop(sprintf(
      paste0("The '%1$s' prior distribution cannot be rescaled using the 'rescale_priors' argument. ",
             "Only one of the following distributions can be rescaled: %2$s"),
      prior[["distribution"]],
      paste0("'", can_be_rescaled, "'", collapse = ", ")
    ), call. = FALSE)

  if (prior[["distribution"]] %in% c("normal", "mnormal")) {
    prior[["parameters"]][["sd"]] <- prior[["parameters"]][["sd"]] * rescale_priors
  } else if (prior[["distribution"]] %in% c("cauchy", "mcauchy", "t", "mt", "invgamma")) {
    prior[["parameters"]][["scale"]] <- prior[["parameters"]][["scale"]] * rescale_priors
  }

  return(prior)
}

### prior specification functions
.check_component_prior_presence <- function(supplied, present, argument, component) {

  if (isTRUE(supplied) && !isTRUE(present)) {
    stop(
      gettextf(
        "The '%1$s' argument can be used only when '%2$s' is specified.",
        argument,
        component
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.check_and_list_priors.brma <- function(
    prior_effect, prior_heterogeneity,
    prior_mods, prior_scale,
    prior_heterogeneity_allocation, prior_bias,
    prior_baserate, prior_lograte,
    rescale_priors,
    prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield,
    data, measure,
    bias_type, steps) {

  ### check input
  if (!missing(rescale_priors))
    BayesTools::check_real(rescale_priors, "rescale_priors", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_unit_information_sd))
    BayesTools::check_real(prior_unit_information_sd, "prior_unit_information_sd", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_informed_field))
    BayesTools::check_char(prior_informed_field, "prior_informed_field", allow_values = c("medicine"), allow_NA = FALSE)
  if (!missing(prior_informed_subfield))
    BayesTools::check_char(prior_informed_subfield, "prior_informed_subfield", allow_NA = FALSE)
  .check_prior_informed_subfield(prior_informed_field, prior_informed_subfield)
  .check_prior_specification_conflict(prior_unit_information_sd, prior_informed_field)
  .check_glmm_nuisance_prior_match(data, prior_baserate, prior_lograte)
  .check_component_prior_presence(
    supplied  = !missing(prior_mods),
    present   = .is_data_mods(data),
    argument  = "prior_mods",
    component = "mods"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_scale),
    present   = .is_data_scale(data),
    argument  = "prior_scale",
    component = "scale"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_heterogeneity_allocation),
    present   = .is_data_multilevel(data),
    argument  = "prior_heterogeneity_allocation",
    component = "cluster"
  )
  if (!missing(steps))
    BayesTools::check_real(steps, "steps", lower = 0, upper = 1, allow_bound = FALSE, check_length = FALSE, allow_NA = FALSE)

  ### set prior distributions
  measure       <- .data_measure(data)
  prior_outcome <- list()
  is_random     <- .is_data_random(data)
  is_random_prior_heterogeneity <- !missing(prior_heterogeneity) &&
    .is_prior_random(prior_heterogeneity)
  random_prior_needs_default_scale <- is_random_prior_heterogeneity &&
    !.assign_prior.random_has_scale(prior_heterogeneity)

  if (is_random_prior_heterogeneity && !is_random) {
    stop(
      "'prior_heterogeneity = BayesTools::prior_random(...)' can be used only when 'random' is specified.",
      call. = FALSE
    )
  }

  prior_outcome[["mu"]] <- .assign_prior.simple(
      prior = prior_effect, parameter = "effect", measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    )
  if (is_random_prior_heterogeneity &&
      (.is_data_scale(data) || random_prior_needs_default_scale)) {
    prior_outcome[["tau"]] <- .assign_prior.simple(
      parameter = "heterogeneity", measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    )
  } else if (!is_random_prior_heterogeneity) {
    prior_outcome[["tau"]] <- .assign_prior.simple(
      prior = prior_heterogeneity, parameter = "heterogeneity", measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    )
  }

  # only for multilevel models
  if (.is_data_multilevel(data)) {
    prior_outcome[["rho"]]   <- .assign_prior.heterogeneity_allocation(prior = prior_heterogeneity_allocation)
    # gamma is an auxiliary parameter for non-centered random-effects parameterization (cluster-level)
    prior_outcome[["gamma"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of OR
  if (.data_outcome_type(data) == "bin") {
    prior_outcome[["pi"]]    <- .assign_prior.baserate(prior = prior_baserate)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of IRR
  if (.data_outcome_type(data) == "pois") {
    prior_outcome[["phi"]] <- .assign_prior.lograte(prior = prior_lograte, data)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for selection / PET / PEESE models
  if (!missing(bias_type)) {
    prior_outcome[["bias"]] <- .assign_prior.bias(
      prior = prior_bias, measure = measure, data = data, prior_unit_information_sd = prior_unit_information_sd,
      bias_type = bias_type, steps = steps)
  }


  if (.is_data_mods(data)) {
    prior_mods <- .assign_prior_list.terms(
      prior_list = prior_mods, prior_intercept = prior_outcome[["mu"]], parameter = "mods",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["mu"]] <- NULL # the prior is forwarded through "mods" to the intercept
  } else {
    prior_mods <- NULL
  }

  prior_random_sd <- prior_outcome[["tau"]]
  if (.is_data_scale(data)) {
    prior_scale <- .assign_prior_list.scale(
      prior_list = prior_scale, prior_intercept = prior_outcome[["tau"]],
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["tau"]] <- NULL # the prior is forwarded through "scale" to the log(intercept)
  } else {
    prior_scale <- NULL
  }

  prior_location <- .assign_prior.location(
    prior_effect = prior_outcome[["mu"]],
    prior_mods   = prior_mods,
    data         = data
  )

  if (is_random) {
    preview_defaults <- BayesTools::prior_random(sd = prior_random_sd)
    preview_prior <- if (is_random_prior_heterogeneity) {
      .assign_prior.random_complete(
        prior    = prior_heterogeneity,
        defaults = preview_defaults
      )
    } else {
      preview_defaults
    }
    random_formula_design <- .compile_prior_random_formula_design(
      data           = data,
      prior_location = prior_location,
      prior_random   = preview_prior
    )
    prior_random <- if (is_random_prior_heterogeneity) {
      if (random_prior_needs_default_scale) {
        random_defaults <- .assign_prior.random(
          prior          = prior_random_sd,
          formula_design = random_formula_design,
          sd_sources     = if (.is_data_scale(data)) {
            .assign_prior.random_scale_sources(data)
          } else {
            NULL
          }
        )
        .assign_prior.random_complete(
          prior    = prior_heterogeneity,
          defaults = random_defaults
        )
      } else {
        prior_heterogeneity
      }
    } else {
      .assign_prior.random(
        prior          = prior_random_sd,
        formula_design = random_formula_design,
        sd_sources     = if (.is_data_scale(data)) {
          .assign_prior.random_scale_sources(data)
        } else {
          NULL
        }
      )
    }
    prior_outcome[["tau"]] <- NULL
  } else {
    prior_random <- NULL
  }

  if (is_random) {
    prior_outcome[["mu"]] <- NULL
    .validate_prior_random(
      data           = data,
      prior_location = prior_location,
      prior_random   = prior_random
    )
  }

  priors <- list(
    outcome  = prior_outcome,
    mods     = prior_mods,
    scale    = prior_scale,
    random   = prior_random,
    location = prior_location
  )
  .check_glmm_no_bias_priors(data, priors)

  return(priors)
}


.check_and_list_priors.RoBMA <- function(
    prior_effect, prior_heterogeneity,
    prior_mods, prior_scale,
    prior_heterogeneity_allocation, prior_bias,
    prior_effect_null, prior_heterogeneity_null,
    prior_mods_null, prior_scale_null,
    prior_heterogeneity_allocation_null, prior_bias_null,
    prior_baserate, prior_lograte,
    rescale_priors,
    prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield,
    data, model_type) {

  ### check input
  if (!missing(rescale_priors))
    BayesTools::check_real(rescale_priors, "rescale_priors", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_unit_information_sd))
    BayesTools::check_real(prior_unit_information_sd, "prior_unit_information_sd", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if (!missing(prior_informed_field))
    BayesTools::check_char(prior_informed_field, "prior_informed_field", allow_values = c("medicine"), allow_NA = FALSE)
  if (!missing(prior_informed_subfield))
    BayesTools::check_char(prior_informed_subfield, "prior_informed_subfield", allow_NA = FALSE)
  .check_prior_informed_subfield(prior_informed_field, prior_informed_subfield)
  .check_prior_specification_conflict(prior_unit_information_sd, prior_informed_field)
  .check_glmm_nuisance_prior_match(data, prior_baserate, prior_lograte)
  .check_component_prior_presence(
    supplied  = !missing(prior_mods),
    present   = .is_data_mods(data),
    argument  = "prior_mods",
    component = "mods"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_mods_null),
    present   = .is_data_mods(data),
    argument  = "prior_mods_null",
    component = "mods"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_scale),
    present   = .is_data_scale(data),
    argument  = "prior_scale",
    component = "scale"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_scale_null),
    present   = .is_data_scale(data),
    argument  = "prior_scale_null",
    component = "scale"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_heterogeneity_allocation),
    present   = .is_data_multilevel(data),
    argument  = "prior_heterogeneity_allocation",
    component = "cluster"
  )
  .check_component_prior_presence(
    supplied  = !missing(prior_heterogeneity_allocation_null),
    present   = .is_data_multilevel(data),
    argument  = "prior_heterogeneity_allocation_null",
    component = "cluster"
  )
  if (!missing(model_type))
    BayesTools::check_char(model_type, "model_type", allow_values = c("2w", "6w", "PP", "PSMA"))

  ### set prior distributions
  measure       <- .data_measure(data)
  prior_outcome <- list()

  prior_outcome[["mu"]] <- .assign_prior.simple_mixture(
    prior = prior_effect, prior_null = prior_effect_null, parameter = "effect", measure = measure,
    data = data, prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors
  )
  prior_outcome[["tau"]] <- .assign_prior.simple_mixture(
    prior = prior_heterogeneity, prior_null = prior_heterogeneity_null, parameter = "heterogeneity", measure = measure,
    data = data, prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors
  )
  # only for multilevel models
  if (.is_data_multilevel(data)) {
    prior_outcome[["rho"]]   <- .assign_prior.heterogeneity_allocation_mixture(
      prior = prior_heterogeneity_allocation, prior_null = prior_heterogeneity_allocation_null)
    # gamma is an auxiliary parameter for non-centered random-effects parameterization (cluster-level)
    prior_outcome[["gamma"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of OR
  if (.data_outcome_type(data) == "bin") {
    prior_outcome[["pi"]]    <- .assign_prior.baserate(prior = prior_baserate)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for glmm models of IRR
  if (.data_outcome_type(data) == "pois") {
    prior_outcome[["phi"]] <- .assign_prior.lograte(prior = prior_lograte, data)
    # theta is an auxiliary parameter for non-centered random-effects parameterization (for estimate-level)
    # (these are marginalized in the normal model)
    prior_outcome[["theta"]] <- prior_factor("normal", parameters = list("mean" = 0, "sd" = 1), contrast = "independent")
  }
  # only for RoBMA
  if (!missing(model_type)) {
    prior_outcome[["bias"]] <- .assign_prior.bias_mixture(
      prior = prior_bias, prior_null = prior_bias_null, measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd, model_type = model_type
    )
  }


  if (.is_data_mods(data)) {
    prior_mods <- .assign_prior_list.terms_mixture(
      prior_list = prior_mods, prior_list_null = prior_mods_null, prior_intercept = prior_outcome[["mu"]], parameter = "mods",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["mu"]] <- NULL # the prior is forwarded through "mods" to the intercept
  } else {
    prior_mods <- NULL
  }

  if (.is_data_scale(data)) {
    prior_scale <- .assign_prior_list.terms_mixture(
      prior_list = prior_scale, prior_list_null = prior_scale_null, prior_intercept = prior_outcome[["tau"]], parameter = "scale",
      measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors)
    prior_outcome[["tau"]] <- NULL # the prior is forwarded through "scale" to the log(intercept)
  } else {
    prior_scale <- NULL
  }

  priors <- list(
    outcome  = prior_outcome,
    mods     = prior_mods,
    scale    = prior_scale
  )
  .check_glmm_no_bias_priors(data, priors)

  return(priors)
}
