.assign_prior.simple                   <- function(
    prior, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {


  ### always use the user-specified prior if possible
  if (!missing(prior)) {

    # allow using FALSE/NULL arguments to set spike(0) prior distributions
    if (is.null(prior) || isFALSE(prior)) {
      return(BayesTools::prior("spike", parameters = list(0)))
    }

    # apply rescaling if specified
    prior <- .rescale_prior_object(prior, rescale_priors)

    # check the specified prior distribution
    prior <- switch(
      parameter,
      "effect"        = .check_prior.effect(prior),
      "heterogeneity" = .check_prior.heterogeneity(prior)
    )

    return(prior)
  }

  ### use informed prior distributions if specified
  if (!missing(prior_informed_field)) {

    if (prior_informed_field == "medicine") {

      # input translation for BayesTools::prior_informed
      if (missing(prior_informed_subfield)) {
        prior_name <- "Cochrane"
      } else {
        prior_name <- prior_informed_subfield
      }

      type <- switch(
        tolower(measure),
        "smd"   = "smd",
        "or"    = "logOR",
        "rr"    = "logRR",
        "rd"    = "RD",
        "hr"    = "logHR"
      )
      if (is.null(type)) {
        stop("Informed medicine priors are only available for measures 'SMD', 'OR', 'RR', 'RD', and 'HR'.", call. = FALSE)
      }

      # simple informed priors from BayesTools package
      prior <- BayesTools::prior_informed(prior_name, parameter = parameter, type = type)
    }

    # apply rescaling if specified
    prior <- .rescale_prior_object(prior, rescale_priors)
    return(prior)
  }

  ### use prior_unit_information_sd otherwise
  if (missing(prior_unit_information_sd)) {
    # if not passed directly, obtain from data
    prior_unit_information_sd <- .get_unit_information_sd(data = data, measure = measure)
  }

  # get the prior standard deviation based on unit information
  prior_sd <- switch(
    parameter,
    "effect"        = prior_unit_information_sd * RoBMA.get_option("default_UISD.effect"),
    "heterogeneity" = prior_unit_information_sd * RoBMA.get_option("default_UISD.heterogeneity")
  )

  # create to corresponding prior
  prior <- switch(
    parameter,
    "effect"        = BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = prior_sd)),
    "heterogeneity" = BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = prior_sd), truncation = list(0, Inf))
  )

  # apply rescaling if specified
  prior <- .rescale_prior_object(prior, rescale_priors)

  return(prior)
}

.assign_prior.baserate                 <- function(prior) {

  # use default settings if prior distribution is not specified
  if (missing(prior) || is.null(prior)) {
    prior <- BayesTools::prior_factor("beta", parameters = list("alpha" = 1, "beta" = 1), contrast = "independent")
    return(prior)
  }

  # check the user specified prior distribution
  prior <- .check_prior.restricted_01(prior, prior_name = "baserate")

  # transform the prior into independent factor prior
  prior <- BayesTools::prior_factor(prior[["distribution"]], parameters = prior[["parameters"]], truncation = prior[["truncation"]], contrast = "independent")

  return(prior)
}
.assign_prior.lograte                  <- function(prior, data) {
  # the log-rate is unbounded, so use a proper data-based default prior
  # rather than an improper flat nuisance-parameter prior

  if (missing(prior) || is.null(prior)) {
    return(.get_unit_information_lograte_prior(
      x1i = data[["outcome"]][["x1i"]],
      x2i = data[["outcome"]][["x2i"]],
      t1i = data[["outcome"]][["t1i"]],
      t2i = data[["outcome"]][["t2i"]])
    )
  }

  # check the user specified prior distribution
  prior <- .check_prior.simple_or_point(prior, prior_name = "prior_lograte")

  # transform the prior into independent factor prior
  prior <- BayesTools::prior_factor(prior[["distribution"]], parameters = prior[["parameters"]], truncation = prior[["truncation"]], contrast = "independent")

  return(prior)
}
.get_PEESE_UISD_ratio                  <- function(measure, data, prior_unit_information_sd) {

  UISD_SMD <- .get_unit_information_sd.known("SMD")
  if (missing(prior_unit_information_sd)) {
    UISD_measure <- .get_unit_information_sd(data = data, measure = measure)
  } else {
    UISD_measure <- prior_unit_information_sd
  }

  return(UISD_SMD / UISD_measure)
}
.assign_prior.bias                     <- function(prior, measure, data, prior_unit_information_sd, bias_type, steps) {

  # assigns either a weight function or PET/PEESE model based on bias_type
  # there is no scaling of the publication bias priors based on rescale_priors
  # (only PEESE prior is rescaled to the match the UISD of the measure)

  if (bias_type == "selmodel") {

    # specify default settings / check the prior distribution
    if (missing(prior) || is.null(prior)) {
      if (missing(steps)) {
        steps <- c(0.025)
      }
      prior <- BayesTools::prior_weightfunction(
        "one-sided",
        steps   = steps,
        weights = BayesTools::wf_cumulative(
          alpha = rep(RoBMA.get_option("default_bias_weightfunction.alpha"), length(steps) + 1)
        )
      )
    } else if (!.prior_is_selection_kernel(prior)) {
      stop(
        "'prior_bias' must be a `prior_weightfunction` or `prior_bias` object",
        call. = FALSE
      )
    } else if (.prior_has_phacking(prior)) {
      .selection_stop_phacking_deferred()
    }

    return(prior)
  }

  if (bias_type == "PET") {

    # specify default settings / check the prior distribution
    if (missing(prior) || is.null(prior)) {
      # PET prior is invariant to different effect size types
      prior <- BayesTools::prior_PET("cauchy", parameters = list("location" = 0, "scale" = RoBMA.get_option("default_bias_PET.scale")), truncation = list(0, Inf))
    } else if (!BayesTools::is.prior.PET(prior)) {
      stop("'prior_bias' must be a `prior_PET` object", call. = FALSE)
    }

    return(prior)
  }

  if (bias_type == "PEESE") {

    # specify default settings / check the prior distribution
    if (missing(prior) || is.null(prior)) {
      ### PEESE prior is inversely related to the specified effect size
      # for SMD we use Cauchy with scale 1
      # for other effect sizes, we re-scale by the corresponding UISD
      UISD_ratio <- .get_PEESE_UISD_ratio(
        measure = measure, data = data, prior_unit_information_sd = prior_unit_information_sd)
      prior <- BayesTools::prior_PEESE("cauchy", parameters = list("location" = 0, "scale" = RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf))
    } else if (!BayesTools::is.prior.PEESE(prior)) {
      stop("'prior_bias' must be a `prior_PEESE` object", call. = FALSE)
    }

    return(prior)
  }

  if (bias_type == "none") {
    # for simplified handling in RoBMA
    return(prior)
  }
}

### assign mixture prior distributions
.assign_prior.simple_mixture                   <- function(
    prior, prior_null, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {

  ### check the alternative prior input
  # prior can be either
  # - unspecified (set default as in brma)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior) || BayesTools::is.prior(prior)) {
    priors_alt <- list(.assign_prior.simple(
      prior = prior, parameter = parameter, measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    ))
  } else if (is.null(prior) || isFALSE(prior)) {
    priors_alt <- list()
  } else if (is.list(prior) && all(sapply(prior, BayesTools::is.prior))) {
    priors_alt <- prior
    for (i in seq_along(priors_alt)) {
      priors_alt[[i]] <- .assign_prior.simple(
        prior = priors_alt[[i]], parameter = parameter, measure = measure,
        data = data, prior_unit_information_sd = prior_unit_information_sd,
        prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
        rescale_priors = rescale_priors
      )
    }
  } else {
    stop(gettextf("The 'prior_%1$s' must be either prior distribution or a list of prior distributions.", parameter), call. = FALSE)
  }

  ### check the null prior input
  # prior can be either
  # - empty (set point prior as default)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior_null)) {
    priors_null <- list(BayesTools::prior("spike", parameters = list(0)))
  } else if (is.null(prior_null) || isFALSE(prior_null)) {
    priors_null <- list()
  } else if (BayesTools::is.prior(prior_null)) {
    priors_null <- list(.assign_prior.simple(
      prior = prior_null, parameter = parameter, measure = measure,
      data = data, prior_unit_information_sd = prior_unit_information_sd,
      prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
      rescale_priors = rescale_priors
    ))
  } else if (is.list(prior_null) && all(sapply(prior_null, BayesTools::is.prior))) {
    priors_null <- prior_null
    for (i in seq_along(priors_null)) {
      priors_null[[i]] <- .assign_prior.simple(
        prior = priors_null[[i]], parameter = parameter, measure = measure,
        data = data, prior_unit_information_sd = prior_unit_information_sd,
        prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
        rescale_priors = rescale_priors
      )
    }
  } else {
    stop(gettextf("The 'prior_%1$s_null' must be either prior distribution or a list of prior distributions.", parameter), call. = FALSE)
  }

  if (length(priors_null) + length(priors_alt) == 0) {
    stop(gettextf("At least one prior distribution needs to be defined for '%1$s'.", parameter), call. = FALSE)
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
.default_prior.bias_alt                         <- function(model_type, measure, data, prior_unit_information_sd) {

  if (model_type == "2w") {
    return(list(
      BayesTools::prior_weightfunction("two-sided", c(0.05),       BayesTools::wf_cumulative(c(1, 1)),    prior_weights = 1/2),
      BayesTools::prior_weightfunction("two-sided", c(0.05, 0.10), BayesTools::wf_cumulative(c(1, 1, 1)), prior_weights = 1/2)
    ))
  }

  if (model_type == "6w") {
    return(list(
      BayesTools::prior_weightfunction("two-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/6),
      BayesTools::prior_weightfunction("two-sided", c(0.05, 0.10),       BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05),      BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.05, 0.5),        BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/6),
      BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05, 0.5), BayesTools::wf_cumulative(c(1, 1, 1, 1)), prior_weights = 1/6)
    ))
  }

  UISD_ratio <- .get_PEESE_UISD_ratio(
    measure = measure, data = data, prior_unit_information_sd = prior_unit_information_sd)

  if (model_type == "PP") {
    return(list(
      BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf), prior_weights = 1/2),
      BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf), prior_weights = 1/2)
    ))
  }

  return(list(
    BayesTools::prior_weightfunction("two-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/12),
    BayesTools::prior_weightfunction("two-sided", c(0.05, 0.10),       BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.05),             BayesTools::wf_cumulative(c(1, 1)),       prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05),      BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.05, 0.5),        BayesTools::wf_cumulative(c(1, 1, 1)),    prior_weights = 1/12),
    BayesTools::prior_weightfunction("one-sided", c(0.025, 0.05, 0.5), BayesTools::wf_cumulative(c(1, 1, 1, 1)), prior_weights = 1/12),
    BayesTools::prior_PET(distribution = "Cauchy",   parameters = list(0, RoBMA.get_option("default_bias_PET.scale")),                truncation = list(0, Inf), prior_weights = 1/4),
    BayesTools::prior_PEESE(distribution = "Cauchy", parameters = list(0, RoBMA.get_option("default_bias_PEESE.scale") * UISD_ratio), truncation = list(0, Inf), prior_weights = 1/4)
  ))
}
.assign_prior.bias_mixture                     <- function(
    prior, prior_null, measure, data, prior_unit_information_sd, model_type) {

  ### use precanned publication bias adjustment models
  # if neither `prior` or `prior_null` input is specified use the `model_type`
  # there is no scaling of the publication bias priors based on rescale_priors
  # (only PEESE prior is rescaled to the match the UISD of the measure)

  if (missing(prior) && missing(prior_null)) {
    priors_null <- list(BayesTools::prior_none(prior_weights = 1))
    priors_alt  <- .default_prior.bias_alt(
      model_type = model_type, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd)
    prior <- BayesTools::prior_mixture(
      prior_list = c(priors_null, priors_alt),
      is_null    = c(TRUE, rep(FALSE, length(priors_alt)))
    )

    return(prior)
  }


  ### use the user specified prior distributions if specified

  ### check the alternative prior input
  # prior can be either
  # - unspecified (set alternative part of PSMA)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior)) {
    priors_alt <- .default_prior.bias_alt(
      model_type = model_type, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd)
  } else if (is.null(prior) || isFALSE(prior)) {
    priors_alt <- list()
  } else if (BayesTools::is.prior(prior)) {
    priors_alt <- list(.assign_prior.bias(
      prior = prior, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(prior)))
  } else if (is.list(prior) && all(sapply(prior, BayesTools::is.prior))) {
    priors_alt <- prior
    for (i in seq_along(priors_alt)) {
      priors_alt[[i]] <- .assign_prior.bias(
        prior = priors_alt[[i]], measure = measure, data = data,
        prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(priors_alt[[i]]))
    }
  } else {
    stop("The 'prior_bias' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  ### check the null prior input
  # prior can be either
  # - empty (set point prior as default)
  # - NULL/FALSE (omit prior distribution)
  # - single prior distribution
  # - list of prior distributions
  # each of these prior distributions needs to satisfy the same properties as in brma
  # if only a single distribution is specified, place it in a list for a simpler treatment
  if (missing(prior_null)) {
    priors_null <- list(BayesTools::prior_none(prior_weights = 1))
  } else if (is.null(prior_null) || isFALSE(prior_null)) {
    priors_null <- list()
  } else if (BayesTools::is.prior(prior_null)) {
    priors_null <- list(.assign_prior.bias(
      prior = prior_null, measure = measure, data = data,
      prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(prior_null)))
  } else if (is.list(prior_null) && all(sapply(prior_null, BayesTools::is.prior))) {
    priors_null <- prior_null
    for (i in seq_along(priors_null)) {
      priors_null[[i]] <- .assign_prior.bias(
        prior = priors_null[[i]], measure = measure, data = data,
        prior_unit_information_sd = prior_unit_information_sd, bias_type = .get_prior_bias_type(priors_null[[i]]))
    }
  } else {
    stop("The 'prior_bias_null' must be either prior distribution or a list of prior distributions.", call. = FALSE)
  }

  if (length(priors_null) + length(priors_alt) == 0) {
    stop("At least one prior distribution needs to be defined for 'bias'.", call. = FALSE)
  }

  if (length(priors_alt) == 0 &&
      length(priors_null) == 1 &&
      BayesTools::is.prior.none(priors_null[[1]])) {
    return(priors_null[[1]])
  }

  ### specify the corresponding mixture prior
  prior <- BayesTools::prior_mixture(
    prior_list = c(priors_null, priors_alt),
    is_null    = c(rep(TRUE, length(priors_null)), rep(FALSE, length(priors_alt)))
  )

  return(prior)
}
