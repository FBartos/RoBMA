.term_default_continuous_prior         <- function(prior) {

  if (!BayesTools::is.prior.factor(prior)) {
    return(prior)
  }

  distribution <- prior[["distribution"]]
  if (distribution %in% c("mnormal", "mcauchy", "mt")) {
    distribution <- substring(distribution, 2)
  }

  return(BayesTools::prior(
    distribution = distribution,
    parameters   = prior[["parameters"]],
    truncation   = prior[["truncation"]]
  ))
}
.term_default_factor_prior             <- function(prior, contrast) {

  if (BayesTools::is.prior.factor(prior)) {
    return(prior)
  }

  distribution <- prior[["distribution"]]
  if (!contrast %in% c("treatment", "independent") &&
      !distribution %in% c("spike", "mnormal", "mcauchy", "mt")) {
    distribution <- paste0("m", distribution)
  }

  return(BayesTools::prior_factor(
    distribution = distribution,
    parameters   = prior[["parameters"]],
    truncation   = prior[["truncation"]],
    contrast     = contrast
  ))
}
.term_default_prior_sd                 <- function(parameter, measure, data, prior_unit_information_sd) {

  if (parameter == "mods") {
    if (is.null(prior_unit_information_sd)) {
      prior_unit_information_sd <- .get_unit_information_sd(data = data, measure = measure)
    }
    return(prior_unit_information_sd * RoBMA.get_option("default_UISD.mods"))
  }

  if (parameter == "scale") {
    return(RoBMA.get_option("default_UISD.scale"))
  }
}
.term_default_continuous_normal_prior  <- function(prior_sd, rescale_priors) {

  prior <- BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = prior_sd))
  return(.rescale_prior_object(prior, rescale_priors))
}
.term_default_factor_normal_prior      <- function(prior_sd, contrast, rescale_priors) {

  if (contrast %in% c("treatment", "independent")) {
    prior <- BayesTools::prior_factor(
      distribution = "normal",
      parameters   = list("mean" = 0, "sd" = prior_sd),
      contrast     = contrast
    )
  } else {
    prior <- BayesTools::prior_factor(
      distribution = "mnormal",
      parameters   = list("mean" = 0, "sd" = prior_sd),
      contrast     = contrast
    )
  }

  return(.rescale_prior_object(prior, rescale_priors))
}
.warn_default_independent_fixed_overspecification <- function(
    formula_design, formula, explicit_prior_names, parameter, contrast) {

  if (!identical(contrast, "independent")) {
    return(invisible(FALSE))
  }

  model_terms      <- formula_design[["model_terms"]]
  model_terms_type <- formula_design[["model_terms_type"]]
  explicit_prior_names <- gsub(
    ":", "__xXx__", explicit_prior_names, fixed = TRUE
  )
  default_factor_terms <- model_terms[
    model_terms_type == "factor" &
      !model_terms %in% explicit_prior_names
  ]
  if (length(default_factor_terms) == 0L) {
    return(invisible(FALSE))
  }

  has_compiled_intercept <- "intercept" %in% model_terms
  term_positions <- match(default_factor_terms, model_terms)
  term_assign <- term_positions - as.integer(has_compiled_intercept)
  selected_assign <- term_assign
  has_formula_intercept <- attr(stats::terms(formula), "intercept") == 1L
  if (has_formula_intercept) {
    selected_assign <- c(0L, selected_assign)
  }

  model_matrix <- formula_design[["model_matrix"]]
  selected <- formula_design[["assign"]] %in% selected_assign
  affected_matrix <- model_matrix[, selected, drop = FALSE]
  affected_rank <- qr(affected_matrix)[["rank"]]
  if (affected_rank == ncol(affected_matrix)) {
    return(invisible(FALSE))
  }

  remedy <- if (has_formula_intercept) {
    paste0(
      "Suppress the intercept with '", parameter,
      " = ~ 0 + ...', use independent coding for at most one factor term, "
    )
  } else {
    "Use independent coding for at most one factor term "
  }
  warning(
    "The '", parameter,
    "' formula is overspecified by default independent factor contrasts: ",
    "the affected design has ", ncol(affected_matrix),
    " columns but rank ", affected_rank, ". ",
    remedy,
    "or set identifiable contrasts with 'prior_", parameter, "'.",
    call. = FALSE
  )

  invisible(TRUE)
}
.assign_prior_list.terms               <- function(
    prior_list, prior_intercept, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors,
    set_default_spike = FALSE, formula_parameter = NULL) {

  # set_default_spike argument is used when calling from .assign_prior_list.terms_mixture
  # to set-up prior distributions under the null hypotheses
  if (is.null(formula_parameter)) {
    formula_parameter <- switch(parameter, "mods" = "mu", "scale" = "log_tau")
  }

  if (parameter == "mods" &&
      attr(stats::terms(attr(data[[parameter]], "formula")), "intercept") == 0) {
    prior_intercept <- BayesTools::prior("spike", parameters = list(0))
  }

  ### check the user-specified priors
  user_default_prior <- NULL
  if (!missing(prior_list) && BayesTools::is.prior(prior_list)) {
    .check_prior.term(prior_list)
    user_default_prior <- prior_list
    prior_list         <- list()
  } else if (!missing(prior_list) && is.null(prior_list)) {
    prior_list <- list()
  } else if (!missing(prior_list)) {
    # check the priors and apply rescaling if specified
    for (i in seq_along(prior_list)) {
      .check_prior.term(prior_list[[i]])
      prior_list[[i]] <- .rescale_prior_object(prior_list[[i]], rescale_priors)
    }
  } else {
    # create an empty placeholder instead
    prior_list <- list()
  }

  if ("intercept" %in% names(prior_list)) {
    prior_intercept          <- prior_list[["intercept"]]
    prior_list[["intercept"]] <- NULL
  }

  ### for scale regression, the intercept cannot be spike(0) since it's used as log(intercept)
  # in exp( log(intercept) + sum(beta_i * x_i) )
  if (parameter == "scale") {
    if (BayesTools::is.prior.point(prior_intercept) && mean(prior_intercept) == 0) {
      stop("Intercept prior distribution for scale regression cannot be set to 0 since the model is parameterized as tau = exp( log(intercept) + sum(beta_i * x_i) ).", call. = FALSE)
    }
  }

  ### add intercept to be beginning of the list (parsed previously)
  prior_list <- c(list(intercept = prior_intercept), prior_list)

  ### add default priors for remaining terms if left unspecified
  # the assignment of missing priors is performed within BayesTools::JAGS_formula
  # (it's easier to check and fill in the formula on assignment when it's being parsed)

  # use user-specified default priors when available; otherwise construct defaults
  if (!is.null(user_default_prior)) {

    default_prior_continuous <- .term_default_continuous_prior(user_default_prior)
    default_prior_factor     <- .term_default_factor_prior(
      prior    = user_default_prior,
      contrast = attr(data, "set_contrast_factor_predictors")
    )

  } else if (set_default_spike) {

    # set_default_spike argument is used when calling from .assign_prior_list.terms_mixture
    # to set-up prior distributions under the null hypotheses
    default_prior_continuous <- BayesTools::prior("spike", parameters = list(0))
    default_prior_factor     <- BayesTools::prior_factor("spike", parameters = list(0), contrast = attr(data, "set_contrast_factor_predictors"))

  } else if (!missing(prior_informed_field)) {

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

      if (parameter == "mods") {
        # treat moderators by scaling them according to 1/2 of the effect
        default_prior_continuous <- BayesTools::prior_informed(prior_name, parameter = "effect", type = type)
        default_prior_continuous <- .rescale_prior_object(default_prior_continuous, RoBMA.get_option("default_informed_priors.mods"))
      } else if (parameter == "scale") {
        # scale parameters use the default scaling (no informed priors exist)
        default_prior_continuous <- BayesTools::prior("normal", parameters = list("mean" = 0, "sd" = RoBMA.get_option("default_informed_priors.scale")))
      }

      # transform the continuous prior distributions into prior distributions for factors
      if (attr(data, "set_contrast_factor_predictors") %in%
          c("treatment", "independent")) {
        default_prior_factor <- BayesTools::prior_factor(
          distribution = default_prior_continuous[["distribution"]],
          parameters   = default_prior_continuous[["parameters"]],
          contrast     = attr(data, "set_contrast_factor_predictors")
        )
      } else {
        default_prior_factor <- BayesTools::prior_factor(
          distribution = paste0("m", default_prior_continuous[["distribution"]]),
          parameters   = default_prior_continuous[["parameters"]],
          contrast     = attr(data, "set_contrast_factor_predictors")
        )
      }

    }

  } else {

    ### use prior_unit_information_sd otherwise
    if (!missing(prior_unit_information_sd)) {
      default_prior_unit_information_sd <- prior_unit_information_sd
    } else {
      default_prior_unit_information_sd <- NULL
    }

    default_prior_sd <- local({
      prior_sd <- NULL
      function() {
        if (is.null(prior_sd)) {
          prior_sd <<- .term_default_prior_sd(
            parameter                 = parameter,
            measure                   = measure,
            data                      = data,
            prior_unit_information_sd = default_prior_unit_information_sd
          )
        }
        return(prior_sd)
      }
    })

    default_prior_continuous <- function() {

      return(.term_default_continuous_normal_prior(
        prior_sd       = default_prior_sd(),
        rescale_priors = rescale_priors
      ))
    }
    default_prior_factor <- function() {

      return(.term_default_factor_normal_prior(
        prior_sd       = default_prior_sd(),
        contrast       = attr(data, "set_contrast_factor_predictors"),
        rescale_priors = rescale_priors
      ))
    }
  }

  # apply rescaling if specified
  if (!is.function(default_prior_continuous)) {
    default_prior_continuous <- .rescale_prior_object(default_prior_continuous, rescale_priors)
  }
  if (!is.function(default_prior_factor)) {
    default_prior_factor <- .rescale_prior_object(default_prior_factor, rescale_priors)
  }

  # list the priors into the default positions
  explicit_prior_names <- names(prior_list)
  prior_list[["__default_continuous"]] <- default_prior_continuous
  prior_list[["__default_factor"]]     <- default_prior_factor

  ### use BayesTools::JAGS_formula to obtain the assigned list of priors
  formula <- attr(data[[parameter]], "formula")
  formula_result <- BayesTools::JAGS_formula(
    parameter  = formula_parameter,
    data       = data[[parameter]],
    formula    = formula,
    prior_list = prior_list
  )
  .warn_default_independent_fixed_overspecification(
    formula_design       = formula_result[["formula_design"]],
    formula              = formula,
    explicit_prior_names = explicit_prior_names,
    parameter            = parameter,
    contrast             = attr(data, "set_contrast_factor_predictors")
  )
  prior_list <- formula_result[["prior_list"]]
  names(prior_list) <- BayesTools::format_parameter_names(
    parameters         = names(prior_list),
    formula_parameters = formula_parameter,
    formula_prefix     = FALSE
  )

  return(prior_list)
}

.assign_prior_list.scale <- function(
    prior_list, prior_intercept, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {

  has_prior_list <- !missing(prior_list) && !is.null(prior_list)
  if (!has_prior_list) {
    prior_list <- list()
  }

  scale_specs <- .data_scale_component_specs(data)
  if (length(scale_specs) == 0L) {
    return(NULL)
  }
  has_prior_unit_information_sd <- !missing(prior_unit_information_sd)
  has_prior_informed_field      <- !missing(prior_informed_field)
  has_prior_informed_subfield   <- !missing(prior_informed_subfield)

  is_component_scale <- inherits(data[["scale"]], "RoBMA_scale_components")
  if (is_component_scale && has_prior_list) {
    prior_list <- .assign_prior_list.scale_match_components(
      prior_list   = prior_list,
      scale_specs  = scale_specs
    )
  }

  assign_one <- function(scale_data, formula_parameter, component_prior_list = prior_list) {
    args <- list(
      prior_list        = component_prior_list,
      prior_intercept   = prior_intercept,
      parameter         = "scale",
      measure           = measure,
      data              = scale_data,
      rescale_priors    = rescale_priors,
      formula_parameter = formula_parameter
    )
    if (has_prior_unit_information_sd) {
      args[["prior_unit_information_sd"]] <- prior_unit_information_sd
    }
    if (has_prior_informed_field) {
      args[["prior_informed_field"]] <- prior_informed_field
    }
    if (has_prior_informed_subfield) {
      args[["prior_informed_subfield"]] <- prior_informed_subfield
    }

    do.call(.assign_prior_list.terms, args)
  }

  if (length(scale_specs) == 1L && identical(scale_specs[[1L]][["parameter"]], "log_tau")) {
    return(assign_one(data, "log_tau"))
  }

  out <- lapply(names(scale_specs), function(component) {
    scale_spec            <- scale_specs[[component]]
    scale_data            <- data
    scale_data[["scale"]] <- scale_spec[["data"]]

    component_prior_list <- if (is_component_scale && has_prior_list) {
      prior_list[[component]]
    } else {
      prior_list
    }

    assign_one(scale_data, scale_spec[["parameter"]], component_prior_list)
  })
  names(out) <- vapply(scale_specs, `[[`, character(1), "parameter")
  class(out) <- c("RoBMA_scale_priors", "list")

  return(out)
}

.assign_prior_list.scale_match_components <- function(prior_list, scale_specs) {

  if (!is.list(prior_list) || BayesTools::is.prior(prior_list)) {
    stop(
      "When 'scale' is a named list of formulas, 'prior_scale' must either ",
      "be omitted or be a named list with one prior list per scale component.",
      call. = FALSE
    )
  }
  if (is.null(names(prior_list)) || any(!nzchar(names(prior_list)))) {
    stop(
      "When 'scale' is a named list of formulas, 'prior_scale' component ",
      "lists must be named.",
      call. = FALSE
    )
  }

  expected <- names(scale_specs)
  supplied <- names(prior_list)
  component_map <- .assign_prior_list.scale_component_map(scale_specs)
  mapped_components <- unname(component_map[supplied])
  unknown_components <- supplied[is.na(mapped_components)]
  missing_components <- setdiff(
    expected,
    mapped_components[!is.na(mapped_components)]
  )
  if (length(missing_components) > 0L || length(unknown_components) > 0L) {
    stop(
      "When 'scale' is a named list of formulas, 'prior_scale' names must ",
      "exactly match the scale components. Missing: ",
      .check_and_list_data.collapse_or_none(missing_components),
      "; unknown: ",
      .check_and_list_data.collapse_or_none(unknown_components),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(supplied) || anyDuplicated(mapped_components)) {
    stop(
      "'prior_scale' component names must uniquely match scale components.",
      call. = FALSE
    )
  }

  prior_list <- prior_list[match(expected, mapped_components)]
  names(prior_list) <- expected
  for (component in expected) {
    if (is.null(prior_list[[component]])) {
      prior_list[[component]] <- list()
    }
    if (!is.list(prior_list[[component]]) ||
        BayesTools::is.prior(prior_list[[component]])) {
      stop(
        "Each 'prior_scale' component must be a list of priors for that ",
        "component's scale formula.",
        call. = FALSE
      )
    }
  }

  prior_list
}

.assign_prior_list.scale_component_map <- function(scale_specs) {

  aliases <- unlist(lapply(scale_specs, function(scale_spec) {
    scale_spec[["aliases"]]
  }), use.names = FALSE)
  components <- rep(names(scale_specs), vapply(scale_specs, function(scale_spec) {
    length(scale_spec[["aliases"]])
  }, integer(1)))

  keep <- !is.na(aliases) & nzchar(aliases)
  aliases    <- aliases[keep]
  components <- components[keep]

  duplicated_aliases <- unique(aliases[duplicated(aliases)])
  if (length(duplicated_aliases) > 0L) {
    components <- components[!aliases %in% duplicated_aliases]
    aliases    <- aliases[!aliases %in% duplicated_aliases]
  }

  stats::setNames(components, aliases)
}

# Remove invalid fixed-effect heterogeneity components for scale regression.
.remove_zero_scale_intercept_prior <- function(prior_intercept) {

  is_zero_point <- function(prior) {
    BayesTools::is.prior.point(prior) && mean(prior) == 0
  }

  if (!BayesTools::is.prior.mixture(prior_intercept)) {
    if (is_zero_point(prior_intercept)) {
      stop(
        "At least one non-zero between-study heterogeneity prior must be defined for scale regression.",
        call. = FALSE
      )
    }
    return(prior_intercept)
  }

  keep <- !vapply(prior_intercept, is_zero_point, logical(1))

  if (all(keep))
    return(prior_intercept)

  warning(
    "Spike(0) prior distribution for between-study heterogeneity tau was removed since scale regression requires non-zero between-study heterogeneity.",
    immediate. = TRUE
  )

  if (!any(keep)) {
    stop(
      "At least one non-zero between-study heterogeneity prior must be defined for scale regression.",
      call. = FALSE
    )
  }

  components <- attr(prior_intercept, "components")
  prior_list <- prior_intercept[keep]
  attributes(prior_list) <- NULL

  return(BayesTools::prior_mixture(
    prior_list = prior_list,
    components = components[keep]
  ))
}

.assign_prior_list.terms_mixture               <- function(
    prior_list, prior_list_null, prior_intercept, parameter, measure, data, prior_unit_information_sd,
    prior_informed_field, prior_informed_subfield, rescale_priors) {

  # in contrast to effect/heterogeneity prior list settings terms currently only allow a single prior for alternative
  # and a single prior for the null hypothesis for each term

  # one issue when re-using .assign_prior_list.terms is that we need to forward prior for the intercept which is
  # an already parsed mixture prior - when merging mixture priors for terms, we ought to not duplicate the intercept

  # we want to allow NULL/FALSE to omit prior distribution for a given component, .assign_prior_list.terms
  # fills all missing prior distributions with default priors
  # as such, we store names of parameters that should be omited and remove them manually
  if (!missing(prior_list) && !BayesTools::is.prior(prior_list) && length(prior_list) > 0 && is.list(prior_list)) {
    par_alt_omit <- names(prior_list)[sapply(prior_list, function(x) is.null(x) || isFALSE(x))]
    # filter out NULL/FALSE values before passing to .assign_prior_list.terms
    prior_list <- prior_list[!sapply(prior_list, function(x) is.null(x) || isFALSE(x))]
  } else {
    par_alt_omit <- NULL
  }
  if (!missing(prior_list_null) && !BayesTools::is.prior(prior_list_null) && length(prior_list_null) > 0 && is.list(prior_list_null)) {
    par_null_omit <- names(prior_list_null)[sapply(prior_list_null, function(x) is.null(x) || isFALSE(x))]
    # filter out NULL/FALSE values before passing to .assign_prior_list.terms
    prior_list_null <- prior_list_null[!sapply(prior_list_null, function(x) is.null(x) || isFALSE(x))]
  } else {
    par_null_omit <- NULL
  }


  if (parameter == "scale") {
    # in the case of scale moderation, spike at zero prior (i.e., fixed-effect model) is not possible since the intercept
    # is on a log scale
    # remove spike(0) priors and rebuild mixture metadata
    prior_intercept <- .remove_zero_scale_intercept_prior(prior_intercept)
  }
  if (parameter == "mods" &&
      attr(stats::terms(attr(data[[parameter]], "formula")), "intercept") == 0) {
    prior_intercept <- BayesTools::prior("spike", parameters = list(0))
  }

  prior_list_alt <- .assign_prior_list.terms(
    prior_list = prior_list, prior_intercept = prior_intercept, parameter = parameter,
    measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors)

  prior_list_null <- .assign_prior_list.terms(
    prior_list = prior_list_null, prior_intercept = prior_intercept, parameter = parameter,
    measure = measure, data = data,  prior_unit_information_sd = prior_unit_information_sd,
    prior_informed_field = prior_informed_field, prior_informed_subfield = prior_informed_subfield,
    rescale_priors = rescale_priors, set_default_spike = TRUE)

  available_terms <- names(prior_list_alt)
  .warn_unmatched_prior_omissions(
    omitted   = par_alt_omit,
    available = available_terms,
    argument  = paste0("prior_", parameter)
  )
  .warn_unmatched_prior_omissions(
    omitted   = par_null_omit,
    available = available_terms,
    argument  = paste0("prior_", parameter, "_null")
  )

  prior_intercept <- prior_list_alt[["intercept"]]

  ### initiate the list of prior mixtures with the intercept
  prior_list <- list(intercept = prior_intercept)
  for (par in names(prior_list_alt)) {

    # skip the already assigned intercept
    if (par == "intercept") {
      next
    }

    # create temporal holder for the parameter
    temp_prior_list  <- list()
    temp_is_null     <- NULL

    if (!par %in% par_null_omit) {
      temp_prior_list <- c(temp_prior_list, list(prior_list_null[[par]]))
      temp_is_null    <- c(temp_is_null, TRUE)
    }

    if (!par %in% par_alt_omit) {
      temp_prior_list <- c(temp_prior_list, list(prior_list_alt[[par]]))
      temp_is_null    <- c(temp_is_null, FALSE)
    }

    if (length(temp_prior_list) == 0) {
      stop(sprintf("At least one prior distribution needs to be defined for '%1$s' parameter in the '%2$s' model.", par, parameter))
    }

    prior_list[[par]] <- BayesTools::prior_mixture(
      prior_list = temp_prior_list,
      is_null    = temp_is_null
    )
  }

  prior_list <- .repair_formula_prior_list(
    prior_list = prior_list,
    parameter  = switch(parameter, "mods" = "mu", "scale" = "log_tau")
  )

  return(prior_list)
}


.warn_unmatched_prior_omissions <- function(omitted, available, argument) {

  if (is.null(omitted) || length(omitted) == 0L) {
    return(invisible(NULL))
  }

  omitted <- omitted[!is.na(omitted) & nzchar(omitted)]
  unknown <- setdiff(omitted, available)
  if (length(unknown) == 0L) {
    return(invisible(NULL))
  }

  warning(
    "Unknown term(s) in '", argument, "' were ignored: ",
    paste0("'", unknown, "'", collapse = ", "),
    ". Available terms are: ",
    paste0("'", available, "'", collapse = ", "),
    ".",
    call. = FALSE
  )

  invisible(NULL)
}


# Restore BayesTools formula metadata lost when formula prior mixtures are rebuilt.
.repair_formula_prior_list <- function(prior_list, parameter) {

  if (is.null(prior_list)) {
    return(prior_list)
  }

  for (i in seq_along(prior_list)) {
    prior_list[[i]] <- .repair_formula_prior(
      prior    = prior_list[[i]],
      parameter = parameter
    )
  }

  return(prior_list)
}
.repair_formula_prior <- function(prior, parameter) {

  attr(prior, "parameter") <- parameter

  formula_attributes <- c(
    "interaction", "interaction_terms", "levels", "level_names",
    "term_components", "factor_terms", "factor_contrasts",
    "factor_design", "factor_cell_names", "multiply_by"
  )

  if (.has_formula_prior_components(prior)) {
    source_attributes <- NULL
    for (i in seq_along(prior)) {
      if (!is.null(attr(prior[[i]], "parameter"))) {
        source_attributes <- attributes(prior[[i]])
        break
      }
    }

    if (!is.null(source_attributes)) {
      for (attribute in formula_attributes) {
        if (is.null(attr(prior, attribute)) && !is.null(source_attributes[[attribute]])) {
          attr(prior, attribute) <- source_attributes[[attribute]]
        }
      }
    }

    for (i in seq_along(prior)) {
      attr(prior[[i]], "parameter") <- parameter
      for (attribute in formula_attributes) {
        if (is.null(attr(prior[[i]], attribute)) && !is.null(attr(prior, attribute))) {
          attr(prior[[i]], attribute) <- attr(prior, attribute)
        }
      }
    }
  }

  return(prior)
}
.has_formula_prior_components <- function(prior) {

  return(BayesTools::is.prior.mixture(prior) ||
    BayesTools::is.prior.spike_and_slab(prior))
}
