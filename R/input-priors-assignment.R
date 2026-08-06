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

.assign_prior.location                 <- function(prior_effect, prior_mods, data) {

  if (!.is_data_random(data)) {
    return(NULL)
  }
  if (!is.null(prior_mods)) {
    return(prior_mods)
  }

  return(list(intercept = prior_effect))
}

.assign_prior.random                   <- function(prior, data, sd_sources = NULL) {

  random_effects <- attr(data[["location"]], "random_effects")
  if (!inherits(random_effects, "BayesTools_random_effects") ||
      length(random_effects[["terms"]]) == 0L) {
    stop("Internal error: missing parsed random-effect object.", call. = FALSE)
  }

  terms      <- random_effects[["terms"]]
  block_info <- .assign_prior.random_block_info(
    terms = terms,
    data  = data[["location"]]
  )
  allocation <- .assign_prior.random_allocation(
    terms           = terms,
    block_info      = block_info,
    sd              = prior,
    sd_sources      = sd_sources,
    allocation_name = "random_total"
  )
  block_sd_sources <- list()
  block_sds        <- list()
  if (is.null(allocation)) {
    if (!is.null(sd_sources)) {
      block_names <- vapply(block_info, `[[`, character(1), "block")
      groups <- .assign_prior.random_component_groups(
        terms       = terms,
        block_names = block_names
      )
      for (group in groups) {
        source <- sd_sources[[group[["label"]]]]
        if (is.null(source)) {
          for (block in group[["blocks"]]) {
            block_sds[[block]] <- prior
          }
          next
        }
        block <- group[["blocks"]][[1L]]
        block_sd_sources[[block]] <- source
      }
    } else {
      block <- block_info[[1L]][["block"]]
      block_sds[[block]] <- prior
    }
  } else if (!is.null(sd_sources)) {
    groups <- .assign_prior.random_component_groups(
      terms       = terms,
      block_names = vapply(block_info, `[[`, character(1), "block")
    )
    block_heterogeneous <- stats::setNames(
      vapply(block_info, `[[`, logical(1), "heterogeneous"),
      vapply(block_info, `[[`, character(1), "block")
    )
    for (group in groups) {
      source <- sd_sources[[group[["label"]]]]
      if (is.null(source)) {
        for (block in group[["blocks"]]) {
          block_sds[[block]] <- prior
        }
        next
      }
      block <- group[["blocks"]][[1L]]
      if (length(group[["blocks"]]) == 1L &&
          !isTRUE(block_heterogeneous[[block]])) {
        block_sd_sources[[block]] <- source
      }
    }
  }
  block_args <- .assign_prior.random_block_args(
    block_info = block_info,
    sds        = block_sds,
    sd_sources = block_sd_sources
  )

  args <- c(
    block_args,
    list(
      allocation = allocation
    )
  )
  args <- args[!vapply(args, is.null, logical(1))]

  do.call(BayesTools::prior_random, args)
}

.assign_prior.random_scale_sources     <- function(data) {

  random_effects <- attr(data[["location"]], "random_effects")
  terms          <- random_effects[["terms"]]
  block_names    <- vapply(terms, `[[`, character(1), "block_name")
  components     <- unique(.assign_prior.random_component_labels(
    terms       = terms,
    block_names = block_names
  ))
  scale_specs    <- .data_scale_component_specs(data)

  if (length(scale_specs) == 1L &&
      identical(scale_specs[[1L]][["source"]], "tau")) {
    sources <- list(BayesTools::random_sd_source("tau", shape = "row"))
    names(sources) <- components
    return(sources)
  }

  sources <- lapply(scale_specs, function(scale_spec) {
    BayesTools::random_sd_source(scale_spec[["source"]], shape = "row")
  })
  names(sources) <- names(scale_specs)

  sources
}

.validate_prior_random                 <- function(data, prior_location,
                                                   prior_random) {

  if (!.is_data_random(data)) {
    return(invisible(TRUE))
  }

  formula_design <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(data = data, parameter = "location"),
    parameter     = "mu",
    data          = data[["location"]],
    prior_list    = prior_location,
    formula_scale = .data_standardize_continuous_predictors(data),
    prior_random  = prior_random
  )[["formula_design"]]

  .validate_prior_random_row_sources(
    data           = data,
    formula_design = formula_design
  )
  .validate_prior_random_scale_sources(
    data           = data,
    formula_design = formula_design
  )

  invisible(TRUE)
}

.validate_prior_random_scale_sources   <- function(data, formula_design) {

  if (!.is_data_scale(data)) {
    return(invisible(TRUE))
  }

  random_effects <- attr(data[["location"]], "random_effects")
  terms          <- random_effects[["terms"]]
  block_names    <- vapply(terms, `[[`, character(1), "block_name")
  components     <- .assign_prior.random_component_labels(
    terms       = terms,
    block_names = block_names
  )
  component_by_block <- stats::setNames(components, block_names)
  expected_sources   <- .assign_prior.random_scale_source_names(data)

  for (term in formula_design[["random_effects"]]) {
    block     <- term[["block_name"]]
    component <- component_by_block[[block]]
    expected  <- if (component %in% names(expected_sources)) {
      expected_sources[[component]]
    } else {
      NULL
    }
    if (is.null(expected)) {
      next
    }
    if (!.validate_prior_random_scale_binding(term[["sd_binding"]], expected)) {
      stop(
        "'prior_heterogeneity = BayesTools::prior_random(...)' with 'scale' ",
        "must bind random component '", component, "' to row-shaped SD source ",
        "'", expected, "' via BayesTools::random_sd_source('", expected,
        "', shape = 'row').",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.validate_prior_random_row_sources     <- function(data, formula_design) {

  expected_sources <- if (.is_data_scale(data)) {
    .assign_prior.random_scale_source_names(data)
  } else {
    character(0)
  }

  for (term in formula_design[["random_effects"]]) {
    binding <- term[["sd_binding"]]
    if (is.null(binding) || is.null(binding[["source"]])) {
      next
    }
    source <- binding[["source"]]
    if (!inherits(source, "random_sd_source") ||
        !identical(source[["shape"]], "row")) {
      next
    }
    if (source[["name"]] %in% expected_sources) {
      next
    }

    stop(
      "Row-shaped SD source '", source[["name"]],
      "' in 'prior_heterogeneity = BayesTools::prior_random(...)' ",
      "is available only for random components targeted by 'scale'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.assign_prior.random_scale_source_names <- function(data) {

  sources <- .assign_prior.random_scale_sources(data)
  vapply(sources, `[[`, character(1), "name")
}

.validate_prior_random_scale_binding   <- function(binding, expected) {

  if (is.null(binding) || is.null(binding[["source"]])) {
    return(FALSE)
  }
  source <- binding[["source"]]

  inherits(source, "random_sd_source") &&
    identical(source[["name"]], expected) &&
    identical(source[["shape"]], "row") &&
    identical(source[["kind"]], "external") &&
    identical(source[["owned"]], FALSE)
}

.assign_prior.random_block_info        <- function(terms, data) {

  lapply(terms, function(term) {
    structure   <- .assign_prior.random_structure(term)
    n_columns   <- .assign_prior.random_n_columns(term, structure, data)
    homogeneous <- .assign_prior.random_homogeneous_sd(term, structure)

    list(
      block         = term[["block_name"]],
      structure     = structure,
      n_columns     = n_columns,
      heterogeneous = !homogeneous && n_columns > 1L,
      needs_cor     = structure == "us" && n_columns > 1L,
      needs_rho     = structure %in% c("cs", "hcs", "ar1", "car", "har") &&
        n_columns > 1L
    )
  })
}

.assign_prior.random_structure         <- function(term) {

  structure <- term[["structure"]]
  if (!is.character(structure) || length(structure) != 1L ||
      is.na(structure) || !nzchar(structure)) {
    stop("Internal error: parsed random-effect term is missing structure.", call. = FALSE)
  }

  tolower(structure)
}

.assign_prior.random_homogeneous_sd    <- function(term, structure) {

  switch(
    structure,
    id   = TRUE,
    diag = identical(term[["hom"]], TRUE),
    cs   = TRUE,
    hcs  = FALSE,
    ar1  = TRUE,
    car  = TRUE,
    har  = FALSE,
    us   = FALSE,
    FALSE
  )
}

.assign_prior.random_n_columns         <- function(term, structure, data) {

  if (structure %in% c("cs", "hcs", "ar1", "har")) {
    variables <- all.vars(term[["term_formula"]])
    values <- lapply(variables, function(variable) {
      if (!variable %in% names(data)) {
        stop("Random-effect variable '", variable, "' is missing.", call. = FALSE)
      }
      .assign_prior.random_index_component(data[[variable]])
    })
    if (length(values) == 1L) {
      return(nlevels(values[[1L]]))
    }
    return(nlevels(do.call(interaction, c(values, list(drop = TRUE)))))
  }
  if (structure == "car") {
    variables <- all.vars(term[["term_formula"]])
    if (length(variables) != 1L || !variables %in% names(data)) {
      stop(
        "CAR random-effect terms require exactly one time variable.",
        call. = FALSE
      )
    }
    return(length(unique(data[[variables]])))
  }

  model_matrix <- stats::model.matrix(term[["term_formula"]], data = data)
  ncol(model_matrix)
}

.assign_prior.random_index_component   <- function(x) {

  if (is.factor(x)) {
    return(x)
  }
  factor(x, levels = sort(unique(x)))
}

.assign_prior.random_block_args        <- function(block_info, sds = list(),
                                                   sd_sources = list()) {

  args <- list()
  for (info in block_info) {
    covariance <- NULL
    if (isTRUE(info[["needs_cor"]])) {
      covariance <- BayesTools::random_covariance(
        cor = BayesTools::prior_lkj(eta = 1)
      )
    } else if (isTRUE(info[["needs_rho"]])) {
      covariance <- BayesTools::random_covariance(
        rho       = BayesTools::prior("beta", list(alpha = 1, beta = 1)),
        rho_scale = "rho"
      )
    }
    sd        <- sds[[info[["block"]]]]
    sd_source <- sd_sources[[info[["block"]]]]
    if (!is.null(covariance) || !is.null(sd) || !is.null(sd_source)) {
      args[[info[["block"]]]] <- BayesTools::random_block(
        sd         = sd,
        covariance = covariance,
        sd_source  = sd_source
      )
    }
  }

  args
}

.assign_prior.random_allocation        <- function(terms, block_info, sd,
                                                   sd_sources = NULL,
                                                   allocation_name) {

  block_names   <- vapply(block_info, `[[`, character(1), "block")
  heterogeneous <- vapply(block_info, `[[`, logical(1), "heterogeneous")
  n_columns     <- vapply(block_info, `[[`, integer(1), "n_columns")

  if (!is.null(sd_sources)) {
    return(.assign_prior.random_component_allocations(
      terms       = terms,
      block_info  = block_info,
      sd_sources  = sd_sources
    ))
  }

  if (length(block_names) == 1L) {
    if (isTRUE(heterogeneous[[1L]])) {
      return(BayesTools::random_variance_allocation(
        terms      = block_names,
        sd         = sd,
        weights    = .assign_prior.random_dirichlet(n_columns[[1L]]),
        target     = "sd_component",
        scale      = "mean_variance"
      ))
    }
    return(NULL)
  }

  groups <- .assign_prior.random_component_groups(terms, block_names)
  if (length(groups) == 1L) {
    group <- groups[[1L]]
    allocation <- BayesTools::random_variance_allocation(
      name      = allocation_name,
      terms     = stats::setNames(group[["blocks"]], group[["child_labels"]]),
      sd        = sd,
      weights   = .assign_prior.random_dirichlet(length(group[["blocks"]]))
    )
    allocations <- list(allocation)
    parent_for_block <- list()
    for (i in seq_along(group[["blocks"]])) {
      parent_for_block[[group[["blocks"]][[i]]]] <- list(
        allocation = allocation_name,
        component  = group[["child_labels"]][[i]]
      )
    }
    for (i in seq_along(block_names)) {
      if (!isTRUE(heterogeneous[[i]])) {
        next
      }
      parent <- parent_for_block[[block_names[[i]]]]
      allocations[[length(allocations) + 1L]] <- BayesTools::random_variance_allocation(
        terms      = block_names[[i]],
        parent     = BayesTools::allocation_ref(
          parent[["allocation"]],
          parent[["component"]]
        ),
        weights    = .assign_prior.random_dirichlet(n_columns[[i]]),
        target     = "sd_component",
        scale      = "mean_variance"
      )
    }
    return(allocations)
  }

  root_terms <- .assign_prior.random_root_terms(groups)
  root_allocation <- BayesTools::random_variance_allocation(
    name      = allocation_name,
    terms     = root_terms,
    sd        = sd,
    weights   = .assign_prior.random_dirichlet(length(root_terms))
  )
  allocations <- list(root_allocation)

  parent_for_block <- .assign_prior.random_parent_map(
    groups          = groups,
    allocation_name = allocation_name
  )
  for (group in groups) {
    if (length(group[["blocks"]]) <= 1L) {
      next
    }
    child_name <- paste0(group[["label"]], "_split")
    allocations[[length(allocations) + 1L]] <- BayesTools::random_variance_allocation(
      name    = child_name,
      terms   = stats::setNames(group[["blocks"]], group[["child_labels"]]),
      parent  = BayesTools::allocation_ref(allocation_name, group[["label"]]),
      weights = .assign_prior.random_dirichlet(length(group[["blocks"]]))
    )
  }

  for (i in seq_along(block_names)) {
    if (!isTRUE(heterogeneous[[i]])) {
      next
    }
    parent <- parent_for_block[[block_names[[i]]]]
    allocations[[length(allocations) + 1L]] <- BayesTools::random_variance_allocation(
      terms      = block_names[[i]],
      parent     = BayesTools::allocation_ref(
        parent[["allocation"]],
        parent[["component"]]
      ),
      weights    = .assign_prior.random_dirichlet(n_columns[[i]]),
      target     = "sd_component",
      scale      = "mean_variance"
    )
  }

  allocations
}

.assign_prior.random_component_allocations <- function(terms, block_info,
                                                       sd_sources) {

  block_names   <- vapply(block_info, `[[`, character(1), "block")
  heterogeneous <- vapply(block_info, `[[`, logical(1), "heterogeneous")
  n_columns     <- vapply(block_info, `[[`, integer(1), "n_columns")
  groups        <- .assign_prior.random_component_groups(terms, block_names)
  allocations   <- list()

  for (group in groups) {
    source <- sd_sources[[group[["label"]]]]
    if (is.null(source)) {
      next
    }

    if (length(group[["blocks"]]) > 1L) {
      allocation_name <- .assign_prior.random_component_allocation_name(
        group  = group,
        groups = groups
      )
      allocations[[allocation_name]] <- BayesTools::random_variance_allocation(
        name      = allocation_name,
        terms     = stats::setNames(group[["blocks"]], group[["child_labels"]]),
        sd_source = source,
        weights   = .assign_prior.random_dirichlet(length(group[["blocks"]]))
      )
    }

    for (block in group[["blocks"]]) {
      index <- match(block, block_names)
      if (!isTRUE(heterogeneous[[index]])) {
        next
      }
      if (length(group[["blocks"]]) > 1L) {
        allocation_name <- .assign_prior.random_component_allocation_name(
          group  = group,
          groups = groups
        )
        parent <- BayesTools::allocation_ref(
          allocation_name,
          group[["child_labels"]][[match(block, group[["blocks"]])]]
        )
        allocations[[length(allocations) + 1L]] <- BayesTools::random_variance_allocation(
          terms      = block,
          parent     = parent,
          weights    = .assign_prior.random_dirichlet(n_columns[[index]]),
          target     = "sd_component",
          scale      = "mean_variance"
        )
      } else {
        allocations[[length(allocations) + 1L]] <- BayesTools::random_variance_allocation(
          terms      = block,
          sd_source  = source,
          weights    = .assign_prior.random_dirichlet(n_columns[[index]]),
          target     = "sd_component",
          scale      = "mean_variance"
        )
      }
    }
  }

  if (length(allocations) == 0L) {
    return(NULL)
  }

  unname(allocations)
}

.assign_prior.random_component_allocation_name <- function(group, groups) {

  if (length(groups) == 1L) {
    return("total")
  }

  group[["label"]]
}

.assign_prior.random_dirichlet         <- function(n) {

  BayesTools::prior("dirichlet", parameters = list(alpha = rep(1, n)))
}

.assign_prior.random_component_groups  <- function(terms, block_names) {

  component_labels <- .assign_prior.random_component_labels(
    terms       = terms,
    block_names = block_names
  )
  child_labels <- vapply(seq_along(terms), function(i) {
    label <- terms[[i]][["component_child_label"]]
    if (is.null(label) || length(label) != 1L ||
        is.na(label) || !nzchar(label)) {
      return(block_names[[i]])
    }
    label
  }, character(1))

  unique_components <- unique(component_labels)
  lapply(unique_components, function(component) {
    index <- which(component_labels == component)
    list(
      label        = component,
      blocks       = block_names[index],
      child_labels = child_labels[index]
    )
  })
}

.assign_prior.random_component_labels  <- function(terms, block_names) {

  vapply(seq_along(terms), function(i) {
    label <- terms[[i]][["component_label"]]
    if (is.null(label) || length(label) != 1L ||
        is.na(label) || !nzchar(label)) {
      return(block_names[[i]])
    }
    label
  }, character(1))
}

.assign_prior.random_root_terms        <- function(groups) {

  terms <- vapply(groups, function(group) {
    if (length(group[["blocks"]]) == 1L) {
      return(group[["blocks"]][[1L]])
    }
    group[["label"]]
  }, character(1))
  stats::setNames(terms, vapply(groups, `[[`, character(1), "label"))
}

.assign_prior.random_parent_map        <- function(groups, allocation_name) {

  out <- list()
  for (group in groups) {
    if (length(group[["blocks"]]) == 1L) {
      out[[group[["blocks"]][[1L]]]] <- list(
        allocation = allocation_name,
        component  = group[["label"]]
      )
      next
    }

    child_name <- paste0(group[["label"]], "_split")
    for (i in seq_along(group[["blocks"]])) {
      out[[group[["blocks"]][[i]]]] <- list(
        allocation = child_name,
        component  = group[["child_labels"]][[i]]
      )
    }
  }

  out
}

.assign_prior.baserate                 <- function(prior) {

  # use default settings if prior distribution is not specified
  if (missing(prior) || is.null(prior)) {
    prior <- BayesTools::prior_factor("beta", parameters = list("alpha" = 1, "beta" = 1), contrast = "independent")
    return(prior)
  }

  # check the user specified prior distribution
  prior <- .check_prior.simple_or_point(
    prior,
    prior_name = "prior_baserate"
  )
  if (BayesTools::is.prior.discrete(prior) &&
      !BayesTools::is.prior.point(prior)) {
    stop(
      "Non-point discrete priors are not supported for 'prior_baserate'; use ",
      "a point or beta prior.",
      call. = FALSE
    )
  }
  if (!BayesTools::is.prior.point(prior) &&
      !identical(prior[["distribution"]], "beta")) {
    stop(
      "Only point and beta priors are supported for 'prior_baserate' because ",
      "post-fit likelihood diagnostics require certified nuisance integration; ",
      "received '", prior[["distribution"]], "'.",
      call. = FALSE
    )
  }

  prior <- .check_prior.restricted_01(
    prior,
    prior_name = "prior_baserate"
  )
  if (BayesTools::is.prior.point(prior) &&
      (mean(prior) <= 0 || mean(prior) >= 1)) {
    stop(
      "Point prior distribution for 'prior_baserate' must be strictly within (0, 1).",
      call. = FALSE
    )
  }

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
  if (BayesTools::is.prior.discrete(prior) &&
      !BayesTools::is.prior.point(prior)) {
    stop(
      "Non-point discrete priors are not supported for 'prior_lograte'; use ",
      "a point or normal prior.",
      call. = FALSE
    )
  }
  if (!BayesTools::is.prior.point(prior) &&
      !identical(prior[["distribution"]], "normal")) {
    stop(
      "Only point and normal priors are supported for 'prior_lograte' because ",
      "post-fit likelihood diagnostics require certified nuisance integration; ",
      "received '", prior[["distribution"]], "'.",
      call. = FALSE
    )
  }

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
