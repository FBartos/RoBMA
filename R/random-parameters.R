# Internal semantic random-parameter extraction.

.brma_random_parameter_bundle <- function(
    object, standardized_coefficients = FALSE, chains = FALSE,
    prior = FALSE, n_prior_samples = 10000L, seed = NULL) {

  fit <- object[["fit"]]
  if (!.is_random(object) || is.null(fit)) {
    return(list(
      samples = matrix(numeric(), nrow = 0L, ncol = 0L),
      specs   = .brma_random_parameter_empty_specs(),
      priors  = list()
    ))
  }

  if (prior) {
    raw_samples <- BayesTools::transform_prior_samples(
      fit           = fit,
      n_samples     = n_prior_samples,
      seed          = seed,
      formula_scale = list()
    )
    extracted <- .brma_random_parameter_extract_fit(
      fit                       = .brma_random_parameter_fit_with_samples(
        fit,
        coda::mcmc.list(coda::mcmc(raw_samples))
      ),
      standardized_coefficients = standardized_coefficients
    )
    return(extracted)
  }

  if (!chains) {
    return(.brma_random_parameter_extract_fit(
      fit                       = fit,
      standardized_coefficients = standardized_coefficients
    ))
  }

  raw_chains <- coda::as.mcmc.list(fit)
  extracted  <- lapply(raw_chains, function(chain) {
    .brma_random_parameter_extract_fit(
      fit = .brma_random_parameter_fit_with_samples(
        fit,
        coda::mcmc.list(chain)
      ),
      standardized_coefficients = standardized_coefficients
    )
  })
  specs <- extracted[[1L]][["specs"]]
  if (any(vapply(extracted, function(x) {
      !identical(x[["specs"]][["parameter"]], specs[["parameter"]])
    }, logical(1)))) {
    stop("Random-parameter columns differ across MCMC chains.", call. = FALSE)
  }

  semantic_chains <- lapply(seq_along(raw_chains), function(i) {
    chain <- raw_chains[[i]]
    coda::mcmc(
      extracted[[i]][["samples"]],
      start = stats::start(chain),
      end   = stats::end(chain),
      thin  = coda::thin(chain)
    )
  })

  list(
    samples = coda::mcmc.list(semantic_chains),
    specs   = specs,
    priors  = extracted[[1L]][["priors"]]
  )
}

.brma_random_parameter_extract_fit <- function(
    fit, standardized_coefficients = FALSE) {

  samples <- BayesTools::JAGS_estimates_table(
    fit                     = fit,
    keep_parameters         = "random_effects",
    random_effects_summary  = "standard",
    random_effects_metadata = TRUE,
    remove_inclusion        = TRUE,
    remove_diagnostics      = TRUE,
    return_samples          = TRUE,
    transform_scaled        = !standardized_coefficients
  )
  samples <- as.matrix(samples)
  priors  <- attr(samples, "prior_list", exact = TRUE)
  if (is.null(priors) || length(priors) == 0L || ncol(samples) == 0L) {
    return(list(
      samples = matrix(numeric(), nrow = nrow(samples), ncol = 0L),
      specs   = .brma_random_parameter_empty_specs(),
      priors  = list()
    ))
  }

  summary_type <- vapply(priors, function(prior) {
    value <- attr(prior, "random_summary", exact = TRUE)
    if (is.null(value)) NA_character_ else value
  }, character(1))
  keep <- summary_type %in% c(
    "sd", "sd_total", "rho", "cor", "var_frac", "var_ratio",
    "sd_multiplier"
  )
  priors <- priors[keep]
  if (length(priors) == 0L) {
    return(list(
      samples = matrix(numeric(), nrow = nrow(samples), ncol = 0L),
      specs   = .brma_random_parameter_empty_specs(),
      priors  = list()
    ))
  }

  labels <- vapply(priors, function(prior) {
    attr(prior, "random_summary_label", exact = TRUE)
  }, character(1))
  formula_parameters <- vapply(priors, function(prior) {
    attr(prior, "parameter", exact = TRUE)
  }, character(1))
  parameters <- paste0("(", formula_parameters, ") ", labels)
  columns    <- match(parameters, colnames(samples))
  if (anyNA(columns)) {
    missing <- parameters[is.na(columns)]
    stop(
      "Semantic random-effect draws are missing expected columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  samples             <- samples[, columns, drop = FALSE]
  names(priors)       <- parameters
  colnames(samples) <- parameters
  specs             <- .brma_random_parameter_specs(priors)
  specs[["source_parameter"]] <- .brma_random_parameter_source_names(
    fit            = fit,
    specs          = specs,
    summary_priors = priors
  )

  list(samples = samples, specs = specs, priors = priors)
}

.brma_random_parameter_specs <- function(priors) {

  attr_string <- function(prior, name) {
    value <- attr(prior, name, exact = TRUE)
    if (is.null(value) || length(value) == 0L) NA_character_ else as.character(value[[1L]])
  }
  attr_string_first <- function(prior, names) {
    values <- vapply(names, function(name) attr_string(prior, name), character(1))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 0L) NA_character_ else values[[1L]]
  }

  data.frame(
    parameter         = names(priors),
    label             = vapply(priors, attr_string, character(1), "random_summary_label"),
    summary_type      = vapply(priors, attr_string, character(1), "random_summary"),
    formula_parameter = vapply(priors, attr_string, character(1), "parameter"),
    block             = vapply(
      priors,
      attr_string_first,
      character(1),
      names = c("random_factor", "random_block")
    ),
    grouping          = vapply(priors, attr_string, character(1), "random_grouping_factor"),
    structure         = vapply(priors, attr_string, character(1), "random_structure"),
    allocation        = vapply(priors, attr_string, character(1), "random_allocation"),
    random_component  = vapply(priors, attr_string, character(1), "random_component"),
    source_parameter  = NA_character_,
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
}

.brma_random_parameter_empty_specs <- function() {

  data.frame(
    parameter         = character(),
    label             = character(),
    summary_type      = character(),
    formula_parameter = character(),
    block             = character(),
    grouping          = character(),
    structure         = character(),
    allocation        = character(),
    random_component  = character(),
    source_parameter  = character(),
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
}

.brma_random_parameter_source_names <- function(
    fit, specs, summary_priors) {

  prior_list     <- attr(fit, "prior_list", exact = TRUE)
  formula_design <- attr(fit, "formula_design", exact = TRUE)
  prior_names    <- names(prior_list)
  out            <- rep(NA_character_, nrow(specs))

  if (length(prior_names) == 0L || nrow(specs) == 0L) {
    return(out)
  }

  for (i in seq_len(nrow(specs))) {
    spec          <- as.list(specs[i, , drop = FALSE])
    summary_prior <- summary_priors[[spec[["parameter"]]]]
    source        <- .brma_random_parameter_source_name(
      spec           = spec,
      summary_prior  = summary_prior,
      prior_list     = prior_list,
      formula_design = formula_design
    )
    if (length(source) == 1L && !is.na(source) && source %in% prior_names) {
      out[i] <- source
    }
  }

  out
}


.brma_random_parameter_source_name <- function(
    spec, summary_prior, prior_list, formula_design) {

  type <- spec[["summary_type"]]
  if (type %in% c("var_frac", "var_ratio", "sd_multiplier")) {
    metadata <- attr(
      summary_prior,
      "random_allocation_metadata",
      exact = TRUE
    )
    return(.brma_random_parameter_unique_name(metadata[["weight_name"]]))
  }

  if (identical(type, "sd_total")) {
    matches <- vapply(prior_list, function(prior) {
      isTRUE(attr(prior, "random_sd_total", exact = TRUE)) &&
        .brma_random_parameter_metadata_matches(
          attr(prior, "parameter", exact = TRUE),
          spec[["formula_parameter"]]
        ) &&
        .brma_random_parameter_metadata_matches(
          attr(prior, "random_allocation", exact = TRUE),
          spec[["allocation"]]
        )
    }, logical(1))
    return(.brma_random_parameter_unique_name(names(prior_list)[matches]))
  }

  term <- .brma_random_parameter_design_term(formula_design, spec)
  if (is.null(term)) {
    return(NA_character_)
  }

  if (identical(type, "sd")) {
    return(.brma_random_parameter_sd_source(term, spec))
  }
  if (identical(type, "rho")) {
    correlation <- term[["correlation"]]
    if (!is.list(correlation) || !identical(correlation[["type"]], "rho")) {
      return(NA_character_)
    }
    if (!is.null(correlation[["sample_fixed"]])) {
      return(NA_character_)
    }
    rho_scale <- correlation[["rho_scale"]]
    if (!is.null(rho_scale) && !identical(rho_scale, "rho")) {
      return(NA_character_)
    }
    candidates <- c(correlation[["rho_name"]], correlation[["sample_name"]])
    candidates <- intersect(unique(candidates), names(prior_list))
    return(.brma_random_parameter_unique_name(candidates))
  }

  NA_character_
}


.brma_random_parameter_design_term <- function(formula_design, spec) {

  if (!is.list(formula_design)) {
    return(NULL)
  }
  designs <- Filter(function(design) {
    .brma_random_parameter_metadata_matches(
      design[["parameter"]],
      spec[["formula_parameter"]]
    )
  }, formula_design)
  terms <- unlist(lapply(designs, `[[`, "random_effects"), recursive = FALSE)
  terms <- Filter(function(term) {
    .brma_random_parameter_metadata_matches(
      term[["block_name"]],
      spec[["block"]]
    ) && .brma_random_parameter_metadata_matches(
      term[["group_label"]],
      spec[["grouping"]]
    )
  }, terms)

  if (length(terms) == 1L) terms[[1L]] else NULL
}


.brma_random_parameter_sd_source <- function(term, spec) {

  sources <- unique(term[["sd_parameter_names"]])
  sources <- sources[!is.na(sources) & nzchar(sources)]
  if (length(sources) == 1L) {
    return(sources)
  }
  if (length(sources) == 0L) {
    return(NA_character_)
  }

  leaves <- term[["sd_leaves"]]
  terms  <- leaves[["leaf_terms"]]
  if (is.null(terms) || is.null(names(terms))) {
    return(NA_character_)
  }
  components <- unname(terms[sources])
  components <- .brma_random_parameter_normalize_components(components, term)
  matches    <- !is.na(components) & components == spec[["random_component"]]
  .brma_random_parameter_unique_name(sources[matches])
}


.brma_random_parameter_normalize_components <- function(components, term) {

  components <- gsub("__xXx__", ":", components, fixed = TRUE)
  components[components == "sd"]          <- "shared"
  components[components == "(Intercept)"] <- "intercept"

  index <- term[["structured_index"]]
  if (is.list(index) && length(index[["name"]]) == 1L &&
      length(index[["label"]]) == 1L &&
      !identical(index[["name"]], index[["label"]])) {
    replace <- components == index[["name"]] |
      startsWith(components, paste0(index[["name"]], "["))
    components[replace] <- paste0(
      index[["label"]],
      substr(
        components[replace],
        nchar(index[["name"]]) + 1L,
        nchar(components[replace])
      )
    )
  }

  components
}


.brma_random_parameter_metadata_matches <- function(value, target) {

  missing_value  <- is.null(value) || length(value) != 1L || is.na(value)
  missing_target <- is.null(target) || length(target) != 1L || is.na(target)
  if (missing_value || missing_target) {
    return(missing_value && missing_target)
  }
  identical(as.character(value), as.character(target))
}


.brma_random_parameter_unique_name <- function(names) {

  names <- unique(names)
  names <- names[!is.na(names) & nzchar(names)]
  if (length(names) == 1L) names else NA_character_
}

.brma_random_parameter_fit_with_samples <- function(fit, samples) {

  fit[["mcmc"]] <- samples
  fit
}

.brma_random_parameter_diagnostic_fit <- function(
    object, parameter, standardized_coefficients = FALSE) {

  selected <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients,
    chains                    = TRUE
  )
  support <- .brma_random_parameter_support(
    selected[["spec"]],
    selected[["prior"]],
    selected[["source_prior"]]
  )
  diagnostic_prior <- BayesTools::prior(
    distribution = "normal",
    parameters   = list(mean = 0, sd = 1),
    truncation   = list(lower = support[1L], upper = support[2L])
  )
  fit <- .brma_random_parameter_fit_with_samples(
    object[["fit"]],
    selected[["samples"]]
  )
  attr(fit, "prior_list") <- stats::setNames(
    list(diagnostic_prior),
    selected[["entry"]][["parameter"]]
  )

  list(
    fit       = fit,
    parameter = selected[["entry"]][["parameter"]],
    label     = selected[["entry"]][["term"]]
  )
}

.brma_random_parameter_select <- function(
    object, parameter, standardized_coefficients = FALSE,
    chains = FALSE, prior = FALSE, n_prior_samples = 10000L,
    seed = NULL) {

  entry <- .brma_parameter_select_entry(
    object    = object,
    parameter = parameter,
    component = "random"
  )
  bundle <- .brma_random_parameter_bundle(
    object                    = object,
    standardized_coefficients = standardized_coefficients,
    chains                    = chains,
    prior                     = prior,
    n_prior_samples           = n_prior_samples,
    seed                      = seed
  )
  index <- match(entry[["parameter"]], bundle[["specs"]][["parameter"]])
  if (is.na(index)) {
    stop(
      "Random-effect quantity '", entry[["parameter"]],
      "' is no longer available in the fitted draws.",
      call. = FALSE
    )
  }

  list(
    entry        = entry,
    spec         = as.list(bundle[["specs"]][index, , drop = FALSE]),
    samples      = if (chains) {
      bundle[["samples"]][, entry[["parameter"]], drop = FALSE]
    } else {
      bundle[["samples"]][, entry[["parameter"]], drop = FALSE]
    },
    prior        = bundle[["priors"]][[entry[["parameter"]]]],
    source_prior = .brma_random_parameter_source_prior(
      object,
      bundle[["specs"]][["source_parameter"]][index]
    )
  )
}

.brma_random_parameter_source_prior <- function(object, source_parameter) {

  if (is.na(source_parameter) || !nzchar(source_parameter)) {
    return(NULL)
  }
  attr(object[["fit"]], "prior_list", exact = TRUE)[[source_parameter]]
}

.brma_random_parameter_support <- function(spec, prior = NULL,
                                           source_prior = NULL) {

  type <- spec[["summary_type"]]
  support <- if (type %in% c("sd", "sd_total", "sd_multiplier")) {
    c(0, Inf)
  } else if (type %in% c("rho", "cor")) {
    c(-1, 1)
  } else if (identical(type, "var_frac")) {
    c(0, 1)
  } else if (identical(type, "var_ratio")) {
    metadata <- if (is.null(prior)) NULL else
      attr(prior, "random_allocation_metadata", exact = TRUE)
    upper <- if (is.null(metadata[["n_targets"]])) Inf else
      as.numeric(metadata[["n_targets"]])
    c(0, upper)
  } else {
    c(-Inf, Inf)
  }

  if (!is.null(source_prior) && !is.null(source_prior[["truncation"]])) {
    truncation <- source_prior[["truncation"]]
    lower      <- truncation[["lower"]]
    upper      <- truncation[["upper"]]
    if (length(lower) == 1L && length(upper) == 1L) {
      support <- c(
        max(support[1L], as.numeric(lower)),
        min(support[2L], as.numeric(upper))
      )
    }
  }

  support
}

.brma_random_parameter_point_test_reason <- function(spec) {

  type   <- spec[["summary_type"]]
  source <- spec[["source_parameter"]]
  label  <- spec[["label"]]
  if (type %in% c("cor", "rho") &&
      (is.na(source) || !nzchar(source))) {
    return(paste0(
      "Point-null Bayes factors are not available for derived pairwise ",
      "correlation '", label, "'."
    ))
  }
  if (identical(type, "sd") &&
      (is.na(source) || !nzchar(source))) {
    return(paste0(
      "Point-null Bayes factors are not available for derived component SD '",
      label, "'."
    ))
  }

  ""
}

.brma_random_parameter_prior_density <- function(samples, support,
                                                 n_points = 4096L) {

  samples <- as.numeric(samples)
  samples <- samples[is.finite(samples)]
  if (length(unique(samples)) < 2L) {
    return(NULL)
  }

  bandwidth <- stats::bw.nrd0(samples)
  reflected <- samples
  if (is.finite(support[1L])) {
    reflected <- c(reflected, 2 * support[1L] - samples)
  }
  if (is.finite(support[2L])) {
    reflected <- c(reflected, 2 * support[2L] - samples)
  }
  args <- list(x = reflected, n = n_points, bw = bandwidth)
  if (is.finite(support[1L])) args[["from"]] <- support[1L]
  if (is.finite(support[2L])) args[["to"]]   <- support[2L]
  density <- do.call(stats::density, args)
  integral <- sum(diff(density[["x"]]) *
    (density[["y"]][-1L] + density[["y"]][-length(density[["y"]])]) / 2)
  if (!is.finite(integral) || integral <= 0) {
    return(NULL)
  }
  density[["y"]] <- density[["y"]] / integral
  out <- list(
    density = list(x = density[["x"]], y = density[["y"]], mass = 1),
    points  = data.frame(x = numeric(), p = numeric()),
    n_grid  = n_points
  )
  class(out) <- c("prior_linear_density", "prior_density")
  attr(out, "support") <- support
  out
}

.brma_random_parameter_mixed_posterior <- function(
    object, parameter, standardized_coefficients = FALSE,
    prior = FALSE, n_prior_samples = 10000L, seed = NULL) {

  selected <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients
  )
  values <- unname(as.numeric(selected[["samples"]][, 1L]))
  attr(values, "sample_ind") <- FALSE
  attr(values, "models_ind") <- rep(1, length(values))
  attr(values, "parameter")  <- selected[["entry"]][["parameter"]]
  attr(values, "prior_list") <- BayesTools::prior_none()
  support <- .brma_random_parameter_support(
    selected[["spec"]],
    selected[["prior"]],
    selected[["source_prior"]]
  )
  attr(values, "posterior_support") <- structure(
    list(
      bounds = support,
      points = numeric(),
      exact  = TRUE,
      source = "model",
      type   = "interval"
    ),
    class = c("BayesTools_posterior_support", "list")
  )

  if (prior) {
    prior_selected <- .brma_random_parameter_select(
      object                    = object,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients,
      prior                     = TRUE,
      n_prior_samples           = n_prior_samples,
      seed                      = seed
    )
    prior_density <- .brma_random_parameter_prior_density(
      prior_selected[["samples"]][, 1L],
      support = support
    )
    if (is.null(prior_density)) {
      stop(
        "A prior-density overlay is not available for fixed random-effect ",
        "quantity '", selected[["entry"]][["term"]], "'.",
        call. = FALSE
      )
    }
    attr(values, "prior_density") <- prior_density
  }

  class(values) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "marginal_posterior.simple",
    "marginal_posterior"
  )
  out <- list(values)
  names(out) <- selected[["entry"]][["parameter"]]
  attr(out, "prior_list") <- stats::setNames(
    list(BayesTools::prior_none()),
    selected[["entry"]][["parameter"]]
  )
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}
