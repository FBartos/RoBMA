# Internal semantic random-parameter extraction.

.brma_random_parameter_io_quantity_map <- function() {

  c(
    sd         = "tau",
    var        = "tau2",
    sd_total   = "tau_total",
    var_total  = "tau2_total",
    sd_common  = "tau_common",
    var_common = "tau2_common",
    cor        = "rho",
    var_prop   = "tau2_prop",
    sd_mult    = "tau_mult",
    var_mult   = "tau2_mult"
  )
}


.brma_random_parameter_io_quantity <- function(quantity) {

  map <- .brma_random_parameter_io_quantity_map()
  out <- unname(map[quantity])
  out[is.na(out)] <- quantity[is.na(out)]
  out
}


# Label parts of random-effect catalog quantities under RoBMA's quantity names
# (tau for sd, rho for cor, ...).
.brma_random_parameter_io_parts <- function(label_parts) {

  lapply(label_parts, function(parts) {
    if (is.null(parts) || is.null(parts[["random"]])) {
      stop(
        "Random-effect catalog quantities have no random-effect label parts. ",
        "Refit the model with the current RoBMA/BayesTools build.",
        call. = FALSE
      )
    }
    parts[["random"]][["quantity"]] <- .brma_random_parameter_io_quantity(
      parts[["random"]][["quantity"]]
    )
    parts
  })
}


# RoBMA names of random-effect catalog quantities, rendered from their label
# parts: the selector ('(mu) study: tau(intercept)') or the label, displayed
# without the formula prefix ('study: tau').
.brma_random_parameter_io_labels <- function(quantities,
                                             type = c("selector", "label")) {

  type <- match.arg(type)
  if (nrow(quantities) == 0L) {
    return(character())
  }

  BayesTools::parameter_labels(
    .brma_random_parameter_io_parts(quantities[["label_parts"]]),
    style          = "table",
    formula_prefix = identical(type, "selector"),
    simplify       = identical(type, "label")
  )
}


# Aliases of one random-effect catalog quantity under RoBMA's quantity names,
# rendered from its label parts in the forms of the BayesTools catalog
# aliases: without the formula prefix, and simplified with and without the
# prefix and the owner (simplified aliases require simplify_names = TRUE).
.brma_random_parameter_io_aliases <- function(quantity) {

  parts <- .brma_random_parameter_io_parts(quantity[["label_parts"]])[[1L]]
  without_owner <- parts
  without_owner[["random"]][["owner"]] <- ""
  render <- function(parts, formula_prefix, simplify) {
    BayesTools::parameter_labels(
      parts,
      style          = "table",
      formula_prefix = formula_prefix,
      simplify       = simplify
    )
  }

  data.frame(
    alias      = c(
      render(parts, FALSE, FALSE),
      render(parts, TRUE, TRUE),
      render(parts, FALSE, TRUE),
      render(without_owner, FALSE, TRUE)
    ),
    simplified = c(FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

.brma_random_parameter_supported_quantities <- function() {

  c(
    "sd", "var", "sd_total", "var_total", "sd_common", "var_common",
    "cor", "var_prop", "var_mult", "sd_mult"
  )
}

.brma_random_parameter_bundle <- function(
    object, standardized_coefficients = FALSE, chains = FALSE,
    prior = FALSE, n_prior_samples = 10000L, seed = NULL,
    selections = NULL) {

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
      standardized_coefficients = standardized_coefficients,
      selections                = selections
    )
    extracted[["raw_samples"]] <- raw_samples
    return(extracted)
  }

  if (!chains) {
    return(.brma_random_parameter_extract_fit(
      fit                       = fit,
      standardized_coefficients = standardized_coefficients,
      selections                = selections
    ))
  }

  raw_chains <- coda::as.mcmc.list(fit)
  extracted  <- lapply(raw_chains, function(chain) {
    .brma_random_parameter_extract_fit(
      fit = .brma_random_parameter_fit_with_samples(
        fit,
        coda::mcmc.list(chain)
      ),
      standardized_coefficients = standardized_coefficients,
      selections                = selections
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

  samples <- coda::mcmc.list(semantic_chains)
  BayesTools::posterior_metadata(samples, "undefined_draws") <-
    BayesTools::posterior_metadata(extracted[[1L]][["samples"]], "undefined_draws")

  list(
    samples = samples,
    specs   = specs,
    priors  = extracted[[1L]][["priors"]]
  )
}

# Draws of a quantity declared as possibly undefined (an original-scale
# correlation with a zero SD; 'undefined_draws' metadata of the extracted
# samples) are left out where the quantity is undefined. Any other missing
# draw is an error. Returns the mask of defined draws.
.brma_random_parameter_defined_draws <- function(values, samples, label) {

  defined <- !is.na(values)
  if (all(defined)) {
    return(defined)
  }
  if (is.null(BayesTools::posterior_metadata(samples, "undefined_draws"))) {
    stop(
      "The draws of random-effect quantity '", label, "' contain missing ",
      "values. Missing draws are accepted only for quantities declared as ",
      "possibly undefined (original-scale random-effect correlations).",
      call. = FALSE
    )
  }
  if (!any(defined)) {
    stop(
      "Random-effect quantity '", label, "' is undefined in every draw.",
      call. = FALSE
    )
  }

  return(defined)
}

# Footnote for results computed over the defined draws only.
.brma_random_parameter_defined_footnote <- function(label, samples,
                                                    posterior_defined,
                                                    prior_defined = NULL) {

  counts <- function(defined, what) {
    if (is.null(defined) || all(defined)) {
      return(NULL)
    }
    paste(sum(defined), "of", length(defined), what, "draws")
  }
  parts <- c(
    counts(posterior_defined, "posterior"),
    counts(prior_defined, "prior")
  )
  if (length(parts) == 0L) {
    return(NULL)
  }
  reason <- BayesTools::posterior_metadata(samples, "undefined_draws")
  condition <- if (identical(unname(reason[[1L]]), "correlation")) {
    "where the correlation is defined, i.e. both SDs are positive."
  } else {
    "where the quantity is defined."
  }

  paste0(
    label, ": computed from ", paste(parts, collapse = " and "), " ",
    condition
  )
}

.brma_random_parameter_extract_fit <- function(
    fit, standardized_coefficients = FALSE, selections = NULL) {

  if (is.null(selections)) {
    catalog    <- BayesTools::parameter_catalog(fit)
    quantities <- catalog[["quantities"]]
    supported  <- .brma_random_parameter_supported_quantities()
    keep <- startsWith(quantities[["role"]], "random_") &
      !quantities[["internal"]] &
      quantities[["status"]] != "unavailable" &
      quantities[["quantity"]] %in% supported
    quantities <- quantities[keep, , drop = FALSE]
    selections <- lapply(seq_len(nrow(quantities)), function(i) {
      BayesTools::parameter_catalog_resolve(
        catalog,
        alias     = quantities[["canonical_name"]][i],
        namespace = quantities[["namespace"]][i]
      )
    })
  } else {
    valid <- is.list(selections) && length(selections) > 0L &&
      all(vapply(
        selections,
        inherits,
        logical(1),
        what = "BayesTools_parameter_selection"
      ))
    if (!valid) {
      stop("Random-parameter selections are invalid.", call. = FALSE)
    }
    quantities <- do.call(rbind, lapply(
      selections,
      `[[`,
      "quantities"
    ))
  }
  n_draws <- sum(vapply(coda::as.mcmc.list(fit), nrow, integer(1)))
  if (nrow(quantities) == 0L) {
    return(list(
      samples = matrix(numeric(), nrow = n_draws, ncol = 0L),
      specs   = .brma_random_parameter_empty_specs(),
      priors  = list()
    ))
  }

  extraction_fit <- fit
  if (standardized_coefficients) {
    attr(extraction_fit, "formula_scale") <- list()
  }
  model_samples <- NULL
  if (length(selections) > 1L) {
    dependencies <- unique(unlist(lapply(
      quantities[["extraction_key"]],
      `[[`,
      "dependencies"
    ), use.names = FALSE))
    model_samples <- if (length(dependencies) > 0L) {
      as.matrix(BayesTools::JAGS_materialize_draws(
        extraction_fit,
        parameters       = dependencies,
        include_internal = TRUE
      ))
    } else {
      matrix(
        numeric(),
        nrow = n_draws,
        ncol = 0L,
        dimnames = list(NULL, character())
      )
    }
  }
  draws <- lapply(selections, function(selection) {
    BayesTools::parameter_draws(
      extraction_fit,
      selection,
      model_samples = model_samples
    )
  })
  samples <- do.call(cbind, lapply(draws, as.matrix))
  parameter_names <- .brma_random_parameter_io_labels(quantities, "selector")
  colnames(samples) <- parameter_names
  # Quantities that can be undefined in some draws (original-scale
  # correlations with a zero SD) keep parameter_draws()' declaration, by
  # column, so that consumers accept only these missing draws.
  undefined <- unlist(lapply(seq_along(draws), function(i) {
    declared <- BayesTools::posterior_metadata(draws[[i]], "undefined_draws")
    if (is.null(declared)) {
      return(NULL)
    }
    stats::setNames(unname(declared[[1L]]), parameter_names[[i]])
  }))
  if (length(undefined) > 0L) {
    BayesTools::posterior_metadata(samples, "undefined_draws") <- undefined
  }
  specs  <- .brma_random_parameter_specs(quantities)
  specs[["display_transform"]] <- I(lapply(selections, function(selection) {
    BayesTools::parameter_transform(extraction_fit, selection)
  }))
  priors <- stats::setNames(
    rep(list(NULL), nrow(quantities)),
    parameter_names
  )

  list(samples = samples, specs = specs, priors = priors)
}

.brma_random_parameter_specs <- function(quantities) {

  keys <- quantities[["extraction_key"]]
  key_string <- function(key, field) {
    value <- key[[field]]
    if (is.null(value) || length(value) != 1L || is.na(value)) "" else
      as.character(value)
  }
  key_number <- function(key, field) {
    value <- key[[field]]
    if (is.null(value) || length(value) != 1L || is.na(value)) NA_real_ else
      as.numeric(value)
  }
  key_logical <- function(key, field) {
    value <- key[[field]]
    is.logical(value) && length(value) == 1L && !is.na(value) && value
  }
  specs <- data.frame(
    parameter          = .brma_random_parameter_io_labels(quantities, "selector"),
    label              = .brma_random_parameter_io_labels(quantities, "label"),
    formula_parameter  = quantities[["formula_parameter"]],
    block              = vapply(keys, key_string, character(1), field = "random_block"),
    grouping           = "",
    structure          = "",
    allocation         = vapply(
      keys,
      key_string,
      character(1),
      field = "allocation_label"
    ),
    random_component   = quantities[["component"]],
    owner_type         = quantities[["owner_type"]],
    owner_name         = quantities[["owner_name"]],
    quantity           = quantities[["quantity"]],
    scale_role         = quantities[["scale_role"]],
    parent_quantity_id = quantities[["parent_quantity_id"]],
    status             = quantities[["status"]],
    allocation_index   = vapply(keys, key_number, numeric(1), field = "index"),
    evaluator          = vapply(keys, key_string, character(1), field = "evaluator"),
    allocation_derived = vapply(
      keys,
      key_logical,
      logical(1),
      field = "allocation_derived"
    ),
    source_type        = quantities[["source_type"]],
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
  specs[["arguments"]]         <- I(quantities[["arguments"]])
  specs[["source_parameter"]]  <- vapply(
    keys, key_string, character(1), field = "source_parameter"
  )
  specs[["source_prior_name"]] <- vapply(
    keys, key_string, character(1), field = "source_prior"
  )
  specs[["source_transform"]]  <- vapply(
    keys, key_string, character(1), field = "source_transform"
  )
  specs[["source_scale"]]      <- vapply(
    keys, key_number, numeric(1), field = "source_scale"
  )
  specs
}

.brma_random_parameter_empty_specs <- function() {

  data.frame(
    parameter         = character(),
    label             = character(),
    formula_parameter = character(),
    block             = character(),
    grouping          = character(),
    structure         = character(),
    allocation        = character(),
    random_component  = character(),
    owner_type        = character(),
    owner_name        = character(),
    quantity          = character(),
    scale_role        = character(),
    parent_quantity_id = character(),
    status            = character(),
    allocation_index  = numeric(),
    evaluator         = character(),
    allocation_derived = logical(),
    arguments         = I(list()),
    source_type       = character(),
    source_parameter  = character(),
    source_prior_name = character(),
    source_transform  = character(),
    source_scale      = numeric(),
    display_transform = I(list()),
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
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
    block_match <- .brma_random_parameter_metadata_matches(
      term[["block_name"]],
      spec[["block"]]
    )
    grouping <- spec[["grouping"]]
    grouping_match <- !is.character(grouping) || length(grouping) != 1L ||
      is.na(grouping) || !nzchar(grouping) ||
      .brma_random_parameter_metadata_matches(term[["group_label"]], grouping)
    block_match && grouping_match
  }, terms)

  if (length(terms) == 1L) terms[[1L]] else NULL
}


.brma_random_parameter_design_allocation <- function(formula_design, spec) {

  allocation_name <- spec[["allocation"]]
  if (!is.list(formula_design) || !is.character(allocation_name) ||
      length(allocation_name) != 1L || is.na(allocation_name) ||
      !nzchar(allocation_name)) {
    return(NULL)
  }
  designs <- Filter(function(design) {
    .brma_random_parameter_metadata_matches(
      design[["parameter"]],
      spec[["formula_parameter"]]
    )
  }, formula_design)
  design_allocations <- unlist(
    lapply(designs, `[[`, "random_allocations"),
    recursive = FALSE
  )
  terms <- unlist(lapply(designs, `[[`, "random_effects"), recursive = FALSE)
  term_allocations <- unlist(lapply(terms, function(random_term) {
    binding <- random_term[["sd_binding"]]
    if (is.null(binding)) list() else binding[["allocations"]]
  }), recursive = FALSE)
  allocations <- c(term_allocations, design_allocations)
  if (length(allocations) == 0L) {
    return(NULL)
  }
  keys <- vapply(allocations, function(allocation) {
    value <- allocation[["label"]]
    if (is.null(value)) "" else value
  }, character(1))
  allocations <- allocations[!duplicated(keys)]
  matches <- vapply(allocations, function(allocation) {
    identical(allocation[["label"]], allocation_name)
  }, logical(1))

  if (sum(matches) == 1L) allocations[[which(matches)]] else NULL
}


.brma_random_parameter_allocation_index <- function(spec, allocation) {

  index <- spec[["allocation_index"]]
  if (!is.numeric(index) || length(index) != 1L || is.na(index) ||
      is.null(allocation)) {
    return(NA_integer_)
  }
  n_targets <- allocation[["n_targets"]]
  if (!is.numeric(n_targets) || length(n_targets) != 1L ||
      is.na(n_targets) || n_targets < 1L) {
    return(NA_integer_)
  }
  if (identical(spec[["quantity"]], "sd_mult") && index > n_targets) {
    index <- index - n_targets
  }
  if (index < 1L || index > n_targets) {
    return(NA_integer_)
  }

  as.integer(index)
}


.brma_random_parameter_normalize_components <- function(components, term) {

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


.brma_random_parameter_fit_with_samples <- function(fit, samples) {

  BayesTools::JAGS_with_draws(fit, samples)
}

.brma_random_parameter_diagnostic_fit <- function(
    object, parameter, standardized_coefficients = FALSE) {

  selected <- .brma_random_parameter_select(
    object                    = object,
    parameter                 = parameter,
    standardized_coefficients = standardized_coefficients,
    chains                    = TRUE
  )
  # The diagnostic density is bounded by the catalog's plotting limits.
  support <- .brma_random_parameter_support(selected, limits = TRUE)
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
    label     = selected[["spec"]][["label"]]
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
    seed                      = seed,
    selections                = list(entry[["selection"]])
  )
  index <- match(entry[["parameter"]], bundle[["specs"]][["parameter"]])
  if (is.na(index)) {
    stop(
      "Random-effect quantity '", entry[["parameter"]],
      "' is no longer available in the fitted draws.",
      call. = FALSE
    )
  }
  spec <- as.list(bundle[["specs"]][index, , drop = FALSE])
  spec[["display_transform"]] <-
    bundle[["specs"]][["display_transform"]][[index]]
  formula_design <- attr(object[["fit"]], "formula_design", exact = TRUE)
  term <- .brma_random_parameter_design_term(formula_design, spec)
  if (!is.null(term)) {
    spec[["grouping"]] <- term[["group_label"]]
    spec[["structure"]] <- term[["structure"]]
  }
  allocation_definition <- .brma_random_parameter_design_allocation(
    formula_design,
    spec
  )
  spec[["allocation_index"]] <- .brma_random_parameter_allocation_index(
    spec,
    allocation_definition
  )

  samples <- bundle[["samples"]][, entry[["parameter"]], drop = FALSE]
  undefined <- BayesTools::posterior_metadata(bundle[["samples"]], "undefined_draws")
  if (entry[["parameter"]] %in% names(undefined)) {
    BayesTools::posterior_metadata(samples, "undefined_draws") <-
      undefined[entry[["parameter"]]]
  }

  list(
    entry        = entry,
    spec         = spec,
    samples      = samples,
    prior        = NULL,
    source_prior = .brma_random_parameter_source_prior(
      object,
      spec
    ),
    allocation_definition = allocation_definition,
    raw_samples            = bundle[["raw_samples"]]
  )
}

.brma_random_parameter_source_prior <- function(object, spec) {

  source_prior_name <- spec[["source_prior_name"]]
  source_prior <- if (is.character(source_prior_name) &&
                        length(source_prior_name) == 1L &&
                        !is.na(source_prior_name) &&
                        nzchar(source_prior_name)) attr(
    object[["fit"]],
    "prior_list",
    exact = TRUE
  )[[source_prior_name]] else NULL
  if (!is.null(source_prior) ||
      !identical(spec[["source_transform"]], "lkj2")) {
    return(source_prior)
  }
  term <- .brma_random_parameter_design_term(
    attr(object[["fit"]], "formula_design", exact = TRUE),
    spec
  )
  eta <- term[["correlation"]][["eta"]]
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta <= 0) {
    return(NULL)
  }

  BayesTools::prior("beta", parameters = list(alpha = eta, beta = eta))
}


.brma_random_parameter_allocation_gate_metadata <- function(selected) {

  allocation <- selected[["allocation_definition"]]
  quantity   <- selected[["spec"]][["quantity"]]
  if (is.null(allocation) ||
      !identical(allocation[["scale"]], "total_variance") ||
      !quantity %in% c("sd_total", "var_total", "var_prop")) {
    return(NULL)
  }
  n_targets <- allocation[["n_targets"]]
  if (!is.numeric(n_targets) || length(n_targets) != 1L ||
      is.na(n_targets) || n_targets != as.integer(n_targets) ||
      n_targets < 1L) {
    stop("Random-effect allocation gate metadata have no valid target count.",
         call. = FALSE)
  }
  n_targets <- as.integer(n_targets)
  component_indicators <- rep(NA_character_, n_targets)
  inclusion <- allocation[["inclusion"]]
  if (is.null(inclusion)) {
    inclusion <- list()
  }
  for (record in inclusion) {
    index     <- record[["index"]]
    indicator <- record[["indicator_name"]]
    if (!is.numeric(index) || length(index) != 1L || is.na(index) ||
        index != as.integer(index) || index < 1L || index > n_targets ||
        !is.character(indicator) || length(indicator) != 1L ||
        is.na(indicator) || !nzchar(indicator)) {
      stop("Random-effect allocation gate metadata are malformed.",
           call. = FALSE)
    }
    component_indicators[[as.integer(index)]] <- indicator
  }
  parent_factors <- allocation[["parent_factors"]]
  if (is.null(parent_factors)) {
    parent_factors <- list()
  }
  parent_indicators <- vapply(parent_factors, function(factor) {
    indicator <- factor[["inclusion_name"]]
    if (is.null(indicator)) {
      return(NA_character_)
    }
    if (!is.character(indicator) || length(indicator) != 1L ||
        is.na(indicator) || !nzchar(indicator)) {
      stop("Random-effect parent-gate metadata are malformed.",
           call. = FALSE)
    }
    indicator
  }, character(1))
  parent_indicators <- unique(parent_indicators[!is.na(parent_indicators)])
  indicators <- c(
    component_indicators[!is.na(component_indicators)],
    parent_indicators
  )
  if (length(indicators) == 0L) {
    return(NULL)
  }
  if (anyDuplicated(indicators)) {
    stop("Random-effect allocation gate metadata contain duplicate indicators.",
         call. = FALSE)
  }

  list(
    quantity             = quantity,
    index                = selected[["spec"]][["allocation_index"]],
    component_indicators = component_indicators,
    parent_indicators    = parent_indicators
  )
}


.brma_random_parameter_allocation_gate_state <- function(metadata,
                                                          raw_samples) {

  if (is.null(metadata)) {
    return(NULL)
  }
  if (!is.matrix(raw_samples)) {
    stop("Random-effect allocation gate samples are unavailable.",
         call. = FALSE)
  }
  n_draws <- nrow(raw_samples)
  component_active <- matrix(
    1,
    nrow = n_draws,
    ncol = length(metadata[["component_indicators"]])
  )
  read_gate <- function(indicator) {
    if (!indicator %in% colnames(raw_samples)) {
      stop(
        "Random-effect allocation samples are missing inclusion indicator '",
        indicator, "'.",
        call. = FALSE
      )
    }
    gate <- raw_samples[, indicator]
    if (any(!is.finite(gate) | !gate %in% c(0, 1))) {
      stop(
        "Random-effect allocation inclusion indicator '", indicator,
        "' is invalid.",
        call. = FALSE
      )
    }
    as.logical(gate)
  }
  for (index in which(!is.na(metadata[["component_indicators"]]))) {
    component_active[, index] <- read_gate(
      metadata[["component_indicators"]][[index]]
    )
  }
  parent_active <- rep(TRUE, n_draws)
  for (indicator in metadata[["parent_indicators"]]) {
    parent_active <- parent_active & read_gate(indicator)
  }
  positive_total <- parent_active & rowSums(component_active) > 0L

  if (metadata[["quantity"]] %in% c("sd_total", "var_total")) {
    return(list(
      defined    = rep(TRUE, n_draws),
      continuous = positive_total,
      point_zero = !positive_total,
      point_one  = rep(FALSE, n_draws)
    ))
  }

  index <- metadata[["index"]]
  if (!is.numeric(index) || length(index) != 1L || is.na(index) ||
      index != as.integer(index) || index < 1L ||
      index > ncol(component_active)) {
    stop("Variance-proportion gate metadata have no valid component index.",
         call. = FALSE)
  }
  index         <- as.integer(index)
  target_active <- component_active[, index]
  other_active <- if (ncol(component_active) == 1L) {
    rep(FALSE, n_draws)
  } else {
    rowSums(component_active[, -index, drop = FALSE]) > 0L
  }

  list(
    defined    = positive_total,
    continuous = positive_total & target_active & other_active,
    point_zero = positive_total & !target_active,
    point_one  = positive_total & target_active & !other_active
  )
}


# 'operation' names what the target is for in the reason returned when there
# is none (e.g. "plots", "point hypotheses").
.brma_random_parameter_density_target <- function(object, parameter,
                                                  operation = "densities") {

  selected  <- .brma_random_parameter_select(object, parameter)
  covariance_update <- BayesTools::random_effects_marginal_update_plan(
    object[["fit"]],
    selected[["entry"]][["selection"]]
  )
  source    <- selected[["spec"]][["source_parameter"]]
  source_type <- selected[["spec"]][["source_type"]]
  type      <- selected[["spec"]][["quantity"]]
  posterior <- as.matrix(.get_posterior_samples(object[["fit"]]))
  conditioning_exclude <- .brma_random_parameter_simplex_exclusions(
    object,
    posterior
  )
  display_transform <- selected[["spec"]][["display_transform"]]
  allocation <- selected[["allocation_definition"]]
  gate_metadata <- .brma_random_parameter_allocation_gate_metadata(selected)
  shared_gate_proportion <- identical(type, "var_prop") &&
    !is.null(gate_metadata) && length(allocation[["inclusion"]]) == 0L &&
    is.character(allocation[["weight_name"]]) &&
    length(allocation[["weight_name"]]) == 1L &&
    !is.na(allocation[["weight_name"]]) && nzchar(allocation[["weight_name"]])
  if (shared_gate_proportion) {
    source            <- allocation[["weight_name"]]
    source_type       <- "identity"
    display_transform <- list(type = "identity")
    selected[["source_prior"]] <- attr(object[["fit"]], "prior_list")[[source]]
  }
  if (source_type %in% c("identity", "one_to_one_transform") &&
      !is.na(source) && nzchar(source) &&
      source %in% colnames(posterior) &&
      !is.null(display_transform)) {
    return(list(
      parameter      = source,
      parameter_spec = list(
        type                 = "primitive",
        target_columns       = source,
        conditioning_exclude = conditioning_exclude,
        covariance_update    = covariance_update
      ),
      display_transform = display_transform
    ))
  }

  if (type %in% c("var_prop", "var_mult", "sd_mult") &&
      identical(selected[["spec"]][["evaluator"]], "allocation") &&
      identical(selected[["spec"]][["source_transform"]], type) &&
      source_type %in% c("identity", "one_to_one_transform") &&
      !is.na(source) && nzchar(source)) {
    metadata  <- selected[["allocation_definition"]]
    index     <- selected[["spec"]][["allocation_index"]]
    n_targets <- metadata[["n_targets"]]
    if (is.numeric(n_targets) && length(n_targets) == 1L &&
        !is.na(n_targets) && n_targets >= 2L && is.numeric(index) &&
        length(index) == 1L && !is.na(index) &&
        index >= 1L && index <= n_targets) {
      n_targets        <- as.integer(n_targets)
      index            <- as.integer(index)
      columns          <- paste0(source, "[", seq_len(n_targets), "]")
      auxiliary_columns <- .iwmde_simplex_auxiliary_columns(
        source,
        n_targets
      )
      allocation_transform <- display_transform
      if (inherits(selected[["source_prior"]], "prior.simplex") &&
          identical(selected[["source_prior"]][["distribution"]], "dirichlet") &&
          !is.null(allocation_transform) &&
          all(c(columns, auxiliary_columns) %in% colnames(posterior))) {
        return(list(
          parameter      = columns[[index]],
          parameter_spec = list(
            type                 = "simplex_pair",
            parameter            = source,
            index                = index,
            n_targets            = n_targets,
            target_columns       = columns,
            auxiliary_columns    = auxiliary_columns,
            conditioning_exclude = columns,
            covariance_update    = covariance_update,
            gate_metadata        = gate_metadata
          ),
          display_transform = allocation_transform
        ))
      }
    }
  }

  if (type %in% c("sd", "sd_total", "sd_common")) {
    target <- .brma_random_parameter_component_density_target(
      object               = object,
      selected             = selected,
      posterior            = posterior,
      conditioning_exclude = conditioning_exclude,
      covariance_update    = covariance_update
    )
    if (!is.null(target)) {
      return(target)
    }
  }

  return(list(
    reason = paste0(
      "qCMDE/IWMDE ", operation, " are not available for random-effect ",
      "quantity '",
      selected[["spec"]][["label"]],
      "' because it has no supported scalar random-component coordinate. ",
      "Use density_method = 'KDE'."
    )
  ))
}


.brma_random_parameter_component_density_target <- function(
    object, selected, posterior, conditioning_exclude, covariance_update) {

  formula_design <- attr(object[["fit"]], "formula_design", exact = TRUE)
  spec <- selected[["spec"]]
  if (!identical(spec[["source_type"]], "composite")) {
    return(NULL)
  }
  if (identical(spec[["evaluator"]], "allocation_sd")) {
    allocation <- selected[["allocation_definition"]]
    if (is.null(allocation) || length(allocation[["inclusion"]]) > 0L) {
      return(NULL)
    }
    factors <- allocation[["parent_factors"]]
  } else if (identical(spec[["evaluator"]], "sd") &&
             isTRUE(spec[["allocation_derived"]])) {
    term <- .brma_random_parameter_design_term(formula_design, spec)
    if (is.null(term) || !.marginalized_random_effect_has_allocation(term)) {
      return(NULL)
    }
    allocation <- term[["sd_binding"]][["allocations"]][[1L]]
    column <- .brma_random_parameter_component_column(term, spec)
    if (is.na(column)) {
      return(NULL)
    }
    factors <- .marginalized_random_effect_allocation_factors(term, column = column)
  } else {
    return(NULL)
  }
  source <- allocation[["source"]]
  if (!is.list(source) || !identical(source[["shape"]], "scalar") ||
      length(factors) == 0L) {
    return(NULL)
  }

  indicators <- unique(unlist(lapply(factors, `[[`, "inclusion_name")))
  gate_metadata <- if (length(indicators) > 0L) list(
    quantity             = "sd_total",
    component_indicators = NA_character_,
    parent_indicators    = indicators
  ) else NULL
  # Inclusion gates determine atoms; continuous rows retain the fitted weights.
  factors <- lapply(factors, function(factor) {

    factor[["inclusion_name"]] <- NULL
    factor
  })
  gate_only <- vapply(factors, function(factor) {

    is.null(factor[["weight_name"]]) && identical(factor[["n_targets"]], 1L)
  }, logical(1))
  factors <- factors[!gate_only]
  source_parameter <- source[["name"]]
  factors          <- lapply(factors, .brma_random_parameter_density_factor)
  if (!is.character(source_parameter) || length(source_parameter) != 1L ||
      is.na(source_parameter) || !nzchar(source_parameter) ||
      !source_parameter %in% colnames(posterior) ||
      any(vapply(factors, is.null, logical(1)))) {
    return(NULL)
  }

  factor_columns <- vapply(factors, function(factor) {
    paste0(factor[["weight_name"]], "[", factor[["index"]], "]")
  }, character(1))
  if (!all(factor_columns %in% colnames(posterior))) {
    return(NULL)
  }

  if (any(!is.finite(posterior[, source_parameter]) |
          posterior[, source_parameter] < 0)) {
    return(NULL)
  }

  if (length(factors) == 0L) {
    return(list(
      parameter      = source_parameter,
      parameter_spec = list(
        type                 = "primitive",
        target_columns       = source_parameter,
        conditioning_exclude = conditioning_exclude,
        covariance_update    = covariance_update,
        gate_metadata        = gate_metadata
      )
    ))
  }

  auxiliary_columns <- unique(unlist(lapply(factors, function(factor) {
    .iwmde_simplex_auxiliary_columns(
      factor[["weight_name"]],
      factor[["n_targets"]]
    )
  }), use.names = FALSE))

  return(list(
    parameter      = source_parameter,
    parameter_spec = list(
      type                 = "random_component_sd",
      source_parameter     = source_parameter,
      factors              = factors,
      target_columns       = source_parameter,
      factor_columns       = factor_columns,
      auxiliary_columns    = auxiliary_columns,
      conditioning_exclude = conditioning_exclude,
      covariance_update    = covariance_update,
      gate_metadata        = gate_metadata
    )
  ))
}


.brma_random_parameter_component_column <- function(term, spec) {

  components <- term[["sd_component_terms"]]
  if (is.null(components)) {
    leaves     <- term[["sd_leaves"]]
    components <- leaves[["leaf_terms"]]
  }
  if (is.null(components)) {
    return(NA_integer_)
  }
  # The SD leaf components are named as BayesTools names formula terms.
  components <- .brma_random_parameter_normalize_components(
    unname(components),
    term
  )
  matches <- which(
    !is.na(components) &
      components == BayesTools::JAGS_parameter_names(spec[["random_component"]])
  )

  if (length(matches) == 1L) as.integer(matches) else NA_integer_
}


.brma_random_parameter_density_factor <- function(factor) {

  fields <- c("weight_name", "index", "scale", "n_targets")
  if (!is.list(factor) || !all(fields %in% names(factor))) {
    return(NULL)
  }
  inclusion_name <- factor[["inclusion_name"]]
  if (is.character(inclusion_name) && length(inclusion_name) == 1L &&
      !is.na(inclusion_name) && nzchar(inclusion_name)) {
    return(NULL)
  }

  out <- factor[fields]
  out[["index"]]     <- as.integer(out[["index"]])
  out[["n_targets"]] <- as.integer(out[["n_targets"]])
  valid <- is.character(out[["weight_name"]]) &&
    length(out[["weight_name"]]) == 1L &&
    !is.na(out[["weight_name"]]) && nzchar(out[["weight_name"]]) &&
    length(out[["index"]]) == 1L && !is.na(out[["index"]]) &&
    length(out[["n_targets"]]) == 1L && !is.na(out[["n_targets"]]) &&
    out[["n_targets"]] >= 2L && out[["index"]] >= 1L &&
    out[["index"]] <= out[["n_targets"]] &&
    out[["scale"]] %in% c("mean_variance", "total_variance")

  if (isTRUE(valid)) out else NULL
}


.brma_random_parameter_simplex_exclusions <- function(object, posterior) {

  prior_list <- attr(object[["fit"]], "prior_list", exact = TRUE)
  if (!is.list(prior_list) || length(prior_list) == 0L) {
    return(character())
  }

  exclusions <- unlist(lapply(names(prior_list), function(parameter) {
    prior <- prior_list[[parameter]]
    if (!inherits(prior, "prior.simplex") ||
        !identical(prior[["distribution"]], "dirichlet")) {
      return(character())
    }
    n_targets <- length(prior[["parameters"]][["alpha"]])
    columns   <- paste0(parameter, "[", seq_len(n_targets), "]")
    if (n_targets < 2L || !all(columns %in% colnames(posterior))) {
      return(character())
    }

    columns[[n_targets]]
  }), use.names = FALSE)

  unique(exclusions)
}

# Exact support of a selected random-effect quantity as declared by the
# parameter catalog (derived by BayesTools from prior provenance), or NULL
# when the catalog cannot derive it.
.brma_random_parameter_catalog_support <- function(selected) {

  quantities <- selected[["entry"]][["selection"]][["quantities"]]
  if (!is.data.frame(quantities) || nrow(quantities) != 1L ||
      !"support" %in% names(quantities)) {
    stop(
      "Random-effect quantity '", selected[["entry"]][["parameter"]],
      "' has no catalog support metadata. Refit the model with the current ",
      "BayesTools version.",
      call. = FALSE
    )
  }

  quantities[["support"]][[1L]]
}

# Bounds of the exact catalog support of a selected random-effect quantity;
# unbounded when the catalog declares no exact support. With 'limits', the
# bounds of a support declared only as plotting limits (not exact) are
# returned as well: they bound plotted densities but exclude no hypothesis.
.brma_random_parameter_support <- function(selected, limits = FALSE) {

  support <- .brma_random_parameter_catalog_support(selected)
  if (is.null(support) || !(isTRUE(support[["exact"]]) || isTRUE(limits))) {
    return(c(-Inf, Inf))
  }

  as.numeric(support[["bounds"]])
}

# Point hypotheses on a random-effect quantity need an exact, regular prior
# ordinate at each value under the BayesTools exactness rule
# (prior_ordinate_status() of the canonical prior density). The status of the
# values (NULL without a prior density); values at prior point masses
# (inclusion and allocation gates, spike components) stop here.
.brma_random_parameter_point_status <- function(selected, prior_density,
                                                values) {

  if (is.null(prior_density)) {
    return(NULL)
  }
  label  <- selected[["spec"]][["label"]]
  status <- BayesTools::prior_ordinate_status(
    prior_density,
    values,
    labels = paste0(label, " = ", format(values, digits = 15L, trim = TRUE))
  )
  point_mass <- which(status[["condition"]] %in% "BayesTools_point_mass_at_null")
  if (length(point_mass) > 0L) {
    value <- format(status[["value"]][[point_mass[[1L]]]], digits = 15L,
                    trim = TRUE)
    stop(
      "Point-null Bayes factors at ", value, " are unavailable for ",
      "random-effect quantity '", label, "' because its prior has a point ",
      "mass at ", value, ". Use a region or directional hypothesis; the ",
      "Component Inclusion table from 'summary(object)' or ",
      "'summary_models(object)' compares the exclusion and inclusion of gated ",
      "components.",
      call. = FALSE
    )
  }

  status
}


# Stops for values without an exact, regular prior ordinate: without a prior
# density, and with the class and message of the refusal of the BayesTools
# exactness rule ('status' from .brma_random_parameter_point_status()).
.brma_random_parameter_check_point_status <- function(selected, status) {

  if (is.null(status)) {
    stop(
      "Point-null Bayes factors are unavailable for random-effect quantity '",
      selected[["spec"]][["label"]], "' because its prior density is ",
      "unavailable. Use a region or directional hypothesis.",
      call. = FALSE
    )
  }
  refused <- which(!status[["eligible"]])
  if (length(refused) > 0L) {
    # The class and message of the refusal of BayesTools' exactness rule.
    stop(structure(
      class = c(
        status[["condition"]][[refused[[1L]]]],
        "BayesTools_hypothesis_ordinate", "error", "condition"
      ),
      list(message = status[["reason"]][[refused[[1L]]]], call = NULL)
    ))
  }

  invisible(status)
}


# Why point hypotheses on a random-effect quantity are unavailable ("" when
# they are available at values that are not prior point masses): the
# quantity needs a canonical prior density whose continuous part BayesTools
# classifies exactly. The classification is taken at a structural interior
# value of the quantity's support, since the support bounds are prior point
# masses or support boundaries for most random-effect quantities.
.brma_random_parameter_point_test_reason <- function(object, selection,
                                                     label) {

  prior_density <- BayesTools::parameter_prior_density(
    object[["fit"]],
    selection
  )
  if (is.null(prior_density)) {
    return(paste0(
      "Point-null Bayes factors are unavailable for random-effect quantity '",
      label, "' because its prior density is unavailable. Use a region or ",
      "directional hypothesis."
    ))
  }
  support <- selection[["quantities"]][["support"]][[1L]]
  value   <- .brma_random_parameter_interior_value(support)
  status  <- BayesTools::prior_ordinate_status(prior_density, value)
  if (identical(status[["condition"]], "BayesTools_inexact_ordinate")) {
    return(paste0(
      "Point-null Bayes factors are unavailable for random-effect quantity '",
      label, "' because its prior ordinate has no exact structural ",
      "classification. Use a region or directional hypothesis."
    ))
  }

  ""
}


# A value inside the support 'support' (a posterior_support_attribute(), or
# NULL for an unbounded support): the midpoint of a bounded support, one unit
# inside a half-bounded support, and 0 otherwise.
.brma_random_parameter_interior_value <- function(support) {

  bounds <- if (is.null(support)) c(-Inf, Inf) else as.numeric(support[["bounds"]])
  if (all(is.finite(bounds))) {
    return(mean(bounds))
  }
  if (is.finite(bounds[[1L]])) {
    return(bounds[[1L]] + 1)
  }
  if (is.finite(bounds[[2L]])) {
    return(bounds[[2L]] - 1)
  }

  0
}


# The draws a random-effect mixed posterior summarizes: when BayesTools left
# out the draws where the quantity is undefined, a mask of the retained share
# of the fitted draws (for the footnote); NULL when every draw is used.
.brma_random_parameter_defined_share <- function(object, values) {

  if (is.null(BayesTools::posterior_metadata(values, "undefined_draws"))) {
    return(NULL)
  }
  n_draws <- sum(vapply(
    coda::as.mcmc.list(object[["fit"]]),
    nrow,
    integer(1)
  ))
  if (length(values) >= n_draws) {
    return(NULL)
  }

  c(rep(TRUE, length(values)), rep(FALSE, n_draws - length(values)))
}


.brma_random_parameter_zero_boundary_alternative <- function(object,
                                                               selected) {

  spec <- selected[["spec"]]
  if (!identical(spec[["quantity"]], "sd") ||
      !isTRUE(spec[["allocation_derived"]])) {
    return(NULL)
  }
  formula_design <- attr(object[["fit"]], "formula_design", exact = TRUE)
  term <- .brma_random_parameter_design_term(formula_design, spec)
  binding <- if (is.null(term)) NULL else term[["sd_binding"]]
  if (is.null(binding) || !isTRUE(binding[["true_allocation"]]) ||
      length(binding[["allocations"]]) != 1L) {
    return(NULL)
  }
  allocation <- binding[["allocations"]][[1L]]
  scale  <- allocation[["scale"]]
  target <- allocation[["target"]]
  if (length(scale) != 1L || is.na(scale) ||
      !scale %in% c("total_variance", "mean_variance") ||
      length(target) != 1L || is.na(target) ||
      !target %in% c("block", "sd_component")) {
    return(NULL)
  }

  component <- if (identical(target, "sd_component")) {
    spec[["random_component"]]
  } else {
    index <- allocation[["index"]]
    component_names <- allocation[["component_names"]]
    if (is.numeric(index) && length(index) == 1L && !is.na(index) &&
        is.character(component_names) && length(component_names) >= index &&
        !is.na(component_names[[index]]) &&
        nzchar(component_names[[index]])) {
      component_names[[index]]
    } else {
      spec[["block"]]
    }
  }
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    return(NULL)
  }
  quantity <- if (identical(scale, "total_variance")) {
    "var_prop"
  } else {
    "var_mult"
  }

  paste0(
    .brma_random_parameter_io_quantity(quantity),
    "(", component, ") = 0"
  )
}


# The mixed posterior of one random-effect quantity with the draw metadata
# of BayesTools::parameter_mixed_posterior(): the catalog support, the
# canonical prior density, the declared atoms (inclusion-gate and
# allocation-gate point masses with masses from the gate states), undefined
# draws and conditioning. 'conditional' conditions on the quantity's
# inclusion event. Returned as the one-element list of as_mixed_posteriors()
# named by the RoBMA parameter name.
.brma_random_parameter_mixed_posterior <- function(
    object, parameter, standardized_coefficients = FALSE,
    conditional = FALSE, selected = NULL) {

  if (is.null(selected)) {
    selected <- .brma_random_parameter_select(
      object                    = object,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients
    )
  }
  fit <- object[["fit"]]
  if (standardized_coefficients) {
    attr(fit, "formula_scale") <- list()
  }
  values <- BayesTools::parameter_mixed_posterior(
    fit         = fit,
    selection   = selected[["entry"]][["selection"]],
    conditional = conditional
  )
  attr(values, "parameter") <- selected[["entry"]][["parameter"]]

  out <- list(values)
  names(out) <- selected[["entry"]][["parameter"]]
  attr(out, "prior_list") <- stats::setNames(
    list(BayesTools::prior_none()),
    selected[["entry"]][["parameter"]]
  )
  attr(out, "random_parameter_label") <- selected[["spec"]][["label"]]
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}
