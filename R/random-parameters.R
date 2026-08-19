# Internal semantic random-parameter extraction.

.brma_random_parameter_supported_quantities <- function() {

  c(
    "sd", "sd_total", "var_total", "sd_common", "var_common",
    "cor", "var_prop", "var_ratio", "sd_ratio"
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

  list(
    samples = coda::mcmc.list(semantic_chains),
    specs   = specs,
    priors  = extracted[[1L]][["priors"]]
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
  columns <- lapply(selections, function(selection) {
    as.matrix(BayesTools::parameter_draws(extraction_fit, selection))
  })
  samples <- do.call(cbind, columns)
  colnames(samples) <- quantities[["canonical_name"]]
  specs  <- .brma_random_parameter_specs(quantities)
  specs[["display_transform"]] <- I(lapply(selections, function(selection) {
    BayesTools::parameter_transform(extraction_fit, selection)
  }))
  priors <- stats::setNames(
    rep(list(NULL), nrow(quantities)),
    quantities[["canonical_name"]]
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
    parameter          = quantities[["canonical_name"]],
    label              = sub("^\\([^)]*\\) ", "", quantities[["canonical_name"]]),
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
    value <- allocation[["weight_name"]]
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
  if (identical(spec[["quantity"]], "sd_ratio") && index > n_targets) {
    index <- index - n_targets
  }
  if (index < 1L || index > n_targets) {
    return(NA_integer_)
  }

  as.integer(index)
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
  support <- .brma_random_parameter_support(
    selected[["spec"]],
    selected[["prior"]],
    selected[["source_prior"]],
    selected[["allocation_definition"]]
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

  list(
    entry        = entry,
    spec         = spec,
    samples      = if (chains) {
      bundle[["samples"]][, entry[["parameter"]], drop = FALSE]
    } else {
      bundle[["samples"]][, entry[["parameter"]], drop = FALSE]
    },
    prior        = NULL,
    source_prior = .brma_random_parameter_source_prior(
      object,
      spec
    ),
    allocation_definition = allocation_definition
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

.brma_random_parameter_exact_prior <- function(selected) {

  type         <- selected[["spec"]][["quantity"]]
  transform    <- selected[["spec"]][["source_transform"]]
  source_prior <- selected[["source_prior"]]
  if (type %in% c("sd", "sd_total", "sd_common", "cor", "sd_ratio") &&
      identical(transform, "identity") &&
      !is.null(source_prior) && BayesTools::is.prior(source_prior)) {
    return(source_prior)
  }
  if (identical(transform, "lkj2") &&
      inherits(source_prior, "prior.simple") &&
      identical(source_prior[["distribution"]], "beta") &&
      identical(source_prior[["parameters"]][["alpha"]], 1) &&
      identical(source_prior[["parameters"]][["beta"]], 1)) {
    return(BayesTools::prior(
      "uniform",
      parameters = list(a = -1, b = 1)
    ))
  }
  if (!identical(type, "var_prop")) {
    return(NULL)
  }

  .brma_random_parameter_allocation_source_prior(selected)
}

.brma_random_parameter_allocation_source_prior <- function(selected) {

  source_prior <- selected[["source_prior"]]
  if (is.null(source_prior) ||
      !inherits(source_prior, "prior.simplex") ||
      !identical(source_prior[["distribution"]], "dirichlet")) {
    return(NULL)
  }

  alpha <- source_prior[["parameters"]][["alpha"]]
  index <- selected[["spec"]][["allocation_index"]]
  if (!is.numeric(alpha) || any(!is.finite(alpha)) || any(alpha <= 0) ||
      length(alpha) < 2L || !is.numeric(index) || length(index) != 1L ||
      is.na(index) || index < 1L || index > length(alpha)) {
    return(NULL)
  }

  BayesTools::prior(
    "beta",
    parameters = list(
      alpha = alpha[[index]],
      beta  = sum(alpha[-index])
    )
  )
}

.brma_random_parameter_density_target <- function(object, parameter) {

  selected  <- .brma_random_parameter_select(object, parameter)
  source    <- selected[["spec"]][["source_parameter"]]
  source_type <- selected[["spec"]][["source_type"]]
  type      <- selected[["spec"]][["quantity"]]
  posterior <- as.matrix(object[["fit"]][["mcmc"]])
  conditioning_exclude <- .brma_random_parameter_simplex_exclusions(
    object,
    posterior
  )
  display_transform <- selected[["spec"]][["display_transform"]]
  if (source_type %in% c("identity", "one_to_one_transform") &&
      !is.na(source) && nzchar(source) &&
      source %in% colnames(posterior) &&
      !is.null(display_transform)) {
    return(list(
      parameter      = source,
      parameter_spec = list(
        type                 = "primitive",
        target_columns       = source,
        conditioning_exclude = conditioning_exclude
      ),
      display_transform = display_transform
    ))
  }

  if (type %in% c("var_prop", "var_ratio", "sd_ratio") &&
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
            prior_density        =
              .brma_random_parameter_allocation_source_prior(selected)
          ),
          display_transform = allocation_transform
        ))
      }
    }
  }

  if (identical(type, "sd")) {
    target <- .brma_random_parameter_component_density_target(
      object               = object,
      selected             = selected,
      posterior            = posterior,
      conditioning_exclude = conditioning_exclude
    )
    if (!is.null(target)) {
      return(target)
    }
  }

  return(list(
    reason = paste0(
      "qCMDE/IWMDE plots are not available for random-effect quantity '",
      selected[["spec"]][["label"]],
      "' because it has no supported scalar random-component coordinate. ",
      "Use density_method = 'KDE'."
    )
  ))
}


.brma_random_parameter_component_density_target <- function(
    object, selected, posterior, conditioning_exclude) {

  formula_design <- attr(object[["fit"]], "formula_design", exact = TRUE)
  spec <- selected[["spec"]]
  if (!identical(spec[["source_type"]], "composite") ||
      !identical(spec[["evaluator"]], "sd") ||
      !isTRUE(spec[["allocation_derived"]])) {
    return(NULL)
  }
  term <- .brma_random_parameter_design_term(
    formula_design,
    spec
  )
  if (is.null(term) || !.marginalized_random_effect_has_allocation(term)) {
    return(NULL)
  }

  allocation <- term[["sd_binding"]][["allocations"]][[1L]]
  source     <- allocation[["source"]]
  column     <- .brma_random_parameter_component_column(
    term,
    spec
  )
  if (is.na(column)) {
    return(NULL)
  }
  factors <- .marginalized_random_effect_allocation_factors(
    term,
    column = column
  )
  if (!is.list(source) || !identical(source[["shape"]], "scalar") ||
      length(factors) == 0L) {
    return(NULL)
  }

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
      conditioning_exclude = conditioning_exclude
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
  components <- .brma_random_parameter_normalize_components(
    unname(components),
    term
  )
  matches <- which(
    !is.na(components) & components == spec[["random_component"]]
  )

  if (length(matches) == 1L) as.integer(matches) else NA_integer_
}


.brma_random_parameter_density_factor <- function(factor) {

  fields <- c("weight_name", "index", "scale", "n_targets")
  if (!is.list(factor) || !all(fields %in% names(factor))) {
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

.brma_random_parameter_support <- function(spec, prior = NULL,
                                           source_prior = NULL,
                                           allocation = NULL) {

  type <- spec[["quantity"]]
  support <- if (type %in% c(
    "sd", "sd_total", "var_total", "sd_common", "var_common"
  )) {
    c(0, Inf)
  } else if (identical(type, "sd_ratio")) {
    scale <- if (is.null(allocation)) NULL else allocation[["scale"]]
    if (is.null(allocation)) {
      c(0, Inf)
    } else if (identical(scale, "mean_variance")) {
      n_targets <- allocation[["n_targets"]]
      if (!is.numeric(n_targets) || length(n_targets) != 1L ||
          is.na(n_targets) || n_targets < 1) {
        stop(
          "SD-ratio metadata are missing a valid allocation target count.",
          call. = FALSE
        )
      }
      c(0, sqrt(n_targets))
    } else if (identical(scale, "total_variance")) {
      c(0, 1)
    } else {
      stop(
        "SD-ratio metadata are missing a canonical allocation scale.",
        call. = FALSE
      )
    }
  } else if (identical(type, "cor")) {
    c(-1, 1)
  } else if (identical(type, "var_prop")) {
    c(0, 1)
  } else if (identical(type, "var_ratio")) {
    upper <- if (is.null(allocation[["n_targets"]])) Inf else
      as.numeric(allocation[["n_targets"]])
    c(0, upper)
  } else {
    c(-Inf, Inf)
  }

  if (!is.null(source_prior) && !is.null(source_prior[["truncation"]])) {
    truncation <- source_prior[["truncation"]]
    lower      <- truncation[["lower"]]
    upper      <- truncation[["upper"]]
    if (length(lower) == 1L && length(upper) == 1L) {
      transform <- spec[["display_transform"]]
      if (is.null(transform) ||
          (identical(transform[["type"]], "square") && lower < 0)) {
        return(support)
      }
      transformed <- BayesTools::parameter_transform_forward(
        c(lower, upper),
        transform
      )
      if (anyNA(transformed)) {
        return(support)
      }
      lower <- min(transformed)
      upper <- max(transformed)
      support <- c(
        max(support[1L], as.numeric(lower)),
        min(support[2L], as.numeric(upper))
      )
    }
  }

  support
}

.brma_random_parameter_point_test_reason <- function(
    spec, prior = NULL, source_prior = NULL, derived = FALSE) {

  type   <- spec[["quantity"]]
  source <- spec[["source_parameter"]]
  label  <- spec[["label"]]
  if (.brma_random_parameter_prior_has_atom(prior) ||
      .brma_random_parameter_prior_has_atom(source_prior)) {
    return(paste0(
      "Point-null Bayes factors are not available for random-effect quantity '",
      label, "' because its induced prior/posterior contains a point mass. ",
      "Use a region or directional hypothesis."
    ))
  }
  if (identical(type, "cor") &&
      (is.na(source) || !nzchar(source)) && !derived) {
    return(paste0(
      "Point-null Bayes factors are not available for derived pairwise ",
      "correlation '", label, "'."
    ))
  }
  if (identical(type, "sd") &&
      (is.na(source) || !nzchar(source)) && !derived) {
    return(paste0(
      "Point-null Bayes factors are not available for derived component SD '",
      label, "'."
    ))
  }

  ""
}

.brma_random_parameter_prior_has_atom <- function(prior) {

  if (is.null(prior) || !BayesTools::is.prior(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.point(prior)) {
    return(TRUE)
  }
  if (BayesTools::is.prior.mixture(prior) ||
      BayesTools::is.prior.spike_and_slab(prior)) {
    return(any(vapply(
      prior,
      .brma_random_parameter_prior_has_atom,
      logical(1)
    )))
  }

  return(FALSE)
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
    prior = FALSE, n_prior_samples = 10000L, seed = NULL,
    selected = NULL, prior_selected = NULL) {

  if (is.null(selected)) {
    selected <- .brma_random_parameter_select(
      object                    = object,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients
    )
  }
  values <- unname(as.numeric(selected[["samples"]][, 1L]))
  attr(values, "sample_ind") <- FALSE
  attr(values, "models_ind") <- rep(1, length(values))
  attr(values, "parameter")  <- selected[["entry"]][["parameter"]]
  attr(values, "prior_list") <- BayesTools::prior_none()
  support <- .brma_random_parameter_support(
    selected[["spec"]],
    selected[["prior"]],
    selected[["source_prior"]],
    selected[["allocation_definition"]]
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
  if (!.brma_random_parameter_prior_has_atom(selected[["prior"]]) &&
      !.brma_random_parameter_prior_has_atom(selected[["source_prior"]])) {
    attr(values, "posterior_atoms") <- BayesTools::posterior_atom_attribute(
      source = "RoBMA semantic random-effect prior"
    )
  }

  target_prior <- BayesTools::prior_none()
  if (prior) {
    target_prior <- .brma_random_parameter_exact_prior(selected)
  }
  if (prior && is.null(target_prior)) {
    if (is.null(prior_selected)) {
      prior_selected <- .brma_random_parameter_select(
        object                    = object,
        parameter                 = parameter,
        standardized_coefficients = standardized_coefficients,
        prior                     = TRUE,
        n_prior_samples           = n_prior_samples,
        seed                      = seed
      )
    }
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
    target_prior <- BayesTools::prior_none()
  }
  attr(values, "prior_list") <- target_prior

  class(values) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "marginal_posterior.simple",
    "marginal_posterior"
  )
  out <- list(values)
  names(out) <- selected[["entry"]][["parameter"]]
  attr(out, "prior_list") <- stats::setNames(
    list(target_prior),
    selected[["entry"]][["parameter"]]
  )
  attr(out, "random_parameter_label") <- selected[["spec"]][["label"]]
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}
