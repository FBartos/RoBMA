# ============================================================================ #
# brma.mv heterogeneity helpers
# ============================================================================ #
#
# Component-aware heterogeneity extraction for brma.mv() models.
#
# ============================================================================ #


.pooled_heterogeneity_brma_mv <- function(object, probs = c(.025, .975),
                                          conditional = FALSE,
                                          component = "all", ...) {

  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(object)) {
    stop("'conditional' predictions are available only for RoBMA objects.",
         call. = FALSE)
  }

  dots <- list(...)
  .check_unused_dots(
    dots    = dots,
    allowed = ".posterior_samples",
    caller  = "pooled_heterogeneity.brma()"
  )

  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )
  components <- .pooled_brma_mv_heterogeneity_components(
    object            = object,
    posterior_samples = posterior_samples
  )

  component <- .normalize_brma_mv_heterogeneity_component(component)
  if (identical(component, "all")) {
    selected <- components
  } else {
    selected <- .select_brma_mv_heterogeneity_components(
      components = components,
      component  = component,
      object     = object
    )
  }

  chain_info <- .brma_samples_chain_info(
    fit       = object[["fit"]],
    n_samples = nrow(posterior_samples)
  )

  out <- lapply(names(selected), function(name) {
    samples <- selected[[name]]
    colnames(samples) <- "tau"
    .new_brma_samples(
      samples  = samples,
      n_chains = chain_info[["n_chains"]],
      n_iter   = chain_info[["n_iter"]],
      title    = if (identical(name, "total")) {
        "Total Pooled Heterogeneity"
      } else {
        paste0("Pooled Heterogeneity (", name, ")")
      },
      probs    = probs,
      data     = NULL
    )
  })
  names(out) <- names(selected)

  out <- .simplify_brma_mv_heterogeneity_output(out)
  if (conditional) {
    out <- .condition_prediction_samples(
      object            = object,
      samples           = out,
      conditional       = conditional,
      parameters        = .conditional_heterogeneity_parameters(object),
      posterior_samples = posterior_samples,
      quiet             = TRUE
    )
  }

  if (is.list(out) && !is.matrix(out)) {
    out <- .new_brma_samples_list(out)
  }

  return(out)
}


.brma_mv_rms_sd_samples <- function(tau_samples) {

  sqrt(rowMeans(tau_samples^2))
}


.pooled_brma_mv_heterogeneity_components <- function(
    object, posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    posterior_samples
  )

  if (!.is_random(object)) {
    tau_result <- .evaluate.brma.pooled_tau(
      fit               = object[["fit"]],
      scale_data        = object[["data"]][["scale"]],
      scale_formula     = if (.is_scale(object)) {
        .create_fit_formula_list(data = object[["data"]], "scale")
      } else {
        NULL
      },
      scale_priors      = object[["priors"]][["scale"]],
      is_scale          = .is_scale(object),
      is_multilevel     = FALSE,
      posterior_samples = posterior_samples,
      fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
      fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
    )
    return(list(tau = tau_result[["tau_total"]]))
  }

  formula_design <- .brma_mv_pooled_formula_design(object)
  blocks <- unique(vapply(
    formula_design[["random_effects"]],
    `[[`,
    character(1),
    "block_name"
  ))
  location_priors <- attr(object[["fit"]], "prior_list")
  if (is.null(location_priors)) {
    location_priors <- formula_design[["prior_list"]]
  }
  if (is.null(location_priors)) {
    location_priors <- object[["priors"]][["location"]]
  }

  out <- lapply(blocks, function(block) {
    # Pass the compiled design directly so the synthetic expanded row is not
    # reconstructed from a raw moderator row.
    random_vcov <- BayesTools::random_effects_marginal_vcov(
      fit               = formula_design,
      parameter         = formula_design[["parameter"]],
      posterior_samples = posterior_samples,
      prior_list        = location_priors,
      blocks            = block,
      diagonal_only     = TRUE
    )
    variance <- .brma_mv_validate_random_marginal_variance_samples(
      random_vcov = random_vcov,
      S           = nrow(posterior_samples),
      K           = 1L
    )
    matrix(sqrt(variance[, 1L]), ncol = 1L)
  })
  names(out) <- blocks

  return(out)
}


.brma_mv_pooled_formula_design <- function(object) {

  formula_design <- if (.is_scale(object)) {
    .predict_known_v_formula_design_with_row_source_values(
      object = object,
      data   = object[["data"]],
      pooled = TRUE
    )
  } else {
    .fitted_formula_design(object, "mu", required = TRUE)
  }

  terms <- formula_design[["random_effects"]]
  if (length(terms) == 0L) {
    stop("Pooled random heterogeneity requires random-effect metadata.",
         call. = FALSE)
  }

  for (i in seq_along(terms)) {
    term         <- terms[[i]]
    model_matrix <- term[["model_matrix"]]
    if (!is.matrix(model_matrix) || !is.numeric(model_matrix) ||
        nrow(model_matrix) < 1L || any(!is.finite(model_matrix))) {
      stop(
        "Pooled random heterogeneity requires a finite fitted random design.",
        call. = FALSE
      )
    }

    pooled_matrix <- matrix(
      colMeans(model_matrix),
      nrow     = 1L,
      ncol     = ncol(model_matrix),
      dimnames = list("pooled", colnames(model_matrix))
    )
    attr(pooled_matrix, "assign") <- attr(model_matrix, "assign")

    group_map <- 1L
    group_covariance <- term[["group_covariance"]]
    if (is.list(group_covariance) &&
        identical(group_covariance[["type"]], "known")) {
      kernel       <- group_covariance[["kernel"]]
      fitted_group <- term[["group_map"]]
      if (!is.matrix(kernel) || nrow(kernel) != ncol(kernel) ||
          any(!is.finite(kernel)) || any(diag(kernel) < 0) ||
          !is.numeric(fitted_group) ||
          length(fitted_group) != nrow(model_matrix) ||
          anyNA(fitted_group) || any(!is.finite(fitted_group)) ||
          any(fitted_group != as.integer(fitted_group)) ||
          any(fitted_group < 1L | fitted_group > nrow(kernel))) {
        stop("Known-R pooled heterogeneity metadata are invalid.",
             call. = FALSE)
      }
      fitted_multiplier <- diag(kernel)[fitted_group]
      pooled_multiplier <- mean(fitted_multiplier)
      positive_groups   <- which(diag(kernel) > 0)
      if (pooled_multiplier == 0) {
        pooled_matrix[] <- 0
      } else if (length(positive_groups) == 0L) {
        stop("Known-R pooled heterogeneity has no positive marginal variance.",
             call. = FALSE)
      } else {
        group_map     <- positive_groups[[1L]]
        pooled_matrix <- pooled_matrix * sqrt(
          pooled_multiplier / diag(kernel)[group_map]
        )
      }
    }

    term[["model_matrix"]] <- pooled_matrix
    term[["group_map"]]    <- group_map
    terms[[i]]              <- term
  }

  formula_design[["random_effects"]] <- terms
  fixed_matrix <- formula_design[["model_matrix"]]
  if (is.matrix(fixed_matrix) && nrow(fixed_matrix) > 0L) {
    formula_design[["model_matrix"]] <- matrix(
      colMeans(fixed_matrix),
      nrow     = 1L,
      ncol     = ncol(fixed_matrix),
      dimnames = list("pooled", colnames(fixed_matrix))
    )
  }
  source_data <- formula_design[["source_data"]]
  if (is.data.frame(source_data) && nrow(source_data) > 0L) {
    formula_design[["source_data"]] <- source_data[1L, , drop = FALSE]
  }

  return(formula_design)
}


.summary_heterogeneity_brma_mv <- function(object, probs = c(.025, .975),
                                           component = "all", ...) {

  dots              <- list(...)
  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )
  components <- .brma_mv_heterogeneity_components(
    object            = object,
    posterior_samples = posterior_samples
  )
  allocation_summaries <- .summary_heterogeneity_brma_mv_allocations(
    object            = object,
    posterior_samples = posterior_samples,
    probs             = probs
  )

  component <- .normalize_brma_mv_heterogeneity_component(component)
  if (identical(component, "all")) {
    replaced_components <- .brma_mv_allocation_replaced_components(
      allocation_summaries
    )
    selected <- components[!names(components) %in% replaced_components]
    out <- lapply(names(selected), function(name) {
      .summary_heterogeneity_brma_mv_one(
        tau_samples = selected[[name]],
        probs       = probs,
        component   = name
      )
    })
    names(out) <- names(selected)
    out <- c(allocation_summaries, out)
  } else {
    allocation_component <- if (identical(component, "total")) {
      NULL
    } else {
      .match_brma_mv_allocation_summary_component(
        component            = component,
        allocation_summaries = allocation_summaries
      )
    }
    if (!is.null(allocation_component)) {
      out <- allocation_summaries[allocation_component]
    } else {
      selected <- .select_brma_mv_heterogeneity_components(
        components            = components,
        component             = component,
        object                = object,
        extra_component_names = names(allocation_summaries)
      )
      out <- lapply(names(selected), function(name) {
        .summary_heterogeneity_brma_mv_one(
          tau_samples = selected[[name]],
          probs       = probs,
          component   = name
        )
      })
      names(out) <- names(selected)
    }
  }

  out <- .simplify_brma_mv_heterogeneity_output(out)
  if (is.list(out) && length(out) > 1L) {
    class(out) <- c("summary_heterogeneity.brma_list", class(out))
  }

  return(out)
}


.brma_mv_allocation_replaced_components <- function(allocation_summaries) {

  replacements <- attr(allocation_summaries, "component_replacements",
                       exact = TRUE)
  if (length(replacements) == 0L) {
    return(character(0))
  }

  unique(unlist(replacements, use.names = FALSE))
}


.summary_heterogeneity_brma_mv_one <- function(tau_samples, probs, component) {

  tau2_samples <- rowMeans(tau_samples^2)
  samples_list <- list(
    tau  = sqrt(tau2_samples),
    tau2 = tau2_samples
  )

  estimates <- BayesTools::ensemble_estimates_table(
    samples    = samples_list,
    parameters = names(samples_list),
    probs      = probs,
    title      = if (identical(component, "total")) {
      "Total Heterogeneity Estimates:"
    } else {
      paste0("Heterogeneity Estimates (", component, "):")
    }
  )

  output <- list(estimates = estimates)
  class(output) <- "summary_heterogeneity.brma"

  return(output)
}


.summary_heterogeneity_brma_mv_allocations <- function(object,
                                                       posterior_samples,
                                                       probs) {

  allocation_samples <- .brma_mv_allocation_sample_lists(
    object            = object,
    posterior_samples = posterior_samples
  )
  if (length(allocation_samples) == 0L) {
    return(list())
  }

  out <- lapply(names(allocation_samples), function(name) {
    samples_list <- allocation_samples[[name]]
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = samples_list,
      parameters = names(samples_list),
      probs      = probs,
      title      = paste0("Heterogeneity Estimates (", name, "):")
    )

    output <- list(estimates = estimates)
    class(output) <- "summary_heterogeneity.brma"

    output
  })
  names(out) <- names(allocation_samples)

  attr(out, "component_aliases") <- attr(
    allocation_samples,
    "component_aliases",
    exact = TRUE
  )
  attr(out, "component_replacements") <- attr(
    allocation_samples,
    "component_replacements",
    exact = TRUE
  )

  return(out)
}


.brma_mv_allocation_sample_lists <- function(object, posterior_samples) {

  if (!.is_random(object)) {
    return(list())
  }

  formula_design <- .fitted_formula_design(object, "mu", required = TRUE)
  allocations    <- formula_design[["random_allocations"]]
  if (length(allocations) == 0L) {
    return(list())
  }

  K              <- nrow(object[["data"]][["outcome"]])
  source_samples <- .brma_mv_scale_source_samples(
    object            = object,
    posterior_samples = posterior_samples
  )
  out <- list()
  aliases <- list()
  replacements <- list()

  for (allocation in allocations) {
    allocation <- .brma_mv_resolve_sd_component_allocation(
      allocation     = allocation,
      formula_design = formula_design
    )
    target <- allocation[["target"]]
    if (!target %in% c("block", "sd_component")) {
      next
    }
    terms <- allocation[["terms"]]
    if (identical(target, "block") && length(terms) <= 1L) {
      next
    }
    if (identical(target, "sd_component") &&
        length(allocation[["leaf_terms"]]) <= 1L) {
      next
    }

    name <- .brma_mv_allocation_summary_name(
      allocation     = allocation,
      formula_design = formula_design
    )
    samples_list <- .brma_mv_allocation_summary_samples(
      allocation        = allocation,
      formula_design    = formula_design,
      posterior_samples = posterior_samples,
      source_samples    = source_samples,
      K                 = K
    )

    out[[length(out) + 1L]] <- samples_list
    names(out)[[length(out)]] <- name
    aliases[[length(aliases) + 1L]] <- .brma_mv_allocation_summary_aliases(
      name       = name,
      allocation = allocation
    )
    names(aliases)[[length(aliases)]] <- name
    replacements[[length(replacements) + 1L]] <- if (identical(target, "sd_component")) {
      terms <- as.character(allocation[["terms"]])
      terms[!is.na(terms) & nzchar(terms)]
    } else {
      character(0)
    }
    names(replacements)[[length(replacements)]] <- name
  }

  if (length(out) == 0L) {
    return(list())
  }

  original_names <- names(out)
  names(out)     <- make.unique(original_names, sep = "_")
  names(aliases) <- names(out)
  names(replacements) <- names(out)

  attr(out, "component_aliases")      <- aliases
  attr(out, "component_replacements") <- replacements

  return(out)
}


.brma_mv_resolve_sd_component_allocation <- function(allocation,
                                                    formula_design) {

  if (!identical(allocation[["target"]], "sd_component") ||
      length(allocation[["leaf_terms"]]) > 0L) {
    return(allocation)
  }

  terms <- as.character(allocation[["terms"]])
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0L) {
    return(allocation)
  }

  random_effects <- formula_design[["random_effects"]]
  for (term in random_effects) {
    if (!identical(term[["block_name"]], terms[[1L]])) {
      next
    }
    term_allocations <- term[["sd_binding"]][["allocations"]]
    for (term_allocation in term_allocations) {
      if (identical(term_allocation[["target"]], "sd_component") &&
          identical(term_allocation[["weight_name"]],
                    allocation[["weight_name"]])) {
        return(term_allocation)
      }
    }
  }

  allocation
}


.brma_mv_allocation_summary_samples <- function(allocation,
                                                formula_design,
                                                posterior_samples,
                                                source_samples,
                                                K) {

  if (identical(allocation[["target"]], "sd_component")) {
    return(.brma_mv_sd_component_allocation_summary_samples(
      allocation        = allocation,
      term              = .brma_mv_sd_component_allocation_term(
        allocation     = allocation,
        formula_design = formula_design
      ),
      posterior_samples = posterior_samples,
      source_samples    = source_samples,
      K                 = K
    ))
  }

  .brma_mv_block_allocation_summary_samples(
    allocation        = allocation,
    posterior_samples = posterior_samples,
    source_samples    = source_samples,
    K                 = K
  )
}


.brma_mv_block_allocation_summary_samples <- function(allocation,
                                                      posterior_samples,
                                                      source_samples,
                                                      K) {

  total <- .random_sd_source_samples(
    source            = allocation[["source"]],
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples
  )
  total <- .expand_brma_mv_heterogeneity_samples(
    samples = total,
    S       = nrow(posterior_samples),
    K       = K
  )
  total_variance <- rowMeans(total^2)

  samples <- list(
    tau  = sqrt(total_variance),
    tau2 = total_variance
  )

  weight_name <- allocation[["weight_name"]]
  terms       <- allocation[["terms"]]
  labels      <- names(terms)
  if (is.null(labels)) {
    labels <- as.character(terms)
  }
  for (i in seq_along(terms)) {
    column <- paste0(weight_name, "[", i, "]")
    if (!column %in% colnames(posterior_samples)) {
      stop("Missing posterior column: ", column, call. = FALSE)
    }
    weights <- posterior_samples[, column]
    if (any(!is.finite(weights) | weights < 0 | weights > 1)) {
      stop(
        "Invalid random-effect variance allocation samples in posterior ",
        "column: ", column,
        call. = FALSE
      )
    }
    samples[[paste0("var_frac(", allocation[["label"]], ": ",
                    labels[[i]], ")")]] <- weights
  }

  samples
}


.brma_mv_sd_component_allocation_summary_samples <- function(allocation,
                                                            term = NULL,
                                                            posterior_samples,
                                                            source_samples,
                                                            K) {

  source_shape <- allocation[["source"]][["shape"]]
  total        <- .random_sd_source_samples(
    source            = allocation[["source"]],
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples
  )
  total <- .expand_brma_mv_heterogeneity_samples(
    samples = total,
    S       = nrow(posterior_samples),
    K       = K
  )
  total_sd <- sqrt(rowMeans(total^2))

  samples <- list()
  samples[[paste0("sd_total(", allocation[["label"]], ")")]] <- total_sd

  weight_name <- allocation[["weight_name"]]
  labels      <- .brma_mv_sd_component_allocation_labels(allocation)
  block       <- .brma_mv_sd_component_allocation_block(allocation)
  leaf_index  <- allocation[["leaf_index_by_column"]]
  n_targets   <- allocation[["n_targets"]]
  scale       <- allocation[["scale"]]

  for (i in seq_along(labels)) {
    column <- paste0(weight_name, "[", i, "]")
    if (!column %in% colnames(posterior_samples)) {
      stop("Missing posterior column: ", column, call. = FALSE)
    }
    weights <- posterior_samples[, column]
    if (any(!is.finite(weights) | weights < 0 | weights > 1)) {
      stop(
        "Invalid random-effect variance allocation samples in posterior ",
        "column: ", column,
        call. = FALSE
      )
    }

    variance_ratio <- .brma_mv_allocation_variance_ratio(
      weights   = weights,
      scale     = scale,
      n_targets = n_targets
    )
    columns <- which(leaf_index == i)
    if (length(columns) == 0L) {
      sd_samples <- total_sd * sqrt(variance_ratio)
    } else {
      rows <- .brma_mv_sd_component_allocation_rows(
        allocation   = allocation,
        term         = term,
        component    = i,
        source_shape = source_shape,
        K            = ncol(total)
      )
      sd_samples <- sqrt(rowMeans(total[, rows, drop = FALSE]^2) *
                           variance_ratio)
    }

    samples[[paste0("sd(", labels[[i]], " | ", block, ")")]] <- sd_samples
  }

  for (i in seq_along(labels)) {
    column <- paste0(weight_name, "[", i, "]")
    weights <- posterior_samples[, column]
    samples[[paste0("var_ratio(", allocation[["label"]], ": ",
                    labels[[i]], ")")]] <- .brma_mv_allocation_variance_ratio(
      weights   = weights,
      scale     = scale,
      n_targets = n_targets
    )
  }

  samples
}


.brma_mv_sd_component_allocation_term <- function(allocation, formula_design) {

  block <- .brma_mv_sd_component_allocation_block(allocation)
  terms <- formula_design[["random_effects"]]
  if (length(terms) == 0L) {
    return(NULL)
  }
  for (term in terms) {
    if (identical(term[["block_name"]], block)) {
      return(term)
    }
  }

  NULL
}


.brma_mv_sd_component_allocation_rows <- function(allocation, term,
                                                  component, source_shape, K) {

  leaf_index <- allocation[["leaf_index_by_column"]]
  columns    <- which(leaf_index == component)

  if (!identical(source_shape, "row")) {
    return(columns)
  }

  model_matrix <- term[["model_matrix"]]
  if (is.null(term) ||
      !is.matrix(model_matrix) ||
      nrow(model_matrix) != K ||
      length(leaf_index) != ncol(model_matrix) ||
      length(columns) == 0L) {
    stop(
      "Cannot summarize row-varying SD-component allocation because its ",
      "design-column allocation does not map to observation-level SD samples.",
      call. = FALSE
    )
  }

  active <- rowSums(model_matrix[, columns, drop = FALSE] != 0) > 0
  rows   <- which(active)
  if (length(rows) == 0L) {
    stop(
      "Cannot summarize row-varying SD-component allocation because no ",
      "observations map to one of its components.",
      call. = FALSE
    )
  }

  return(rows)
}


.brma_mv_sd_component_allocation_labels <- function(allocation) {

  labels <- allocation[["leaf_terms"]]
  if (length(labels) == 0L) {
    labels <- paste0("component[", seq_len(allocation[["n_targets"]]), "]")
  }
  labels <- unname(as.character(labels))
  labels[is.na(labels) | !nzchar(labels)] <- paste0(
    "component[",
    which(is.na(labels) | !nzchar(labels)),
    "]"
  )

  labels
}


.brma_mv_sd_component_allocation_block <- function(allocation) {

  terms <- as.character(allocation[["terms"]])
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) > 0L) {
    return(terms[[1L]])
  }

  "random"
}


.brma_mv_allocation_variance_ratio <- function(weights, scale, n_targets) {

  if (identical(scale, "mean_variance")) {
    return(n_targets * weights)
  }
  if (identical(scale, "total_variance")) {
    return(weights)
  }

  stop("Random-effect allocation factor scale is unsupported.",
       call. = FALSE)
}


.brma_mv_allocation_summary_name <- function(allocation, formula_design) {

  if (identical(allocation[["target"]], "sd_component")) {
    return(.brma_mv_sd_component_allocation_block(allocation))
  }

  group_labels <- .brma_mv_allocation_group_labels(
    allocation     = allocation,
    formula_design = formula_design
  )
  nested_name <- .brma_mv_nested_allocation_name(group_labels)
  if (!is.null(nested_name)) {
    return(nested_name)
  }

  labels <- allocation[["component_labels"]]
  if (length(labels) == 0L) {
    labels <- names(allocation[["terms"]])
  }
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) > 0L) {
    return(paste(labels, collapse = "+"))
  }

  allocation[["label"]]
}


.brma_mv_allocation_group_labels <- function(allocation, formula_design) {

  random_effects <- formula_design[["random_effects"]]
  if (length(random_effects) == 0L) {
    return(character(0))
  }

  group_labels <- vapply(random_effects, `[[`, character(1), "group_label")
  names(group_labels) <- vapply(random_effects, `[[`, character(1),
                                "block_name")

  terms <- as.character(allocation[["terms"]])
  group_labels[terms[terms %in% names(group_labels)]]
}


.brma_mv_nested_allocation_name <- function(group_labels) {

  group_labels <- group_labels[!is.na(group_labels) & nzchar(group_labels)]
  if (length(group_labels) <= 1L) {
    return(NULL)
  }

  groups <- lapply(group_labels, function(label) {
    strsplit(label, ":", fixed = TRUE)[[1L]]
  })
  group_lengths <- vapply(groups, length, integer(1))
  order         <- order(group_lengths)
  groups        <- groups[order]

  variables <- groups[[1L]]
  for (i in seq.int(2L, length(groups))) {
    if (!all(variables %in% groups[[i]])) {
      return(NULL)
    }
    variables <- c(variables, setdiff(groups[[i]], variables))
  }

  paste(variables, collapse = "/")
}


.brma_mv_allocation_summary_aliases <- function(name, allocation) {

  aliases <- c(
    name,
    allocation[["label"]],
    allocation[["total_name"]],
    allocation[["weight_name"]]
  )
  if (identical(allocation[["target"]], "sd_component")) {
    aliases <- c(
      aliases,
      as.character(allocation[["terms"]]),
      as.character(allocation[["leaf_terms"]])
    )
  }
  aliases <- aliases[!is.na(aliases) & nzchar(aliases)]

  unique(aliases)
}


.match_brma_mv_allocation_summary_component <- function(component,
                                                       allocation_summaries) {

  if (length(allocation_summaries) == 0L) {
    return(NULL)
  }
  if (component %in% names(allocation_summaries)) {
    return(component)
  }

  aliases <- attr(allocation_summaries, "component_aliases", exact = TRUE)
  matches <- names(allocation_summaries)[vapply(names(allocation_summaries),
                                                function(name) {
    component %in% aliases[[name]]
  }, logical(1))]

  if (length(matches) == 1L) {
    return(matches)
  }
  if (length(matches) > 1L) {
    stop(
      "Heterogeneity component '", component, "' is ambiguous. Available ",
      "components: ", paste(names(allocation_summaries), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  NULL
}


.brma_mv_heterogeneity_components <- function(object,
                                              posterior_samples = NULL) {

  posterior_samples <- .get_posterior_samples(object[["fit"]],
                                              posterior_samples)
  K <- nrow(object[["data"]][["outcome"]])

  if (.is_random(object)) {
    return(.brma_mv_random_heterogeneity_components(
      object            = object,
      posterior_samples = posterior_samples,
      K                 = K
    ))
  }

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) {
      .create_fit_formula_list(data = object[["data"]], "scale")
    } else {
      NULL
    },
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = FALSE,
    K                 = K,
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )
  out <- list(tau = .expand_brma_mv_heterogeneity_samples(
    samples = tau_result[["tau_total"]],
    S       = nrow(posterior_samples),
    K       = K
  ))

  return(out)
}


.brma_mv_random_heterogeneity_components <- function(object,
                                                     posterior_samples,
                                                     K,
                                                     data = object[["data"]],
                                                     new_levels = NULL) {

  formula_design <- .fitted_formula_design(object, "mu", required = TRUE)
  terms          <- formula_design[["random_effects"]]
  if (length(terms) == 0L) {
    stop("brma.mv heterogeneity components require random-effect metadata.",
         call. = FALSE)
  }

  pure_intercepts <- vapply(
    terms,
    .brma_mv_random_term_is_pure_intercept,
    logical(1)
  )
  source_samples <- NULL
  if (.is_scale(object) && any(pure_intercepts)) {
    source_samples <- .brma_mv_scale_source_samples(
      object            = object,
      posterior_samples = posterior_samples,
      data              = data
    )
  }

  out <- lapply(seq_along(terms), function(i) {
    term <- terms[[i]]
    if (pure_intercepts[[i]]) {
      samples <- tryCatch(
        .random_effect_term_sd_samples(
          term              = term,
          posterior_samples = posterior_samples,
          K                 = K,
          source_samples    = source_samples
        ),
        error = function(e) {
          .brma_mv_random_block_row_sd_samples(
            object            = object,
            posterior_samples = posterior_samples,
            block             = term[["block_name"]],
            K                 = K,
            data              = data,
            new_levels        = new_levels,
            original_error    = e
          )
        }
      )
    } else {
      samples <- .brma_mv_random_block_row_sd_samples(
        object            = object,
        posterior_samples = posterior_samples,
        block             = term[["block_name"]],
        K                 = K,
        data              = data,
        new_levels        = new_levels,
        original_error    = NULL
      )
    }
    .expand_brma_mv_heterogeneity_samples(
      samples = samples,
      S       = nrow(posterior_samples),
      K       = K
    )
  })
  names(out) <- vapply(terms, `[[`, character(1), "block_name")

  return(out)
}


.brma_mv_random_term_is_pure_intercept <- function(term) {

  if (.random_effect_term_has_known_group_covariance(term) ||
      .marginalized_random_effect_has_row_multiplier(term)) {
    return(FALSE)
  }

  term_formula <- term[["term_formula"]]
  if (!inherits(term_formula, "formula")) {
    return(FALSE)
  }
  formula_terms <- stats::terms(term_formula)

  identical(attr(formula_terms, "intercept"), 1L) &&
    length(attr(formula_terms, "term.labels")) == 0L
}


.random_effect_term_sd_samples <- function(term, posterior_samples, K,
                                           source_samples = NULL) {

  sd_samples <- .marginalized_random_effect_sd_samples(
    term              = term,
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples
  )
  variance <- .marginalized_random_effect_variance_samples(
    term       = term,
    sd_samples = sd_samples,
    K          = K
  )

  return(sqrt(variance))
}


.brma_mv_scale_source_samples <- function(object, posterior_samples,
                                          data = object[["data"]]) {

  if (!.is_scale(object)) {
    return(NULL)
  }

  scale_specs <- .data_scale_component_specs(data)
  if (length(scale_specs) == 0L) {
    return(NULL)
  }

  scale_samples <- .evaluate.brma.scale_terms(
    fit               = object[["fit"]],
    data              = data,
    priors            = object[["priors"]],
    posterior_samples = posterior_samples,
    as_list           = TRUE
  )
  names(scale_samples) <- vapply(scale_specs, `[[`, character(1), "source")

  return(scale_samples)
}


.brma_mv_random_block_row_sd_samples <- function(
    object, posterior_samples, block, K, data = object[["data"]],
    new_levels = NULL, original_error = NULL) {

  random_vcov <- tryCatch(
    .brma_mv_random_effects_marginal_vcov(
      object            = object,
      posterior_samples = posterior_samples,
      blocks            = block,
      data              = data,
      new_levels        = new_levels,
      diagonal_only     = TRUE
    ),
    error = function(e) {
      direct_message <- if (is.null(original_error)) {
        "the random design requires row-marginal covariance evaluation"
      } else {
        paste0("direct SD evaluation failed: ", conditionMessage(original_error))
      }
      stop(
        "Cannot evaluate heterogeneity component '", block, "': ",
        direct_message, "; marginal covariance evaluation failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  variance <- .brma_mv_validate_random_marginal_variance_samples(
    random_vcov = random_vcov,
    S           = nrow(posterior_samples),
    K           = K
  )

  return(sqrt(variance))
}


# Validate the BayesTools diagonal marginal-variance return schema.
.brma_mv_validate_random_marginal_variance_samples <- function(random_vcov,
                                                               S, K) {

  metadata <- random_vcov[["metadata"]]
  variance <- random_vcov[["samples"]]
  valid_metadata <- is.list(metadata) &&
    identical(metadata[["representation"]], "diagonal") &&
    identical(metadata[["quantity"]], "variance") &&
    isTRUE(metadata[["diagonal_only"]]) &&
    identical(metadata[["dense"]], FALSE) &&
    identical(metadata[["n_draws"]], S) &&
    identical(metadata[["n_rows"]], K) &&
    identical(metadata[["row_order"]], seq_len(K)) &&
    is.character(metadata[["row_names"]]) &&
    length(metadata[["row_names"]]) == K

  if (!valid_metadata || !is.matrix(variance) || !is.numeric(variance) ||
      !identical(dim(variance), c(S, K)) || any(!is.finite(variance)) ||
      any(variance < 0) ||
      !identical(colnames(variance), metadata[["row_names"]])) {
    stop(
      "Random-effect marginal variance samples have an invalid diagonal ",
      "variance schema.",
      call. = FALSE
    )
  }

  return(variance)
}


.expand_brma_mv_heterogeneity_samples <- function(samples, S, K) {

  if (!is.matrix(samples) || nrow(samples) != S) {
    stop("Heterogeneity samples must be a posterior-draw matrix.",
         call. = FALSE)
  }
  if (ncol(samples) == K) {
    return(samples)
  }
  if (ncol(samples) == 1L) {
    return(matrix(samples[, 1L], nrow = S, ncol = K))
  }

  stop("Heterogeneity samples must be scalar or row-wise.",
       call. = FALSE)
}


.select_brma_mv_heterogeneity_components <- function(components, component,
                                                     object,
                                                     extra_component_names = character()) {

  component <- .normalize_brma_mv_heterogeneity_component(component)

  if (identical(component, "all")) {
    return(components)
  }

  if (identical(component, "total")) {
    out <- list(total = .total_brma_mv_heterogeneity_samples(components))
    return(out)
  }

  selected <- .match_brma_mv_heterogeneity_component(
    component             = component,
    components            = components,
    object                = object,
    extra_component_names = extra_component_names
  )

  return(components[selected])
}


.normalize_brma_mv_heterogeneity_component <- function(component) {

  if (is.null(component)) {
    return("all")
  }
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    stop("'component' must be a single component name.", call. = FALSE)
  }

  return(component)
}


.match_brma_mv_heterogeneity_component <- function(component, components,
                                                   object,
                                                   extra_component_names = character()) {

  component_names <- names(components)
  if (component %in% component_names) {
    return(component)
  }

  aliases <- .brma_mv_heterogeneity_component_aliases(object, component_names)
  matches <- component_names[vapply(component_names, function(name) {
    component %in% aliases[[name]]
  }, logical(1))]

  if (length(matches) == 1L) {
    return(matches)
  }
  if (length(matches) > 1L) {
    stop(
      "Heterogeneity component '", component, "' is ambiguous. Available ",
      "components: ",
      paste(c(component_names, extra_component_names, "all", "total"),
            collapse = ", "), ".",
      call. = FALSE
    )
  }

  stop(
    "Unknown heterogeneity component '", component, "'. Available components: ",
    paste(c(component_names, extra_component_names, "all", "total"),
          collapse = ", "), ".",
    call. = FALSE
  )
}


.brma_mv_heterogeneity_component_aliases <- function(object, component_names) {

  out <- stats::setNames(as.list(component_names), component_names)

  random_terms <- .data_random_effect_terms(object[["data"]])
  for (term in random_terms) {
    block <- term[["block_name"]]
    if (!block %in% component_names) {
      next
    }
    out[[block]] <- unique(c(
      out[[block]],
      term[["component"]],
      term[["component_label"]],
      term[["component_child_label"]]
    ))
  }

  scale_specs <- .data_scale_component_specs(object[["data"]])
  for (scale_spec in scale_specs) {
    name <- scale_spec[["name"]]
    if (!name %in% component_names) {
      next
    }
    out[[name]] <- unique(c(out[[name]], scale_spec[["aliases"]]))
  }

  lapply(out, function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    unique(x)
  })
}


.total_brma_mv_heterogeneity_samples <- function(components) {

  variance <- Reduce(`+`, lapply(components, function(samples) samples^2))
  sqrt(variance)
}


.simplify_brma_mv_heterogeneity_output <- function(out) {

  if (length(out) == 1L) {
    return(out[[1L]])
  }

  return(out)
}


.check_univariate_heterogeneity_component <- function(component) {

  component <- .normalize_brma_mv_heterogeneity_component(component)
  if (!component %in% c("all", "total", "tau")) {
    stop(
      "Unknown heterogeneity component '", component,
      "'. Available components: tau, all, total.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
