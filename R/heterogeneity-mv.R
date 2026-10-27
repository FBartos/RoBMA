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
  .warn_unused_dots(
    dots    = dots,
    allowed = ".posterior_samples",
    caller  = "pooled_heterogeneity.brma()"
  )

  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )
  component <- .normalize_brma_mv_heterogeneity_component(component)
  components <- .pooled_brma_mv_heterogeneity_components(
    object            = object,
    posterior_samples = posterior_samples,
    covariance        = identical(component, "total")
  )
  component_quantities <- attr(components, "quantities", exact = TRUE)

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
    quantity <- if (identical(name, "total")) {
      if (length(components) > 1L) {
        "sd_total"
      } else {
        unname(component_quantities[[1L]])
      }
    } else {
      unname(component_quantities[[name]])
    }
    colnames(samples) <- .brma_random_parameter_io_quantity(quantity)
    .new_brma_samples(
      samples   = samples,
      n_chains  = chain_info[["n_chains"]],
      n_iter    = chain_info[["n_iter"]],
      title     = if (identical(name, "total")) {
        "Total Pooled Heterogeneity"
      } else {
        paste0("Pooled Heterogeneity (", name, ")")
      },
      component = name,
      probs     = probs,
      data      = NULL
    )
  })
  names(out) <- names(selected)

  if (conditional) {
    inclusion_parameters <- .random_component_conditioning_parameters(
      object     = object,
      components = names(out),
      fallback   = .conditional_heterogeneity_parameters(object)
    )
    out <- lapply(names(out), function(name) {
      .condition_prediction_samples(
        object            = object,
        samples           = out[[name]],
        conditional       = TRUE,
        parameters        = inclusion_parameters[[name]],
        posterior_samples = posterior_samples,
        quiet             = TRUE
      )
    })
    names(out) <- names(selected)
  }

  out <- .simplify_brma_mv_heterogeneity_output(out)

  if (is.list(out) && !is.matrix(out)) {
    out <- .new_brma_samples_list(out)
  }

  return(out)
}


.brma_mv_rms_sd_samples <- function(sd_samples) {

  sqrt(rowMeans(sd_samples^2))
}


.pooled_brma_mv_heterogeneity_components <- function(
    object, posterior_samples = NULL, covariance = FALSE) {

  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    posterior_samples
  )

  if (!.is_random(object)) {
    tau_result <- .evaluate.brma.pooled_tau(
      fit               = .object_formula_fit(object),
      scale_data        = object[["data"]][["scale"]],
      scale_formula     = if (.is_scale(object)) {
        .create_fit_formula_list(data = object[["data"]], "scale")
      } else {
        NULL
      },
      scale_priors      = object[["priors"]][["scale"]],
      is_scale          = .is_scale(object),
      is_multilevel     = FALSE,
      is_random         = .is_random(object),
      posterior_samples = posterior_samples,
      fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
      fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
    )
    out <- list(estimate = tau_result[["tau_total"]])
    attr(out, "quantities") <- c(estimate = "sd")
    return(out)
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
  allocation_components <- if (covariance) {
    list()
  } else {
    .brma_mv_pooled_allocation_block_samples(
      object            = object,
      formula_design    = formula_design,
      posterior_samples = posterior_samples
    )
  }

  covariance_blocks <- setdiff(blocks, names(allocation_components))
  covariance_components <- list()
  if (length(covariance_blocks) > 0L) {
    factors <- BayesTools::random_effects_marginal_factor_states(
      fit               = formula_design,
      parameter         = formula_design[["parameter"]],
      posterior_samples = posterior_samples,
      prior_list        = location_priors,
      blocks            = covariance_blocks,
      row_blocks        = list(1L)
    )
    variance <- BayesTools::random_effects_marginal_factor_diagonal(
      factors,
      by_block = TRUE
    )
    covariance_components <- lapply(variance, function(samples) {
      if (any(!is.finite(samples)) || any(samples < 0)) {
        stop("Pooled heterogeneity requires finite nonnegative random variance.",
             call. = FALSE)
      }
      matrix(sqrt(samples[, 1L]), ncol = 1L)
    })
  }

  out <- lapply(blocks, function(block) {
    if (block %in% names(allocation_components)) {
      return(allocation_components[[block]][["samples"]])
    }
    covariance_components[[block]]
  })
  names(out) <- blocks
  quantities <- stats::setNames(rep("sd", length(blocks)), blocks)
  for (block in intersect(blocks, names(allocation_components))) {
    quantities[[block]] <- allocation_components[[block]][["quantity"]]
  }
  if (covariance) {
    for (term in formula_design[["random_effects"]]) {
      if (term[["n_columns"]] > 1L) {
        quantities[[term[["block_name"]]]] <- "sd_total"
      }
    }
  }
  attr(out, "quantities") <- quantities

  return(out)
}


.brma_mv_pooled_allocation_block_samples <- function(
    object, formula_design, posterior_samples) {

  allocations <- formula_design[["random_allocations"]]
  if (length(allocations) == 0L) {
    return(list())
  }
  K   <- nrow(object[["data"]][["outcome"]])
  out <- list()

  for (allocation in allocations) {
    allocation <- .brma_mv_resolve_sd_component_allocation(
      allocation     = allocation,
      formula_design = formula_design
    )
    if (!identical(allocation[["target"]], "sd_component") ||
        length(allocation[["leaf_terms"]]) <= 1L) {
      next
    }
    block <- unique(as.character(allocation[["terms"]]))
    block <- block[!is.na(block) & nzchar(block)]
    if (length(block) != 1L) {
      stop(
        "Pooled SD-component allocation metadata must identify one random ",
        "block.",
        call. = FALSE
      )
    }
    quantity <- .brma_mv_allocation_aggregate_quantities(allocation)[["sd"]]

    source <- allocation[["source"]]
    if (identical(source[["shape"]], "row")) {
      evaluator <- source[["source"]][["values"]]
      inputs <- source[["source"]][["inputs"]]
      if (!is.function(evaluator)) {
        stop("Pooled SD-component allocation requires its prepared row-source evaluator.",
             call. = FALSE)
      }
      missing <- setdiff(inputs, colnames(posterior_samples))
      if (length(missing) > 0L) {
        stop("Pooled SD-component allocation is unavailable; missing source coordinates: ",
             paste(missing, collapse = ", "), ".", call. = FALSE)
      }
      samples <- vapply(seq_len(nrow(posterior_samples)), function(i) {
        evaluator(parameters = as.list(stats::setNames(
          as.numeric(posterior_samples[i, inputs, drop = FALSE]), inputs)),
          data = formula_design[["source_data"]], n_rows = 1L)
      }, numeric(1L))
      samples <- matrix(samples, ncol = 1L)
      for (factor in allocation[["parent_factors"]]) {
        samples <- samples * .random_allocation_factor_samples(
          factor = factor, posterior_samples = posterior_samples
        )
      }
      if (any(!is.finite(samples)) || any(samples < 0)) {
        stop("Pooled heterogeneity requires finite nonnegative allocation SDs.",
             call. = FALSE)
      }
    } else if (!is.null(object[["fit"]])) {
      catalog    <- BayesTools::parameter_catalog(object[["fit"]])
      quantities <- catalog[["quantities"]]
      keys       <- quantities[["extraction_key"]]
      rows <- which(
        !quantities[["internal"]] &
          quantities[["quantity"]] == quantity &
          vapply(keys, function(key) {
            identical(key[["allocation_label"]], allocation[["label"]])
          }, logical(1))
      )
      if (length(rows) != 1L) {
        stop(
          "Pooled SD-component allocation metadata do not identify exactly ",
          "one aggregate SD quantity.",
          call. = FALSE
        )
      }
      selection <- BayesTools::parameter_catalog_resolve(
        catalog   = catalog,
        alias     = quantities[["canonical_name"]][rows],
        namespace = quantities[["namespace"]][rows]
      )
      draws <- BayesTools::parameter_draws(
        object[["fit"]],
        selection,
        model_samples = posterior_samples
      )
      samples <- matrix(as.numeric(draws[[1L]][, 1L]), ncol = 1L)
    } else {
      source_samples <- .brma_mv_scale_source_samples(
        object = object, posterior_samples = posterior_samples
      )
      name   <- source[["name"]]
      available <- !is.null(source_samples) && name %in% names(source_samples)
      if (!available && identical(source[["shape"]], "scalar")) {
        available <- name %in% colnames(posterior_samples)
      }
      if (!available && identical(source[["shape"]], "row")) {
        available <- all(
          paste0(name, "[", seq_len(K), "]") %in% colnames(posterior_samples)
        )
      }
      if (!available) {
        next
      }
      samples <- .brma_mv_allocation_observed_sources(
        allocation        = allocation,
        posterior_samples = posterior_samples,
        K                 = K,
        source_samples    = source_samples
      )[[name]]
      samples <- matrix(
        sqrt(rowMeans(samples^2)),
        ncol = 1L
      )
    }

    out[[block]] <- list(samples = samples, quantity = quantity)
  }

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

    term[["model_matrix"]] <- pooled_matrix
    term[["group_map"]]    <- 1L
    terms[[i]]              <- term
  }

  formula_design[["random_effects"]] <- terms
  formula_design <- .brma_mv_remove_known_group_covariance(formula_design)
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

  component <- .normalize_brma_mv_heterogeneity_component(component)
  dots              <- list(...)
  posterior_samples <- .get_posterior_samples(
    object[["fit"]],
    dots[[".posterior_samples"]]
  )
  components <- .brma_mv_heterogeneity_components(
    object                         = object,
    posterior_samples              = posterior_samples,
    include_known_group_covariance = FALSE
  )
  unavailable <- .brma_mv_original_scale_unavailable(object)
  allocation_summaries <- .summary_heterogeneity_brma_mv_allocations(
    object            = object,
    posterior_samples = posterior_samples,
    probs             = probs,
    component         = component,
    unavailable       = unavailable
  )
  correlation_samples <- .brma_mv_correlation_sample_lists(
    object            = object,
    posterior_samples = posterior_samples
  )
  if (length(components) == 1L) {
    correlation_samples[["total"]] <- correlation_samples[[names(components)]]
  }

  if (identical(component, "all")) {
    replaced_components <- .brma_mv_allocation_replaced_components(
      allocation_summaries
    )
    selected <- components[!names(components) %in% replaced_components]
    out <- lapply(names(selected), function(name) {
      .summary_heterogeneity_brma_mv_one(
        sd_samples          = selected[[name]],
        probs               = probs,
        component           = name,
        aggregate           = FALSE,
        correlation_samples = correlation_samples[[name]],
        unavailable         = unavailable
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
        extra_component_names = attr(allocation_summaries, "available_components",
                                     exact = TRUE)
      )
      out <- lapply(names(selected), function(name) {
        .summary_heterogeneity_brma_mv_one(
          sd_samples          = selected[[name]],
          probs               = probs,
          component           = name,
          aggregate           = identical(name, "total") &&
            length(components) > 1L,
          correlation_samples = correlation_samples[[name]],
          unavailable         = unavailable
        )
      })
      names(out) <- names(selected)
    }
  }

  out <- .simplify_brma_mv_heterogeneity_output(out)
  if (is.list(out) && length(out) > 1L &&
      !inherits(out, "summary_heterogeneity.brma")) {
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


.summary_heterogeneity_brma_mv_one <- function(sd_samples, probs, component,
                                               aggregate = FALSE,
                                               correlation_samples = NULL,
                                               unavailable = character()) {

  var_samples <- rowMeans(sd_samples^2)
  sd_name     <- if (aggregate) "tau_total" else "tau"
  var_name    <- if (aggregate) "tau2_total" else "tau2"
  samples_list <- list(
    sqrt(var_samples),
    var_samples
  )
  names(samples_list) <- c(sd_name, var_name)
  samples_list <- c(samples_list, correlation_samples)

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
  estimates <- .summary_heterogeneity_mark_unavailable(
    estimates, samples_list, unavailable
  )

  output <- list(estimates = estimates, component = component)
  class(output) <- "summary_heterogeneity.brma"

  return(output)
}


# The random-effect quantities that summary() reports as NA on the original
# scale because the fitted random block cannot represent them, by catalog name,
# with the note of the row under RoBMA's argument. BayesTools owns the rule: it
# is read from the table summary() is built from, so both tables agree. The set
# depends on the design alone, so the table is built on a few draws of the fit,
# not on all of them.
.brma_mv_original_scale_unavailable <- function(object, n_draws = 3L) {

  fit <- object[["fit"]]
  if (!.is_random(object) || is.null(fit) ||
      length(attr(fit, "formula_scale", exact = TRUE)) == 0L) {
    return(character())
  }

  draws <- as.matrix(.get_posterior_samples(fit))
  keep  <- unique(round(seq(1, nrow(draws), length.out = min(nrow(draws), n_draws))))
  object[["fit"]] <- BayesTools::JAGS_with_draws(
    fit, coda::mcmc.list(coda::mcmc(draws[keep, , drop = FALSE]))
  )
  table <- .summary_estimates_table(
    object                   = object,
    probs                    = c(.025, .975),
    include_mcmc_diagnostics = FALSE,
    is_robma                 = .is_robust_RoBMA(object),
    keep_parameters          = "random",
    transform_scaled         = TRUE,
    random_effects_summary   = "full",
    formula_prefix           = FALSE
  )
  table      <- .summary_translate_scale_advice(table)
  footnotes  <- attr(table, "footnotes", exact = TRUE)
  parameters <- attr(table, "parameters", exact = TRUE)
  if (length(footnotes) == 0L) {
    return(character())
  }

  # A row is unavailable when it is NA with the BayesTools note of an NA row
  # ("<row>: NA: <reason>").
  notes <- stats::setNames(rep(NA_character_, nrow(table)), parameters)
  for (row in seq_len(nrow(table))) {
    prefix    <- paste0(rownames(table)[row], ": NA: ")
    note      <- footnotes[names(footnotes) == rownames(table)[row]]
    note      <- note[startsWith(note, prefix)]
    if (length(note) > 0L && is.na(table[row, "Mean"])) {
      notes[[row]] <- substring(note[[1L]], nchar(prefix) - 3L)
    }
  }

  return(notes[!is.na(notes)])
}


# Rows of a heterogeneity table whose quantity is unavailable on the original
# scale are NA, with the note summary() gives the same row. The draws of those
# quantities carry their catalog name.
.summary_heterogeneity_mark_unavailable <- function(estimates, samples,
                                                    unavailable) {

  if (length(unavailable) == 0L) {
    return(estimates)
  }
  names_catalog <- vapply(samples, function(values) {
    name <- attr(values, "catalog_name", exact = TRUE)
    if (is.null(name)) NA_character_ else name
  }, character(1))
  marked <- which(names_catalog %in% names(unavailable))
  if (length(marked) == 0L) {
    return(estimates)
  }

  labels <- names(samples)[marked]
  rows   <- match(labels, rownames(estimates))
  rows   <- rows[!is.na(rows)]
  estimates[rows, ] <- NA_real_
  notes <- stats::setNames(
    paste0(labels, ": ", unavailable[names_catalog[marked]]),
    labels
  )
  attr(estimates, "footnotes") <- c(
    attr(estimates, "footnotes", exact = TRUE), notes
  )

  return(estimates)
}


.brma_mv_correlation_sample_lists <- function(object, posterior_samples) {

  if (!.is_random(object) || is.null(object[["fit"]])) {
    return(list())
  }

  catalog    <- BayesTools::parameter_catalog(object[["fit"]])
  quantities <- catalog[["quantities"]]
  selected <- quantities[
    quantities[["namespace"]] == "mu" &
      quantities[["owner_type"]] == "random_block" &
      quantities[["quantity"]] == "cor" &
      !quantities[["internal"]] &
      quantities[["status"]] != "unavailable",
    ,
    drop = FALSE
  ]
  out <- list()
  for (i in seq_len(nrow(selected))) {
    selection <- BayesTools::parameter_catalog_resolve(
      catalog   = catalog,
      alias     = selected[["canonical_name"]][i],
      namespace = selected[["namespace"]][i]
    )
    draws <- BayesTools::parameter_draws(
      object[["fit"]],
      selection,
      model_samples = posterior_samples
    )
    block <- selected[["extraction_key"]][[i]][["random_block"]]
    label <- .brma_random_parameter_io_labels(
      selected[i, , drop = FALSE],
      "label"
    )
    values <- .brma_mv_catalog_draw_values(draws)
    attr(values, "catalog_name") <- selected[["canonical_name"]][i]
    out[[block]][[label]] <- values
  }

  return(out)
}


# The draws of one catalog quantity as a plain vector. Quantities that can be
# undefined in some draws (original-scale correlations with a zero SD) keep
# the 'undefined_draws' declaration of parameter_draws(), so summaries use the
# defined draws and footnote their share instead of rejecting the NA draws.
.brma_mv_catalog_draw_values <- function(draws) {

  values    <- as.numeric(draws[[1L]][, 1L])
  undefined <- BayesTools::posterior_metadata(draws, "undefined_draws")
  if (!is.null(undefined)) {
    BayesTools::posterior_metadata(values, "undefined_draws") <- unname(undefined[[1L]])
  }

  return(values)
}


.summary_heterogeneity_brma_mv_allocations <- function(object,
                                                       posterior_samples,
                                                       probs,
                                                       component = "all",
                                                       unavailable = character()) {

  if (identical(component, "total")) {
    return(list())
  }

  allocation_samples <- .brma_mv_allocation_sample_lists(
    object            = object,
    posterior_samples = posterior_samples,
    component         = component
  )
  available_components <- attr(allocation_samples, "available_components",
                               exact = TRUE)
  if (is.null(available_components)) {
    available_components <- names(allocation_samples)
  }
  selected <- if (identical(component, "all")) {
    names(allocation_samples)
  } else {
    .match_brma_mv_allocation_summary_component(component, allocation_samples)
  }

  out <- lapply(selected, function(name) {
    samples_list <- allocation_samples[[name]]
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = samples_list,
      parameters = names(samples_list),
      probs      = probs,
      title      = paste0("Heterogeneity Estimates (", name, "):")
    )
    estimates <- .summary_heterogeneity_mark_unavailable(
      estimates, samples_list, unavailable
    )

    output <- list(estimates = estimates, component = name)
    class(output) <- "summary_heterogeneity.brma"

    output
  })
  names(out) <- selected
  attr(out, "available_components") <- available_components

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


.brma_mv_allocation_sample_lists <- function(object, posterior_samples,
                                              component = "all") {

  if (!.is_random(object)) {
    return(list())
  }

  if (!is.null(object[["fit"]])) {
    return(.brma_mv_allocation_catalog_sample_lists(
      object            = object,
      posterior_samples = posterior_samples,
      component         = component
    ))
  }

  formula_design <- .fitted_formula_design(object, "mu", required = TRUE)
  allocations    <- formula_design[["random_allocations"]]
  if (length(allocations) == 0L) {
    return(list())
  }

  K       <- nrow(object[["data"]][["outcome"]])
  entries <- list()

  for (allocation in allocations) {
    allocation <- .brma_mv_resolve_sd_component_allocation(
      allocation     = allocation,
      formula_design = formula_design
    )
    target <- allocation[["target"]]
    if (!target %in% c("block", "sd_component")) {
      next
    }
    if (identical(target, "block") && length(allocation[["terms"]]) <= 1L) {
      next
    }
    if (identical(target, "sd_component") &&
        length(allocation[["leaf_terms"]]) <= 1L) {
      next
    }

    entries[[length(entries) + 1L]] <- list(
      allocation = allocation,
      name       = .brma_mv_allocation_summary_name(
        allocation     = allocation,
        formula_design = formula_design
      )
    )
  }

  plan           <- .brma_mv_allocation_summary_plan(entries, component)
  source_samples <- NULL
  out            <- list()

  for (i in plan[["selected"]]) {
    allocation <- entries[[i]][["allocation"]]
    if (is.null(source_samples)) {
      source_samples <- .brma_mv_scale_source_samples(
        object = object, posterior_samples = posterior_samples
      )
    }
    out[[length(out) + 1L]] <- .brma_mv_allocation_summary_samples(
      allocation        = allocation,
      formula_design    = formula_design,
      posterior_samples = posterior_samples,
      source_samples    = .brma_mv_allocation_observed_sources(
        allocation        = allocation,
        posterior_samples = posterior_samples,
        K                 = K,
        source_samples    = source_samples
      ),
      K                 = K
    )
  }

  available <- c(unlist(plan[["aliases"]], use.names = FALSE), plan[["names"]])

  return(.brma_mv_allocation_summary_result(out, plan, entries, available))
}


# The summary names of the allocations are made unique over all of them, in
# declaration order, before 'component' selects among them, so that a selector
# matches the names of the complete result. An exact name selects that
# allocation alone; otherwise every allocation declaring the alias is kept.
.brma_mv_allocation_summary_plan <- function(entries, component) {

  summary_names <- vapply(entries, `[[`, character(1), "name")
  unique_names  <- make.unique(summary_names, sep = "_")
  aliases <- lapply(seq_along(entries), function(i) {
    .brma_mv_allocation_summary_aliases(summary_names[[i]], entries[[i]][["allocation"]])
  })
  names(aliases) <- unique_names

  selected <- if (identical(component, "all")) {
    seq_along(entries)
  } else if (component %in% unique_names) {
    match(component, unique_names)
  } else {
    which(vapply(aliases, function(x) component %in% x, logical(1)))
  }

  list(names = unique_names, aliases = aliases, selected = selected)
}


.brma_mv_allocation_summary_result <- function(samples, plan, entries,
                                               available) {

  attr(samples, "available_components") <- unique(available)
  if (length(samples) == 0L) {
    return(samples)
  }

  selected <- plan[["selected"]]
  names(samples) <- plan[["names"]][selected]
  attr(samples, "component_aliases")      <- plan[["aliases"]][selected]
  attr(samples, "component_replacements") <- stats::setNames(
    lapply(entries[selected], function(entry) {
      .brma_mv_allocation_summary_replacements(entry[["allocation"]])
    }),
    plan[["names"]][selected]
  )

  return(samples)
}


# The source samples with the SD source of 'allocation' replaced by its observed
# values: evaluated on the posterior draws and scaled by each declared parent
# factor (weights and inclusion gates), as its summaries and its pooled SD use it.
.brma_mv_allocation_observed_sources <- function(allocation, posterior_samples,
                                                 K, source_samples = NULL) {

  observed <- .random_sd_source_samples(
    source            = allocation[["source"]],
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples
  )
  for (factor in allocation[["parent_factors"]]) {
    observed <- .recycle_rows(
      observed,
      .random_allocation_factor_samples(factor, posterior_samples)
    )
  }
  if (is.null(source_samples)) {
    source_samples <- list()
  }
  source_samples[[allocation[["source"]][["name"]]]] <- observed

  return(source_samples)
}


# The random blocks that an SD-component allocation replaces in the complete
# summary.
.brma_mv_allocation_summary_replacements <- function(allocation) {

  if (!identical(allocation[["target"]], "sd_component")) {
    return(character(0))
  }
  terms <- as.character(allocation[["terms"]])

  terms[!is.na(terms) & nzchar(terms)]
}


.brma_mv_allocation_catalog_sample_lists <- function(object,
                                                      posterior_samples,
                                                      component = "all") {

  formula_design <- .fitted_formula_design(object, "mu", required = TRUE)
  allocations    <- formula_design[["random_allocations"]]
  if (length(allocations) == 0L) {
    return(list())
  }

  catalog    <- BayesTools::parameter_catalog(object[["fit"]])
  quantities <- catalog[["quantities"]]
  quantities <- quantities[
    startsWith(quantities[["role"]], "random_") &
      !quantities[["internal"]] &
      quantities[["status"]] != "unavailable",
    ,
    drop = FALSE
  ]
  keys <- quantities[["extraction_key"]]
  entries <- list()
  available_components <- character()

  for (allocation in allocations) {
    allocation <- .brma_mv_resolve_sd_component_allocation(
      allocation     = allocation,
      formula_design = formula_design
    )
    target <- allocation[["target"]]
    if (!target %in% c("block", "sd_component")) {
      next
    }
    if (identical(target, "block") &&
        length(allocation[["terms"]]) <= 1L) {
      next
    }
    if (identical(target, "sd_component") &&
        length(allocation[["leaf_terms"]]) <= 1L) {
      next
    }
    name <- .brma_mv_allocation_summary_name(
      allocation = allocation, formula_design = formula_design
    )
    available_components <- c(
      available_components,
      .brma_mv_allocation_summary_aliases(name, allocation)
    )

    allocation_rows <- vapply(keys, function(key) {
      identical(key[["allocation_label"]], allocation[["label"]])
    }, logical(1))
    component_rows <- rep(FALSE, nrow(quantities))
    if (identical(target, "sd_component")) {
      block <- .brma_mv_sd_component_allocation_block(allocation)
      block_rows <- quantities[["owner_type"]] == "random_block" &
        vapply(keys, function(key) {
          identical(key[["random_block"]], block)
        }, logical(1))
      allocation_derived <- vapply(keys, function(key) {
        isTRUE(key[["allocation_derived"]])
      }, logical(1))
      component_rows <- block_rows & (
        quantities[["quantity"]] == "cor" |
          (quantities[["quantity"]] %in% c("sd", "var", "sd_mult") &
            allocation_derived)
      )
    }
    selected <- quantities[allocation_rows | component_rows, , drop = FALSE]
    if (nrow(selected) == 0L) {
      next
    }
    order_group <- match(
      selected[["quantity"]],
      c(
        "sd_total", "sd_common", "var_total", "var_common", "sd", "var",
        "cor", "sd_mult", "var_prop", "var_mult", "inclusion"
      )
    )

    entries[[length(entries) + 1L]] <- list(
      allocation = allocation,
      name       = name,
      rows       = selected[order(order_group, na.last = TRUE), , drop = FALSE]
    )
  }

  plan           <- .brma_mv_allocation_summary_plan(entries, component)
  source_samples <- NULL
  out            <- list()

  for (i in plan[["selected"]]) {
    allocation <- entries[[i]][["allocation"]]
    target     <- allocation[["target"]]
    selected   <- entries[[i]][["rows"]]

    samples_list <- lapply(seq_len(nrow(selected)), function(j) {
      selection <- BayesTools::parameter_catalog_resolve(
        catalog   = catalog,
        alias     = selected[["canonical_name"]][j],
        namespace = selected[["namespace"]][j]
      )
      draws <- BayesTools::parameter_draws(
        object[["fit"]],
        selection,
        model_samples = posterior_samples
      )
      values <- .brma_mv_catalog_draw_values(draws)
      attr(values, "catalog_name") <- selected[["canonical_name"]][j]
      values
    })
    names(samples_list) <- .brma_random_parameter_io_labels(selected, "label")
    if (identical(target, "sd_component") &&
        identical(allocation[["source"]][["shape"]], "row")) {
      if (is.null(source_samples)) {
        source_samples <- .brma_mv_scale_source_samples(
          object = object, posterior_samples = posterior_samples
        )
      }
      observed_sources <- .brma_mv_allocation_observed_sources(
        allocation        = allocation,
        posterior_samples = posterior_samples,
        K                 = nrow(object[["data"]][["outcome"]]),
        source_samples    = source_samples
      )
      block_rows <- which(selected[["owner_type"]] == "random_block")
      block_owner <- if (length(block_rows) > 0L) {
        selected[["label_parts"]][[block_rows[1L]]][["random"]][["owner"]]
      } else {
        .brma_mv_allocation_public_name(allocation)
      }
      observed <- .brma_mv_sd_component_allocation_summary_samples(
        allocation = allocation,
        term = .brma_mv_sd_component_allocation_term(allocation, formula_design),
        block_owner = block_owner, posterior_samples = posterior_samples,
        source_samples = observed_sources,
        K = nrow(object[["data"]][["outcome"]])
      )
      samples_list <- c(observed[!names(observed) %in% names(samples_list)],
        samples_list)
    }

    out[[length(out) + 1L]] <- samples_list
  }

  available <- c(available_components, plan[["names"]])

  return(.brma_mv_allocation_summary_result(out, plan, entries, available))
}



.brma_mv_resolve_sd_component_allocation <- function(allocation,
                                                    formula_design) {

  if (!identical(allocation[["target"]], "sd_component")) {
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
    blocks <- unique(vapply(
      formula_design[["random_effects"]],
      `[[`,
      character(1),
      "block_name"
    ))
    term <- .brma_mv_sd_component_allocation_term(allocation, formula_design)
    block_owner <- term[["block_name"]]
    if (is.null(block_owner) || !nzchar(block_owner)) {
      block_owner <- term[["group_label"]]
    } else if (identical(block_owner, term[["component_label"]]) &&
               !is.null(term[["component"]]) && nzchar(term[["component"]])) {
      block_owner <- term[["component"]]
    }
    if (length(blocks) <= 1L && !isTRUE(term[["component_visible"]])) {
      block_owner <- ""
    }
    return(.brma_mv_sd_component_allocation_summary_samples(
      allocation        = allocation,
      term              = term,
      block_owner       = block_owner,
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

  weight_name <- allocation[["weight_name"]]
  terms       <- allocation[["terms"]]
  labels      <- names(terms)
  if (is.null(labels)) {
    labels <- as.character(terms)
  }
  component_names <- allocation[["component_names"]]
  if (is.character(component_names) &&
      length(component_names) == length(labels) &&
      !anyNA(component_names) && all(nzchar(component_names))) {
    labels <- unname(component_names)
  }
  component_weights <- matrix(NA_real_, nrow(posterior_samples), length(terms))
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
    component_weights[, i] <- weights
  }

  inclusion <- allocation[["inclusion"]]
  for (gate in inclusion) {
    index <- gate[["index"]]
    BayesTools::check_int(index, "allocation inclusion index", lower = 1,
      upper = length(terms), allow_NULL = FALSE)
    component_weights[, index] <- component_weights[, index] *
      .random_allocation_factor_samples(
        list(inclusion_name = gate[["indicator_name"]]), posterior_samples
      )[, 1L]
  }
  active_weight <- rowSums(component_weights)
  if (length(inclusion) > 0L) {
    total_variance <- total_variance * active_weight
  }

  parent_active <- rep(TRUE, nrow(posterior_samples))
  parent_gated <- FALSE
  for (factor in allocation[["parent_factors"]]) {
    if (!is.null(factor[["inclusion_name"]])) {
      parent_gated <- TRUE
      parent_active <- parent_active & .random_allocation_factor_samples(
        list(inclusion_name = factor[["inclusion_name"]]), posterior_samples
      )[, 1L] == 1
    }
  }
  gated <- length(inclusion) > 0L || parent_gated
  if (gated) {
    active <- parent_active & active_weight > 0
    component_weights[active, ] <- component_weights[active, , drop = FALSE] /
      active_weight[active]
    component_weights[!active, ] <- NA_real_
  }

  aggregate_quantities <- .brma_mv_allocation_aggregate_quantities(allocation)
  allocation_owner <- .brma_mv_allocation_public_name(allocation)
  samples <- list()
  samples[[.brma_mv_allocation_parameter_name(
    allocation_owner,
    .brma_random_parameter_io_quantity(aggregate_quantities[["sd"]])
  )]] <- sqrt(total_variance)
  samples[[.brma_mv_allocation_parameter_name(
    allocation_owner,
    .brma_random_parameter_io_quantity(aggregate_quantities[["var"]])
  )]] <- total_variance
  for (i in seq_along(terms)) {
    share <- component_weights[, i]
    if (gated) {
      BayesTools::posterior_metadata(share, "undefined_draws") <- "allocation_active"
    }
    samples[[.brma_mv_allocation_parameter_name(
      allocation_owner, paste0("tau2_prop(", labels[[i]], ")")
    )]] <- share
  }

  samples
}


.brma_mv_allocation_aggregate_quantities <- function(allocation) {

  if (identical(allocation[["scale"]], "mean_variance")) {
    return(c(sd = "sd_common", var = "var_common"))
  }

  c(sd = "sd_total", var = "var_total")
}


.brma_mv_sd_component_allocation_summary_samples <- function(allocation,
                                                            term = NULL,
                                                            block_owner = NULL,
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
  block <- .brma_mv_sd_component_allocation_block(allocation)
  if (is.null(block_owner)) {
    block_owner <- block
  }

  aggregate_quantities <- .brma_mv_allocation_aggregate_quantities(allocation)
  allocation_owner <- .brma_mv_allocation_public_name(allocation)
  if (!nzchar(allocation_owner)) {
    allocation_owner <- block_owner
  }
  samples <- list()
  samples[[.brma_mv_allocation_parameter_name(
    allocation_owner,
    .brma_random_parameter_io_quantity(aggregate_quantities[["sd"]])
  )]] <- total_sd
  samples[[.brma_mv_allocation_parameter_name(
    allocation_owner,
    .brma_random_parameter_io_quantity(aggregate_quantities[["var"]])
  )]] <- total_sd^2

  weight_name <- allocation[["weight_name"]]
  labels      <- .brma_mv_sd_component_allocation_labels(allocation)
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

    variance_multiplier <- .brma_mv_allocation_variance_multiplier(
      weights   = weights,
      scale     = scale,
      n_targets = n_targets
    )
    if (!identical(source_shape, "row")) {
      sd_samples <- total_sd * sqrt(variance_multiplier)
    } else {
      rows <- .brma_mv_sd_component_allocation_rows(
        allocation   = allocation,
        term         = term,
        component    = i,
        source_shape = source_shape,
        K            = ncol(total)
      )
      sd_samples <- sqrt(rowMeans(total[, rows, drop = FALSE]^2) *
                           variance_multiplier)
    }

    samples[[.brma_mv_allocation_parameter_name(
      block_owner,
      paste0("tau(", labels[[i]], ")")
    )]] <- sd_samples
  }

  for (i in seq_along(labels)) {
    sd_name <- .brma_mv_allocation_parameter_name(
      block_owner,
      paste0("tau(", labels[[i]], ")")
    )
    samples[[.brma_mv_allocation_parameter_name(
      block_owner,
      paste0("tau2(", labels[[i]], ")")
    )]] <- samples[[sd_name]]^2
  }

  for (i in seq_along(labels)) {
    column <- paste0(weight_name, "[", i, "]")
    weights <- posterior_samples[, column]
    variance_multiplier <- .brma_mv_allocation_variance_multiplier(
      weights   = weights,
      scale     = scale,
      n_targets = n_targets
    )
    samples[[.brma_mv_allocation_parameter_name(
      allocation_owner,
      paste0("tau_mult(", labels[[i]], ")")
    )]] <- sqrt(variance_multiplier)
  }

  for (i in seq_along(labels)) {
    column <- paste0(weight_name, "[", i, "]")
    weights <- posterior_samples[, column]
    samples[[.brma_mv_allocation_parameter_name(
      allocation_owner,
      paste0("tau2_mult(", labels[[i]], ")")
    )]] <- .brma_mv_allocation_variance_multiplier(
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
      "Row-varying SD-component allocation summaries are unavailable because ",
      "the design-column allocation does not map to observation-level SD samples. ",
      "Inspect the fitted random-effect design and allocation metadata.",
      call. = FALSE
    )
  }

  active <- rowSums(model_matrix[, columns, drop = FALSE] != 0) > 0
  rows   <- which(active)
  if (length(rows) == 0L) {
    stop(
      "Row-varying SD-component allocation summaries are unavailable because no ",
      "observations map to one of the components. ",
      "Inspect the fitted random-effect design and allocation metadata.",
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


.brma_mv_allocation_variance_multiplier <- function(weights, scale, n_targets) {

  if (identical(scale, "mean_variance")) {
    return(n_targets * weights)
  }
  if (identical(scale, "total_variance")) {
    return(weights)
  }

  stop("Random-effect allocation factor scale is unsupported.",
       call. = FALSE)
}


.brma_mv_allocation_public_name <- function(allocation) {

  display_name <- allocation[["display_name"]]
  if (is.character(display_name) && length(display_name) == 1L &&
      !is.na(display_name)) {
    return(display_name)
  }

  allocation[["label"]]
}


.brma_mv_allocation_parameter_name <- function(owner, quantity) {

  if (!nzchar(owner)) {
    return(quantity)
  }

  paste0(owner, ": ", quantity)
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
    .brma_mv_allocation_public_name(allocation)
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


.brma_mv_heterogeneity_components <- function(
    object, posterior_samples = NULL,
    include_known_group_covariance = FALSE) {

  posterior_samples <- .get_posterior_samples(object[["fit"]],
                                              posterior_samples)
  K <- nrow(object[["data"]][["outcome"]])

  if (.is_random(object)) {
    return(.brma_mv_random_heterogeneity_components(
      object                         = object,
      posterior_samples              = posterior_samples,
      K                              = K,
      include_known_group_covariance = include_known_group_covariance
    ))
  }

  tau_result <- .evaluate.brma.tau(
    fit               = .object_formula_fit(object),
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (.is_scale(object)) {
      .create_fit_formula_list(data = object[["data"]], "scale")
    } else {
      NULL
    },
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = .is_scale(object),
    is_multilevel     = FALSE,
    is_random         = .is_random(object),
    K                 = K,
    posterior_samples = posterior_samples,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )
  out <- list(estimate = .expand_brma_mv_heterogeneity_samples(
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
                                                     new_levels = NULL,
                                                     include_known_group_covariance = TRUE) {

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

  block_names     <- vapply(terms, `[[`, character(1), "block_name")
  out             <- vector("list", length(terms))
  fallback_blocks <- character()
  fallback_errors <- list()
  names(out)      <- block_names
  sd_evaluator    <- .marginalized_random_sd_evaluator(
    fit   = object[["fit"]],
    terms = terms[pure_intercepts]
  )

  for (i in seq_along(terms)) {
    term <- terms[[i]]
    if (pure_intercepts[[i]]) {
      direct <- tryCatch(
        .random_effect_term_sd_samples(
          term              = term,
          posterior_samples = posterior_samples,
          K                 = K,
          source_samples    = source_samples,
          sd_evaluator      = sd_evaluator
        ),
        error = identity
      )
      if (!inherits(direct, "error")) {
        out[[i]] <- .expand_brma_mv_heterogeneity_samples(
          samples = direct,
          S       = nrow(posterior_samples),
          K       = K
        )
        next
      }
      fallback_errors[[term[["block_name"]]]] <- direct
    }
    fallback_blocks <- c(fallback_blocks, term[["block_name"]])
  }

  if (length(fallback_blocks) > 0L) {
    fallback <- .brma_mv_random_block_row_sd_samples_batch(
      object                         = object,
      posterior_samples              = posterior_samples,
      blocks                          = fallback_blocks,
      K                               = K,
      data                            = data,
      new_levels                      = new_levels,
      original_errors                 = fallback_errors,
      include_known_group_covariance = include_known_group_covariance
    )
    for (block in fallback_blocks) {
      out[[block]] <- .expand_brma_mv_heterogeneity_samples(
        samples = fallback[[block]],
        S       = nrow(posterior_samples),
        K       = K
      )
    }
  }

  return(out)
}


.brma_mv_random_block_row_sd_samples_batch <- function(
    object, posterior_samples, blocks, K, data = object[["data"]],
    new_levels = NULL, original_errors = list(),
    include_known_group_covariance = TRUE) {

  if (is.null(new_levels)) {
    variance <- tryCatch(
      .brma_mv_random_effects_marginal_diagonal_by_block(
        object                         = object,
        posterior_samples              = posterior_samples,
        blocks                          = blocks,
        data                            = data,
        include_known_group_covariance = include_known_group_covariance
      ),
      error = identity
    )
    valid <- is.list(variance) && identical(names(variance), blocks) &&
      all(vapply(variance, function(samples) {
        is.numeric(samples) && is.matrix(samples) &&
          identical(dim(samples), c(nrow(posterior_samples), K)) &&
          all(is.finite(samples)) && all(samples >= 0)
      }, logical(1)))
    if (isTRUE(valid)) {
      return(lapply(variance, sqrt))
    }
    # Preserve the established block-specific error if batch compilation is
    # structurally unavailable.
  }

  out <- lapply(blocks, function(block) {
    .brma_mv_random_block_row_sd_samples(
      object                         = object,
      posterior_samples              = posterior_samples,
      block                          = block,
      K                              = K,
      data                           = data,
      new_levels                     = new_levels,
      original_error                 = original_errors[[block]],
      include_known_group_covariance = include_known_group_covariance
    )
  })
  names(out) <- blocks

  return(out)
}


.brma_mv_random_term_is_pure_intercept <- function(term) {

  if (.random_effect_term_has_known_group_covariance(term) ||
      .marginalized_random_effect_has_row_multiplier(term) ||
      !inherits(term[["term_formula"]], "formula")) {
    return(FALSE)
  }

  return(.random_effect_term_is_intercept_only(term))
}


.random_effect_term_sd_samples <- function(term, posterior_samples, K,
                                           source_samples = NULL,
                                           sd_evaluator = NULL) {

  sd_samples <- .marginalized_random_effect_sd_samples(
    term              = term,
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples,
    sd_evaluator      = sd_evaluator
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
    fit               = .object_formula_fit(object),
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
    new_levels = NULL, original_error = NULL,
    include_known_group_covariance = TRUE) {

  random_vcov <- tryCatch(
    .brma_mv_random_effects_marginal_vcov(
      object                         = object,
      posterior_samples              = posterior_samples,
      blocks                          = block,
      data                            = data,
      new_levels                      = new_levels,
      diagonal_only                   = TRUE,
      include_known_group_covariance = include_known_group_covariance
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
      paste(unique(c(component_names, extra_component_names, "all", "total")),
            collapse = ", "), ".",
      call. = FALSE
    )
  }

  stop(
    "Unknown heterogeneity component '", component, "'. Available components: ",
    paste(unique(c(component_names, extra_component_names, "all", "total")),
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

  if (!is.list(components) || length(components) == 0L ||
      any(vapply(components, function(samples) !is.numeric(samples) ||
                   length(samples) == 0L, logical(1L)))) {
    stop("Heterogeneity components must contain numeric samples.", call. = FALSE)
  }
  variance <- Reduce(`+`, lapply(components, function(samples) samples^2))
  if (any(!is.finite(variance)) || any(variance < 0)) {
    stop("Total heterogeneity requires finite nonnegative variance.",
         call. = FALSE)
  }
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
