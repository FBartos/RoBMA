# ============================================================================ #
# Random-component inclusion metadata
# ============================================================================ #


.random_gated_root_allocations <- function(allocations) {

  if (length(allocations) == 0L) {
    return(list())
  }

  gated <- vapply(allocations, function(allocation) {
    is.null(allocation[["parent"]]) &&
      identical(allocation[["target"]], "block") &&
      length(allocation[["inclusion"]]) > 0L
  }, logical(1))

  allocations[gated]
}


.random_component_inclusion_map <- function(object) {

  formula_design <- .fitted_formula_design(object, "mu", required = FALSE)
  allocations    <- formula_design[["random_allocations"]]
  if (length(allocations) == 0L) {
    return(list())
  }

  allocation_labels <- vapply(allocations, `[[`, character(1), "label")
  allocations <- stats::setNames(allocations, allocation_labels)
  gated_roots <- .random_gated_root_allocations(allocations)
  out <- list()

  add_component <- function(allocation_label, component_id, indicator,
                            visited = character()) {

    visit_key <- paste(allocation_label, component_id, sep = "::")
    if (visit_key %in% visited) {
      stop("Random-effect allocation metadata contain a parent cycle.",
           call. = FALSE)
    }
    allocation <- allocations[[allocation_label]]
    if (is.null(allocation)) {
      stop(
        "Random-effect allocation metadata reference unknown parent '",
        allocation_label, "'.",
        call. = FALSE
      )
    }
    component_index <- match(component_id, names(allocation[["terms"]]))
    if (is.na(component_index)) {
      stop(
        "Random-effect allocation metadata reference unknown component '",
        component_id, "' in allocation '", allocation_label, "'.",
        call. = FALSE
      )
    }
    aliases <- c(
      component_id,
      unname(allocation[["terms"]][[component_id]]),
      allocation[["component_names"]][[component_index]]
    )
    aliases <- unique(aliases[!is.na(aliases) & nzchar(aliases)])
    for (alias in aliases) {
      out[[alias]] <<- unique(c(out[[alias]], indicator))
    }

    visited <- c(visited, visit_key)
    for (child in allocations) {
      parent <- child[["parent"]]
      if (is.null(parent) ||
          !identical(parent[["allocation"]], allocation_label) ||
          !identical(parent[["component"]], component_id)) {
        next
      }
      for (child_component in names(child[["terms"]])) {
        add_component(
          allocation_label = child[["label"]],
          component_id     = child_component,
          indicator        = indicator,
          visited          = visited
        )
      }
    }
  }

  for (allocation in gated_roots) {
    for (component_id in names(allocation[["inclusion"]])) {
      indicator <- allocation[["inclusion"]][[component_id]][[
        "indicator_name"
      ]]
      add_component(
        allocation_label = allocation[["label"]],
        component_id     = component_id,
        indicator        = indicator
      )
    }
  }

  out
}


# Pair each random-inclusion gate with the SD quantity it switches on, keyed by
# the gate's backend indicator parameter and valued by the public SD draw name.
# Summaries and `as_draws()` both report gates through this one mapping.
.random_inclusion_sd_names <- function(object) {

  quantities <- BayesTools::parameter_catalog(object[["fit"]])[["quantities"]]
  keys       <- quantities[["extraction_key"]]
  gate_rows  <- which(quantities[["role"]] == "random_inclusion")
  if (length(gate_rows) == 0L) {
    return(stats::setNames(character(), character()))
  }

  indicators    <- vapply(keys[gate_rows], `[[`, character(1), "source_parameter")
  sd_rows       <- which(quantities[["quantity"]] %in%
                           c("sd", "sd_total", "sd_common"))
  component_map <- .random_component_inclusion_map(object)
  sd_indicators <- lapply(sd_rows, function(row) {

    key <- keys[[row]]
    if (identical(quantities[["quantity"]][[row]], "sd")) {
      owners <- c(key[["random_block"]], quantities[["owner_name"]][[row]])
      return(unique(unlist(component_map[intersect(owners, names(component_map))],
                           use.names = FALSE)))
    }
    intersect(key[["dependencies"]], indicators)
  })

  paired <- vapply(seq_along(gate_rows), function(i) {

    candidates <- sd_rows[vapply(sd_indicators, identical, logical(1),
                                 indicators[[i]])]
    totals <- candidates[quantities[["quantity"]][candidates] %in%
                           c("sd_total", "sd_common")]
    if (length(totals) == 1L) {
      candidates <- totals
    }
    if (length(candidates) != 1L) {
      return(NA_character_)
    }
    .brma_random_parameter_io_name(
      sub("^\\([^)]*\\) ", "", quantities[["display_label"]][[candidates]]),
      quantities[["quantity"]][[candidates]]
    )
  }, character(1))

  out <- stats::setNames(paired, indicators)
  out[!is.na(out)]
}


# Pair each allocation-SD mixture prior with the aggregate SD it is a prior on,
# keyed by the mixture's backend indicator and valued by the public SD draw
# name. This indicator selects which mixture component allocated the total SD,
# so it takes one value per component rather than a 0/1 inclusion state. It
# sits on the shared SD budget rather than on a single random-effect component,
# so it is not a `random_inclusion` catalog row and is resolved separately from
# `.random_inclusion_sd_names()`.
.random_slab_sd_names <- function(object) {

  prior_list <- attr(object[["fit"]], "prior_list", exact = TRUE)
  slab       <- .random_slab_prior_parameters(prior_list)
  if (length(slab) == 0L) {
    return(stats::setNames(character(), character()))
  }

  quantities <- BayesTools::parameter_catalog(object[["fit"]])[["quantities"]]
  keys       <- quantities[["extraction_key"]]
  total_rows <- which(quantities[["quantity"]] %in% c("sd_total", "sd_common"))
  sd_rows    <- which(quantities[["quantity"]] == "sd")

  owning_row <- function(parameter) {

    for (rows in list(total_rows, sd_rows)) {
      candidates <- rows[vapply(rows, function(row) {
        parameter %in% keys[[row]][["dependencies"]]
      }, logical(1))]
      if (length(candidates) == 1L) {
        return(candidates)
      }
    }
    NA_integer_
  }

  paired <- vapply(slab, function(parameter) {

    row <- owning_row(parameter)
    if (is.na(row)) {
      return(NA_character_)
    }
    .brma_random_parameter_io_name(
      sub("^\\([^)]*\\) ", "", quantities[["display_label"]][[row]]),
      quantities[["quantity"]][[row]]
    )
  }, character(1))

  out <- stats::setNames(unname(paired), paste0(slab, "_indicator"))
  out[!is.na(out)]
}


.random_component_conditioning_parameters <- function(object, components,
                                                        fallback = NULL) {

  inclusion_parameters <- .random_component_inclusion_map(object)
  if (length(inclusion_parameters) == 0L && !is.null(fallback)) {
    inclusion_parameters <- stats::setNames(
      rep(list(fallback), length(components)),
      components
    )
  }
  all_parameters <- unique(unlist(inclusion_parameters, use.names = FALSE))

  out <- lapply(components, function(component) {
    parameters <- if (identical(component, "total")) {
      all_parameters
    } else {
      inclusion_parameters[[component]]
    }
    if (is.null(parameters) || length(parameters) == 0L) {
      stop(
        "No inclusion gate is available for random-effect component '",
        component, "'.",
        call. = FALSE
      )
    }
    parameters
  })
  names(out) <- components

  out
}


.random_allocation_inclusion_prior <- function(prior_list, indicator_name,
                                               required = FALSE) {

  if (!is.character(indicator_name) || length(indicator_name) != 1L ||
      is.na(indicator_name) || !nzchar(indicator_name)) {
    if (required) {
      stop("Random-effect inclusion metadata have no valid indicator name.",
           call. = FALSE)
    }
    return(NULL)
  }

  indicators <- .random_allocation_inclusion_indicators(prior_list)
  matches    <- indicators == indicator_name
  if (sum(matches) != 1L) {
    if (required) {
      stop(
        "Random-effect inclusion indicator '", indicator_name,
        "' has no unique prior probability.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  prior_list[[which(matches)]]
}


.random_allocation_inclusion_indicators <- function(prior_list) {

  indicators <- vapply(prior_list, function(prior) {
    indicator <- attr(prior, "random_allocation_indicator", exact = TRUE)
    if (is.character(indicator) && length(indicator) == 1L &&
        !is.na(indicator) && nzchar(indicator)) {
      indicator
    } else {
      ""
    }
  }, character(1))

  indicators
}


.random_slab_prior_parameters <- function(prior_list) {

  parameters <- names(prior_list)
  parameters[vapply(prior_list, function(prior) {
    isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)) &&
      BayesTools::is.prior.mixture(prior)
  }, logical(1))]
}
