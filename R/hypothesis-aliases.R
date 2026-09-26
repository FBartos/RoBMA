.hypothesis_brma_select_parameter <- function(object, hypothesis,
                                              component, metadata = NULL) {

  component          <- .parameter_component_normalize(component)
  if (is.null(metadata)) {
    metadata <- .brma_parameter_catalog_metadata(object)
  }
  ast                <- .hypothesis_brma_ast(hypothesis, metadata[["catalog"]])
  occurrences        <- BayesTools::hypothesis_symbols(
    ast,
    occurrences = TRUE
  )
  resolver_component <- if (identical(component, "auto") ||
                            any(!is.na(occurrences[["level"]]))) {
    NULL
  } else {
    component
  }
  resolved <- tryCatch(
    BayesTools::hypothesis_resolve(
      ast       = ast,
      catalog   = metadata[["catalog"]],
      component = resolver_component,
      simplify_names = TRUE
    ),
    BayesTools_parameter_ambiguous = function(error) {
      ambiguity <- .brma_parameter_catalog_group_ambiguity(
        entries      = metadata[["entries"]],
        quantity_ids = error[["candidates"]][["quantity_id"]],
        component    = component
      )
      if (!is.null(ambiguity[["component"]])) {
        return(BayesTools::hypothesis_resolve(
          ast       = ast,
          catalog   = metadata[["catalog"]],
          component = ambiguity[["component"]],
          simplify_names = TRUE
        ))
      }
      candidates <- ambiguity[["candidates"]]
      namespace  <- .hypothesis_brma_component_namespace(
        candidates = candidates,
        quantities = error[["candidates"]],
        component  = component
      )
      if (!is.null(namespace)) {
        return(BayesTools::hypothesis_resolve(
          ast       = ast,
          catalog   = metadata[["catalog"]],
          namespace = namespace,
          component = resolver_component,
          simplify_names = TRUE
        ))
      }
      if (length(unique(candidates[["parameter"]])) > 1L) {
        .hypothesis_brma_stop_multiple_parameters(candidates)
      }
      stop(error)
    }
  )
  quantity_ids <- unique(resolved[["occurrences"]][["quantity_id"]])
  entries <- .brma_parameter_catalog_entries_for_quantities(
    entries      = metadata[["entries"]],
    quantity_ids = quantity_ids
  )
  covered <- vapply(quantity_ids, function(quantity_id) {
    any(entries[["quantity_id"]] == quantity_id) ||
      any(vapply(
        entries[["member_quantity_ids"]],
        function(members) quantity_id %in% members,
        logical(1)
      ))
  }, logical(1))
  if (length(quantity_ids) == 0L || !all(covered) || nrow(entries) == 0L) {
    coefficients <- .brma_contrast_coefficient_quantities(
      metadata, quantity_ids[!covered]
    )
    if (nrow(coefficients) > 0L) {
      .brma_stop_contrast_coefficient(
        metadata    = metadata,
        selector    = coefficients[["canonical_name"]][[1L]],
        coefficient = coefficients[1L, , drop = FALSE]
      )
    }
    stop(
      "Resolved hypothesis metadata are unavailable. Refit the model with ",
      "the current RoBMA/BayesTools build.",
      call. = FALSE
    )
  }
  if (!identical(component, "auto")) {
    compatible <- entries[["component"]] == component
    if (!any(compatible)) {
      selected_components <- unique(entries[["component"]])
      if (length(selected_components) == 1L) {
        .parameter_component_check_compatible(
          component          = component,
          selected_component = selected_components,
          argument           = "hypothesis"
        )
      }
      stop(
        "The hypothesis does not resolve to component = '", component, "'.",
        call. = FALSE
      )
    }
    entries <- entries[compatible, , drop = FALSE]
  }
  if (nrow(entries) > 1L) {
    .hypothesis_brma_stop_multiple_parameters(entries)
  }

  entry   <- entries[1L, , drop = FALSE]
  aliases <- as.list(rep(entry[["parameter"]], length(entry[["aliases"]][[1L]])))
  names(aliases) <- entry[["aliases"]][[1L]]
  aliases[[entry[["parameter"]]]] <- entry[["parameter"]]
  return(list(
    parameter  = entry[["parameter"]],
    aliases    = aliases,
    component  = entry[["component"]],
    entry      = as.list(entry[1L, setdiff(names(entry), "aliases"), drop = FALSE]),
    resolution = resolved
  ))
}


# A statement whose references name several model parameters is ambiguous
# ("RoBMA_hypothesis_ambiguous" with the parent of the hypothesis refusals,
# .hypothesis_refusal()); 'component' selects one of them.
.hypothesis_brma_stop_multiple_parameters <- function(entries) {

  .hypothesis_stop(.hypothesis_refusal(
    paste0(
      "Hypothesis references multiple model parameters (",
      paste(unique(entries[["parameter"]]), collapse = ", "),
      "). Set 'component' to 'mods'/'location', 'scale', or 'bias'."
    ),
    "ambiguous"
  ))
}


# Level references are resolved without the requested component: the catalog
# component of a level quantity is its level. When a level reference names
# levels of several parameters (e.g. the alias of a factor term of both the
# location and the scale formula), an explicit component restricts them to
# its own parameter, as the resolver restricts unleveled references: the
# namespace of that parameter's candidate quantities, in which the
# hypothesis is resolved again. NULL when the component is "auto" or does not
# select one parameter with one namespace; stops when no candidate belongs to
# the component. 'candidates' are the catalog entries of the ambiguous
# quantities 'quantities'.
.hypothesis_brma_component_namespace <- function(candidates, quantities,
                                                 component) {

  if (identical(component, "auto") ||
      length(unique(candidates[["parameter"]])) < 2L) {
    return(NULL)
  }
  selected <- candidates[candidates[["component"]] == component, , drop = FALSE]
  if (nrow(selected) == 0L) {
    stop(
      "The hypothesis does not resolve to component = '", component, "'.",
      call. = FALSE
    )
  }
  if (nrow(selected) != 1L) {
    return(NULL)
  }
  members   <- c(
    selected[["quantity_id"]],
    unlist(selected[["member_quantity_ids"]], use.names = FALSE)
  )
  namespace <- unique(quantities[["namespace"]][
    quantities[["quantity_id"]] %in% members
  ])
  if (length(namespace) != 1L) {
    return(NULL)
  }

  namespace
}


.hypothesis_brma_check_supported_component <- function(component) {

  if (identical(component, "bias")) {
    stop("Hypothesis tests for publication-bias parameters are not supported.",
         call. = FALSE)
  }

  invisible(TRUE)
}


.hypothesis_brma_alias_label <- function(aliases, parameter) {

  keep <- vapply(aliases, identical, logical(1), y = parameter)
  candidates <- names(aliases)[keep]
  candidates <- candidates[nzchar(candidates)]
  candidates <- setdiff(candidates, parameter)
  if (length(candidates) > 0L) {
    return(candidates[1L])
  }

  return(parameter)
}


.hypothesis_brma_rewrite <- function(hypothesis, aliases, parameter) {

  ast <- .hypothesis_brma_ast(hypothesis)
  mapping <- unlist(aliases, use.names = TRUE)
  mapping <- mapping[mapping == parameter & names(mapping) != mapping]
  statements <- BayesTools::hypothesis_render(ast)
  if (length(statements) == 1L) {
    roots   <- BayesTools::hypothesis_symbols(ast)
    mapping <- mapping[names(mapping) %in% roots]
    return(BayesTools::hypothesis_rewrite(ast, mapping))
  }

  rewritten <- vapply(statements, function(statement) {

    statement_ast <- BayesTools::hypothesis_parse(statement)
    roots <- BayesTools::hypothesis_symbols(statement_ast)
    statement_mapping <- mapping[names(mapping) %in% roots]
    BayesTools::hypothesis_render(
      BayesTools::hypothesis_rewrite(statement_ast, statement_mapping)
    )
  }, character(1))
  return(BayesTools::hypothesis_parse(unname(rewritten)))
}


.hypothesis_brma_symbol_roots <- function(hypothesis) {

  return(BayesTools::hypothesis_symbols(.hypothesis_brma_ast(hypothesis)))
}


# Hypotheses are parsed against the fitted parameter catalog, whose aliases
# include every rendered table label, level cell and contrast coefficient
# selector.
.hypothesis_brma_ast <- function(hypothesis, catalog = NULL) {

  if (inherits(hypothesis, "BayesTools_hypothesis_ast")) {
    return(hypothesis)
  }
  return(BayesTools::hypothesis_parse(
    hypothesis     = hypothesis,
    catalog        = catalog,
    simplify_names = !is.null(catalog)
  ))
}
