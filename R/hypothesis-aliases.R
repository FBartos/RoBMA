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
  # 'group_component' is the component with which an ambiguity between a
  # factor term and its contrast coefficients is resolved to the term
  # (.brma_parameter_catalog_group_ambiguity()).
  resolve <- function(resolver_component, group_component = component) {
    tryCatch(
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
          component    = group_component
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
        # Candidates of no several model parameters (e.g. contrast
        # coefficients without parameter entries): BayesTools' ambiguity,
        # with the classes of an ambiguous statement first.
        class(error) <- unique(c(.hypothesis_ambiguous_class(), class(error)))
        stop(error)
      }
    )
  }
  # BayesTools' refusals of the statement's references are statement errors:
  # BayesTools' condition, with its classes and fields, after
  # "RoBMA_hypothesis_statement" (a statement without parameter symbols,
  # "BayesTools_hypothesis_no_parameters"; a level of another component,
  # "BayesTools_hypothesis_component_mismatch"; a contrast selector of a
  # level, "BayesTools_selector_unavailable"; an unknown reference). A
  # reference that is unknown within an explicit 'component' (whose aliases
  # it is not) but names quantities without the component is resolved
  # without it: a name of a quantity of the component (e.g. its display
  # label, such as '(mu) intercept' for component = "mods") is accepted; a
  # name known only outside the component ('outside') is refused by the
  # checks below as a component mismatch (or as the target it is), also next
  # to references of the component. Every reference of the statement counts,
  # not only the unknown name that BayesTools reports first, so that a
  # display label of the component next to a name of another component is a
  # mismatch in either order. A reference unknown in every component is the
  # unknown name refused.
  outside  <- FALSE
  resolved <- tryCatch(
    resolve(resolver_component),
    BayesTools_parameter_not_found = function(error) {
      if (!is.null(resolver_component) && NROW(metadata[["entries"]]) > 0L) {
        unrestricted <- tryCatch(
          resolve(NULL, group_component = "auto"),
          BayesTools_parameter_not_found = .hypothesis_stop_statement_condition
        )
        inside <- vapply(
          unrestricted[["occurrences"]][["quantity_id"]],
          function(quantity_id) {
            own <- .brma_parameter_catalog_entries_for_quantities(
              entries      = metadata[["entries"]],
              quantity_ids = quantity_id
            )
            any(own[["component"]] == component)
          },
          logical(1)
        )
        outside <<- !all(inside)
        return(unrestricted)
      }
      .hypothesis_stop_statement_condition(error)
    },
    BayesTools_selector_unavailable          = .hypothesis_stop_statement_condition,
    BayesTools_hypothesis_no_parameters      = .hypothesis_stop_statement_condition,
    BayesTools_hypothesis_component_mismatch = .hypothesis_stop_statement_condition
  )
  # A parameter catalog without RoBMA entries is metadata of an older build.
  if (NROW(metadata[["entries"]]) == 0L) {
    .stop_refit_required(
      "Resolved hypothesis metadata are unavailable. Refit the model with ",
      "the current RoBMA/BayesTools build."
    )
  }
  quantity_ids <- unique(resolved[["occurrences"]][["quantity_id"]])
  if (length(quantity_ids) == 0L) {
    stop("Internal error: the hypothesis resolved to no catalog quantities.",
         call. = FALSE)
  }
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
  if (!all(covered)) {
    .hypothesis_brma_stop_uncovered(
      object     = object,
      metadata   = metadata,
      resolution = resolved,
      uncovered  = quantity_ids[!covered]
    )
  }
  # A statement whose parameters do not belong to the requested component is
  # a component mismatch of the statement, as for plot() selections; so is a
  # statement with a reference known only outside the component.
  if (!identical(component, "auto")) {
    compatible <- entries[["component"]] == component
    if (!any(compatible) || (outside && !all(compatible))) {
      selected_components <- unique(entries[["component"]][!compatible])
      if (length(selected_components) == 1L) {
        .parameter_component_check_compatible(
          component          = component,
          selected_component = selected_components,
          argument           = "hypothesis",
          class              = "RoBMA_hypothesis_statement"
        )
      }
      .stop_component_mismatch(
        paste0("The hypothesis does not resolve to component = '", component, "'."),
        class = "RoBMA_hypothesis_statement"
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
  # The roots with which the statement references the selected parameter:
  # every resolved occurrence of its quantity or of its levels (its aliases,
  # and display labels such as '(mu) intercept', 'exp(intercept)' or the
  # root '(mu) g' of '(mu) g[a]'), which .hypothesis_brma_rewrite() maps to
  # the parameter.
  occurrences <- resolved[["occurrences"]]
  own_ids     <- c(entry[["quantity_id"]],
                   unlist(entry[["member_quantity_ids"]], use.names = FALSE))
  roots       <- unique(occurrences[["parameter"]][
    occurrences[["quantity_id"]] %in% own_ids
  ])
  return(list(
    parameter  = entry[["parameter"]],
    aliases    = aliases,
    roots      = roots,
    component  = entry[["component"]],
    entry      = as.list(entry[1L, setdiff(names(entry), "aliases"), drop = FALSE]),
    resolution = resolved
  ))
}


# Resolved catalog quantities without a RoBMA parameter entry ('uncovered',
# quantity ids of 'resolution') are no hypothesis targets of a current fit,
# so a refit cannot help: contrast coefficients of formula factor terms are
# statements to restate on factor levels ("RoBMA_hypothesis_statement");
# quantities of the publication-bias prior (e.g. the weight-function
# coordinates 'omega[1]' and the mixture indicator 'bias_indicator') have the
# publication-bias refusal; other quantities (e.g. inclusion indicators,
# 'inclusion(<component>)' of variance allocations, and latent cluster
# effects) are refused as targets ("RoBMA_hypothesis_target"), named as the
# statement references them (.brma_stop_uncovered_quantity(), which the
# parameter selection of plot() shares).
.hypothesis_brma_stop_uncovered <- function(object, metadata, resolution,
                                            uncovered) {

  occurrences <- resolution[["occurrences"]]
  symbol      <- occurrences[["symbol"]][
    match(uncovered[[1L]], occurrences[["quantity_id"]])
  ]

  .brma_stop_uncovered_quantity(
    object       = object,
    metadata     = metadata,
    quantity_ids = uncovered,
    symbol       = symbol,
    hypothesis   = TRUE
  )
}


# The catalog quantities among 'quantity_ids' of the publication-bias prior:
# those whose fitted coordinates are all nodes that BayesTools monitors for
# that prior (BayesTools::JAGS_to_monitor(), e.g. 'omega', 'bias_indicator',
# 'PET', and 'PEESE').
.hypothesis_brma_bias_quantity_ids <- function(object, metadata,
                                               quantity_ids) {

  prior <- object[["priors"]][["outcome"]][["bias"]]
  if (is.null(prior) || length(quantity_ids) == 0L) {
    return(character())
  }
  monitors    <- BayesTools::JAGS_to_monitor(list(bias = prior))
  coordinates <- BayesTools::parameter_coordinates(object[["fit"]])
  bias_coordinates <- coordinates[["coordinate_name"]][
    coordinates[["monitor_name"]] %in% monitors
  ]
  quantities  <- metadata[["catalog"]][["quantities"]]
  rows        <- quantities[
    quantities[["quantity_id"]] %in% quantity_ids,
    ,
    drop = FALSE
  ]
  bias <- vapply(rows[["extraction_key"]], function(key) {
    dependencies <- if (is.list(key)) key[["dependencies"]]
    length(dependencies) > 0L && all(dependencies %in% bias_coordinates)
  }, logical(1))

  rows[["quantity_id"]][bias]
}


# A statement whose references name several model parameters is ambiguous
# ("RoBMA_hypothesis_ambiguous" with the parent "RoBMA_hypothesis_statement",
# .hypothesis_ambiguous_class()); 'component' selects one of them.
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
    .stop_component_mismatch(
      paste0("The hypothesis does not resolve to component = '", component, "'."),
      class = "RoBMA_hypothesis_statement"
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


# Publication-bias parameters are no hypothesis targets
# ("RoBMA_hypothesis_target").
.hypothesis_brma_check_supported_component <- function(component) {

  if (identical(component, "bias")) {
    .hypothesis_stop(.hypothesis_refusal(
      "Hypothesis tests for publication-bias parameters are not supported.",
      "target"
    ))
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


# The statement written against 'parameter': the roots that name it (its
# 'aliases', and the 'resolved_roots' of its quantities, which include
# display labels that are no alias of the parameter) are rewritten to it.
.hypothesis_brma_rewrite <- function(hypothesis, aliases, parameter,
                                     resolved_roots = character()) {

  ast <- .hypothesis_brma_ast(hypothesis)
  mapping <- unlist(aliases, use.names = TRUE)
  mapping <- mapping[mapping == parameter & names(mapping) != mapping]
  labels  <- setdiff(resolved_roots, c(names(mapping), parameter))
  if (length(labels) > 0L) {
    mapping <- c(mapping, stats::setNames(rep(parameter, length(labels)), labels))
  }
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
# selector. A contrast-coefficient selector of a level coordinate (e.g.
# 'g{1}' of a treatment factor), which BayesTools refuses naming the level
# form, is a statement to restate on that level, as the selectors of
# contrast coefficients are (.brma_stop_contrast_coefficient()).
.hypothesis_brma_ast <- function(hypothesis, catalog = NULL) {

  if (inherits(hypothesis, "BayesTools_hypothesis_ast")) {
    return(hypothesis)
  }
  return(tryCatch(
    BayesTools::hypothesis_parse(
      hypothesis     = hypothesis,
      catalog        = catalog,
      simplify_names = !is.null(catalog)
    ),
    BayesTools_selector_unavailable = .hypothesis_stop_statement_condition
  ))
}


# Re-raises a BayesTools refusal of a hypothesis statement (its message,
# classes, and fields) as a statement error, with "RoBMA_hypothesis_statement"
# before its classes.
.hypothesis_stop_statement_condition <- function(condition) {

  class(condition) <- unique(c("RoBMA_hypothesis_statement", class(condition)))
  stop(condition)
}
