.hypothesis_brma_select_parameter <- function(object, hypothesis,
                                              component) {

  component <- .parameter_component_normalize(component)
  metadata  <- .brma_parameter_catalog_metadata(object)
  ast       <- .hypothesis_brma_ast(hypothesis)
  resolved  <- BayesTools::hypothesis_resolve(
    ast       = ast,
    catalog   = metadata[["catalog"]],
    component = if (identical(component, "auto")) NULL else component
  )
  quantity_ids <- unique(resolved[["occurrences"]][["quantity_id"]])
  entries <- metadata[["entries"]][
    metadata[["entries"]][["quantity_id"]] %in% quantity_ids,
    ,
    drop = FALSE
  ]
  if (nrow(entries) != length(quantity_ids)) {
    stop(
      "Resolved hypothesis metadata are unavailable. Refit the model with ",
      "the current RoBMA/BayesTools build.",
      call. = FALSE
    )
  }
  if (length(quantity_ids) > 1L) {
    stop(
      "Hypothesis references multiple model parameters (",
      paste(unique(entries[["parameter"]]), collapse = ", "),
      "). Set 'component' to 'mods'/'location', 'scale', or 'bias'.",
      call. = FALSE
    )
  }

  entry   <- entries[1L, , drop = FALSE]
  aliases <- as.list(rep(entry[["parameter"]], length(entry[["aliases"]][[1L]])))
  names(aliases) <- entry[["aliases"]][[1L]]
  aliases[[entry[["parameter"]]]] <- entry[["parameter"]]
  return(list(
    parameter  = entry[["parameter"]],
    aliases    = aliases,
    component  = entry[["component"]],
    resolution = resolved
  ))
}


.hypothesis_brma_check_supported_component <- function(component) {

  if (identical(component, "bias")) {
    stop("Hypothesis tests for publication-bias parameters are not supported.",
         call. = FALSE)
  }

  invisible(TRUE)
}


.hypothesis_brma_aliases <- function(object) {

  catalog <- .brma_parameter_catalog(object)
  .hypothesis_brma_aliases_for_catalog_parameter(
    catalog   = catalog,
    parameter = unique(catalog[["parameter"]])
  )
}


.hypothesis_brma_aliases_for_catalog_parameter <- function(catalog,
                                                           parameter) {

  keep <- catalog[["parameter"]] %in% parameter
  rows <- catalog[keep, , drop = FALSE]

  aliases <- as.list(rows[["parameter"]])
  names(aliases) <- rows[["alias"]]

  for (name in parameter) {
    aliases[[name]] <- name
  }

  return(aliases)
}


.hypothesis_brma_aliases_for_parameter <- function(aliases, parameter) {

  keep <- vapply(aliases, identical, logical(1), y = parameter)
  aliases <- aliases[keep]
  aliases[[parameter]] <- parameter

  return(aliases)
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


.hypothesis_brma_available_aliases <- function(aliases) {

  available <- unique(names(aliases))
  available <- available[nzchar(available)]
  available <- sort(available)

  paste0("'", available, "'", collapse = ", ")
}


.hypothesis_brma_rewrite <- function(hypothesis, aliases, parameter) {

  ast <- .hypothesis_brma_ast(hypothesis)
  mapping <- unlist(aliases, use.names = TRUE)
  mapping <- mapping[mapping == parameter & names(mapping) != mapping]
  roots   <- BayesTools::hypothesis_symbols(ast)
  mapping <- mapping[names(mapping) %in% roots]
  return(BayesTools::hypothesis_rewrite(ast, mapping))
}


.hypothesis_brma_symbol_roots <- function(hypothesis) {

  return(BayesTools::hypothesis_symbols(.hypothesis_brma_ast(hypothesis)))
}


.hypothesis_brma_ast <- function(hypothesis) {

  if (inherits(hypothesis, "BayesTools_hypothesis_ast")) {
    return(hypothesis)
  }
  return(BayesTools::hypothesis_parse(hypothesis))
}
