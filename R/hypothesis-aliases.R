.hypothesis_brma_select_parameter <- function(object, hypothesis,
                                              component) {

  component <- .parameter_component_normalize(component)
  catalog   <- .brma_parameter_catalog(object)

  roots <- .hypothesis_brma_symbol_roots(hypothesis)
  roots <- unique(roots[nzchar(roots)])
  if (length(roots) == 0L) {
    stop("Hypothesis must reference a model parameter.", call. = FALSE)
  }

  matches <- lapply(roots, function(root) {
    rows <- catalog[catalog[["alias"]] == root, , drop = FALSE]
    if (!identical(component, "auto")) {
      rows <- rows[rows[["component"]] == component, , drop = FALSE]
    }
    rows
  })
  matches <- do.call(rbind, matches)
  if (is.null(matches) || nrow(matches) == 0L) {
    stop(
      "Could not infer a model parameter from the hypothesis. Available ",
      "quantities are: ", .brma_parameter_available(catalog, component),
      ".",
      call. = FALSE
    )
  }

  known <- matches[["parameter"]]
  known <- unique(known)

  if (length(known) == 1L) {
    aliases <- .hypothesis_brma_aliases_for_catalog_parameter(catalog, known)
    component <- unique(matches[["component"]][matches[["parameter"]] == known])
    return(list(
      parameter = known,
      aliases   = aliases,
      component = component[1L]
    ))
  }

  if (length(known) > 1L) {
    stop(
      "Hypothesis references multiple model parameters (",
      paste(known, collapse = ", "),
      "). Set 'component' to 'mods'/'location', 'scale', or 'bias'.",
      call. = FALSE
    )
  }
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

  rewritten <- vapply(hypothesis, function(text) {
    text <- .hypothesis_brma_normalize_level_references(text)
    symbols <- .hypothesis_brma_symbols_normalized(text)
    symbols <- symbols[order(nchar(symbols), decreasing = TRUE)]

    for (symbol in symbols) {
      ref <- .hypothesis_brma_level_ref(symbol)
      if (!is.null(ref)) {
        root <- ref[["parameter"]]
        if (!is.null(aliases[[root]]) && identical(aliases[[root]], parameter)) {
          replacement <- paste0(aliases[[root]], "[", ref[["level"]], "]")
          text <- gsub(
            paste0("`", symbol, "`"),
            paste0("`", replacement, "`"),
            text,
            fixed = TRUE
          )
        }
      } else if (!is.null(aliases[[symbol]]) &&
                 identical(aliases[[symbol]], parameter)) {
        text <- .hypothesis_brma_replace_symbol(
          text        = text,
          symbol      = symbol,
          replacement = aliases[[symbol]]
        )
      }
    }

    .hypothesis_brma_unquote_level_references(text)
  }, character(1), USE.NAMES = FALSE)

  return(rewritten)
}


.hypothesis_brma_unquote_level_references <- function(text) {

  gsub("`([^`]+\\[[^`]+\\])`", "\\1", text, perl = TRUE)
}


.hypothesis_brma_replace_symbol <- function(text, symbol, replacement) {

  n_delimiters <- nchar(text, type = "bytes") -
    nchar(gsub("`", "", text, fixed = TRUE), type = "bytes")
  pieces <- strsplit(text, "`", fixed = TRUE)[[1L]]
  expected_length <- n_delimiters + 1L
  if (length(pieces) < expected_length) {
    pieces <- c(pieces, rep("", expected_length - length(pieces)))
  }
  if (length(pieces) == 0L) {
    return(text)
  }

  for (i in seq_along(pieces)) {
    if (i %% 2L == 0L) {
      if (identical(pieces[[i]], symbol)) {
        pieces[[i]] <- replacement
      }
      next
    }

    pieces[[i]] <- gsub(
      paste0("(?<![A-Za-z0-9._])", BayesTools::JAGS_regex_escape(symbol),
             "(?![A-Za-z0-9._])"),
      replacement,
      pieces[[i]],
      perl = TRUE
    )
  }

  paste(pieces, collapse = "`")
}


.hypothesis_brma_normalize_level_references <- function(text) {

  return(BayesTools::hypothesis_normalize_level_references(text))
}


.hypothesis_brma_symbols <- function(hypothesis) {

  text <- .hypothesis_brma_normalize_level_references(hypothesis)
  .hypothesis_brma_symbols_normalized(text)
}


.hypothesis_brma_symbols_normalized <- function(text) {

  sides <- unlist(strsplit(text, "\\s+[Vv][Ss]\\s+", perl = TRUE),
                  use.names = FALSE)

  symbols <- character()
  for (side in sides) {
    parsed <- try(parse(text = side, keep.source = FALSE), silent = TRUE)
    if (inherits(parsed, "try-error")) {
      next
    }
    for (expr in parsed) {
      symbols <- c(symbols, all.names(expr, functions = TRUE, unique = TRUE))
    }
  }

  blocked <- c(
    "(", "+", "-", "*", "/", "^", "<", "<=", ">", ">=", "=", "==", "!=",
    "&", "|", "!", "abs", "exp", "log", "sqrt", "plogis", "qlogis",
    "Inf", "NaN", "NA", "TRUE", "FALSE"
  )

  unique(setdiff(symbols, blocked))
}


.hypothesis_brma_symbol_roots <- function(hypothesis) {

  symbols <- unique(unlist(
    lapply(hypothesis, .hypothesis_brma_symbols),
    use.names = FALSE
  ))
  roots <- vapply(symbols, function(symbol) {
    ref <- .hypothesis_brma_level_ref(symbol)
    if (is.null(ref)) symbol else ref[["parameter"]]
  }, character(1))

  return(roots)
}


.hypothesis_brma_level_ref <- function(symbol) {

  ref <- BayesTools::hypothesis_parse_level_reference(symbol)
  if (nrow(ref) != 1L || !isTRUE(ref[["direct"]][[1L]])) {
    return(NULL)
  }

  list(
    parameter = ref[["parameter"]][[1L]],
    level     = ref[["level"]][[1L]]
  )
}
