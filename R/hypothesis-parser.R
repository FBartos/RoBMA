.hypothesis_brma_point_refs <- function(hypothesis, parameter,
                                        require_direct = TRUE) {

  refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  if (nrow(refs) == 0L) {
    return(data.frame(
      symbol           = character(),
      level            = character(),
      value            = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  if (isTRUE(require_direct) && any(!refs[["direct"]])) {
    unsupported <- refs[["hypothesis"]][!refs[["direct"]]][1L]
    stop(
      "Point-null hypotheses require a direct parameter or level reference; ",
      "unsupported point expression in: '",
      unsupported, "'.",
      call. = FALSE
    )
  }

  refs <- refs[refs[["direct"]] & refs[["parameter"]] == parameter, , drop = FALSE]
  if (nrow(refs) == 0L) {
    return(data.frame(
      symbol           = character(),
      level            = character(),
      value            = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- unique(refs[, c("symbol", "level", "value"), drop = FALSE])
  rownames(out) <- NULL
  return(out)
}


.hypothesis_brma_level_contrast_candidate <- function(hypothesis,
                                                       parameter) {

  refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  if (nrow(refs) == 0L || !any(!refs[["direct"]])) {
    return(FALSE)
  }

  occurrences <- BayesTools::hypothesis_symbols(
    hypothesis,
    occurrences = TRUE
  )
  if (nrow(occurrences) == 0L ||
      any(is.na(occurrences[["level"]])) ||
      any(occurrences[["parameter"]] != parameter)) {
    return(FALSE)
  }

  return(length(unique(occurrences[["level"]])) == 2L)
}
