.hypothesis_brma_point_refs <- function(hypothesis, parameter,
                                        require_direct = TRUE) {

  if (inherits(hypothesis, "BayesTools_hypothesis_ast")) {
    hypothesis <- BayesTools::hypothesis_render(hypothesis)
  }
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
      "qCMDE/IWMDE ordinates require direct parameter or ",
      "level point hypotheses; unsupported point expression in: '",
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


.hypothesis_brma_has_compound_point <- function(hypothesis) {

  refs <- BayesTools::hypothesis_parse_point_reference(
    hypothesis     = BayesTools::hypothesis_render(hypothesis),
    allow_compound = TRUE
  )
  nrow(refs) > 0L && any(!refs[["direct"]])
}
