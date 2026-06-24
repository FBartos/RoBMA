# ============================================================================ #
# brma-mv-known-r.R
# ============================================================================ #
#
# RoBMA-facing adapter from metafor-style R/Rscale arguments to BayesTools
# known group covariance metadata for formula random effects.
#
# ============================================================================ #


.brma_mv_has_R <- function(R) {

  !is.null(R)
}


.brma_mv_make_group_covariance <- function(R, Rscale) {

  if (!.brma_mv_has_R(R)) {
    return(NULL)
  }

  .brma_mv_require_group_covariance_api()

  entries <- .brma_mv_R_entries(R)
  scales  <- .brma_mv_normalize_Rscale(
    Rscale  = Rscale,
    R_names = names(entries)
  )

  group_covariance <- Map(
    function(covariance, scale) {
      BayesTools::random_group_covariance(
        covariance = covariance,
        scale      = scale
      )
    },
    entries,
    scales
  )

  entry_names <- names(entries)
  if (length(group_covariance) == 1L && !nzchar(entry_names[[1L]])) {
    return(group_covariance[[1L]])
  }

  names(group_covariance) <- entry_names
  return(group_covariance)
}


.brma_mv_require_group_covariance_api <- function() {

  has_constructor <- exists(
    "random_group_covariance",
    envir    = asNamespace("BayesTools"),
    inherits = FALSE
  )
  has_parser_arg <- "group_covariance" %in%
    names(formals(BayesTools::random_effects_formula))

  if (!has_constructor || !has_parser_arg) {
    stop(
      "Known random-effect group covariance requires a BayesTools build with ",
      "random_group_covariance() and random_effects_formula(..., ",
      "group_covariance = ...).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.brma_mv_group_covariance_error_message <- function(message) {

  gsub("'group_covariance'", "'R'", message, fixed = TRUE)
}


.brma_mv_R_entries <- function(R) {

  if (is.matrix(R) || is.data.frame(R)) {
    out <- list(R)
    names(out) <- ""
    return(out)
  }

  if (!is.list(R) || length(R) == 0L) {
    stop(
      "'R' must be a matrix or a non-empty list of matrices.",
      call. = FALSE
    )
  }

  R_names <- names(R)
  if (is.null(R_names)) {
    R_names <- rep("", length(R))
  }
  R_names[is.na(R_names)] <- ""

  if (length(R) > 1L && any(!nzchar(R_names))) {
    stop(
      "'R' must be fully named when multiple covariance matrices are supplied.",
      call. = FALSE
    )
  }
  if (anyDuplicated(R_names[nzchar(R_names)])) {
    stop("'R' names must be unique.", call. = FALSE)
  }

  names(R) <- R_names
  return(R)
}


.brma_mv_normalize_Rscale <- function(Rscale, R_names) {

  n_R <- length(R_names)
  if (n_R < 1L) {
    stop("Internal error: missing R entries.", call. = FALSE)
  }

  if (is.list(Rscale) && !is.data.frame(Rscale)) {
    return(.brma_mv_normalize_Rscale_list(Rscale, R_names))
  }

  Rscale_names <- names(Rscale)
  if (length(Rscale) == 1L && (is.null(Rscale_names) || !nzchar(Rscale_names[[1L]]))) {
    out <- rep(.brma_mv_normalize_Rscale_value(Rscale), n_R)
    names(out) <- R_names
    return(out)
  }

  if (is.null(Rscale_names) || any(!nzchar(Rscale_names))) {
    stop(
      "Multi-entry 'Rscale' must be named to match supplied 'R' entries.",
      call. = FALSE
    )
  }

  .brma_mv_check_Rscale_names(Rscale_names, R_names)
  out <- vapply(
    Rscale[R_names],
    .brma_mv_normalize_Rscale_value,
    character(1)
  )
  names(out) <- R_names
  return(out)
}


.brma_mv_normalize_Rscale_list <- function(Rscale, R_names) {

  n_R <- length(R_names)
  if (length(Rscale) == 0L) {
    stop("'Rscale' must not be an empty list.", call. = FALSE)
  }

  Rscale_names <- names(Rscale)
  if (is.null(Rscale_names)) {
    Rscale_names <- rep("", length(Rscale))
  }
  Rscale_names[is.na(Rscale_names)] <- ""

  if (length(Rscale) == 1L && !nzchar(Rscale_names[[1L]])) {
    if (n_R != 1L) {
      stop(
        "Unnamed list 'Rscale' is only supported with a single 'R' matrix.",
        call. = FALSE
      )
    }
    out <- .brma_mv_normalize_Rscale_value(Rscale[[1L]])
    names(out) <- R_names
    return(out)
  }

  if (any(!nzchar(Rscale_names))) {
    stop("'Rscale' list entries must be named.", call. = FALSE)
  }

  .brma_mv_check_Rscale_names(Rscale_names, R_names)
  out <- vapply(
    Rscale[R_names],
    .brma_mv_normalize_Rscale_value,
    character(1)
  )
  names(out) <- R_names
  return(out)
}


.brma_mv_check_Rscale_names <- function(Rscale_names, R_names) {

  if (anyDuplicated(Rscale_names)) {
    stop("'Rscale' names must be unique.", call. = FALSE)
  }
  if (any(!nzchar(R_names))) {
    stop(
      "Named 'Rscale' cannot be matched when 'R' is unnamed.",
      call. = FALSE
    )
  }

  missing_names <- setdiff(R_names, Rscale_names)
  unknown_names <- setdiff(Rscale_names, R_names)
  if (length(missing_names) > 0L || length(unknown_names) > 0L) {
    details <- character()
    if (length(missing_names) > 0L) {
      details <- c(details, paste0(
        "missing: ", paste(missing_names, collapse = ", ")
      ))
    }
    if (length(unknown_names) > 0L) {
      details <- c(details, paste0(
        "unknown: ", paste(unknown_names, collapse = ", ")
      ))
    }
    stop(
      "'Rscale' names must match supplied 'R' names (",
      paste(details, collapse = "; "),
      ").",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.brma_mv_normalize_Rscale_value <- function(value) {

  if (is.list(value) && !is.data.frame(value)) {
    if (length(value) != 1L) {
      stop("'Rscale' entries must be scalar.", call. = FALSE)
    }
    value <- value[[1L]]
  }

  if (is.logical(value)) {
    if (length(value) != 1L || is.na(value)) {
      stop("'Rscale' logical aliases must be TRUE or FALSE.", call. = FALSE)
    }
    return(if (isTRUE(value)) "cor" else "none")
  }

  if (is.numeric(value)) {
    if (length(value) != 1L || is.na(value) || !is.finite(value) ||
        !(value %in% 0:3)) {
      stop(
        "'Rscale' numeric aliases must be one of 0, 1, 2, or 3.",
        call. = FALSE
      )
    }
    return(c("none", "cor", "cor0", "cov0")[[as.integer(value) + 1L]])
  }

  if (is.character(value)) {
    if (length(value) != 1L || is.na(value) ||
        !value %in% c("cor", "none", "cor0", "cov0")) {
      stop(
        "'Rscale' must be one of 'cor', 'none', 'cor0', or 'cov0'.",
        call. = FALSE
      )
    }
    return(value)
  }

  stop(
    "'Rscale' must be a character, logical, or numeric alias.",
    call. = FALSE
  )
}
