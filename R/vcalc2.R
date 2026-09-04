#' Sampling covariance with retained construction metadata
#'
#' @description
#' `vcalc2()` is a metadata-preserving wrapper around [metafor::vcalc()]. It
#' returns the same sampling covariance while retaining enough construction
#' information for RoBMA to use a certified diagonal-plus-factor
#' representation when one is available. Calls outside the supported factor
#' contract remain ordinary dense sampling covariances.
#'
#' The optimized representation is currently available when `type`, `obs`, and
#' a numeric two-element `rho` are supplied without `grp1`, `grp2`, `time1`,
#' `time2`, `w1`, `w2`, or `rvars`, and when `nearpd = FALSE`. This covers the
#' common within- and between-construct correlation specification. The
#' covariance returned by `metafor::vcalc()` remains authoritative.
#'
#' @param vi,cluster,subgroup,obs,type,time1,time2,grp1,grp2,w1,w2,data,rho,phi,rvars
#'   Arguments passed unchanged to [metafor::vcalc()].
#' @param checkpd,nearpd,sparse Logical arguments passed unchanged to
#'   [metafor::vcalc()].
#' @param ... Additional arguments passed unchanged to [metafor::vcalc()].
#'
#' @return The object returned by [metafor::vcalc()] with retained RoBMA
#'   construction metadata.
#'
#' @examples \dontrun{
#' if (requireNamespace("metafor", quietly = TRUE)) {
#'   dat <- data.frame(
#'     vi      = c(0.04, 0.05, 0.03),
#'     study   = c(1, 1, 1),
#'     outcome = c("a", "a", "b"),
#'     effect  = 1:3
#'   )
#'   V <- vcalc2(
#'     vi, cluster = study, type = outcome, obs = effect,
#'     rho = c(0.6, 0.3), data = dat
#'   )
#' }
#' }
#'
#' @export
vcalc2 <- function(vi, cluster, subgroup, obs, type, time1, time2,
                   grp1, grp2, w1, w2, data, rho, phi, rvars,
                   checkpd = TRUE, nearpd = FALSE, sparse = FALSE, ...) {

  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop(
      "Package 'metafor' is required for vcalc2(). Install it with ",
      "install.packages(\"metafor\").",
      call. = FALSE
    )
  }

  caller       <- parent.frame()
  matched_call <- match.call(expand.dots = TRUE)
  evaluated    <- .vcalc2_evaluate_supported_call(matched_call, caller)
  if (is.null(evaluated)) {
    metafor_call <- matched_call
    metafor_call[[1L]] <- quote(metafor::vcalc)
    V <- eval(metafor_call, envir = caller)
  } else {
    V <- do.call(metafor::vcalc, evaluated[["metafor_arguments"]])
  }

  metadata <- .vcalc2_metadata(
    V            = V,
    matched_call = matched_call,
    evaluated    = evaluated
  )
  attr(V, "RoBMA_vcalc_metadata") <- metadata
  V
}


.vcalc2_evaluate_supported_call <- function(matched_call, caller) {

  supplied <- names(matched_call)[-1L]
  required <- c("vi", "cluster", "type", "obs", "rho")
  unsupported <- c(
    "grp1", "grp2", "time1", "time2", "w1", "w2", "rvars", "phi"
  )
  if (!all(required %in% supplied) || any(unsupported %in% supplied)) {
    return(NULL)
  }

  data <- if ("data" %in% supplied) {
    eval(matched_call[["data"]], envir = caller)
  } else {
    NULL
  }
  if (is.null(data)) {
    data <- caller
  } else if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  evaluate <- function(name) {
    eval(matched_call[[name]], envir = data, enclos = caller)
  }

  vi       <- evaluate("vi")
  cluster  <- evaluate("cluster")
  type     <- evaluate("type")
  obs      <- evaluate("obs")
  rho      <- evaluate("rho")
  subgroup <- if ("subgroup" %in% supplied) evaluate("subgroup") else NULL

  control_defaults <- list(checkpd = TRUE, nearpd = FALSE, sparse = FALSE)
  controls <- lapply(names(control_defaults), function(name) {
    if (name %in% supplied) {
      eval(matched_call[[name]], envir = caller)
    } else {
      control_defaults[[name]]
    }
  })
  names(controls) <- names(control_defaults)
  formal_names <- setdiff(names(formals(vcalc2)), "...")
  dot_names    <- setdiff(supplied, formal_names)
  dots <- lapply(dot_names, function(name) {
    eval(matched_call[[name]], envir = caller)
  })
  names(dots) <- dot_names

  metafor_arguments <- list(
    vi      = vi,
    cluster = cluster,
    type    = type,
    obs     = obs,
    rho     = rho
  )
  if (!is.null(subgroup)) {
    metafor_arguments[["subgroup"]] <- subgroup
  }
  metafor_arguments <- c(metafor_arguments, controls, dots)

  list(
    vi                = vi,
    cluster           = cluster,
    subgroup          = subgroup,
    type              = type,
    obs               = obs,
    rho               = rho,
    nearpd            = isTRUE(controls[["nearpd"]]) ||
      isTRUE(dots[["nearPD"]]),
    metafor_arguments = metafor_arguments
  )
}


.vcalc2_metadata <- function(V, matched_call, evaluated) {

  metadata <- list(
    version       = 1L,
    call          = matched_call,
    vi            = NULL,
    cluster       = NULL,
    subgroup      = NULL,
    type          = NULL,
    obs           = NULL,
    rho           = NULL,
    factor        = NULL,
    factor_status = "unsupported"
  )
  if (is.null(evaluated)) {
    return(structure(metadata, class = c("RoBMA_vcalc_metadata", "list")))
  }

  metadata[c("vi", "cluster", "subgroup", "type", "obs", "rho")] <-
    evaluated[c("vi", "cluster", "subgroup", "type", "obs", "rho")]
  vi       <- evaluated[["vi"]]
  cluster  <- evaluated[["cluster"]]
  subgroup <- evaluated[["subgroup"]]
  type     <- evaluated[["type"]]
  obs      <- evaluated[["obs"]]
  rho      <- evaluated[["rho"]]

  if (isTRUE(evaluated[["nearpd"]]) || !is.numeric(rho) ||
      !is.null(dim(rho)) || length(rho) != 2L || anyNA(rho) ||
      any(!is.finite(rho))) {
    return(structure(metadata, class = c("RoBMA_vcalc_metadata", "list")))
  }
  if (length(vi) == 1L && length(cluster) > 1L) {
    vi <- rep(vi, length(cluster))
  }
  K <- length(vi)
  if (!is.numeric(vi) || anyNA(vi) || any(!is.finite(vi)) || any(vi < 0) ||
      length(cluster) != K || length(type) != K || length(obs) != K ||
      (!is.null(subgroup) && length(subgroup) != K)) {
    return(structure(metadata, class = c("RoBMA_vcalc_metadata", "list")))
  }

  if (!is.null(subgroup)) {
    cluster <- paste0(cluster, ".", subgroup)
  }

  factor <- .vcalc2_scalar_type_obs_factor(
    vi      = as.numeric(vi),
    cluster = cluster,
    type    = type,
    obs     = obs,
    rho     = as.numeric(rho)
  )
  if (is.null(factor) || !.vcalc2_factor_matches(V, factor)) {
    return(structure(metadata, class = c("RoBMA_vcalc_metadata", "list")))
  }

  metadata[["factor"]]        <- factor
  metadata[["factor_status"]] <- "certified"
  structure(metadata, class = c("RoBMA_vcalc_metadata", "list"))
}


.vcalc2_scalar_type_obs_factor <- function(vi, cluster, type, obs, rho) {

  within  <- rho[[1L]]
  between <- rho[[2L]]
  residual <- 1 - within
  if (residual < 0) {
    return(NULL)
  }

  K                <- length(vi)
  diagonal         <- numeric(K)
  loading_columns  <- list()
  sampling_sei     <- sqrt(vi)
  sampling_blocks  <- split(seq_len(K), cluster)

  add_loading <- function(rows, values) {
    column       <- numeric(K)
    column[rows] <- values
    loading_columns[[length(loading_columns) + 1L]] <<- column
  }

  for (rows in sampling_blocks) {
    types <- unique(as.character(type[rows]))
    type_covariance <- matrix(
      between,
      nrow = length(types),
      ncol = length(types)
    )
    diag(type_covariance) <- within

    if (all(type_covariance == 0)) {
      type_loading <- matrix(0, nrow = length(types), ncol = 0L)
    } else {
      type_loading <- tryCatch(
        t(chol(type_covariance)),
        error = function(e) NULL
      )
      if (is.null(type_loading)) {
        return(NULL)
      }
    }
    for (column in seq_len(ncol(type_loading))) {
      add_loading(
        rows,
        sampling_sei[rows] * type_loading[
          match(as.character(type[rows]), types), column
        ]
      )
    }

    cells <- split(rows, list(type[rows], obs[rows]), drop = TRUE)
    for (cell_rows in cells) {
      if (length(cell_rows) == 1L) {
        diagonal[cell_rows] <- diagonal[cell_rows] + residual * vi[cell_rows]
      } else if (residual > 0) {
        add_loading(
          cell_rows,
          sqrt(residual) * sampling_sei[cell_rows]
        )
      }
    }
  }

  loading <- if (length(loading_columns) == 0L) {
    matrix(0, nrow = K, ncol = 0L)
  } else {
    do.call(cbind, loading_columns)
  }
  known_v_factor(diagonal = diagonal, loading = loading)
}


.vcalc2_factor_matches <- function(V, factor) {

  if (!is.matrix(V) && !inherits(V, "Matrix")) {
    return(FALSE)
  }
  covariance <- as.matrix(V)
  expected <- .known_v_factor_covariance(
    factor[["diagonal"]], factor[["loading"]]
  )
  if (!identical(dim(covariance), dim(expected))) {
    return(FALSE)
  }
  if (anyNA(covariance) || any(!is.finite(covariance)) ||
      anyNA(expected) || any(!is.finite(expected))) {
    return(FALSE)
  }
  scale <- max(1, abs(covariance), abs(expected))
  max(abs(covariance - expected)) <= 100 * .Machine$double.eps * scale
}
