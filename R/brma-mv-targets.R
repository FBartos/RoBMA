# ============================================================================ #
# brma-mv-targets.R
# ============================================================================ #
#
# Central registry of brma.mv() diagnostic target semantics.
#
# ============================================================================ #


.brma_mv_diagnostic_target_table <- function() {

  data.frame(
    method = c(
      "logLik(unit = 'estimate')",
      "add_loo()/loo()",
      "add_waic()/waic()",
      "rstudent()/LOO-PIT residuals",
      "dfbetas()",
      "covratio()",
      "add_marglik()/bridge_sampler()",
      "hatvalues()",
      "rstandard(conditioning_depth = 'marginal')",
      "vif()",
      "dffits()",
      "cooks.distance()",
      "cluster/block LOO"
    ),
    status = c(
      rep("implemented", 12),
      "deferred"
    ),
    target = c(
      rep("estimate-unit log-score", 4),
      rep("estimate-unit PSIS influence", 2),
      "full joint fitted likelihood",
      rep("marginal GLS covariance", 3),
      rep("fixed-location fitted-value influence", 2),
      "joint dependency-block deletion"
    ),
    known_v_semantics = c(
      rep("correlated V uses p(y_i | y_-i, theta) within dependency blocks", 6),
      "uses the full joint fitted known-V likelihood",
      rep("uses V + ZGZ' marginal covariance", 3),
      rep("conditional estimate deletion; reported target fixed_location_fitted_value", 2),
      "deferred for brma.mv()"
    ),
    stringsAsFactors = FALSE
  )
}


.brma_mv_target_row <- function(method) {

  table <- .brma_mv_diagnostic_target_table()
  row   <- table[table[["method"]] == method, , drop = FALSE]

  if (nrow(row) != 1L) {
    stop("Unknown brma.mv() diagnostic target method: ", method,
         call. = FALSE)
  }

  return(row)
}


.brma_mv_target_metadata <- function(method) {

  row <- .brma_mv_target_row(method)

  return(list(
    method              = row[["method"]],
    status              = row[["status"]],
    deletion_target     = "estimate",
    known_v_semantics   = row[["known_v_semantics"]],
    reported_target     = if (method %in% c("dffits()", "cooks.distance()")) {
      "fixed_location_fitted_value"
    } else {
      row[["target"]]
    }
  ))
}


.brma_mv_known_v_backend_metadata <- function(object_or_data) {

  data <- if (is.list(object_or_data) && !is.null(object_or_data[["data"]]) &&
              (inherits(object_or_data, "brma") ||
               inherits(object_or_data, "brma.mv"))) {
    object_or_data[["data"]]
  } else {
    object_or_data
  }

  if (!.is_data_known_v(data)) {
    return(list(
      known_v                              = FALSE,
      known_v_parameterization            = NULL,
      known_v_parameterization_requested  = NULL
    ))
  }

  known_V <- .data_known_v_data(data)
  list(
    known_v                              = TRUE,
    known_v_parameterization            = known_V[["parameterization"]],
    known_v_parameterization_requested  = known_V[["parameterization_requested"]]
  )
}


.brma_mv_known_v_backend_label <- function(metadata) {

  if (!isTRUE(metadata[["known_v"]])) {
    return(NULL)
  }

  label <- paste0("Known-V backend: ", metadata[["known_v_parameterization"]])
  if (!identical(
    metadata[["known_v_parameterization_requested"]],
    metadata[["known_v_parameterization"]]
  )) {
    label <- paste0(
      label,
      " (requested: ",
      metadata[["known_v_parameterization_requested"]],
      ")"
    )
  }

  label
}


.brma_mv_marglik_target_metadata <- function(object) {

  c(
    .brma_mv_target_metadata("add_marglik()/bridge_sampler()"),
    .brma_mv_known_v_backend_metadata(object)
  )
}


.brma_mv_attach_target_metadata <- function(x, method) {

  attr(x, "brma_mv_target") <- .brma_mv_target_metadata(method)
  return(x)
}


.brma_mv_attach_marglik_target_metadata <- function(x, object) {

  if (!inherits(object, "brma.mv")) {
    return(x)
  }

  attr(x, "RoBMA_target") <- .brma_mv_marglik_target_metadata(object)
  return(x)
}
