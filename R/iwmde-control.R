# ============================================================================ #
# IWMDE Controls and Reliability Policy
# ============================================================================ #

.density_control_normalize <- function(density_method, density_control = NULL,
                                       allow_normal = FALSE,
                                       purpose = c("density", "ordinate")) {

  density_method <- .density_method_normalize(
    density_method = density_method,
    allow_normal   = allow_normal
  )
  purpose <- match.arg(purpose)
  allowed_names <- c(
    "n_points", "max_samples", "initial_samples", "target_relative_mcse",
    "normalization_points", "normalization_prob", "display_grid"
  )
  defaults <- list(
    n_points             = 100L,
    max_samples          = if (!identical(purpose, "density")) {
      Inf
    } else if (identical(density_method, "IWMDE")) {
      1000L
    } else {
      500L
    },
    initial_samples      = 500L,
    target_relative_mcse = .05,
    normalization_points = NULL,
    normalization_prob   = .999,
    display_grid         = "adaptive"
  )

  if (is.null(density_control)) {
    return(defaults)
  }
  if (!is.list(density_control)) {
    stop("'density_control' must be a named list.", call. = FALSE)
  }
  if (length(density_control) == 0L) {
    return(defaults)
  }

  control_names <- names(density_control)
  if (is.null(control_names) || any(!nzchar(control_names))) {
    stop("'density_control' must be a fully named list.", call. = FALSE)
  }
  duplicated_names <- unique(control_names[duplicated(control_names)])
  if (length(duplicated_names) > 0L) {
    stop(
      "'density_control' contains duplicate setting(s): ",
      paste0("'", duplicated_names, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  unknown_names <- setdiff(control_names, allowed_names)
  if (length(unknown_names) > 0L) {
    stop(
      "'density_control' contains unrecognized setting(s): ",
      paste0("'", unknown_names, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (density_method %in% c("KDE", "normal")) {
    stop(
      "'density_control' is only used when 'density_method' is ",
      "'qCMDE' or 'IWMDE'.",
      call. = FALSE
    )
  }

  for (name in control_names) {
    defaults[[name]] <- density_control[[name]]
  }

  BayesTools::check_int(
    defaults[["n_points"]],
    "density_control$n_points",
    lower = 20
  )
  .iwmde_check_max_samples(defaults[["max_samples"]])
  BayesTools::check_int(
    defaults[["initial_samples"]],
    "density_control$initial_samples",
    lower = 20
  )
  BayesTools::check_real(
    defaults[["target_relative_mcse"]],
    "density_control$target_relative_mcse",
    lower        = 0,
    upper        = 1,
    check_length = 1,
    allow_NA     = FALSE
  )
  if (defaults[["target_relative_mcse"]] <= 0) {
    stop(
      "'density_control$target_relative_mcse' must be higher than 0.",
      call. = FALSE
    )
  }
  if (is.finite(defaults[["max_samples"]]) &&
      defaults[["initial_samples"]] > defaults[["max_samples"]]) {
    if (all(c("initial_samples", "max_samples") %in% control_names)) {
      stop(
        "'density_control$initial_samples' cannot exceed ",
        "'density_control$max_samples'.",
        call. = FALSE
      )
    }
    defaults[["initial_samples"]] <- as.integer(defaults[["max_samples"]])
  }
  BayesTools::check_int(
    defaults[["normalization_points"]],
    "density_control$normalization_points",
    lower        = 20,
    allow_NULL   = TRUE,
    check_length = 1
  )
  BayesTools::check_real(
    defaults[["normalization_prob"]],
    "density_control$normalization_prob",
    lower        = 0,
    upper        = 1,
    check_length = 1,
    allow_NA     = FALSE
  )
  if (defaults[["normalization_prob"]] <= 0) {
    stop(
      "'density_control$normalization_prob' must be higher than 0.",
      call. = FALSE
    )
  }
  defaults[["display_grid"]] <- .iwmde_normalize_display_grid(
    defaults[["display_grid"]]
  )

  return(defaults)
}


# Resolve public density controls while preserving the private ordinate marker.
.iwmde_density_control_resolve <- function(
    density_method, density_control = NULL, allow_normal = FALSE,
    purpose = c("density", "ordinate")) {

  purpose <- match.arg(purpose)
  ordinate_marker <- is.list(density_control) &&
    identical(density_control[["display_grid"]], "ordinate")
  if (ordinate_marker) {
    density_control[["display_grid"]] <- "adaptive"
  }
  control <- .density_control_normalize(
    density_method  = density_method,
    density_control = density_control,
    allow_normal    = allow_normal,
    purpose         = purpose
  )
  if (ordinate_marker) {
    control[["display_grid"]] <- "ordinate"
  }

  return(control)
}


.iwmde_check_max_samples <- function(max_samples) {

  if (length(max_samples) != 1L || is.na(max_samples) ||
      (!is.finite(max_samples) && !identical(as.numeric(max_samples), Inf)) ||
      (is.finite(max_samples) &&
       (max_samples < 20 || max_samples != as.integer(max_samples)))) {
    stop(
      "'density_control$max_samples' must be an integer at least 20 or Inf.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.iwmde_bf_max_relative_mcse <- function() {

  return(.25)
}


.iwmde_bf_warning_relative_mcse <- function() {

  return(.05)
}


.iwmde_bf_min_ess <- function() {

  return(20)
}


.iwmde_bf_warning_min_ess <- function() {

  return(100)
}


.iwmde_bf_max_weight_share <- function() {

  return(.50)
}


.iwmde_bf_warning_weight_share <- function() {

  return(.20)
}


.iwmde_bf_warning_min_finite_terms <- function() {

  return(100)
}


.iwmde_bf_min_finite_terms <- function() {

  return(20)
}


.iwmde_bf_mass_warning_tolerance <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(.025)
  }
  if (identical(estimator, "iwmde")) {
    return(.05)
  }

  return(Inf)
}


.iwmde_bf_mass_fail_tolerance <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(.05)
  }
  if (identical(estimator, "iwmde")) {
    return(.10)
  }

  return(0)
}


.iwmde_density_min_estimator_rows <- function() {

  return(300)
}


.iwmde_density_warning_min_estimator_rows <- function() {

  return(500)
}


.iwmde_density_max_relative_mcse <- function() {

  return(.25)
}


.iwmde_density_warning_relative_mcse <- function() {

  return(.10)
}


.iwmde_density_min_ess <- function(estimator_rows = NA_real_) {

  estimator_rows <- as.numeric(estimator_rows)[1L]
  if (is.finite(estimator_rows) && estimator_rows > 0) {
    return(min(50, max(4, .25 * estimator_rows)))
  }

  return(50)
}


.iwmde_density_warning_min_ess <- function(estimator_rows = NA_real_) {

  estimator_rows <- as.numeric(estimator_rows)[1L]
  if (is.finite(estimator_rows) && estimator_rows > 0) {
    return(min(100, max(20, .50 * estimator_rows)))
  }

  return(100)
}


.iwmde_density_max_weight_share <- function() {

  return(.25)
}


.iwmde_density_warning_weight_share <- function() {

  return(.10)
}


.iwmde_quadrature_warning_tolerance <- function() {

  return(.025)
}


.iwmde_quadrature_fail_tolerance <- function() {

  return(.05)
}


# Whether likelihood-aware density estimation supports this model and estimator.
.iwmde_density_method_supported <- function(object, density_method) {

  density_method <- .density_method_normalize(density_method)
  return(!(
    inherits(object, "brma.glmm") && identical(density_method, "IWMDE")
  ))
}


# Fail before an uncertified estimator can enter any public result.
.iwmde_check_density_method_supported <- function(object, density_method) {

  if (!.iwmde_density_method_supported(object, density_method)) {
    stop(
      "IWMDE density estimation is unavailable for binomial and ",
      "Poisson GLMMs because their high-dimensional conditional weights do ",
      "not meet the bridge-sampling certification tolerance. Use ",
      "density_method = 'qCMDE'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.iwmde_point_ordinate_supported <- function(object, density_method) {

  return(.iwmde_density_method_supported(object, density_method))
}


.iwmde_check_point_ordinate_supported <- function(object, density_method) {

  .iwmde_check_density_method_supported(object, density_method)
}


.iwmde_check_context_density_method_supported <- function(context,
                                                           density_method) {

  density_method <- .density_method_normalize(density_method)
  if (is.null(context[["data"]])) {
    return(invisible(TRUE))
  }
  outcome_type <- .data_outcome_type(context[["data"]])
  if (outcome_type %in% c("bin", "pois") &&
      identical(density_method, "IWMDE")) {
    stop(
      "IWMDE density estimation is unavailable for binomial and ",
      "Poisson GLMMs because their high-dimensional conditional weights do ",
      "not meet the bridge-sampling certification tolerance. Use ",
      "density_method = 'qCMDE'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
