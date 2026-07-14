# ============================================================================ #
# qCMDE Normalization and Refinement
# ============================================================================ #

.iwmde_qcmde_normalizer_plan <- function(normalization_grid, transform) {

  grid_sequence <- .iwmde_qcmde_grid_sequence(
    normalization_grid = normalization_grid,
    transform          = transform
  )
  if (length(grid_sequence) == 0L) {
    grid_sequence <- list(normalization_grid)
  }

  all_z <- sort(unique(unlist(lapply(grid_sequence, `[[`, "z"),
                              use.names = FALSE)))
  all_x <- .iwmde_from_internal(all_z, transform)
  keep  <- is.finite(all_z) & is.finite(all_x)
  all_z <- all_z[keep]
  all_x <- all_x[keep]

  if (length(all_z) < 2L || any(diff(all_z) <= 0)) {
    all_z <- as.numeric(normalization_grid[["z"]])
    all_x <- as.numeric(normalization_grid[["x"]])
  }

  all_grid <- list(
    x            = all_x,
    z            = all_z,
    log_jacobian = .iwmde_log_jacobian(all_z, transform)
  )

  grid_sequence <- lapply(grid_sequence, function(grid) {
    index <- match(grid[["z"]], all_grid[["z"]])
    keep  <- !is.na(index)
    grid[["x"]]            <- grid[["x"]][keep]
    grid[["z"]]            <- grid[["z"]][keep]
    grid[["log_jacobian"]] <- grid[["log_jacobian"]][keep]
    grid[["all_index"]]    <- index[keep]
    grid
  })
  grid_sequence <- grid_sequence[
    vapply(grid_sequence, function(grid) length(grid[["z"]]) >= 2L, logical(1))
  ]
  if (length(grid_sequence) == 0L) {
    grid_sequence <- list(c(all_grid, list(all_index = seq_along(all_z))))
  }

  return(list(
    grid_sequence      = grid_sequence,
    all_grid           = all_grid,
    pilot_grid         = grid_sequence[[1L]],
    final_grid         = grid_sequence[[length(grid_sequence)]],
    pilot_index        = 1L,
    final_index        = length(grid_sequence),
    validation_index   = length(grid_sequence),
    n_refinement_steps = length(grid_sequence) - 1L
  ))
}


.iwmde_qcmde_grid_sequence <- function(normalization_grid, transform,
                                       max_refinement_steps = 3L) {

  z <- sort(unique(as.numeric(normalization_grid[["z"]])))
  z <- z[is.finite(z)]
  if (length(z) < 2L || any(diff(z) <= 0)) {
    return(list())
  }

  width <- diff(range(z))
  if (!is.finite(width) || width <= 0) {
    return(list())
  }

  n_base <- length(z)
  out    <- vector("list", max_refinement_steps + 1L)
  out[[1L]] <- .iwmde_qcmde_grid_from_z(z, transform)

  for (step in seq_len(max_refinement_steps)) {
    extension <- width * .iwmde_qcmde_extension_fraction(step)
    n_points  <- .iwmde_qcmde_refinement_points(
      n_base = n_base,
      step   = step
    )
    z_step <- seq(
      min(z) - extension,
      max(z) + extension,
      length.out = n_points
    )
    z_step <- sort(unique(c(z, z_step)))
    out[[step + 1L]] <- .iwmde_qcmde_grid_from_z(z_step, transform)
  }

  out <- out[!vapply(out, is.null, logical(1))]

  return(out)
}


.iwmde_qcmde_grid_from_z <- function(z, transform) {

  x    <- .iwmde_from_internal(z, transform)
  keep <- is.finite(z) & is.finite(x)
  z    <- z[keep]
  x    <- x[keep]
  if (length(z) < 2L || any(diff(z) <= 0)) {
    return(NULL)
  }

  return(list(
    x            = x,
    z            = z,
    log_jacobian = .iwmde_log_jacobian(z, transform)
  ))
}


.iwmde_qcmde_extension_fraction <- function(step) {

  fractions <- c(.25, .50, 1)
  step      <- min(length(fractions), max(1L, as.integer(step[[1L]])))

  return(fractions[[step]])
}


.iwmde_qcmde_refinement_points <- function(n_base, step) {

  multiplier <- c(1.5, 2, 2.5)
  step       <- min(length(multiplier), max(1L, as.integer(step[[1L]])))

  return(max(2L, as.integer(ceiling(n_base * multiplier[[step]]))))
}


.iwmde_qcmde_select_refinement <- function(log_q_display,
                                           log_normalizer_sequence,
                                           active_mass, denominator) {

  n_sequence <- length(log_normalizer_sequence)
  if (n_sequence <= 1L) {
    return(list(
      pilot_index        = 1L,
      final_index        = 1L,
      validation_index   = 1L,
      n_refinement_steps = 0L
    ))
  }

  if (n_sequence == 2L) {
    return(list(
      pilot_index        = 1L,
      final_index        = 2L,
      validation_index   = 2L,
      n_refinement_steps = 1L
    ))
  }

  for (index in seq.int(2L, n_sequence - 1L)) {
    candidate_y <- .iwmde_qcmde_density_from_normalizer(
      log_q_display  = log_q_display,
      log_normalizer = log_normalizer_sequence[[index]],
      active_mass    = active_mass,
      denominator    = denominator
    )
    validation_y <- .iwmde_qcmde_density_from_normalizer(
      log_q_display  = log_q_display,
      log_normalizer = log_normalizer_sequence[[index + 1L]],
      active_mass    = active_mass,
      denominator    = denominator
    )
    change <- .iwmde_qcmde_ordinate_change(
      pilot_y = candidate_y,
      final_y = validation_y
    )
    max_change <- .iwmde_max_or_na(change[["relative"]])
    if (is.finite(max_change) &&
        max_change <= .iwmde_qcmde_refinement_target()) {
      return(list(
        pilot_index        = index - 1L,
        final_index        = index,
        validation_index   = index + 1L,
        n_refinement_steps = index - 1L
      ))
    }
  }

  return(list(
    pilot_index        = max(1L, n_sequence - 2L),
    final_index        = n_sequence - 1L,
    validation_index   = n_sequence,
    n_refinement_steps = n_sequence - 2L
  ))
}


.iwmde_qcmde_density_from_normalizer <- function(log_q_display,
                                                 log_normalizer,
                                                 active_mass,
                                                 denominator) {

  keep_rows <- is.finite(log_normalizer)
  return(.iwmde_qcmde_pilot_density(
    log_q_display  = log_q_display,
    log_normalizer = log_normalizer,
    keep_rows      = keep_rows,
    active_mass    = active_mass,
    denominator    = denominator
  ))
}


.iwmde_qcmde_refinement_target <- function() {

  return(.025)
}


.iwmde_qcmde_pilot_density <- function(log_q_display, log_normalizer,
                                       keep_rows, active_mass,
                                       denominator) {

  if (!any(keep_rows)) {
    return(rep(NA_real_, nrow(log_q_display)))
  }

  density_terms <- .iwmde_density_aggregate(
    log_terms = sweep(
      log_q_display[, keep_rows, drop = FALSE],
      2L,
      log_normalizer[keep_rows],
      "-"
    ),
    active_mass = active_mass,
    denominator = denominator
  )

  return(density_terms[["y"]])
}


.iwmde_qcmde_normalizer_change <- function(initial_log_normalizer,
                                           final_log_normalizer) {

  finite <- is.finite(initial_log_normalizer) &
    is.finite(final_log_normalizer)
  if (!any(finite)) {
    return(list(max = NA_real_, p95 = NA_real_, median = NA_real_))
  }

  relative_change <- abs(expm1(
    final_log_normalizer[finite] - initial_log_normalizer[finite]
  ))

  return(list(
    max    = max(relative_change),
    p95    = stats::quantile(relative_change, .95, names = FALSE, type = 8),
    median = stats::median(relative_change)
  ))
}


.iwmde_qcmde_ordinate_change <- function(pilot_y, final_y) {

  relative_change <- rep(NA_real_, length(final_y))
  log_change      <- rep(NA_real_, length(final_y))

  positive_pilot <- is.finite(pilot_y) & pilot_y > 0
  finite_final   <- is.finite(final_y)
  relative_rows  <- positive_pilot & finite_final
  if (any(relative_rows)) {
    relative_change[relative_rows] <- abs(
      final_y[relative_rows] / pilot_y[relative_rows] - 1
    )
  }

  zero_to_positive <- is.finite(pilot_y) & pilot_y == 0 &
    finite_final & final_y > 0
  relative_change[zero_to_positive] <- Inf

  positive_final <- finite_final & final_y > 0
  log_rows       <- positive_pilot & positive_final
  if (any(log_rows)) {
    log_change[log_rows] <- abs(log(final_y[log_rows]) - log(pilot_y[log_rows]))
  }

  return(list(relative = relative_change, log = log_change))
}


.iwmde_normalization_density <- function(log_q_norm, log_normalizer,
                                         log_jacobian, normalization_grid,
                                         active_mass, denominator) {

  y <- numeric(length(normalization_grid))
  for (g in seq_along(normalization_grid)) {
    log_terms <- log_q_norm[g, ] + log_jacobian[g] - log_normalizer
    finite    <- is.finite(log_terms)
    if (any(finite)) {
      max_term         <- max(log_terms[finite])
      scaled_terms     <- exp(log_terms[finite] - max_term)
      y[g]             <- active_mass * exp(max_term) *
        sum(scaled_terms) / denominator
    }
  }

  return(y)
}


.iwmde_row_drop_fraction <- function(n_candidate_rows, n_normalized_rows) {

  n_candidate_rows  <- as.numeric(n_candidate_rows)[1]
  n_normalized_rows <- as.numeric(n_normalized_rows)[1]
  if (!is.finite(n_candidate_rows) || n_candidate_rows <= 0 ||
      !is.finite(n_normalized_rows) || n_normalized_rows < 0) {
    return(NA_real_)
  }

  return(min(1, max(0, 1 - n_normalized_rows / n_candidate_rows)))
}


.iwmde_log_trapz_columns <- function(x, log_y) {

  log_y <- as.matrix(log_y)
  out   <- rep(-Inf, ncol(log_y))
  keep  <- colSums(is.finite(log_y)) >= 2L
  if (!any(keep)) {
    return(out)
  }

  log_y_keep <- log_y[, keep, drop = FALSE]
  max_log    <- apply(log_y_keep, 2L, max, na.rm = TRUE)
  y          <- exp(sweep(log_y_keep, 2L, max_log, "-"))
  y[!is.finite(y)] <- 0
  area       <- .iwmde_trapz_columns(x = x, y = y)
  valid      <- is.finite(area) & area > 0
  out[which(keep)[valid]] <- log(area[valid]) + max_log[valid]

  return(out)
}
