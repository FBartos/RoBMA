# ============================================================================ #
# qCMDE Normalization
# ============================================================================ #
#
# Every qCMDE row is normalized by the trapezoid rule on the chart z. One pass
#
#   1. evaluates the display values, which also reports each row's Gaussian
#      likelihood kernel when the likelihood route has one;
#   2. places the normalization range at the rows' conditional quantiles at
#      `normalization_prob` (R/iwmde-qcmde-range.R);
#   3. evaluates a uniform lattice of `normalization_points` nodes over that
#      range and extends it by whole lattice steps where estimated rows need
#      it; and
#   4. evaluates the midpoints of the final lattice.
#
# The nodes and midpoints form the selected nested grid (2N - 1 points); the
# node grid over the same range validates it. The two differ only in their
# spacing, so their disagreement measures discretization alone, and it
# overstates the error of the finer selected grid. The mass the range leaves
# out of each row is reported separately as truncation.

.iwmde_qcmde_normalization_pass <- function(context, parameter, display_grid,
                                            normalization_grid, transform,
                                            normalization_prob, row_states,
                                            replacement, estimator_rows) {

  n_states         <- length(row_states)
  n_points         <- length(normalization_grid[["z"]])
  base_z           <- range(normalization_grid[["z"]])
  ordinary_rows    <- .iwmde_qcmde_ordinary_rows(row_states)
  conditional_rows <- setdiff(seq_len(n_states), ordinary_rows)
  conditional_normalizers <- lapply(conditional_rows, function(row) {
    tryCatch(.iwmde_retained_location_normalizer(row_states[[row]]), error = function(e) {
      .iwmde_stop_construction_failure("q_grid_cmde", parameter, estimator_rows[[row]],
        stage = "retained-location conditional normalization", detail = conditionMessage(e))
    })
  })
  conditional_log_mass <- vapply(conditional_normalizers, `[[`, numeric(1L), "log_normalizer")
  quadrature_changes <- vapply(conditional_normalizers, `[[`, numeric(1L), "relative_error")

  record_quadrature <- function(log_q) {

    change <- attr(log_q, "max_quadrature_relative_change", exact = TRUE)
    if (is.numeric(change) && length(change) == 1L && !is.na(change)) {
      quadrature_changes <<- c(quadrature_changes, change)
    }

    invisible(log_q)
  }
  evaluate <- function(evaluation_context, values) {

    log_q <- tryCatch(
      .iwmde_log_q_grid(
        context     = evaluation_context,
        parameter   = parameter,
        values      = values,
        row_states  = row_states,
        replacement = replacement
      ),
      error = function(e) {
        if (inherits(e, "iwmde_construction_error")) {
          stop(e)
        }
        .iwmde_stop_construction_failure(
          estimator = "q_grid_cmde",
          parameter = parameter,
          rows      = estimator_rows,
          stage     = "joint-density grid evaluation",
          detail    = conditionMessage(e)
        )
      }
    )
    .iwmde_validate_log_grid(
      log_q_grid = log_q,
      estimator  = "q_grid_cmde",
      parameter  = parameter,
      rows       = estimator_rows,
      n_values   = length(values),
      stage      = "joint-density grid evaluation"
    )

    return(record_quadrature(log_q))
  }
  # Interpolation grids of fixed selection normalizers and covariance sweeps
  # answer exactly the values they are built for; values outside that set are
  # evaluated directly.
  direct <- context
  direct[["normalizer_grid"]] <- NULL
  direct[["covariance_grid"]] <- NULL
  with_grids <- function(values) {

    out <- context
    out[["normalizer_grid"]] <- .selection_normalizer_grid(out, values,
      vapply(row_states, `[[`, integer(1L), "row_index"))
    out[["covariance_grid"]] <- .selection_covariance_grid(out, values,
      vapply(row_states, `[[`, integer(1L), "row_index"),
      row_states, replacement, parameter)
    out
  }
  planned_values <- function(lattice) {

    .iwmde_qcmde_lattice_points(
      lattice,
      seq(0, lattice[["n_initial"]] - 1L, by = .5)
    )[["x"]]
  }
  uses_grids <- function(evaluation_context) {

    !is.null(evaluation_context[["normalizer_grid"]]) ||
      !is.null(evaluation_context[["covariance_grid"]])
  }

  # 1. Display values, and the rows' Gaussian kernels.
  base_lattice  <- .iwmde_qcmde_lattice(base_z, n_points, transform)
  grid_context  <- with_grids(c(display_grid, planned_values(base_lattice)))
  log_q_display <- evaluate(grid_context, display_grid)
  kernel        <- attr(log_q_display, "gaussian_kernel", exact = TRUE)

  # 2. The range: the rows' conditional quantiles.
  laws <- .iwmde_qcmde_row_laws(
    context     = context,
    row_states  = row_states,
    replacement = replacement,
    kernel      = kernel,
    rows        = ordinary_rows
  )
  intervals <- .iwmde_qcmde_law_intervals(laws, normalization_prob)
  laws      <- intervals[["laws"]]
  z_range   <- .iwmde_qcmde_initial_range(
    laws      = laws,
    intervals = intervals[["intervals"]],
    base_z    = base_z,
    transform = transform
  )
  lattice <- .iwmde_qcmde_lattice(z_range, n_points, transform)
  if (!identical(z_range, base_z) && uses_grids(grid_context)) {
    grid_context  <- with_grids(c(display_grid, planned_values(lattice)))
    log_q_display <- evaluate(grid_context, display_grid)
  }

  # 3. Nodes, extended where estimated rows need it.
  extension <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = which(laws[["kind"]] %in% "estimate"),
    target           = (1 - normalization_prob) / 2,
    evaluate_inside  = function(values) evaluate(grid_context, values),
    evaluate_outside = function(values) {
      log_q <- .iwmde_qcmde_try_log_q(direct, parameter, values, row_states,
                                      replacement)
      if (!is.null(log_q)) {
        record_quadrature(log_q)
      }
      log_q
    }
  )
  if (is.null(extension)) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows,
      stage     = "conditional-density normalization",
      detail    = "the normalization range has no valid grid"
    )
  }
  nodes <- extension[["nodes"]]

  # 4. Midpoints of the final lattice.
  midpoint_index <- nodes[["index"]][-length(nodes[["index"]])] + .5
  midpoints      <- .iwmde_qcmde_lattice_points(lattice, midpoint_index)
  if (!all(midpoints[["valid"]])) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows,
      stage     = "conditional-density normalization",
      detail    = "the nested validation grid has an unrepresentable value"
    )
  }
  planned  <- midpoint_index < n_points - 1L & midpoint_index > 0
  midpoints[["log_q"]] <- matrix(NA_real_, length(midpoint_index), n_states)
  if (any(planned)) {
    midpoints[["log_q"]][planned, ] <- evaluate(grid_context,
                                                midpoints[["x"]][planned])
  }
  if (any(!planned)) {
    midpoints[["log_q"]][!planned, ] <- evaluate(direct,
                                                 midpoints[["x"]][!planned])
  }

  initial <- nodes[["index"]] >= 0 & nodes[["index"]] <= n_points - 1L
  grids   <- list(
    initial = .iwmde_qcmde_subset_grid(nodes, initial),
    nodes   = .iwmde_qcmde_subset_grid(nodes, rep(TRUE, length(initial))),
    nested  = .iwmde_qcmde_merge_grids(nodes, midpoints)
  )
  if (any(diff(grids[["nested"]][["z"]]) <= 0)) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows,
      stage     = "conditional-density normalization",
      detail    = "the nested validation grid is not increasing"
    )
  }
  for (name in names(grids)) {
    log_normalizer <- numeric(n_states)
    if (length(ordinary_rows) > 0L) {
      log_normalizer[ordinary_rows] <- .iwmde_log_trapz_columns(
        x     = grids[[name]][["z"]],
        log_y = grids[[name]][["log_q"]][, ordinary_rows, drop = FALSE] +
          grids[[name]][["log_jacobian"]]
      )
    }
    log_normalizer[conditional_rows] <- conditional_log_mass
    grids[[name]][["log_normalizer"]] <- log_normalizer
  }

  return(list(
    log_q_display     = log_q_display,
    initial           = grids[["initial"]],
    nodes             = grids[["nodes"]],
    nested            = grids[["nested"]],
    laws              = laws,
    log_mass          = intervals[["log_mass"]],
    extension_steps   = extension[["steps"]],
    extension_passes  = extension[["passes"]],
    conditional_normalization = list(rows = conditional_rows,
      methods = vapply(conditional_normalizers, `[[`, character(1), "method")),
    quadrature_change = if (length(quadrature_changes) > 0L) {
      max(quadrature_changes)
    } else {
      NA_real_
    },
    normalizer_interpolation = .selection_normalizer_grid_diagnostics(
      grid_context[["normalizer_grid"]]
    ),
    covariance_interpolation = .selection_covariance_grid_diagnostics(
      grid_context[["covariance_grid"]]
    )
  ))
}


.iwmde_qcmde_subset_grid <- function(points, keep) {

  return(list(
    x            = points[["x"]][keep],
    z            = points[["z"]][keep],
    log_jacobian = points[["log_jacobian"]][keep],
    log_q        = points[["log_q"]][keep, , drop = FALSE]
  ))
}


.iwmde_qcmde_merge_grids <- function(nodes, midpoints) {

  index <- order(c(nodes[["index"]], midpoints[["index"]]))

  return(list(
    x            = c(nodes[["x"]], midpoints[["x"]])[index],
    z            = c(nodes[["z"]], midpoints[["z"]])[index],
    log_jacobian = c(nodes[["log_jacobian"]], midpoints[["log_jacobian"]])[index],
    log_q        = rbind(nodes[["log_q"]], midpoints[["log_q"]])[index, , drop = FALSE]
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


# Quadrature tolerance of the retained-location conditional normalizers.
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


.iwmde_log_trapz_columns <- function(x, log_y) {

  log_y <- as.matrix(log_y)
  out   <- rep(-Inf, ncol(log_y))
  keep  <- colSums(is.finite(log_y)) >= 2L
  if (!any(keep)) {
    return(out)
  }

  # Dropping no column needs no copy of the grid, which carries one column per
  # retained posterior row.
  log_y_keep <- if (all(keep)) log_y else log_y[, keep, drop = FALSE]
  max_log    <- apply(log_y_keep, 2L, max, na.rm = TRUE)
  # The column sweep written directly on the matrix: same subtraction, one
  # allocation instead of sweep()'s permuted copy.
  y          <- exp(log_y_keep - rep(max_log, each = nrow(log_y_keep)))
  # Every retained column was divided by its own finite maximum, so the scaled
  # weights lie in [0, 1] and only a missing log density can leave one
  # non-finite. Testing for that costs no mask over the whole grid.
  if (anyNA(y)) {
    y[!is.finite(y)] <- 0
  }
  area       <- .iwmde_trapz_columns(x = x, y = y)
  valid      <- is.finite(area) & area > 0
  out[which(keep)[valid]] <- log(area[valid]) + max_log[valid]

  return(out)
}
