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
# When no row can report a kernel every row is an estimate row, the range is the
# draws-based one, and the display values, the nodes and the midpoints of the
# initial lattice are one call of the joint density (which evaluates each value
# independently of the others in its call, so the values are those of separate
# calls). Contexts with call-owned interpolation grids keep the display values
# alone in the first call, because the grids place their anchors by it, and
# declare the nodes, the midpoints and a reserve of extension steps on each
# side as their queries, so that the extension the tails need is requested from
# the grids too.
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

    return(record_quadrature(.iwmde_qcmde_evaluate_grid(
      evaluation_context, parameter, values, row_states, replacement,
      estimator_rows
    )))
  }
  # Interpolation grids of fixed selection normalizers and covariance sweeps
  # answer exactly the values they are built for; values outside that set are
  # evaluated directly. The set holds the display values, the nodes and
  # midpoints of the lattice, and a reserve of extension steps on each side
  # (.iwmde_qcmde_extension_reserve()) for the extension the tails need.
  direct <- context
  direct[["normalizer_grid"]] <- NULL
  direct[["covariance_grid"]] <- NULL
  reserve <- .iwmde_qcmde_extension_reserve(n_points)
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
      seq(-reserve, lattice[["n_initial"]] - 1L + reserve, by = .5)
    )[["x"]]
  }
  uses_grids <- function(evaluation_context) {

    !is.null(evaluation_context[["normalizer_grid"]]) ||
      !is.null(evaluation_context[["covariance_grid"]])
  }

  # 1. Display values, and the rows' Gaussian kernels. When no row can report a
  # kernel, every row is an estimate row: the range is the draws-based one, the
  # nodes of its lattice and the midpoints between them are known before the
  # display values are evaluated, and all go in one call. Call-owned grids place
  # their anchors by the values of their first request, so a context that uses
  # them keeps the display values alone in it.
  base_lattice  <- .iwmde_qcmde_lattice(base_z, n_points, transform)
  grid_context  <- with_grids(c(display_grid, planned_values(base_lattice)))
  log_q_display <- NULL
  initial_log_q <- NULL
  initial_midpoints <- NULL
  if (!uses_grids(grid_context) &&
      .iwmde_qcmde_estimate_only(context, row_states, ordinary_rows)) {
    base_nodes     <- .iwmde_qcmde_lattice_points(base_lattice,
                                                  seq.int(0L, n_points - 1L))
    base_midpoints <- .iwmde_qcmde_lattice_points(base_lattice,
                                                  seq.int(0L, n_points - 2L) + .5)
    if (all(base_nodes[["valid"]]) && !any(diff(base_nodes[["z"]]) <= 0) &&
        all(base_midpoints[["valid"]])) {
      joint <- evaluate(
        grid_context,
        c(display_grid, base_nodes[["x"]], base_midpoints[["x"]])
      )
      n_display     <- length(display_grid)
      log_q_display <- joint[seq_len(n_display), , drop = FALSE]
      initial_log_q <- joint[n_display + seq_len(n_points), , drop = FALSE]
      initial_midpoints <- list(
        x     = base_midpoints[["x"]],
        log_q = joint[n_display + n_points + seq_len(n_points - 1L), ,
                      drop = FALSE]
      )
      kernel <- attr(joint, "gaussian_kernel", exact = TRUE)
    }
  }
  if (is.null(log_q_display)) {
    log_q_display <- evaluate(grid_context, display_grid)
    kernel        <- attr(log_q_display, "gaussian_kernel", exact = TRUE)
  }

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
  if (!identical(z_range, base_z)) {
    # The nodes and midpoints evaluated ahead belong to the draws-based range.
    initial_log_q     <- NULL
    initial_midpoints <- NULL
    if (uses_grids(grid_context)) {
      grid_context  <- with_grids(c(display_grid, planned_values(lattice)))
      log_q_display <- evaluate(grid_context, display_grid)
    }
  }
  reserved <- if (uses_grids(grid_context)) planned_values(lattice) else numeric()

  # 3. Nodes, extended where estimated rows need it. The extension values in
  # the reserve are requests of the call-owned grids; the rest (and every value
  # without grids) is evaluated directly, and values the joint density cannot
  # evaluate are left missing, which ends the extension there.
  evaluate_directly <- function(values) {

    log_q <- .iwmde_qcmde_try_log_q(direct, parameter, values, row_states,
                                    replacement)
    if (!is.null(log_q)) {
      record_quadrature(log_q)
    }
    log_q
  }
  evaluate_extension <- function(values) {

    inside <- values %in% reserved
    if (!any(inside)) {
      return(evaluate_directly(values))
    }
    out <- matrix(NA_real_, length(values), n_states)
    requested <- tryCatch(evaluate(grid_context, values[inside]),
                          error = function(e) NULL)
    if (is.null(requested)) {
      inside <- rep(FALSE, length(values))
    } else {
      out[inside, ] <- requested
    }
    if (any(!inside)) {
      other <- evaluate_directly(values[!inside])
      if (!is.null(other)) {
        out[!inside, ] <- other
      }
    }

    out
  }
  extension <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = which(laws[["kind"]] %in% "estimate"),
    target           = (1 - normalization_prob) / 2,
    evaluate_inside  = function(values) evaluate(grid_context, values),
    evaluate_outside = evaluate_extension,
    initial_log_q    = initial_log_q
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
  midpoints[["log_q"]] <- matrix(NA_real_, length(midpoint_index), n_states)
  ahead <- if (is.null(initial_midpoints)) {
    rep(NA_integer_, length(midpoint_index))
  } else {
    match(midpoints[["x"]], initial_midpoints[["x"]])
  }
  if (any(!is.na(ahead))) {
    midpoints[["log_q"]][!is.na(ahead), ] <-
      initial_midpoints[["log_q"]][ahead[!is.na(ahead)], , drop = FALSE]
  }
  # The midpoints of the initial lattice are requests of the call-owned grids.
  # The midpoints beyond it are too while their nodes were in the reserve, and
  # a request of the grid that fails there is evaluated without it, like the
  # extension nodes; every other midpoint is evaluated directly.
  inner     <- midpoint_index > 0 & midpoint_index < n_points - 1L
  in_grid   <- is.na(ahead) & midpoints[["x"]] %in% reserved
  planned   <- in_grid & inner
  band      <- in_grid & !inner
  remaining <- is.na(ahead) & !in_grid
  if (any(planned)) {
    midpoints[["log_q"]][planned, ] <- evaluate(grid_context,
                                                midpoints[["x"]][planned])
  }
  if (any(band)) {
    values <- midpoints[["x"]][band]
    midpoints[["log_q"]][band, ] <- tryCatch(
      evaluate(grid_context, values),
      error = function(e) evaluate(direct, values)
    )
  }
  if (any(remaining)) {
    midpoints[["log_q"]][remaining, ] <- evaluate(direct,
                                                  midpoints[["x"]][remaining])
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


# Steps of the lattice, on each side, that the call-owned grids reserve for the
# extension of estimated rows: one slice of the initial range, which is what a
# side with an undecided row first evaluates and where the tails of the usual
# rows end.
.iwmde_qcmde_extension_reserve <- function(n_points) {

  return(as.integer(max(1L, ceiling((n_points - 1L) / 8))))
}


# Whether no row of the batch can report a Gaussian kernel, so that every
# ordinary row is an estimate row. Only the normal-location routes of rows
# without a weight function report one; a fitted binomial or Poisson model, or
# a normal model whose rows all carry a weight function, reports none. An
# unknown model reports nothing certain.
.iwmde_qcmde_estimate_only <- function(context, row_states, rows) {

  data <- context[["data"]]
  if (is.null(data) || length(rows) == 0L) {
    return(FALSE)
  }
  if (.data_outcome_type(data) %in% c("bin", "pois")) {
    return(TRUE)
  }
  if (!identical(.data_outcome_type(data), "norm")) {
    return(FALSE)
  }

  return(all(vapply(row_states[rows], function(state) {
    isTRUE(state[["active_setup"]][["is_weightfunction"]])
  }, logical(1L))))
}


# ---------------------------------------------------------------------------- #
# Session cache of the statement-independent normalization
# ---------------------------------------------------------------------------- #
#
# The normalization of a fitted model's rows (the pass above and what
# .iwmde_density_grid() derives from its grids without the display values)
# depends on the fit, the target and the control, not on the value a statement
# tests. Repeated calls on one fit and target, such as the statements of a
# hypothesis test made one at a time, reuse it; only the display values are
# evaluated again. The keys are content hashes (the fit's draws, data and priors
# among them), so entries of different fits never meet. The cache is an
# environment of the package: it never travels with a saved fit, leaves with the
# session, and is dropped when the package code is reloaded, so a development
# session never reads a normalization the previous code computed. It is bounded
# by the entries it holds and their size.


.iwmde_qcmde_normalization_cache <- new.env(parent = emptyenv())


# Entries kept, oldest first out.
.iwmde_qcmde_cache_limit <- function() 8L


# The largest entry, in bytes, that is kept; a larger one is recomputed.
.iwmde_qcmde_cache_max_bytes <- function() 64 * 1024^2


# The cache of a fitted model's normalizations: the environment of the package.
.iwmde_qcmde_cache_env <- function(context) {

  .iwmde_qcmde_normalization_cache
}


# Drop every cached normalization. A test that replaces a function the
# normalization is computed with, between two requests on one fit, drops them
# too.
.iwmde_qcmde_cache_clear <- function() {

  rm(list = ls(.iwmde_qcmde_normalization_cache, all.names = TRUE),
     envir = .iwmde_qcmde_normalization_cache)

  invisible(NULL)
}


# The cache, the key of the normalization a plan executes, and whether the key
# holds the display values, for `.iwmde_density_grid()`; NULL when the
# normalization is not cached. The key names everything the normalization
# depends on: the fit's content, the target and its conditioning, the rows, the
# replacement, the range and its chart, the control, and the memory bound of
# the call-owned grids. Contexts with
# call-owned interpolation grids also place the display values in the grids'
# queries, so the key holds them and the cached entry keeps their density.
.iwmde_plan_normalization_cache <- function(context, plan, display_grid) {

  cache <- .iwmde_qcmde_cache_env(context)
  if (is.null(cache) || !identical(plan[["method"]], "q_grid_cmde")) {
    return(NULL)
  }
  # The call-owned grids exist for known-V joint selection models only.
  display_in_key <- .is_data_joint_selection(context[["data"]]) &&
    .is_data_known_v(context[["data"]])
  control <- plan[["control"]]
  key <- .iwmde_hash("iwmde_qcmde_normalization", .iwmde_compact_nulls(list(
    schema_version     = .iwmde_schema_version(),
    algorithm_version  = .iwmde_algorithm_version(),
    source_fingerprint = plan[["source_fingerprint"]],
    parameter          = plan[["target"]][["parameter"]],
    target_key         = plan[["target"]][["target_key"]],
    execution_spec     = plan[["execution_spec"]],
    control            = .iwmde_density_control_provenance(control),
    row_budget         = plan[["row_budget"]],
    rows               = .iwmde_plan_rows_provenance(plan),
    estimator_rows     = plan[["rows"]][["estimator_rows"]],
    replacement        = .iwmde_hash("iwmde_replacement", plan[["replacement"]]),
    support            = plan[["support"]][["support"]],
    transform          = plan[["support"]][["transform"]],
    normalization_grid = .iwmde_hash(
      "iwmde_normalization_grid",
      plan[["grids"]][["normalization_grid"]]
    ),
    covariance_bytes   = .known_v_covariance_max_bytes(),
    display            = if (display_in_key) {
      .iwmde_hash("iwmde_display_grid", display_grid)
    }
  )))

  list(cache = cache, key = key, display_in_key = display_in_key)
}


.iwmde_qcmde_cache_get <- function(cache) {

  if (is.null(cache) || !exists(cache[["key"]], envir = cache[["cache"]],
                                inherits = FALSE)) {
    return(NULL)
  }

  get(cache[["key"]], envir = cache[["cache"]], inherits = FALSE)
}


.iwmde_qcmde_cache_set <- function(cache, value) {

  if (is.null(cache) ||
      as.numeric(utils::object.size(value)) > .iwmde_qcmde_cache_max_bytes()) {
    return(invisible(FALSE))
  }
  env   <- cache[["cache"]]
  order <- c(setdiff(env[[".order"]], cache[["key"]]), cache[["key"]])
  assign(cache[["key"]], value, envir = env)
  if (length(order) > .iwmde_qcmde_cache_limit()) {
    evicted <- order[seq_len(length(order) - .iwmde_qcmde_cache_limit())]
    rm(list = evicted, envir = env)
    order <- setdiff(order, evicted)
  }
  assign(".order", order, envir = env)

  invisible(TRUE)
}


# Evaluate values of the joint density of the rows, as a normalization pass
# does: a failure is a construction failure of the target, and the values are
# validated.
.iwmde_qcmde_evaluate_grid <- function(context, parameter, values, row_states,
                                       replacement, estimator_rows) {

  log_q <- tryCatch(
    .iwmde_log_q_grid(
      context     = context,
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

  return(log_q)
}


# The largest of the quadrature changes reported, NA when none was.
.iwmde_qcmde_max_change <- function(changes) {

  changes <- changes[!is.na(changes)]

  if (length(changes) > 0L) max(changes) else NA_real_
}


# What of a normalization pass does not depend on the display values: the
# normalizers of the three grids, the normalization integrals of the pilot and
# the selected grid, the row truncation, and the descriptions of the grids.
# The checks that a normalizer is finite are the pass's own, so an entry is only
# built from a pass that passed them.
.iwmde_qcmde_normalization_summary <- function(pass, parameter, estimator_rows,
                                               active_mass, n_candidate_rows) {

  pilot_grid      <- pass[["initial"]]
  final_grid      <- pass[["nested"]]
  validation_grid <- pass[["nodes"]]
  final_log_normalizer      <- final_grid[["log_normalizer"]]
  validation_log_normalizer <- validation_grid[["log_normalizer"]]
  final_finite      <- is.finite(final_log_normalizer)
  validation_finite <- is.finite(validation_log_normalizer)
  if (any(!final_finite)) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows[!final_finite],
      stage     = "conditional-density normalization",
      detail    = "no finite positive normalizer was obtained on the normalization grid"
    )
  }
  if (any(!validation_finite)) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows[!validation_finite],
      stage     = "conditional-density normalization validation",
      detail    = "the validation grid did not produce a finite positive normalizer"
    )
  }
  norm_y_initial <- .iwmde_normalization_density(
    log_q_norm         = pilot_grid[["log_q"]],
    log_normalizer     = final_log_normalizer,
    log_jacobian       = pilot_grid[["log_jacobian"]],
    normalization_grid = pilot_grid[["z"]],
    active_mass        = active_mass,
    denominator        = n_candidate_rows
  )
  norm_y_final <- .iwmde_normalization_density(
    log_q_norm         = final_grid[["log_q"]],
    log_normalizer     = final_log_normalizer,
    log_jacobian       = final_grid[["log_jacobian"]],
    normalization_grid = final_grid[["z"]],
    active_mass        = active_mass,
    denominator        = n_candidate_rows
  )
  truncation <- .iwmde_qcmde_truncation_summary(.iwmde_qcmde_row_truncation(
    laws      = pass[["laws"]],
    log_mass  = pass[["log_mass"]],
    x_range   = range(final_grid[["x"]]),
    estimates = .iwmde_qcmde_grid_tail_estimates(final_grid)
  ))

  list(
    pilot_log_normalizer         = pilot_grid[["log_normalizer"]],
    log_normalizer               = final_log_normalizer,
    validation_log_normalizer    = validation_log_normalizer,
    pilot_normalization_integral = .iwmde_trapz(pilot_grid[["z"]], norm_y_initial),
    final_normalization_integral = .iwmde_trapz(final_grid[["z"]], norm_y_final),
    truncation                   = truncation,
    normalization_points         = length(final_grid[["x"]]),
    normalization_range          = range(final_grid[["x"]]),
    normalization_initial_points = length(pilot_grid[["x"]]),
    normalization_initial_range  = range(pilot_grid[["x"]]),
    conditional_normalization    = pass[["conditional_normalization"]],
    normalizer_interpolation     = pass[["normalizer_interpolation"]],
    covariance_interpolation     = pass[["covariance_interpolation"]],
    quadrature_change            = pass[["quadrature_change"]],
    n_refinement_steps           = pass[["extension_passes"]]
  )
}


# The normalization of a set of rows for a display grid: the display values and
# the summary. A cached summary is reused, evaluating only the display values
# (or, where they are part of the key, taking their density from the entry);
# otherwise the pass runs and, when it is cacheable, its summary is kept.
.iwmde_qcmde_normalization <- function(context, parameter, display_grid,
                                       normalization_grid, transform,
                                       normalization_prob, row_states,
                                       replacement, estimator_rows,
                                       active_mass, n_candidate_rows,
                                       cache = NULL) {

  entry <- .iwmde_qcmde_cache_get(cache)
  if (!is.null(entry)) {
    log_q_display     <- entry[["log_q_display"]]
    quadrature_change <- entry[["summary"]][["quadrature_change"]]
    if (is.null(log_q_display)) {
      evaluated <- .iwmde_qcmde_evaluate_grid(
        context, parameter, display_grid, row_states, replacement,
        estimator_rows
      )
      quadrature_change <- .iwmde_qcmde_max_change(c(
        quadrature_change,
        attr(evaluated, "max_quadrature_relative_change", exact = TRUE)
      ))
      log_q_display <- matrix(as.numeric(evaluated), nrow(evaluated),
                              ncol(evaluated))
    }

    return(list(log_q_display = log_q_display, summary = entry[["summary"]],
                quadrature_change = quadrature_change))
  }

  pass <- .iwmde_qcmde_normalization_pass(
    context            = context,
    parameter          = parameter,
    display_grid       = display_grid,
    normalization_grid = normalization_grid,
    transform          = transform,
    normalization_prob = normalization_prob,
    row_states         = row_states,
    replacement        = replacement,
    estimator_rows     = estimator_rows
  )
  summary <- .iwmde_qcmde_normalization_summary(
    pass, parameter, estimator_rows, active_mass, n_candidate_rows
  )
  # The display values as a plain matrix: the evaluations that carry them
  # differ in the attributes they leave on it, which nothing reads.
  log_q_display <- pass[["log_q_display"]]
  log_q_display <- matrix(as.numeric(log_q_display), nrow(log_q_display),
                          ncol(log_q_display))
  .iwmde_qcmde_cache_set(cache, list(
    summary       = summary,
    log_q_display = if (isTRUE(cache[["display_in_key"]])) log_q_display
  ))

  list(log_q_display = log_q_display, summary = summary,
       quadrature_change = summary[["quadrature_change"]])
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
