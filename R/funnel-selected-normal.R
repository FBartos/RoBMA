# ============================================================================ #
# funnel-selected-normal.R
# ============================================================================ #

# .funnel_model_averaged_quantile
# ---------------------------------------------------------------------------- #
#
# Invert the model-averaged posterior-row CDF at one SE value.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_quantile <- function(se, p, setup, effect_direction) {

  se_setup <- .funnel_model_averaged_se_setup(
    se               = se,
    setup            = setup,
    effect_direction = effect_direction
  )
  total_sd <- se_setup[["total_sd"]]
  location <- se_setup[["location"]]
  if (all(total_sd == 0)) {
    return(unname(stats::quantile(location, probs = p, names = FALSE, type = 1)))
  }

  spread <- total_sd
  lower  <- min(location - 10 * spread, na.rm = TRUE)
  upper  <- max(location + 10 * spread, na.rm = TRUE)

  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop(
      "Funnel quantiles could not be computed from a valid bracketed CDF.",
      call. = FALSE
    )
  }

  obj_fun <- function(q) {
    .funnel_model_averaged_cdf_precomputed(q, se_setup) - p
  }

  lower_value <- obj_fun(lower)
  upper_value <- obj_fun(upper)
  step        <- max(spread, na.rm = TRUE)
  if (!is.finite(lower_value) || !is.finite(upper_value)) {
    return(.funnel_grid_quantile_precomputed(p, lower, upper, se_setup))
  }
  for (i in seq_len(25)) {
    if (lower_value <= 0 && upper_value >= 0) {
      break
    }
    if (lower_value > 0) {
      lower       <- lower - step
      lower_value <- obj_fun(lower)
      if (!is.finite(lower_value)) {
        return(.funnel_grid_quantile_precomputed(p, lower, upper, se_setup))
      }
    }
    if (upper_value < 0) {
      upper       <- upper + step
      upper_value <- obj_fun(upper)
      if (!is.finite(upper_value)) {
        return(.funnel_grid_quantile_precomputed(p, lower, upper, se_setup))
      }
    }
    step <- step * 2
  }

  if (lower_value > 0 || upper_value < 0) {
    return(.funnel_grid_quantile_precomputed(p, lower, upper, se_setup))
  }

  tolerance <- max(
    .Machine$double.xmin,
    .Machine$double.eps * max(abs(lower), abs(upper), step)
  )
  out <- tryCatch(
    stats::uniroot(
      obj_fun,
      interval = c(lower, upper),
      tol      = tolerance
    )[["root"]],
    error = function(e) NA_real_
  )

  if (is.na(out)) {
    out <- .funnel_grid_quantile_precomputed(p, lower, upper, se_setup)
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_se_setup
# ---------------------------------------------------------------------------- #
#
# Precompute row quantities that stay constant while inverting quantiles for
# one standard error.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_se_setup <- function(se, setup, effect_direction) {

  total_sd      <- .root_sum_squares(se, setup[["tau"]])
  location      <- .funnel_row_location(se, setup, effect_direction)
  zero_sd       <- total_sd == 0
  normal_rows   <- !setup[["is_weightfunction"]] & !zero_sd
  weighted_rows <- setup[["is_weightfunction"]] & !zero_sd
  selection     <- setup[["selection"]]

  out <- list(
    se               = se,
    total_sd         = total_sd,
    location         = location,
    zero_sd          = zero_sd,
    normal_rows      = normal_rows,
    weighted_rows    = weighted_rows,
    setup            = setup,
    effect_direction = effect_direction
  )

  if (any(weighted_rows) && !is.null(selection)) {
    .selection_require_step_evaluable(selection, ".funnel_model_averaged_cdf()")
  }

  if (any(weighted_rows) && !is.null(selection) &&
      se > 0 && isFALSE(selection[["has_phack"]])) {
    rows <- which(weighted_rows)
    out[["step_selection"]] <- .funnel_step_selection_precompute(
      rows      = rows,
      se        = se,
      location  = location,
      total_sd  = total_sd,
      selection = selection
    )
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .funnel_step_selection_precompute
# ---------------------------------------------------------------------------- #
#
# Precompute selected-normal denominators for one standard error.
#
# ---------------------------------------------------------------------------- #
.funnel_step_selection_precompute <- function(rows, se, location, total_sd,
                                              selection) {

  sign          <- selection[["sign"]]
  omega         <- selection[["omega"]][rows, , drop = FALSE]
  mean          <- sign * location[rows]
  sd            <- total_sd[rows]
  n_bins        <- selection[["n_bins"]]
  use_telescope <- isTRUE(selection[["telescope_probabilities"]]) && n_bins > 1L

  if (use_telescope) {
    boundary <- selection[["z_lower"]][seq_len(n_bins - 1L)] * se
    tails    <- vapply(
      boundary,
      stats::pnorm,
      numeric(length(rows)),
      mean       = mean,
      sd         = sd,
      lower.tail = FALSE
    )
    tails      <- matrix(tails, nrow = length(rows), ncol = n_bins - 1L)
    omega_diff <- omega[, seq_len(n_bins - 1L), drop = FALSE] -
      omega[, seq_len(n_bins - 1L) + 1L, drop = FALSE]
    denom <- omega[, n_bins] + rowSums(omega_diff * tails)

    use_telescope <- all(is.finite(denom) & denom > 0)
  }

  if (!use_telescope) {
    denom <- rep(0, length(rows))

    for (b in seq_len(n_bins)) {
      lower <- selection[["z_lower"]][b] * se
      upper <- selection[["z_upper"]][b] * se
      denom <- denom + omega[, b] * .selection_interval_prob_vec(lower, upper, mean, sd)
    }
  }

  return(list(
    rows                    = rows,
    se                      = se,
    sign                    = sign,
    omega                   = omega,
    omega_diff              = if (use_telescope) omega_diff else NULL,
    mean                    = mean,
    sd                      = sd,
    denom                   = denom,
    boundary                = if (use_telescope) boundary else NULL,
    boundary_tail           = if (use_telescope) tails else NULL,
    telescope_probabilities = use_telescope,
    z_lower                 = selection[["z_lower"]],
    z_upper                 = selection[["z_upper"]],
    n_bins                  = n_bins
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_cdf
# ---------------------------------------------------------------------------- #
#
# Average CDF values across posterior rows at one x-coordinate and SE.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_cdf <- function(q, se, setup, effect_direction) {

  se_setup <- .funnel_model_averaged_se_setup(
    se               = se,
    setup            = setup,
    effect_direction = effect_direction
  )

  return(.funnel_model_averaged_cdf_precomputed(q, se_setup))
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_cdf_precomputed
# ---------------------------------------------------------------------------- #
#
# Average CDF values using quantities precomputed for one standard error.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_cdf_precomputed <- function(q, se_setup) {

  S          <- length(se_setup[["location"]])
  total_sd   <- se_setup[["total_sd"]]
  location   <- se_setup[["location"]]
  cdf_values <- rep(NA_real_, S)

  zero_sd <- se_setup[["zero_sd"]]
  if (any(zero_sd)) {
    cdf_values[zero_sd] <- as.numeric(q >= location[zero_sd])
  }

  normal_rows <- se_setup[["normal_rows"]]
  if (any(normal_rows)) {
    cdf_values[normal_rows] <- stats::pnorm(
      q,
      mean = location[normal_rows],
      sd   = total_sd[normal_rows]
    )
  }

  weighted_rows <- se_setup[["weighted_rows"]]
  if (any(weighted_rows)) {
    if (!is.null(se_setup[["step_selection"]])) {
      step_selection <- se_setup[["step_selection"]]
      cdf_values[step_selection[["rows"]]] <-
        .funnel_step_selected_cdf_precomputed(q, step_selection)
    } else {
      rows <- which(weighted_rows)
      cdf_values[rows] <- .funnel_selected_cdf(
        q                = q,
        rows             = rows,
        se               = se_setup[["se"]],
        total_sd         = total_sd,
        setup            = se_setup[["setup"]],
        effect_direction = se_setup[["effect_direction"]]
      )
    }
  }

  cdf_values <- .plot_validate_cdf(cdf_values, "Model-averaged funnel")
  return(mean(cdf_values))
}


# ---------------------------------------------------------------------------- #
# .funnel_step_selected_cdf_precomputed
# ---------------------------------------------------------------------------- #
#
# Selected-normal CDF for step-selection rows with precomputed denominators.
#
# ---------------------------------------------------------------------------- #
.funnel_step_selected_cdf_precomputed <- function(q, step_selection) {

  q_signed <- step_selection[["sign"]] * q

  if (isTRUE(step_selection[["telescope_probabilities"]])) {
    q_tail <- stats::pnorm(
      q_signed,
      mean       = step_selection[["mean"]],
      sd         = step_selection[["sd"]],
      lower.tail = FALSE
    )
    q_cdf <- stats::pnorm(
      q_signed,
      mean = step_selection[["mean"]],
      sd   = step_selection[["sd"]]
    )

    if (step_selection[["sign"]] == 1L) {
      mass <- step_selection[["omega"]][, step_selection[["n_bins"]]] * q_cdf

      for (b in seq_len(step_selection[["n_bins"]] - 1L)) {
        if (q_signed > step_selection[["boundary"]][b]) {
          mass <- mass + step_selection[["omega_diff"]][, b] *
            (step_selection[["boundary_tail"]][, b] - q_tail)
        }
      }
    } else {
      mass <- step_selection[["omega"]][, step_selection[["n_bins"]]] * q_tail

      for (b in seq_len(step_selection[["n_bins"]] - 1L)) {
        tail_part <- if (q_signed >= step_selection[["boundary"]][b]) {
          q_tail
        } else {
          step_selection[["boundary_tail"]][, b]
        }
        mass <- mass + step_selection[["omega_diff"]][, b] *
          tail_part
      }
    }

    out <- mass / step_selection[["denom"]]
    if (all(is.finite(out)) && all(out >= 0 & out <= 1)) {
      return(out)
    }
  }

  mass     <- rep(0, length(step_selection[["rows"]]))

  for (b in seq_len(step_selection[["n_bins"]])) {
    lower <- step_selection[["z_lower"]][b] * step_selection[["se"]]
    upper <- step_selection[["z_upper"]][b] * step_selection[["se"]]

    if (step_selection[["sign"]] == 1L) {
      selected_mass <- .selection_interval_prob_vec(
        lower,
        min(upper, q_signed),
        step_selection[["mean"]],
        step_selection[["sd"]]
      )
    } else {
      selected_mass <- .selection_interval_prob_vec(
        max(lower, q_signed),
        upper,
        step_selection[["mean"]],
        step_selection[["sd"]]
      )
    }

    mass <- mass + step_selection[["omega"]][, b] * selected_mass
  }

  return(.plot_validate_cdf(
    mass / step_selection[["denom"]],
    "Selected-normal funnel",
    operation_count = 4L * step_selection[["n_bins"]] + 4L
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_row_location
# ---------------------------------------------------------------------------- #
#
# Row-specific expected observed location on the original effect-size scale.
#
# ---------------------------------------------------------------------------- #
.funnel_row_location <- function(se, setup, effect_direction) {

  direction <- ifelse(effect_direction == "negative", -1, 1)

  return(
    setup[["mu"]] +
      direction * setup[["PET"]] * se +
      direction * setup[["PEESE"]] * se^2
  )
}


# ---------------------------------------------------------------------------- #
# .funnel_selected_cdf
# ---------------------------------------------------------------------------- #
#
# Selected-normal CDF for selection-model posterior rows.
#
# ---------------------------------------------------------------------------- #
.funnel_selected_cdf <- function(q, rows, se, total_sd, setup, effect_direction) {

  selection <- setup[["selection"]]

  if (is.null(selection)) {
    return(stats::pnorm(q, mean = setup[["mu"]][rows], sd = total_sd[rows]))
  }

  .selection_require_step_evaluable(selection, ".funnel_selected_cdf()")

  if (se <= 0) {
    return(.funnel_selected_cdf_zero_se(
      q        = q,
      rows     = rows,
      total_sd = total_sd,
      setup    = setup
    ))
  }

  local_context <- BayesTools::selection_context_subset_rows(
    context = selection,
    rows    = rows
  )
  cdf           <- .selection_step_cdf_matrix(
    q                 = q,
    mean              = matrix(setup[["mu"]][rows], ncol = 1),
    sd                = matrix(total_sd[rows], ncol = 1),
    sei               = se,
    selection_context = local_context,
    lower.tail        = TRUE
  )

  return(as.numeric(cdf[, 1]))
}


# ---------------------------------------------------------------------------- #
# .funnel_selected_cdf_zero_se
# ---------------------------------------------------------------------------- #
#
# Limit of the selected-normal CDF as the standard error approaches zero.
#
# Selection bins are defined on signed z = sign * y / se. At se -> 0, finite
# z-bin intervals collapse to zero mass. Only the two tail intervals remain:
# the upper z tail maps to signed effects above 0, and the lower z tail maps to
# signed effects below 0. This explicit limit avoids evaluating the native
# selected-normal kernel at the singular se = 0 boundary.
#
# ---------------------------------------------------------------------------- #
.funnel_selected_cdf_zero_se <- function(q, rows, total_sd, setup) {

  selection <- setup[["selection"]]
  sign      <- selection[["sign"]]
  omega     <- selection[["omega"]][rows, , drop = FALSE]
  mean      <- sign * setup[["mu"]][rows]
  sd        <- total_sd[rows]
  q_signed  <- sign * q
  S         <- length(rows)

  denom <- rep(0, S)
  mass  <- rep(0, S)

  for (b in seq_len(selection[["n_bins"]])) {
    z_lower <- selection[["z_lower"]][b]
    z_upper <- selection[["z_upper"]][b]

    if (is.infinite(z_lower) && z_lower < 0 &&
        is.infinite(z_upper) && z_upper > 0) {
      lower <- -Inf
      upper <-  Inf
    } else if (is.infinite(z_upper) && z_upper > 0) {
      lower <- 0
      upper <- Inf
    } else if (is.infinite(z_lower) && z_lower < 0) {
      lower <- -Inf
      upper <- 0
    } else {
      next
    }

    bin_mass <- .selection_interval_prob_vec(lower, upper, mean, sd)
    denom    <- denom + omega[, b] * bin_mass

    if (sign == 1L) {
      selected_mass <- .selection_interval_prob_vec(lower, min(upper, q_signed), mean, sd)
    } else {
      selected_mass <- .selection_interval_prob_vec(max(lower, q_signed), upper, mean, sd)
    }
    mass <- mass + omega[, b] * selected_mass
  }

  out <- mass / denom
  return(.plot_validate_cdf(
    out,
    "Zero-SE selected-normal funnel",
    operation_count = 4L * selection[["n_bins"]] + 4L
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_grid_quantile
# ---------------------------------------------------------------------------- #
#
# Fallback inversion for discontinuous CDFs from degenerate posterior rows.
#
# ---------------------------------------------------------------------------- #
.funnel_grid_quantile <- function(p, lower, upper, se, setup, effect_direction) {

  se_setup <- .funnel_model_averaged_se_setup(
    se               = se,
    setup            = setup,
    effect_direction = effect_direction
  )

  return(.funnel_grid_quantile_precomputed(p, lower, upper, se_setup))
}


# ---------------------------------------------------------------------------- #
# .funnel_grid_quantile_precomputed
# ---------------------------------------------------------------------------- #
#
# Fallback inversion using a precomputed standard-error context.
#
# ---------------------------------------------------------------------------- #
.funnel_grid_quantile_precomputed <- function(p, lower, upper, se_setup) {

  grid <- seq(lower, upper, length.out = 1000)
  cdf  <- vapply(
    grid,
    .funnel_model_averaged_cdf_precomputed,
    numeric(1),
    se_setup = se_setup
  )
  if (any(!is.finite(cdf))) {
    stop(
      "Funnel quantiles could not be computed from a valid bracketed CDF.",
      call. = FALSE
    )
  }
  index <- which(cdf >= p)[1]

  if (is.na(index)) {
    stop(
      "Funnel quantiles could not be computed from a valid bracketed CDF.",
      call. = FALSE
    )
  }

  return(grid[index])
}
