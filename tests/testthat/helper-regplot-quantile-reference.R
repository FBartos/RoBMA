# ---------------------------------------------------------------------------- #
# Independent R references for regression-plot mixture quantiles.
# ---------------------------------------------------------------------------- #
#
# Deterministic quantiles of posterior mixtures of normal distributions.
#
# ---------------------------------------------------------------------------- #
.has_native_regplot_mixture <- function() {

  return(is.loaded("RoBMA_plot_normal_mixture_quantiles", PACKAGE = "RoBMA"))
}

.has_native_regplot_selection_mixture <- function() {

  return(is.loaded("RoBMA_plot_selnorm_mixture_quantiles", PACKAGE = "RoBMA"))
}


.regplot_mixture_interval_quantiles_r <- function(mean_samples, sd_samples, probs) {

  lower <- numeric(ncol(mean_samples))
  upper <- numeric(ncol(mean_samples))

  for (i in seq_len(ncol(mean_samples))) {
    cdf_fun <- function(q) {
      .regplot_normal_mixture_cdf(
        q    = q,
        mean = mean_samples[, i],
        sd   = sd_samples[, i]
      )
    }

    lower[i] <- .regplot_mixture_quantile(
      probs[1], mean_samples[, i], sd_samples[, i], cdf_fun,
      full_support = TRUE
    )
    upper[i] <- .regplot_mixture_quantile(
      probs[2], mean_samples[, i], sd_samples[, i], cdf_fun,
      full_support = TRUE
    )
  }

  return(list(lower = lower, upper = upper))
}


# ---------------------------------------------------------------------------- #
# .regplot_selection_mixture_interval_quantiles
# ---------------------------------------------------------------------------- #
#
# Deterministic quantiles of observed-effect mixtures with selection branches.
#
# ---------------------------------------------------------------------------- #
.regplot_selection_mixture_interval_quantiles_r <- function(mean_samples,
                                                            sd_samples, se,
                                                            probs, setup,
                                                            effect_direction) {

  lower <- numeric(ncol(mean_samples))
  upper <- numeric(ncol(mean_samples))
  full_support <- .plot_selection_mixture_has_full_support(
    selection_context = setup[["selection"]],
    selected_rows      = setup[["is_weightfunction"]]
  )

  for (i in seq_len(ncol(mean_samples))) {
    se_i <- if (length(se) == 1L) se else se[[i]]
    cdf_fun <- function(q) {
      .regplot_selection_mixture_cdf(
        q                = q,
        mean             = mean_samples[, i],
        sd               = sd_samples[, i],
        se               = se_i,
        setup            = setup,
        effect_direction = effect_direction
      )
    }

    lower[i] <- .regplot_mixture_quantile(
      probs[1], mean_samples[, i], sd_samples[, i], cdf_fun,
      full_support = full_support
    )
    upper[i] <- .regplot_mixture_quantile(
      probs[2], mean_samples[, i], sd_samples[, i], cdf_fun,
      full_support = full_support
    )
  }

  return(list(lower = lower, upper = upper))
}


# ---------------------------------------------------------------------------- #
# .regplot_mixture_quantile
# ---------------------------------------------------------------------------- #
#
# Invert a posterior-mixture CDF.
#
# ---------------------------------------------------------------------------- #
.regplot_mixture_quantile <- function(p, mean, sd, cdf_fun, full_support) {

  if (all(sd == 0)) {
    return(unname(stats::quantile(mean, probs = p, names = FALSE, type = 1)))
  }

  spread <- sd
  lower  <- min(mean - 10 * spread, na.rm = TRUE)
  upper  <- max(mean + 10 * spread, na.rm = TRUE)

  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop(
      "Regression-plot quantiles could not be computed from a valid bracketed CDF.",
      call. = FALSE
    )
  }

  obj_fun     <- function(q) cdf_fun(q) - p
  lower_value <- obj_fun(lower)
  upper_value <- obj_fun(upper)
  step        <- max(spread, na.rm = TRUE)
  if (!is.finite(lower_value) || !is.finite(upper_value)) {
    return(.regplot_bracketed_quantile(p, lower, upper, cdf_fun))
  }

  for (i in seq_len(25)) {
    if (lower_value < 0 && upper_value >= 0) {
      break
    }
    if (lower_value >= 0) {
      lower       <- lower - step
      lower_value <- obj_fun(lower)
      if (!is.finite(lower_value)) {
        return(.regplot_bracketed_quantile(p, lower, upper, cdf_fun))
      }
    }
    if (upper_value < 0) {
      upper       <- upper + step
      upper_value <- obj_fun(upper)
      if (!is.finite(upper_value)) {
        return(.regplot_bracketed_quantile(p, lower, upper, cdf_fun))
      }
    }
    step <- step * 2
  }

  if (all(spread > 0) && full_support) {
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
    if (is.finite(out)) {
      return(out)
    }
  }

  return(.regplot_bracketed_quantile(p, lower, upper, cdf_fun))
}


# ---------------------------------------------------------------------------- #
# .regplot_normal_mixture_cdf
# ---------------------------------------------------------------------------- #
.regplot_normal_mixture_cdf <- function(q, mean, sd) {

  cdf_values <- rep(NA_real_, length(mean))
  zero_sd    <- sd == 0

  if (any(zero_sd)) {
    cdf_values[zero_sd] <- as.numeric(q >= mean[zero_sd])
  }
  if (any(!zero_sd)) {
    cdf_values[!zero_sd] <- stats::pnorm(
      q,
      mean = mean[!zero_sd],
      sd   = sd[!zero_sd]
    )
  }

  cdf_values <- .plot_validate_cdf(
    cdf_values,
    "Regression-plot normal mixture"
  )
  return(base::mean(cdf_values))
}


# ---------------------------------------------------------------------------- #
# .regplot_selection_mixture_cdf
# ---------------------------------------------------------------------------- #
.regplot_selection_mixture_cdf <- function(q, mean, sd, se, setup,
                                           effect_direction) {

  cdf_values <- rep(NA_real_, length(mean))
  zero_sd    <- sd == 0

  if (any(zero_sd)) {
    cdf_values[zero_sd] <- as.numeric(q >= mean[zero_sd])
  }

  normal_rows <- !setup[["is_weightfunction"]] & !zero_sd
  if (any(normal_rows)) {
    cdf_values[normal_rows] <- stats::pnorm(
      q,
      mean = mean[normal_rows],
      sd   = sd[normal_rows]
    )
  }

  selected_rows <- setup[["is_weightfunction"]] & !zero_sd
  if (any(selected_rows)) {
    setup[["mu"]] <- mean
    rows <- which(selected_rows)
    cdf_values[rows] <- .funnel_selected_cdf(
      q                = q,
      rows             = rows,
      se               = se,
      total_sd         = sd,
      setup            = setup,
      effect_direction = effect_direction
    )
  }

  cdf_values <- .plot_validate_cdf(
    cdf_values,
    "Regression-plot selected-normal mixture"
  )
  return(base::mean(cdf_values))
}


# ---------------------------------------------------------------------------- #
# .regplot_selection_setup
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# .regplot_bracketed_quantile
# ---------------------------------------------------------------------------- #
#
# Return inf{x: F(x) >= p} within a valid finite CDF bracket.
#
# ---------------------------------------------------------------------------- #
.regplot_bracketed_quantile <- function(p, lower, upper, cdf_fun) {

  lower_cdf <- cdf_fun(lower)
  upper_cdf <- cdf_fun(upper)
  if (!is.finite(lower_cdf) || !is.finite(upper_cdf) ||
      lower_cdf >= p || upper_cdf < p) {
    stop(
      "Regression-plot quantiles could not be computed from a valid bracketed CDF.",
      call. = FALSE
    )
  }
  repeat {
    mid <- 0.5 * lower + 0.5 * upper
    if (mid <= lower || mid >= upper) {
      break
    }
    mid_cdf <- cdf_fun(mid)
    if (!is.finite(mid_cdf)) {
      stop(
        "Regression-plot quantiles could not be computed from a valid bracketed CDF.",
        call. = FALSE
      )
    }
    if (mid_cdf >= p) {
      upper <- mid
    } else {
      lower <- mid
    }
  }

  return(upper)
}
