# ============================================================================ #
# funnel-data.R
# ============================================================================ #

# .funnel_data_outcome
# ---------------------------------------------------------------------------- #
#
# Generate funnel plot data for outcome mode (no mods AND no scale).
#
# Shows observed effect sizes against the fitted sampling distribution funnel.
#
# @param x                      brma object
# @param sampling_heterogeneity logical; incorporate tau into funnel
# @param sampling_bias          logical; incorporate bias into funnel
# @param max_samples            maximum posterior samples for biased contours
# @param dots                   list of graphical parameters
#
# @return list with funnel plot data components
#
# ---------------------------------------------------------------------------- #
.funnel_data_outcome <- function(x, sampling_heterogeneity, sampling_bias,
                                 max_samples, dots) {
  # get model characteristics
  is_weightfunction <- .is_weightfunction(x)
  is_PET            <- .is_PET(x)
  is_PEESE          <- .is_PEESE(x)
  effect_direction  <- .effect_direction(x)
  is_glmm            <- .outcome_type(x) %in% c("bin", "pois")

  if (is_glmm) {
    warning(
      "GLMM outcome-mode funnel contours are a descriptive normal ",
      "effect-size approximation; they are not coverage intervals from ",
      "the fitted discrete likelihood.",
      call. = FALSE
    )
  }

  # get observed effect sizes (yi) - these go on the x-axis
  yi <- .outcome_data_yi(x)

  # get pooled effect estimate - this is the funnel center
  mu      <- pooled_effect(x)
  mu_mean <- summary(mu)["mu", "Mean"]

  # get standard errors
  se <- .outcome_data_sei(x)
  K  <- length(se)

  # get tau for incorporating heterogeneity into sampling distribution
  if (sampling_heterogeneity) {
    tau <- .get_funnel_tau(x)
  } else {
    tau <- 0
  }

  # compute standard-error plotting range
  se_range <- pretty(c(0, max(se)))

  # apply custom y-axis limits if provided
  if (!is.null(dots[["ylim"]])) {
    se_range <- pretty(dots[["ylim"]])
  }

  # set axis labels - mark the GLMM effect-size approximation explicitly
  default_xlab <- if (is_glmm) {
    "Observed Effect Size (Descriptive Normal Approximation)"
  } else {
    "Observed Effect Size"
  }
  xlab <- if (!is.null(dots[["xlab"]])) dots[["xlab"]] else default_xlab
  ylab <- if (!is.null(dots[["ylab"]])) dots[["ylab"]] else "Standard Error"

  # set data for reference line
  if (!is.null(dots[["refline"]])) {
    df_refline <- data.frame(
      x = rep(dots[["refline"]], 2),
      y = range(se_range)
    )
  } else {
    # use model center
    df_refline <- NULL
  }

  # determine y-axis plot limit before computing clipped contours
  if (!is.null(dots[["ylim"]])) ylim <- dots[["ylim"]] else ylim <- range(se_range)

  # compute funnel region STRICTLY within ylim to avoid artifacts
  # generate sequence within ylim
  se_sequence_clipped <- seq(min(ylim), max(ylim), length.out = 100)

  # compute contours for clipped sequence
  ci_offsets_clipped <- .get_funnel_quantiles(
    x                      = x,
    se_sequence            = se_sequence_clipped,
    mu                     = mu_mean,
    tau                    = tau,
    sampling_heterogeneity = sampling_heterogeneity,
    sampling_bias          = sampling_bias,
    is_weightfunction      = is_weightfunction,
    is_PET                 = is_PET,
    is_PEESE               = is_PEESE,
    effect_direction       = effect_direction,
    max_samples            = max_samples
  )
  ci_left  <- ci_offsets_clipped$lower
  ci_right <- ci_offsets_clipped$upper

  # create data frames for plotting
  x_range <- pretty(c(ci_left, ci_right, yi))

  # apply custom x-axis limits after funnel computation
  if (!is.null(dots[["xlim"]])) {
    x_range <- pretty(dots[["xlim"]])
  }

  if (!is.null(dots[["xlim"]])) xlim <- dots[["xlim"]] else xlim <- range(x_range)

  # Clamp CIs to xlim to prevent drawing outside plot area (for POLYGON)
  ci_left_clipped  <- pmax(ci_left, min(xlim))
  ci_right_clipped <- pmin(ci_right, max(xlim))

  # Compute reference line for plotting (if not user-specified)
  if (is.null(df_refline)) {
      df_refline <- data.frame(
        x = ci_offsets_clipped$mid,
        y = se_sequence_clipped
      )
      # Clip refline to xlim just in case
      df_refline <- .clip_line_x(df_refline$x, df_refline$y, xlim)
  }

  # create output data structures
  # Points: observed data
  df_points <- data.frame(
    x = yi,
    y = se
  )
  # Funnel: use CLIPPED data for filled polygon to respect limits
  df_funnel <- data.frame(
    x = c(rev(ci_left_clipped), ci_right_clipped),
    y = c(rev(se_sequence_clipped), se_sequence_clipped)
  )
  # Edges: use CLIPPED lines using exact intersection
  edge1 <- .clip_line_x(ci_left, se_sequence_clipped, xlim)
  edge2 <- .clip_line_x(ci_right, se_sequence_clipped, xlim)

  df_funnel_edge1 <- data.frame(
    x = edge1$x,
    y = edge1$y
  )
  df_funnel_edge2 <- data.frame(
    x = edge2$x,
    y = edge2$y
  )
  df_background <- data.frame(
    x = c(min(xlim), max(xlim), max(xlim), min(xlim)),
    y = c(min(ylim), min(ylim), max(ylim), max(ylim))
  )

  out <- list(
    points       = df_points,
    funnel       = df_funnel,
    funnel_edge1 = df_funnel_edge1,
    funnel_edge2 = df_funnel_edge2,
    background   = df_background,
    x_range      = x_range, # Keep original ticks for axis drawing if needed
    y_range      = se_range,
    xlim         = xlim,
    ylim         = ylim,
    xlab         = xlab,
    ylab         = ylab,
    refline      = df_refline
  )
  if (is_glmm) {
    out[["approximation"]] <- list(
      method                              = "normal_effect_size_approximation",
      descriptive                         = TRUE,
      fitted_discrete_likelihood_coverage = FALSE
    )
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .funnel_data_residual
# ---------------------------------------------------------------------------- #
#
# Generate funnel plot data for residual mode (mods OR scale).
#
# @param x                  brma object
# @param type               character; type of residuals
# @param unit               character; output unit
# @param conditioning_depth character; conditioning depth
# @param dots               list of graphical parameters
#
# @return list with funnel plot data components
#
# ---------------------------------------------------------------------------- #
.funnel_data_residual <- function(x, type, unit, conditioning_depth, dots) {

  # get residuals based on type
  # rstandard and rstudent return data.frame with resid, se, z columns
  if (type == "rstandard") {
    res_obj <- rstandard.brma(
      model              = x,
      unit               = unit,
      conditioning_depth = conditioning_depth
    )
    res <- res_obj$resid
    se  <- res_obj$se
  } else if (type == "LOO-PIT" || type == "rstudent") {
    res_obj <- rstudent.brma(x, unit = unit)
    res <- res_obj$resid
    se  <- res_obj$se
  } else if (type == "outcome") {
    # raw outcome residuals
    res_obj <- residuals.brma(
      object             = x,
      type               = "outcome",
      unit               = unit,
      conditioning_depth = conditioning_depth
    )
    if (is.data.frame(res_obj)) {
      res <- res_obj$resid
      se  <- res_obj$se
    } else {
      res <- res_obj
      se  <- .outcome_data_sei(x)
    }
  }
  K <- length(se)

  # Compute a descriptive sampling-error reference region. Residual funnels do
  # not incorporate tau and the bounds are not posterior-predictive intervals.
  # se_range determines the axis ticks
  se_range    <- pretty(c(0, max(se)))

  # determine plot limits
  # For xlim: pretty range of residuals + quantiles
  # For ylim: defaults to se_range range
  x_init_range <- pretty(c(stats::qnorm(0.025) * max(se_range), stats::qnorm(0.975) * max(se_range), res))

  if (!is.null(dots[["xlim"]])) xlim <- dots[["xlim"]] else xlim <- range(x_init_range)
  if (!is.null(dots[["ylim"]])) ylim <- dots[["ylim"]] else ylim <- range(se_range)

  # generate sequence strictly within ylim for clean polygons
  se_sequence_clipped <- seq(min(ylim), max(ylim), length.out = 100)

  # Reference bounds use qnorm(p) * se.
  ci_left  <- stats::qnorm(0.025) * se_sequence_clipped
  ci_right <- stats::qnorm(0.975) * se_sequence_clipped

  # Clamp CIs to xlim to prevent drawing outside plot area (for POLYGON)
  ci_left_clipped  <- pmax(ci_left, min(xlim))
  ci_right_clipped <- pmin(ci_right, max(xlim))

  # set axis labels
  xlab <- if (!is.null(dots[["xlab"]])) dots[["xlab"]] else "Residual"
  ylab <- if (!is.null(dots[["ylab"]])) dots[["ylab"]] else "Standard Error"

  # set data for reference line (residuals centered at 0 or user specified)
  if (!is.null(dots[["refline"]])) {
    refline_x <- dots[["refline"]]
  } else {
    refline_x <- 0
  }
  df_refline <- data.frame(
    x = rep(refline_x, 2),
    y = range(se_range)
  )

  # create output data structures
  df_points <- data.frame(
    x                  = res,
    y                  = se,
    unit               = rep(unit, K),
    conditioning_depth = rep(conditioning_depth, K)
  )
  if (exists("res_obj") && is.data.frame(res_obj) && "cluster" %in% names(res_obj)) {
    df_points[["cluster"]]     <- res_obj[["cluster"]]
    df_points[["n_estimates"]] <- res_obj[["n_estimates"]]
  }
  # Funnel: use CLIPPED data for filled polygon
  df_funnel <- data.frame(
    x = c(rev(ci_left_clipped), ci_right_clipped),
    y = c(rev(se_sequence_clipped), se_sequence_clipped)
  )
  # Edges: use CLIPPED lines using exact intersection
  edge1 <- .clip_line_x(ci_left, se_sequence_clipped, xlim)
  edge2 <- .clip_line_x(ci_right, se_sequence_clipped, xlim)

  df_funnel_edge1 <- data.frame(
    x = edge1$x,
    y = edge1$y
  )
  df_funnel_edge2 <- data.frame(
    x = edge2$x,
    y = edge2$y
  )
  df_background <- data.frame(
    x = c(min(xlim), max(xlim), max(xlim), min(xlim)),
    y = c(min(ylim), min(ylim), max(ylim), max(ylim))
  )

  return(list(
    points       = df_points,
    funnel       = df_funnel,
    funnel_edge1 = df_funnel_edge1,
    funnel_edge2 = df_funnel_edge2,
    background   = df_background,
    x_range      = x_init_range, # ticks
    y_range      = se_range,     # ticks
    xlim         = xlim,
    ylim         = ylim,
    xlab         = xlab,
    ylab         = ylab,
    refline      = df_refline
  ))
}
