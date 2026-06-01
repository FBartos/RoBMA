# ---------------------------------------------------------------------------- #
# Plotting helpers
# ---------------------------------------------------------------------------- #
.regplot_has_shade <- function(shade) {
  return(isTRUE(shade))
}

.regplot_palette <- function(groups, default_col) {

  n <- length(groups)
  if (n == 1L) {
    return(stats::setNames(default_col[1], groups))
  }

  if (length(default_col) >= n && length(unique(default_col)) > 1L) {
    cols <- default_col[seq_len(n)]
  } else {
    cols <- grDevices::hcl.colors(n, palette = "Dark 3")
  }

  return(stats::setNames(cols, groups))
}

.regplot_dodge_x <- function(x, group_id, n_groups, dodge_width) {

  if (n_groups <= 1L) {
    return(x)
  }

  offset <- (group_id - (n_groups + 1) / 2) * dodge_width / n_groups
  return(x + offset)
}

.regplot_jitter_values <- function(n, amount) {

  if (amount <= 0) {
    return(rep(0, n))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  set.seed(42)
  return(stats::runif(n, -amount, amount))
}

.regplot_band_edges <- function(df_band) {

  out <- do.call(rbind, lapply(split(df_band, df_band[["group"]]), function(df_group) {
    n <- nrow(df_group) / 2
    data.frame(
      x        = rep(df_group[["x"]][seq_len(n)], 2),
      y        = c(df_group[["lower"]][seq_len(n)], df_group[["upper"]][seq_len(n)]),
      edge     = rep(c("lower", "upper"), each = n),
      group    = rep(df_group[["group"]][1], 2 * n),
      group_id = rep(df_group[["group_id"]][1], 2 * n),
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  return(out)
}

.regplot_draw_continuous_band_base <- function(df_band, fill_col, alpha,
                                               border_col, shade, lwd) {

  if (is.null(df_band)) {
    return(invisible(NULL))
  }

  for (group in unique(df_band[["group"]])) {
    df_group <- df_band[df_band[["group"]] == group, , drop = FALSE]
    n        <- nrow(df_group) / 2

    if (.regplot_has_shade(shade)) {
      graphics::polygon(
        df_group[["x"]],
        df_group[["y"]],
        col    = grDevices::adjustcolor(fill_col[group], alpha.f = alpha),
        border = NA
      )
    }

    graphics::lines(
      df_group[["x"]][seq_len(n)],
      df_group[["lower"]][seq_len(n)],
      col = border_col[group],
      lwd = max(1, lwd / 2)
    )
    graphics::lines(
      df_group[["x"]][seq_len(n)],
      df_group[["upper"]][seq_len(n)],
      col = border_col[group],
      lwd = max(1, lwd / 2)
    )
  }

  return(invisible(NULL))
}

.regplot_draw_categorical_band_base <- function(df_band, fill_col, alpha,
                                                border_col, shade, width,
                                                n_groups, dodge_width, lwd) {

  if (is.null(df_band)) {
    return(invisible(NULL))
  }

  x <- .regplot_dodge_x(df_band[["x"]], df_band[["group_id"]], n_groups, dodge_width)
  w <- width / max(1, n_groups)

  for (i in seq_len(nrow(df_band))) {
    group <- df_band[["group"]][i]
    graphics::rect(
      xleft   = x[i] - w / 2,
      ybottom = df_band[["lower"]][i],
      xright  = x[i] + w / 2,
      ytop    = df_band[["upper"]][i],
      col     = if (.regplot_has_shade(shade)) grDevices::adjustcolor(fill_col[group], alpha.f = alpha) else NA,
      border  = border_col[group],
      lwd     = max(1, lwd / 2)
    )
    graphics::segments(
      x0  = x[i] - w / 2,
      y0  = df_band[["middle"]][i],
      x1  = x[i] + w / 2,
      y1  = df_band[["middle"]][i],
      col = border_col[group],
      lwd = max(1, lwd / 2)
    )
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
.regplot_plot_base <- function(data, dots, legend = TRUE) {

  pch       <- dots[["pch"]]
  col       <- dots[["col"]]
  bg        <- dots[["bg"]]
  lcol      <- dots[["lcol"]]
  lwd       <- dots[["lwd"]]
  col_ci    <- dots[["col.ci"]]
  col_pi    <- dots[["col.pi"]]
  col_si    <- dots[["col.si"]]
  alpha_ci  <- dots[["alpha.ci"]]
  alpha_pi  <- dots[["alpha.pi"]]
  alpha_si  <- dots[["alpha.si"]]
  shade     <- dots[["shade"]]
  main      <- dots[["main"]]
  jitter_am <- dots[["jitter"]]
  box_width <- dots[["box.width"]]
  dodge_w   <- dots[["dodge.width"]]
  las       <- dots[["las"]]

  df_points  <- data$points
  df_pred    <- data$pred
  df_ci      <- data$ci
  df_pi      <- data$pi
  df_si      <- data$si
  df_refline <- data$refline
  mod_type   <- data$mod_type
  groups     <- data$groups
  n_groups   <- length(groups)
  has_by     <- n_groups > 1L

  line_cols <- .regplot_palette(groups, lcol)
  point_cols <- if (has_by) line_cols else stats::setNames(rep(col, n_groups), groups)
  point_bgs  <- if (has_by) {
    stats::setNames(grDevices::adjustcolor(line_cols, alpha.f = 0.45), groups)
  } else {
    stats::setNames(rep(bg, n_groups), groups)
  }

  if (has_by) {
    ci_cols <- pi_cols <- si_cols <- line_cols
  } else {
    ci_cols <- stats::setNames(rep(col_ci, n_groups), groups)
    pi_cols <- stats::setNames(rep(col_pi, n_groups), groups)
    si_cols <- stats::setNames(rep(col_si, n_groups), groups)
  }

  if (mod_type == "categorical") {
    graphics::plot(
      NA, NA,
      xlim = data$xlim,
      ylim = data$ylim,
      xlab = data$xlab,
      ylab = data$ylab,
      main = main,
      type = "n",
      xaxt = "n",
      las  = las
    )
    graphics::axis(
      1,
      at     = seq_along(data$levels),
      labels = data$levels,
      las    = las
    )
  } else {
    graphics::plot(
      NA, NA,
      xlim = data$xlim,
      ylim = data$ylim,
      xlab = data$xlab,
      ylab = data$ylab,
      main = main,
      type = "n",
      las  = las
    )
  }

  if (!is.null(df_refline)) {
    graphics::abline(h = df_refline$y, lty = "dashed", col = "gray50")
  }

  if (!is.null(df_si)) {
    if (mod_type == "continuous") {
      .regplot_draw_continuous_band_base(df_si, si_cols, alpha_si, si_cols, shade, lwd)
    } else {
      .regplot_draw_categorical_band_base(
        df_si, si_cols, alpha_si, si_cols, shade,
        box_width * 1.4, n_groups, dodge_w, lwd
      )
    }
  }

  if (!is.null(df_pi)) {
    if (mod_type == "continuous") {
      .regplot_draw_continuous_band_base(df_pi, pi_cols, alpha_pi, pi_cols, shade, lwd)
    } else {
      .regplot_draw_categorical_band_base(
        df_pi, pi_cols, alpha_pi, pi_cols, shade,
        box_width * 1.15, n_groups, dodge_w, lwd
      )
    }
  }

  if (!is.null(df_ci)) {
    if (mod_type == "continuous") {
      .regplot_draw_continuous_band_base(df_ci, ci_cols, alpha_ci, ci_cols, shade, lwd)
    } else {
      .regplot_draw_categorical_band_base(
        df_ci, ci_cols, alpha_ci, ci_cols, shade,
        box_width, n_groups, dodge_w, lwd
      )
    }
  }

  if (!is.null(df_pred)) {
    if (mod_type == "continuous") {
      for (group in groups) {
        df_group <- df_pred[df_pred[["group"]] == group, , drop = FALSE]
        graphics::lines(df_group$x, df_group$y, col = line_cols[group], lwd = lwd)
      }
    } else {
      x_pred <- .regplot_dodge_x(df_pred$x, df_pred$group_id, n_groups, dodge_w)
      graphics::points(
        x_pred,
        df_pred$y,
        pch = 18,
        col = line_cols[df_pred$group],
        cex = 1.5
      )
    }
  }

  if (mod_type == "categorical") {
    x_points <- .regplot_dodge_x(df_points$x, df_points$group_id, n_groups, dodge_w)
    x_points <- x_points + .regplot_jitter_values(
      n      = nrow(df_points),
      amount = jitter_am / max(1, n_groups)
    )
    graphics::points(
      x_points,
      df_points$y,
      pch = pch,
      col = point_cols[df_points$group],
      bg  = point_bgs[df_points$group],
      cex = df_points$size
    )
  } else {
    graphics::points(
      df_points$x,
      df_points$y,
      pch = pch,
      col = point_cols[df_points$group],
      bg  = point_bgs[df_points$group],
      cex = df_points$size
    )
  }

  if (has_by && legend) {
    graphics::legend(
      "topright",
      legend = groups,
      col    = line_cols,
      pt.bg  = point_bgs,
      pch    = 21,
      lty    = 1,
      lwd    = lwd,
      bty    = "n"
    )
  }

  return(invisible(NULL))
}
