.regplot_add_continuous_band_ggplot <- function(out, df_band, fill_col, alpha,
                                                border_col, shade, has_by) {

  if (is.null(df_band)) {
    return(out)
  }

  if (.regplot_has_shade(shade)) {
    if (has_by) {
      out <- out +
        ggplot2::geom_polygon(
          data    = df_band,
          mapping = ggplot2::aes(
            x     = .data[["x"]],
            y     = .data[["y"]],
            group = .data[["group"]],
            fill  = .data[["group"]]
          ),
          alpha   = alpha
        )
    } else {
      out <- out +
        ggplot2::geom_polygon(
          data    = df_band,
          mapping = ggplot2::aes(
            x     = .data[["x"]],
            y     = .data[["y"]],
            group = .data[["group"]]
          ),
          fill    = fill_col[1],
          alpha   = alpha
        )
    }
  }

  df_edges <- .regplot_band_edges(df_band)
  if (has_by) {
    out <- out +
      ggplot2::geom_line(
        data    = df_edges,
        mapping = ggplot2::aes(
          x      = .data[["x"]],
          y      = .data[["y"]],
          group  = interaction(.data[["group"]], .data[["edge"]]),
          colour = .data[["group"]]
        ),
        linewidth = 0.35
      )
  } else {
    out <- out +
      ggplot2::geom_line(
        data    = df_edges,
        mapping = ggplot2::aes(
          x     = .data[["x"]],
          y     = .data[["y"]],
          group = .data[["edge"]]
        ),
        colour  = border_col[1],
        linewidth = 0.35
      )
  }

  return(out)
}

.regplot_add_categorical_band_ggplot <- function(out, df_band, fill_col, alpha,
                                                 border_col, shade, has_by,
                                                 width, n_groups, dodge_width) {

  if (is.null(df_band)) {
    return(out)
  }

  df_band[["x_plot"]] <- .regplot_dodge_x(
    df_band[["x"]],
    df_band[["group_id"]],
    n_groups,
    dodge_width
  )
  width <- width / max(1, n_groups)

  if (has_by && .regplot_has_shade(shade)) {
    out <- out +
      ggplot2::geom_crossbar(
        data    = df_band,
        mapping = ggplot2::aes(
          x      = .data[["x_plot"]],
          y      = .data[["middle"]],
          ymin   = .data[["lower"]],
          ymax   = .data[["upper"]],
          colour = .data[["group"]],
          fill   = .data[["group"]]
        ),
        width   = width,
        alpha   = alpha,
        linewidth = 0.35
      )
  } else if (has_by) {
    out <- out +
      ggplot2::geom_crossbar(
        data    = df_band,
        mapping = ggplot2::aes(
          x      = .data[["x_plot"]],
          y      = .data[["middle"]],
          ymin   = .data[["lower"]],
          ymax   = .data[["upper"]],
          colour = .data[["group"]]
        ),
        width   = width,
        fill    = NA,
        alpha   = 1,
        linewidth = 0.35
      )
  } else {
    out <- out +
      ggplot2::geom_crossbar(
        data    = df_band,
        mapping = ggplot2::aes(
          x    = .data[["x_plot"]],
          y    = .data[["middle"]],
          ymin = .data[["lower"]],
          ymax = .data[["upper"]]
        ),
        width   = width,
        colour  = border_col[1],
        fill    = if (.regplot_has_shade(shade)) fill_col[1] else NA,
        alpha   = if (.regplot_has_shade(shade)) alpha else 1,
        linewidth = 0.35
      )
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .regplot_plot_base
# ---------------------------------------------------------------------------- #
#
# Create regression plot using base R graphics.
#
# @param data   list of plot data from .regplot_data
# @param dots   list of graphical parameters
# @param legend logical; show legend
#
# @return NULL invisibly
#

# ---------------------------------------------------------------------------- #
# .regplot_plot_ggplot
# ---------------------------------------------------------------------------- #
#
# Create regression plot using ggplot2.
#
# @param data   list of plot data from .regplot_data
# @param dots   list of graphical parameters
# @param legend logical; show legend
#
# @return ggplot object
#
# ---------------------------------------------------------------------------- #
.regplot_plot_ggplot <- function(data, dots, legend = TRUE) {

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
  if (has_by) {
    ci_cols <- pi_cols <- si_cols <- line_cols
  } else {
    ci_cols <- stats::setNames(rep(col_ci, n_groups), groups)
    pi_cols <- stats::setNames(rep(col_pi, n_groups), groups)
    si_cols <- stats::setNames(rep(col_si, n_groups), groups)
  }

  out <- ggplot2::ggplot()

  if (!is.null(df_refline)) {
    out <- out +
      ggplot2::geom_hline(
        yintercept = df_refline$y,
        linetype   = "dashed",
        colour     = "gray50"
      )
  }

  if (!is.null(df_si)) {
    if (mod_type == "continuous") {
      out <- .regplot_add_continuous_band_ggplot(out, df_si, si_cols, alpha_si, si_cols, shade, has_by)
    } else {
      out <- .regplot_add_categorical_band_ggplot(
        out, df_si, si_cols, alpha_si, si_cols, shade, has_by,
        box_width * 1.4, n_groups, dodge_w
      )
    }
  }

  if (!is.null(df_pi)) {
    if (mod_type == "continuous") {
      out <- .regplot_add_continuous_band_ggplot(out, df_pi, pi_cols, alpha_pi, pi_cols, shade, has_by)
    } else {
      out <- .regplot_add_categorical_band_ggplot(
        out, df_pi, pi_cols, alpha_pi, pi_cols, shade, has_by,
        box_width * 1.15, n_groups, dodge_w
      )
    }
  }

  if (!is.null(df_ci)) {
    if (mod_type == "continuous") {
      out <- .regplot_add_continuous_band_ggplot(out, df_ci, ci_cols, alpha_ci, ci_cols, shade, has_by)
    } else {
      out <- .regplot_add_categorical_band_ggplot(
        out, df_ci, ci_cols, alpha_ci, ci_cols, shade, has_by,
        box_width, n_groups, dodge_w
      )
    }
  }

  if (!is.null(df_pred)) {
    if (mod_type == "continuous") {
      if (has_by) {
        out <- out +
          ggplot2::geom_line(
            data    = df_pred,
            mapping = ggplot2::aes(
              x      = .data[["x"]],
              y      = .data[["y"]],
              colour = .data[["group"]],
              group  = .data[["group"]]
            ),
            linewidth = lwd / 2
          )
      } else {
        out <- out +
          ggplot2::geom_line(
            data    = df_pred,
            mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
            colour  = lcol[1],
            linewidth = lwd / 2
          )
      }
    } else {
      df_pred[["x_plot"]] <- .regplot_dodge_x(df_pred$x, df_pred$group_id, n_groups, dodge_w)
      if (has_by) {
        out <- out +
          ggplot2::geom_point(
            data    = df_pred,
            mapping = ggplot2::aes(
              x      = .data[["x_plot"]],
              y      = .data[["y"]],
              colour = .data[["group"]]
            ),
            shape   = 18,
            size    = 4
          )
      } else {
        out <- out +
          ggplot2::geom_point(
            data    = df_pred,
            mapping = ggplot2::aes(x = .data[["x_plot"]], y = .data[["y"]]),
            shape   = 18,
            colour  = lcol[1],
            size    = 4
          )
      }
    }
  }

  if (mod_type == "categorical") {
    df_points[["x_plot"]] <- .regplot_dodge_x(df_points$x, df_points$group_id, n_groups, dodge_w)
    df_points[["x_plot"]] <- df_points[["x_plot"]] + .regplot_jitter_values(
      n      = nrow(df_points),
      amount = jitter_am / max(1, n_groups)
    )
  } else {
    df_points[["x_plot"]] <- df_points[["x"]]
  }

  if (has_by) {
    out <- out +
      ggplot2::geom_point(
        data    = df_points,
        mapping = ggplot2::aes(
          x      = .data[["x_plot"]],
          y      = .data[["y"]],
          size   = .data[["size"]],
          colour = .data[["group"]],
          fill   = .data[["group"]]
        ),
        shape   = pch
      )
  } else {
    out <- out +
      ggplot2::geom_point(
        data    = df_points,
        mapping = ggplot2::aes(
          x    = .data[["x_plot"]],
          y    = .data[["y"]],
          size = .data[["size"]]
        ),
        shape   = pch,
        colour  = col,
        fill    = bg
      )
  }

  out <- out +
    ggplot2::scale_size_identity()

  if (mod_type == "categorical") {
    out <- out +
      ggplot2::scale_x_continuous(
        breaks = seq_along(data$levels),
        labels = data$levels,
        limits = data$xlim,
        name   = data$xlab
      )
  } else {
    out <- out +
      ggplot2::scale_x_continuous(
        limits = data$xlim,
        name   = data$xlab
      )
  }

  out <- out +
    ggplot2::scale_y_continuous(
      limits = data$ylim,
      name   = data$ylab
    )

  if (has_by) {
    out <- out +
      ggplot2::scale_colour_manual(values = line_cols, name = data$by_name) +
      ggplot2::scale_fill_manual(values = line_cols, name = data$by_name)
  }

  if (!is.null(main)) {
    out <- out + ggplot2::ggtitle(main)
  }

  out <- out + ggplot2::guides(size = "none")
  if (has_by && !legend) {
    out <- out + ggplot2::guides(colour = "none", fill = "none")
  }

  return(out)
}

