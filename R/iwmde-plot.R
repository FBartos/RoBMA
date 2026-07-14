# ============================================================================ #
# IWMDE Plotting and Plot Data Helpers
# ============================================================================ #

.iwmde_kde <- function(values, xlim, n_points, mass = 1) {

  values <- values[is.finite(values)]
  if (length(unique(values)) < 2L) {
    return(NULL)
  }

  density <- stats::density(values, n = n_points, from = xlim[1], to = xlim[2])
  return(list(x = density[["x"]], y = density[["y"]] * mass))
}


.iwmde_histogram <- function(values, xlim, mass = 1) {

  values <- values[is.finite(values)]
  breaks <- seq(xlim[1], xlim[2], length.out = 18)
  mids   <- (utils::head(breaks, -1L) + utils::tail(breaks, -1L)) / 2

  if (length(values) == 0L) {
    return(list(
      mids    = mids,
      density = rep(0, length(mids)),
      breaks  = breaks
    ))
  }

  values_in_window <- values[values >= xlim[1] & values <= xlim[2]]
  if (length(values_in_window) == 0L) {
    return(list(
      mids    = mids,
      density = rep(0, length(mids)),
      breaks  = breaks
    ))
  }

  hist_data <- graphics::hist(
    values_in_window,
    breaks         = breaks,
    plot           = FALSE,
    include.lowest = TRUE
  )
  widths <- diff(hist_data[["breaks"]])
  return(list(
    mids    = hist_data[["mids"]],
    density = mass * hist_data[["counts"]] / (length(values) * widths),
    breaks  = hist_data[["breaks"]]
  ))
}


.iwmde_trapz <- function(x, y) {

  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }

  return(sum(.iwmde_trapz_weights(x) * y))
}


.iwmde_trapz_columns <- function(x, y) {

  y <- as.matrix(y)
  if (nrow(y) != length(x) || length(x) < 2L) {
    return(rep(NA_real_, ncol(y)))
  }

  return(as.numeric(crossprod(.iwmde_trapz_weights(x), y)))
}


.iwmde_trapz_weights <- function(x) {

  n <- length(x)
  if (n < 2L) {
    return(rep(0, n))
  }

  dx      <- diff(x)
  weights <- numeric(n)
  weights[-n] <- weights[-n] + dx / 2
  weights[-1L] <- weights[-1L] + dx / 2

  return(weights)
}
