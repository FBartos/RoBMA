# ============================================================================ #
# IWMDE Support Transforms and Grids
# ============================================================================ #

.iwmde_plot_range <- function(values, support) {

  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(NA_real_, NA_real_))
  }

  probs <- stats::quantile(values, c(.005, .995), names = FALSE, type = 8)
  if (!all(is.finite(probs)) || probs[1] == probs[2]) {
    probs <- range(values)
  }
  width <- diff(probs)
  if (!is.finite(width) || width <= 0) {
    width <- max(1, abs(probs[1]))
  }

  xlim <- probs + c(-.20, .20) * width
  xlim[1] <- max(xlim[1], support[1])
  xlim[2] <- min(xlim[2], support[2])

  if (xlim[1] == xlim[2]) {
    xlim <- xlim + c(-.5, .5)
    xlim[1] <- max(xlim[1], support[1])
    xlim[2] <- min(xlim[2], support[2])
  }

  xlim <- .iwmde_open_finite_support(xlim, support)

  return(xlim)
}


.iwmde_open_finite_support <- function(xlim, support) {

  offset <- .iwmde_interior_offsets(xlim, support)

  if (is.finite(support[1]) && xlim[1] <= support[1]) {
    xlim[1] <- support[1] + offset[1L]
  }
  if (is.finite(support[2]) && xlim[2] >= support[2]) {
    xlim[2] <- support[2] - offset[2L]
  }
  if (xlim[1] >= xlim[2] && all(is.finite(support))) {
    midpoint <- mean(support)
    radius   <- diff(support) / 4
    xlim     <- midpoint + c(-radius, radius)
  }

  return(xlim)
}


.iwmde_interior_offsets <- function(xlim, support) {

  width <- diff(xlim)
  if (!is.finite(width) || width <= 0) {
    width <- diff(support)
  }
  if (!is.finite(width) || width <= 0) {
    finite <- c(xlim, support)
    finite <- finite[is.finite(finite)]
    scale  <- if (length(finite) > 0L) max(abs(finite)) else 0
    width  <- .Machine$double.eps * scale
  }
  if (!is.finite(width) || width <= 0) {
    width <- .Machine$double.xmin * .Machine$double.eps
  }

  boundary_scale <- abs(support)
  boundary_scale[!is.finite(boundary_scale)] <- 0

  return(pmax(
    sqrt(.Machine$double.eps) * width,
    .Machine$double.eps * boundary_scale,
    .Machine$double.xmin * .Machine$double.eps
  ))
}


.iwmde_include_plot_values <- function(xlim, values, support) {

  values <- as.numeric(values)
  values <- values[
    is.finite(values) &
      values >= support[1] &
      values <= support[2]
  ]
  if (length(values) == 0L) {
    return(xlim)
  }

  xlim[1] <- min(xlim[1], values)
  xlim[2] <- max(xlim[2], values)

  return(.iwmde_open_finite_support(xlim, support))
}


.iwmde_parameter_transform <- function(support) {

  if (all(is.finite(support))) {
    return(list(
      type  = "logit",
      lower = support[1],
      upper = support[2]
    ))
  }
  if (is.finite(support[1])) {
    return(list(
      type  = "log",
      lower = support[1],
      upper = Inf
    ))
  }
  if (is.finite(support[2])) {
    return(list(
      type  = "reverse_log",
      lower = -Inf,
      upper = support[2]
    ))
  }

  return(list(
    type  = "identity",
    lower = -Inf,
    upper = Inf
  ))
}


.iwmde_to_internal <- function(x, transform) {

  if (identical(transform[["type"]], "identity")) {
    return(x)
  }
  if (identical(transform[["type"]], "log")) {
    return(log(x - transform[["lower"]]))
  }
  if (identical(transform[["type"]], "reverse_log")) {
    return(log(transform[["upper"]] - x))
  }
  if (identical(transform[["type"]], "logit")) {
    width <- transform[["upper"]] - transform[["lower"]]
    p     <- (x - transform[["lower"]]) / width
    return(stats::qlogis(p))
  }

  stop("Unknown IWMDE transform.", call. = FALSE)
}


.iwmde_from_internal <- function(z, transform) {

  if (identical(transform[["type"]], "identity")) {
    return(z)
  }
  if (identical(transform[["type"]], "log")) {
    return(transform[["lower"]] + exp(z))
  }
  if (identical(transform[["type"]], "reverse_log")) {
    return(transform[["upper"]] - exp(z))
  }
  if (identical(transform[["type"]], "logit")) {
    return(transform[["lower"]] +
      (transform[["upper"]] - transform[["lower"]]) * stats::plogis(z))
  }

  stop("Unknown IWMDE transform.", call. = FALSE)
}


.iwmde_log_jacobian <- function(z, transform) {

  if (identical(transform[["type"]], "identity")) {
    return(rep(0, length(z)))
  }
  if (identical(transform[["type"]], "log") ||
      identical(transform[["type"]], "reverse_log")) {
    return(z)
  }
  if (identical(transform[["type"]], "logit")) {
    p <- stats::plogis(z)
    return(log(transform[["upper"]] - transform[["lower"]]) +
      log(p) + log1p(-p))
  }

  stop("Unknown IWMDE transform.", call. = FALSE)
}


.iwmde_normalization_grid <- function(values, display_grid, support, transform,
                                      normalization_points,
                                      normalization_prob) {

  values <- values[is.finite(values)]
  values <- values[values > support[1] & values < support[2]]
  if (length(values) < 2L) {
    return(NULL)
  }

  z_values <- .iwmde_to_internal(values, transform)
  z_values <- z_values[is.finite(z_values)]
  if (length(z_values) < 2L) {
    return(NULL)
  }

  tail_prob <- (1 - normalization_prob) / 2
  probs     <- stats::quantile(
    z_values,
    probs = c(tail_prob, 1 - tail_prob),
    names = FALSE,
    type  = 8
  )
  if (!all(is.finite(probs)) || probs[1] == probs[2]) {
    probs <- range(z_values)
  }

  z_range <- range(probs)
  width   <- diff(probs)
  if (!is.finite(width) || width <= 0) {
    width <- diff(range(z_values))
  }
  if (!is.finite(width) || width <= 0) {
    width <- 1
  }

  z_range <- z_range + c(-.10, .10) * width
  if (z_range[1] >= z_range[2]) {
    z_range <- mean(z_range) + c(-.5, .5)
  }

  z <- seq(z_range[1], z_range[2], length.out = normalization_points)
  x <- .iwmde_from_internal(z, transform)
  offset <- .iwmde_interior_offsets(range(x), support)
  if (is.finite(support[1])) {
    x[x <= support[1]] <- support[1] + offset[1L]
  }
  if (is.finite(support[2])) {
    x[x >= support[2]] <- support[2] - offset[2L]
  }
  z <- .iwmde_to_internal(x, transform)
  if (any(!is.finite(x)) || any(!is.finite(z))) {
    return(NULL)
  }

  return(list(
    x            = x,
    z            = z,
    log_jacobian = .iwmde_log_jacobian(z, transform)
  ))
}


.iwmde_display_grid <- function(xlim, n_points, transform, values = NULL,
                                grid_method = "uniform") {

  grid_method <- .iwmde_normalize_display_grid(grid_method)
  if (identical(grid_method, "uniform") || is.null(values)) {
    return(seq(xlim[1], xlim[2], length.out = n_points))
  }

  grid <- .iwmde_adaptive_display_grid(
    xlim      = xlim,
    n_points  = n_points,
    transform = transform,
    values    = values
  )
  if (is.null(grid)) {
    return(seq(xlim[1], xlim[2], length.out = n_points))
  }

  return(grid)
}


.iwmde_include_display_values <- function(grid, values, xlim, support) {

  values <- as.numeric(values)
  values <- values[
    is.finite(values) &
      values >= xlim[1] &
      values <= xlim[2] &
      values >= support[1] &
      values <= support[2]
  ]
  if (length(values) == 0L) {
    return(grid)
  }

  for (value in sort(unique(values))) {
    if (any(grid == value)) {
      next
    }
    grid[which.min(abs(grid - value))] <- value
  }

  return(sort(unique(grid)))
}


.iwmde_adaptive_display_grid <- function(xlim, n_points, transform, values) {

  values <- values[is.finite(values)]
  values <- values[values >= xlim[1] & values <= xlim[2]]
  if (length(unique(values)) < 4L || n_points < 20L) {
    return(NULL)
  }

  n_inner    <- n_points - 2L
  n_uniform  <- max(2L, floor(.20 * n_inner))
  n_quantile <- max(2L, floor(.30 * n_inner))
  n_height   <- max(2L, floor(.25 * n_inner))
  n_curvature <- n_inner - n_uniform - n_quantile - n_height

  uniform_grid <- seq(xlim[1], xlim[2], length.out = n_uniform + 2L)
  uniform_grid <- uniform_grid[-c(1L, length(uniform_grid))]

  quantile_grid <- stats::quantile(
    values,
    probs = seq(0, 1, length.out = n_quantile + 2L)[-c(1L, n_quantile + 2L)],
    names = FALSE,
    type  = 8
  )
  curvature_grid <- .iwmde_curvature_display_grid(
    xlim        = xlim,
    transform   = transform,
    values      = values,
    n_height    = n_height,
    n_curvature = n_curvature
  )

  grid <- .iwmde_complete_display_grid(
    grid     = c(xlim[1], uniform_grid, quantile_grid, curvature_grid, xlim[2]),
    xlim     = xlim,
    n_points = n_points
  )
  if (length(grid) != n_points || any(!is.finite(grid)) ||
      any(diff(grid) <= 0)) {
    return(NULL)
  }

  return(grid)
}


.iwmde_curvature_display_grid <- function(xlim, transform, values,
                                           n_height, n_curvature) {

  if (n_height + n_curvature <= 0L) {
    return(numeric(0))
  }

  z_values <- .iwmde_to_internal(values, transform)
  z_xlim   <- range(.iwmde_to_internal(xlim, transform))
  z_values <- z_values[is.finite(z_values)]
  if (length(unique(z_values)) < 4L || any(!is.finite(z_xlim)) ||
      z_xlim[1] == z_xlim[2]) {
    return(numeric(0))
  }

  pilot <- stats::density(
    z_values,
    from = z_xlim[1],
    to   = z_xlim[2],
    n    = max(256L, 4L * (n_height + n_curvature))
  )
  if (length(pilot[["x"]]) < 4L || any(!is.finite(pilot[["y"]]))) {
    return(numeric(0))
  }

  dz <- mean(diff(pilot[["x"]]))
  if (!is.finite(dz) || dz <= 0) {
    return(numeric(0))
  }

  density_score <- .iwmde_rescale_positive(pilot[["y"]])
  second_diff   <- diff(pilot[["y"]], differences = 2L) / dz^2
  curvature     <- c(0, abs(second_diff), 0)
  curvature_score <- .iwmde_rescale_positive(curvature) * sqrt(density_score)

  z_grid <- c(
    .iwmde_weighted_display_grid(pilot[["x"]], density_score, n_height),
    .iwmde_weighted_display_grid(pilot[["x"]], curvature_score, n_curvature)
  )
  x_grid <- .iwmde_from_internal(z_grid, transform)
  x_grid <- x_grid[is.finite(x_grid) & x_grid >= xlim[1] & x_grid <= xlim[2]]

  return(x_grid)
}


.iwmde_rescale_positive <- function(x) {

  x[!is.finite(x) | x < 0] <- 0
  maximum <- max(x)
  if (!is.finite(maximum) || maximum <= 0) {
    return(rep(0, length(x)))
  }

  return(x / maximum)
}


.iwmde_weighted_display_grid <- function(grid, weights, n_points) {

  if (n_points <= 0L) {
    return(numeric(0))
  }

  keep <- is.finite(grid) & is.finite(weights) & weights > 0
  grid <- grid[keep]
  weights <- weights[keep]
  if (length(unique(grid)) < 2L || sum(weights) <= 0) {
    return(numeric(0))
  }

  order_index <- order(grid)
  grid        <- grid[order_index]
  weights     <- weights[order_index]
  cdf         <- cumsum(weights) / sum(weights)
  probs       <- seq(0, 1, length.out = n_points + 2L)
  probs       <- probs[-c(1L, length(probs))]

  return(stats::approx(cdf, grid, xout = probs, ties = "ordered",
                       rule = 2)[["y"]])
}


.iwmde_complete_display_grid <- function(grid, xlim, n_points) {

  grid <- sort(unique(as.numeric(grid)))
  grid <- grid[is.finite(grid) & grid >= xlim[1] & grid <= xlim[2]]
  grid <- sort(unique(c(xlim[1], grid, xlim[2])))

  while (length(grid) < n_points) {
    gaps <- diff(grid)
    if (length(gaps) == 0L || max(gaps) <= 0) {
      break
    }
    index <- which.max(gaps)
    grid <- sort(unique(c(grid, mean(grid[c(index, index + 1L)]))))
  }

  if (length(grid) > n_points) {
    keep <- unique(round(seq(1, length(grid), length.out = n_points)))
    grid <- grid[keep]
  }

  return(grid)
}
