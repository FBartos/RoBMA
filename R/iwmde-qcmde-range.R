# ============================================================================ #
# qCMDE Normalization Range
# ============================================================================ #
#
# Every qCMDE row density is the full conditional of the target coordinate
# given that row's other draws, normalized by a quadrature over one shared
# grid. The grid therefore has to cover each row's own conditional mass: its
# range is the union over rows of every row's central interval at
# `normalization_prob`. How the interval is found, and how well the mass left
# outside the final range (the row's truncation) is known, depends on the row:
#
# - "exact": the likelihood factor is a Gaussian kernel in the target and the
#   prior factor a (truncated) normal density, so the row is a truncated normal
#   with closed-form quantiles and truncation.
# - "bound": a Gaussian kernel times another proper scalar prior. Outside the
#   radius r (in kernel standard deviations) the kernel is at most
#   exp(-r^2 / 2), so the omitted mass is at most exp(-r^2 / 2) times the prior
#   mass there, relative to a lower bound on the row's total mass.
# - "estimate": any other row. Its interval starts at the draws-based range and
#   each side is extended until the exponential-envelope tail estimate
#   q(b) / |d log q / dz| at the endpoint is below the side's target.
#
# Rows with a retained-location conditional are normalized over their full
# support elsewhere and leave no truncation.


.iwmde_qcmde_truncation_statuses <- function() {

  return(c("exact", "bound", "estimate"))
}


# Classify the rows in `rows` from the Gaussian kernels the log-density route
# reported (NULL when it reported none). Unclassified rows are estimates.
.iwmde_qcmde_row_laws <- function(context, row_states, replacement, kernel,
                                  rows = seq_along(row_states)) {

  n    <- length(row_states)
  laws <- list(
    kind   = rep("estimate", n),
    mean   = rep(NA_real_, n),
    sd     = rep(NA_real_, n),
    lower  = rep(NA_real_, n),
    upper  = rep(NA_real_, n),
    center = rep(NA_real_, n),
    scale  = rep(NA_real_, n),
    prior  = vector("list", n)
  )
  laws[["kind"]][setdiff(seq_len(n), rows)] <- NA_character_
  if (!is.list(kernel) || length(rows) == 0L) {
    return(laws)
  }

  for (row in rows) {
    current   <- kernel[["current"]][[row]]
    linear    <- kernel[["linear"]][[row]]
    quadratic <- kernel[["quadratic"]][[row]]
    route     <- kernel[["prior_route"]][[row]]
    if (!is.finite(current) || !is.finite(linear) || !is.finite(quadratic) ||
        quadratic <= 0 || is.na(route)) {
      next
    }
    center <- current + linear / quadratic
    scale  <- 1 / sqrt(quadratic)
    if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
      next
    }

    form <- if (identical(route, "focal")) {
      .iwmde_qcmde_focal_prior_form(row_states[[row]][["focal_prior"]])
    } else if (identical(route, "linear")) {
      .iwmde_qcmde_linear_prior_form(context, row_states[[row]], replacement,
                                     current)
    } else {
      NULL
    }
    if (is.null(form)) {
      next
    }

    if (identical(form[["kind"]], "normal")) {
      # Product of the kernel and the prior's Gaussian factors in the target.
      precision <- quadratic + sum(1 / form[["sd"]]^2)
      mean      <- (quadratic * center + sum(form[["mean"]] / form[["sd"]]^2)) /
        precision
      sd        <- 1 / sqrt(precision)
      if (!is.finite(mean) || !is.finite(sd) || sd <= 0 ||
          !(form[["lower"]] < form[["upper"]])) {
        next
      }
      laws[["kind"]][[row]]  <- "exact"
      laws[["mean"]][[row]]  <- mean
      laws[["sd"]][[row]]    <- sd
      laws[["lower"]][[row]] <- form[["lower"]]
      laws[["upper"]][[row]] <- form[["upper"]]
    } else if (identical(form[["kind"]], "proper")) {
      laws[["kind"]][[row]]   <- "bound"
      laws[["center"]][[row]] <- center
      laws[["scale"]][[row]]  <- scale
      laws[["prior"]][[row]]  <- form[["prior"]]
    }
  }

  return(laws)
}


# The prior factor of a focal-prior row, as it is evaluated by
# .iwmde_focal_log_prior_values(): a normal prior.simple is a truncated normal
# density; another proper scalar prior keeps its own distribution.
.iwmde_qcmde_focal_prior_form <- function(prior) {

  normal <- .iwmde_qcmde_normal_prior_parameters(prior)
  if (!is.null(normal)) {
    return(c(list(kind = "normal"), normal))
  }
  if (.iwmde_qcmde_prior_is_bounded_proper(prior)) {
    return(list(kind = "proper", prior = prior))
  }

  return(NULL)
}


# The prior factor of a linear-target row moves every active coordinate along
# its coefficient (.iwmde_predictor_linear_log_prior_delta()). With normal
# coordinate priors each factor is a Gaussian in the target; any other prior
# leaves the row an estimate.
.iwmde_qcmde_linear_prior_form <- function(context, state, replacement,
                                           current) {

  linear <- tryCatch(
    .iwmde_linear_replacement_state(context, state, replacement),
    error = function(e) NULL
  )
  if (!isTRUE(linear[["valid"]]) ||
      length(linear[["active_columns"]]) == 0L ||
      !isTRUE(all.equal(linear[["current"]], current))) {
    return(NULL)
  }

  row   <- state[["row"]]
  mean  <- numeric()
  sd    <- numeric()
  lower <- -Inf
  upper <- Inf
  for (column in linear[["active_columns"]]) {
    prior <- .iwmde_focal_prior(context, column, row)
    if (is.null(prior) || BayesTools::is.prior.none(prior)) {
      next
    }
    normal <- .iwmde_qcmde_normal_prior_parameters(prior)
    if (is.null(normal)) {
      return(NULL)
    }
    coefficient <- linear[["coefficients"]][[column]]
    value       <- as.numeric(row[[column]])
    if (!is.finite(coefficient) || !is.finite(value)) {
      return(NULL)
    }
    if (coefficient == 0) {
      next
    }
    # The coordinate is value + (v - current) * coefficient.
    mean   <- c(mean, current + (normal[["mean"]] - value) / coefficient)
    sd     <- c(sd, normal[["sd"]] / abs(coefficient))
    limits <- sort(current + (c(normal[["lower"]], normal[["upper"]]) - value) /
      coefficient)
    lower  <- max(lower, limits[[1L]])
    upper  <- min(upper, limits[[2L]])
  }

  return(list(kind = "normal", mean = mean, sd = sd, lower = lower,
              upper = upper))
}


.iwmde_qcmde_normal_prior_parameters <- function(prior) {

  if (!inherits(prior, "prior.simple") ||
      !identical(prior[["distribution"]], "normal")) {
    return(NULL)
  }
  mean  <- as.numeric(prior[["parameters"]][["mean"]])
  sd    <- as.numeric(prior[["parameters"]][["sd"]])
  lower <- prior[["truncation"]][["lower"]]
  upper <- prior[["truncation"]][["upper"]]
  lower <- if (is.null(lower)) -Inf else as.numeric(lower)
  upper <- if (is.null(upper)) Inf else as.numeric(upper)
  if (length(mean) != 1L || length(sd) != 1L || !is.finite(mean) ||
      !is.finite(sd) || sd <= 0 || length(lower) != 1L ||
      length(upper) != 1L || is.na(lower) || is.na(upper) ||
      !(lower < upper)) {
    return(NULL)
  }

  return(list(mean = mean, sd = sd, lower = lower, upper = upper))
}


# A proper continuous scalar prior whose distribution function BayesTools
# evaluates, as .iwmde_focal_log_prior_values() evaluates its density.
.iwmde_qcmde_prior_is_bounded_proper <- function(prior) {

  if (is.null(prior) || !inherits(prior, "prior.simple") ||
      BayesTools::is.prior.none(prior) ||
      BayesTools::is.prior.point(prior) ||
      BayesTools::is.prior.discrete(prior) ||
      BayesTools::is.prior.factor(prior) ||
      !is.null(attr(prior, "multiply_by", exact = TRUE))) {
    return(FALSE)
  }
  probe <- tryCatch(
    BayesTools::cdf(prior, c(-1, 0, 1)),
    error = function(e) NULL
  )

  return(is.numeric(probe) && length(probe) == 3L && all(is.finite(probe)))
}


# Probability that a normal variable falls in (lower, upper), on the log scale,
# taken from the tail on the interval's side of the mean. Vectorized.
.iwmde_qcmde_log_normal_mass <- function(lower, upper, mean, sd) {

  n     <- max(length(lower), length(upper), length(mean), length(sd))
  lower <- rep_len(lower, n)
  upper <- rep_len(upper, n)
  a     <- (lower - rep_len(mean, n)) / rep_len(sd, n)
  b     <- (upper - rep_len(mean, n)) / rep_len(sd, n)
  out   <- rep(-Inf, n)
  valid <- !is.na(a) & !is.na(b) & a < b

  right <- valid & a > 0
  if (any(right)) {
    log_a <- stats::pnorm(a[right], lower.tail = FALSE, log.p = TRUE)
    log_b <- stats::pnorm(b[right], lower.tail = FALSE, log.p = TRUE)
    out[right] <- log_a + log1p(-exp(log_b - log_a))
  }
  left <- valid & !right
  if (any(left)) {
    log_a <- stats::pnorm(a[left], log.p = TRUE)
    log_b <- stats::pnorm(b[left], log.p = TRUE)
    out[left] <- log_b + log1p(-exp(log_a - log_b))
  }

  return(out)
}


# Lower-tail quantile of a normal truncated to (lower, upper): the q with
# P(X <= q) = probability. Vectorized over rows.
.iwmde_qcmde_truncated_normal_lower_quantile <- function(probability, mean, sd,
                                                         lower, upper) {

  n   <- length(mean)
  a   <- (lower - mean) / sd
  b   <- (upper - mean) / sd
  out <- rep(NA_real_, n)
  log_probability <- log(probability)

  right <- a > 0
  if (any(right)) {
    # P(X > q) = S(a) - probability * (S(a) - S(b)) with the survival S.
    log_sa <- stats::pnorm(a[right], lower.tail = FALSE, log.p = TRUE)
    log_sb <- stats::pnorm(b[right], lower.tail = FALSE, log.p = TRUE)
    target <- log_sa + log1p(probability * expm1(log_sb - log_sa))
    out[right] <- stats::qnorm(target, lower.tail = FALSE, log.p = TRUE)
  }
  left <- !right
  if (any(left)) {
    # P(X < q) = F(a) + probability * (F(b) - F(a)).
    log_fa   <- stats::pnorm(a[left], log.p = TRUE)
    log_mass <- .iwmde_qcmde_log_normal_mass(a[left], b[left], 0, 1)
    shifted  <- log_probability + log_mass
    top      <- pmax(log_fa, shifted)
    target   <- top + log(exp(log_fa - top) + exp(shifted - top))
    out[left] <- stats::qnorm(target, log.p = TRUE)
  }
  out <- mean + sd * out

  return(pmin(pmax(out, lower), upper))
}


.iwmde_qcmde_truncated_normal_interval <- function(tail, mean, sd, lower,
                                                   upper) {

  return(cbind(
    .iwmde_qcmde_truncated_normal_lower_quantile(tail, mean, sd, lower, upper),
    -.iwmde_qcmde_truncated_normal_lower_quantile(tail, -mean, sd, -upper,
                                                 -lower)
  ))
}


# Log lower bound on the mass M = integral of prior(v) exp(-((v - c) / s)^2 / 2)
# dv of a "bound" row: for every u, M >= exp(-u^2 / 2) P(|V - c| <= u s) under
# the prior. The bound is maximized over a fixed ladder of u.
.iwmde_qcmde_log_mass_lower_bound <- function(prior, center, scale) {

  ladder      <- seq(.25, 64, by = .25)
  radius      <- outer(scale, ladder)
  lower_value <- as.vector(center - radius)
  upper_value <- as.vector(center + radius)
  lower_cdf   <- BayesTools::cdf(prior, lower_value)
  upper_cdf   <- BayesTools::cdf(prior, upper_value)
  lower_ccdf  <- BayesTools::ccdf(prior, lower_value)
  upper_ccdf  <- BayesTools::ccdf(prior, upper_value)
  # Difference the tail that is small on this side of the prior.
  mass <- ifelse(lower_cdf > .5, lower_ccdf - upper_ccdf,
                 upper_cdf - lower_cdf)
  mass <- matrix(pmax(mass, 0), nrow = length(center))
  log_bound <- sweep(log(mass), 2L, ladder^2 / 2, "-")

  return(apply(log_bound, 1L, max))
}


# The central interval of every exact or bound row at `probability`, on the x
# scale, as a two-column matrix (NA for other rows). Bound rows also return the
# log mass bound their truncation bound divides by; a bound row whose mass
# bound vanishes has no finite radius and becomes an estimate.
.iwmde_qcmde_law_intervals <- function(laws, probability) {

  n         <- length(laws[["kind"]])
  intervals <- matrix(NA_real_, nrow = n, ncol = 2L)
  log_mass  <- rep(NA_real_, n)
  tail      <- (1 - probability) / 2

  exact <- which(laws[["kind"]] %in% "exact")
  if (length(exact) > 0L) {
    intervals[exact, ] <- .iwmde_qcmde_truncated_normal_interval(
      tail  = tail,
      mean  = laws[["mean"]][exact],
      sd    = laws[["sd"]][exact],
      lower = laws[["lower"]][exact],
      upper = laws[["upper"]][exact]
    )
  }

  bound <- which(laws[["kind"]] %in% "bound")
  for (prior in unique(laws[["prior"]][bound])) {
    rows <- bound[vapply(laws[["prior"]][bound], identical, logical(1L), prior)]
    log_mass[rows] <- .iwmde_qcmde_log_mass_lower_bound(
      prior  = prior,
      center = laws[["center"]][rows],
      scale  = laws[["scale"]][rows]
    )
  }
  unusable <- bound[!is.finite(log_mass[bound])]
  laws[["kind"]][unusable] <- "estimate"
  bound <- setdiff(bound, unusable)
  if (length(bound) > 0L) {
    # exp(-r^2 / 2) / M <= 1 - probability; M <= 1 for a proper prior.
    radius <- sqrt(2 * pmax(0, -log1p(-probability) - log_mass[bound]))
    intervals[bound, 1L] <- laws[["center"]][bound] -
      radius * laws[["scale"]][bound]
    intervals[bound, 2L] <- laws[["center"]][bound] +
      radius * laws[["scale"]][bound]
  }

  return(list(laws = laws, intervals = intervals, log_mass = log_mass))
}


# The x interval of one row, moved inside the open support, on the z chart.
.iwmde_qcmde_z_interval <- function(interval, transform) {

  support  <- c(transform[["lower"]], transform[["upper"]])
  interval <- c(max(interval[1L], support[1L]), min(interval[2L], support[2L]))
  if (any(!is.finite(interval)) || !(interval[1L] < interval[2L])) {
    return(c(NA_real_, NA_real_))
  }
  interval <- .iwmde_open_finite_support(interval, support)
  z        <- range(.iwmde_to_internal(interval, transform))
  if (any(!is.finite(z)) || !(z[1L] < z[2L])) {
    return(c(NA_real_, NA_real_))
  }

  return(z)
}


# The starting normalization range on the z chart: the union of the exact and
# bound rows' intervals, and the draws-based range when any row is an estimate.
.iwmde_qcmde_initial_range <- function(laws, intervals, base_z, transform) {

  kind   <- laws[["kind"]]
  ranges <- if (any(kind %in% "estimate") || !any(!is.na(kind))) {
    list(base_z)
  } else {
    list()
  }
  closed <- which(kind %in% c("exact", "bound"))
  for (row in closed) {
    z <- .iwmde_qcmde_z_interval(intervals[row, ], transform)
    if (all(is.finite(z))) {
      ranges[[length(ranges) + 1L]] <- z
    }
  }
  if (length(ranges) == 0L) {
    return(base_z)
  }

  return(range(unlist(ranges, use.names = FALSE)))
}


# ---------------------------------------------------------------------------- #
# Lattice of normalization nodes
# ---------------------------------------------------------------------------- #
#
# Nodes sit at z = origin + index * step. The index is an integer for the
# nodes of the selected grid and a half-integer for the midpoints of its
# nested validation grid, so every value is computed by the same expression
# wherever it is needed and repeated requests match exactly.

.iwmde_qcmde_lattice <- function(z_range, n_points, transform) {

  n_points <- as.integer(n_points)
  step     <- diff(z_range) / (n_points - 1L)

  return(list(
    origin    = z_range[[1L]],
    step      = step,
    n_initial = n_points,
    transform = transform
  ))
}


.iwmde_qcmde_lattice_points <- function(lattice, index) {

  transform <- lattice[["transform"]]
  support   <- c(transform[["lower"]], transform[["upper"]])
  x <- .iwmde_from_internal(lattice[["origin"]] + index * lattice[["step"]],
                            transform)
  z <- .iwmde_to_internal(x, transform)
  log_jacobian <- .iwmde_log_jacobian(z, transform)
  valid <- is.finite(x) & x > support[1L] & x < support[2L] &
    is.finite(z) & is.finite(log_jacobian)

  return(list(
    index        = index,
    x            = x,
    z            = z,
    log_jacobian = log_jacobian,
    valid        = valid
  ))
}


# Per-row tail estimates at the two ends of a grid: q(b) / |d log q / dz| with
# the slope of the end step, relative to the row's grid normalizer. A row whose
# density does not decrease towards an end has no estimate there (Inf).
.iwmde_qcmde_tail_estimates <- function(z, log_density, log_normalizer) {

  n     <- length(z)
  ends  <- list(lower = c(1L, 2L), upper = c(n, n - 1L))
  out   <- list()
  for (side in names(ends)) {
    edge  <- log_density[ends[[side]][[1L]], ]
    inner <- log_density[ends[[side]][[2L]], ]
    width <- abs(z[ends[[side]][[2L]]] - z[ends[[side]][[1L]]])
    slope <- (inner - edge) / width
    estimate <- rep(Inf, length(edge))
    estimate[edge == -Inf] <- 0
    decaying <- is.finite(edge) & is.finite(slope) & slope > 0
    estimate[decaying] <- exp(edge[decaying] - log_normalizer[decaying]) /
      slope[decaying]
    estimate[!is.finite(log_normalizer)] <- NA_real_
    out[[side]]                  <- estimate
    out[[paste0(side, "_slope")]] <- slope
  }

  return(out)
}


# Largest number of nodes the extension of estimated rows may reach, relative
# to the requested nodes; a row that has not decayed by then keeps its reported
# estimate.
.iwmde_qcmde_extension_limit <- function() 8L


.iwmde_qcmde_extension_max_iterations <- function() 40L


# Evaluate the nodes of the initial lattice and extend its sides, one lattice
# step at a time, until every estimated row's tail estimate is below `target`
# on both sides. `evaluate_inside()` evaluates values of the planned lattice;
# `evaluate_outside()` evaluates new values and returns NULL when they cannot
# be evaluated, which closes that side.
.iwmde_qcmde_extend_lattice <- function(lattice, estimate_rows, target,
                                        evaluate_inside, evaluate_outside) {

  n_initial <- lattice[["n_initial"]]
  nodes     <- .iwmde_qcmde_lattice_points(lattice, seq.int(0L, n_initial - 1L))
  if (!all(nodes[["valid"]]) || any(diff(nodes[["z"]]) <= 0)) {
    return(NULL)
  }
  log_q <- evaluate_inside(nodes[["x"]])
  nodes[["log_q"]] <- log_q
  open     <- c(lower = TRUE, upper = TRUE)
  steps    <- c(lower = 0L, upper = 0L)
  limit    <- .iwmde_qcmde_extension_limit() * (n_initial - 1L) + 1L
  n_passes <- 0L

  if (length(estimate_rows) == 0L) {
    return(list(nodes = nodes, steps = steps, passes = n_passes,
                open = open))
  }

  for (iteration in seq_len(.iwmde_qcmde_extension_max_iterations())) {
    log_density <- nodes[["log_q"]][, estimate_rows, drop = FALSE] +
      nodes[["log_jacobian"]]
    log_normalizer <- .iwmde_log_trapz_columns(nodes[["z"]], log_density)
    tails <- .iwmde_qcmde_tail_estimates(nodes[["z"]], log_density,
                                         log_normalizer)
    added <- FALSE
    for (side in c("lower", "upper")) {
      estimate <- tails[[side]]
      failing  <- !is.na(estimate) & estimate > target
      if (!open[[side]] || !any(failing)) {
        next
      }
      n_nodes <- length(nodes[["index"]])
      budget  <- limit - n_nodes
      if (budget <= 0L) {
        open[[side]] <- FALSE
        next
      }
      slope    <- tails[[paste0(side, "_slope")]][failing]
      distance <- log(estimate[failing] / target) / slope
      k <- if (all(is.finite(distance) & distance > 0)) {
        ceiling(max(distance) / lattice[["step"]])
      } else {
        n_nodes - 1L
      }
      k <- as.integer(max(1L, min(k, n_nodes - 1L, budget)))
      index <- if (identical(side, "lower")) {
        min(nodes[["index"]]) - rev(seq_len(k))
      } else {
        max(nodes[["index"]]) + seq_len(k)
      }
      extension <- .iwmde_qcmde_side_extension(lattice, index, side,
                                               evaluate_outside)
      if (length(extension[["index"]]) < k) {
        open[[side]] <- FALSE
      }
      if (length(extension[["index"]]) == 0L) {
        next
      }
      extended <- .iwmde_qcmde_bind_nodes(nodes, extension, side)
      if (any(diff(extended[["z"]]) <= 0)) {
        open[[side]] <- FALSE
        next
      }
      nodes <- extended
      steps[[side]] <- steps[[side]] + length(extension[["index"]])
      added <- TRUE
    }
    if (!added) {
      break
    }
    n_passes <- n_passes + 1L
  }

  return(list(nodes = nodes, steps = steps, passes = n_passes, open = open))
}


# Evaluate new outer nodes, keeping the valid run that continues the grid.
.iwmde_qcmde_side_extension <- function(lattice, index, side,
                                        evaluate_outside) {

  points <- .iwmde_qcmde_lattice_points(lattice, index)
  # Keep the nodes nearest the current grid up to the first unusable one.
  order_out <- if (identical(side, "lower")) rev(seq_along(index)) else
    seq_along(index)
  usable <- cumsum(!points[["valid"]][order_out]) == 0L
  keep   <- sort(order_out[usable])
  empty  <- list(index = integer(), x = numeric(), z = numeric(),
                 log_jacobian = numeric(), log_q = NULL)
  if (length(keep) == 0L) {
    return(empty)
  }
  points <- lapply(points[c("index", "x", "z", "log_jacobian")], `[`, keep)
  log_q  <- evaluate_outside(points[["x"]])
  if (is.null(log_q)) {
    return(empty)
  }
  # Values the joint density cannot evaluate end the extension there.
  bad <- rowSums(is.na(log_q) | log_q == Inf) > 0L
  run <- if (identical(side, "lower")) {
    rev(cumsum(rev(bad)) == 0L)
  } else {
    cumsum(bad) == 0L
  }
  if (!any(run)) {
    return(empty)
  }
  points <- lapply(points, `[`, run)
  points[["log_q"]] <- log_q[run, , drop = FALSE]

  return(points)
}


.iwmde_qcmde_bind_nodes <- function(nodes, extension, side) {

  fields <- c("index", "x", "z", "log_jacobian")
  if (identical(side, "lower")) {
    for (field in fields) {
      nodes[[field]] <- c(extension[[field]], nodes[[field]])
    }
    nodes[["log_q"]] <- rbind(extension[["log_q"]], nodes[["log_q"]])
    nodes[["valid"]] <- c(rep(TRUE, length(extension[["index"]])),
                          nodes[["valid"]])
  } else {
    for (field in fields) {
      nodes[[field]] <- c(nodes[[field]], extension[[field]])
    }
    nodes[["log_q"]] <- rbind(nodes[["log_q"]], extension[["log_q"]])
    nodes[["valid"]] <- c(nodes[["valid"]],
                          rep(TRUE, length(extension[["index"]])))
  }

  return(nodes)
}


# Rows normalized by the shared grid, as opposed to rows whose retained-location
# conditional has its own full-support normalizer.
.iwmde_qcmde_ordinary_rows <- function(row_states) {

  which(!vapply(row_states, function(state) {
    !is.null(state[["conditioning_transform"]])
  }, logical(1L)))
}


# Evaluate new values of the joint log density, or NULL when the route cannot.
.iwmde_qcmde_try_log_q <- function(context, parameter, values, row_states,
                                   replacement) {

  log_q <- tryCatch(
    .iwmde_log_q_grid(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = row_states,
      replacement = replacement
    ),
    error = function(e) NULL
  )
  if (!is.numeric(log_q) || !is.matrix(log_q) ||
      !identical(dim(log_q), c(length(values), length(row_states)))) {
    return(NULL)
  }

  return(log_q)
}


# ---------------------------------------------------------------------------- #
# Truncation
# ---------------------------------------------------------------------------- #

# Mass of every row outside the final x range: exact for truncated-normal
# rows, a bound for Gaussian-kernel rows with another proper prior, and the
# endpoint tail estimate `estimates` for the rest. Conditional rows are
# normalized over their full support and leave none.
.iwmde_qcmde_row_truncation <- function(laws, log_mass, x_range, estimates) {

  n      <- length(laws[["kind"]])
  value  <- rep(0, n)
  status <- rep("exact", n)
  lower  <- min(x_range)
  upper  <- max(x_range)

  exact <- which(laws[["kind"]] %in% "exact")
  if (length(exact) > 0L) {
    mean <- laws[["mean"]][exact]
    sd   <- laws[["sd"]][exact]
    law_lower <- laws[["lower"]][exact]
    law_upper <- laws[["upper"]][exact]
    log_total <- .iwmde_qcmde_log_normal_mass(law_lower, law_upper, mean, sd)
    log_below <- .iwmde_qcmde_log_normal_mass(law_lower,
      pmin(lower, law_upper), mean, sd)
    log_above <- .iwmde_qcmde_log_normal_mass(pmax(upper, law_lower),
      law_upper, mean, sd)
    value[exact] <- exp(log_below - log_total) + exp(log_above - log_total)
  }

  bound <- which(laws[["kind"]] %in% "bound")
  for (row in bound) {
    center <- laws[["center"]][[row]]
    scale  <- laws[["scale"]][[row]]
    prior  <- laws[["prior"]][[row]]
    # Outside (lower, upper) the kernel is at most exp(-r^2 / 2) on each side.
    below <- exp(-max(0, (center - lower) / scale)^2 / 2) *
      BayesTools::cdf(prior, lower)
    above <- exp(-max(0, (upper - center) / scale)^2 / 2) *
      BayesTools::ccdf(prior, upper)
    value[[row]]  <- (below + above) / exp(log_mass[[row]])
    status[[row]] <- "bound"
  }

  estimate <- which(laws[["kind"]] %in% "estimate")
  if (length(estimate) > 0L) {
    value[estimate]  <- estimates[estimate]
    status[estimate] <- "estimate"
  }
  value <- pmin(value, 1)

  return(list(value = value, status = status))
}


# The largest row truncation, the weakest status among the rows, and the
# relative bound it places on a normalized ordinate: a row normalizer that
# misses the fraction t of its mass overstates that row's density by the
# factor 1 / (1 - t).
.iwmde_qcmde_truncation_summary <- function(truncation) {

  value  <- truncation[["value"]]
  status <- truncation[["status"]]
  if (length(value) == 0L) {
    return(list(max = 0, status = "exact", ordinate_bound = 0))
  }
  statuses <- .iwmde_qcmde_truncation_statuses()
  weakest  <- statuses[[max(match(status, statuses))]]
  largest  <- if (anyNA(value)) NA_real_ else max(value)
  bound    <- if (is.na(largest)) {
    NA_real_
  } else if (largest >= 1) {
    Inf
  } else {
    largest / (1 - largest)
  }

  return(list(max = largest, status = weakest, ordinate_bound = bound))
}
