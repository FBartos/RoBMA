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
# - "bound": a Gaussian kernel times another proper scalar prior. Beyond a
#   distance r (in kernel standard deviations) from the kernel center the
#   kernel is at most exp(-r^2 / 2), so the omitted mass on each side is at most
#   exp(-r^2 / 2) times the prior mass there, relative to a lower bound on the
#   row's total mass. Each end sits where its side's bound reaches
#   (1 - normalization_prob) / 2.
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
  for (prior in unique(laws[["prior"]][bound])) {
    rows <- bound[vapply(laws[["prior"]][bound], identical, logical(1L), prior)]
    intervals[rows, ] <- .iwmde_qcmde_bound_interval(
      prior       = prior,
      center      = laws[["center"]][rows],
      scale       = laws[["scale"]][rows],
      log_mass    = log_mass[rows],
      probability = probability
    )
  }

  return(list(laws = laws, intervals = intervals, log_mass = log_mass))
}


# The tail bound of one side of a "bound" row, on the log scale: beyond b the
# kernel is at most exp(-max(0, (c - b) / s)^2 / 2) (mirrored above), so the
# omitted mass is at most that times the prior mass beyond b, relative to the
# lower bound M on the row's mass. Increasing in b below, decreasing above.
.iwmde_qcmde_bound_side_log <- function(prior, value, center, scale, log_mass,
                                        side) {

  if (identical(side, "lower")) {
    distance <- pmax(0, (center - value) / scale)
    tail     <- BayesTools::cdf(prior, value)
  } else {
    distance <- pmax(0, (value - center) / scale)
    tail     <- BayesTools::ccdf(prior, value)
  }

  return(-distance^2 / 2 + log(pmax(tail, 0)) - log_mass)
}


# The interval of "bound" rows: each end is where that side's own tail bound
# reaches (1 - probability) / 2. The symmetric radius
# r^2 = 2 (-log(1 - probability) - log M) brackets both ends (widened until each
# side meets its target) and a bisection places them, so an end stops where the
# prior's own mass runs out instead of crossing a finite support boundary.
.iwmde_qcmde_bound_interval <- function(prior, center, scale, log_mass,
                                        probability) {

  log_target <- log((1 - probability) / 2)
  radius     <- sqrt(2 * pmax(0, -log1p(-probability) - log_mass))
  side_log   <- function(value, side) {
    .iwmde_qcmde_bound_side_log(prior, value, center, scale, log_mass, side)
  }
  lower <- center - radius * scale
  upper <- center + radius * scale
  for (iteration in seq_len(60L)) {
    lower_open <- side_log(lower, "lower") > log_target
    upper_open <- side_log(upper, "upper") > log_target
    if (!any(lower_open | upper_open)) {
      break
    }
    lower[lower_open] <- center[lower_open] -
      2 * (center[lower_open] - lower[lower_open])
    upper[upper_open] <- center[upper_open] +
      2 * (upper[upper_open] - center[upper_open])
  }

  # The lower end is the largest value whose lower bound meets the target and
  # the upper end the smallest whose upper bound does, both within the bracket.
  out <- matrix(NA_real_, nrow = length(center), ncol = 2L)
  for (side in c("lower", "upper")) {
    meeting <- if (identical(side, "lower")) lower else upper
    failing <- if (identical(side, "lower")) upper else lower
    reached <- side_log(failing, side) <= log_target
    for (iteration in seq_len(80L)) {
      middle <- (meeting + failing) / 2
      meets  <- side_log(middle, side) <= log_target
      meeting[meets]  <- middle[meets]
      failing[!meets] <- middle[!meets]
    }
    meeting[reached] <- if (identical(side, "lower")) upper[reached] else
      lower[reached]
    out[, if (identical(side, "lower")) 1L else 2L] <- meeting
  }

  return(out)
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


# The tail estimate q(b) / |d log q / dz| of one end of a grid: `edge` and
# `inner` are the log densities at the end and at its inner neighbour, `width`
# the step between them, and `log_normalizer` the row's grid normalizer, all
# vectors over rows. A row whose density does not decrease towards the end has
# no estimate (Inf), a row that vanishes there has none to make (0), and a row
# without a finite normalizer is left out (NA).
.iwmde_qcmde_side_tail <- function(edge, inner, width, log_normalizer) {

  log_normalizer <- rep_len(log_normalizer, length(edge))
  slope          <- (inner - edge) / width
  estimate       <- rep(Inf, length(edge))
  estimate[edge == -Inf] <- 0
  decaying <- is.finite(edge) & is.finite(slope) & slope > 0
  estimate[decaying] <- exp(edge[decaying] - log_normalizer[decaying]) /
    slope[decaying]
  estimate[!is.finite(log_normalizer)] <- NA_real_

  return(list(estimate = estimate, slope = slope))
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
    tail  <- .iwmde_qcmde_side_tail(edge, inner, width, log_normalizer)
    out[[side]]                   <- tail[["estimate"]]
    out[[paste0(side, "_slope")]] <- tail[["slope"]]
  }

  return(out)
}


# How far, in widths of the initial range, each side may be extended for
# estimated rows; a row that has not decayed by then keeps its reported
# estimate.
.iwmde_qcmde_extension_limit <- function() 4L


.iwmde_qcmde_extension_max_iterations <- function() 40L


# Evaluate the nodes of the initial lattice and extend its sides until every
# estimated row's tail estimate is below `target` on both sides.
# `evaluate_inside()` evaluates values of the planned lattice, unless the caller
# already evaluated them (`initial_log_q`, one row per node); `evaluate_outside()`
# evaluates new values and returns NULL when they cannot be evaluated, which
# closes that side; it may return rows of missing values for the values it could
# not evaluate, which end the extension there.
#
# Each side grows by lattice steps, in as few evaluations as the rows allow: the
# tails of rows that decay at the end predict the distance the side needs
# (.iwmde_qcmde_extension_distance()), and a row that does not decay there has
# no prediction, so a side with such a row is extended by a slice of the
# initial range that grows with the steps already taken. Whatever the size of
# the slice, the side keeps only the shortest run of new nodes at which the rows
# that failed meet the target (.iwmde_qcmde_extension_prefix()); the nodes
# beyond it were evaluated in the same call and are dropped, and the next pass
# re-checks every row against the new grid.
.iwmde_qcmde_extend_lattice <- function(lattice, estimate_rows, target,
                                        evaluate_inside, evaluate_outside,
                                        initial_log_q = NULL) {

  n_initial <- lattice[["n_initial"]]
  nodes     <- .iwmde_qcmde_lattice_points(lattice, seq.int(0L, n_initial - 1L))
  if (!all(nodes[["valid"]]) || any(diff(nodes[["z"]]) <= 0)) {
    return(NULL)
  }
  log_q <- if (is.null(initial_log_q)) {
    evaluate_inside(nodes[["x"]])
  } else {
    initial_log_q
  }
  nodes[["log_q"]] <- log_q
  open     <- c(lower = TRUE, upper = TRUE)
  steps    <- c(lower = 0L, upper = 0L)
  limit    <- .iwmde_qcmde_extension_limit() * (n_initial - 1L)
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
      failing  <- which(!is.na(estimate) & estimate > target)
      if (!open[[side]] || length(failing) == 0L) {
        next
      }
      n_nodes <- length(nodes[["index"]])
      budget  <- limit - steps[[side]]
      if (budget <= 0L) {
        open[[side]] <- FALSE
        next
      }
      k <- .iwmde_qcmde_extension_size(
        distance  = .iwmde_qcmde_extension_distance(
          log_density = log_density[, failing, drop = FALSE],
          side        = side,
          step        = lattice[["step"]],
          excess      = log(estimate[failing] / target),
          slope       = tails[[paste0(side, "_slope")]][failing]
        ),
        step      = lattice[["step"]],
        n_initial = n_initial,
        taken     = steps[[side]],
        n_nodes   = n_nodes,
        budget    = budget
      )
      index <- if (identical(side, "lower")) {
        min(nodes[["index"]]) - rev(seq_len(k))
      } else {
        max(nodes[["index"]]) + seq_len(k)
      }
      extension <- .iwmde_qcmde_side_extension(lattice, index, side,
                                               evaluate_outside)
      n_valid <- length(extension[["index"]])
      if (n_valid == 0L) {
        open[[side]] <- FALSE
        next
      }
      prefix <- .iwmde_qcmde_extension_prefix(
        nodes          = nodes,
        extension      = extension,
        side           = side,
        rows           = estimate_rows[failing],
        log_normalizer = log_normalizer[failing],
        target         = target
      )
      if (n_valid < k && !prefix[["met"]]) {
        open[[side]] <- FALSE
      }
      if (prefix[["n"]] < n_valid) {
        kept <- if (identical(side, "lower")) {
          seq.int(n_valid - prefix[["n"]] + 1L, n_valid)
        } else {
          seq_len(prefix[["n"]])
        }
        extension <- .iwmde_qcmde_subset_extension(extension, kept)
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


# Distance beyond an end at which every failing row that decays there reaches
# the target: `excess` is log(estimate / target). A log density that is
# concave at the end (a Gaussian-like tail) follows its local quadratic, any
# other decaying tail its exponential envelope. A row that does not decay at
# the end has no distance; `undecided` reports whether a failing row has none.
# The next pass re-checks the estimates.
.iwmde_qcmde_extension_distance <- function(log_density, side, step, excess,
                                            slope) {

  n         <- nrow(log_density)
  decaying  <- is.finite(slope) & slope > 0
  undecided <- any(!decaying)
  if (!any(decaying)) {
    return(list(distance = 0, undecided = undecided))
  }
  inner <- if (identical(side, "lower")) seq_len(min(3L, n)) else
    rev(seq.int(max(1L, n - 2L), n))
  curvature <- if (length(inner) == 3L) {
    -(log_density[inner[[1L]], ] - 2 * log_density[inner[[2L]], ] +
        log_density[inner[[3L]], ]) / step^2
  } else {
    rep(NA_real_, ncol(log_density))
  }
  curvature <- curvature[decaying]
  excess    <- excess[decaying]
  slope     <- slope[decaying]
  distance  <- excess / slope
  concave   <- is.finite(curvature) & curvature > 0
  if (any(concave)) {
    edge_slope <- slope[concave] + curvature[concave] * step / 2
    distance[concave] <- (sqrt(edge_slope^2 + 2 * curvature[concave] *
      excess[concave]) - edge_slope) / curvature[concave]
  }
  distance <- distance[is.finite(distance)]

  return(list(distance = if (length(distance) > 0L) max(distance) else 0,
              undecided = undecided || length(distance) == 0L))
}


# Lattice steps to evaluate at once on one side: the predicted distance of the
# rows that decay at the end, and at least one slice of the initial range for a
# side with a row that does not, growing with the steps already taken so that a
# tail that never decays reaches the limit in a few calls. Never beyond the
# width of the current grid or the remaining budget.
.iwmde_qcmde_extension_size <- function(distance, step, n_initial, taken,
                                        n_nodes, budget) {

  predicted <- if (distance[["distance"]] > 0) {
    ceiling(distance[["distance"]] / step)
  } else {
    0
  }
  if (!is.finite(predicted)) {
    predicted <- n_nodes - 1L
  }
  size <- if (isTRUE(distance[["undecided"]])) {
    slice <- max(1L, ceiling((n_initial - 1L) / 8))
    max(predicted, slice, taken)
  } else {
    max(1L, predicted)
  }

  return(as.integer(max(1L, min(size, n_nodes - 1L, budget))))
}


# The shortest run of the new nodes `extension` (from
# .iwmde_qcmde_side_extension(), ascending in the lattice index) at which every
# row in `rows` meets the target. The rows are the ones that failed at the
# current end; the normalizer of each grows with the trapezoids of the nodes
# it adds. `met` is FALSE when even the whole extension does not do it.
.iwmde_qcmde_extension_prefix <- function(nodes, extension, side, rows,
                                          log_normalizer, target) {

  n_new   <- length(extension[["index"]])
  outward <- if (identical(side, "lower")) rev(seq_len(n_new)) else
    seq_len(n_new)
  edge    <- if (identical(side, "lower")) 1L else length(nodes[["index"]])
  previous_log_density <- nodes[["log_q"]][edge, rows] +
    nodes[["log_jacobian"]][[edge]]
  previous_z      <- nodes[["z"]][[edge]]
  previous_height <- exp(previous_log_density - log_normalizer)
  area            <- rep(1, length(rows))

  for (j in seq_len(n_new)) {
    node        <- outward[[j]]
    log_density <- extension[["log_q"]][node, rows] +
      extension[["log_jacobian"]][[node]]
    width  <- abs(extension[["z"]][[node]] - previous_z)
    height <- exp(log_density - log_normalizer)
    # Areas in units of the normalizer of the grid without the new nodes.
    area <- area + width * (previous_height + height) / 2
    tail <- .iwmde_qcmde_side_tail(log_density, previous_log_density, width,
                                   log_normalizer + log(area))
    if (all(is.na(tail[["estimate"]]) | tail[["estimate"]] <= target)) {
      return(list(n = j, met = TRUE))
    }
    previous_log_density <- log_density
    previous_z           <- extension[["z"]][[node]]
    previous_height      <- height
  }

  return(list(n = n_new, met = FALSE))
}


.iwmde_qcmde_subset_extension <- function(extension, keep) {

  for (field in c("index", "x", "z", "log_jacobian")) {
    extension[[field]] <- extension[[field]][keep]
  }
  extension[["log_q"]] <- extension[["log_q"]][keep, , drop = FALSE]

  return(extension)
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


# Both ends' tail estimates of every row of a normalized grid, summed.
.iwmde_qcmde_grid_tail_estimates <- function(grid) {

  tails <- .iwmde_qcmde_tail_estimates(
    z              = grid[["z"]],
    log_density    = grid[["log_q"]] + grid[["log_jacobian"]],
    log_normalizer = grid[["log_normalizer"]]
  )

  return(tails[["lower"]] + tails[["upper"]])
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
