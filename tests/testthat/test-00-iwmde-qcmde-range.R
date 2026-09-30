context("(00) qCMDE per-row normalization range")

# References: the truncated-normal quantiles and interval masses, and the
# Gaussian-kernel x Cauchy-prior masses, were computed with mpmath at 60
# significant digits (bisection on the normal distribution function and
# adaptive quadrature over the whole line).

.qcmde_range_state <- function(prior) {

  list(focal_prior = prior, use_focal_prior_delta = TRUE)
}


.qcmde_range_kernel <- function(center, scale, current = center,
                                route = "focal") {

  quadratic <- 1 / scale^2
  list(
    current     = current,
    linear      = (center - current) * quadratic,
    quadratic   = rep(quadratic, length.out = length(current)),
    prior_route = rep(route, length(current))
  )
}


test_that("truncated-normal intervals and masses match mpmath references", {

  cases <- list(
    list(tail = 1e-10, mean = 0, sd = 1, lower = 2, upper = Inf,
         expected = c(2.0000000000421369229, 6.9189502390077374593)),
    list(tail = 5e-4, mean = .3, sd = .2, lower = -Inf, upper = -1,
         expected = c(-1.212095000777065775, -1.0000150473776711931)),
    list(tail = 5e-9, mean = -.25, sd = 1.5, lower = -1, upper = 3,
         expected = c(-0.99999998559218414401, 2.9999998670473625082)),
    list(tail = 5e-4, mean = .1, sd = .05, lower = 0, upper = Inf,
         expected = c(0.00044846569708682588706, 0.26484978652029922883))
  )
  for (case in cases) {
    interval <- .iwmde_qcmde_truncated_normal_interval(
      case[["tail"]], case[["mean"]], case[["sd"]], case[["lower"]],
      case[["upper"]]
    )
    expect_lt(max(abs(interval[1L, ] - case[["expected"]])), 1e-12)
  }

  masses <- .iwmde_qcmde_log_normal_mass(
    lower = c(40, -3, -1, 8),
    upper = c(41, -2.5, 1e6, 9),
    mean  = 0,
    sd    = 1
  )
  expected <- c(-804.60844201375378817, -5.3267647240790622482,
                -0.17275377902344988953, -35.013618593437148117)
  expect_lt(max(abs(masses / expected - 1)), 1e-12)
  expect_identical(.iwmde_qcmde_log_normal_mass(1, 1, 0, 1), -Inf)
})


test_that("a Gaussian kernel with a normal prior is an exact truncated normal", {

  prior  <- BayesTools::prior("normal", list(mean = .2, sd = .5),
                              truncation = list(lower = 0, upper = Inf))
  center <- c(.35, -.1)
  scale  <- c(.08, .3)
  states <- lapply(center, function(value) .qcmde_range_state(prior))
  kernel <- .qcmde_range_kernel(center, scale, current = c(.3, .05))

  laws <- .iwmde_qcmde_row_laws(list(), states, list(type = "scalar"), kernel)
  precision <- 1 / scale^2 + 1 / .5^2
  expect_identical(laws[["kind"]], c("exact", "exact"))
  expect_equal(laws[["sd"]], 1 / sqrt(precision), tolerance = 1e-14)
  expect_equal(laws[["mean"]], (center / scale^2 + .2 / .5^2) / precision,
               tolerance = 1e-14)
  expect_identical(laws[["lower"]], c(0, 0))
  expect_identical(laws[["upper"]], c(Inf, Inf))

  probability <- .999
  intervals <- .iwmde_qcmde_law_intervals(laws, probability)
  truncation <- vapply(seq_along(center), function(row) {
    .iwmde_qcmde_row_truncation(
      laws, intervals[["log_mass"]], intervals[["intervals"]][row, ],
      estimates = rep(NA_real_, 2L)
    )[["value"]][[row]]
  }, numeric(1L))
  # Each row's own central interval leaves exactly 1 - probability outside.
  expect_equal(truncation, rep(1 - probability, 2L), tolerance = 1e-10)

  # Direct normal-distribution arithmetic for a wider range.
  x_range <- c(.05, .9)
  value <- .iwmde_qcmde_row_truncation(laws, intervals[["log_mass"]], x_range,
                                       rep(NA_real_, 2L))[["value"]]
  total <- 1 - stats::pnorm(0, laws[["mean"]], laws[["sd"]])
  outside <- stats::pnorm(x_range[[1L]], laws[["mean"]], laws[["sd"]]) -
    stats::pnorm(0, laws[["mean"]], laws[["sd"]]) +
    stats::pnorm(x_range[[2L]], laws[["mean"]], laws[["sd"]],
                 lower.tail = FALSE)
  expect_equal(value, outside / total, tolerance = 1e-12)
})


test_that("a Gaussian kernel with a Cauchy prior keeps a valid mass bound", {

  prior  <- BayesTools::prior("cauchy", list(location = 0, scale = .5))
  center <- .8
  scale  <- .15
  laws   <- .iwmde_qcmde_row_laws(list(), list(.qcmde_range_state(prior)),
                                  list(type = "scalar"),
                                  .qcmde_range_kernel(center, scale))
  expect_identical(laws[["kind"]], "bound")
  expect_equal(laws[["center"]], center)
  expect_equal(laws[["scale"]], scale)

  # M = integral of prior(v) exp(-((v - center) / scale)^2 / 2) dv.
  mass <- 0.070483749361896878824
  intervals <- .iwmde_qcmde_law_intervals(laws, .999)
  expect_lte(intervals[["log_mass"]], log(mass))
  expect_gt(intervals[["log_mass"]], log(mass) - log(2.5))

  outside <- c(`3` = 0.0039233162877938912104,
               `4.5` = 0.000012254571437277958002)
  for (radius in c(3, 4.5)) {
    bound <- .iwmde_qcmde_row_truncation(
      laws, intervals[["log_mass"]], center + c(-1, 1) * radius * scale,
      estimates = NA_real_
    )
    expect_identical(bound[["status"]], "bound")
    expect_gte(bound[["value"]], outside[[as.character(radius)]])
  }
  own <- .iwmde_qcmde_row_truncation(laws, intervals[["log_mass"]],
                                     intervals[["intervals"]][1L, ],
                                     estimates = NA_real_)
  expect_lte(own[["value"]], 1 - .999 + 1e-12)
})


test_that("a bound row at a finite support boundary ends where its side bound meets the target", {

  # A Gaussian kernel times a half-Cauchy prior on (0, Inf), with the kernel
  # centre near and below the boundary. The symmetric radius c -/+ r s crosses
  # 0; each end must instead sit where its own side's bound reaches
  # (1 - p) / 2, which stays inside the support. The reference masses come
  # from adaptive quadrature of prior(v) exp(-((v - c) / s)^2 / 2).
  prior <- BayesTools::prior("cauchy", list(location = 0, scale = .5),
                             truncation = list(lower = 0, upper = Inf))
  probability <- .999
  target      <- (1 - probability) / 2
  for (center in c(.05, -.1)) {
    scale <- .2
    laws  <- .iwmde_qcmde_row_laws(list(), list(.qcmde_range_state(prior)),
                                   list(type = "scalar"),
                                   .qcmde_range_kernel(center, scale))
    expect_identical(laws[["kind"]], "bound")
    intervals <- .iwmde_qcmde_law_intervals(laws, probability)
    interval  <- intervals[["intervals"]][1L, ]
    radius    <- sqrt(2 * (-log1p(-probability) - intervals[["log_mass"]]))
    expect_lt(center - radius * scale, 0)
    expect_gt(interval[[1L]], 1e-5)
    expect_lt(interval[[1L]], interval[[2L]])

    density <- function(v) {
      2 * stats::dcauchy(v, 0, .5) * exp(-((v - center) / scale)^2 / 2)
    }
    mass  <- function(lower, upper) {
      stats::integrate(density, lower, upper, rel.tol = 1e-12,
                       subdivisions = 1000L)$value
    }
    total <- mass(0, Inf)
    below <- mass(0, interval[[1L]]) / total
    above <- mass(interval[[2L]], Inf) / total
    expect_lte(below, target)
    expect_lte(above, target)
    # The lower end is not needlessly close to the boundary: the side bound
    # is within the looseness of the mass bound of the true lower mass.
    expect_gt(below, target / 20)

    # The reported bound over the row's own interval keeps the guarantee.
    own <- .iwmde_qcmde_row_truncation(laws, intervals[["log_mass"]], interval,
                                       estimates = NA_real_)
    expect_identical(own[["status"]], "bound")
    expect_gte(own[["value"]], below + above)
    expect_lte(own[["value"]], 1 - probability + 1e-12)
  }
})


test_that("rows without a Gaussian kernel or a supported prior are estimates", {

  normal <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  point  <- BayesTools::prior("point", list(location = 0))
  states <- list(.qcmde_range_state(normal), .qcmde_range_state(point),
                 .qcmde_range_state(normal), .qcmde_range_state(normal))
  kernel <- .qcmde_range_kernel(rep(0, 4L), rep(1, 4L))
  kernel[["quadratic"]][[3L]]   <- -1
  kernel[["prior_route"]][[4L]] <- "generic"
  laws <- .iwmde_qcmde_row_laws(list(), states, list(type = "scalar"), kernel,
                                rows = 1:4)
  expect_identical(laws[["kind"]], c("exact", "estimate", "estimate",
                                     "estimate"))

  none <- .iwmde_qcmde_row_laws(list(), states, list(type = "scalar"), NULL,
                                rows = 2:3)
  expect_identical(none[["kind"]], c(NA, "estimate", "estimate", NA))
})


test_that("the initial range adds the draws range only for estimated rows", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  prior     <- BayesTools::prior("normal", list(mean = 0, sd = 10))
  states    <- list(.qcmde_range_state(prior), .qcmde_range_state(prior))
  kernel    <- .qcmde_range_kernel(c(-1, 2), c(.1, .2))
  laws      <- .iwmde_qcmde_row_laws(list(), states, list(type = "scalar"),
                                     kernel)
  intervals <- .iwmde_qcmde_law_intervals(laws, .99)
  union     <- range(intervals[["intervals"]])

  exact <- .iwmde_qcmde_initial_range(laws, intervals[["intervals"]],
                                      base_z = c(-5, 5), transform = transform)
  expect_equal(exact, union)

  laws[["kind"]][[2L]] <- "estimate"
  mixed <- .iwmde_qcmde_initial_range(laws, intervals[["intervals"]],
                                      base_z = c(-5, 5), transform = transform)
  expect_equal(mixed, c(-5, 5))

  # On a log chart the interval is mapped through the chart.
  log_transform <- .iwmde_parameter_transform(c(0, Inf))
  positive <- list(kind = "exact", mean = 1, sd = .1, lower = 0, upper = Inf)
  interval <- .iwmde_qcmde_law_intervals(positive, .99)[["intervals"]]
  expect_equal(.iwmde_qcmde_initial_range(positive, interval, c(-9, 9),
                                          log_transform),
               log(interval[1L, ]))
})


test_that("the envelope tail estimate is exact for an exponential tail", {

  z <- seq(-2, 2.5, length.out = 46L)
  log_density <- cbind(-abs(z - .3) / .5)
  # The density integrates to 1 over the whole line.
  tails <- .iwmde_qcmde_tail_estimates(z, log_density, log_normalizer = 0)
  expect_equal(tails[["lower"]], .5 * exp(-2.3 / .5), tolerance = 1e-12)
  expect_equal(tails[["upper"]], .5 * exp(-2.2 / .5), tolerance = 1e-12)

  flat <- .iwmde_qcmde_tail_estimates(z, cbind(rep(0, length(z))), 0)
  expect_identical(flat[["lower"]], Inf)
  empty <- .iwmde_qcmde_tail_estimates(z, cbind(c(-Inf, rep(0, 45L))), 0)
  expect_identical(empty[["lower"]], 0)
})


test_that("estimated rows extend the grid by lattice steps to the target", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  laplace   <- function(values) {
    cbind(-abs(values - .3) / .5, -abs(values + .2) / .25)
  }
  evaluated <- numeric()
  evaluate  <- function(values) {
    evaluated <<- c(evaluated, values)
    laplace(values)
  }
  lattice <- .iwmde_qcmde_lattice(c(-1, 1), 21L, transform)
  target  <- 1e-6
  result  <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1:2,
    target           = target,
    evaluate_inside  = evaluate,
    evaluate_outside = evaluate
  )
  nodes <- result[["nodes"]]
  # The nodes stay one uniform lattice that contains the initial grid, and
  # every value was evaluated once.
  expect_equal(diff(nodes[["z"]]), rep(.1, length(nodes[["z"]]) - 1L),
               tolerance = 1e-12)
  expect_true(all(round(seq(-1, 1, by = .1), 12) %in% round(nodes[["z"]], 12)))
  expect_equal(sort(evaluated), nodes[["x"]])
  expect_gt(result[["steps"]][["lower"]], 0L)
  expect_gt(result[["steps"]][["upper"]], 0L)

  # The analytic tail mass of both rows is at the target on each side.
  lower <- min(nodes[["z"]])
  upper <- max(nodes[["z"]])
  mass  <- cbind(.5 * exp(-(c(.3, -.2) - lower) / c(.5, .25)),
                 .5 * exp(-(upper - c(.3, -.2)) / c(.5, .25)))
  expect_lte(max(mass), target)
  expect_gt(max(mass), target * exp(-.1 / .25))
})


test_that("an extension stops at values the joint density cannot evaluate", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  flat      <- function(values) cbind(rep(0, length(values)))
  lattice   <- .iwmde_qcmde_lattice(c(0, 1), 11L, transform)
  limited <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1L,
    target           = 1e-3,
    evaluate_inside  = flat,
    evaluate_outside = function(values) {
      out <- flat(values)
      out[values < -.25 | values > 1.35] <- NA_real_
      out
    }
  )
  expect_equal(range(limited[["nodes"]][["z"]]), c(-.2, 1.3),
               tolerance = 1e-12)
  expect_false(any(limited[["open"]]))

  # A row that never decays stops at the extension limit on each side.
  capped <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1L,
    target           = 1e-3,
    evaluate_inside  = flat,
    evaluate_outside = flat
  )
  limit <- .iwmde_qcmde_extension_limit()
  expect_equal(range(capped[["nodes"]][["z"]]), c(-limit, 1 + limit),
               tolerance = 1e-12)
  expect_identical(length(capped[["nodes"]][["z"]]), (2L * limit + 1L) * 10L + 1L)
})


test_that("a Gaussian tail extends close to its target in few passes", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  normal    <- function(values) cbind(stats::dnorm(values, log = TRUE))
  lattice   <- .iwmde_qcmde_lattice(c(-1, 1), 41L, transform)
  target    <- 5e-7
  result    <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1L,
    target           = target,
    evaluate_inside  = normal,
    evaluate_outside = normal
  )
  ends <- range(result[["nodes"]][["z"]])
  # The envelope overstates a Gaussian tail, so the true tails are below the
  # target; the local quadratic keeps the overshoot within a few steps.
  expect_lte(stats::pnorm(ends[[1L]]), target)
  expect_lte(stats::pnorm(ends[[2L]], lower.tail = FALSE), target)
  needed <- stats::qnorm(target, lower.tail = FALSE)
  expect_lt(max(abs(ends) - needed), .3)
  expect_lte(result[["passes"]], 4L)
})


.qcmde_flat_tail_rows <- function(half_width = 3, scale = .5, gaussian = TRUE) {

  # Row 1 is Gaussian. Row 2 is flat on (-half_width, half_width) with Laplace
  # tails outside, so it does not decay at the ends of a range inside the flat
  # part and has no predicted distance there.
  function(values) {
    flat <- -pmax(0, abs(values) - half_width) / scale
    if (gaussian) cbind(stats::dnorm(values, log = TRUE), flat) else cbind(flat)
  }
}


# The plain-R reference of the extension: the smallest symmetric extension of
# [-1, 1], in steps of .1, at which the tail estimate q(b) / |slope| / (trapezoid
# normalizer) of every row meets the target at both ends.
.qcmde_symmetric_extension_needed <- function(rows, target, max_steps = 100L) {

  worst <- vapply(0:max_steps, function(steps) {
    z           <- seq(-1 - .1 * steps, 1 + .1 * steps, by = .1)
    log_density <- rows(z)
    n           <- length(z)
    max(vapply(seq_len(ncol(log_density)), function(row) {
      y      <- exp(log_density[, row])
      area   <- sum(diff(z) * (y[-1L] + y[-n]) / 2)
      slopes <- c((log_density[2L, row] - log_density[1L, row]) / .1,
                  (log_density[n - 1L, row] - log_density[n, row]) / .1)
      max(ifelse(slopes > 0, c(y[1L], y[n]) / slopes / area, Inf))
    }, numeric(1)))
  }, numeric(1))

  min(which(worst <= target)) - 1L
}


test_that("a row that does not decay at the initial end extends the grid only as far as its tail needs", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  rows      <- .qcmde_flat_tail_rows()
  evaluated <- list()
  evaluate  <- function(values) {
    evaluated[[length(evaluated) + 1L]] <<- values
    rows(values)
  }
  lattice <- .iwmde_qcmde_lattice(c(-1, 1), 21L, transform)
  target  <- 1e-6
  needed  <- .qcmde_symmetric_extension_needed(rows, target)
  # The flat row needs the extension to leave its flat part: more than three
  # widths of the initial range, near the limit of four.
  expect_gt(needed, 60L)
  expect_lt(needed, 80L)

  result <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1:2,
    target           = target,
    evaluate_inside  = evaluate,
    evaluate_outside = evaluate
  )
  nodes <- result[["nodes"]]
  # Each side stops within one step of the smallest extension that meets the
  # target, and short of the extension limit that a jump by the width of the
  # grid would have overshot.
  expect_lte(abs(result[["steps"]][["lower"]] - needed), 1L)
  expect_lte(abs(result[["steps"]][["upper"]] - needed), 1L)
  expect_true(all(result[["open"]]))
  expect_equal(diff(range(nodes[["z"]])), .1 * (20L + sum(result[["steps"]])),
               tolerance = 1e-12)
  log_density <- nodes[["log_q"]] + nodes[["log_jacobian"]]
  tails <- .iwmde_qcmde_tail_estimates(
    nodes[["z"]], log_density,
    .iwmde_log_trapz_columns(nodes[["z"]], log_density)
  )
  expect_lte(max(tails[["lower"]], tails[["upper"]]), target)

  # No value is evaluated twice, and every kept node was evaluated; the values
  # of a slice beyond the step that meets the target are dropped.
  values  <- unlist(evaluated, use.names = FALSE)
  expect_length(unique(values), length(values))
  expect_true(all(nodes[["x"]] %in% values))
  expect_lt(length(values) - length(nodes[["x"]]), 21L)
})


test_that("a side whose rows do not decay grows its slices from an eighth of the range", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  rows      <- .qcmde_flat_tail_rows(gaussian = FALSE)
  evaluated <- list()
  evaluate  <- function(values) {
    evaluated[[length(evaluated) + 1L]] <<- values
    rows(values)
  }
  lattice <- .iwmde_qcmde_lattice(c(-1, 1), 21L, transform)
  target  <- 1e-6
  needed  <- .qcmde_symmetric_extension_needed(rows, target)

  result <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1L,
    target           = target,
    evaluate_inside  = evaluate,
    evaluate_outside = evaluate
  )
  expect_lte(abs(result[["steps"]][["lower"]] - needed), 1L)
  expect_lte(abs(result[["steps"]][["upper"]] - needed), 1L)

  # The calls of one side (every second call after the initial nodes): a slice
  # of ceiling(20 / 8) = 3 steps, then the steps already taken until the row
  # leaves its flat part (after 24 steps the end is at 3.4), and then the
  # predicted distance of the decaying tail, one step beyond the need at most.
  sizes <- vapply(evaluated, length, integer(1))
  expect_identical(sizes[[1L]], 21L)
  lower <- sizes[seq(2L, length(sizes), by = 2L)]
  expect_identical(lower[1:4], c(3L, 3L, 6L, 12L))
  expect_identical(length(lower), 5L)
  expect_gte(lower[[5L]], needed - 24L)
  expect_lte(lower[[5L]] - (needed - 24L), 2L)
})


test_that("a flat tail reaches the extension limit in a few calls", {

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  flat      <- function(values) cbind(rep(0, length(values)))
  calls     <- 0L
  lattice   <- .iwmde_qcmde_lattice(c(0, 1), 41L, transform)
  result <- .iwmde_qcmde_extend_lattice(
    lattice          = lattice,
    estimate_rows    = 1L,
    target           = 1e-3,
    evaluate_inside  = flat,
    evaluate_outside = function(values) {
      calls <<- calls + 1L
      flat(values)
    }
  )
  limit <- .iwmde_qcmde_extension_limit()
  expect_identical(unname(result[["steps"]]), rep(limit * 40L, 2L))
  # Slices of 5, 5, 10, 20, 40 and 80 steps reach the limit of 160 steps.
  expect_identical(calls, 2L * 6L)
})


test_that("the slice of an extension follows the predicted distance and the growth of an undecided side", {

  size <- function(distance, undecided, taken = 0L, n_nodes = 101L,
                   budget = 400L) {
    .iwmde_qcmde_extension_size(
      distance  = list(distance = distance, undecided = undecided),
      step      = .1,
      n_initial = 101L,
      taken     = taken,
      n_nodes   = n_nodes,
      budget    = budget
    )
  }
  # A side whose failing rows all decay takes their predicted distance.
  expect_identical(size(.72, FALSE), 8L)
  expect_identical(size(0, FALSE), 1L)
  # An undecided side takes an eighth of the range, doubling with the steps taken.
  expect_identical(size(0, TRUE), 13L)
  expect_identical(size(0, TRUE, taken = 13L), 13L)
  expect_identical(size(0, TRUE, taken = 52L), 52L)
  # Never beyond the width of the current grid or the remaining budget.
  expect_identical(size(30, TRUE), 100L)
  expect_identical(size(0, TRUE, taken = 300L, n_nodes = 401L), 300L)
  expect_identical(size(0, TRUE, taken = 300L, budget = 40L), 40L)
  expect_identical(size(Inf, FALSE, n_nodes = 21L), 20L)

  # The distance of a mixed side comes from its decaying rows alone.
  log_density <- cbind(-abs(seq(-1, 1, length.out = 5L)) / .5,
                       rep(0, 5L))
  distance <- .iwmde_qcmde_extension_distance(
    log_density = log_density,
    side        = "upper",
    step        = .5,
    excess      = c(log(1e3), Inf),
    slope       = c(2, 0)
  )
  expect_true(distance[["undecided"]])
  expect_equal(distance[["distance"]], log(1e3) / 2, tolerance = 1e-12)
  none <- .iwmde_qcmde_extension_distance(log_density[, 2L, drop = FALSE],
                                          "upper", .5, Inf, 0)
  expect_identical(none, list(distance = 0, undecided = TRUE))
})


test_that("the prefix of an extension is the shortest run that meets the target", {

  # A Laplace row with unit rate on [-2, 2] (the normalizer is about 2).
  z       <- seq(-2, 2, by = .5)
  laplace <- function(values) cbind(-abs(values))
  nodes   <- list(
    index        = seq_along(z) - 1L,
    x            = z,
    z            = z,
    log_jacobian = rep(0, length(z)),
    log_q        = laplace(z)
  )
  run <- function(values, index) {
    list(index = index, x = values, z = values,
         log_jacobian = rep(0, length(values)), log_q = laplace(values))
  }
  upper  <- run(2 + .5 * seq_len(12L), max(nodes[["index"]]) + seq_len(12L))
  lower  <- run(-2 - .5 * rev(seq_len(12L)),
                min(nodes[["index"]]) - rev(seq_len(12L)))
  log_normalizer <- .iwmde_log_trapz_columns(z, laplace(z))
  target <- 1e-3
  prefix <- .iwmde_qcmde_extension_prefix(nodes, upper, "upper", 1L,
                                          log_normalizer, target)
  expect_true(prefix[["met"]])

  # The plain-R reference: extend one node at a time and test the definition,
  # the density at the end over its slope and the trapezoid normalizer.
  reference <- function(n) {
    grid <- seq(-2, 2 + .5 * n, by = .5)
    y    <- exp(-abs(grid))
    area <- sum(diff(grid) * (y[-1L] + y[-length(y)]) / 2)
    y[[length(y)]] / 1 / area
  }
  first <- min(which(vapply(1:12, reference, numeric(1)) <= target))
  expect_identical(prefix[["n"]], first)
  expect_lt(prefix[["n"]], 12L)

  # A target no run reaches keeps the whole extension and reports it.
  none <- .iwmde_qcmde_extension_prefix(nodes, upper, "upper", 1L,
                                        log_normalizer, 1e-12)
  expect_false(none[["met"]])
  expect_identical(none[["n"]], 12L)

  # The lower side of the symmetric row reads its run from the end nearest the
  # grid and gives the same prefix.
  expect_identical(
    .iwmde_qcmde_extension_prefix(nodes, lower, "lower", 1L, log_normalizer,
                                  target),
    prefix
  )
  kept <- .iwmde_qcmde_subset_extension(lower, seq.int(12L - prefix[["n"]] + 1L,
                                                       12L))
  expect_identical(kept[["index"]], min(nodes[["index"]]) -
    rev(seq_len(prefix[["n"]])))
  expect_identical(nrow(kept[["log_q"]]), prefix[["n"]])
})


test_that("IWMDE keeps the draws-based range for its normalization-mass check", {

  # IWMDE integrates the marginal density to check its weights, so the range
  # stays the central `normalization_prob` quantile range of the draws, widened
  # by 10% on each side; the per-row qCMDE range does not apply to it.
  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  values    <- stats::qnorm(stats::ppoints(400L), .2, .5)
  grid <- .iwmde_normalization_grid(values, numeric(), c(-Inf, Inf), transform,
                                    normalization_points = 50L,
                                    normalization_prob   = .999)
  quantiles <- stats::quantile(values, c(5e-4, 1 - 5e-4), names = FALSE,
                               type = 8)
  expect_equal(range(grid[["x"]]), quantiles + c(-.1, .1) * diff(quantiles),
               tolerance = 1e-12)

  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {
      matrix(stats::dnorm(values, .2, .5, log = TRUE), nrow = length(values),
             ncol = length(row_states))
    },
    .package = "RoBMA"
  )
  rows   <- seq_len(40L)
  states <- lapply(values[rows], function(value) {
    list(baseline_log_q = stats::dnorm(value, .2, .5, log = TRUE))
  })
  iwmde <- .iwmde_density_iwmde(
    context            = list(),
    parameter          = "mu",
    display_grid       = 0,
    row_states         = states,
    active_rows        = rows,
    active_values      = values[rows],
    proposal_weight    = list(
      log_weight = stats::dnorm(values[rows], .2, .5, log = TRUE),
      method     = "oracle_posterior"
    ),
    active_mass        = 1,
    replacement        = list(type = "scalar"),
    normalization_grid = grid
  )
  expect_equal(iwmde[["normalization_range"]], range(grid[["x"]]))
  expect_identical(iwmde[["normalization_points"]], 50L)
  expect_null(iwmde[["normalization_truncation"]])
  expect_equal(iwmde[["support_grid_normalization_integral"]],
               diff(stats::pnorm(range(grid[["x"]]), .2, .5)),
               tolerance = 1e-3)
})
