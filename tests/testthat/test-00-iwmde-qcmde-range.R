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
