context("Selection certificate counterexamples")
skip_on_cran()

# Blocks of coincident rows whose accepted deterministic normalizer was wrong
# under the geometric transition rule, or that the acceptance rule must keep
# right. Each reference is computed here from independent integration
# (composite Gauss-Legendre or nested stats::integrate), never from the kernel.
# A block either carries a deterministic result (relative MCSE exactly zero)
# within the relative tolerance of the reference, or reaches the randomized
# fallback; the expected route is stated for every block. The shift-amplitude
# certificate defers the blocks that fed95b63 accepted wrongly; the refinement
# tier (the average over four shifts per axis at the finest affordable rung)
# certifies them, within 1e-5 of the reference, except the forest of 371 rows,
# whose rows narrow the tier gate beyond its edges and which therefore reaches
# the fallback.

.test_cert_tolerance <- .005

.test_cert_check_case <- function(name, case, native) {

  reference <- case[["reference"]](case[["block"]])
  if (!is.null(case[["pinned"]])) {
    # a changed reference must be visible
    expect_lt(abs(reference - case[["pinned"]]), 1e-8, label = paste(name, "reference"))
  }
  actual <- native(case[["block"]])
  route  <- if (.test_cert_deterministic(actual)) "deterministic" else "fallback"
  if (route == "deterministic") {
    expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference),
              .test_cert_tolerance, label = paste(name, "relative error"))
  }
  if (route == "fallback") {
    # a randomized result is held to the same tolerance (the fixed design makes
    # it reproducible)
    expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], reference),
              .test_cert_tolerance, label = paste(name, "relative error of the fallback"))
  }
  expect_identical(route, case[["route"]], label = paste(name, "route"))
  invisible(actual)
}

.test_cert_parallel_block <- function(loading, sd, mean, count, omega, z_lower = c(0, -Inf),
                                      sei = 1) {

  .test_cert_block(list(.test_cert_parallel_group(count, loading, sd, mean, sei)), omega, z_lower)
}

test_that("parallel coincident rows on two axes keep their normalizer within the tolerance", {

  golden <- (1 + sqrt(5)) / 2
  scale  <- sqrt(1 + golden^2)
  cases  <- list(
    # Sol, review 33: eight rows with the selection deep in a tail
    sol8_tail = list(pinned = -12.495579111516, route = "deterministic",
      block = .test_cert_parallel_block(c(.224, .975), .344, -4.08, 8L, c(1, 0))),
    # an 80-row block with weights (1, .2)
    sol80_w2 = list(pinned = -6.013912874415, route = "deterministic",
      block = .test_cert_parallel_block(c(.42, .91), .53, -1.65, 80L, c(1, .2))),
    # the documented 80-row case that the geometric rule accepted with error 1.04e-2
    fed80 = list(pinned = -9.012420977204, route = "deterministic",
      block = .test_cert_parallel_block(c(.496441, .542314), .246872, -2.125758, 80L, c(1, .01))),
    # Opus golden-ratio blocks, always routed to the fallback by the geometric rule
    golden20 = list(pinned = -11.228950361197, route = "deterministic",
      block = .test_cert_parallel_block(c(1, golden), .11136 * scale, -4 * scale, 20L, c(1, 0))),
    golden40 = list(pinned = -11.224743078254, route = "deterministic",
      block = .test_cert_parallel_block(c(1, golden), .11136 * scale, -4 * scale, 40L, c(1, .5))),
    # a block that an aliasing-free rule resolves early
    qmc10 = list(pinned = -6.920829350705, route = "deterministic",
      block = .test_cert_parallel_block(c(sqrt(.5), sqrt(.2)), .3, -4 * sqrt(.7), 10L, c(1, .5)))
  )
  for (name in names(cases)) {
    cases[[name]][["reference"]] <- .test_cert_reference_parallel
    .test_cert_check_case(name, cases[[name]], .test_cert_native_factor)
  }
})

test_that("a rounded mean of the documented 80-row case keeps the contract", {

  # The geometric rule's contract said that rounding the mean to -2.13 reaches
  # the fallback; the refinement tier now certifies this block as well.
  case <- list(route = "deterministic", reference = .test_cert_reference_parallel,
    block = .test_cert_parallel_block(c(.496441, .542314), .246872, -2.13, 80L, c(1, .01)))
  .test_cert_check_case("fed80_rounded", case, .test_cert_native_factor)
})

test_that("independent coincident groups on separate axes keep their normalizer within the tolerance", {

  band   <- .test_cert_band_weights(.97)
  cuts   <- .test_cert_band_cuts()
  group  <- function(count, loading, sd, mean, sei, columns, rank) {
    full <- numeric(rank)
    full[columns] <- loading
    .test_cert_parallel_group(count, full, sd, mean, sei)
  }

  # Sol, review 35: two independent two-axis groups (327 rows, four axes)
  reduced <- .test_cert_block(list(
    group(114L, c(.9732175711714035, .2298859699138557), .6328158871331484,
          -1.469713798724115, .463680501201846, 1:2, 4L),
    group(213L, c(.2063287166395255, .9784827339763775), .7306215995525465,
          .9045836022123694, .8706529675272152, 3:4, 4L)), band, cuts)
  .test_cert_check_case("sol35_reduced", list(block = reduced, route = "deterministic",
    pinned = .738439020239 + 3.776697872506, reference = .test_cert_reference_parallel),
    .test_cert_native_factor)

  # four independent two-axis groups (594 rows, eight axes)
  table <- data.frame(count = c(235L, 44L, 124L, 191L), a1 = c(.826, .733, .577, .035),
    a2 = c(.564, .680, .817, .999), sd = c(.558, .557, .679, .660),
    center = c(-1.456, 1.887, -.399, -.354), sei = c(.310, .061, .127, 1.078))
  eight <- .test_cert_block(lapply(seq_len(4L), function(i) {
    row <- table[i, ]
    group(row$count, c(row$a1, row$a2), row$sd, -row$center, row$sei, 2L * i - 1:0, 8L)
  }), band, cuts)
  .test_cert_check_case("sol35_eight_axes", list(block = eight, route = "deterministic",
    pinned = 1.320954983677 + .015317152283 + .322823247972 + 4.353098088774,
    reference = .test_cert_reference_parallel), .test_cert_native_factor)
})

test_that("two crossing coincident groups keep their normalizer within the tolerance", {

  # An adversarial search found 92 + 32 rows on two crossing directions where the
  # geometric rule accepted order 63 with relative error 7.66e-3.
  block <- .test_cert_block(list(
    list(count = 92L, loading = c(1.59892823963330177, .46248906601712753),
         sd = 1.17935060431781880, mean = -.12799436126613484, sei = 1),
    list(count = 32L, loading = c(-.37056447402550319, .48747887335740198),
         sd = .23232799104534663, mean = -2.32439955742265258, sei = 1)),
    c(1, .03633489870625161))
  .test_cert_check_case("crossing124", list(block = block, route = "deterministic",
    pinned = -21.277252217891, reference = .test_cert_reference_two_groups),
    .test_cert_native_factor)
})

test_that("a depth-one forest with cancelling leaf aliases keeps its normalizer within the tolerance", {

  # Astra's 371-row forest: a shared axis (loading .001) and two leaf axes of
  # opposite orientation. Leaf amplitudes that are shifted together nearly
  # cancel (0.0006) and accepted order 31 with relative error -1.02e-2; the
  # amplitudes of separate leaf shifts sum to 0.034.
  block <- .test_cert_block(list(
    list(count = 115L, loading = c(.001, 1, 0), sd = .249373398318475, mean = 1.93906024750322, sei = 1),
    list(count = 256L, loading = c(.001, 0, -1), sd = .103488784220421, mean = 1.91516049392521, sei = 1)),
    c(1, .01))
  # The refinement tier's gate narrows with the 371 rows on the shared axis and
  # refuses the rung (see the refinement tests): the block reaches the fallback.
  .test_cert_check_case("forest371", list(block = block, route = "fallback",
    pinned = -.15663896170926442, reference = .test_cert_reference_forest),
    .test_cert_native_factor)
})

test_that("coincident rank-one clusters keep their normalizer within the tolerance", {

  # Rank-one clusters from the research battery where the geometric rule accepted
  # relative errors of 9.7e-3, 2.0e-2 and 1.8e-2 (orders 31, 127 and 63).
  cluster <- function(count, loading, sd, sei, t0, cut, w2) {
    .test_cert_block(list(.test_cert_parallel_group(count, loading, sd, cut * sei - loading * t0, sei)),
      c(1, w2), c(cut, -Inf))
  }
  cases <- list(
    cluster180 = list(route = "deterministic",
      block = cluster(180L, .49190045987952141, .34623272304241315, 1.1742373526276058,
                      -.93498589762020856, 1.645, .01)),
    cluster300 = list(route = "deterministic",
      block = cluster(300L, 1.3060798299269314, .51812101191729221, .55784746236513971,
                      3.0003575886366889, 1.645, .1)),
    cluster85 = list(route = "deterministic",
      block = cluster(85L, 2.1997086546647093, 1.0857977532550365, .44680967672203892,
                      2.1270223344909027, 0, .05))
  )
  for (name in names(cases)) {
    cases[[name]][["reference"]] <- .test_cert_reference_parallel
    .test_cert_check_case(name, cases[[name]], .test_cert_native_cluster)
  }
})

test_that("a rung that coincides with the aliasing of the tilted density is detected by its amplitude", {

  # Sol 8 rows: the order-31 rung estimates the normalizer 7.9e-3 low, and the
  # neighbouring rungs 15, 21 and 31 sit near zero crossings of their errors, so
  # their successive change is small. The shift-amplitude certificate of the
  # rung itself is not: its summed amplitude covers the actual error and exceeds
  # half the tolerance.
  block     <- .test_cert_parallel_block(c(.224, .975), .344, -4.08, 8L, c(1, 0))
  reference <- .test_cert_reference_parallel(block)
  for (order in c(21L, 31L)) {
    certificate <- .test_cert_amplitudes(block, order)
    error       <- abs(expm1(certificate[["log_integral"]] - reference))
    expect_lt(error, sum(certificate[["amplitude"]]))
    expect_gt(sum(certificate[["amplitude"]]), .test_cert_tolerance / 2)
  }
})

test_that("captured Assink factor states are accepted at the order the certificate recorded", {

  # Cost guard. Six states of the three factor blocks of the Assink CCI fit
  # (22, 5 and 16 rows; ranks 3, 2 and 4) with the order that the shift-amplitude
  # certificate accepted (research record assink-cert.rds) and the plain-rule
  # normalizers at every order. A deterministic acceptance at the recorded order
  # returns that rung's value; a higher order would show as another value.
  states <- readRDS(testthat::test_path("fixtures", "selection-certificate-assink-states.rds"))
  rules  <- .selection_joint_factor_quadrature_rules()
  for (state in states) {
    K      <- state[["K"]]
    qmc    <- BayesTools::selection_qmc_design(
      dimensions = 2L * state[["rank"]], points = 1024L, scrambles = 8L, seed = 1L)
    actual <- .Call("RoBMA_selnorm_factor_step_loglik_batch", rep(0, K), matrix(state[["mean"]], 1L),
      matrix(state[["residual_sd"]], 1L), matrix(state[["loading"]], 1L), state[["sei"]],
      matrix(state[["omega"]], 1L), state[["z_lower"]], state[["z_upper"]], state[["obs_bin"]],
      state[["sign"]], TRUE, 1L, rules[["nodes"]], rules[["log_weights"]], as.double(rules[["orders"]]),
      as.double(qmc), 256L, 1024L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
    info <- sprintf("block %d state %d", state[["block"]], state[["state"]])
    expect_true(.test_cert_deterministic(actual), info = info)
    accepted <- which.min(abs(actual[["log_normalizer"]] - state[["recorded_values"]]))
    expect_identical(as.integer(rules[["orders"]][accepted]), as.integer(state[["recorded_order"]]),
                     info = info)
    expect_equal(actual[["log_normalizer"]], state[["recorded_values"]][accepted],
                 tolerance = 1e-9, info = info)
    expect_lt(.test_cert_relative_error(actual[["log_normalizer"]], state[["reference"]]),
              .test_cert_tolerance, label = info)
  }
})
