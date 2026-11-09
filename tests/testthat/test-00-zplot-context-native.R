# The native batch (src/selnorm/zplot-context.cc) ports the R controller
# .zplot_context_projection() statement by statement; these tests pin the two
# to identical results, errors and refusals, and the batch to thread-count
# invariance and ordered replay inside the joint zplot loop.

.zplot_native_test_selection <- function(sei, cutoffs, weights, direction) {

  prior <- BayesTools::prior_weightfunction("one-sided",
    stats::pnorm(cutoffs, lower.tail = FALSE), BayesTools::wf_fixed(weights),
    model = BayesTools::selection_model(weight_rule = "product", group = "paper"))
  context <- .selection_spec(list(outcome = list(bias = prior)), rep(0, length(sei)), sei, direction)
  context$omega <- matrix(context$fixed_omega[1L, ], 1L)
  context$alpha <- 0
  context$phack_kind <- 0L
  context$kernel_mode <- SELKERNEL_STEP
  context$vector_rule <- 0L
  context$use_normal <- FALSE
  BayesTools::selection_context_validate(context, n_samples = 1L)
}


# One fixture per status and per reachable refinement step of the ladder,
# found by a randomized search over the R controller and by review (each
# comment names the steps the native trace records; ZplotTraceCode in
# zplot-context.h). Pruning routing is not reachable: its allocation keeps the
# pruning error below the routing threshold.
.zplot_native_test_fixtures <- list(
  # ok-point-absolute: route compression (eight halvings, the last to zero)
  list(
    cutoffs = c(2.7073688422096893, 2.5156306988792494, 1.522085587028414),
    weights = c(1, 0.67502524894375393, 0.20437466374088981, 0.063700387783312107),
    direction = "positive",
    sei = c(0.243497576341975, 0.20383368450977943),
    covariance = c(0.0077905055288595287, 0.00099697707429221107, 0.00099697707429221107,
      0.0022804408084654212),
    factor = c(0, 0),
    mean = c(-0.016321328318342997, -0.15736476617395681),
    tolerance = 0.0095661360668959864,
    grid = 41,
    context_means = c(0.051629034384735088, 0.12530872485333505, -0.03610120369623919, 0.27549200724393885,
      -0.12994916115724406, 0.097788129583496752),
    absolute = 1e-07
  ),
  # ok-continuous: omission allowance halving
  list(
    cutoffs = 1.96,
    weights = c(1, 1e-04),
    direction = "positive",
    sei = c(1, 1),
    covariance = c(1, 0.1, 0.1, 1),
    factor = c(0, 0),
    mean = c(-1, -1),
    tolerance = 0.1,
    grid = 41,
    z = c(-2, 0),
    context_means = NULL,
    absolute = NULL
  ),
  # ok-point-absolute
  list(
    cutoffs = c(2.7073688422096893, 2.5156306988792494, 1.522085587028414),
    weights = c(1, 0.67502524894375393, 0.20437466374088981, 0.063700387783312107),
    direction = "positive",
    sei = c(0.243497576341975, 0.20383368450977943),
    covariance = c(0.0077905055288595287, 0.00099697707429221107, 0.00099697707429221107,
      0.0022804408084654212),
    factor = c(0, 0),
    mean = c(-0.016321328318342997, -0.15736476617395681),
    tolerance = 0.0095661360668959864,
    grid = 41,
    context_means = c(0.051629034384735088, 0.12530872485333505, -0.03610120369623919, 0.27549200724393885,
      -0.12994916115724406, 0.097788129583496752),
    absolute = 0.00019843652922397337
  ),
  # ok-uniform
  list(
    cutoffs = 2.0834214928327128,
    weights = c(1, 1),
    direction = "positive",
    sei = c(1.1057117523169455, 0.70770654483608808),
    covariance = c(1.2340146109475532, 0.065647298183371158, 0.065647298183371158, 0.33289824178764987
      ),
    factor = c(0, 0),
    mean = c(2.0697022027987568, 0.62378863733408185),
    tolerance = 0.033247505800204587,
    grid = 81,
    context_means = NULL,
    absolute = NULL
  ),
  # ok-continuous: route inner
  list(
    cutoffs = c(1.7858512910315767, 0.98434114316478372),
    weights = c(1, 0.68249933431040666, 0.043856536904050517),
    direction = "negative",
    sei = c(0.17424745146589066, 1.8249270203470578),
    covariance = c(0.029763287975451309, 0.086594644260001075, 0.086594644260001075, 2.7021663776077696
      ),
    factor = c(0, 0),
    mean = c(0.081593784466890251, 2.4027718880828632),
    tolerance = 0.00028108202925614182,
    grid = 81,
    context_means = NULL,
    absolute = NULL
  ),
  # null
  list(
    cutoffs = c(1.7458741059526801, 1.7367047689622268),
    weights = c(1, 0.59278764351732893, 0.0086651965236613279),
    direction = "positive",
    sei = c(0.43381405421000124, 0.77436830156400505),
    covariance = c(-0.11116340237786237, 0.19054646592531996, 0.19054646592531996, 0.6774981717006302
      ),
    factor = c(0.46607838743547519, 0.48283273056300652),
    mean = c(0.36206224664075565, 0.043265982782622572),
    tolerance = 0.0013020227485876737,
    grid = 81,
    context_means = NULL,
    absolute = NULL
  ),
  # error
  list(
    cutoffs = 1.3609706454444677,
    weights = c(1, 0.0010302596544285272),
    direction = "positive",
    sei = c(0.18622174335067482, 0.28365108104507425),
    covariance = c(0.0066640418462526208, -0.0077450423408460464, -0.0077450423408460464,
      0.033190391823351366),
    factor = c(0.099897735945640789, 0.17927616046226097),
    mean = c(NA, 0.23579064813493272),
    tolerance = 2.3099341553207714e-05,
    grid = 81,
    context_means = NULL,
    absolute = NULL
  ),
  # ok-point: resolution inner, own normalizer, route normalizer
  list(
    cutoffs = c(2.767848196788691, 2.1952514951117337, 1.9072400843724608),
    weights = c(1, 0.03736417946859602, 8.5670031909373718e-06, 2.8043677137402575e-06
      ),
    direction = "positive",
    sei = c(0.30353715291723782, 1.1221395716925358),
    covariance = c(0.006760051851716146, 0.0092567268996465463, 0.0092567268996465463, 0.059641573622465929
      ),
    factor = c(0, 0),
    mean = c(0.28553234916624914, 4.0025590362726744),
    tolerance = 5.416883940969056e-05,
    grid = 81,
    context_means = c(0.024988078364047359, 0.70547472465551442, -0.14444503797689329, 0.96991950246216962,
      0.23171046697891939, 0.084058979793193392),
    absolute = NULL
  ),
  # ok-continuous: resolution context, route outer
  list(
    cutoffs = 1.7634531810181215,
    weights = c(1, 0.011958018657704102),
    direction = "positive",
    sei = c(1.1256546202638633, 0.14461667434327746),
    covariance = c(0.93634714473871217, 0.0038728741638097961, 0.0038728741638097961, 0.017915956249708081
      ),
    factor = c(1.0820340379343232, 0.038403492716783463),
    mean = c(0.080733078058379107, 0.057129735389922873),
    tolerance = 0.0040395381602232987,
    grid = 41,
    context_means = NULL,
    absolute = NULL
  ),
  # ok-continuous: route cdf
  list(
    cutoffs = 1.3269231696613133,
    weights = c(1, 5.4976276847508648e-06),
    direction = "positive",
    sei = c(1.7224406777897796, 0.15915531155258916),
    covariance = c(2.6021172377876622, 0.050502757155361608, 0.050502757155361608, 0.024502915177519048
      ),
    factor = c(0, 0),
    mean = c(-0.90253675021894064, -0.019451908382400608),
    tolerance = 0.00017698030938607088,
    grid = 81,
    context_means = NULL,
    absolute = NULL
  ),
  # null: resolution inner, selection inner, selection denominator, own normalizer, route covariance
  list(
    cutoffs = c(2.8622508143307641, 2.6100954705616459, 2.5472676382632926),
    weights = c(1, 3.8390178849804782e-07, 1.1931289848081423e-07, 3.9114827286973549e-08
      ),
    direction = "positive",
    sei = c(0.84717377405758421, 1.5975466611053655),
    covariance = c(0.42514980011793274, 0.033100102932744686, 0.033100102932744686, 0.12209572971013377
      ),
    factor = c(0, 0),
    mean = c(0.25277433419942802, -3.1231142723832113),
    tolerance = 3.9321708131807455e-07,
    grid = 81,
    context_means = c(-0.31554292254361277, 1.30392164112145, -0.18446141332253058, 1.2059609292277993,
      0.58964508276721128, 1.8122613720931509, 0.025315174250746416, -0.18546422533462575
      ),
    absolute = NULL
  ),
  # ok-continuous: resolution context, selection outer
  list(
    cutoffs = c(1.6178404465317726, 1.2139540510252118),
    weights = c(1, 0.52291121609652547, 2.1025261148339732e-05),
    direction = "positive",
    sei = c(0.21956220700762635, 0.33591917812722555),
    covariance = c(0.03372155391884734, 0.016594087077130542, 0.016594087077130542, 0.094943118493394887
      ),
    factor = c(0.32432360101371954, 0.14428380498761176),
    mean = c(0.43859677923243873, 0.36390721333935422),
    tolerance = 0.0093815171664404298,
    grid = 81,
    context_means = NULL,
    absolute = NULL
  )
)


.zplot_native_test_case <- function(fixture) {

  k <- length(fixture$sei)
  context_means <- if (is.null(fixture$context_means)) NULL else matrix(fixture$context_means, ncol = k)
  range <- if (is.null(fixture$z)) c(-4, 6) else fixture$z
  list(z = seq(range[[1L]], range[[2L]], length.out = fixture$grid),
    mean = if (is.null(context_means)) fixture$mean else colMeans(context_means),
    covariance = matrix(fixture$covariance, k), factor = matrix(fixture$factor, 1L), sei = fixture$sei,
    selection = .zplot_native_test_selection(fixture$sei, fixture$cutoffs, fixture$weights, fixture$direction),
    control = set_selection_likelihood_control(relative_tolerance = fixture$tolerance),
    context_means = context_means, absolute = fixture$absolute)
}


.zplot_native_test_reference <- function(case) {

  tryCatch(.zplot_context_projection(case$z, case$mean, case$covariance, case$factor, case$sei,
    case$selection, case$control, case$absolute, context_means = case$context_means),
    error = function(e) e)
}


.zplot_native_test_batch <- function(case) {

  k <- length(case$sei)
  .zplot_context_projection_batch(case$z, matrix(case$mean, 1L), array(case$covariance, c(k, k, 1L)),
    array(case$factor, c(1L, k, 1L)), case$sei, case$selection, case$selection$omega, case$control,
    case$absolute, context_means = case$context_means)
}


.zplot_native_test_same <- function(reference, batch, draw) {

  native <- tryCatch(.zplot_context_projection_result(batch, draw), error = function(e) e)
  if (inherits(reference, "error")) {
    expect_s3_class(native, "error")
    expect_identical(conditionMessage(native), conditionMessage(reference))
  } else if (.test_fused_multiply_add() && !is.null(reference)) {
    # The R controller rounds a * b + c twice and the native batch once. Measured
    # gaps on Apple silicon: 1.8e-15 relative on the density, 3.1e-11 absolute on
    # the error diagnostics, whose small sums cancel; everything else stays exact.
    errors <- grep("_error$", names(reference), value = TRUE)
    exact  <- setdiff(names(reference), c("density", errors))
    expect_identical(names(native), names(reference))
    expect_identical(native[exact], reference[exact])
    gap <- abs(native[["density"]] - reference[["density"]])
    scale <- pmax(abs(native[["density"]]), abs(reference[["density"]]))
    expect_lte(max(gap[gap > 0] / scale[gap > 0], 0), 16 * .Machine$double.eps)
    expect_lte(max(abs(unlist(native[errors]) - unlist(reference[errors]))), 1e-10)
  } else {
    expect_identical(native, reference)
  }
}


test_that("the native context projection reproduces the R controller on every ladder branch", {

  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)
  RoBMA.options(native_threads = 1L)
  statuses <- character()
  steps <- integer()
  for (fixture in .zplot_native_test_fixtures) {
    case <- .zplot_native_test_case(fixture)
    reference <- .zplot_native_test_reference(case)
    batch <- .zplot_native_test_batch(case)
    .zplot_native_test_same(reference, batch, 1L)
    statuses <- c(statuses, if (inherits(reference, "error")) "error" else if (is.null(reference)) "null" else
      if (is.null(reference$quadrature)) "uniform" else if (!is.null(case$context_means)) "point" else "continuous")
    steps <- c(steps, batch$trace[[1L]])
  }
  expect_setequal(statuses, c("error", "null", "uniform", "point", "continuous"))
  expect_true(all(c(1:11, 13L, 14L) %in% steps))
  expect_identical(sum(steps == 11L), 8L)
})


test_that("the native context projection keeps the controller's argument checks", {

  case <- .zplot_native_test_case(.zplot_native_test_fixtures[[5L]])
  k <- length(case$sei)
  run <- function(case) {
    list(reference = .zplot_native_test_reference(case), batch = .zplot_native_test_batch(case))
  }
  # A missing weight is an error; a nonpositive one decides first.
  for (omega in list(c(1, NA, .2), c(0, NA, .2))) {
    modified <- case
    modified$selection$omega <- matrix(omega[seq_len(ncol(case$selection$omega))], 1L)
    out <- run(modified)
    .zplot_native_test_same(out$reference, out$batch, 1L)
  }
  # The tolerance is checked where the controller checks it.
  modified <- case
  modified$control$relative_tolerance <- 0
  out <- run(modified)
  .zplot_native_test_same(out$reference, out$batch, 1L)
  uniform <- .zplot_native_test_case(.zplot_native_test_fixtures[[4L]])
  uniform$control$relative_tolerance <- 0
  out <- run(uniform)
  .zplot_native_test_same(out$reference, out$batch, 1L)
  expect_error(.zplot_context_projection_batch(case$z, matrix(case$mean, 1L),
    array(case$covariance, c(k, k, 1L)), array(case$factor, c(1L, k, 1L)), case$sei, case$selection,
    case$selection$omega, case$control, absolute_tolerance = -1),
    "absolute tolerance must be finite and positive")
})


test_that("the batched context projection is invariant to the native thread count", {

  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)
  case <- .zplot_native_test_case(.zplot_native_test_fixtures[[12L]])
  k <- length(case$sei)
  withr::local_seed(20261105L)
  S <- 40L
  means <- matrix(case$mean, S, k, byrow = TRUE) + matrix(stats::rnorm(S * k, 0, .1), S, k)
  covariance <- array(case$covariance, c(k, k, S)) * rep(stats::runif(S, .8, 1.2), each = k * k)
  factor <- array(case$factor, c(1L, k, S)) * rep(stats::runif(S, .5, 1.5), each = k)
  factor[, , 7L] <- 0
  omega <- matrix(case$selection$omega[1L, ], S, ncol(case$selection$omega), byrow = TRUE)
  omega[, -1L] <- omega[, -1L] * stats::runif(S, .6, 1.4)
  omega[5L, ] <- omega[5L, 1L]
  omega[9L, 2L] <- 0
  means[12L, 1L] <- Inf
  context <- case$selection
  context$omega <- omega[1L, , drop = FALSE]
  batches <- lapply(c(1L, 2L, 3L, 8L), function(threads) {
    RoBMA.options(native_threads = threads)
    .zplot_context_projection_batch(case$z, means, covariance, factor, case$sei, context, omega, case$control)
  })
  for (batch in batches[-1L]) expect_identical(batch, batches[[1L]])
  expect_true(all(c(0L, 1L, 2L) %in% batches[[1L]]$status))
  expect_true(batches[[1L]]$uniform[[5L]])
  for (draw in seq_len(S)) {
    context$omega <- omega[draw, , drop = FALSE]
    reference <- tryCatch(.zplot_context_projection(case$z, means[draw, ], covariance[, , draw],
      matrix(factor[, , draw], 1L), case$sei, context, case$control), error = function(e) e)
    .zplot_native_test_same(reference, batches[[2L]], draw)
  }
})


test_that("the native route follows the default matrix product only", {

  withr::local_options(matprod = "default")
  expect_true(.zplot_context_native_available())
  withr::local_options(matprod = "internal")
  expect_false(.zplot_context_native_available())
})


test_that("parallel uniform-weight draws stay interruptible", {

  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)
  RoBMA.options(native_threads = 2L)
  k <- 40L
  S <- 4L
  sei <- rep(1, k)
  selection <- .zplot_native_test_selection(sei, 1.96, c(1, 1), "positive")
  omega <- matrix(selection$omega[1L, ], S, ncol(selection$omega), byrow = TRUE)
  z <- seq(-4, 4, length.out = 400000L)
  means <- matrix(.1, S, k)
  covariance <- array(diag(k), c(k, k, S))
  factor <- array(0, c(1L, k, S))
  control <- set_selection_likelihood_control()
  run <- function() .zplot_context_projection_batch(z, means, covariance, factor, sei, selection, omega, control)
  # The uncancelled batch, which also warms every code path the limited run takes.
  full <- system.time(run())[["elapsed"]]
  skip_if(full < .5, "The uniform-weight workload is too small to time an interrupt.")
  # The grid work polls on the main thread, so an elapsed-time limit stops the
  # batch well before the remaining draws would have finished.
  setTimeLimit(elapsed = .05, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  started <- proc.time()[["elapsed"]]
  expect_error(run(), "time limit")
  interrupted <- proc.time()[["elapsed"]] - started
  setTimeLimit(elapsed = Inf)
  expect_lt(interrupted, full / 2)
})
