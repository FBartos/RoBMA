test_that("exact selection caching preserves targets under controls and memory pressure", {

  old_limit <- RoBMA.get_option("selection.cache_max_bytes")
  on.exit({
    RoBMA.options(selection.cache_max_bytes = old_limit)
    selection_cache_info(clear = TRUE)
  }, add = TRUE)
  info <- function(clear = FALSE) {

    result <- selection_cache_info(clear = clear)
    result[result[["component"]] == "total", , drop = FALSE]
  }
  rho <- .4
  weight <- .2
  sei <- c(.2, .3, .4)
  covariance <- rho * outer(sei, sei)
  diag(covariance) <- sei^2
  quadrature <- .selection_joint_cluster_quadrature_rules(c(7L, 15L, 31L, 63L))
  arguments <- list(rep(0, 3L), matrix(rep(0, 3L), 1L),
    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L), sei,
    matrix(c(1, weight), 1L), c(0, -Inf), c(Inf, 0), rep(1L, 3L),
    1L, TRUE, 1L,
    as.double(BayesTools::selection_qmc_design(6L, 512L, 8L, seed = 173L)),
    512L, 8L, .005, TRUE, 0L, quadrature)
  evaluate <- function(values) {

    do.call(.Call, c(list("RoBMA_selnorm_mnorm_step_loglik_batch"), values,
      list(PACKAGE = "RoBMA")))
  }
  diagnostic_error <- function(result) {

    diagnostic <- result[["integration_diagnostics"]][1L, ]
    diagnostic[["covariance_width"]] + 2 * diagnostic[["quadrature_change"]] +
      diagnostic[["tail_bound"]]
  }
  # Expansion of product weights into centered Gaussian orthant probabilities:
  # P_i=1/2, P_ij=1/4+asin(rho)/(2*pi), P_123=1/8+3*asin(rho)/(4*pi).
  pair <- .25 + asin(rho) / (2 * pi)
  triple <- .125 + 3 * asin(rho) / (4 * pi)
  reference <- log(weight^3 + 3 * weight^2 * (1 - weight) / 2 +
    3 * weight * (1 - weight)^2 * pair + (1 - weight)^3 * triple)

  changed_mean <- changed_omega <- changed_covariance <- changed_controls <- tightened <- arguments
  changed_mean[[2L]][1L, ] <- c(.1, -.1, .2)
  changed_omega[[5L]][1L, 2L] <- .6
  changed_covariance[[3L]][1L, c(1L, 4L, 6L)] <-
    changed_covariance[[3L]][1L, c(1L, 4L, 6L)] + .01
  changed_controls[[18L]] <- .selection_joint_cluster_quadrature_rules(c(15L, 31L, 63L))
  tightened[[15L]] <- 1e-8
  cases <- list(arguments, changed_mean, changed_omega, changed_covariance,
    changed_controls, arguments, tightened)
  RoBMA.options(selection.cache_max_bytes = 0)
  uncached <- lapply(cases, evaluate)
  expect_lt(abs(uncached[[1L]][["log_normalizer"]] - reference), 5e-4)
  numerator <- mvtnorm::dmvnorm(rep(0, 3L), sigma = covariance, log = TRUE)
  expect_equal(uncached[[1L]][["log_density"]] + uncached[[1L]][["log_normalizer"]],
    numerator, tolerance = 1e-12)
  expect_equal(info()[["allocated_bytes"]], 0)

  # A large permitted cap is not a request to allocate that capacity.
  RoBMA.options(selection.cache_max_bytes = 2^31)
  expect_equal(info()[["capacity_bytes"]], 2^31)
  expect_equal(info()[["allocated_bytes"]], 0)
  first <- evaluate(arguments)
  expect_identical(first, uncached[[1L]])
  expect_gt(info()[["allocated_bytes"]], 0)
  expect_lt(info()[["allocated_bytes"]], 2^31)
  expect_identical(evaluate(arguments), first)
  expect_gt(info()[["hits"]], 0)

  # Copies may be changed and their original values revisited. Rule/tolerance
  # changes must not restore a previously cached accuracy result.
  cached <- lapply(cases, evaluate)
  expect_identical(cached, uncached)
  expect_lte(diagnostic_error(cached[[7L]]), 1e-8)
  RoBMA.options(selection.cache_max_bytes = 1)
  expect_identical(evaluate(arguments), uncached[[1L]])
  expect_equal(info()[["allocated_bytes"]], 0)

  RoBMA.options(selection.cache_max_bytes = 64 * 1024)
  pressure_arguments <- arguments
  for (index in seq_len(400L)) {
    pressure_arguments[[2L]][] <- index * 1e-4
    last <- evaluate(pressure_arguments)
  }
  pressure <- info()
  expect_true(is.finite(last[["log_normalizer"]]))
  expect_lte(diagnostic_error(last), .005)
  expect_gt(pressure[["evictions"]], 0)
  expect_lte(pressure[["allocated_bytes"]], pressure[["capacity_bytes"]])
  expect_lte(pressure[["peak_bytes"]], pressure[["capacity_bytes"]])
  expect_identical(evaluate(arguments), uncached[[1L]])
  cleared <- info(clear = TRUE)
  expect_equal(cleared[["allocated_bytes"]], 0)
  expect_equal(cleared[["entries"]], 0)
})

test_that("selection runtime settings accept exact axes and remain frozen for workers", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(selection.sampler = "coarse_corrected",
    selection.cache_max_bytes = 65536,
    selection.coarse_grid = c(log_weight = 0, variance = .02, mean = .03),
    selection.coarse_max_rules = 3)
  expected <- list(cache_max_bytes = 65536, cache_budget_request = 65536,
    available_memory_bytes = NA_real_, sampler = "coarse_corrected",
    coarse_grid = c(mean = .03, variance = .02, log_weight = 0),
    coarse_max_rules = 3L)
  expect_identical(.selection_runtime_settings(), expected)
  setup <- unserialize(serialize(.selection_runtime_setup(), NULL))
  expect_identical(attr(setup, "selection_runtime"), expected)

  RoBMA.options(selection.sampler = "default", selection.cache_max_bytes = 0,
    selection.coarse_grid = c(mean = 0, variance = 0, log_weight = 0))
  context <- list(phase = "start", role = "local", chains = 10L,
    processes = 1L, process_id = 1L, process_chains = 10L, parallel = FALSE)
  expect_invisible(setup(context))
  native <- .Call("RoBMA_selnorm_sampler_control", NULL, NULL, FALSE, PACKAGE = "RoBMA")
  expect_true(native[["enabled"]])
  expect_equal(selection_cache_info()[["capacity_bytes"]], c(65536, NA, NA))
  context$parallel <- TRUE
  context$processes <- 3L
  context$role <- "coordinator"
  context$process_id <- 0L
  context$process_chains <- 0L
  setup(context)
  expect_equal(selection_cache_info()[["capacity_bytes"]][1L], 0)
  shares <- vapply(seq_len(3L), function(worker) {
    context$role <- "worker"
    context$process_id <- worker
    context$process_chains <- c(4L, 3L, 3L)[worker]
    setup(context)
    selection_cache_info()[["capacity_bytes"]][1L]
  }, numeric(1L))
  expect_equal(shares, floor(65536 / 10) * c(4, 3, 3))
  expect_lte(sum(shares), 65536)
  context$phase <- "finish"
  setup(context)
  expect_equal(selection_cache_info()[["capacity_bytes"]][1L], 65536)
  # The transported callback configures the native runtime without mutating
  # process-local R options, which may still be defaults in a fresh worker.
  expect_identical(RoBMA.get_option("selection.sampler"), "default")
  expect_identical(attr(setup, "selection_runtime"), expected)
})

test_that("invalid selection option updates leave R and native settings unchanged", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(selection.sampler = "coarse_corrected",
    selection.cache_max_bytes = 65536,
    selection.coarse_grid = c(mean = .05, variance = .05, log_weight = .01),
    selection.coarse_max_rules = 3L)
  before_options <- RoBMA.options()
  before_cache <- selection_cache_info(component = "all")
  before_sampler <- .Call("RoBMA_selnorm_sampler_control", NULL, NULL, FALSE, PACKAGE = "RoBMA")
  grid_message <- paste0("Option 'selection.coarse_grid' must contain finite nonnegative ",
    "'mean', 'variance', and 'log_weight' steps.")
  cases <- list(
    list(args = list(selection.cache_max_bytes = 0, selection.sampler = "coarse"),
      message = "Option 'selection.sampler' must be 'default' or 'coarse_corrected'."),
    list(args = list(selection.sampler = "default",
      selection.coarse_grid = c(mean = .1, variance = .1, variance = .1)), message = grid_message),
    list(args = list(selection.cache_max_bytes = 0,
      selection.coarse_grid = c(mean = .1, variance = Inf, log_weight = .1)), message = grid_message),
    list(args = list(selection.cache_max_bytes = 0, selection.coarse_max_rules = 2L),
      message = "Option 'selection.coarse_max_rules' must be an integer >= 3."),
    list(args = list(selection.sampler = "default", selection.cache_max_bytes = .5),
      message = paste0("Option 'selection.cache_max_bytes' must be a nonnegative ",
        "whole number of bytes representable exactly in R."))
  )
  for (case in cases) {
    error <- tryCatch(do.call(RoBMA.options, case$args), error = identity)
    expect_s3_class(error, "error")
    expect_identical(conditionMessage(error), case$message)
    expect_null(conditionCall(error))
    expect_identical(RoBMA.options(), before_options)
    expect_identical(selection_cache_info(component = "all"), before_cache)
    expect_identical(.Call("RoBMA_selnorm_sampler_control", NULL, NULL, FALSE, PACKAGE = "RoBMA"),
      before_sampler)
  }
})

test_that("automatic shared cache budgets use currently available RAM once", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  available <- 8 * 1024^3
  testthat::local_mocked_bindings(
    ps_system_memory = function() list(total = 64 * 1024^3, avail = available),
    .package = "ps")
  RoBMA.options(selection.cache_max_bytes = "auto")
  setup <- .selection_runtime_setup()
  expect_equal(attr(setup, "selection_runtime")$cache_max_bytes, 2 * 1024^3)
  available <- 48 * 1024^3
  expect_equal(.selection_runtime_settings()$cache_max_bytes, 4 * 1024^3)
  RoBMA.options(selection.sampler = "coarse_corrected")
  expect_equal(selection_cache_info()$capacity_bytes[1L], 2 * 1024^3)
  # Already captured fits retain their resolved total as available RAM changes.
  expect_equal(attr(setup, "selection_runtime")$cache_max_bytes, 2 * 1024^3)
  available <- 0
  expect_equal(.selection_runtime_settings()$cache_max_bytes, 0)
  available <- NA_real_
  warning <- NULL
  settings <- withCallingHandlers(.selection_runtime_settings(), warning = function(w) {
    warning <<- w
    invokeRestart("muffleWarning")
  })
  expect_identical(conditionMessage(warning), paste0(
    "Automatic selection cache memory budget is unavailable; caching is disabled. ",
    "Set 'selection.cache_max_bytes' to a numeric byte budget."))
  expect_null(conditionCall(warning))
  expect_equal(settings$cache_max_bytes, 0)
  RoBMA.options(selection.cache_max_bytes = 1024)
  expect_equal(.selection_runtime_settings()$cache_max_bytes, 1024)
  available <- 48 * 1024^3
})
