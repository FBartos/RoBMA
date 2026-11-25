probe_state("test-00-selection-cache.R")
test_that("selection cache inspection initializes its configured budget", {

  old <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old), add = TRUE)
  RoBMA.options(selection.cache_max_bytes = 65536)
  .Call("RoBMA_selnorm_cache_control", 0, 3L, PACKAGE = "RoBMA")
  RoBMA.private[["selection_runtime_initialized"]] <- FALSE
  info <- selection_cache_info()
  expect_identical(info$capacity_bytes[1L], 65536)
  expect_true(RoBMA.private[["selection_runtime_initialized"]])
  expect_identical(selection_cache_info(), info)
})


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
  # Exercise many admissions after the cache is full; allocator/eviction bugs
  # appear on repeated turnover, not on the first cache fill.
  for (index in seq_len(4000L)) {
    pressure_arguments[[2L]][] <- index * 1e-4
    last <- evaluate(pressure_arguments)
  }
  pressure <- info()
  expect_true(is.finite(last[["log_normalizer"]]))
  expect_lte(diagnostic_error(last), .005)
  expect_gt(pressure[["evictions"]], 0)
  expect_equal(pressure[["allocation_failures"]], 0)
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
  RoBMA.options(selection.sampler = "slice.coarse",
    selection.cache_max_bytes = 65536,
    selection.coarse_grid = c(log_weight = 0, variance = .02, mean = .03),
    selection.coarse_max_rules = 3, selection.coarse_min_acceptance = .8)
  expected <- list(cache_max_bytes = 65536, cache_budget_request = 65536,
    available_memory_bytes = NA_real_, sampler = "slice.coarse",
    coarse_grid = c(mean = .03, variance = .02, log_weight = 0),
    coarse_max_rules = 3L, coarse_min_acceptance = .8)
  expect_identical(.selection_runtime_settings(), expected)
  setup <- unserialize(serialize(.selection_runtime_setup(), NULL))
  expect_identical(attr(setup, "selection_runtime"), expected)

  RoBMA.options(selection.sampler = "slice", selection.cache_max_bytes = 0,
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
  expect_identical(RoBMA.get_option("selection.sampler"), "slice")
  expect_identical(attr(setup, "selection_runtime"), expected)
})

test_that("invalid selection option updates leave R and native settings unchanged", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(selection.sampler = "slice.coarse",
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
      message = "Option 'selection.sampler' must be 'slice.coarse', 'slice' or 'jags.default'."),
    list(args = list(selection.sampler = "slice",
      selection.coarse_grid = c(mean = .1, variance = .1, variance = .1)), message = grid_message),
    list(args = list(selection.cache_max_bytes = 0,
      selection.coarse_grid = c(mean = .1, variance = Inf, log_weight = .1)), message = grid_message),
    list(args = list(selection.cache_max_bytes = 0, selection.coarse_max_rules = 2L),
      message = "Option 'selection.coarse_max_rules' must be an integer >= 3."),
    list(args = list(selection.sampler = "slice", selection.cache_max_bytes = .5),
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

# The warnings an expression raises, muffled, with its value.
.sampler_option_warnings <- function(expr) {

  caught <- list()
  value  <- withCallingHandlers(expr, warning = function(w) {
    caught[[length(caught) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = caught)
}

test_that("selection.sampler takes three values and maps the deprecated ones with one warning", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  expect_identical(.RoBMA_option_schema[["selection.sampler"]][["default"]], "slice.coarse")

  for (value in c("slice.coarse", "slice", "jags.default")) {
    result <- .sampler_option_warnings(RoBMA.options(selection.sampler = value))
    expect_length(result$warnings, 0L)
    expect_identical(result$value[["selection.sampler"]], value)
    expect_identical(RoBMA.get_option("selection.sampler"), value)
  }

  deprecated <- list(
    coarse_corrected = list(value = "slice.coarse", message = paste0(
      "'selection.sampler = \"coarse_corrected\"' is deprecated; use \"slice.coarse\".")),
    default          = list(value = "jags.default", message = paste0(
      "'selection.sampler = \"default\"' is deprecated and now means \"jags.default\" (only the samplers JAGS ",
      "assigns, as in RoBMA 4.1.1); use \"jags.default\", or \"slice\" to keep the gated and factor slice samplers."))
  )
  for (alias in names(deprecated)) {
    RoBMA.options(selection.sampler = "slice")
    result <- .sampler_option_warnings(RoBMA.options(selection.sampler = alias, selection.coarse_max_rules = 3L))
    expect_length(result$warnings, 1L)
    condition <- result$warnings[[1L]]
    expect_identical(class(condition), c("RoBMA_deprecated_option_value", "RoBMA_deprecated", "deprecatedWarning",
                                         "warning", "condition"), info = alias)
    expect_identical(conditionMessage(condition), deprecated[[alias]][["message"]], info = alias)
    expect_null(conditionCall(condition))
    # the option stores and returns the new value, and fits capture it
    expect_identical(result$value[["selection.sampler"]], deprecated[[alias]][["value"]], info = alias)
    expect_identical(RoBMA.get_option("selection.sampler"), deprecated[[alias]][["value"]], info = alias)
    expect_identical(.selection_runtime_settings()[["sampler"]], deprecated[[alias]][["value"]], info = alias)
  }

  # a named value is stored and captured as the bare value
  for (value in c("slice.coarse", "slice", "jags.default")) {
    RoBMA.options(selection.sampler = setdiff(c("slice", "jags.default"), value)[[1L]])
    result <- .sampler_option_warnings(RoBMA.options(selection.sampler = c(stock = value)))
    expect_length(result$warnings, 0L)
    expect_identical(result$value[["selection.sampler"]], value, info = value)
    expect_identical(RoBMA.get_option("selection.sampler"), value, info = value)
    expect_identical(.selection_runtime_settings()[["sampler"]], value, info = value)
  }
  result <- .sampler_option_warnings(RoBMA.options(selection.sampler = c(stock = "default")))
  expect_length(result$warnings, 1L)
  expect_identical(RoBMA.get_option("selection.sampler"), "jags.default")

  # the deprecation warning waits for the other options of the call: a call
  # that fails on a later option warns about nothing and stores nothing
  RoBMA.options(selection.sampler = "slice")
  before <- RoBMA.options()
  result <- .sampler_option_warnings(tryCatch(RoBMA.options(selection.sampler = "default", max_cores = 0),
                                              error = identity))
  expect_s3_class(result$value, "error")
  expect_identical(conditionMessage(result$value), "Option 'max_cores' must be an integer >= 1.")
  expect_length(result$warnings, 0L)
  expect_identical(RoBMA.options(), before)

  before <- RoBMA.options()
  for (value in list("coarse", "Slice", "jags", NA_character_, c("slice", "slice.coarse"), 1, TRUE, NULL)) {
    expect_error(do.call(RoBMA.options, list(selection.sampler = value)),
                 "Option 'selection.sampler' must be 'slice.coarse', 'slice' or 'jags.default'.", fixed = TRUE)
  }
  expect_identical(RoBMA.options(), before)
})

test_that("frozen sampler settings of earlier versions keep their meaning without a warning", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  native_enabled <- function() {
    .Call("RoBMA_selnorm_sampler_control", NULL, NULL, FALSE, PACKAGE = "RoBMA")[["enabled"]]
  }

  # settings frozen in fits of RoBMA 4.1.1 and of this version
  corrected <- c(coarse_corrected = TRUE, default = FALSE,
                 slice.coarse = TRUE, slice = FALSE, jags.default = FALSE)
  for (sampler in names(corrected)) {
    expect_identical(.selection_sampler_corrected(sampler), corrected[[sampler]], info = sampler)
    settings <- .selection_runtime_settings()
    settings[["sampler"]] <- sampler
    if (sampler %in% c("coarse_corrected", "default")) {
      settings[["coarse_min_acceptance"]] <- NULL
    }
    # start from the other state, so that the frozen value sets it
    .Call("RoBMA_selnorm_sampler_control", !corrected[[sampler]], NULL, FALSE, PACKAGE = "RoBMA")
    setup <- .selection_runtime_setup()
    environment(setup)[["settings"]] <- settings
    result <- .sampler_option_warnings(setup(list(phase = "start", role = "local", chains = 1L, processes = 1L,
                                                  process_id = 1L, process_chains = 1L, parallel = FALSE)))
    expect_length(result$warnings, 0L)
    expect_identical(native_enabled(), corrected[[sampler]], info = sampler)
  }
  expect_false(.selection_sampler_corrected(NULL))
  expect_false(.selection_sampler_corrected(NA_character_))
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
  RoBMA.options(selection.sampler = "slice.coarse")
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
