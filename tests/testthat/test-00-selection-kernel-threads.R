# ============================================================================ #
# test-00-selection-kernel-threads.R
# ============================================================================ #

context("Selection kernel threading")
skip_on_cran()

source(testthat::test_path("helper-selection-kernel.R"))


test_that("threaded selection kernels return identical values at any thread count", {

  .with_kernel_thread_worker({

    expect_true(.has_native_selnorm_kernel())
    expect_true(.has_native_selnorm_log_norm_delta())

    inputs <- .kernel_thread_inputs(S = 6000L)
    expect_identical(inputs[["spec"]][["n_bins"]], 4L)
    spec   <- inputs[["spec"]]
    mean   <- inputs[["mean"]]
    sd     <- inputs[["sd"]]
    values <- seq(-.20, .40, length.out = 7L)
    basis  <- matrix(
      seq(.5, 1.5, length.out = inputs[["S"]] * inputs[["K"]]),
      nrow = inputs[["S"]]
    )
    current <- seq(-.05, .25, length.out = inputs[["S"]])

    current_log_norm <- .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = inputs[["sei"]],
      omega          = inputs[["omega"]],
      selection_spec = spec,
      kernel_mode    = inputs[["kernel_mode"]]
    )

    .expect_thread_invariant("loglik_matrix", function() {
      .selnorm_kernel_loglik_matrix(
        yi             = inputs[["yi"]],
        mu_num         = mean,
        sigma_num      = sd,
        sei            = inputs[["sei"]],
        omega          = inputs[["omega"]],
        selection_spec = spec,
        kernel_mode    = inputs[["kernel_mode"]]
      )
    }, rows = inputs[["S"]], work_per_row = inputs[["K"]] * inputs[["spec"]][["n_bins"]])
    .expect_thread_invariant("loglik_row_sum", function() {
      .selnorm_kernel_loglik_row_sum(
        yi             = inputs[["yi"]],
        mu_num         = mean,
        sigma_num      = sd,
        sei            = inputs[["sei"]],
        omega          = inputs[["omega"]],
        selection_spec = spec,
        kernel_mode    = inputs[["kernel_mode"]]
      )
    }, rows = inputs[["S"]], work_per_row = inputs[["K"]] * inputs[["spec"]][["n_bins"]])
    .expect_thread_invariant("log_norm_matrix", function() {
      .selnorm_kernel_log_norm_matrix(
        mean           = mean,
        sd             = sd,
        sei            = inputs[["sei"]],
        omega          = inputs[["omega"]],
        selection_spec = spec,
        kernel_mode    = inputs[["kernel_mode"]]
      )
    }, rows = inputs[["S"]], work_per_row = inputs[["K"]] * inputs[["spec"]][["n_bins"]])
    .expect_thread_invariant("log_norm_delta_grid", function() {
      .selnorm_kernel_log_norm_delta_grid(
        mean             = mean,
        sd               = sd,
        basis            = basis,
        current_log_norm = current_log_norm,
        current          = current,
        values           = values,
        sei              = inputs[["sei"]],
        omega            = inputs[["omega"]],
        selection_spec   = spec,
        kernel_mode      = inputs[["kernel_mode"]]
      )
    }, rows = inputs[["S"]], work_per_row = length(values) * inputs[["K"]] * inputs[["spec"]][["n_bins"]])
    .expect_thread_invariant("cdf_matrix", function() {
      .selnorm_kernel_cdf_matrix(
        q              = inputs[["yi"]],
        mean           = mean,
        sd             = sd,
        sei            = inputs[["sei"]],
        omega          = inputs[["omega"]],
        selection_spec = spec,
        kernel_mode    = inputs[["kernel_mode"]]
      )
    }, rows = inputs[["S"]], work_per_row = inputs[["K"]] * inputs[["spec"]][["n_bins"]])
    .expect_thread_invariant("moments_matrix", function() {
      .selnorm_kernel_moments_matrix(
        mean           = mean,
        sd             = sd,
        sei            = inputs[["sei"]],
        omega          = inputs[["omega"]],
        selection_spec = spec,
        kernel_mode    = inputs[["kernel_mode"]]
      )
    }, rows = inputs[["S"]], work_per_row = inputs[["K"]] * inputs[["spec"]][["n_bins"]])
  })
})


test_that("threaded zcurve kernels return identical values at any thread count", {

  .with_kernel_thread_worker({

    expect_true(.has_native_zplot_density(selection = TRUE))
    expect_true(.has_native_zplot_threshold())

    inputs     <- .kernel_thread_inputs(S = 6000L)
    expect_identical(inputs[["spec"]][["n_bins"]], 4L)
    selection  <- .kernel_thread_selection_context(inputs)
    z_sequence <- seq(-4, 4, length.out = 41L)

    .expect_thread_invariant("zcurve_normal_density_matrix", function() {
      .zplot_normal_density_matrix(
        z_sequence = z_sequence,
        mean       = inputs[["mean"]],
        sd         = inputs[["sd"]],
        sei        = inputs[["sei"]]
      )
    }, rows = inputs[["S"]], work_per_row = length(z_sequence) * inputs[["K"]])
    .expect_thread_invariant("zcurve_density_matrix (fitted)", function() {
      .zplot_selnorm_density_matrix(
        z_sequence        = z_sequence,
        mean              = inputs[["mean"]],
        sd                = inputs[["sd"]],
        sei               = inputs[["sei"]],
        selection_context = selection,
        extrapolate       = FALSE
      )
    }, rows = inputs[["S"]], work_per_row = length(z_sequence) * inputs[["K"]])
    .expect_thread_invariant("zcurve_density_matrix (extrapolated)", function() {
      .zplot_selnorm_density_matrix(
        z_sequence        = z_sequence,
        mean              = inputs[["mean"]],
        sd                = inputs[["sd"]],
        sei               = inputs[["sei"]],
        selection_context = selection,
        extrapolate       = TRUE
      )
    }, rows = inputs[["S"]], work_per_row = length(z_sequence) * inputs[["K"]])
    .expect_thread_invariant("zcurve_density_matrix (pair)", function() {
      .zplot_selnorm_density_pair(
        z_sequence        = z_sequence,
        mean              = inputs[["mean"]],
        sd                = inputs[["sd"]],
        sei               = inputs[["sei"]],
        selection_context = selection
      )
    }, rows = inputs[["S"]], work_per_row = length(z_sequence) * inputs[["K"]])
    .expect_thread_invariant("zcurve_threshold_summary", function() {
      .zplot_selnorm_threshold_summary(
        z_threshold       = 1.96,
        mean              = inputs[["mean"]],
        sd                = inputs[["sd"]],
        sei               = inputs[["sei"]],
        selection_context = selection,
        extrapolate       = FALSE
      )
    }, rows = inputs[["S"]], work_per_row = inputs[["K"]] * inputs[["spec"]][["n_bins"]])
  })
})


# ---------------------------------------------------------------------------- #
# Work-based row schedule
# ---------------------------------------------------------------------------- #

.row_schedule <- function(rows, work_per_row) {

  return(.Call("RoBMA_selnorm_row_schedule", as.numeric(rows),
               as.numeric(work_per_row), PACKAGE = "RoBMA"))
}


.current_native_thread_budget <- function() {

  current <- .native_threads_configure(1L)
  if (!is.null(current)) {
    .native_threads_configure(current)
  }
  current
}


test_that("native thread scopes restore their process budget", {

  skip_if_not(is.loaded("RoBMA_selnorm_set_native_threads", PACKAGE = "RoBMA"))

  previous_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous_options), add = TRUE)
  RoBMA.options(native_threads = NA_integer_, max_cores = 3L)
  .native_threads_configure(5L)

  serial_object   <- list(fit_control = list(parallel = FALSE))
  parallel_object <- list(fit_control = list(parallel = TRUE))

  expect_identical(
    .with_native_threads(serial_object, .current_native_thread_budget()),
    1L
  )
  expect_identical(.current_native_thread_budget(), 5L)

  expect_error(
    .with_native_threads(serial_object, stop("scope failure", call. = FALSE)),
    "scope failure",
    fixed = TRUE
  )
  expect_identical(.current_native_thread_budget(), 5L)

  expect_identical(
    .with_native_threads(parallel_object, .current_native_thread_budget()),
    3L
  )
  expect_identical(.current_native_thread_budget(), 5L)
})


test_that("object-facing post-fit entry points inherit a parallel fit budget", {

  skip_if_not(is.loaded("RoBMA_selnorm_set_native_threads", PACKAGE = "RoBMA"))

  previous_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous_options), add = TRUE)
  RoBMA.options(native_threads = NA_integer_, max_cores = 3L)
  .native_threads_configure(5L)
  object <- structure(
    list(fit_control = list(parallel = TRUE), fit = list()),
    class = "brma"
  )

  testthat::local_mocked_bindings(
    .check_log_lik_target_available = function(...) invisible(TRUE),
    .log_lik_estimate.brma = function(...) .current_native_thread_budget(),
    .predict_brma_context = function(...) list(),
    .predict_brma_from_context = function(...) .current_native_thread_budget(),
    .check_legacy_level_arg = function(...) invisible(TRUE),
    .normalize_funnel_max_samples = function(x) x,
    .set_dots_funnel = function(dots) list(as_data = TRUE),
    .is_mods = function(...) FALSE,
    .is_scale = function(...) FALSE,
    .funnel_common_heterogeneity = function(...) list(common = TRUE),
    .funnel_data_outcome = function(...) .current_native_thread_budget(),
    .get_posterior_samples = function(...) matrix(1, nrow = 2L, ncol = 1L),
    .thin_sample_rows = function(...) NULL,
    .zplot_requires_selection_marginal = function(...) TRUE,
    .zplot_selection_marginal = function(...) .current_native_thread_budget(),
    .package = "RoBMA"
  )

  expect_identical(.log_lik.brma(object), 3L)
  expect_identical(predict.brma(object), 3L)
  expect_identical(funnel.brma(object, as_data = TRUE), 3L)
  expect_identical(
    .zplot_density_pair(object, 0, max_samples = Inf),
    3L
  )
  expect_identical(.current_native_thread_budget(), 5L)
})


test_that("a fresh RoBMA worker starts with one native thread", {

  skip_on_cran()
  cluster <- parallel::makePSOCKcluster(1L, rscript_args = "--vanilla")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(paths) {
    .libPaths(paths)
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  }, .libPaths())
  parallel::clusterEvalQ(cluster, loadNamespace("RoBMA"))

  initial <- parallel::clusterEvalQ(cluster, {
    configure <- getFromNamespace(".native_threads_configure", "RoBMA")
    previous  <- configure(1L)
    configure(previous)
    previous
  })
  expect_identical(initial, list(1L))
})


test_that("the row schedule keeps small batches serial and sizes large regions", {

  skip_if_not(is.loaded("RoBMA_selnorm_row_schedule", PACKAGE = "RoBMA"))

  # Every shape below the threshold runs the entry point's own serial loop, so
  # no thread budget can make it slower than one thread.
  small <- list(
    c(rows = 16L, work = 5 * 4),        # 16 rows, 5 observations, 4 bins
    c(rows = 64L, work = 5 * 4),
    c(rows = 16L, work = 200 * 4),
    c(rows = 500L, work = 12 * 4),
    c(rows = 64L, work = 20 * 5 * 4),   # 20 grid points
    c(rows = 16L, work = 20 * 12)       # zcurve grid, no bins
  )
  # The shapes the post-fit batches spend their time in do thread, and one
  # parallel region covers enough rows that the region overhead is amortized.
  large <- list(
    c(rows = 20000L, work = 12 * 4),
    c(rows = 60000L, work = 48 * 4),
    c(rows = 4000L, work = 20 * 48 * 4),
    c(rows = 20000L, work = 350 * 5)
  )

  # OpenMP reads its available thread cap during DLL initialization. Query
  # schedules in a fresh two-thread worker even when this test process is serial.
  schedules <- (function() {
    withr::local_envvar(c(OMP_NUM_THREADS = "2", OMP_THREAD_LIMIT = "2",
      OMP_DYNAMIC = "FALSE"))
    cluster <- parallel::makePSOCKcluster(1L, rscript_args = "--vanilla")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(paths)
      Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    }, .libPaths())
    parallel::clusterEvalQ(cluster, loadNamespace("RoBMA"))
    parallel::clusterCall(cluster, function(small, large) {
      row_schedule <- function(rows, work_per_row) {
        .Call("RoBMA_selnorm_row_schedule", as.numeric(rows),
              as.numeric(work_per_row), PACKAGE = "RoBMA")
      }
      query_shapes <- function(shapes) {
        lapply(shapes, function(shape) row_schedule(shape[["rows"]], shape[["work"]]))
      }
      RoBMA::RoBMA.options(native_threads = 31L)
      small_schedules <- query_shapes(small)
      large_schedules <- query_shapes(large)
      RoBMA::RoBMA.options(native_threads = 8L)
      eight_schedule <- row_schedule(60000L, 48 * 4)
      RoBMA::RoBMA.options(native_threads = 1L)
      one_schedule <- row_schedule(60000L, 48 * 4)
      list(openmp = .Call("RoBMA_native_unwind_probe", "openmp", NULL, 1L,
                         PACKAGE = "RoBMA"),
           small = small_schedules, large = large_schedules,
           eight = eight_schedule, one = one_schedule)
    }, small, large)[[1L]]
  })()

  skip_if_not(isTRUE(schedules[["openmp"]]), "This build has no OpenMP support.")

  for (i in seq_along(small)) {
    shape <- small[[i]]
    schedule <- schedules[["small"]][[i]]
    expect_identical(schedule[["threads"]], 1L,
      info = paste("rows", shape[["rows"]], "work", shape[["work"]]))
  }
  for (i in seq_along(large)) {
    shape <- large[[i]]
    schedule <- schedules[["large"]][[i]]
    expect_gt(schedule[["threads"]], 1L)
    expect_gte(schedule[["chunk_rows"]], schedule[["threads"]])
    # Either the region carries its full share of work or it covers the batch.
    expect_true(
      schedule[["chunk_rows"]] == shape[["rows"]] ||
        as.numeric(schedule[["chunk_rows"]]) * shape[["work"]] >=
          schedule[["threads"]] * 1e6,
      info = paste("rows", shape[["rows"]], "work", shape[["work"]])
    )
  }

  # A smaller budget is never exceeded, and one thread stays one thread.
  expect_lte(schedules[["eight"]][["threads"]], 8L)
  expect_identical(schedules[["one"]][["threads"]], 1L)
})


test_that("a batch below the work gate is identical at every thread budget", {

  skip_if_not(.has_native_selnorm_kernel())

  # 20 rows and 5 observations are below the work gate, so this batch keeps the
  # serial loop; the values must still be identical at every budget.
  inputs    <- .kernel_thread_inputs(S = 20L, K = 5L)
  selection <- .kernel_thread_selection_context(inputs)

  .expect_thread_invariant("small loglik_row_sum", function() {
    .selnorm_kernel_loglik_row_sum(
      yi             = inputs[["yi"]],
      mu_num         = inputs[["mean"]],
      sigma_num      = inputs[["sd"]],
      sei            = inputs[["sei"]],
      omega          = inputs[["omega"]],
      selection_spec = inputs[["spec"]],
      kernel_mode    = inputs[["kernel_mode"]]
    )
  }, threads = c(1L, 2L, 8L, 31L))

  expect_identical(
    .row_schedule(20L, length(inputs[["sei"]]) * inputs[["spec"]][["n_bins"]])[["threads"]],
    1L
  )
})
