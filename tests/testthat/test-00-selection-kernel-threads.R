# ============================================================================ #
# test-00-selection-kernel-threads.R
# ============================================================================ #

context("Selection kernel threading")
skip_on_cran()

source(testthat::test_path("helper-selection-kernel.R"))


# One synthetic batch with enough posterior rows to pass the parallel gate of
# robma_batch_threads() and with row-wise selection weights, kernel modes and
# scales, so the row arguments the kernels resolve really vary by row.
.kernel_thread_inputs <- function(S = 96L, K = 5L) {

  set.seed(20260917)
  yi   <- c(.12, -.05, .31, .22, -.18)[seq_len(K)]
  sei  <- c(.10, .18, .14, .22, .09)[seq_len(K)]
  spec <- .test_step_spec(yi, sei)

  mean <- matrix(stats::rnorm(S * K, mean = .15, sd = .10), nrow = S)
  sd   <- matrix(stats::runif(S * K, .08, .30), nrow = S)
  # A normal row every fourth draw: the selected and unselected branches of the
  # kernels must both be visited inside the parallel region.
  kernel_mode <- rep(SELKERNEL_STEP, S)
  kernel_mode[seq.int(4L, S, by = 4L)] <- SELKERNEL_NORMAL
  omega <- cbind(
    rep(1, S),
    stats::runif(S, .40, .95),
    stats::runif(S, .15, .60),
    stats::runif(S, .05, .40)
  )
  omega <- omega[, seq_len(spec[["n_bins"]]), drop = FALSE]
  omega[, 1L] <- 1

  return(list(
    yi          = yi,
    sei         = sei,
    spec        = spec,
    S           = S,
    K           = K,
    mean        = mean,
    sd          = sd,
    omega       = omega,
    kernel_mode = kernel_mode
  ))
}


.kernel_thread_selection_context <- function(inputs) {

  context <- inputs[["spec"]]
  context[["omega"]]       <- inputs[["omega"]]
  context[["alpha"]]       <- 0
  context[["phack_kind"]]  <- 0L
  context[["kernel_mode"]] <- inputs[["kernel_mode"]]
  context[["vector_rule"]] <- 0L
  context[["use_normal"]]  <- inputs[["kernel_mode"]] == SELKERNEL_NORMAL

  return(BayesTools::selection_context_validate(
    context,
    n_samples = inputs[["S"]]
  ))
}


# Evaluate one entry point at several thread budgets and require the values to
# be identical, not merely close: the row partition must never change a result.
.expect_thread_invariant <- function(label, evaluate,
                                     threads = c(1L, 2L, 8L)) {

  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)

  reference <- NULL
  for (count in threads) {
    RoBMA.options(native_threads = count)
    values <- unlist(evaluate(), use.names = FALSE)
    if (is.null(reference)) {
      expect_true(
        any(is.finite(values)),
        info = paste0(label, ": no finite value to compare")
      )
      reference <- values
    } else {
      expect_identical(
        values,
        reference,
        info = paste0(label, " at ", count, " threads")
      )
    }
  }

  return(invisible(reference))
}


test_that("threaded selection kernels return identical values at any thread count", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_selnorm_log_norm_delta())

  inputs <- .kernel_thread_inputs()
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
  })
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
  })
  .expect_thread_invariant("log_norm_matrix", function() {
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = inputs[["sei"]],
      omega          = inputs[["omega"]],
      selection_spec = spec,
      kernel_mode    = inputs[["kernel_mode"]]
    )
  })
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
  })
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
  })
  .expect_thread_invariant("moments_matrix", function() {
    .selnorm_kernel_moments_matrix(
      mean           = mean,
      sd             = sd,
      sei            = inputs[["sei"]],
      omega          = inputs[["omega"]],
      selection_spec = spec,
      kernel_mode    = inputs[["kernel_mode"]]
    )
  })
})


test_that("threaded zcurve kernels return identical values at any thread count", {

  skip_if_not(.has_native_zplot_density(selection = TRUE))
  skip_if_not(.has_native_zplot_threshold())

  inputs     <- .kernel_thread_inputs()
  selection  <- .kernel_thread_selection_context(inputs)
  z_sequence <- seq(-4, 4, length.out = 41L)

  .expect_thread_invariant("zcurve_normal_density_matrix", function() {
    .zplot_normal_density_matrix(
      z_sequence = z_sequence,
      mean       = inputs[["mean"]],
      sd         = inputs[["sd"]],
      sei        = inputs[["sei"]]
    )
  })
  .expect_thread_invariant("zcurve_density_matrix (fitted)", function() {
    .zplot_selnorm_density_matrix(
      z_sequence        = z_sequence,
      mean              = inputs[["mean"]],
      sd                = inputs[["sd"]],
      sei               = inputs[["sei"]],
      selection_context = selection,
      extrapolate       = FALSE
    )
  })
  .expect_thread_invariant("zcurve_density_matrix (extrapolated)", function() {
    .zplot_selnorm_density_matrix(
      z_sequence        = z_sequence,
      mean              = inputs[["mean"]],
      sd                = inputs[["sd"]],
      sei               = inputs[["sei"]],
      selection_context = selection,
      extrapolate       = TRUE
    )
  })
  .expect_thread_invariant("zcurve_density_matrix (pair)", function() {
    .zplot_selnorm_density_pair(
      z_sequence        = z_sequence,
      mean              = inputs[["mean"]],
      sd                = inputs[["sd"]],
      sei               = inputs[["sei"]],
      selection_context = selection
    )
  })
  .expect_thread_invariant("zcurve_threshold_summary", function() {
    .zplot_selnorm_threshold_summary(
      z_threshold       = 1.96,
      mean              = inputs[["mean"]],
      sd                = inputs[["sd"]],
      sei               = inputs[["sei"]],
      selection_context = selection,
      extrapolate       = FALSE
    )
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
  parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
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

  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)
  RoBMA.options(native_threads = 31L)

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
  for (shape in small) {
    schedule <- .row_schedule(shape[["rows"]], shape[["work"]])
    expect_identical(schedule[["threads"]], 1L,
      info = paste("rows", shape[["rows"]], "work", shape[["work"]]))
  }

  # The shapes the post-fit batches spend their time in do thread, and one
  # parallel region covers enough rows that the region overhead is amortized.
  large <- list(
    c(rows = 20000L, work = 12 * 4),
    c(rows = 60000L, work = 48 * 4),
    c(rows = 4000L, work = 20 * 48 * 4),
    c(rows = 20000L, work = 350 * 5)
  )
  for (shape in large) {
    schedule <- .row_schedule(shape[["rows"]], shape[["work"]])
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
  RoBMA.options(native_threads = 8L)
  expect_lte(.row_schedule(60000L, 48 * 4)[["threads"]], 8L)
  RoBMA.options(native_threads = 1L)
  expect_identical(.row_schedule(60000L, 48 * 4)[["threads"]], 1L)
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
