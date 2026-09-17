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
