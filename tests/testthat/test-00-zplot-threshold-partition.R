source(testthat::test_path("helper-selection-kernel.R"))


test_that("positive threshold partitions agree with independent interval masses", {

  sei <- c(.1, 1, 3)
  means <- matrix(c(-12, -3, 0, .5, 4, 12), 6L, 3L)
  scales <- matrix(rep(c(.3, 1, 2), each = 6L), 6L)
  omega <- rbind(c(1, .7, .35, .2), c(.2, .8, .3, 1),
    c(.3, .3, .1, .1), c(1, 1, 1, 1), c(1, 1e-12, .7, .01),
    c(.1, .2, .9, .5))
  for (direction in c("positive", "negative")) {
    selection <- .test_step_spec(rep(0, 3L), sei, direction)
    selection$omega <- omega
    selection$alpha <- 0
    selection$phack_kind <- 0L
    selection$kernel_mode <- SELKERNEL_STEP
    selection$telescope_probabilities <- TRUE
    log_normalizer <- .test_step_log_norm_reference(means, scales, sei, omega, selection)
    for (threshold in c(0, 1e-8, .7, stats::qnorm(.975), 6, Inf)) {
      reference <- matrix(0, 6L, 3L)
      for (s in seq_len(6L)) for (k in seq_len(3L)) {
        tail_mass <- vapply(seq_len(ncol(omega)), function(bin) {
          bounds <- sort(selection$sign *
            c(selection$z_lower[bin], selection$z_upper[bin]) * sei[k])
          lower <- .test_interval_log_prob(bounds[1L],
            min(bounds[2L], -threshold * sei[k]), means[s, k], scales[s, k])
          upper <- .test_interval_log_prob(max(bounds[1L], threshold * sei[k]),
            bounds[2L], means[s, k], scales[s, k])
          log(omega[s, bin]) + .test_logsumexp(c(lower, upper))
        }, numeric(1))
        reference[s, k] <- exp(.test_logsumexp(tail_mass) - log_normalizer[s, k])
      }
      observed <- .zplot_selnorm_threshold_summary(
        threshold, means, scales, sei, selection, FALSE)
      expect_equal(observed$EDR, rowMeans(reference), tolerance = 1e-12)
      expect_equal(observed$weights, rowMeans(exp(-log_normalizer)), tolerance = 1e-12)
    }
  }
})


test_that("threshold partitions retain zero-weight and collapsed-score gates", {

  selection <- .test_step_spec(0, 1)
  selection$alpha <- 0
  selection$phack_kind <- 0L
  selection$kernel_mode <- SELKERNEL_STEP
  selection$telescope_probabilities <- TRUE
  selection$omega <- matrix(c(1, .7, .35, .2), 1L)
  collapsed <- .zplot_selnorm_threshold_summary(
    1.96, matrix(-1e308), matrix(1), 1, selection, FALSE)
  expect_true(is.na(collapsed$EDR))
  expect_true(is.na(collapsed$weights))
  infinite_score <- .zplot_selnorm_threshold_summary(
    1.96, matrix(0), matrix(1e-310), 1, selection, FALSE)
  expect_true(is.na(infinite_score$EDR))
  expect_true(is.na(infinite_score$weights))

  # Equal adjacent weights remove their boundaries before the score check.
  selection$omega[,] <- .2
  coalesced <- .zplot_selnorm_threshold_summary(
    1.96, matrix(-1e308), matrix(1), 1, selection, FALSE)
  expect_identical(coalesced, list(EDR = 1, weights = 5))

  for (weights in list(c(1, 0, 0, 0), c(0, 0, 0, 1), c(1, .7, 0, .2))) {
    selection$omega <- matrix(weights, 1L)
    for (mean in c(-1e308, -4, 0, 4, 1e308)) {
      observed <- .zplot_selnorm_threshold_summary(
        1.96, matrix(mean), matrix(1), 1, selection, FALSE)
      reference <- selection
      reference$telescope_probabilities <- FALSE
      expected <- .zplot_selnorm_threshold_summary(
        1.96, matrix(mean), matrix(1), 1, reference, FALSE)
      expect_equal(observed, expected, tolerance = 1e-12)
    }
  }
})


test_that("repeated threshold inputs retain observation order and exact keys", {

  sei <- c(.1, 1, .1, 2, 1, .1, 2, 1)
  means <- rbind(c(.2, .3, .2, -.4, .3, .5, -.4, .3),
    c(-1, 2, -1, -4, 2, -1, -4, 2))
  scales <- rbind(c(.5, 1, .5, 2, 1, .5, 3, 1),
    c(.3, 1, .3, 4, 2, .3, 4, 1))
  selection <- .test_step_spec(rep(0, length(sei)), sei)
  selection$omega <- rbind(c(1, .7, .35, .2), c(1, 0, .8, .2))
  selection$alpha <- 0
  selection$phack_kind <- 0L
  for (mode in list(c(SELKERNEL_STEP, SELKERNEL_NORMAL), rep(SELKERNEL_STEP, 2L))) {
    selection$kernel_mode <- mode
    for (extrapolate in c(FALSE, TRUE)) {
      total_probability <- total_weight <- numeric(nrow(means))
      for (k in seq_along(sei)) {
        value <- .zplot_selnorm_threshold_summary(1.96,
          means[, k, drop = FALSE], scales[, k, drop = FALSE], sei[k], selection, extrapolate)
        total_probability <- total_probability + value$EDR * if (extrapolate) value$weights else 1
        total_weight <- total_weight + value$weights
      }
      expected <- list(EDR = total_probability / if (extrapolate) total_weight else length(sei),
        weights = total_weight / length(sei))
      observed <- .zplot_selnorm_threshold_summary(1.96, means, scales, sei, selection, extrapolate)
      expect_identical(observed, expected)
    }
  }
})


test_that("repeated threshold inputs do not mask invalid cells", {

  sei <- rep(.5, 6L)
  means <- matrix(.2, 4L, 6L)
  scales <- matrix(.7, 4L, 6L)
  selection <- .test_step_spec(rep(0, length(sei)), sei)
  selection$omega <- matrix(rep(c(1, .7, .35, .2), each = 4L), 4L)
  selection$alpha <- 0
  selection$phack_kind <- 0L
  selection$kernel_mode <- SELKERNEL_STEP
  means[1L, 2L] <- NA_real_
  means[2L, 3L] <- Inf
  scales[3L, 4L] <- 0
  scales[4L, 5L] <- NaN
  for (extrapolate in c(FALSE, TRUE)) {
    observed <- .zplot_selnorm_threshold_summary(1.96, means, scales, sei, selection, extrapolate)
    expect_true(all(is.na(observed$EDR)))
    expect_true(all(is.na(observed$weights)))
  }
  means[,] <- .2
  scales[,] <- .7
  for (value in c(NA_real_, Inf, 0)) {
    sei[3L] <- value
    observed <- .zplot_selnorm_threshold_summary(1.96, means, scales, sei, selection, FALSE)
    expect_true(all(is.na(observed$EDR)))
    expect_true(all(is.na(observed$weights)))
  }
})


test_that("threshold reuse is identical with one or two native threads", {

  skip_on_cran()
  .with_kernel_thread_worker({
    S <- 64L
    K <- 200L
    sei <- rep(c(.1, .3, .5, 1), length.out = K)
    means <- matrix(seq(-.5, .5, length.out = S), S, K)
    scales <- matrix(rep(sei + .3, each = S), S)
    selection <- .test_step_spec(rep(0, K), sei)
    selection$omega <- cbind(1, seq(.05, .8, length.out = S), 0, .2)
    selection$alpha <- 0
    selection$phack_kind <- 0L
    selection$kernel_mode <- rep(c(SELKERNEL_STEP, SELKERNEL_NORMAL), length.out = S)
    for (extrapolate in c(FALSE, TRUE)) {
      .expect_thread_invariant("repeated threshold cells", function() {
        .zplot_selnorm_threshold_summary(1.96, means, scales, sei, selection, extrapolate)
      }, threads = c(1L, 2L), rows = S, work_per_row = K * selection$n_bins)
    }
  })
})
