context("Selection kernel")
skip_on_cran()

test_that("selection-only prior_bias wrapper is delegated to BayesTools priors", {

  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  object <- bselmodel(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    measure                   = "SMD",
    prior_bias                = BayesTools::prior_bias(selection = selection),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  priors <- .create_fit_priors(object[["data"]], object[["priors"]])

  expect_true("bias" %in% names(priors))
  expect_true(BayesTools::is_prior_bias(priors[["bias"]]))
  expect_equal(BayesTools::JAGS_to_monitor(priors), c("mu", "tau", "omega"))
})

test_that("selection-only prior_bias wrapper is visible to zplot thresholds", {

  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  object <- bselmodel(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    measure                   = "SMD",
    prior_bias                = BayesTools::prior_bias(selection = selection),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  expect_equal(.zplot_threshold(object[["priors"]]), stats::qnorm(1 - .025))
})

test_that("zplot extrapolation colors deduplicate user line arguments", {

  zplot_object <- structure(list(), class = "zplot_brma")
  density_data <- data.frame(
    x     = c(0, 1),
    y     = c(1, 1),
    y_lCI = c(.9, .9),
    y_uCI = c(1.1, 1.1)
  )
  observed_col <- character()

  testthat::local_mocked_bindings(
    lines.zplot_brma = function(x, plot_type = "base",
                                probs = c(.025, .975), max_samples = 10000,
                                plot_ci = TRUE, extrapolate = FALSE,
                                from = -6, to = 6, by = 0.05,
                                length.out = NULL, col = "black",
                                as_data = FALSE, ...) {

      if (isTRUE(extrapolate)) {
        observed_col <<- c(observed_col, col)
      }
      return(density_data)
    },
    hist.zplot_brma = function(...) {

      graphics::plot(c(0, 1), c(0, 1), type = "n")
      invisible(NULL)
    },
    .package = "RoBMA"
  )

  file <- tempfile(fileext = ".png")
  grDevices::png(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)

  expect_silent(plot.zplot_brma(
    zplot_object,
    plot_fit = FALSE,
    plot_ci  = FALSE,
    col      = "red"
  ))
  expect_equal(observed_col, "red")

  observed_col <- character()
  expect_silent(plot.zplot_brma(
    zplot_object,
    plot_fit           = FALSE,
    plot_ci            = FALSE,
    dots_extrapolation = list(col = "green")
  ))
  expect_equal(observed_col, "green")
})

test_that("native zplot density matches selected-normal R reference", {

  skip_if_not(.has_native_zplot_density(selection = TRUE))

  S        <- 4L
  K        <- 3L
  sei      <- c(.07, .13, .29)
  mu       <- matrix(c(
    -.20, .05, .30,
     .10, .12, .18,
     .40, .00, -.15,
     .05, -.25, .20
  ), nrow = S, byrow = TRUE)
  tau      <- matrix(c(
    .05, .07, .10,
    .08, .06, .11,
    .03, .09, .12,
    .04, .10, .15
  ), nrow = S, byrow = TRUE)
  z_values <- c(
    -3,
    0,
    stats::qnorm(.05, lower.tail = FALSE),
    stats::qnorm(.025, lower.tail = FALSE),
    3
  )
  spec     <- .test_step_spec(yi = c(.10, .25, .40), sei = sei)

  selection <- spec
  selection[["omega"]]       <- matrix(c(
    1, .7, .2, .1,
    .2, .4, .8, 1.5,
    1.5, .2, 2.0, .05,
    1, 1, 1, 1
  ), nrow = S, byrow = TRUE)
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- c(
    SELKERNEL_STEP,
    SELKERNEL_NORMAL,
    SELKERNEL_STEP,
    SELKERNEL_NORMAL
  )
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  reference_density <- function(extrapolate) {

    sei_mat     <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
    total_sd    <- sqrt(tau^2 + sei_mat^2)
    densities   <- matrix(NA_real_, nrow = S, ncol = length(z_values))
    const_rows  <- if (extrapolate) which(!selection[["use_normal"]]) else integer(0)
    weight_rows <- if (extrapolate) integer(0) else which(!selection[["use_normal"]])

    log_norm <- NULL
    if (length(const_rows) > 0L) {
      log_norm <- .selection_step_log_norm_matrix(
        mean              = mu[const_rows, , drop = FALSE],
        sd                = total_sd[const_rows, , drop = FALSE],
        sei               = sei,
        selection_context = BayesTools::selection_context_subset_rows(selection, const_rows)
      )
    }

    for (ell in seq_along(z_values)) {
      y_seq <- z_values[ell] * sei
      dens  <- stats::dnorm(
        matrix(y_seq, nrow = S, ncol = K, byrow = TRUE),
        mean = mu,
        sd   = total_sd
      ) * sei_mat

      if (length(const_rows) > 0L) {
        dens[const_rows, ] <- dens[const_rows, , drop = FALSE] / exp(log_norm)
      }
      if (length(weight_rows) > 0L) {
        selection_weight <- BayesTools::selection_context_subset_rows(
          selection,
          weight_rows
        )
        log_norm_weight <- .selection_step_log_norm_matrix(
          mean              = mu[weight_rows, , drop = FALSE],
          sd                = total_sd[weight_rows, , drop = FALSE],
          sei               = sei,
          selection_context = selection_weight
        )
        bin <- .selection_step_bin_from_z(
          selection[["sign"]] * z_values[ell],
          selection[["p_cuts"]]
        )
        log_pdf <- stats::dnorm(
          matrix(y_seq, nrow = length(weight_rows), ncol = K, byrow = TRUE),
          mean = mu[weight_rows, , drop = FALSE],
          sd   = total_sd[weight_rows, , drop = FALSE],
          log  = TRUE
        ) - log_norm_weight + log(selection[["omega"]][weight_rows, bin])
        dens[weight_rows, ] <- exp(log_pdf) * sei_mat[weight_rows, , drop = FALSE]
      }

      densities[, ell] <- rowMeans(dens)
    }

    return(densities)
  }

  for (extrapolate in c(FALSE, TRUE)) {
    expect_equal(
      .zplot_density_vectorized(
        z_sequence       = z_values,
        mu_samples       = mu,
        tau_within       = tau,
        sei              = sei,
        selection        = selection,
        extrapolate      = extrapolate,
        effect_direction = "positive"
      ),
      reference_density(extrapolate),
      tolerance = 1e-12
    )
  }
})


test_that("zplot threshold separates fitted, extrapolated, and missing-study targets", {

  skip_if_not(.has_native_selnorm_kernel())

  S         <- 4L
  K         <- 3L
  sei       <- c(.07, .13, .29)
  mu        <- matrix(c(
    -.20, .05, .30,
     .10, .12, .18,
     .40, .00, -.15,
     .05, -.25, .20
  ), nrow = S, byrow = TRUE)
  tau       <- matrix(c(
    .05, .07, .10,
    .08, .06, .11,
    .03, .09, .12,
    .04, .10, .15
  ), nrow = S, byrow = TRUE)
  threshold <- stats::qnorm(.025, lower.tail = FALSE)
  spec      <- .test_step_spec(yi = c(.10, .25, .40), sei = sei)

  selection <- spec
  selection[["omega"]]       <- matrix(c(
    1, .7, .2, .1,
    .2, .4, .8, 1.5,
    1.5, .2, 2.0, .05,
    1, 1, 1, 1
  ), nrow = S, byrow = TRUE)
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- c(
    SELKERNEL_STEP,
    SELKERNEL_NORMAL,
    SELKERNEL_STEP,
    SELKERNEL_NORMAL
  )
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau^2 + sei_mat^2)
  q_upper  <- threshold * sei
  q_lower  <- -threshold * sei

  unselected_thresholds <- stats::pnorm(
    matrix(q_upper, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu,
    sd         = total_sd,
    lower.tail = FALSE
  ) + stats::pnorm(
    matrix(q_lower, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu,
    sd         = total_sd,
    lower.tail = TRUE
  )

  step_rows      <- which(!selection[["use_normal"]])
  selection_step <- BayesTools::selection_context_subset_rows(selection, step_rows)
  log_norm       <- .selection_step_log_norm_matrix(
    mean              = mu[step_rows, , drop = FALSE],
    sd                = total_sd[step_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_step
  )
  expected_weights <- matrix(1, nrow = S, ncol = K)
  expected_weights[step_rows, ] <- exp(-log_norm)

  selected_thresholds <- unselected_thresholds
  selected_thresholds[step_rows, ] <- .selection_step_cdf_matrix(
    q                 = q_upper,
    mean              = mu[step_rows, , drop = FALSE],
    sd                = total_sd[step_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_step,
    lower.tail        = FALSE
  ) + .selection_step_cdf_matrix(
    q                 = q_lower,
    mean              = mu[step_rows, , drop = FALSE],
    sd                = total_sd[step_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_step,
    lower.tail        = TRUE
  )

  extrapolated <- .zplot_threshold_vectorized(
    z_threshold      = threshold,
    mu_samples       = mu,
    tau_within       = tau,
    sei              = sei,
    selection        = selection,
    extrapolate      = TRUE,
    effect_direction = "positive"
  )
  fitted <- .zplot_threshold_vectorized(
    z_threshold      = threshold,
    mu_samples       = mu,
    tau_within       = tau,
    sei              = sei,
    selection        = selection,
    extrapolate      = FALSE,
    effect_direction = "positive"
  )

  expect_equal(
    extrapolated[["EDR"]],
    rowSums(unselected_thresholds * expected_weights) / rowSums(expected_weights),
    tolerance = 1e-12
  )
  expect_equal(extrapolated[["weights"]], rowMeans(expected_weights), tolerance = 1e-12)
  expect_equal(fitted[["EDR"]], rowMeans(selected_thresholds), tolerance = 1e-12)
  expect_equal(fitted[["weights"]], rowMeans(expected_weights), tolerance = 1e-12)
})


test_that("native regplot selection intervals match R reference", {

  skip_if_not(.has_native_regplot_selection_mixture())

  set.seed(131)
  S            <- 15L
  K            <- 4L
  se           <- .14
  probs        <- c(.025, .975)
  mean_samples <- matrix(rnorm(S * K, 0, .35), nrow = S, ncol = K)
  sd_samples   <- matrix(runif(S * K, .02, .30), nrow = S, ncol = K)
  sd_samples[c(3, 11)] <- sqrt(.Machine$double.eps) / 10
  spec         <- .test_step_spec(yi = seq_len(K) / 10, sei = rep(se, K))

  selection <- spec
  selection[["omega"]]       <- matrix(
    runif(S * spec[["n_bins"]], .02, 3),
    nrow = S,
    ncol = spec[["n_bins"]]
  )
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection[["kernel_mode"]][seq(2, S, by = 4)] <- SELKERNEL_NORMAL
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  setup <- list(
    is_weightfunction = !selection[["use_normal"]],
    selection         = selection
  )

  native <- .regplot_selnorm_mixture_interval_quantiles(
    mean_samples      = mean_samples,
    sd_samples        = sd_samples,
    se                = se,
    probs             = probs,
    selection_context = selection
  )
  ref <- .regplot_selection_mixture_interval_quantiles_r(
    mean_samples     = mean_samples,
    sd_samples       = sd_samples,
    se               = se,
    probs            = probs,
    setup            = setup,
    effect_direction = "positive"
  )

  expect_equal(native[["lower"]], ref[["lower"]], tolerance = 1e-5)
  expect_equal(native[["upper"]], ref[["upper"]], tolerance = 1e-5)
})

test_that("plot quantiles retain tiny positive normal scales", {

  tiny_sd <- sqrt(.Machine$double.eps) / 100
  probs   <- c(.025, .975)

  regplot_r <- .regplot_mixture_interval_quantiles_r(
    mean_samples = matrix(0, nrow = 1L),
    sd_samples   = matrix(tiny_sd, nrow = 1L),
    probs        = probs
  )
  expect_equal(
    unlist(regplot_r, use.names = FALSE) / tiny_sd,
    stats::qnorm(probs),
    tolerance = 1e-7
  )

  if (.has_native_regplot_mixture()) {
    regplot_native <- .regplot_mixture_interval_quantiles(
      mean_samples = matrix(0, nrow = 1L),
      sd_samples   = matrix(tiny_sd, nrow = 1L),
      probs        = probs
    )
    expect_equal(
      unlist(regplot_native, use.names = FALSE) / tiny_sd,
      stats::qnorm(probs),
      tolerance = 1e-7
    )
  }

  funnel_setup <- list(
    mu                = 0,
    tau               = tiny_sd,
    PET               = 0,
    PEESE             = 0,
    is_weightfunction = FALSE,
    selection         = NULL
  )
  expect_equal(
    vapply(
      probs,
      .funnel_model_averaged_quantile,
      numeric(1L),
      se               = 0,
      setup            = funnel_setup,
      effect_direction = "positive"
    ) / tiny_sd,
    stats::qnorm(probs),
    tolerance = 1e-7
  )
})

test_that("native weighted selected-normal summary matches matrix reductions", {

  skip_if_not(.has_native_selnorm_kernel())

  set.seed(141)
  S    <- 17L
  K    <- 5L
  yi   <- -seq_len(K) / 20
  sei  <- seq(.06, .21, length.out = K)
  spec <- .test_step_spec(yi = yi, sei = sei, effect_direction = "negative")
  mean <- matrix(rnorm(S * K, 0, .25), nrow = S, ncol = K)
  sd   <- matrix(runif(S * K, .04, .35), nrow = S, ncol = K)
  psis <- matrix(rexp(S * K), nrow = S, ncol = K)
  psis <- sweep(psis, 2, colSums(psis), "/")

  selection <- spec
  selection[["omega"]]       <- matrix(
    runif(S * spec[["n_bins"]], .03, 2),
    nrow = S,
    ncol = spec[["n_bins"]]
  )
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection[["kernel_mode"]][seq(3, S, by = 5)] <- SELKERNEL_NORMAL
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  native <- .selection_step_weighted_summary(
    yi                = yi,
    mean              = mean,
    sd                = sd,
    sei               = sei,
    psis_weights      = psis,
    selection_context = selection
  )
  cdf_lower <- .selection_step_cdf_matrix(
    yi, mean, sd, sei, selection, lower.tail = TRUE
  )
  cdf_upper <- .selection_step_cdf_matrix(
    yi, mean, sd, sei, selection, lower.tail = FALSE
  )
  moments <- .selection_step_moments_matrix(mean, sd, sei, selection)

  expect_equal(
    native[["log_lower"]],
    log(colSums(psis * cdf_lower)),
    tolerance = 1e-12
  )
  expect_equal(
    native[["log_upper"]],
    log(colSums(psis * cdf_upper)),
    tolerance = 1e-12
  )
  expect_equal(
    exp(native[["log_lower"]]) + exp(native[["log_upper"]]),
    rep(1, K),
    tolerance = 1e-15
  )
  expected_mean <- colSums(psis * moments[["mean"]])
  centered      <- sweep(moments[["mean"]], 2L, expected_mean, "-")
  expected_variance <- colSums(psis * (
    moments[["second"]] - moments[["mean"]]^2 + centered^2
  ))
  expect_equal(native[["mean"]], expected_mean, tolerance = 1e-12)
  expect_equal(native[["variance"]], expected_variance, tolerance = 1e-12)
})


test_that("selected-normal PIT preserves finite affine contrasts", {

  skip_if_not(.has_native_selnorm_kernel())

  q_base    <- 0x1.d5f4b3ac79505p+560
  mean_base <- 0x1.d5f4b3ac79506p+560
  sd        <- 0x1p+508
  sei       <- 0x1.6ec9d2937021cp-463
  expect_identical(q_base / sei, mean_base / sei)

  for (effect_direction in c("positive", "negative")) {
    direction <- if (effect_direction == "positive") 1 else -1
    yi        <- direction * q_base
    mean      <- direction * mean_base
    expected  <- (yi - mean) / sd
    spec      <- .test_step_spec(
      yi               = yi,
      sei              = sei,
      effect_direction = effect_direction
    )
    selection <- spec
    selection[["omega"]]       <- matrix(1, nrow = 1L, ncol = spec[["n_bins"]])
    selection[["alpha"]]       <- 0
    selection[["phack_kind"]]  <- 0L
    selection[["kernel_mode"]] <- SELKERNEL_STEP
    selection[["use_normal"]]  <- FALSE

    summary <- .selection_step_weighted_summary(
      yi                = yi,
      mean              = matrix(mean, nrow = 1L),
      sd                = matrix(sd, nrow = 1L),
      sei               = sei,
      psis_weights      = matrix(7, nrow = 1L),
      selection_context = selection
    )

    expect_equal(summary[["log_lower"]], stats::pnorm(
      expected,
      log.p = TRUE
    ))
    expect_equal(summary[["log_upper"]], stats::pnorm(
      expected,
      lower.tail = FALSE,
      log.p      = TRUE
    ))
    expect_equal(
      .loo_pit_z_from_log_probabilities(
        summary[["log_lower"]],
        summary[["log_upper"]]
      ),
      expected
    )
  }
})


test_that("selected-normal CDF handles extreme finite scales and rejects invalid ones", {

  skip_if_not(.has_native_selnorm_kernel())

  spec  <- .test_step_spec(0, 1)
  omega <- matrix(1, nrow = 1L, ncol = spec[["n_bins"]])
  tiny  <- .Machine$double.xmin * .Machine$double.eps
  cdf   <- .selnorm_kernel_cdf_matrix(
    q              = 4,
    mean           = matrix(4, nrow = 1L),
    sd             = matrix(1, nrow = 1L),
    sei            = tiny,
    omega          = omega,
    selection_spec = spec
  )
  log_norm <- .selnorm_kernel_log_norm_matrix(
    mean           = matrix(4, nrow = 1L),
    sd             = matrix(1, nrow = 1L),
    sei            = tiny,
    omega          = omega,
    selection_spec = spec
  )

  expect_identical(cdf[1L, 1L], .5)
  expect_identical(log_norm[1L, 1L], 0)

  evaluate <- function(q = 0, sd = 1, sei = 1) {
    .selnorm_kernel_cdf_matrix(
      q              = q,
      mean           = matrix(0, nrow = 1L),
      sd             = matrix(sd, nrow = 1L),
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )[1L, 1L]
  }
  expect_true(is.nan(evaluate(q = NaN)))
  expect_true(is.nan(evaluate(sd = Inf)))
  expect_true(is.nan(evaluate(sei = Inf)))
})


test_that("normal scale composition avoids intermediate overflow", {

  component <- .Machine$double.xmax / 2
  expected  <- component * sqrt(2)
  tiny      <- .Machine$double.xmin / 2

  expect_true(is.infinite(sqrt(component^2 + component^2)))
  expect_equal(.root_sum_squares(component, component), expected)
  expect_identical(sqrt(tiny^2 + tiny^2), 0)
  expect_gt(.root_sum_squares(tiny, tiny), 0)
  expect_identical(.root_sum_squares(0, 0), 0)
  expect_identical(.root_sum_squares(Inf, 1), Inf)
  expect_true(is.na(.root_sum_squares(NA_real_, 1)))

  matrix_result <- .root_sum_squares(
    matrix(c(3, 5), nrow = 1L),
    matrix(c(4, 12), nrow = 1L)
  )
  expect_identical(matrix_result, matrix(c(5, 13), nrow = 1L))
})


test_that("zero-weight selection intervals produce exact CDF plateaus", {

  skip_if_not(.has_native_selnorm_kernel())

  spec        <- .test_step_spec(0, 1)
  lower       <- spec[["z_lower"]][2L]
  upper       <- spec[["z_upper"]][2L]
  query       <- lower + c(1 / 3, 2 / 3) * (upper - lower)
  omega       <- matrix(c(1, 0, .35, .2), nrow = 1L)
  mean        <- matrix(0, nrow = 1L, ncol = 2L)
  sd          <- matrix(1, nrow = 1L, ncol = 2L)
  lower_tail  <- .selnorm_kernel_cdf_matrix(
    q              = query,
    mean           = mean,
    sd             = sd,
    sei            = rep(1, 2L),
    omega          = omega,
    selection_spec = spec
  )
  upper_tail <- .selnorm_kernel_cdf_matrix(
    q              = query,
    mean           = mean,
    sd             = sd,
    sei            = rep(1, 2L),
    omega          = omega,
    selection_spec = spec,
    lower.tail     = FALSE
  )

  expect_identical(lower_tail[1L, 1L], lower_tail[1L, 2L])
  expect_identical(upper_tail[1L, 1L], upper_tail[1L, 2L])
  expect_equal(lower_tail + upper_tail, matrix(1, nrow = 1L, ncol = 2L))
})


test_that("zero-weight selected-normal rows do not enter mixture evaluation", {

  skip_if_not(.has_native_selnorm_kernel())

  spec      <- .test_step_spec(yi = 0, sei = 1)
  selection <- spec
  selection[["omega"]]       <- matrix(1, nrow = 2L, ncol = spec[["n_bins"]])
  selection[["alpha"]]       <- c(0, 0)
  selection[["phack_kind"]]  <- c(0L, 0L)
  selection[["kernel_mode"]] <- c(SELKERNEL_STEP, SELKERNEL_STEP)
  selection[["use_normal"]]  <- c(FALSE, FALSE)

  summary <- .selection_step_weighted_summary(
    yi                = 0,
    mean              = matrix(c(0, 0), ncol = 1L),
    sd                = matrix(c(1, 0), ncol = 1L),
    sei               = 1,
    psis_weights      = matrix(c(1, 0), ncol = 1L),
    selection_context = selection
  )

  expect_equal(summary[["log_lower"]], stats::pnorm(0, log.p = TRUE))
  expect_equal(summary[["log_upper"]], stats::pnorm(
    0,
    lower.tail = FALSE,
    log.p      = TRUE
  ))
  expect_equal(summary[["mean"]], 0)
  expect_equal(summary[["variance"]], 1)
})


test_that("selected-normal variance is centered before reducing moments", {

  skip_if_not(.has_native_selnorm_kernel())

  spec      <- .test_step_spec(yi = 1e12, sei = 1)
  selection <- spec
  selection[["omega"]]       <- matrix(1, nrow = 1L, ncol = spec[["n_bins"]])
  selection[["alpha"]]       <- 0
  selection[["phack_kind"]]  <- 0L
  selection[["kernel_mode"]] <- SELKERNEL_STEP
  selection[["use_normal"]]  <- FALSE

  summary <- .selection_step_weighted_summary(
    yi                = 1e12,
    mean              = matrix(1e12, nrow = 1L),
    sd                = matrix(1, nrow = 1L),
    sei               = 1,
    psis_weights      = matrix(1, nrow = 1L),
    selection_context = selection
  )

  expect_equal(summary[["mean"]], 1e12)
  expect_equal(summary[["variance"]], 1)
})


test_that("selected-normal log tails survive probability underflow", {

  skip_if_not(.has_native_selnorm_kernel())

  yi        <- c(-100, 100)
  sei       <- c(1, 1)
  S         <- 2L
  spec      <- .test_step_spec(yi = yi, sei = sei)
  selection <- spec
  selection[["omega"]]       <- matrix(1, nrow = S, ncol = spec[["n_bins"]])
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection[["use_normal"]]  <- rep(FALSE, S)

  summary <- .selection_step_weighted_summary(
    yi                = yi,
    mean              = matrix(0, nrow = S, ncol = length(yi)),
    sd                = matrix(1, nrow = S, ncol = length(yi)),
    sei               = sei,
    psis_weights      = matrix(0.5, nrow = S, ncol = length(yi)),
    selection_context = selection
  )

  expect_true(all(is.finite(summary[["log_lower"]])))
  expect_true(all(is.finite(summary[["log_upper"]])))
  expect_equal(
    .loo_pit_z_from_log_probabilities(
      summary[["log_lower"]],
      summary[["log_upper"]]
    ),
    yi,
    tolerance = 1e-10
  )
})


test_that("p-hacking priors are deferred before RoBMA model construction", {

  phacking <- BayesTools::prior_phacking(form = "linear")
  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  expect_false("prior_phacking" %in% getNamespaceExports("RoBMA"))
  expect_false("prior_bias" %in% getNamespaceExports("RoBMA"))
  expect_error(
    .use_plot_prior_list_dispatch(phacking),
    "not enabled"
  )
  expect_error(
    bselmodel(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_bias                = phacking,
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    bselmodel(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_bias                = BayesTools::prior_bias(
        selection = selection,
        phacking  = phacking
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    RoBMA(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_bias                = phacking,
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    RoBMA(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_bias                = BayesTools::prior_bias(
        selection = selection,
        phacking  = phacking
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    .selection_spec(
      priors           = list(outcome = list(bias = phacking)),
      yi               = c(.1, .2),
      sei              = c(.1, .1),
      effect_direction = "positive"
    ),
    "not enabled"
  )
  expect_error(
    .selection_spec(
      priors = list(outcome = list(
        bias = BayesTools::prior_bias(selection = selection, phacking = phacking)
      )),
      yi               = c(.1, .2),
      sei              = c(.1, .1),
      effect_direction = "positive"
    ),
    "not enabled"
  )
})
