context("Selection kernel")
skip_on_cran()

test_that("selected-normal native boundary handles structural zero omega", {

  skip_if_not(.has_native_selnorm_kernel())

  yi     <- c(.10, .20)
  sei    <- c(.10, .10)
  spec   <- .test_step_spec(yi, sei)
  mu     <- matrix(0, nrow = 1, ncol = 2)
  sigma  <- matrix(.20, nrow = 1, ncol = 2)
  omega_empty_bin <- matrix(c(1, 0, .5, .25), nrow = 1)
  omega_obs_bin   <- matrix(c(1, .5, 0, .25), nrow = 1)

  empty_bin_log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega_empty_bin,
    selection_spec = spec
  )
  obs_bin_log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega_obs_bin,
    selection_spec = spec
  )

  expect_true(all(is.finite(empty_bin_log_lik)))
  expect_equal(obs_bin_log_lik[1, 1], -Inf)
  expect_true(is.finite(obs_bin_log_lik[1, 2]))
})

test_that("selected-normal native row arguments reject incompatible lengths", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- c(.10, .20)
  sei   <- c(.10, .10)
  spec  <- .test_step_spec(yi, sei)
  mean  <- matrix(0, nrow = 2, ncol = 2)
  sd    <- matrix(.20, nrow = 2, ncol = 2)
  omega <- matrix(1, nrow = 2, ncol = spec[["n_bins"]])

  expect_error(
    .selnorm_kernel_loglik_matrix(
      yi             = yi,
      mu_num         = mean,
      sigma_num      = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = c(0, 0, 0)
    ),
    "length 1 or"
  )
  expect_error(
    .selnorm_kernel_loglik_matrix(
      yi             = yi,
      mu_num         = mean,
      sigma_num      = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      phack_kind     = c(0L, 0L, 0L)
    ),
    "length 1 or"
  )
  expect_error(
    .selnorm_kernel_loglik_matrix(
      yi             = yi,
      mu_num         = mean,
      sigma_num      = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      kernel_mode    = rep(SELKERNEL_STEP, 3L)
    ),
    "length 1 or"
  )
})

test_that("active p-hacking kernels error on unsupported non-likelihood paths", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- .10
  sei   <- .10
  spec  <- .test_step_spec(yi, sei)
  mean  <- matrix(c(0, .02), nrow = 2, ncol = 1)
  sd    <- matrix(.20, nrow = 2, ncol = 1)
  omega <- matrix(1, nrow = 2, ncol = spec[["n_bins"]])
  alpha <- c(.50, 0)
  phack <- c(1L, 0L)
  mode  <- c(SELKERNEL_STEP_PHACK_POWER, SELKERNEL_STEP)

  log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mean,
    sigma_num      = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec,
    alpha          = alpha,
    phack_kind     = phack,
    kernel_mode    = mode
  )
  expect_true(all(is.finite(log_lik)))

  expect_error(
    .selnorm_kernel_cdf_matrix(
      q              = yi,
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )
  expect_error(
    .selnorm_kernel_moments_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )
  expect_error(
    .selnorm_kernel_rng_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )
  expect_error(
    .selnorm_kernel_weighted_summary(
      yi             = yi,
      mean           = mean,
      sd             = sd,
      sei            = sei,
      psis_weights   = matrix(c(.5, .5), nrow = 2, ncol = 1),
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )

  context <- spec
  context[["omega"]]          <- omega
  context[["alpha"]]          <- alpha
  context[["phack_kind"]]     <- phack
  context[["kernel_mode"]]    <- mode
  context[["bias_indicator"]] <- c(1L, 1L)
  context[["use_normal"]]     <- c(FALSE, FALSE)

  setup <- list(
    mu        = c(0, .02),
    selection = context
  )

  expect_error(
    .funnel_selected_cdf(
      q                = yi,
      rows             = 1:2,
      se               = 0,
      total_sd         = c(.20, .20),
      setup            = setup,
      effect_direction = "positive"
    ),
    "active p-hacking"
  )
})

test_that("funnel selected-normal CDF is finite at zero standard error", {

  skip_if_not(.has_native_selnorm_kernel())

  spec <- .test_step_spec(yi = c(.10, .25), sei = c(.10, .15))

  selection_context <- spec
  selection_context[["omega"]]          <- matrix(c(1, .70, .35, .20), nrow = 1)
  selection_context[["alpha"]]          <- 0
  selection_context[["phack_kind"]]     <- 0L
  selection_context[["kernel_mode"]]    <- SELKERNEL_STEP
  selection_context[["bias_indicator"]] <- 1L
  selection_context[["use_normal"]]     <- FALSE

  setup <- list(
    mu                = 0,
    tau               = .20,
    PET               = 0,
    PEESE             = 0,
    is_weightfunction = TRUE,
    selection         = selection_context
  )

  out <- .funnel_model_averaged_cdf(
    q                = 0,
    se               = 0,
    setup            = setup,
    effect_direction = "positive"
  )

  expected <- (.20 * stats::pnorm(0, mean = 0, sd = .20)) /
    (stats::pnorm(0, mean = 0, sd = .20, lower.tail = FALSE) +
      .20 * stats::pnorm(0, mean = 0, sd = .20))

  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("funnel and regplot selected CDF paths match non-telescoped fallback", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_regplot_selection_mixture())

  set.seed(109)
  S     <- 11L
  K     <- 3L
  yi    <- c(.08, .20, .42)
  sei   <- c(.09, .14, .22)
  spec  <- .test_step_spec(yi, sei)
  omega <- matrix(runif(S * spec[["n_bins"]], .10, 1.80), nrow = S)

  selection <- spec
  selection[["omega"]]       <- omega
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection[["use_normal"]]  <- rep(FALSE, S)
  selection[["has_phack"]]   <- FALSE

  selection_fallback <- selection
  selection_fallback[["telescope_probabilities"]] <- FALSE

  setup <- list(
    mu                = seq(-.12, .20, length.out = S),
    tau               = seq(.03, .18, length.out = S),
    PET               = rep(0, S),
    PEESE             = rep(0, S),
    is_weightfunction = rep(TRUE, S),
    selection         = selection
  )
  setup_fallback <- setup
  setup_fallback[["selection"]] <- selection_fallback

  expect_equal(
    .funnel_model_averaged_cdf(.14, .13, setup, "positive"),
    .funnel_model_averaged_cdf(.14, .13, setup_fallback, "positive"),
    tolerance = 1e-12
  )

  mean_samples <- matrix(rnorm(S * K, .05, .18), nrow = S)
  sd_samples   <- matrix(runif(S * K, .08, .35), nrow = S)

  expect_equal(
    .regplot_selnorm_mixture_interval_quantiles(
      mean_samples      = mean_samples,
      sd_samples        = sd_samples,
      se                = .13,
      probs             = c(.10, .90),
      selection_context = selection
    ),
    .regplot_selnorm_mixture_interval_quantiles(
      mean_samples      = mean_samples,
      sd_samples        = sd_samples,
      se                = .13,
      probs             = c(.10, .90),
      selection_context = selection_fallback
    ),
    tolerance = 1e-8
  )
})

test_that("native funnel model-averaged quantiles match R fallback", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_funnel_model_averaged_quantiles(list(
    selection = list()
  )))

  set.seed(127)
  S     <- 19L
  yi    <- c(.08, .20, .42)
  sei   <- c(.09, .14, .22)
  omega <- matrix(runif(S * 4L, .10, 1.80), nrow = S)

  for (effect_direction in c("positive", "negative")) {
    sign <- if (effect_direction == "positive") 1 else -1
    spec <- .test_step_spec(sign * yi, sei, effect_direction = effect_direction)

    for (telescope_probabilities in c(TRUE, FALSE)) {
      selection <- spec
      selection[["omega"]] <- omega[
        , seq_len(spec[["n_bins"]]), drop = FALSE
      ]
      selection[["alpha"]]       <- rep(0, S)
      selection[["phack_kind"]]  <- rep(0L, S)
      selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
      selection[["kernel_mode"]][seq(3, S, by = 5)] <- SELKERNEL_NORMAL
      selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL
      selection[["has_phack"]]  <- FALSE
      selection[["telescope_probabilities"]] <- telescope_probabilities

      setup <- list(
        mu                = sign * seq(-.12, .20, length.out = S),
        tau               = seq(.00, .18, length.out = S),
        PET               = rnorm(S, 0, .03),
        PEESE             = rnorm(S, 0, .02),
        is_weightfunction = !selection[["use_normal"]],
        selection         = selection
      )
      se_sequence <- c(0, .05, .13, .27)

      native <- .funnel_model_averaged_quantiles_native(
        se_sequence      = se_sequence,
        setup            = setup,
        effect_direction = effect_direction
      )
      ref <- list(
        lower = vapply(
          se_sequence,
          .funnel_model_averaged_quantile,
          numeric(1),
          p                = .025,
          setup            = setup,
          effect_direction = effect_direction
        ),
        upper = vapply(
          se_sequence,
          .funnel_model_averaged_quantile,
          numeric(1),
          p                = .975,
          setup            = setup,
          effect_direction = effect_direction
        ),
        mid = vapply(
          se_sequence,
          .funnel_model_averaged_quantile,
          numeric(1),
          p                = .5,
          setup            = setup,
          effect_direction = effect_direction
        )
      )

      expect_equal(native[["lower"]], ref[["lower"]], tolerance = 1e-5)
      expect_equal(native[["upper"]], ref[["upper"]], tolerance = 1e-5)
      expect_equal(native[["mid"]],   ref[["mid"]],   tolerance = 1e-5)
    }
  }
})

test_that("native funnel model-averaged quantiles reject active p-hacking", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_funnel_model_averaged_quantiles(list(
    selection = list()
  )))

  S    <- 3L
  spec <- .test_step_spec(c(.08, .20), c(.09, .14))

  selection <- spec
  selection[["omega"]]       <- matrix(
    runif(S * spec[["n_bins"]], .10, 1.80),
    nrow = S
  )
  selection[["alpha"]]       <- rep(.20, S)
  selection[["phack_kind"]]  <- rep(1L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP_PHACK_POWER, S)
  selection[["use_normal"]] <- rep(FALSE, S)
  selection[["has_phack"]]  <- TRUE

  setup <- list(
    mu                = seq(-.12, .20, length.out = S),
    tau               = seq(.00, .18, length.out = S),
    PET               = rep(0, S),
    PEESE             = rep(0, S),
    is_weightfunction = rep(TRUE, S),
    selection         = selection
  )

  expect_error(
    .funnel_model_averaged_quantiles_native(
      se_sequence      = c(.05, .13),
      setup            = setup,
      effect_direction = "positive"
    ),
    "active p-hacking"
  )
})

test_that("funnel max_samples uses a global nested simple random sample", {

  bias_indicator <- c(rep(1L, 5000), rep(2L, 3000), rep(3L, 2000))

  expect_equal(.normalize_max_samples(10), 10L)
  expect_equal(.normalize_max_samples(Inf), Inf)
  expect_error(.normalize_max_samples(9), "at least 10")
  expect_error(.normalize_max_samples(10.5), "positive integer or Inf")
  expect_equal(.normalize_funnel_max_samples(10), 10L)
  expect_equal(.normalize_funnel_max_samples(Inf), Inf)
  expect_error(.normalize_funnel_max_samples(9), "at least 10")
  expect_error(.normalize_funnel_max_samples(10.5), "positive integer or Inf")

  plain_rows <- .thin_sample_rows(10000, 1000)
  expect_equal(plain_rows, round(seq(from = 1, to = 10000, length.out = 1000)))
  expect_true(max(diff(plain_rows)) - min(diff(plain_rows)) <= 1)
  expect_null(.thin_sample_rows(1000, 1000))
  expect_null(.thin_sample_rows(100, Inf))

  rows <- .nested_srs_rows(seq_along(bias_indicator), 1000)
  expect_equal(rows, .nested_srs_rows(seq_along(bias_indicator), 1000))
  expect_equal(length(rows), 1000)
  expect_true(all(
    rows %in% .nested_srs_rows(seq_along(bias_indicator), 2000)
  ))
  expect_equal(
    as.integer(table(factor(bias_indicator[rows], levels = 1:3))),
    c(485L, 323L, 192L)
  )
  expect_equal(
    .nested_srs_rows(seq_along(bias_indicator), Inf),
    seq_along(bias_indicator)
  )
})
