context("Selection probability numerics")

source(testthat::test_path("common-functions.R"))

skip_on_cran()


test_that("native zero-SE funnel retains remote-tail interval mass", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_funnel_model_averaged_quantiles(list(
    selection = list()
  )))

  spec      <- .test_step_spec(c(.10, .25), c(.10, .15))
  selection <- spec
  selection[["omega"]]                   <- matrix(c(1, 0, 0, 0), nrow = 1)
  selection[["alpha"]]                   <- 0
  selection[["phack_kind"]]              <- 0L
  selection[["kernel_mode"]]             <- SELKERNEL_STEP
  selection[["use_normal"]]              <- FALSE
  selection[["has_phack"]]               <- FALSE
  selection[["telescope_probabilities"]] <- FALSE

  setup <- list(
    mu                = -10,
    tau               = 1,
    PET               = 0,
    PEESE             = 0,
    is_weightfunction = TRUE,
    selection         = selection
  )
  actual <- .funnel_model_averaged_quantiles_native(
    se_sequence      = 0,
    setup            = setup,
    effect_direction = "positive"
  )

  probabilities      <- c(.025, .975, .5)
  log_survival_zero  <- stats::pnorm(
    0,
    mean       = -10,
    sd         = 1,
    lower.tail = FALSE,
    log.p      = TRUE
  )
  expected           <- stats::qnorm(
    log_survival_zero + log1p(-probabilities),
    mean       = -10,
    sd         = 1,
    lower.tail = FALSE,
    log.p      = TRUE
  )

  expect_equal(
    unname(unlist(actual, use.names = FALSE)),
    expected,
    tolerance = 1e-6
  )
})


test_that("native funnel CDF uses structurally matched tail masses", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_funnel_model_averaged_quantiles(list(
    selection = list()
  )))

  make_setup <- function(mu, tau, lower_weight) {

    prior <- BayesTools::prior_weightfunction(
      side    = "one-sided",
      steps   = .025,
      weights = BayesTools::wf_fixed(c(1, lower_weight))
    )
    selection <- .selection_spec(
      priors           = list(outcome = list(bias = prior)),
      yi               = .20,
      sei              = .10,
      effect_direction = "positive",
      signed_data      = FALSE
    )
    selection[["omega"]]                   <- matrix(c(1, lower_weight), nrow = 1L)
    selection[["alpha"]]                   <- 0
    selection[["phack_kind"]]              <- 0L
    selection[["kernel_mode"]]             <- SELKERNEL_STEP
    selection[["use_normal"]]              <- FALSE
    selection[["has_phack"]]               <- FALSE
    selection[["telescope_probabilities"]] <- TRUE

    return(list(
      mu                = mu,
      tau               = tau,
      PET               = 0,
      PEESE             = 0,
      is_weightfunction = TRUE,
      selection         = selection
    ))
  }

  selected_quantile <- function(p, mu, tau, se, lower_weight) {

    sd       <- sqrt(tau^2 + se^2)
    cut      <- stats::qnorm(.975) * se
    cut_prob <- stats::pnorm(cut, mean = mu, sd = sd)
    total    <- lower_weight * cut_prob + 1 - cut_prob
    target   <- p * total

    if (target <= lower_weight * cut_prob) {
      return(stats::qnorm(target / lower_weight, mean = mu, sd = sd))
    }

    return(stats::qnorm(
      cut_prob + target - lower_weight * cut_prob,
      mean = mu,
      sd   = sd
    ))
  }

  cases <- list(
    list(
      mu           = 0.077382373995749,
      tau          = 0.360189673816798,
      se           = 0,
      lower_weight = 0.500377518771078
    ),
    list(
      mu           = 0.0182293024488417,
      tau          = 0.313440158972955,
      se           = 0.00808080808080808,
      lower_weight = 0.477068312576876
    )
  )

  for (case in cases) {
    setup <- make_setup(
      mu           = case[["mu"]],
      tau          = case[["tau"]],
      lower_weight = case[["lower_weight"]]
    )
    actual <- .funnel_model_averaged_quantiles_native(
      se_sequence      = case[["se"]],
      setup            = setup,
      effect_direction = "positive"
    )
    expected <- vapply(c(.025, .975, .5), function(p) {
      selected_quantile(
        p            = p,
        mu           = case[["mu"]],
        tau          = case[["tau"]],
        se           = case[["se"]],
        lower_weight = case[["lower_weight"]]
      )
    }, numeric(1))

    expect_equal(
      unname(unlist(actual, use.names = FALSE)),
      expected,
      tolerance = 1e-10
    )
  }
})


test_that("native central interval probabilities are not projected", {

  skip_if_not(.has_native_selnorm_kernel())

  native_log_probability <- function(mean, sd, lower, upper) {

    return(.Call(
      "RoBMA_selnorm_kernel_log_norm_matrix",
      matrix(mean, nrow = 1L, ncol = 1L),
      matrix(sd, nrow = 1L, ncol = 1L),
      as.numeric(1),
      matrix(c(0, 1, 0), nrow = 1L),
      as.numeric(0),
      as.integer(0),
      as.integer(SELKERNEL_STEP),
      as.numeric(c(upper, lower, -Inf)),
      as.numeric(c(Inf, upper, lower)),
      as.integer(1),
      as.integer(0),
      as.numeric(c(0, 0)),
      as.numeric(c(0, 0)),
      as.numeric(c(-Inf, Inf)),
      as.integer(1),
      as.integer(0),
      FALSE,
      PACKAGE = "RoBMA"
    ))
  }

  mean  <- .125
  sd    <- .8
  lower <- -.75
  upper <- 1.25

  actual   <- as.numeric(native_log_probability(mean, sd, lower, upper))
  expected <- log(
    stats::pnorm(upper, mean = mean, sd = sd) -
      stats::pnorm(lower, mean = mean, sd = sd)
  )

  expect_identical(actual, expected)
  expect_true(is.nan(native_log_probability(NaN, sd, lower, upper)))
  expect_true(is.nan(native_log_probability(mean, Inf, lower, upper)))
})


test_that("selected-normal likelihood preserves representable extreme scores", {

  skip_if_not(.has_native_selnorm_kernel())

  spec   <- .test_step_spec(0, 1)
  omega  <- matrix(1, nrow = 1L, ncol = spec[["n_bins"]])
  actual <- .selnorm_kernel_loglik_matrix(
    yi             = 1e308,
    mu_num         = matrix(-1e308, nrow = 1L),
    sigma_num      = matrix(1e308, nrow = 1L),
    sei            = 1,
    omega          = omega,
    selection_spec = spec
  )
  expected <- stats::dnorm(2, log = TRUE) - log(1e308)

  expect_equal(actual[1L, 1L], expected, tolerance = 1e-12)
})


test_that("selected-normal fails closed for score-collapsed STEP bins", {

  skip_if_not(.has_native_selnorm_kernel())

  spec   <- .test_step_spec(0, 1)
  mean   <- matrix(1, nrow = 1L)
  sd     <- matrix(1, nrow = 1L)
  sei    <- 1e-310
  omega  <- matrix(c(0, 1, 0, 0), nrow = 1L)
  actual <- .selnorm_kernel_log_norm_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )
  moments <- .selnorm_kernel_moments_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )

  expect_true(is.nan(actual[1L, 1L]))
  expect_true(all(is.nan(unlist(moments))))
})


test_that("selected-normal RNG avoids division-first score coordinates", {

  skip_if_not(.has_native_selnorm_kernel())

  S              <- 32L
  spec           <- .test_step_spec(0, 1)
  omega          <- matrix(rep(c(1, 0, 0, 0), each = S), nrow = S)
  expected_lower <- spec[["z_lower"]][1L] * 1e-310

  set.seed(8041)
  actual <- .selnorm_kernel_rng_matrix(
    mean           = matrix(1, nrow = S),
    sd             = matrix(1, nrow = S),
    sei            = 1e-310,
    omega          = omega,
    selection_spec = spec
  )

  expect_true(all(is.finite(actual)))
  expect_true(all(actual[, 1L] >= expected_lower))
})


test_that("selected-normal threshold and z-density use affine scores", {

  skip_if_not(.has_native_selnorm_kernel())
  skip_if_not(.has_native_zplot_density(selection = TRUE))
  skip_if_not(.has_native_zplot_threshold())

  spec       <- .test_step_spec(0, 1)
  omega      <- matrix(1, nrow = 1L, ncol = spec[["n_bins"]])
  mean       <- matrix(-1e308, nrow = 1L)
  sd         <- matrix(1e308, nrow = 1L)
  sei        <- 1e308
  expected_p <- stats::pnorm(-2, mean = -1, sd = 1) +
    stats::pnorm(2, mean = -1, sd = 1, lower.tail = FALSE)
  expected_d <- stats::dnorm(2, mean = -1, sd = 1)

  for (kernel_mode in c(SELKERNEL_NORMAL, SELKERNEL_STEP)) {
    selection <- spec
    selection[["omega"]]       <- omega
    selection[["alpha"]]       <- 0
    selection[["phack_kind"]]  <- 0L
    selection[["kernel_mode"]] <- kernel_mode

    threshold <- .zplot_selnorm_threshold_summary(
      z_threshold       = 2,
      mean              = mean,
      sd                = sd,
      sei               = sei,
      selection_context = selection,
      extrapolate       = FALSE
    )
    density <- .zplot_selnorm_density_matrix(
      z_sequence        = 2,
      mean              = mean,
      sd                = sd,
      sei               = sei,
      selection_context = selection,
      extrapolate       = FALSE
    )

    expect_equal(threshold[["EDR"]], expected_p, tolerance = 1e-12)
    expect_equal(density[1L, 1L], expected_d, tolerance = 1e-12)
  }
})
