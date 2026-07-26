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
