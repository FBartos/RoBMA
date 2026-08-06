context("Options and brma_samples metadata")

test_that("RoBMA options expose only public options", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  expect_true("silent" %in% names(old_options))
  expect_false("RoBMA_version" %in% names(old_options))
  expect_false("module_location" %in% names(old_options))

  updated <- RoBMA.options(silent = TRUE)
  expect_true(updated[["silent"]])
  expect_true(RoBMA.get_option("silent"))

  expect_error(RoBMA.options(RoBMA_version = "bad"), "Unmatched or ambiguous option")
  expect_error(RoBMA.get_option("RoBMA_version"), "Unmatched or ambiguous option")
  expect_error(RoBMA.options(TRUE), "All options must be named")
  expect_error(RoBMA.options(silent = 1), "must be TRUE or FALSE")
  expect_error(RoBMA.options(max_cores = 0), "must be an integer")
  expect_error(RoBMA.options(default_bias_PET.scale = expression(stop("boom"))), "finite number")

  updated <- RoBMA.options(max_cores = 1L)
  expect_equal(updated[["max_cores"]], 1L)
})

test_that("fitting constructors inherit silent option when omitted", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  norm_args <- list(
    yi        = c(0.1, 0.2, 0.3),
    sei       = c(0.2, 0.2, 0.2),
    measure   = "GEN",
    only_data = TRUE
  )
  glmm_args <- list(
    ai        = c(1, 2, 3),
    bi        = c(9, 8, 7),
    ci        = c(2, 1, 2),
    di        = c(8, 9, 8),
    measure   = "OR",
    only_data = TRUE
  )

  constructors <- list(
    brma       = list(fun = brma,       args = norm_args),
    RoBMA      = list(fun = RoBMA,      args = norm_args),
    BMA        = list(fun = BMA,        args = norm_args),
    bselmodel  = list(fun = bselmodel,  args = norm_args),
    bPET       = list(fun = bPET,       args = norm_args),
    bPEESE     = list(fun = bPEESE,     args = norm_args),
    brma.glmm  = list(fun = brma.glmm,  args = glmm_args),
    BMA.glmm   = list(fun = BMA.glmm,   args = glmm_args)
  )

  constructor_silent <- function(constructor) {

    object <- do.call(constructor[["fun"]], constructor[["args"]])

    return(object[["fit_control"]][["silent"]])
  }

  RoBMA.options(silent = TRUE)
  for (constructor_name in names(constructors)) {
    expect_true(
      constructor_silent(constructors[[constructor_name]]),
      info = constructor_name
    )
  }

  RoBMA.options(silent = FALSE)
  for (constructor_name in names(constructors)) {
    expect_false(
      constructor_silent(constructors[[constructor_name]]),
      info = constructor_name
    )
  }

  explicit_args             <- norm_args
  explicit_args[["silent"]] <- TRUE
  expect_true(constructor_silent(list(fun = BMA, args = explicit_args)))

  RoBMA.options(silent = TRUE)
  explicit_args[["silent"]] <- FALSE
  expect_false(constructor_silent(list(fun = RoBMA, args = explicit_args)))
})


test_that("convergence checks expose BayesTools routing controls", {

  checks <- set_convergence_checks()
  autofit <- set_autofit_control()

  expect_true(all(c(
    "max_Rhat", "min_ESS", "max_error", "max_SD_error",
    "check_indicators", "monitor", "allow_not_assessable"
  ) %in% names(checks)))
  expect_true(checks[["check_indicators"]])
  expect_true(autofit[["check_indicators"]])
  expect_null(checks[["monitor"]])
  expect_false(checks[["allow_not_assessable"]])

  routed <- set_convergence_checks(
    check_indicators     = FALSE,
    monitor              = c("mu", "tau"),
    allow_not_assessable = TRUE
  )
  expect_false(routed[["check_indicators"]])
  expect_false(set_autofit_control(
    check_indicators = FALSE
  )[["check_indicators"]])
  expect_identical(routed[["monitor"]], c("mu", "tau"))
  expect_true(routed[["allow_not_assessable"]])
  expect_false("remove_failed" %in% names(checks))
  expect_false("balance_probability" %in% names(checks))
  expect_error(set_convergence_checks(remove_failed = TRUE), "unused argument")
  expect_error(set_convergence_checks(balance_probability = FALSE), "unused argument")
})

test_that("control updates distinguish missing and explicit indicator settings", {

  old_autofit <- set_autofit_control()
  old_autofit[["check_indicators"]] <- NULL
  expect_true(RoBMA:::.update_autofit_control(
    old_autofit,
    list()
  )[["check_indicators"]])
  expect_false(RoBMA:::.update_autofit_control(
    set_autofit_control(check_indicators = FALSE),
    list()
  )[["check_indicators"]])
  expect_false(RoBMA:::.update_autofit_control(
    set_autofit_control(),
    list(check_indicators = FALSE)
  )[["check_indicators"]])

  old_checks <- set_convergence_checks()
  old_checks[["check_indicators"]] <- NULL
  expect_true(RoBMA:::.update_convergence_checks(
    old_checks,
    list()
  )[["check_indicators"]])
  expect_false(RoBMA:::.update_convergence_checks(
    set_convergence_checks(check_indicators = FALSE),
    list()
  )[["check_indicators"]])
  expect_false(RoBMA:::.update_convergence_checks(
    set_convergence_checks(),
    list(check_indicators = FALSE)
  )[["check_indicators"]])
})


test_that("brma_samples metadata must match sample rows", {

  expect_error(
    RoBMA:::.new_brma_samples(
      samples  = matrix(seq_len(3), ncol = 1),
      n_chains = 2,
      n_iter   = 2,
      title    = "Invalid"
    ),
    "Invalid brma_samples chain metadata"
  )

  samples <- RoBMA:::.new_brma_samples(
    samples  = matrix(seq_len(3), ncol = 1),
    n_chains = 1,
    n_iter   = 3,
    title    = "Valid"
  )
  attr(samples, "nchains") <- 2
  attr(samples, "niter")   <- 2

  expect_error(
    RoBMA:::.brma_samples_to_mcmc.list(samples),
    "Invalid brma_samples chain metadata"
  )
})


test_that("brma_samples chain info preserves only balanced chains", {

  balanced_fit <- list(
    mcmc = list(
      matrix(seq_len(4), nrow = 2),
      matrix(seq_len(4), nrow = 2)
    )
  )
  unequal_fit <- list(
    mcmc = list(
      matrix(seq_len(2), nrow = 1),
      matrix(seq_len(4), nrow = 2)
    )
  )

  expect_equal(
    RoBMA:::.brma_samples_chain_info(balanced_fit, n_samples = 4),
    list(n_chains = 2L, n_iter = 2L)
  )
  expect_equal(
    RoBMA:::.brma_samples_chain_info(unequal_fit, n_samples = 3),
    list(n_chains = 1L, n_iter = 3L)
  )
})


test_that("posterior model indicators must be exact integers", {

  samples <- matrix(
    c(1, 2),
    ncol     = 1L,
    dimnames = list(NULL, "mu_indicator")
  )
  expect_identical(
    .extract_posterior_indicator(
      posterior_samples = samples,
      parameter         = "mu",
      prior             = list(NULL, NULL)
    ),
    c(1L, 2L)
  )

  samples[1L, 1L] <- 1 + .Machine$double.eps
  expect_error(
    .extract_posterior_indicator(
      posterior_samples = samples,
      parameter         = "mu",
      prior             = list(NULL, NULL)
    ),
    "integer-valued"
  )
})
