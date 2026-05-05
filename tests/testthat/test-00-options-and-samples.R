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


test_that("convergence checks expose only active thresholds", {

  checks <- set_convergence_checks()

  expect_true(all(c("max_Rhat", "min_ESS", "max_error", "max_SD_error") %in% names(checks)))
  expect_false("remove_failed" %in% names(checks))
  expect_false("balance_probability" %in% names(checks))
  expect_error(set_convergence_checks(remove_failed = TRUE), "unused argument")
  expect_error(set_convergence_checks(balance_probability = FALSE), "unused argument")
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
