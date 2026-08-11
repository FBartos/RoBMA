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

test_that("targeted convergence controls retain product-space indicator checks", {

  set.seed(45)
  chain_1 <- cbind(
    mu = stats::rnorm(200),
    mu_indicator = 0,
    mu_inclusion = 0
  )
  chain_2 <- cbind(
    mu = stats::rnorm(200),
    mu_indicator = 1,
    mu_inclusion = 1
  )
  fit <- list(
    mcmc = coda::mcmc.list(
      coda::mcmc(chain_1),
      coda::mcmc(chain_2)
    ),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"
  attr(fit, "prior_list") <- list(mu = BayesTools::prior_spike_and_slab(
    BayesTools::prior("normal", list(0, 1)),
    prior_inclusion = BayesTools::prior("beta", list(1, 1))
  ))

  checked <- RoBMA:::.recheck_brma_fit(list(
    fit    = fit,
    priors = list(),
    convergence_checks = set_convergence_checks(
      max_Rhat         = 1.05,
      min_ESS          = NULL,
      max_error        = NULL,
      max_SD_error     = NULL,
      check_indicators = TRUE,
      monitor          = "mu"
    )
  ))

  expect_false(checked[["converged"]])
  expect_match(checked[["warnings"]], "indicator|not assessable|R-hat")
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


test_that("control updates distinguish absent and explicitly null settings", {

  old_autofit <- set_autofit_control(
    max_Rhat     = 1.01,
    min_ESS      = 1000,
    max_error    = 0.01,
    max_SD_error = 0.02,
    max_time     = list(time = 2, unit = "hours"),
    restarts     = 4,
    max_extend   = 5,
    monitor      = c("mu", "tau")
  )
  expect_equal(
    RoBMA:::.update_autofit_control(old_autofit, NULL),
    old_autofit
  )
  expect_equal(
    RoBMA:::.update_autofit_control(old_autofit, list()),
    old_autofit
  )

  cleared_autofit <- RoBMA:::.update_autofit_control(
    old_autofit,
    list(
      max_Rhat     = NULL,
      min_ESS      = NULL,
      max_error    = NULL,
      max_SD_error = NULL,
      max_time     = NULL,
      restarts     = NULL,
      max_extend   = NULL,
      monitor      = NULL
    )
  )
  for (setting in c(
    "max_Rhat", "min_ESS", "max_error", "max_SD_error", "max_time",
    "restarts", "max_extend", "monitor"
  )) {
    expect_null(cleared_autofit[[setting]], info = setting)
  }
  expect_identical(
    cleared_autofit[["sample_extend"]],
    old_autofit[["sample_extend"]]
  )

  old_checks <- set_convergence_checks(
    max_Rhat     = 1.01,
    min_ESS      = 1000,
    max_error    = 0.01,
    max_SD_error = 0.02,
    monitor      = c("mu", "tau")
  )
  expect_equal(
    RoBMA:::.update_convergence_checks(old_checks, NULL),
    old_checks
  )
  expect_equal(
    RoBMA:::.update_convergence_checks(old_checks, list()),
    old_checks
  )

  cleared_checks <- RoBMA:::.update_convergence_checks(
    old_checks,
    list(
      max_Rhat     = NULL,
      min_ESS      = NULL,
      max_error    = NULL,
      max_SD_error = NULL,
      monitor      = NULL
    )
  )
  for (setting in c(
    "max_Rhat", "min_ESS", "max_error", "max_SD_error", "monitor"
  )) {
    expect_null(cleared_checks[[setting]], info = setting)
  }
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


test_that("nested brma_samples results preserve their summary-table structure", {

  study <- .new_brma_samples(
    samples  = matrix(
      1:8,
      nrow     = 4L,
      dimnames = list(NULL, c("study 1", "study 2"))
    ),
    n_chains = 1L,
    n_iter   = 4L,
    title    = "Study effects"
  )
  estimate <- .new_brma_samples(
    samples  = matrix(
      9:16,
      nrow     = 4L,
      dimnames = list(NULL, c("estimate 1", "estimate 2"))
    ),
    n_chains = 1L,
    n_iter   = 4L,
    title    = "Estimate effects"
  )
  result <- .new_brma_samples_list(list(
    location = list(study = study, estimate = estimate),
    scale    = list(study = study)
  ))

  long_table <- as.data.frame(result)
  tables     <- as.data.frame(result, format = "list")

  expect_s3_class(result, "brma_samples_list")
  expect_s3_class(long_table, "data.frame")
  expect_setequal(
    long_table[["component"]],
    c("location/study", "location/estimate", "scale/study")
  )
  expect_true(all(long_table[["parameter"]] %in% c(
    "study 1", "study 2", "estimate 1", "estimate 2"
  )))
  expect_identical(
    names(long_table),
    c("component", "parameter", "Mean", "Median", "CI_0.025", "CI_0.975")
  )
  expect_identical(data.frame(result), data.frame(long_table))
  expect_named(tables, c("location", "scale"))
  expect_named(tables[["location"]], c("study", "estimate"))
  expect_identical(
    tables[["location"]][["study"]],
    as.data.frame(study)
  )
  expect_identical(
    tables[["location"]][["estimate"]],
    as.data.frame(estimate)
  )
  expect_identical(
    tables[["scale"]][["study"]],
    as.data.frame(study)
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
