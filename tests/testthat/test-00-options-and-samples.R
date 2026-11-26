context("Options and brma_samples metadata")

test_that("RoBMA options expose only public options", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)

  expect_true("silent" %in% names(old_options))
  expect_false(old_options[["silent"]])
  expect_equal(old_options[["default_lograte.sd"]], 1)
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
  expect_error(RoBMA.options(default_lograte.sd = 0), "must be > 0")
  expect_error(RoBMA.options(default_bias_PET.scale = expression(stop("boom"))), "finite number")

  updated <- RoBMA.options(max_cores = 1L)
  expect_equal(updated[["max_cores"]], 1L)
})

test_that("worker output is opt-in and reaches fitting and extension", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  log_path <- tempfile("RoBMA worker ", fileext = ".log")
  for (value in list(NULL, TRUE, NA_character_, c("a", "b"))) {
    expect_error(RoBMA.options(jags.worker_output = value),
      "Option 'jags.worker_output' must be a file path or an empty string.", fixed = TRUE)
  }
  RoBMA.options(jags.worker_output = "")
  object <- brma(yi = c(.1, .2, .3), sei = c(.2, .2, .2),
    measure = "SMD", only_priors = TRUE, silent = TRUE)
  received <- list()
  backend <- function(..., worker_output) {
    received[length(received) + 1L] <<- list(worker_output)
    simpleError("Worker output forwarding probe")
  }
  testthat::local_mocked_bindings(JAGS_fit = backend, JAGS_extend = backend,
    .package = "BayesTools")
  .fit(object)
  expect_null(received[[1L]])
  RoBMA.options(jags.worker_output = log_path)
  expect_identical(RoBMA.get_option("jags.worker_output"), log_path)
  .fit(object)
  object$fit <- list("extension probe")
  .fit(object, extend = TRUE)
  expect_identical(received[2:3], list(log_path, log_path))
  RoBMA.options(jags.worker_output = "")
  .fit(object, extend = TRUE)
  expect_null(received[[4L]])
  expect_false(file.exists(log_path))
})


test_that("brma-family constructors inherit the silent option and model averaging is quiet", {

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
  mv_args <- list(
    yi        = norm_args[["yi"]],
    V         = diag(norm_args[["sei"]]^2),
    measure   = "GEN",
    only_data = TRUE
  )

  constructors <- list(
    brma       = list(fun = brma,       args = norm_args),
    RoBMA      = list(fun = RoBMA,      args = norm_args),
    BMA        = list(fun = BMA,        args = norm_args),
    bselmodel  = list(fun = bselmodel,  args = norm_args),
    bPET       = list(fun = bPET,       args = norm_args),
    bPEESE     = list(fun = bPEESE,     args = norm_args),
    brma.mv    = list(fun = brma.mv,    args = mv_args),
    bselmodel.mv = list(fun = bselmodel.mv, args = mv_args),
    BMA.mv     = list(fun = BMA.mv,     args = mv_args),
    RoBMA.mv   = list(fun = RoBMA.mv,   args = mv_args),
    bPET.mv    = list(fun = bPET.mv,    args = mv_args),
    bPEESE.mv  = list(fun = bPEESE.mv,  args = mv_args),
    brma.glmm  = list(fun = brma.glmm,  args = glmm_args),
    BMA.glmm   = list(fun = BMA.glmm,   args = glmm_args)
  )

  constructor_silent <- function(constructor) {

    object <- do.call(constructor[["fun"]], constructor[["args"]])

    return(object[["fit_control"]][["silent"]])
  }

  # as in 4.0.0, the model-averaging constructors default to silent = TRUE
  # whatever the option, the brma family follows the option
  model_averaging <- c("RoBMA", "BMA", "BMA.glmm", "BMA.mv", "RoBMA.mv")
  expect_false(RoBMA.get_option("silent"))
  for (value in c(FALSE, TRUE)) {
    RoBMA.options(silent = value)
    for (constructor_name in names(constructors)) {
      expected <- if (constructor_name %in% model_averaging) TRUE else value
      expect_identical(
        constructor_silent(constructors[[constructor_name]]),
        expected,
        info = paste(constructor_name, "option", value)
      )
    }
  }

  RoBMA.options(silent = FALSE)
  for (constructor_name in c("BMA.norm", "BMA", "RoBMA", "BMA.glmm", "BMA.mv", "RoBMA.mv")) {
    expect_true(
      eval(formals(get(constructor_name))[["silent"]]),
      info = constructor_name
    )
  }

  explicit_args             <- norm_args
  explicit_args[["silent"]] <- TRUE
  expect_true(constructor_silent(list(fun = BMA, args = explicit_args)))
  expect_true(constructor_silent(list(fun = brma, args = explicit_args)))

  RoBMA.options(silent = TRUE)
  explicit_args[["silent"]] <- FALSE
  expect_false(constructor_silent(list(fun = RoBMA, args = explicit_args)))
  expect_false(constructor_silent(list(fun = brma, args = explicit_args)))
})


test_that("automatic extensions default to sample_extend = 'auto' and keep a fixed number", {

  expect_identical(formals(set_autofit_control)[["sample_extend"]], "auto")
  expect_identical(set_autofit_control()[["sample_extend"]], "auto")
  expect_identical(set_autofit_control(sample_extend = 2000)[["sample_extend"]], 2000)
  # a named "auto" is stored as the plain string, which BayesTools runs as "auto"
  expect_identical(set_autofit_control(sample_extend = c(step = "auto"))[["sample_extend"]], "auto")
  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "SMD",
                 only_priors = TRUE, silent = TRUE)
  expect_identical(object[["autofit_control"]][["sample_extend"]], "auto")

  # validated by BayesTools: a positive whole number or "auto"
  for (value in list("Auto", c("auto", "auto"), NA_character_, 0, 1.5, TRUE, NULL)) {
    expect_error(set_autofit_control(sample_extend = value), "sample_extend", info = deparse(value))
  }
  # a stored fixed step stays when other settings change
  expect_identical(
    .update_autofit_control(set_autofit_control(sample_extend = 7), list(max_extend = 3))[["sample_extend"]],
    7
  )
})


test_that("the readiness settings of burn-in are validated and survive control updates", {

  burnin_names <- c("burnin_max_Rhat", "burnin_extend", "burnin_max_extend")
  defaults <- set_autofit_control()
  expect_identical(utils::tail(names(formals(set_autofit_control)), 3L), burnin_names)
  expect_identical(defaults[burnin_names],
    list(burnin_max_Rhat = "auto", burnin_extend = 500, burnin_max_extend = 10))
  # an explicit NULL turns the check off and is kept as a named element
  off <- set_autofit_control(burnin_max_Rhat = NULL)
  expect_true("burnin_max_Rhat" %in% names(off))
  expect_null(off[["burnin_max_Rhat"]])
  # a number fixes the threshold
  expect_identical(set_autofit_control(burnin_max_Rhat = 1.8)[["burnin_max_Rhat"]], 1.8)

  # validated by BayesTools
  expect_error(set_autofit_control(burnin_extend = 50),
    "Checking 'autofit_control':\n\tThe 'burnin_extend' must be equal or higher than 100.", fixed = TRUE)
  expect_error(set_autofit_control(burnin_max_Rhat = 0.5), "'burnin_max_Rhat' must be equal or higher than 1.")
  expect_error(set_autofit_control(burnin_max_Rhat = Inf), "'burnin_max_Rhat' argument must contain only finite values.")
  expect_error(set_autofit_control(burnin_max_Rhat = "Auto"),
    "Checking 'autofit_control':\n\tThe 'Auto' values are not recognized by the 'burnin_max_Rhat' argument.",
    fixed = TRUE)
  expect_error(set_autofit_control(burnin_max_Rhat = c("auto", "auto")),
    "The 'burnin_max_Rhat' argument must have length '1'.", fixed = TRUE)
  expect_error(set_autofit_control(burnin_max_extend = -1), "'burnin_max_extend' must be equal or higher than 0.")
  expect_error(set_autofit_control(burnin_extend = NULL), "The 'burnin_extend' argument cannot be NULL.")

  # the convergence checks do not get them
  expect_false(any(grepl("^burnin_", names(formals(set_convergence_checks)))))
  expect_false(any(grepl("^burnin_", names(set_convergence_checks()))))
  expect_false(any(grepl("^burnin_", names(.update_convergence_checks(set_convergence_checks(), list())))))

  # update(): a supplied value or NULL replaces the stored one, an omitted one is kept
  stored <- set_autofit_control(burnin_max_Rhat = 1.3, burnin_extend = 300, burnin_max_extend = 4)
  expect_identical(.update_autofit_control(stored, NULL)[burnin_names], stored[burnin_names])
  expect_identical(.update_autofit_control(stored, list(burnin_extend = 200))[burnin_names],
    list(burnin_max_Rhat = 1.3, burnin_extend = 200, burnin_max_extend = 4))
  cleared <- .update_autofit_control(stored, list(burnin_max_Rhat = NULL))
  expect_true("burnin_max_Rhat" %in% names(cleared))
  expect_null(cleared[["burnin_max_Rhat"]])
  later <- .update_autofit_control(cleared, list(max_extend = 2))
  expect_true("burnin_max_Rhat" %in% names(later))
  expect_null(later[["burnin_max_Rhat"]])
  # "auto" is kept when omitted and replaced by a number, and the reverse
  expect_identical(.update_autofit_control(defaults, list(max_extend = 2))[["burnin_max_Rhat"]], "auto")
  expect_identical(.update_autofit_control(defaults, list(burnin_max_Rhat = 1.3))[["burnin_max_Rhat"]], 1.3)
  expect_identical(.update_autofit_control(stored, list(burnin_max_Rhat = "auto"))[["burnin_max_Rhat"]], "auto")
  # a control stored before the settings existed takes their defaults
  legacy <- set_autofit_control(max_extend = 3)
  legacy[burnin_names] <- NULL
  expect_identical(.update_autofit_control(legacy, list())[burnin_names], defaults[burnin_names])

  # the constructors store the control as given
  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "SMD", only_priors = TRUE, silent = TRUE,
    autofit_control = set_autofit_control(burnin_max_Rhat = NULL, burnin_extend = 200))
  expect_true("burnin_max_Rhat" %in% names(object[["autofit_control"]]))
  expect_null(object[["autofit_control"]][["burnin_max_Rhat"]])
  expect_identical(object[["autofit_control"]][["burnin_extend"]], 200)
})


test_that("convergence checks expose BayesTools routing controls", {

  checks <- set_convergence_checks()
  autofit <- set_autofit_control()

  expect_true(all(c(
    "max_Rhat", "min_ESS", "max_error", "max_SD_error",
    "check_indicators", "monitor", "allow_not_assessable"
  ) %in% names(checks)))
  # the product-space indicator checks are opt-in, as in 4.0.0
  expect_false(checks[["check_indicators"]])
  expect_false(autofit[["check_indicators"]])
  expect_null(checks[["monitor"]])
  expect_false(checks[["allow_not_assessable"]])

  routed <- set_convergence_checks(
    check_indicators     = TRUE,
    monitor              = c("mu", "tau"),
    allow_not_assessable = TRUE
  )
  expect_true(routed[["check_indicators"]])
  expect_true(set_autofit_control(
    check_indicators = TRUE
  )[["check_indicators"]])
  expect_identical(routed[["monitor"]], c("mu", "tau"))
  expect_true(routed[["allow_not_assessable"]])
  expect_identical(names(checks), c(
    "max_Rhat", "min_ESS", "max_error", "max_SD_error",
    "check_indicators", "monitor", "allow_not_assessable"
  ))
  expect_error(
    set_convergence_checks(allow_not_assessable = NA),
    "Checking 'convergence_checks':\n\tThe 'allow_not_assessable' argument cannot contain NA/NaN values.",
    fixed = TRUE
  )
  # the setting alone is a valid convergence-check list
  only_allow <- RoBMA:::.check_and_list_convergence_checks(list(allow_not_assessable = TRUE))
  expect_true(only_allow[["allow_not_assessable"]])
  expect_false(only_allow[["check_indicators"]])
  expect_false(RoBMA:::.check_and_list_convergence_checks(list(max_Rhat = 1.05))[["allow_not_assessable"]])
  expect_false("remove_failed" %in% names(checks))
  expect_false("balance_probability" %in% names(checks))
  expect_error(set_convergence_checks(remove_failed = TRUE), "unused argument")
  expect_error(set_convergence_checks(balance_probability = FALSE), "unused argument")
})

test_that("the automatic-fitting control has no allow_not_assessable", {

  # removed after 4.1.1: automatic fitting never extends a model for
  # parameters that cannot be assessed; set_convergence_checks() keeps it
  expect_identical(names(set_autofit_control()), c(
    "max_Rhat", "min_ESS", "max_error", "max_SD_error", "max_time",
    "sample_extend", "restarts", "max_extend", "check_indicators", "monitor",
    "burnin_max_Rhat", "burnin_extend", "burnin_max_extend"
  ))
  expect_false("allow_not_assessable" %in% names(formals(set_autofit_control)))
  expect_error(set_autofit_control(allow_not_assessable = TRUE), "unused argument")
  expect_error(
    RoBMA(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "SMD", only_priors = TRUE, silent = TRUE,
          autofit_control = list(sample_extend = "auto", allow_not_assessable = TRUE)),
    "The 'allow_not_assessable' objects are not recognized by the 'autofit_control' argument.",
    fixed = TRUE
  )

  # a control stored by RoBMA 4.1.1 is rebuilt without it, whatever its value
  for (allow in c(FALSE, TRUE)) {
    stored <- c(set_autofit_control(sample_extend = 1000), list(allow_not_assessable = allow))
    expect_identical(.update_autofit_control(stored, NULL), set_autofit_control(sample_extend = 1000))
    expect_identical(
      .update_autofit_control(stored, list(max_extend = 3)),
      set_autofit_control(sample_extend = 1000, max_extend = 3)
    )
  }
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
    )
  )
  class(fit) <- "BayesTools_fit"
  attr(fit, "prior_list") <- list(mu = BayesTools::prior_spike_and_slab(
    BayesTools::prior("normal", list(0, 1)),
    prior_inclusion = BayesTools::prior("beta", list(1, 1))
  ))
  fit <- as_bayestools_fit(fit)

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

test_that("convergence rechecks leave structural weights to BayesTools", {

  # A one-sided weight function fixes the most significant bin (omega[1]) at
  # one; BayesTools reports it as a structural constant from the prior list.
  set.seed(46)
  make_chain <- function() {
    cbind(
      mu         = stats::rnorm(200),
      "omega[1]" = 1,
      "omega[2]" = stats::runif(200, 0.5, 1),
      "omega[3]" = stats::runif(200, 0.2, 0.5)
    )
  }
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(make_chain()), coda::mcmc(make_chain()))
  )
  class(fit) <- "BayesTools_fit"
  weightfunction <- BayesTools::prior_weightfunction(
    "one-sided", c(.025, .05), BayesTools::wf_cumulative(c(1, 1, 1))
  )
  attr(fit, "prior_list") <- list(
    mu    = BayesTools::prior("normal", list(0, 1)),
    omega = weightfunction
  )
  fit <- as_bayestools_fit(fit)

  checked <- RoBMA:::.recheck_brma_fit(list(
    fit    = fit,
    priors = list(outcome = list(bias = weightfunction)),
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL
    )
  ))
  diagnostics <- attr(checked[["converged"]], "diagnostics", exact = TRUE)
  states <- stats::setNames(diagnostics[["state"]], diagnostics[["parameter"]])

  expect_identical(
    unname(states[c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,1]")]),
    c("structural_constant", "assessable", "assessable")
  )
  expect_true(isTRUE(as.logical(checked[["converged"]])))
  expect_length(checked[["warnings"]], 0L)
})

test_that("control updates distinguish missing and explicit indicator settings", {

  old_autofit <- set_autofit_control()
  old_autofit[["check_indicators"]] <- NULL
  expect_false(RoBMA:::.update_autofit_control(
    old_autofit,
    list()
  )[["check_indicators"]])
  expect_true(RoBMA:::.update_autofit_control(
    set_autofit_control(check_indicators = TRUE),
    list()
  )[["check_indicators"]])
  expect_true(RoBMA:::.update_autofit_control(
    set_autofit_control(),
    list(check_indicators = TRUE)
  )[["check_indicators"]])
  expect_false(RoBMA:::.update_autofit_control(
    set_autofit_control(check_indicators = TRUE),
    list(check_indicators = FALSE)
  )[["check_indicators"]])

  old_checks <- set_convergence_checks()
  old_checks[["check_indicators"]] <- NULL
  expect_false(RoBMA:::.update_convergence_checks(
    old_checks,
    list()
  )[["check_indicators"]])
  expect_true(RoBMA:::.update_convergence_checks(
    set_convergence_checks(check_indicators = TRUE),
    list()
  )[["check_indicators"]])
  expect_true(RoBMA:::.update_convergence_checks(
    set_convergence_checks(),
    list(check_indicators = TRUE)
  )[["check_indicators"]])
  expect_false(RoBMA:::.update_convergence_checks(
    set_convergence_checks(check_indicators = TRUE),
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
  expected_location_study <- as.data.frame(summary(study))
  expected_location_study[["component"]] <- "location/study"
  expected_location_estimate <- as.data.frame(summary(estimate))
  expected_location_estimate[["component"]] <- "location/estimate"
  expected_scale_study <- as.data.frame(summary(study))
  expected_scale_study[["component"]] <- "scale/study"

  expect_identical(
    tables[["location"]][["study"]],
    expected_location_study
  )
  expect_identical(
    tables[["location"]][["estimate"]],
    expected_location_estimate
  )
  expect_identical(
    tables[["scale"]][["study"]],
    expected_scale_study
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
