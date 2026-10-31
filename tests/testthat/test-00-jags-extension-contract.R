context("JAGS extension contract")


test_that(".fit forwards only supported BayesTools extension controls", {

  extend_formals <- formals(BayesTools::JAGS_extend)
  expect_false("seed" %in% names(extend_formals))

  forwarded <- NULL
  mock_extend <- function() NULL
  formals(mock_extend) <- extend_formals
  body(mock_extend) <- quote({
    forwarded <<- list(
      fit             = fit,
      autofit_control = autofit_control,
      parallel        = parallel,
      cores           = cores,
      silent          = silent,
      worker_output   = worker_output,
      runtime_cache   = runtime_cache
    )
    return(fit)
  })

  stored_fit <- structure(
    list(existing_chain = TRUE),
    class = c("BayesTools_fit", "list")
  )
  attr(stored_fit, "prior_list") <- list()
  object <- list(
    fit                = stored_fit,
    fit_control        = list(
      parallel = FALSE,
      cores    = NULL,
      silent   = TRUE,
      seed     = 719L
    ),
    autofit_control    = list(sample_extend = 20L),
    convergence_checks = list(
      max_Rhat             = 1.05,
      min_ESS              = 10,
      max_error            = 0.10,
      max_SD_error         = 0.10,
      check_indicators     = FALSE,
      monitor              = NULL,
      allow_not_assessable = FALSE
    ),
    data               = list(),
    priors             = list(outcome = list())
  )

  testthat::local_mocked_bindings(
    .create_fit_priors = function(...) list(),
    .create_fit_data = function(...) stop("fresh-fit data builder was called"),
    .create_jags_formula_args = function(...) stop("fresh-fit formula builder was called"),
    .create_model_syntax = function(...) stop("fresh-fit syntax builder was called"),
    .package = "RoBMA"
  )
  checked <- NULL
  testthat::local_mocked_bindings(
    JAGS_extend = mock_extend,
    JAGS_check_convergence = function(...) {
      checked <<- list(...)
      TRUE
    },
    .package = "BayesTools"
  )

  expect_silent(result <- .fit(object, extend = TRUE))

  expect_identical(
    names(forwarded),
    c("fit", "autofit_control", "parallel", "cores", "silent",
      "worker_output", "runtime_cache")
  )
  expect_identical(forwarded[["fit"]], stored_fit)
  expect_identical(forwarded[["autofit_control"]], object[["autofit_control"]])
  expect_false("seed" %in% names(forwarded))
  expect_null(forwarded[["runtime_cache"]])
  expected_output <- RoBMA.get_option("jags.worker_output")
  if (!nzchar(expected_output)) expected_output <- NULL
  expect_identical(forwarded[["worker_output"]], expected_output)
  expect_true(result[["has_posterior"]])
  # BayesTools classifies parameters by the fit's coordinate roles; RoBMA
  # passes no prior list and excludes no parameters by name.
  expect_false("prior_list" %in% names(checked))
  expect_false("add_parameters" %in% names(checked))
})


test_that("update.brma rejects reseeding a chain continuation", {

  object <- brma.norm(
    yi        = c(0.10, 0.20, 0.05),
    sei       = c(0.05, 0.06, 0.07),
    only_data = TRUE,
    silent    = TRUE
  )
  object[["fit"]] <- structure(
    list(has_posterior = TRUE),
    class = c("BayesTools_fit", "list")
  )

  expect_error(
    update(object, sample_extend = 10, seed = 719),
    "Unused argument in update.brma(): 'seed'", fixed = TRUE
  )
})


test_that("fresh and empty-fit attempts still build their inputs and retain error priors", {

  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "SMD",
    only_priors = TRUE)
  marker_priors <- list(mu = BayesTools::prior("point", list(0)))
  called <- character()
  testthat::local_mocked_bindings(
    .create_fit_priors = function(...) marker_priors,
    .create_fit_data = function(...) { called <<- c(called, "data"); list(y = 1) },
    .create_jags_formula_args = function(...) { called <<- c(called, "formula"); list() },
    .create_model_syntax = function(...) { called <<- c(called, "syntax"); "model {}" },
    .package = "RoBMA")
  forwarded <- NULL
  testthat::local_mocked_bindings(JAGS_fit = function(...) {
    forwarded <<- list(...)
    simpleError("Fresh fit diagnostic.")
  }, .package = "BayesTools")
  for (extend in c(FALSE, TRUE)) {
    called <- character()
    object[["fit"]] <- if (extend) list() else list(existing_chain = TRUE)
    result <- .fit(object, extend = extend)
    expect_identical(called, c("data", "formula", "syntax"))
    expect_identical(result[["errors"]], "Fresh fit diagnostic.")
    expect_false(result[["has_posterior"]])
    expect_identical(attr(result, "prior_list"), marker_priors)
    expect_false("add_parameters" %in% names(forwarded))
    expect_identical(forwarded[["seed"]], object[["fit_control"]][["seed"]])
  }
})


test_that("selection extension forwards its fitted runtime cache and preserves error priors", {

  object <- bselmodel(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "SMD",
    only_priors = TRUE)
  object[["fit"]] <- list(existing_chain = TRUE)
  stored_fit <- object[["fit"]]
  marker_cache <- list(cache_marker = TRUE)
  marker_priors <- list(mu = BayesTools::prior("point", list(0)))
  testthat::local_mocked_bindings(
    .create_fit_priors = function(...) marker_priors,
    .create_fit_data = function(...) stop("extension built data"),
    .create_jags_formula_args = function(...) stop("extension built formulas"),
    .create_model_syntax = function(...) stop("extension built syntax"),
    .selection_cache_runtime = function(fit) {
      expect_identical(fit, stored_fit)
      marker_cache
    }, .package = "RoBMA")
  received <- NULL
  testthat::local_mocked_bindings(JAGS_extend = function(...) {
    received <<- list(...)
    simpleError("Extension diagnostic.")
  }, .package = "BayesTools")
  result <- .fit(object, extend = TRUE)
  expect_identical(received[["runtime_cache"]], marker_cache)
  expect_false("seed" %in% names(received))
  expect_identical(result[["errors"]], "Extension diagnostic.")
  expect_identical(attr(result, "prior_list"), marker_priors)
})
