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
      silent          = silent
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
    .create_fit_data = function(...) list(),
    .create_jags_formula_args = function(...) {
      list(
        formula_list                          = NULL,
        formula_data_list                     = NULL,
        formula_prior_list                    = NULL,
        formula_scale_list                    = NULL,
        formula_random_prior_list             = NULL,
        formula_random_effects_compile_list   = NULL,
        add_parameters                        = NULL
      )
    },
    .create_model_syntax = function(...) "model {}",
    .convergence_structural_parameters = function(...) character(),
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    JAGS_extend = mock_extend,
    JAGS_check_convergence = function(...) TRUE,
    .package = "BayesTools"
  )

  expect_silent(result <- .fit(object, extend = TRUE))

  expect_identical(
    names(forwarded),
    c("fit", "autofit_control", "parallel", "cores", "silent")
  )
  expect_identical(forwarded[["fit"]], stored_fit)
  expect_identical(forwarded[["autofit_control"]], object[["autofit_control"]])
  expect_false("seed" %in% names(forwarded))
  expect_true(result[["has_posterior"]])
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
    "cannot reseed an existing JAGS fit"
  )
})
