context("Fitting errors")


test_that(".stop_fit_errors preserves stored fitting diagnostics", {

  expect_true(.stop_fit_errors(list(
    errors        = NULL,
    has_posterior = TRUE
  )))
  expect_error(
    .stop_fit_errors(list(
      errors        = "Original JAGS fitting diagnostic.",
      has_posterior = FALSE
    )),
    "^Original JAGS fitting diagnostic[.]$"
  )
  expect_error(
    .stop_fit_errors(list(
      errors        = NULL,
      has_posterior = FALSE
    )),
    "failed before posterior samples"
  )
})


test_that("constructors stop on fitting errors before creating summaries", {

  summary_called <- FALSE
  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      return(list(
        errors        = "Original constructor fitting diagnostic.",
        has_posterior = FALSE
      ))
    },
    .object_summary = function(object) {

      summary_called <<- TRUE
      stop("Summary should not be created.", call. = FALSE)
    },
    .package = "RoBMA"
  )

  expect_error(
    brma.norm(
      yi                        = c(0.10, 0.20),
      sei                       = c(0.05, 0.06),
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      silent                    = TRUE
    ),
    "^Original constructor fitting diagnostic[.]$"
  )
  expect_false(summary_called)

  expect_error(
    brma.mv(
      yi                        = c(0.10, 0.20),
      V                         = diag(c(0.05^2, 0.06^2)),
      random                    = NULL,
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      silent                    = TRUE
    ),
    "^Original constructor fitting diagnostic[.]$"
  )
  expect_false(summary_called)
})


test_that("all public fitting constructors preserve original fitting errors", {

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      return(list(
        errors        = "Original package-wide fitting diagnostic.",
        has_posterior = FALSE
      ))
    },
    .object_summary = function(object) {

      stop("Summary should not be created.", call. = FALSE)
    },
    .package = "RoBMA"
  )

  normal_args <- list(
    yi                        = c(0.10, 0.20),
    sei                       = c(0.05, 0.06),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    silent                    = TRUE
  )
  glmm_args <- list(
    ai      = c(1, 2),
    bi      = c(9, 8),
    ci      = c(2, 1),
    di      = c(8, 9),
    measure = "OR",
    silent  = TRUE
  )
  mv_args <- list(
    yi                        = normal_args[["yi"]],
    V                         = diag(normal_args[["sei"]]^2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    silent                    = TRUE
  )
  constructors <- list(
    BMA.norm   = function() do.call(BMA.norm, normal_args),
    bPET       = function() do.call(bPET, normal_args),
    bPET.mv    = function() do.call(bPET.mv, mv_args),
    bPEESE     = function() do.call(bPEESE, normal_args),
    bPEESE.mv  = function() do.call(bPEESE.mv, mv_args),
    bselmodel  = function() do.call(bselmodel, normal_args),
    RoBMA      = function() do.call(RoBMA, normal_args),
    brma.glmm  = function() do.call(brma.glmm, glmm_args),
    BMA.glmm   = function() do.call(BMA.glmm, glmm_args)
  )

  for (constructor in constructors) {
    expect_error(
      constructor(),
      "^Original package-wide fitting diagnostic[.]$"
    )
  }
})
