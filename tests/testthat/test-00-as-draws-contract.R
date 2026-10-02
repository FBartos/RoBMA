context("BayesTools draw materialization contract")


test_that("brma draw conversion delegates structural geometry to BayesTools", {

  values <- matrix(
    c(.25, 0, 11, .25, 0, 12, .25, 0, 13),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("mu", "tau", "sampling_z[1]"))
  )
  materialized <- coda::mcmc.list(coda::mcmc(
    values,
    start = 7,
    thin  = 3
  ))
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  attr(fit, "prior_list") <- list()
  object <- structure(list(fit = fit), class = "brma")
  checked <- NULL
  include_internal <- NULL
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(fit, requires) {
      checked <<- requires
      invisible(TRUE)
    },
    JAGS_ordered_parameter_spec = function(...) list(),
    JAGS_materialize_draws = function(fit, parameters = NULL,
                                      include_internal = FALSE) {
      include_internal <<- include_internal
      materialized
    },
    .package = "BayesTools"
  )

  public <- .brma_to_mcmc.list(object)
  raw    <- .brma_to_mcmc.list(object, include_auxiliary = TRUE)

  expect_setequal(checked, c("parameter_map", "draw_geometry"))
  expect_false(include_internal)
  expect_identical(
    colnames(as.matrix(public[[1L]])),
    c("mu", "tau")
  )
  expect_identical(as.numeric(public[[1L]][, "mu"]), rep(.25, 3L))
  expect_identical(as.numeric(public[[1L]][, "tau"]), rep(0, 3L))
  expect_identical(coda::mcpar(public[[1L]]), c(7, 13, 3))
  expect_identical(raw, materialized)
})


test_that("zero-coordinate fits retain draws and chain timing", {

  materialized <- coda::mcmc.list(coda::mcmc(
    matrix(numeric(), nrow = 4L, ncol = 0L),
    start = 5,
    thin  = 2
  ))
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  object <- structure(list(fit = fit), class = "brma")
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(...) invisible(TRUE),
    JAGS_ordered_parameter_spec = function(...) list(),
    JAGS_materialize_draws = function(...) materialized,
    .package = "BayesTools"
  )

  observed <- .brma_to_mcmc.list(object)
  expect_identical(dim(observed[[1L]]), c(4L, 0L))
  expect_identical(coda::mcpar(observed[[1L]]), c(5, 11, 2))

  draws <- as_draws_matrix(object)
  expect_identical(posterior::ndraws(draws), 4L)
  expect_identical(posterior::nchains(draws), 1L)
  expect_length(posterior::variables(draws), 0L)
})


test_that("ordered extraction has one restricted source read and preserves chain geometry", {

  values <- matrix(c(1, 99, 96, 2, 98, 95, 3, 97, 94), 3, byrow = TRUE,
    dimnames = list(NULL, c("total", "increment", "increment2")))
  main <- coda::mcmc.list(coda::mcmc(values, start = 7, thin = 3))
  source <- coda::mcmc.list(coda::mcmc(cbind(total = 1:3, gamma = c(.2, .3, .4),
    total2 = 4:6, gamma2 = c(.5, .6, .7)), start = 7, thin = 3))
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  attr(fit, "prior_list") <- list()
  object <- structure(list(fit = fit), class = "brma")
  calls <- list()
  conflict <- FALSE
  source_error <- NULL
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(...) invisible(TRUE),
    parameter_coordinates = function(...) data.frame(
      coordinate_name = c("total", "gamma", "total2", "gamma2", "optional"),
      monitor_status = c("sampled", "sampled", "sampled", "sampled", "unmonitored")),
    JAGS_materialize_draws = function(fit, parameters = NULL, include_internal = FALSE) {
      calls[[length(calls) + 1L]] <<- list(parameters = parameters, internal = include_internal)
      if (include_internal) source else main
    },
    JAGS_ordered_parameter_spec = function(fit, draws = NULL) {
      spec <- list(coefficient_names = "increment",
        source_coordinates = c("total", "gamma", "optional"),
        total_auxiliary_coordinates = character())
      spec2 <- list(coefficient_names = "increment2",
        source_coordinates = c("total2", "gamma", "gamma2", "optional"),
        total_auxiliary_coordinates = character())
      if (!is.null(draws)) {
        if (!is.null(source_error)) {
          stop(source_error)
        }
        spec[["sampled_parameters"]] <- list(
          values = cbind(total = as.numeric(draws[, "total"]) + as.integer(conflict),
            share = as.numeric(draws[, "gamma"])), allocation_names = "share")
        spec2[["sampled_parameters"]] <- list(
          values = cbind(total2 = as.numeric(draws[, "total2"]),
            shared_share = as.numeric(draws[, "gamma"]),
            share2 = as.numeric(draws[, "gamma2"])),
          allocation_names = c("shared_share", "share2"))
      }
      list(term = spec, term2 = spec2)
    }, .package = "BayesTools"
  )
  extracted <- .brma_to_mcmc.list(object)
  expect_identical(calls, list(list(parameters = NULL, internal = FALSE),
    list(parameters = c("total", "gamma", "total2", "gamma2"), internal = TRUE)))
  expect_identical(colnames(extracted[[1L]]), c("total", "share", "total2", "shared_share", "share2"))
  expect_identical(as.numeric(extracted[[1L]][, "share"]), c(.2, .3, .4))
  expect_identical(as.numeric(extracted[[1L]][, "shared_share"]), c(.2, .3, .4))
  expect_identical(as.numeric(extracted[[1L]][, "share2"]), c(.5, .6, .7))
  expect_identical(as.numeric(extracted[[1L]][, "total2"]), c(4, 5, 6))
  expect_identical(coda::mcpar(extracted[[1L]]), c(7, 13, 3))
  conflict <- TRUE
  error <- tryCatch(.brma_to_mcmc.list(object), error = identity)
  expect_s3_class(error, "RoBMA_ordered_draws_conflict")
  expect_s3_class(error, "BayesTools_refit_required")
  expect_identical(conditionMessage(error), paste0(
    "Ordered sampled parameter 'total' conflicts with existing fitted draws. ",
    "Refit the model with the current RoBMA/BayesTools build."))
  conflict <- FALSE
  source_error <- simpleError("Unrelated accessor failure.")
  error <- tryCatch(.brma_to_mcmc.list(object), error = identity)
  expect_identical(error, source_error)
  source_error <- NULL
  source <- coda::mcmc.list(coda::mcmc(as.matrix(source[[1L]]), start = 8, thin = 3))
  expect_error(.brma_to_mcmc.list(object),
    "Ordered source draws disagree with the fitted draw geometry.", fixed = TRUE)
})
