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
    JAGS_materialize_draws = function(fit, parameters = NULL,
                                      include_internal = FALSE) {
      include_internal <<- include_internal
      materialized
    },
    .package = "BayesTools"
  )

  public <- .brma_to_mcmc.list(object)
  raw    <- .brma_to_mcmc.list(object, include_auxiliary = TRUE)

  expect_setequal(checked, c("parameter_registry", "draw_geometry"))
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
