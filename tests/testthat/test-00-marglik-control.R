test_that("marginal-likelihood parallel settings inherit and respect max_cores", {

  old_max_cores <- RoBMA.get_option("max_cores")
  on.exit(RoBMA.options(max_cores = old_max_cores), add = TRUE)
  RoBMA.options(max_cores = 4L)

  object <- list(fit_control = list(parallel = TRUE, cores = 3L))
  expect_identical(
    .marglik_parallel_control(object),
    list(parallel = TRUE, cores = 3L)
  )
  expect_identical(
    .marglik_parallel_control(object, cores = 8L),
    list(parallel = TRUE, cores = 4L)
  )
  expect_identical(
    .marglik_parallel_control(object, parallel = FALSE, cores = 4L),
    list(parallel = FALSE, cores = 1L)
  )

  serial_object <- list(fit_control = list(parallel = FALSE, cores = 3L))
  expect_identical(
    .marglik_parallel_control(serial_object),
    list(parallel = FALSE, cores = 1L)
  )
  expect_identical(
    .marglik_parallel_control(serial_object, parallel = TRUE),
    list(parallel = TRUE, cores = 3L)
  )

  expect_error(
    .marglik_parallel_control(object, parallel = NA),
    "parallel"
  )
  expect_error(
    .marglik_parallel_control(object, cores = 0L),
    "cores"
  )
})


test_that("parallel marginal likelihoods load fitted bridge packages", {

  fit <- structure(list(), required_packages = c("RoBMA", "customPackage"))

  expect_null(.marglik_bridge_packages(fit, cores = 1L))
  expect_identical(
    .marglik_bridge_packages(fit, cores = 2L),
    c("BayesTools", "RoBMA", "customPackage")
  )
})


test_that("add_marglik forwards bridge and parallel controls", {

  old_max_cores <- RoBMA.get_option("max_cores")
  on.exit(RoBMA.options(max_cores = old_max_cores), add = TRUE)
  RoBMA.options(max_cores = 4L)

  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .check_marglik_available = function(object, caller) NULL,
    .marglik = function(object, cores, repetitions, method, maxiter, silent){
      seen$controls <- list(
        cores       = cores,
        repetitions = repetitions,
        method      = method,
        maxiter     = maxiter,
        silent      = silent
      )
      list(logml = 0)
    },
    .brma_mv_attach_marglik_target_metadata = function(marglik, object) marglik
  )
  object <- structure(
    list(fit_control = list(parallel = TRUE, cores = 3L)),
    class = c("brma.norm", "brma")
  )

  inherited <- add_marglik(object)
  expect_identical(
    seen$controls,
    list(
      cores       = 3L,
      repetitions = 1L,
      method      = "normal",
      maxiter     = 10000L,
      silent      = TRUE
    )
  )
  expect_identical(inherited[["marglik"]][["logml"]], 0)

  explicit <- add_marglik(
    object,
    parallel    = TRUE,
    cores       = 8L,
    repetitions = 4L,
    method      = "warp3",
    maxiter     = 2500L,
    silent      = FALSE
  )
  expect_identical(
    seen$controls,
    list(
      cores       = 4L,
      repetitions = 4L,
      method      = "warp3",
      maxiter     = 2500L,
      silent      = FALSE
    )
  )
  expect_identical(explicit[["marglik"]][["logml"]], 0)
})
