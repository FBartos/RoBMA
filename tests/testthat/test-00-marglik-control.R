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


test_that("add_marglik forwards inherited and explicit parallel controls", {

  old_max_cores <- RoBMA.get_option("max_cores")
  on.exit(RoBMA.options(max_cores = old_max_cores), add = TRUE)
  RoBMA.options(max_cores = 4L)

  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .check_marglik_available = function(object, caller) NULL,
    .marglik = function(object, cores){
      seen$cores <- cores
      list(logml = 0)
    },
    .brma_mv_attach_marglik_target_metadata = function(marglik, object) marglik
  )
  object <- structure(
    list(fit_control = list(parallel = TRUE, cores = 3L)),
    class = c("brma.norm", "brma")
  )

  inherited <- add_marglik(object)
  expect_identical(seen$cores, 3L)
  expect_identical(inherited[["marglik"]][["logml"]], 0)

  explicit <- add_marglik(object, parallel = TRUE, cores = 8L)
  expect_identical(seen$cores, 4L)
  expect_identical(explicit[["marglik"]][["logml"]], 0)
})
