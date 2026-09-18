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


test_that("the cached bridge sample layout reproduces the general assembly", {

  # A bridge fixes the one-row sample matrix's columns once and fills the
  # values by position. Every state must give the matrix the general path
  # assembles from the parameter list, the bridge row metadata and the nodes.
  parameters <- list(mu = 0.5, beta = c(1, 2, 3), empty = numeric(),
                     tau = 0.25)
  attr(parameters, "posterior_samples") <- matrix(
    c(0.5, 7, 9), nrow = 1L,
    dimnames = list(NULL, c("mu", "extra_a", "extra_b"))
  )
  context <- structure(
    list(nodes = c(tau = 0.25, extra_a = 7, node_c = 4)),
    class = c("BayesTools_bridge_context", "list")
  )

  reference <- .marglik_bridge_posterior_samples(parameters, context)
  cache     <- new.env(parent = emptyenv())
  expect_identical(
    .marglik_bridge_posterior_samples(parameters, context, cache),
    reference
  )

  # a second state of the same shape replays the layout
  other <- parameters
  other[["beta"]] <- c(-1, -2, -3)
  attr(other, "posterior_samples") <- matrix(
    c(0.5, 70, 90), nrow = 1L,
    dimnames = list(NULL, c("mu", "extra_a", "extra_b"))
  )
  other_context <- context
  other_context[["nodes"]] <- c(tau = 0.25, extra_a = 70, node_c = 40)
  expect_identical(
    .marglik_bridge_posterior_samples(other, other_context, cache),
    .marglik_bridge_posterior_samples(other, other_context)
  )

  # another shape rebuilds the layout instead of replaying a stale one
  changed <- parameters
  changed[["beta"]] <- c(1, 2)
  expect_identical(
    .marglik_bridge_posterior_samples(changed, context, cache),
    .marglik_bridge_posterior_samples(changed, context)
  )

  # a context without nodes and a call without one keep their own layouts
  expect_identical(
    .marglik_bridge_posterior_samples(parameters, NULL, cache),
    .marglik_bridge_posterior_samples(parameters, NULL)
  )
  expect_identical(
    .marglik_bridge_posterior_samples(parameters, context, cache),
    reference
  )
})
