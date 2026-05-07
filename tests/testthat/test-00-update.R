context("update.brma")

mock_update_object <- function() {

  object <- brma.norm(
    yi        = c(0.10, 0.20, 0.05),
    sei       = c(0.05, 0.06, 0.07),
    only_data = TRUE,
    silent    = TRUE
  )

  object[["priors"]]       <- list(outcome = list())
  object[["fit"]]          <- structure(list(has_posterior = TRUE), class = "BayesTools_fit")
  object[["summary"]]      <- data.frame(Mean = 1, row.names = "mu")
  object[["coefficients"]] <- c(mu = 1)

  return(object)
}

test_that("update.brma updates slab labels without fitting", {

  object <- brma.norm(
    yi        = c(0.10, 0.20, 0.05),
    sei       = c(0.05, 0.06, 0.07),
    only_data = TRUE,
    silent    = TRUE
  )

  updated <- update(object, slab = c("A", "B", NA))

  expect_equal(updated[["data"]][["outcome"]][["slab"]], c("A", "B", NA))
  expect_true(attr(updated[["data"]], "slab"))
  expect_null(updated[["fit"]])
})

test_that("update.brma rejects structural update arguments", {

  object <- mock_update_object()

  expect_error(update(object, formula. = ~ x), "Updating formulas")
  expect_error(update(object, evaluate = FALSE), "evaluate = FALSE")
  expect_error(update(object, sample = 10), "Unused argument")
  expect_error(update(object, slab = c("A", "B")), "same as the fitted data")
})

test_that("update.brma extends one sample chunk and recomputes stored caches", {

  object <- mock_update_object()
  object[["loo"]]     <- list(estimate = "old-loo")
  object[["waic"]]    <- list(estimate = "old-waic")
  object[["marglik"]] <- "old-marglik"

  calls <- list()

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      calls[["extend"]]          <<- extend
      calls[["sample_extend"]]   <<- object[["autofit_control"]][["sample_extend"]]
      calls[["max_extend"]]      <<- object[["autofit_control"]][["max_extend"]]
      object[["fit"]][["mock"]] <- "extended"
      return(object[["fit"]])
    },
    .object_summary = function(object) {

      return(data.frame(Mean = 2, row.names = "mu"))
    },
    .object_coefficients = function(object) {

      return(c(mu = 2))
    },
    add_loo = function(object, unit = "estimate", ...) {

      calls[["loo_units"]] <<- c(calls[["loo_units"]], unit)
      if (is.null(object[["loo"]])) {
        object[["loo"]] <- list()
      }
      object[["loo"]][[unit]] <- paste0("new-loo-", unit)
      return(object)
    },
    add_waic = function(object, unit = "estimate", ...) {

      calls[["waic_units"]] <<- c(calls[["waic_units"]], unit)
      if (is.null(object[["waic"]])) {
        object[["waic"]] <- list()
      }
      object[["waic"]][[unit]] <- paste0("new-waic-", unit)
      return(object)
    },
    add_marglik = function(object, ...) {

      calls[["marglik"]] <<- TRUE
      object[["marglik"]] <- "new-marglik"
      return(object)
    },
    .package = "RoBMA"
  )

  updated <- update(
    object,
    sample_extend = 7,
    autofit_control = set_autofit_control(sample_extend = 3, max_extend = 4)
  )

  expect_true(calls[["extend"]])
  expect_equal(calls[["sample_extend"]], 7)
  expect_equal(calls[["max_extend"]], 1)
  expect_equal(updated[["autofit_control"]][["sample_extend"]], 7)
  expect_equal(updated[["autofit_control"]][["max_extend"]], 4)
  expect_equal(updated[["fit"]][["mock"]], "extended")
  expect_equal(updated[["summary"]][["Mean"]], 2)
  expect_equal(updated[["coefficients"]], c(mu = 2))
  expect_equal(updated[["loo"]][["estimate"]], "new-loo-estimate")
  expect_equal(updated[["waic"]][["estimate"]], "new-waic-estimate")
  expect_equal(updated[["marglik"]], "new-marglik")
  expect_equal(calls[["loo_units"]], "estimate")
  expect_equal(calls[["waic_units"]], "estimate")
  expect_true(calls[["marglik"]])
})

test_that("update.brma can drop stored caches after extension", {

  object <- mock_update_object()
  object[["loo"]]     <- list(estimate = "old-loo")
  object[["waic"]]    <- list(estimate = "old-waic")
  object[["marglik"]] <- "old-marglik"

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      object[["fit"]][["mock"]] <- "extended"
      return(object[["fit"]])
    },
    .object_summary = function(object) {

      return(data.frame(Mean = 2, row.names = "mu"))
    },
    .object_coefficients = function(object) {

      return(c(mu = 2))
    },
    .package = "RoBMA"
  )

  updated <- expect_warning(
    update(object, sample_extend = 1, recompute = "drop"),
    "Dropping cached"
  )

  expect_null(updated[["loo"]])
  expect_null(updated[["waic"]])
  expect_null(updated[["marglik"]])
})

test_that("update.brma rejects adversarial extension inputs before fitting", {

  object <- mock_update_object()

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      stop(".fit should not be called", call. = FALSE)
    },
    .package = "RoBMA"
  )

  expect_error(update(object, sample_extend = 0), "sample_extend")
  expect_error(update(object, sample_extend = 1.5), "sample_extend")
  expect_error(update(object, sample_extend = c(1, 2)), "sample_extend")
  expect_error(update(object, sample_extend = NA_integer_), "sample_extend")
  expect_error(update(object, sample_extend = 1, recompute = "none"), "all")
})

test_that("update.brma rejects extension when no fit is stored", {

  object <- mock_update_object()
  object[["fit"]] <- NULL

  expect_error(
    update(object, sample_extend = 1),
    "does not contain a fitted model"
  )
})

test_that("update.brma recomputes all named cached units", {

  object <- mock_update_object()
  object[["loo"]]  <- list(estimate = "old-loo-estimate", cluster = "old-loo-cluster")
  object[["waic"]] <- list(estimate = "old-waic-estimate", cluster = "old-waic-cluster")

  calls <- list()

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      return(object[["fit"]])
    },
    .object_summary = function(object) {

      return(data.frame(Mean = 2, row.names = "mu"))
    },
    .object_coefficients = function(object) {

      return(c(mu = 2))
    },
    add_loo = function(object, unit = "estimate", ...) {

      calls[["loo_units"]] <<- c(calls[["loo_units"]], unit)
      if (is.null(object[["loo"]])) {
        object[["loo"]] <- list()
      }
      object[["loo"]][[unit]] <- paste0("new-loo-", unit)
      return(object)
    },
    add_waic = function(object, unit = "estimate", ...) {

      calls[["waic_units"]] <<- c(calls[["waic_units"]], unit)
      if (is.null(object[["waic"]])) {
        object[["waic"]] <- list()
      }
      object[["waic"]][[unit]] <- paste0("new-waic-", unit)
      return(object)
    },
    .package = "RoBMA"
  )

  updated <- update(object, sample_extend = 1)

  expect_equal(calls[["loo_units"]], c("estimate", "cluster"))
  expect_equal(calls[["waic_units"]], c("estimate", "cluster"))
  expect_equal(updated[["loo"]][["cluster"]], "new-loo-cluster")
  expect_equal(updated[["waic"]][["cluster"]], "new-waic-cluster")
})

test_that("update.brma rejects malformed multi-unit cache names", {

  object <- mock_update_object()
  object[["loo"]] <- list("old-loo-1", "old-loo-2")

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      return(object[["fit"]])
    },
    .object_summary = function(object) {

      return(data.frame(Mean = 2, row.names = "mu"))
    },
    .object_coefficients = function(object) {

      return(c(mu = 2))
    },
    .package = "RoBMA"
  )

  expect_error(
    update(object, sample_extend = 1),
    "must be named by unit"
  )
})

test_that("update.brma does not recompute marginal likelihood for RoBMA objects", {

  object <- mock_update_object()
  class(object) <- c("RoBMA", "brma")
  object[["marglik"]] <- "stale-marglik"

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE) {

      return(object[["fit"]])
    },
    .object_summary = function(object) {

      return(data.frame(Mean = 2, row.names = "mu"))
    },
    .object_coefficients = function(object) {

      return(c(mu = 2))
    },
    add_marglik = function(object, ...) {

      stop("add_marglik should not be called", call. = FALSE)
    },
    .package = "RoBMA"
  )

  updated <- update(object, sample_extend = 1)

  expect_null(updated[["marglik"]])
  expect_s3_class(updated, "RoBMA")
})
