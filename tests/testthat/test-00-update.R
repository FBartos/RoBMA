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
    .fit = function(object, extend = FALSE, adapt = NULL) {

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

test_that("update.brma extends by 'auto' or a number, and not at all without sample_extend", {

  object <- mock_update_object()
  expect_identical(object[["autofit_control"]][["sample_extend"]], "auto")
  # a fit of the earlier layout stores a fixed step
  object[["autofit_control"]][["sample_extend"]] <- 7

  calls <- list()
  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE, adapt = NULL) {
      calls[[length(calls) + 1L]] <<- object[["autofit_control"]][c("sample_extend", "max_extend")]
      return(object[["fit"]])
    },
    .object_summary = function(object) data.frame(Mean = 2, row.names = "mu"),
    .object_coefficients = function(object) c(mu = 2),
    .package = "RoBMA"
  )

  auto   <- update(object, sample_extend = "auto")
  fixed  <- update(object, sample_extend = 30)
  absent <- update(object, autofit_control = list(max_extend = 3))

  # one extension each, sized by BayesTools for "auto"
  expect_equal(calls, list(
    list(sample_extend = "auto", max_extend = 1),
    list(sample_extend = 30, max_extend = 1)
  ))
  expect_identical(auto[["autofit_control"]][["sample_extend"]], "auto")
  expect_identical(fixed[["autofit_control"]][["sample_extend"]], 30)
  # without sample_extend, update() does not extend and keeps the stored step
  expect_identical(absent[["autofit_control"]][["sample_extend"]], 7)
  expect_identical(absent[["autofit_control"]][["max_extend"]], 3)

  for (value in list("Auto", "", c("auto", "auto"), NA_character_)) {
    expect_error(update(object, sample_extend = value), "sample_extend", info = deparse(value))
  }
  expect_length(calls, 2L)
})

test_that("update.brma can drop stored caches after extension", {

  object <- mock_update_object()
  object[["loo"]]     <- list(estimate = "old-loo")
  object[["waic"]]    <- list(estimate = "old-waic")
  object[["marglik"]] <- "old-marglik"

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE, adapt = NULL) {

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
    .fit = function(object, extend = FALSE, adapt = NULL) {

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


test_that("update.brma preserves extension fitting errors", {

  object         <- mock_update_object()
  summary_called <- FALSE

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE, adapt = NULL) {

      return(list(
        errors        = "Original extension fitting diagnostic.",
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
    update(object, sample_extend = 1),
    "^Original extension fitting diagnostic[.]$"
  )
  expect_false(summary_called)
})

test_that("update.brma recomputes all named cached units", {

  object <- mock_update_object()
  object[["loo"]]  <- list(estimate = "old-loo-estimate", cluster = "old-loo-cluster")
  object[["waic"]] <- list(estimate = "old-waic-estimate", cluster = "old-waic-cluster")

  calls <- list()

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE, adapt = NULL) {

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
    .fit = function(object, extend = FALSE, adapt = NULL) {

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
    .fit = function(object, extend = FALSE, adapt = NULL) {

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

test_that("update.brma applies an 'adapt' override to its extension only", {

  object <- mock_update_object()
  object[["fit_control"]][["adapt"]] <- 5000
  forwarded <- list()

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE, adapt = NULL) {

      forwarded[[length(forwarded) + 1L]] <<- adapt
      return(object[["fit"]])
    },
    .object_summary = function(object) data.frame(Mean = 2, row.names = "mu"),
    .object_coefficients = function(object) c(mu = 2),
    .package = "RoBMA"
  )

  omitted  <- update(object, sample_extend = 10)
  explicit <- update(object, sample_extend = 10, adapt = 5000)
  override <- update(object, sample_extend = 10, adapt = 20000)
  after    <- update(override, sample_extend = 10)

  # an omitted value is the initial setting; an override is passed for its
  # extension and leaves the stored initial setting unchanged
  expect_identical(forwarded, list(5000, 5000, 20000, 5000))
  expect_identical(omitted, explicit)
  for (updated in list(omitted, override, after)) {
    expect_identical(updated[["fit_control"]], object[["fit_control"]])
  }
})

test_that("update.brma validates 'adapt' as the fitting functions do", {

  object <- mock_update_object()

  testthat::local_mocked_bindings(
    .fit = function(object, extend = FALSE, adapt = NULL) {

      stop(".fit should not be called", call. = FALSE)
    },
    .package = "RoBMA"
  )

  for (adapt in list(49, 0, 50.5, NA, Inf, "100", c(100, 200), numeric(0), .Machine$integer.max + 1)) {
    fitting_error <- tryCatch(
      brma.norm(yi = c(0.1, 0.2), sei = c(0.05, 0.06), measure = "GEN", adapt = adapt, only_data = TRUE),
      error = conditionMessage
    )
    expect_true(is.character(fitting_error), info = deparse(adapt))
    expect_identical(
      tryCatch(update(object, sample_extend = 1, adapt = adapt), error = conditionMessage),
      fitting_error,
      info = deparse(adapt)
    )
  }
  expect_error(
    update(object, adapt = 1000),
    "The 'adapt' argument applies only to an extension; supply 'sample_extend'.",
    fixed = TRUE
  )
})

# A tiny normal model fitted with a non-default adaptation, and the
# settings that reach the direct JAGS driver while extending it.
.adaptation_fit <- function(..., constructor = brma.norm) {

  suppressWarnings(constructor(
    yi = c(0.10, 0.25, 0.05, 0.32, 0.18), sei = c(0.08, 0.10, 0.07, 0.12, 0.09),
    measure = "GEN", prior_unit_information_sd = 1, chains = 2, adapt = 300, burnin = 100, sample = 100,
    seed = 11, silent = TRUE, ...
  ))
}

.record_jags_runs <- function(env = parent.frame()) {

  calls <- list()
  original <- getFromNamespace(".bt_rjags_run", "BayesTools")
  testthat::local_mocked_bindings(
    .bt_rjags_run = function(model, data, inits, monitor, adapt, burnin, sample,
                            thin, cl, silent, live_model = NULL, probability_plan = NULL) {
      calls[[length(calls) + 1L]] <<- list(
        adapt = adapt, burnin = burnin, parallel = !is.null(cl), probability_plan = probability_plan
      )
      original(model, data, inits, monitor, adapt, burnin, sample, thin, cl, silent,
               live_model = live_model, probability_plan = probability_plan)
    },
    .package = "BayesTools",
    .env = env
  )
  function() calls
}

.draws_by_chain <- function(object) {

  lapply(coda::as.mcmc.list(object[["fit"]]), function(chain) {
    matrix(as.numeric(chain), nrow = nrow(chain), dimnames = list(NULL, colnames(chain)))
  })
}

test_that("extensions inherit the initial adaptation, preserve the draws, and add no burn-in", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(sampler.probabilities = character())
  fit   <- .adaptation_fit()
  calls <- .record_jags_runs()

  omitted  <- update(fit, sample_extend = 20)
  explicit <- update(fit, sample_extend = 20, adapt = 300)
  override <- update(fit, sample_extend = 20, adapt = 1000)
  after    <- update(override, sample_extend = 20)

  # the arguments that reached the driver
  expect_identical(
    lapply(calls(), `[[`, "adapt"),
    list(300L, 300L, 1000L, 300L)
  )
  expect_identical(lapply(calls(), `[[`, "burnin"), rep(list(0L), 4))
  expect_identical(lapply(calls(), `[[`, "parallel"), rep(list(FALSE), 4))
  expect_identical(lapply(calls(), `[[`, "probability_plan"), rep(list(NULL), 4))

  # an explicit value equal to the initial one gives identical draws, an
  # override other draws; no extension changes the stored setting
  expect_identical(.draws_by_chain(omitted), .draws_by_chain(explicit))
  expect_false(identical(.draws_by_chain(omitted), .draws_by_chain(override)))
  expect_identical(fit[["fit_control"]][["adapt"]], 300)
  for (updated in list(omitted, override, after)) {
    expect_identical(updated[["fit_control"]], fit[["fit_control"]])
  }
  expect_identical(
    attr(after[["fit"]], "adaptation"),
    list(initial = 300L, extensions = data.frame(adapt = c(1000L, 300L), rebuilt = c(TRUE, TRUE)))
  )

  # every chain keeps its draws and gains sample_extend new ones
  old_draws <- .draws_by_chain(fit)
  for (updated in list(omitted, override)) {
    new_draws <- .draws_by_chain(updated)
    expect_identical(vapply(new_draws, nrow, integer(1)), vapply(old_draws, nrow, integer(1)) + 20L)
    expect_identical(lapply(new_draws, function(draws) draws[seq_len(100), , drop = FALSE]), old_draws)
  }
  expect_identical(
    lapply(.draws_by_chain(after), function(draws) draws[seq_len(120), , drop = FALSE]),
    .draws_by_chain(override)
  )

  # a saved and reloaded fit extends as the original
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(fit, saved)
  expect_identical(.draws_by_chain(update(readRDS(saved), sample_extend = 20)), .draws_by_chain(omitted))
  expect_identical(calls()[[5L]][["adapt"]], 300L)
})

test_that("extensions retain the probability plan captured at the initial fit", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  RoBMA.options(selection.sampler = "slice", sampler.joint = "off",
                sampler.probabilities = "mu_indicator")
  calls <- .record_jags_runs()
  fit <- .adaptation_fit(constructor = RoBMA, model_type = "PP")
  plan <- attr(fit[["fit"]], "sampler_declarations", exact = TRUE)[["probabilities"]]
  expect_identical(plan[["node"]], "mu_indicator")
  before <- fit[["fit"]][["sampler_probabilities"]][["mu_indicator"]]
  expect_s3_class(before, "mcmc.list")

  RoBMA.options(sampler.probabilities = character())
  extended <- update(fit, sample_extend = 20)
  expect_identical(lapply(calls(), `[[`, "probability_plan"), rep(list(plan), 2L))

  after <- extended[["fit"]][["sampler_probabilities"]][["mu_indicator"]]
  expect_s3_class(after, "mcmc.list")
  expect_identical(vapply(after, nrow, integer(1)), c(120L, 120L))
  expect_identical(
    lapply(after, function(chain) as.matrix(chain)[seq_len(100L), , drop = FALSE]),
    lapply(before, as.matrix)
  )
})

test_that("automatic extensions of a fit use its adaptation", {

  calls <- .record_jags_runs()
  fit <- .adaptation_fit(
    autofit = TRUE,
    autofit_control = set_autofit_control(min_ESS = 1e6, sample_extend = 10, max_extend = 2)
  )

  backend <- calls()
  expect_equal(backend[[1L]][c("adapt", "burnin")], list(adapt = 300, burnin = 100))
  expect_identical(lapply(backend[-1L], `[[`, "adapt"), list(300L, 300L))
  expect_identical(lapply(backend[-1L], `[[`, "burnin"), list(0L, 0L))
  # serial automatic extensions continue the sampler adapted by the initial run
  expect_identical(
    attr(fit[["fit"]], "adaptation"),
    list(initial = 300L, extensions = data.frame(adapt = c(300L, 300L), rebuilt = c(FALSE, FALSE)))
  )
  expect_identical(vapply(.draws_by_chain(fit), nrow, integer(1)), c(120L, 120L))
})

test_that("update.brma sizes an 'auto' extension from the convergence check of the stored fit", {

  fit <- .adaptation_fit()
  # 100 draws per chain: twice the draws is below the minimum extension of 1000
  extended <- update(fit, sample_extend = "auto")
  expect_identical(vapply(.draws_by_chain(extended), nrow, integer(1)), c(1100L, 1100L))
  expect_identical(
    lapply(.draws_by_chain(extended), function(draws) draws[seq_len(100), , drop = FALSE]),
    .draws_by_chain(fit)
  )
  expect_identical(extended[["autofit_control"]][["sample_extend"]], "auto")

  # an unreachable ESS target is capped at twice the draws per chain, and
  # criteria that hold give the minimum
  capped <- update(extended, sample_extend = "auto", autofit_control = list(min_ESS = 1e6))
  expect_identical(vapply(.draws_by_chain(capped), nrow, integer(1)), c(3300L, 3300L))
  met <- update(extended, sample_extend = "auto", autofit_control = list(max_Rhat = NULL, min_ESS = NULL))
  expect_identical(vapply(.draws_by_chain(met), nrow, integer(1)), c(2100L, 2100L))

  # a number extends by that number
  expect_identical(vapply(.draws_by_chain(update(fit, sample_extend = 20)), nrow, integer(1)), c(120L, 120L))
})
