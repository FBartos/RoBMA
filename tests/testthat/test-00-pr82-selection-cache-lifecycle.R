.pr82_cache_restore <- function() {

  options <- RoBMA.options()
  settings <- .selection_runtime_settings(
    capacity_bytes = selection_cache_info()[["capacity_bytes"]][[1L]])
  initialized <- RoBMA.private[["selection_runtime_initialized"]]
  loading <- RoBMA.private[["selection_runtime_loading"]]
  function() {
    RoBMA.private[["selection_runtime_loading"]] <- TRUE
    do.call(RoBMA.options, options)
    .selection_runtime_configure(settings)
    RoBMA.private[["selection_runtime_initialized"]] <- initialized
    RoBMA.private[["selection_runtime_loading"]] <- loading
  }
}

test_that("namespace loading does not probe RAM or warn when automatic probing fails", {

  restore <- .pr82_cache_restore()
  on.exit(restore(), add = TRUE)
  calls <- 0L
  testthat::local_mocked_bindings(ps_system_memory = function() {
    calls <<- calls + 1L
    stop("Unavailable memory probe")
  }, .package = "ps")
  previous_options <- get0(".RoBMA.options", envir = .GlobalEnv, inherits = FALSE)
  had_options <- exists(".RoBMA.options", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_options) assign(".RoBMA.options", previous_options, envir = .GlobalEnv)
    else rm(list = ".RoBMA.options", envir = .GlobalEnv)
  }, add = TRUE)
  assign(".RoBMA.options", list(selection.cache_max_bytes = "auto"), envir = .GlobalEnv)
  expect_warning(.onLoad(RoBMA.private[["lib_name"]], "RoBMA"), NA)
  expect_identical(calls, 0L)
  expect_false(RoBMA.private[["selection_runtime_initialized"]])
  expect_identical(RoBMA.get_option("selection.cache_max_bytes"), "auto")
  expect_equal(selection_cache_info()[["capacity_bytes"]][[1L]], 0)
  expect_warning(.selection_runtime_ensure(), "Automatic selection cache memory budget is unavailable", fixed = TRUE)
  expect_warning(.selection_runtime_ensure(), NA)
  expect_identical(calls, 1L)
})

test_that("saved selection contexts initialize a serial process cache exactly once", {

  restore <- .pr82_cache_restore()
  on.exit(restore(), add = TRUE)
  object <- bselmodel(yi = c(-.1, .2), sei = c(.3, .4), measure = "SMD",
                       prior_bias = BayesTools::prior_weightfunction(
                         "one-sided", .5, BayesTools::wf_fixed(c(1, .4))),
                       only_priors = TRUE, silent = TRUE)
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(list(object = object, samples = cbind(mu = c(0, .1), tau = c(.2, .3))), saved)
  restored <- readRDS(saved)
  calls <- 0L
  testthat::local_mocked_bindings(ps_system_memory = function() {
    calls <<- calls + 1L
    list(avail = 16 * 1024^2)
  }, .package = "ps")
  RoBMA.private[["selection.cache_max_bytes"]] <- "auto"
  .selection_cache_configure(0)
  RoBMA.private[["selection_runtime_initialized"]] <- FALSE
  score <- function() .log_lik_from_posterior_samples_sum(
    restored[["object"]][["fit"]], restored[["samples"]],
    restored[["object"]][["data"]], restored[["object"]][["priors"]]
  )
  first <- score()
  expect_equal(selection_cache_info()[["capacity_bytes"]][[1L]], 4 * 1024^2)
  expect_true(RoBMA.private[["selection_runtime_initialized"]])
  second <- score()
  expect_identical(first, second)
  expect_true(all(is.finite(first)))
  expect_identical(calls, 1L)
})

test_that("lazy initialization preserves explicitly assigned worker and disabled budgets", {

  restore <- .pr82_cache_restore()
  on.exit(restore(), add = TRUE)
  testthat::local_mocked_bindings(ps_system_memory = function() stop("Unexpected probe"),
                                 .package = "ps")
  settings <- .selection_runtime_settings(capacity_bytes = 1024^2)
  .selection_runtime_configure(settings,
    list(phase = "start", role = "worker", parallel = TRUE, chains = 4L, process_chains = 1L))
  expect_silent(.selection_runtime_ensure())
  expect_equal(selection_cache_info()[["capacity_bytes"]][[1L]], 1024^2 / 4)
  RoBMA.options(selection.cache_max_bytes = 0)
  expect_silent(.selection_runtime_ensure())
  expect_equal(selection_cache_info()[["capacity_bytes"]][[1L]], 0)
})
