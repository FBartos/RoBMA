test_that("cache retention is opt-in and removal only drops the object's reference", {

  old_options <- RoBMA.options()
  on.exit(do.call(RoBMA.options, old_options), add = TRUE)
  expect_false(.RoBMA_option_schema[["selection.cache_retain"]][["default"]])
  RoBMA.options(selection.cache_retain = FALSE)
  expect_null(.selection_cache_runtime())
  error <- tryCatch(RoBMA.options(selection.cache_retain = 1), error = identity)
  expect_identical(conditionMessage(error),
    "Option 'selection.cache_retain' must be TRUE or FALSE.")
  expect_null(conditionCall(error))

  fit <- list(draws = matrix(seq_len(6L), ncol = 2L))
  attr(fit, "runtime_state") <- list(list(payload = as.raw(1:16)))
  object <- structure(list(fit = fit, summary = "unchanged"), class = "brma")
  before <- selection_cache_info()
  cleared <- remove_selection_cache(object)
  expect_null(attr(cleared[["fit"]], "runtime_state", exact = TRUE))
  expect_identical(cleared[["fit"]][["draws"]], object[["fit"]][["draws"]])
  expect_identical(cleared[["summary"]], object[["summary"]])
  expect_identical(attr(object[["fit"]], "runtime_state"), attr(fit, "runtime_state"))
  expect_identical(selection_cache_info(), before)
  expect_identical(remove_selection_cache(cleared), cleared)
  error <- tryCatch(remove_selection_cache(list()), error = identity)
  expect_identical(conditionMessage(error), "'object' must be a 'brma' object.")
  expect_null(conditionCall(error))
})

test_that("retained cache callbacks round-trip and respect the current retention option", {

  old_options <- RoBMA.options()
  on.exit({
    do.call(RoBMA.options, old_options)
    selection_cache_info(clear = TRUE)
  }, add = TRUE)
  current_build <- "test-build"
  testthat::local_mocked_bindings(.selection_cache_build = function() current_build)
  RoBMA.options(selection.cache_max_bytes = 1024^2, selection.cache_retain = TRUE)
  callback <- .selection_cache_runtime()
  # A serialized callback contains settings, never retained bytes or an old
  # build identity that could validate stale snapshots after an upgrade.
  callback <- unserialize(serialize(callback, NULL))
  context <- list(phase = "restore")
  callback(context)
  sei <- c(.2, .3, .4)
  covariance <- .4 * outer(sei, sei)
  diag(covariance) <- sei^2
  quadrature <- .selection_joint_cluster_quadrature_rules(c(7L, 15L, 31L, 63L))
  arguments <- list(rep(0, 3L), matrix(rep(0, 3L), 1L),
    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L), sei,
    matrix(c(1, .2), 1L), c(0, -Inf), c(Inf, 0), rep(1L, 3L),
    1L, TRUE, 1L,
    as.double(BayesTools::selection_qmc_design(6L, 512L, 8L, seed = 173L)),
    512L, 8L, .005, TRUE, 0L, quadrature)
  evaluate <- function() {

    do.call(.Call, c(list("RoBMA_selnorm_mnorm_step_loglik_batch"), arguments,
      list(PACKAGE = "RoBMA")))
  }
  first <- evaluate()
  expect_gt(selection_cache_info()$entries[1L], 0)
  snapshot <- callback(list(phase = "capture"))
  expect_true(is.raw(snapshot[["payload"]]))
  expect_gt(length(snapshot[["payload"]]), 0)
  expect_equal(selection_cache_info()$allocated_bytes[1L], 0)
  snapshot <- unserialize(serialize(snapshot, NULL))
  set.seed(16)
  rng <- .Random.seed
  callback(context, list(snapshot))
  expect_identical(.Random.seed, rng)
  expect_identical(evaluate(), first)
  expect_equal(selection_cache_info()$hits[1L], 1)

  saved_fit <- structure(list(), runtime_state = list(snapshot))
  RoBMA.options(selection.cache_retain = FALSE)
  callback <- .selection_cache_runtime(saved_fit)
  callback(context, list(snapshot))
  expect_identical(evaluate(), first)
  expect_equal(selection_cache_info()$hits[1L], 1)
  expect_null(callback(list(phase = "capture")))
  expect_gt(selection_cache_info()$entries[1L], 0)

  current_build <- "different-build"
  warning <- NULL
  withCallingHandlers(callback(context, list(snapshot)), warning = function(w) {
    warning <<- w
    invokeRestart("muffleWarning")
  })
  expect_identical(conditionMessage(warning), paste0(
    "Retained selection cache is unavailable with mismatched or missing build metadata. ",
    "Extending with an empty cache for those entries."))
  expect_null(conditionCall(warning))
  expect_equal(selection_cache_info()$entries[1L], 0)
  expect_identical(evaluate(), first)

  current_build <- "test-build"
  corrupted <- snapshot
  last <- length(corrupted[["payload"]])
  corrupted[["payload"]][last] <- as.raw(bitwXor(as.integer(corrupted[["payload"]][last]), 1L))
  warning <- NULL
  withCallingHandlers(callback(context, list(corrupted)), warning = function(w) {
    warning <<- w
    invokeRestart("muffleWarning")
  })
  expect_identical(conditionMessage(warning), paste0(
    "Retained selection cache is unavailable: Selection cache snapshot checksum does not match. ",
    "Extending with an empty cache."))
  expect_null(conditionCall(warning))
  expect_equal(selection_cache_info()$entries[1L], 0)
  expect_identical(evaluate(), first)
  .selection_cache_configure(0)
  callback(context, list(snapshot))
  expect_equal(selection_cache_info()$allocated_bytes[1L], 0)
  expect_identical(evaluate(), first)
})
