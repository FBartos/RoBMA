test_that("scalar Normal intervals use the appropriate shared tail", {

  reference <- readRDS(testthat::test_path("fixtures", "normal-tail-reference.rds"))
  points <- reference[["intervals"]]
  expect_identical(nrow(points), 32L)
  for (i in seq_len(nrow(points))) {
    lower <- as.numeric(points[["lower_hex"]][[i]])
    upper <- as.numeric(points[["upper_hex"]][[i]])
    mean <- as.numeric(points[["mean_hex"]][[i]])
    sd <- as.numeric(points[["sd_hex"]][[i]])
    actual <- .iwmde_log_normal_interval_prob(lower, upper, mean, sd)
    expect_identical(actual, .iwmde_qcmde_log_normal_mass(lower, upper, mean, sd))
    oracle <- reference[["oracles"]][[i]]
    allowance <- 64 * .Machine$double.eps * max(1, abs(oracle))
    if (startsWith(points[["case"]][[i]], "legacy169-")) {
      allowance <- allowance + .ordered_normal_interval_arithmetic_allowance(lower, upper, mean, sd)[["log_error"]]
    }
    expect_lte(abs(actual - oracle), allowance)
    if (i <= 16L) {
      z <- as.numeric(points[["z_hex"]][[i]])
      expect_identical(actual, stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
    }
  }
  expect_identical(.iwmde_normal_interval_conditioning_scale(-Inf, Inf, 0, 1), 0)
  expect_identical(.iwmde_normal_interval_conditioning_scale(0, Inf, 1, 1e-320), Inf)
  expect_identical(.iwmde_normal_interval_conditioning_scale(1, 1, 0, 1), Inf)
})


test_that("scalar Normal interval edge and shape behavior remains unchanged", {

  reference <- readRDS(testthat::test_path("fixtures", "normal-tail-reference.rds"))[["edges"]]
  expect_length(reference[["inputs"]], 14L)
  capture <- function(arguments) {
    warnings <- list()
    value <- tryCatch(withCallingHandlers(do.call(.iwmde_log_normal_interval_prob, arguments), warning = function(warning) {
      warnings[[length(warnings) + 1L]] <<- list(message = conditionMessage(warning), class = class(warning))
    }), error = function(error) structure(list(message = conditionMessage(error), class = class(error)), class = "captured_error"))
    list(value = if (inherits(value, "captured_error")) NULL else value,
      error = if (inherits(value, "captured_error")) unclass(value) else NULL, warnings = warnings)
  }
  for (name in names(reference[["inputs"]])) {
    actual <- capture(reference[["inputs"]][[name]])
    expected <- reference[["results"]][[name]][["original"]]
    if (is.null(expected[["error"]])) {
      expect_identical(actual, expected)
    } else {
      # Base-R condition wording is localized by gettext; preserve the error
      # behavior and family without prescribing the captured English text.
      expect_false(is.null(actual[["error"]]))
      expect_null(actual[["value"]])
      expect_identical(actual[["error"]][["class"]], expected[["error"]][["class"]])
      expect_identical(actual[["warnings"]], expected[["warnings"]])
    }
  }
})
