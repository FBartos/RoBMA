#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) > 0L) args[[1L]] else "all"
filters <- if (length(args) > 1L) args[-1L] else character()

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "test-profile.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

Sys.setenv(AGENT = "1")

source(file.path("tests", "testthat", "common-functions.R"))

set_full_test_env <- function() {

  Sys.setenv(
    ROBMA_TEST_EXTENDED         = "TRUE",
    ROBMA_TEST_FULL_VISUALS    = "TRUE",
    ROBMA_TEST_FULL_DIAGNOSTICS = "TRUE"
  )

  return(invisible(TRUE))
}


run_tests <- function(filter = NULL) {

  test_args <- list(
    reporter        = "llm",
    stop_on_failure = TRUE
  )
  if (!is.null(filter) && nzchar(filter)) {
    test_args[["filter"]] <- filter
  }

  do.call(devtools::test, test_args)
}


validate_fit_cache <- function() {

  expected <- fit_catalog()[["name"]]
  available <- list_fits(validate = TRUE)
  missing <- setdiff(expected, available)

  if (length(missing) > 0L) {
    stop(
      "Cached fit generation incomplete. Missing or stale fits: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  message("Validated ", length(available), " cached fits.")

  return(invisible(available))
}


run_release_profile <- function() {

  set_full_test_env()
  clean_cached_fits()
  run_tests("01-")
  validate_fit_cache()
  Sys.setenv(
    ROBMA_TEST_FORCE_REFIT = "FALSE",
    ROBMA_TEST_SKIP_REFIT  = "TRUE"
  )
  run_tests()
  devtools::check(error_on = "warning")

  return(invisible(TRUE))
}


if (identical(profile, "all")) {
  clean_cached_fits()
  set_full_test_env()
  run_tests()
} else if (profile %in% c("release", "pre-release", "pre_release")) {
  run_release_profile()
} else if (profile %in% c("cache", "01-cache")) {
  clean_cached_fits()
  run_tests("01-")
} else if (identical(profile, "quick")) {
  run_tests("00-|02-")
} else if (identical(profile, "filter")) {
  if (length(filters) == 0L) {
    stop("Profile 'filter' requires at least one testthat filter.", call. = FALSE)
  }
  run_tests(paste(filters, collapse = "|"))
} else {
  run_tests(profile)
}
