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

if (identical(profile, "all")) {
  clean_cached_fits()
  Sys.setenv(
    ROBMA_TEST_EXTENDED         = "TRUE",
    ROBMA_TEST_FULL_VISUALS    = "TRUE",
    ROBMA_TEST_FULL_DIAGNOSTICS = "TRUE"
  )
  run_tests()
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
