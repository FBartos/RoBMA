#!/usr/bin/env Rscript

args       <- commandArgs(trailingOnly = TRUE)
profile    <- if (length(args) > 0L) args[[1L]] else "standard"
clean      <- "--clean" %in% args
extra_args <- setdiff(args[-1L], "--clean")

cmd      <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "test-profile.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
script_path  <- normalizePath(script_path, mustWork = TRUE)
setwd(project_root)

Sys.setenv(AGENT = "1")

run_subprofile <- function(name) {

  rscript <- file.path(R.home("bin"), "Rscript.exe")
  status  <- system2(
    rscript,
    args = c(shQuote(script_path), name, "--clean")
  )

  if (!identical(status, 0L)) {
    stop("Test profile failed: ", name, call. = FALSE)
  }

  return(invisible(TRUE))
}

if (profile %in% c("release", "pre-release", "pre_release")) {
  run_subprofile("standard")
  run_subprofile("certification")
  devtools::check(error_on = "warning")
  quit(save = "no", status = 0L)
}

if (identical(profile, "all")) {
  profile <- "standard"
  clean   <- TRUE
}

filter <- NULL
if (profile %in% c("cache", "01-cache")) {
  profile <- "standard"
  filter  <- "01-"
  clean   <- TRUE
} else if (identical(profile, "quick")) {
  profile <- "standard"
  filter  <- "00-|02-"
} else if (identical(profile, "filter")) {
  if (length(extra_args) == 0L) {
    stop("Profile 'filter' requires at least one testthat filter.", call. = FALSE)
  }
  profile <- "standard"
  filter  <- paste(extra_args, collapse = "|")
} else if (profile %in% c("standard", "certification") &&
           length(extra_args) > 0L) {
  filter <- paste(extra_args, collapse = "|")
} else if (!profile %in% c("standard", "certification")) {
  filter  <- profile
  profile <- "standard"
}

Sys.setenv(ROBMA_TEST_PROFILE = profile)
source(file.path("tests", "testthat", "common-functions.R"))

if (clean) {
  clean_cached_fits()
}

started <- proc.time()[["elapsed"]]
run_tests <- function(filter = NULL) {

  test_args <- list(
    reporter        = "llm",
    stop_on_failure = TRUE
  )
  if (!is.null(filter)) {
    test_args[["filter"]] <- filter
  }

  do.call(devtools::test, test_args)
}

validate_fit_cache <- function() {

  expected  <- active_fit_catalog()[["name"]]
  available <- list_fits(validate = TRUE)
  missing   <- setdiff(expected, available)
  if (length(missing) > 0L) {
    stop(
      "Cached fit generation incomplete. Missing or stale fits: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  message("Validated ", length(available), " ", profile, " cached fits.")

  return(invisible(available))
}

if (is.null(filter)) {
  run_tests("01-")
  validate_fit_cache()
  run_tests("00-|02-|03-")
} else {
  run_tests(filter)
  if (identical(filter, "01-")) {
    validate_fit_cache()
  }
}

elapsed <- proc.time()[["elapsed"]] - started
message(profile, " profile completed in ", round(elapsed, 1), " seconds.")

if (identical(profile, "standard") && is.null(filter) && elapsed > 15 * 60) {
  stop(
    "The standard test profile exceeded its 15-minute runtime budget (",
    round(elapsed, 1), " seconds).",
    call. = FALSE
  )
}
