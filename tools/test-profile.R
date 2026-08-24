#!/usr/bin/env Rscript

args       <- commandArgs(trailingOnly = TRUE)
profile    <- if (length(args) > 0L) args[[1L]] else "standard"
clean      <- "--clean" %in% args
list_only  <- "--list" %in% args
worker_arg <- grep("^--case-worker=", args, value = TRUE)
extra_args <- setdiff(args[-1L], c("--clean", "--list", worker_arg))

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

if (Sys.getenv("AGENT") == "" && Sys.getenv("ROBMA_TEST_REPORTER") == "") {
  Sys.setenv(AGENT = "1")
}

message("Running RoBMA test profile: ", profile)
message(
  "ROBMA_TEST_FILES_DIR: ",
  Sys.getenv("ROBMA_TEST_FILES_DIR", unset = "<package default>")
)
message("ROBMA_TEST_REPORTER: ", Sys.getenv("ROBMA_TEST_REPORTER", unset = "<unset>"))
message("AGENT: ", Sys.getenv("AGENT", unset = "<unset>"))

profile_started <- Sys.time()

finish_profile <- function(status = 0L) {

  elapsed <- difftime(Sys.time(), profile_started, units = "secs")
  message("Elapsed seconds: ", round(as.numeric(elapsed), 1))
  quit(save = "no", status = status)
}


.rscript <- function() {

  command <- Sys.which("Rscript")
  if (!nzchar(command)) {
    command <- file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
  }

  return(command)
}


run_subprofile <- function(name, clean = FALSE) {

  subprofile_args <- c(shQuote(script_path), name)
  if (clean) {
    subprofile_args <- c(subprofile_args, "--clean")
  }
  status <- system2(.rscript(), args = subprofile_args)
  if (!identical(status, 0L)) {
    stop("Test profile failed: ", name, call. = FALSE)
  }

  return(invisible(TRUE))
}


run_tests <- function(filter = NULL) {

  reporter <- Sys.getenv("ROBMA_TEST_REPORTER", unset = "llm")
  stop_on_failure <- Sys.getenv("ROBMA_TEST_STOP_ON_FAILURE", unset = "")
  stop_on_failure <- if (!nzchar(stop_on_failure)) {
    TRUE
  } else {
    is_true_env("ROBMA_TEST_STOP_ON_FAILURE")
  }
  test_args <- list(
    reporter        = if (identical(reporter, "llm") &&
                           is_true_env("ROBMA_TEST_QUIET_SKIPS")) {
      quiet_llm_reporter()
    } else {
      reporter
    },
    stop_on_failure = stop_on_failure
  )
  if (!is.null(filter)) {
    test_args[["filter"]] <- filter
  }

  message(
    "testthat filter: ",
    if (is.null(filter)) "<all>" else filter
  )
  message("testthat reporter: ", reporter)

  return(do.call(devtools::test, test_args))
}


validate_fit_cache <- function(label) {

  expected  <- active_fit_catalog()[["name"]]
  available <- list_fits(validate = TRUE)
  missing   <- setdiff(expected, available)
  if (length(missing) > 0L) {
    stop(
      "The ", label, " cache is missing or stale: ",
      paste(missing, collapse = ", "),
      ". Run the matching cache refresh before this profile.",
      call. = FALSE
    )
  }
  message("Validated ", length(available), " ", label, " cached fits.")

  return(invisible(available))
}


set_active_fits <- function(fit_names) {

  value <- if (length(fit_names) == 0L) {
    "__none__"
  } else {
    paste(fit_names, collapse = ",")
  }
  Sys.setenv(ROBMA_TEST_ACTIVE_FITS = value)

  return(invisible(TRUE))
}


source_profile_helpers <- function(profile_name, active_fits = NULL) {

  Sys.setenv(ROBMA_TEST_PROFILE = profile_name)
  if (!is.null(active_fits)) {
    set_active_fits(active_fits)
  } else {
    Sys.unsetenv("ROBMA_TEST_ACTIVE_FITS")
  }
  source(file.path("tests", "testthat", "common-functions.R"), local = .GlobalEnv)
  source(
    file.path("tests", "testthat", "helper-profile-cases.R"),
    local = .GlobalEnv
  )

  message("Cached-fit profile: ", profile_name)
  message("Cached-fit directory: ", test_files_dir)

  return(invisible(TRUE))
}


run_certification_worker <- function(name) {

  source_profile_helpers("certification")
  fit_names <- certification_case_fit_names(name)
  set_active_fits(fit_names)
  case    <- certification_case(name)
  catalog <- fit_catalog()

  if (clean && length(fit_names) > 0L) {
    clean_cached_fits(fit_names)
  }

  started <- proc.time()[["elapsed"]]
  for (source_file in case[["fit_sources"]]) {
    source_fit_names <- fit_names[
      fit_names %in% catalog[["name"]][catalog[["source_file"]] == source_file]
    ]
    if (length(source_fit_names) == 0L) {
      next
    }

    set_active_fits(source_fit_names)
    source_filter <- sub("^test-", "", source_file)
    source_filter <- sub("\\.[Rr]$", "", source_filter)
    run_tests(source_filter)
    validate_fit_cache(paste0(
      "certification case '", name, "' source '", source_file, "'"
    ))
  }

  set_active_fits(fit_names)
  results <- run_tests(case[["test_filter"]])
  validate_certification_evidence(
    results        = results,
    required_tests = case[["required_tests"]],
    case_name      = name
  )

  elapsed <- proc.time()[["elapsed"]] - started
  message(
    "Certification case '", name, "' completed in ",
    round(elapsed, 1), " seconds."
  )

  return(invisible(TRUE))
}


run_certification_case <- function(name, clean = FALSE) {

  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("Package 'processx' is required to enforce certification timeouts.",
         call. = FALSE)
  }

  case_args <- c(
    script_path,
    "certification",
    paste0("--case-worker=", name)
  )
  if (clean) {
    case_args <- c(case_args, "--clean")
  }

  message("Starting certification case '", name, "' (one-hour limit).")
  result <- processx::run(
    command         = .rscript(),
    args            = case_args,
    timeout         = CERTIFICATION_CASE_TIMEOUT_SECONDS,
    echo            = TRUE,
    echo_cmd        = TRUE,
    spinner         = FALSE,
    error_on_status = FALSE
  )
  if (isTRUE(result[["timeout"]])) {
    stop(
      "Certification case '", name, "' exceeded its one-hour limit.",
      call. = FALSE
    )
  }
  if (!identical(result[["status"]], 0L)) {
    stop(
      "Certification case '", name, "' failed with status ",
      result[["status"]], ".",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


if (profile %in% c("release", "pre-release", "pre_release")) {
  run_subprofile("refresh-standard", clean = clean)
  run_subprofile("standard")
  run_subprofile("certification", clean = clean)
  devtools::check(error_on = "warning")
  finish_profile()
}

if (profile %in% c("cache", "01-cache", "refresh-standard")) {
  source_profile_helpers("standard")
  if (clean) {
    clean_cached_fits()
  }
  run_tests("01-")
  validate_fit_cache("standard")
  finish_profile()
}

if (identical(profile, "certification")) {
  source_profile_helpers("certification")
  if (list_only) {
    cases <- certification_cases()
    for (name in names(cases)) {
      cat(name, "\t", cases[[name]][["description"]], "\n", sep = "")
    }
    finish_profile()
  }
  if (length(worker_arg) > 1L) {
    stop("Only one certification worker may be requested.", call. = FALSE)
  }
  if (length(worker_arg) == 1L) {
    name <- sub("^--case-worker=", "", worker_arg)
    run_certification_worker(name)
    finish_profile()
  }

  requested <- if (length(extra_args) == 0L) {
    certification_case_names()
  } else {
    extra_args
  }
  unknown <- setdiff(requested, certification_case_names())
  if (length(unknown) > 0L) {
    stop(
      "Unknown certification case: ", paste(unknown, collapse = ", "),
      ". Use '--list' to list available cases.",
      call. = FALSE
    )
  }
  for (name in requested) {
    run_certification_case(name, clean = clean)
  }
  finish_profile()
}

if (identical(profile, "standard")) {
  if (clean) {
    stop(
      "The standard profile never deletes or creates fits. Use ",
      "'refresh-standard --clean' first.",
      call. = FALSE
    )
  }
  source_profile_helpers("standard")
  if (length(extra_args) > 0L) {
    run_tests(paste(extra_args, collapse = "|"))
    finish_profile()
  }

  validate_fit_cache("standard")
  started <- proc.time()[["elapsed"]]
  run_tests("00-|02-|03-")
  elapsed <- proc.time()[["elapsed"]] - started
  message("Standard profile completed in ", round(elapsed, 1), " seconds.")
  if (elapsed > 15 * 60) {
    stop(
      "The standard test profile exceeded its 15-minute runtime budget (",
      round(elapsed, 1), " seconds).",
      call. = FALSE
    )
  }
  finish_profile()
}

if (identical(profile, "quick")) {
  source_profile_helpers("standard")
  run_tests("00-|02-")
  finish_profile()
}

if (identical(profile, "filter")) {
  if (length(extra_args) == 0L) {
    stop("Profile 'filter' requires at least one testthat filter.", call. = FALSE)
  }
  source_profile_helpers("standard")
  run_tests(paste(extra_args, collapse = "|"))
  finish_profile()
}

stop(
  "Unknown profile: ", profile,
  ". Use standard, refresh-standard, certification, quick, filter, or release.",
  call. = FALSE
)
