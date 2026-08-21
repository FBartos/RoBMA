.robma_test_project_root <- local({

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    ofile <- frame[["ofile"]]
    if (is.null(ofile)) "" else ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  candidates <- c(
    dirname(normalizePath(ofiles, winslash = "/", mustWork = FALSE)),
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  )

  for (candidate in unique(candidates)) {
    repeat {
      description <- file.path(candidate, "DESCRIPTION")
      if (file.exists(description)) {
        package <- tryCatch(
          read.dcf(description, fields = "Package")[[1L]],
          error = function(error) ""
        )
        if (identical(package, "RoBMA")) {
          return(candidate)
        }
      }

      parent <- dirname(candidate)
      if (identical(parent, candidate)) {
        break
      }
      candidate <- parent
    }
  }

  stop(
    "Cannot locate the RoBMA package root from the sourced file or working directory.",
    call. = FALSE
  )
})

if (!exists("quiet_llm_reporter", mode = "function")) {
  source(
    file.path(
      .robma_test_project_root,
      "tests",
      "testthat",
      "common-functions.R"
    ),
    local = TRUE
  )
}
if (!exists("review_test_snapshots", mode = "function") ||
    !exists("test_scenarios", mode = "function")) {
  source(
    file.path(
      .robma_test_project_root,
      "tests",
      "scenarios",
      "helper-scenarios.R"
    ),
    local = TRUE
  )
}


.robma_test_validate_flag <- function(value, name) {

  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }

  return(value)
}


.robma_test_rscript <- function() {

  command <- Sys.which("Rscript")
  if (!nzchar(command)) {
    command <- file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
  }

  return(command)
}


.run_robma_test_profile <- function(profile, clean = FALSE, filter = NULL) {

  profile_script <- file.path(
    .robma_test_project_root,
    "tools",
    "test-profile.R"
  )
  arguments <- c(shQuote(profile_script), profile)
  if (clean) {
    arguments <- c(arguments, "--clean")
  }
  if (!is.null(filter)) {
    arguments <- c(arguments, shQuote(filter))
  }

  message("Running RoBMA test profile: ", profile)
  status <- system2(.robma_test_rscript(), args = arguments)
  if (status != 0L) {
    stop("RoBMA test profile failed: ", profile, call. = FALSE)
  }

  return(invisible(TRUE))
}


# Run ordinary and certification tests from an interactive development session.
# A filter selects standard testthat files; an unfiltered run also executes all
# independently bounded certification cases.
test_tests <- function(filter = NULL, reporter = "llm", refit = FALSE,
                       update = FALSE, update_timings = FALSE,
                       regenerate = FALSE,
                       load_package = TRUE, stop_on_failure = FALSE,
                       root = file.path(
                         .robma_test_project_root, "tests", "testthat"
                       )) {

  refit           <- .robma_test_validate_flag(refit, "refit")
  update          <- .robma_test_validate_flag(update, "update")
  update_timings  <- .robma_test_validate_flag(
    update_timings,
    "update_timings"
  )
  regenerate      <- .robma_test_validate_flag(regenerate, "regenerate")
  load_package    <- .robma_test_validate_flag(load_package, "load_package")
  stop_on_failure <- .robma_test_validate_flag(
    stop_on_failure,
    "stop_on_failure"
  )
  if (!is.null(filter) &&
      (!is.character(filter) || length(filter) != 1L || is.na(filter))) {
    stop("'filter' must be NULL or one regular expression.", call. = FALSE)
  }
  if (!is.character(reporter) || length(reporter) != 1L ||
      is.na(reporter) || !nzchar(reporter)) {
    stop("'reporter' must be one reporter name.", call. = FALSE)
  }
  if (update_timings) {
    stop(
      "Ordinary tests have no timing baselines; the standard profile ",
      "enforces its 15-minute total budget.",
      call. = FALSE
    )
  }

  refit  <- refit || regenerate
  update <- update || regenerate
  root         <- normalizePath(root, winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(
    file.path(root, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  if (!identical(project_root, .robma_test_project_root)) {
    stop("'root' must be RoBMA's 'tests/testthat' directory.", call. = FALSE)
  }

  if (load_package) {
    devtools::load_all(project_root, quiet = TRUE)
  }

  environment_names <- c(
    "AGENT",
    "NOT_CRAN",
    "ROBMA_TEST_REPORTER",
    "ROBMA_TEST_STOP_ON_FAILURE",
    "ROBMA_TEST_QUIET_SKIPS",
    "ROBMA_TEST_FULL_VISUALS",
    "ROBMA_TEST_FULL_DIAGNOSTICS",
    "ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS",
    "ROBMA_SCENARIO_REGENERATE",
    "ROBMA_SCENARIO_REFIT",
    "ROBMA_SCENARIO_UPDATE",
    "ROBMA_SCENARIO_UPDATE_TIMINGS",
    "ROBMA_SCENARIO_RUNNER"
  )
  old_environment <- Sys.getenv(environment_names, unset = NA_character_)
  on.exit({
    for (i in seq_along(environment_names)) {
      name  <- environment_names[[i]]
      value <- old_environment[[i]]
      if (is.na(value)) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(value), name))
      }
    }
  }, add = TRUE)
  on.exit({
    if (interactive()) {
      try(review_test_snapshots(root = root))
    }
  }, add = TRUE)

  Sys.unsetenv("AGENT")
  Sys.setenv(
    NOT_CRAN                          = "true",
    ROBMA_TEST_REPORTER               = reporter,
    ROBMA_TEST_STOP_ON_FAILURE        = if (stop_on_failure) "TRUE" else "FALSE",
    ROBMA_TEST_QUIET_SKIPS            = "TRUE",
    ROBMA_TEST_FULL_VISUALS           = "TRUE",
    ROBMA_TEST_FULL_DIAGNOSTICS       = "TRUE",
    ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS = if (update) "TRUE" else "FALSE",
    ROBMA_SCENARIO_REGENERATE         = "FALSE",
    ROBMA_SCENARIO_REFIT              = "FALSE",
    ROBMA_SCENARIO_UPDATE             = "FALSE",
    ROBMA_SCENARIO_UPDATE_TIMINGS     = "FALSE",
    ROBMA_SCENARIO_RUNNER             = "TRUE"
  )

  results <- list()
  if (is.null(filter)) {
    results[["refresh_standard"]] <- .run_robma_test_profile(
      "refresh-standard",
      clean = refit
    )
    results[["standard"]] <- .run_robma_test_profile("standard")
    results[["certification"]] <- .run_robma_test_profile(
      "certification",
      clean = refit
    )
  } else {
    if (refit) {
      results[["refresh_standard"]] <- .run_robma_test_profile(
        "refresh-standard",
        clean = TRUE
      )
    }
    results[["tests"]] <- .run_robma_test_profile(
      "filter",
      filter = filter
    )
  }

  return(invisible(results))
}
