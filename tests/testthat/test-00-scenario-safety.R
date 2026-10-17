context("Scenario artifact and run safety")

source(testthat::test_path("..", "scenarios", "helper-scenarios.R"))


test_that("scenario and text names reserve the candidate suffix before evaluation", {

  root <- tempfile("robma-scenario-safety-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  bytes <- function(path) {
    if (!file.exists(path)) return(NULL)
    readBin(path, "raw", n = file.size(path))
  }
  scenario_start("unit", root = root, create_missing = TRUE,
                 update = FALSE, show_output = FALSE)
  scenario_text("summary", print("LOCKED PRIMARY MARKER"))
  primary <- file.path(root, "results", "unit", "summary.txt")
  original <- bytes(primary)
  evaluated <- FALSE
  for (name in c("summary.new", "summary.NEW", "summary.NeW")) {
    expect_error(scenario_text(name, {
      evaluated <- TRUE
      print("LOCKED OTHER DISTINCT MARKER")
    }), "reserved")
  }
  expect_false(evaluated)
  expect_identical(bytes(primary), original)
  expect_false(file.exists(file.path(root, "results", "unit", "summary.new.txt")))

  # Both the matching branch and a new primary would otherwise unlink the
  # candidate path. A reserved logical name must reach neither branch.
  for (create in c(FALSE, TRUE)) {
    scenario_start("unit", root = root, create_missing = create,
                   update = FALSE, show_output = FALSE)
    locked <- file.path(root, "results", "unit", "blocked.new.txt")
    candidate <- .scenario_candidate_path(locked)
    .scenario_write_lines('[1] "old candidate"', candidate)
    if (!create) .scenario_write_lines('[1] "same"', locked)
    before <- bytes(candidate)
    expect_error(scenario_text("blocked.new", print("same")), "reserved")
    expect_identical(bytes(candidate), before)
    if (create) expect_false(file.exists(locked))
    if (!create) unlink(locked)
  }

  timing_candidate <- file.path(root, "timings", "unit.new.tsv")
  .scenario_write_lines("CURRENT LAST RUN", timing_candidate)
  before <- bytes(timing_candidate)
  for (name in c("unit.new", "unit.NEW")) {
    expect_error(scenario_start(name, root = root), "reserved")
  }
  expect_identical(bytes(timing_candidate), before)
  expect_identical(.scenario_validate_name("fit.new", "fit"), "fit.new")
  expect_identical(.scenario_validate_name("summary.new.details", "text"),
                   "summary.new.details")
})


test_that("ordinary text artifacts preserve locked bytes and candidate behavior", {

  root <- tempfile("robma-scenario-text-control-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  scenario_start("unit", root = root, create_missing = TRUE,
                 show_output = FALSE)
  scenario_text("summary", print("LOCKED PRIMARY MARKER"))
  scenario_text("details", print("LOCKED OTHER DISTINCT MARKER"))
  paths <- file.path(root, "results", "unit", c("summary.txt", "details.txt"))
  before <- tools::md5sum(paths)
  expect_failure(scenario_text("summary", print("changed")), "candidate is cached")
  expect_identical(tools::md5sum(paths), before)
  expect_identical(readLines(.scenario_candidate_path(paths[[1L]])),
                   '[1] "changed"')
  scenario_text("summary", print("LOCKED PRIMARY MARKER"))
  expect_false(file.exists(.scenario_candidate_path(paths[[1L]])))
  expect_identical(tools::md5sum(paths), before)
})


test_that("plot aliases fail before evaluating or writing a second artifact", {

  skip_if_not_installed("vdiffr")
  root <- tempfile("robma-scenario-plot-alias-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env <- environment(scenario_plot)
  old <- mget(c(".scenario_is_interactive", ".scenario_snapshot_context"),
              envir = helper_env, inherits = FALSE)
  on.exit(list2env(old, helper_env), add = TRUE)
  assign(".scenario_is_interactive", function() TRUE, helper_env)
  assign(".scenario_snapshot_context", function() NULL, helper_env)
  scenario_start("unit", root = root, create_missing = TRUE,
                 show_output = FALSE)
  scenario_plot("A_b", graphics::plot(1:3))
  path <- file.path(root, "_snaps", "unit", .scenario_plot_filename("A_b"))
  before <- tools::md5sum(path)
  evaluated <- FALSE
  for (name in c("a-b", "a.b", "A--b")) {
    expect_error(scenario_plot(name, {
      evaluated <- TRUE
      graphics::plot(1:4)
    }), "same.*filename|filename.*already")
  }
  expect_false(evaluated)
  expect_identical(tools::md5sum(path), before)
  expect_false(file.exists(.scenario_candidate_path(path)))
  scenario_plot("A_b", graphics::plot(1:3))
  expect_identical(tools::md5sum(path), before)
  scenario_start("unit", root = root, create_missing = TRUE,
                 update = TRUE, show_output = FALSE)
  scenario_plot("a-b", graphics::plot(1:3))
  expect_true(file.exists(path))
})


.scenario_safety_timing_run <- function(outcome, strict = FALSE,
                                        update_timings = TRUE) {

  project <- tempfile("robma-scenario-run-safety-")
  root <- file.path(project, "tests", "scenarios")
  dir.create(root, recursive = TRUE)
  dir.create(file.path(project, "tests", "testthat"))
  file.copy(file.path(.scenario_helpers_dir, "helper-scenarios.R"), root)
  file.copy(.scenario_visual_writer, file.path(project, "tests", "testthat"))
  scenario_start("unit", root = root, show_output = FALSE)
  provenance <- .scenario_timing_provenance()
  baseline <- data.frame(
    type = c("text", "time"), name = c("first", "obsolete"),
    elapsed = c(8, 8), memory_gb = NA_real_,
    r_version = provenance[["r_version"]], platform = provenance[["platform"]],
    stringsAsFactors = FALSE
  )
  expected <- file.path(root, "timings", "unit.tsv")
  .scenario_write_timings(baseline, expected)
  before <- readBin(expected, "raw", n = file.size(expected))
  action <- switch(outcome,
    failure = 'testthat::expect_true(FALSE, info = "Deliberate synthetic scenario failure")',
    skip = 'testthat::skip("Deliberate synthetic rendering skip")',
    error = 'stop("Deliberate synthetic scenario error")',
    interrupt = 'stop(structure(list(message = "Deliberate synthetic interruption"), class = c("interrupt", "condition")))',
    success = 'testthat::succeed()'
  )
  .scenario_write_lines(c(
    paste0('scenario_start("unit", root = ', deparse(root),
           ', show_output = FALSE)'),
    'testthat::test_that("synthetic scenario run status", {',
    paste0('  .scenario_register_timing("text", "first", ',
           if (update_timings) "9" else "7", ')'),
    paste0("  ", action),
    '})'
  ), file.path(root, "test-unit.R"))
  result <- NULL
  error <- NULL
  warnings <- character()
  console <- capture.output(withCallingHandlers(
    result <- tryCatch(
      test_scenario(root = root, load_package = FALSE, reporter = "silent",
                    update_timings = update_timings, stop_on_failure = strict),
      error = function(e) { error <<- e; NULL },
      interrupt = function(e) { error <<- e; NULL }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  return(list(
    project = project, before = before,
    after = readBin(expected, "raw", n = file.size(expected)),
    baseline = .scenario_read_timings(expected),
    last = .scenario_read_timings(.scenario_candidate_path(expected)),
    result = if (is.null(result)) NULL else as.data.frame(result),
    error = error, warnings = warnings, console = console
  ))
}


test_that("failed and skipped managed runs merge without losing unmeasured timings", {

  for (outcome in c("failure", "skip", "error", "interrupt")) {
    run <- .scenario_safety_timing_run(outcome)
    on.exit(unlink(run[["project"]], recursive = TRUE), add = TRUE)
    expect_setequal(run[["baseline"]][["name"]], c("first", "obsolete"))
    expect_equal(run[["baseline"]][["elapsed"]][
      run[["baseline"]][["name"]] == "first"], 9, info = outcome)
    expect_equal(run[["baseline"]][["elapsed"]][
      run[["baseline"]][["name"]] == "obsolete"], 8, info = outcome)
    expect_identical(run[["last"]][["name"]], "first", info = outcome)
    expect_equal(run[["last"]][["elapsed"]], 9, info = outcome)
    if (outcome == "failure") expect_equal(run[["result"]][["failed"]], 1L)
    if (outcome == "skip") expect_true(run[["result"]][["skipped"]])
    if (outcome == "error") expect_true(run[["result"]][["error"]])
    if (outcome == "interrupt") expect_s3_class(run[["error"]], "interrupt")
  }
})


test_that("strict failure preserves baseline and complete success may replace it", {

  strict <- .scenario_safety_timing_run("failure", strict = TRUE)
  on.exit(unlink(strict[["project"]], recursive = TRUE), add = TRUE)
  expect_s3_class(strict[["error"]], "error")
  expect_identical(strict[["after"]], strict[["before"]])
  expect_identical(strict[["last"]][["name"]], "first")
  expect_equal(strict[["last"]][["elapsed"]], 9)

  complete <- .scenario_safety_timing_run("success")
  on.exit(unlink(complete[["project"]], recursive = TRUE), add = TRUE)
  expect_identical(complete[["baseline"]][["name"]], "first")
  expect_equal(complete[["baseline"]][["elapsed"]], 9)
  expect_equal(complete[["result"]][["failed"]], 0L)
  expect_false(complete[["result"]][["skipped"]])

  partial <- .scenario_safety_timing_run("success", update_timings = FALSE)
  on.exit(unlink(partial[["project"]], recursive = TRUE), add = TRUE)
  expect_setequal(partial[["baseline"]][["name"]], c("first", "obsolete"))
  expect_equal(partial[["baseline"]][["elapsed"]][
    partial[["baseline"]][["name"]] == "first"], 7)
  expect_equal(partial[["baseline"]][["elapsed"]][
    partial[["baseline"]][["name"]] == "obsolete"], 8)
})


test_that("direct timing warnings compare with the immutable run-start baseline", {

  helper_env <- environment(scenario_start)
  old_interactive <- get(".scenario_is_interactive", helper_env)
  on.exit(assign(".scenario_is_interactive", old_interactive, helper_env),
          add = TRUE)
  assign(".scenario_is_interactive", function() TRUE, helper_env)
  withr::local_envvar(ROBMA_SCENARIO_RUNNER = NA_character_)
  for (case in c("accepted", "best", "unchanged")) {
    root <- tempfile("robma-scenario-baseline-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    scenario_start("unit", root = root, show_output = FALSE)
    .scenario_register_timing("text", "a", 10)
    .scenario_register_timing("text", "b", 10)
    scenario_start("unit", root = root, show_output = FALSE,
                   update_timings = case == "accepted")
    current <- switch(case, accepted = c(12, 12), best = c(8, 14),
                      unchanged = c(10, 10))
    warnings <- character()
    withCallingHandlers({
      .scenario_register_timing("text", "a", current[[1L]])
      .scenario_register_timing("text", "b", current[[2L]])
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    issues <- .scenario_timing_issues(final = TRUE)
    average <- issues[["text"]][issues[["id"]] == "average"]
    if (case == "accepted") {
      expect_length(average, 1L)
      expect_match(average, "average timing regression: 20%", fixed = TRUE)
      expect_true(any(grepl("average timing regression: 20%", warnings,
                           fixed = TRUE)))
    } else {
      expect_length(average, 0L)
      expect_false(any(grepl("average timing regression", warnings, fixed = TRUE)))
    }
    baseline <- .scenario_read_timings(file.path(root, "timings/unit.tsv"))
    expect_equal(baseline[["elapsed"]],
                 switch(case, accepted = c(12, 12), best = c(8, 10),
                        unchanged = c(10, 10)))
    last <- .scenario_read_timings(file.path(root, "timings/unit.new.tsv"))
    expect_equal(last[["elapsed"]], current)
  }
})


test_that("default scenario loading prepares the release DLL without test helpers", {

  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  skip_if_not(file.exists(file.path(package_root, "tools/optimized-dll.R")),
              "Developer release loader requires a source checkout.")
  project <- tempfile("robma-interactive-loader-")
  root <- file.path(project, "tests", "scenarios")
  dir.create(root, recursive = TRUE)
  dir.create(file.path(project, "tests", "testthat"))
  dir.create(file.path(project, "R"))
  dir.create(file.path(project, "tools"))
  on.exit(unlink(project, recursive = TRUE), add = TRUE)
  for (file in c("tools/optimized-dll.R", "tests/scenarios/helper-scenarios.R",
                  "tests/testthat/helper-visual-writer.R")) {
    file.copy(file.path(package_root, file), file.path(project, file))
  }
  writeLines(c("Package: RoBMA", "Version: 0.0.0.9000",
               "Title: Synthetic Scenario Loader Fixture",
               "Description: A minimal scenario loader fixture without native code.",
               "Author: RoBMA tests", "Maintainer: RoBMA tests <tests@example.org>",
               "License: MIT"), file.path(project, "DESCRIPTION"))
  writeLines("export(fixture_marker)", file.path(project, "NAMESPACE"))
  writeLines('fixture_marker <- function() TRUE', file.path(project, "R", "marker.R"))
  writeLines('stop("PACKAGE TEST HELPERS MUST NOT LOAD")',
             file.path(project, "tests/testthat/helper-must-not-load.R"))
  writeLines(c(paste0('scenario_start("unit", root = ', deparse(root), ')'),
               'testthat::test_that("synthetic loader smoke", testthat::succeed())'),
             file.path(root, "test-unit.R"))
  script <- file.path(project, "loader-smoke.R")
  writeLines(c(paste0('source(', deparse(file.path(root, "helper-scenarios.R")), ')'),
               paste0('result <- test_scenario(root = ', deparse(root), ', reporter = "silent")'),
               'stopifnot(all(as.data.frame(result)$failed == 0L), !any(as.data.frame(result)$error))',
               'cat("LOADED_NAMESPACE", getNamespaceInfo(asNamespace("RoBMA"), "path"), "\\n")'), script)
  executable <- file.path(R.home("bin"),
                          if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  withr::local_envvar(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  for (debug in c("FALSE", "TRUE")) {
    withr::local_envvar(ROBMA_DEBUG_DLL = debug)
    run <- processx::run(executable, c("--vanilla", script),
                         error_on_status = FALSE, timeout = 60000)
    cat("\nSynthetic scenario loader exit ", run$status, "\n", run$stdout, run$stderr, sep = "")
    expect_equal(run$status, 0L)
    expect_match(run$stdout, normalizePath(project, winslash = "/"), fixed = TRUE)
    expect_match(run$stderr, if (debug == "TRUE") "debug opt-out" else "release", fixed = TRUE)
    expect_false(grepl("PACKAGE TEST HELPERS MUST NOT LOAD", run$stderr, fixed = TRUE))
  }
})
