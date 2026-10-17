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
