context("Maintainer scenario helpers")

source(testthat::test_path("..", "scenarios", "helper-scenarios.R"))


.scenario_test_root <- function() {

  root <- tempfile("robma-scenario-test-")
  dir.create(root, recursive = TRUE)
  return(root)
}


test_that("scenario_fit caches fits until regeneration is requested", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  counter <- 0L

  scenario_start("unit", root = root, create_missing = TRUE)
  first <- scenario_fit("model", {
    counter <- counter + 1L
    list(value = counter)
  })
  second <- scenario_fit("model", {
    counter <- counter + 1L
    list(value = counter)
  })

  expect_identical(counter, 1L)
  expect_identical(first, second)
  expect_true(file.exists(file.path(root, "cache", "unit", "model.rds")))

  REGENERATE_SCENARIO_FILES <- TRUE
  scenario_start("unit", root = root, create_missing = TRUE)
  third <- scenario_fit("model", {
    counter <- counter + 1L
    list(value = counter)
  })

  expect_identical(counter, 2L)
  expect_identical(third[["value"]], 2L)
  expect_identical(readRDS(file.path(
    root, "cache", "unit", "model.rds"
  )), third)
})


test_that("scenario_text adds, compares, and regenerates tracked output", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root, create_missing = TRUE, width = 80L)
  value <- scenario_text("summary", {
    print(data.frame(value = 1:2))
    invisible(42L)
  })
  path <- file.path(root, "results", "unit", "summary.txt")

  expect_identical(value, 42L)
  expect_true(file.exists(path))
  expect_equal(readLines(path, warn = FALSE), c("  value", "1     1", "2     2"))

  scenario_text("summary", {
    print(data.frame(value = 1:2))
    invisible(42L)
  })
  expect_failure(
    scenario_text("summary", print("changed")),
    "Scenario text.*changed"
  )

  scenario_start(
    "unit",
    root           = root,
    regenerate     = TRUE,
    create_missing = TRUE
  )
  scenario_text("summary", print("changed"))
  expect_equal(readLines(path, warn = FALSE), "[1] \"changed\"")
})


test_that("scenario_text rejects missing baselines when creation is disabled", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root, create_missing = FALSE)
  expect_failure(
    scenario_text("missing", print("output")),
    "Missing locked scenario text"
  )
  expect_false(file.exists(file.path(
    root, "results", "unit", "missing.txt"
  )))
})


test_that("scenario_plot delegates plot comparison to vdiffr", {

  skip_if_not_installed("vdiffr")
  skip_if_no_vdiffr_snapshots()
  scenario_start(
    "00-scenario-helpers",
    root           = testthat::test_path(),
    create_missing = TRUE
  )

  scenario_plot("simple-base-plot", {
    graphics::plot(1:3, 1:3, xlab = "x", ylab = "y")
  })
})
