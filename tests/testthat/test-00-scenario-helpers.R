context("Maintainer scenario helpers")

source(testthat::test_path("..", "scenarios", "helper-scenarios.R"))


.scenario_test_root <- function() {

  root <- tempfile("robma-scenario-test-")
  dir.create(root, recursive = TRUE)
  return(root)
}


test_that("scenario_fit invalidates changed calls and supports regeneration", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  counter <- 0L

  scenario_start("unit", root = root, create_missing = TRUE)
  first <- scenario_fit("model.one", {
    counter <- counter + 1L
    list(value = counter)
  })
  second <- scenario_fit("model.one", {
    counter <- counter + 1L
    list(value = counter)
  })

  expect_identical(counter, 1L)
  expect_identical(first, second)
  expect_true(file.exists(file.path(root, "cache", "unit", "model.one.rds")))
  expect_true(file.exists(file.path(
    root, "cache", "unit", "model.one.call.txt"
  )))

  third <- expect_message(
    scenario_fit("model.one", {
      counter <- counter + 1L
      list(value = counter, changed = TRUE)
    }),
    "Fit call changed.*refitting"
  )
  fourth <- scenario_fit("model.one", {
    counter <- counter + 1L
    list(value = counter, changed = TRUE)
  })

  expect_identical(counter, 2L)
  expect_true(third[["changed"]])
  expect_identical(third, fourth)

  REGENERATE_SCENARIO_FILES <- TRUE
  scenario_start("unit", root = root, create_missing = TRUE)
  fifth <- scenario_fit("model.one", {
    counter <- counter + 1L
    list(value = counter, changed = TRUE)
  })

  expect_identical(counter, 3L)
  expect_identical(fifth[["value"]], 3L)
  expect_identical(readRDS(file.path(
    root, "cache", "unit", "model.one.rds"
  )), fifth)
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


test_that("scenario_text replays captured output interactively", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env      <- environment(scenario_text)
  old_interactive <- get(".scenario_is_interactive", envir = helper_env)
  on.exit(assign(
    ".scenario_is_interactive",
    old_interactive,
    envir = helper_env
  ), add = TRUE)
  assign(
    ".scenario_is_interactive",
    function() TRUE,
    envir = helper_env
  )

  scenario_start("unit", root = root, create_missing = TRUE)
  expect_output(
    scenario_text("visible", print("interactive output")),
    '[1] "interactive output"',
    fixed = TRUE
  )
})


test_that("scenario plot evaluation restores base graphics parameters", {

  state <- new.env(parent = emptyenv())
  .with_temp_plot_device({
    state[["before"]] <- graphics::par("mar")
    .scenario_evaluate_plot(quote({
      graphics::par(mar = c(2, 1, 2, 1))
    }), environment())
    state[["after_success"]] <- graphics::par("mar")

    expect_error(
      .scenario_evaluate_plot(quote({
        graphics::par(mar = c(1, 1, 1, 1))
        stop("plot failed")
      }), environment()),
      "plot failed"
    )
    state[["after_error"]] <- graphics::par("mar")
  })

  expect_equal(state[["after_success"]], state[["before"]])
  expect_equal(state[["after_error"]], state[["before"]])
})


test_that("scenario_plot accepts mixed-case names and draws once", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env   <- environment(scenario_plot)
  helper_names <- c(
    ".scenario_is_interactive",
    ".scenario_snapshot_context"
  )
  old_helpers <- mget(helper_names, envir = helper_env, inherits = FALSE)
  on.exit(list2env(old_helpers, envir = helper_env), add = TRUE)
  assign(".scenario_is_interactive", function() TRUE, envir = helper_env)
  assign(".scenario_snapshot_context", function() NULL, envir = helper_env)
  state <- new.env(parent = emptyenv())
  state[["draws"]] <- 0L

  scenario_start("unit", root = root, create_missing = TRUE)
  .with_temp_plot_device({
    graphics::plot(1:3, 1:3)
    graphics::par(new = TRUE)
    scenario_plot("mu_BF_comparison", {
      state[["entry_new"]] <- graphics::par("new")
      state[["draws"]] <- state[["draws"]] + 1L
      graphics::plot(1:3, 1:3)
      graphics::par(new = TRUE)
    })
    state[["exit_new"]] <- graphics::par("new")
  })

  expect_identical(state[["draws"]], 1L)
  expect_false(state[["entry_new"]])
  expect_false(state[["exit_new"]])
})


test_that("scenario_plot delegates plot comparison to vdiffr", {

  skip_if_not_installed("vdiffr")
  skip_if_no_vdiffr_snapshots()
  scenario_start(
    "00-scenario-helpers",
    root           = testthat::test_path(),
    create_missing = TRUE
  )
  helper_env      <- environment(scenario_plot)
  old_interactive <- get(".scenario_is_interactive", envir = helper_env)
  on.exit(assign(
    ".scenario_is_interactive",
    old_interactive,
    envir = helper_env
  ), add = TRUE)
  assign(
    ".scenario_is_interactive",
    function() TRUE,
    envir = helper_env
  )
  draws <- new.env(parent = emptyenv())
  draws[["count"]] <- 0L

  .with_temp_plot_device({
    scenario_plot("simple-base-plot", {
      draws[["count"]] <- draws[["count"]] + 1L
      graphics::plot(1:3, 1:3, xlab = "x", ylab = "y")
    })
  })
  expect_gte(draws[["count"]], 2L)
})
