context("Maintainer scenario helpers")

source(testthat::test_path("..", "scenarios", "helper-scenarios.R"))

.scenario_test_environment_names <- c(
  "ROBMA_SCENARIO_REGENERATE",
  "ROBMA_SCENARIO_REFIT",
  "ROBMA_SCENARIO_UPDATE",
  "ROBMA_SCENARIO_UPDATE_TIMINGS",
  "ROBMA_SCENARIO_RUNNER"
)
.scenario_test_old_environment <- Sys.getenv(
  .scenario_test_environment_names,
  unset = NA_character_
)
teardown({
  for (i in seq_along(.scenario_test_environment_names)) {
    name  <- .scenario_test_environment_names[[i]]
    value <- .scenario_test_old_environment[[i]]
    if (is.na(value)) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, stats::setNames(list(value), name))
    }
  }
})
Sys.setenv(
  ROBMA_SCENARIO_REGENERATE     = "FALSE",
  ROBMA_SCENARIO_REFIT          = "FALSE",
  ROBMA_SCENARIO_UPDATE         = "FALSE",
  ROBMA_SCENARIO_UPDATE_TIMINGS = "FALSE",
  ROBMA_SCENARIO_RUNNER         = "TRUE"
)


.scenario_test_root <- function() {

  root <- tempfile("robma-scenario-test-")
  dir.create(root, recursive = TRUE)
  return(root)
}


test_that("plural scenario runner retains the established interface", {

  expect_identical(test_scenarios, test_scenario)
  expect_identical(formals(test_scenarios), formals(test_scenario))
})


test_that("nested helper sourcing reuses the runner state", {

  helper_env <- environment(scenario_start)
  nested_env <- new.env(parent = helper_env)
  source(
    normalizePath(
      testthat::test_path("..", "scenarios", "helper-scenarios.R"),
      winslash = "/",
      mustWork = TRUE
    ),
    local = nested_env
  )

  expect_false(exists(
    ".scenario_state", envir = nested_env, inherits = FALSE
  ))
  expect_identical(
    get(".scenario_state", envir = nested_env, inherits = TRUE),
    .scenario_state
  )
})


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

  scenario_start(
    "unit",
    root           = root,
    regenerate     = TRUE,
    create_missing = TRUE
  )
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


test_that("scenario defaults distinguish direct and managed interactive runs", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env      <- environment(scenario_start)
  old_interactive <- get(
    ".scenario_is_interactive",
    envir    = helper_env,
    inherits = FALSE
  )
  old_snapshot_context <- get(
    ".scenario_snapshot_context",
    envir    = helper_env,
    inherits = FALSE
  )
  on.exit(assign(
    ".scenario_is_interactive",
    old_interactive,
    envir = helper_env
  ), add = TRUE)
  on.exit(assign(
    ".scenario_snapshot_context",
    old_snapshot_context,
    envir = helper_env
  ), add = TRUE)
  assign(".scenario_is_interactive", function() TRUE, envir = helper_env)
  assign(".scenario_snapshot_context", function() NULL, envir = helper_env)

  environment_names <- c(
    "ROBMA_SCENARIO_REGENERATE",
    "ROBMA_SCENARIO_REFIT",
    "ROBMA_SCENARIO_UPDATE",
    "ROBMA_SCENARIO_UPDATE_TIMINGS",
    "ROBMA_SCENARIO_RUNNER"
  )
  old_environment <- Sys.getenv(environment_names, unset = NA_character_)
  on.exit({
    for (i in seq_along(environment_names)) {
      name <- environment_names[[i]]
      if (is.na(old_environment[[i]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old_environment[[i]]), name))
      }
    }
  }, add = TRUE)
  Sys.unsetenv(environment_names)

  scenario_start(
    "unit",
    root        = root,
    update      = TRUE,
    show_output = FALSE
  )
  scenario_text("result", "old")
  .with_temp_plot_device({
    scenario_plot("direct-plot", graphics::plot(1:3, 1:3))
  })

  direct <- scenario_start("unit", root = root)
  expect_false(direct[["refit"]])
  expect_true(direct[["show_output"]])
  expect_false(direct[["update"]])
  expect_false(direct[["update_timings"]])
  expect_true(direct[["create_missing"]])
  expect_output(
    scenario_text("new-result", "new baseline"),
    '[1] "new baseline"',
    fixed = TRUE
  )
  expect_true(file.exists(file.path(
    root, "results", "unit", "new-result.txt"
  )))
  direct_output <- capture.output({
    expect_failure(scenario_text("result", "new"), "candidate is cached")
  })
  expect_match(paste(direct_output, collapse = "\n"), '[1] "new"',
               fixed = TRUE)
  path <- file.path(root, "results", "unit", "result.txt")
  expect_equal(readLines(path, warn = FALSE), '[1] "old"')
  expect_true(file.exists(.scenario_candidate_path(path)))
  .with_temp_plot_device({
    expect_failure(
      scenario_plot("direct-plot", graphics::plot(1:4, 1:4)),
      "Scenario plot.*changed"
    )
  })
  plot_path <- file.path(root, "_snaps", "unit", "direct-plot.svg")
  expect_true(file.exists(plot_path))
  expect_true(file.exists(.scenario_candidate_path(plot_path)))

  Sys.setenv(ROBMA_SCENARIO_RUNNER = "TRUE")
  managed <- scenario_start("unit", root = root)
  expect_false(managed[["refit"]])
  expect_false(managed[["show_output"]])
  expect_false(managed[["update"]])
  expect_false(managed[["update_timings"]])
  expect_false(managed[["create_missing"]])
  expect_failure(scenario_text("result", "managed"), "candidate is cached")
  expect_equal(readLines(path, warn = FALSE), '[1] "old"')
})


test_that("interactive timing backfills but requires consent for slowdowns", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env      <- environment(scenario_start)
  old_interactive <- get(
    ".scenario_is_interactive", envir = helper_env, inherits = FALSE
  )
  on.exit(assign(
    ".scenario_is_interactive",
    old_interactive,
    envir = helper_env
  ), add = TRUE)
  old_runner <- Sys.getenv("ROBMA_SCENARIO_RUNNER", unset = NA_character_)
  on.exit({
    if (is.na(old_runner)) {
      Sys.unsetenv("ROBMA_SCENARIO_RUNNER")
    } else {
      Sys.setenv(ROBMA_SCENARIO_RUNNER = old_runner)
    }
  }, add = TRUE)
  Sys.unsetenv("ROBMA_SCENARIO_RUNNER")
  assign(".scenario_is_interactive", function() TRUE, envir = helper_env)

  scenario_start("unit", root = root)
  .scenario_register_timing("text", "summary", 10)
  timing_path    <- file.path(root, "timings", "unit.tsv")
  candidate_path <- file.path(root, "timings", "unit.new.tsv")
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 10)

  expect_warning(
    .scenario_register_timing("text", "summary", 12),
    "average timing regression: 20%"
  )
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 10)
  expect_equal(.scenario_read_timings(candidate_path)[["elapsed"]], 12)

  scenario_start("unit", root = root, update_timings = TRUE)
  expect_warning(
    .scenario_register_timing("text", "summary", 12),
    "average timing regression: 20%"
  )
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 12)
  expect_false(file.exists(candidate_path))
})


test_that("scenario_time returns values and records successful computations", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env <- environment(scenario_time)
  old_clock  <- get(".scenario_clock", envir = helper_env, inherits = FALSE)
  on.exit(assign(".scenario_clock", old_clock, envir = helper_env), add = TRUE)
  ticks <- c(10, 13, 20, 24, 30)
  assign(
    ".scenario_clock",
    function() {

      value <- ticks[[1L]]
      ticks <<- ticks[-1L]
      value
    },
    envir = helper_env
  )
  counter <- 0L

  scenario_start("unit", root = root, update_timings = TRUE)
  visible <- withVisible(scenario_time("visible", {
    counter <- counter + 1L
    42L
  }))
  hidden <- withVisible(scenario_time("hidden", invisible("value")))
  expect_error(scenario_time("failed", stop("calculation failed")), "calculation failed")

  expect_identical(counter, 1L)
  expect_identical(visible, list(value = 42L, visible = TRUE))
  expect_identical(hidden, list(value = "value", visible = FALSE))
  current <- .scenario_current_timings()
  expect_equal(current[c("type", "name", "elapsed")], data.frame(
    type    = c("time", "time"),
    name    = c("hidden", "visible"),
    elapsed = c(4, 3)
  ))
  expect_false("failed" %in% current[["name"]])
  expect_length(list.files(file.path(root, "results", "unit")), 0L)
  expect_length(list.files(file.path(root, "_snaps", "unit")), 0L)
  expect_no_warning(.scenario_finalize_timing())
  expect_equal(
    .scenario_read_timings(file.path(root, "timings", "unit.tsv"))[
      c("type", "name", "elapsed")
    ],
    current[c("type", "name", "elapsed")]
  )
})


test_that("scenario_fit supports optional targeted cache versions", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  counter <- 0L

  scenario_start("unit", root = root, create_missing = TRUE)
  first <- scenario_fit("versioned", {
    counter <- counter + 1L
    counter
  })
  second <- scenario_fit("versioned", {
    counter <- counter + 1L
    counter
  })
  third <- expect_message(
    scenario_fit("versioned", {
      counter <- counter + 1L
      counter
    }, cache_version = 1L),
    "Fit call changed.*refitting"
  )
  fourth <- scenario_fit("versioned", {
    counter <- counter + 1L
    counter
  }, cache_version = 1L)

  expect_identical(first, 1L)
  expect_identical(second, 1L)
  expect_identical(third, 2L)
  expect_identical(fourth, 2L)
  expect_error(
    scenario_fit("invalid-version", 1L, cache_version = 0L),
    "positive integer"
  )
})


test_that("scenario refitting and output updates are independent", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  counter <- 0L

  scenario_start("unit", root = root, create_missing = TRUE)
  first <- scenario_fit("model", {
    counter <- counter + 1L
    counter
  })
  scenario_text("result", "old")
  path <- file.path(root, "results", "unit", "result.txt")

  scenario_start(
    "unit",
    root           = root,
    refit          = TRUE,
    update         = FALSE,
    create_missing = TRUE
  )
  second <- scenario_fit("model", {
    counter <- counter + 1L
    counter
  })
  expect_failure(scenario_text("result", "new"), "candidate is cached")
  expect_equal(readLines(path, warn = FALSE), '[1] "old"')

  scenario_start(
    "unit",
    root           = root,
    refit          = FALSE,
    update         = TRUE,
    create_missing = TRUE
  )
  third <- scenario_fit("model", {
    counter <- counter + 1L
    counter
  })
  scenario_text("result", "new")

  expect_identical(first, 1L)
  expect_identical(second, 2L)
  expect_identical(third, 2L)
  expect_equal(readLines(path, warn = FALSE), '[1] "new"')
})


test_that("scenario_fit does not test cached production timing", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env <- environment(scenario_fit)
  old_clock  <- get(".scenario_clock", envir = helper_env, inherits = FALSE)
  on.exit(assign(".scenario_clock", old_clock, envir = helper_env), add = TRUE)
  ticks <- c(10, 13)
  assign(
    ".scenario_clock",
    function() {

      value <- ticks[[1L]]
      ticks <<- ticks[-1L]
      value
    },
    envir = helper_env
  )
  counter <- 0L

  scenario_start(
    "unit",
    root           = root,
    update_timings = TRUE,
    create_missing = TRUE
  )
  first <- scenario_fit("timed", {
    counter <- counter + 1L
    counter
  })
  expect_no_warning(.scenario_finalize_timing())
  timing_path <- file.path(root, "timings", "unit.tsv")
  fit_metadata <- readRDS(file.path(
    root, "cache", "unit", "timed.timing.rds"
  ))
  fit_metadata[["elapsed"]] <- 30
  saveRDS(fit_metadata, file.path(
    root, "cache", "unit", "timed.timing.rds"
  ))

  scenario_start("unit", root = root, create_missing = TRUE)
  second <- scenario_fit("timed", {
    counter <- counter + 1L
    counter
  })
  .scenario_register_timing("text", "summary", 2)
  cached_timing <- .scenario_current_timings()
  expect_no_warning(.scenario_finalize_timing())

  expect_identical(first, 1L)
  expect_identical(second, 1L)
  expect_identical(counter, 1L)
  expect_equal(fit_metadata[["elapsed"]], 30)
  expect_equal(cached_timing[c("type", "name", "elapsed")], data.frame(
    type = "text", name = "summary", elapsed = 2
  ))
  expect_equal(
    .scenario_read_timings(timing_path)[c("type", "name", "elapsed")],
    data.frame(
      type    = c("fit", "text"),
      name    = c("timed", "summary"),
      elapsed = c(3, 2)
    )
  )
})


test_that("scenario_fit automatically stores post-fit phase timings", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env <- environment(scenario_fit)
  old_clock  <- get(".scenario_clock", envir = helper_env, inherits = FALSE)
  on.exit(assign(".scenario_clock", old_clock, envir = helper_env), add = TRUE)
  ticks <- c(0, 2, 5, 6, 10, 12)
  assign(
    ".scenario_clock",
    function() {

      value <- ticks[[1L]]
      ticks <<- ticks[-1L]
      value
    },
    envir = helper_env
  )
  counter <- 0L
  original_add_loo <- function(object) {

    object[["loo"]] <- TRUE
    object
  }
  original_add_marglik <- function(object) {

    object[["marglik"]] <- TRUE
    object
  }
  add_loo     <- original_add_loo
  add_marglik <- original_add_marglik

  scenario_start(
    "unit",
    root           = root,
    update_timings = TRUE,
    create_missing = TRUE
  )
  first <- scenario_fit("phased", {
    counter <- counter + 1L
    fit <- list(value = counter)
    fit <- add_loo(fit)
    fit <- add_marglik(fit)
    return(fit)
  })
  current <- .scenario_current_timings()
  expect_no_warning(.scenario_finalize_timing())
  metadata <- readRDS(file.path(
    root, "cache", "unit", "phased.timing.rds"
  ))

  expect_identical(add_loo, original_add_loo)
  expect_identical(add_marglik, original_add_marglik)
  expect_identical(metadata[["version"]], 2L)
  expect_equal(metadata[["elapsed"]], 12)
  expect_equal(metadata[["phases"]], c(model = 5, loo = 3, marglik = 4))
  expect_equal(current[c("type", "name", "elapsed")], data.frame(
    type = c("fit", "fit_model", "fit_loo", "fit_marglik"),
    name = rep("phased", 4L),
    elapsed = c(12, 5, 3, 4)
  ))

  scenario_start("unit", root = root, create_missing = TRUE)
  second <- scenario_fit("phased", {
    counter <- counter + 1L
    fit <- list(value = counter)
    fit <- add_loo(fit)
    fit <- add_marglik(fit)
    return(fit)
  })
  cached <- .scenario_current_timings()
  expect_no_warning(.scenario_finalize_timing())

  expect_identical(counter, 1L)
  expect_identical(second, first)
  expect_equal(cached, .scenario_empty_timings())
})


test_that("cached fits do not require timing metadata", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  counter <- 0L

  scenario_start("unit", root = root, create_missing = TRUE)
  scenario_fit("old-cache", {
    counter <- counter + 1L
    counter
  })
  unlink(file.path(root, "cache", "unit", "old-cache.timing.rds"))

  scenario_start("unit", root = root, update_timings = TRUE)
  scenario_fit("old-cache", {
    counter <- counter + 1L
    counter
  })
  .scenario_register_timing("text", "summary", 2)
  expect_no_warning(.scenario_finalize_timing())

  expect_identical(counter, 1L)
  available <- .scenario_read_timings(file.path(root, "timings", "unit.tsv"))
  expect_equal(available[c("type", "name", "elapsed")], data.frame(
    type = "text", name = "summary", elapsed = 2
  ))
})


test_that("scenario timings backfill and retain the fastest baseline", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 10)
  expect_no_warning(.scenario_finalize_timing())
  timing_path <- file.path(root, "timings", "unit.tsv")
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 10)

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 12.1)
  .scenario_register_timing("text", "summary", 10)
  expect_warning(
    .scenario_finalize_timing(),
    "+2.1 s (+21%; 10.0 s -> 12.1 s) fit/model",
    fixed = TRUE
  )
  backfilled <- .scenario_read_timings(timing_path)
  expect_equal(backfilled[["elapsed"]], c(10, 10))
  expect_true(file.exists(file.path(root, "timings", "unit.new.tsv")))

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 8)
  .scenario_register_timing("text", "summary", 8)
  expect_no_warning(.scenario_finalize_timing())
  improved <- .scenario_read_timings(timing_path)
  expect_equal(improved[["elapsed"]], c(8, 8))
  expect_false(file.exists(file.path(root, "timings", "unit.new.tsv")))

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 10)
  .scenario_register_timing("text", "summary", 8)
  warning_message <- character()
  withCallingHandlers(
    .scenario_finalize_timing(),
    warning = function(warning) {

      warning_message <<- conditionMessage(warning)
      invokeRestart("muffleWarning")
    }
  )
  expect_match(
    warning_message,
    "+2.0 s (+25%; 8.0 s -> 10.0 s) fit/model",
    fixed = TRUE
  )
  expect_match(warning_message, "average timing regression: 12%", fixed = TRUE)
  retained <- .scenario_read_timings(timing_path)
  expect_equal(retained[["elapsed"]], c(8, 8))
  expect_true(file.exists(file.path(root, "timings", "unit.new.tsv")))

  scenario_start("unit", root = root, update_timings = TRUE)
  .scenario_register_timing("fit", "model", 10)
  .scenario_register_timing("text", "summary", 8)
  expect_warning(
    .scenario_finalize_timing(),
    "average timing regression: 12%"
  )
  accepted <- .scenario_read_timings(timing_path)
  expect_equal(accepted[["elapsed"]], c(10, 8))
  expect_false(file.exists(file.path(root, "timings", "unit.new.tsv")))

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 10)
  .scenario_register_timing("text", "summary", 9.6)
  warning_message <- character()
  withCallingHandlers(
    .scenario_finalize_timing(),
    warning = function(warning) {

      warning_message <<- conditionMessage(warning)
      invokeRestart("muffleWarning")
    }
  )
  expect_match(warning_message, "average timing regression: 10%", fixed = TRUE)
  expect_false(grepl("text/summary:", warning_message, fixed = TRUE))

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 10.4)
  .scenario_register_timing("text", "summary", 8.4)
  expect_no_warning(.scenario_finalize_timing())
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], c(10, 8))
  expect_true(file.exists(file.path(root, "timings", "unit.new.tsv")))
})


test_that("scenario timings below three quarters of a second are not assessed", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root)
  .scenario_register_timing("text", "quick", 0.06)
  expect_no_warning(.scenario_finalize_timing())

  scenario_start("unit", root = root)
  .scenario_register_timing("text", "quick", 0.74)
  expect_no_warning(.scenario_finalize_timing())
  expect_equal(
    .scenario_read_timings(file.path(root, "timings", "unit.tsv"))[["elapsed"]],
    0.06
  )
  expect_equal(
    .scenario_read_timings(file.path(root, "timings", "unit.new.tsv"))[["elapsed"]],
    0.74
  )

  scenario_start("unit", root = root)
  .scenario_register_timing("text", "quick", 0.75)
  expect_warning(
    .scenario_finalize_timing(),
    "+0.7 s (+1150%; 0.1 s -> 0.8 s) text/quick",
    fixed = TRUE
  )
})


test_that("split fit timing averages exclude the redundant total", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 12)
  .scenario_register_timing("fit_model", "model", 5)
  .scenario_register_timing("fit_loo", "model", 3)
  .scenario_register_timing("fit_marglik", "model", 4)
  expect_no_warning(.scenario_finalize_timing())

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 18)
  .scenario_register_timing("fit_model", "model", 5.5)
  .scenario_register_timing("fit_loo", "model", 3.3)
  .scenario_register_timing("fit_marglik", "model", 4.4)
  warning_message <- character()
  withCallingHandlers(
    .scenario_finalize_timing(),
    warning = function(warning) {

      warning_message <<- conditionMessage(warning)
      invokeRestart("muffleWarning")
    }
  )

  expect_match(
    warning_message,
    "+6.0 s (+50%; 12.0 s -> 18.0 s) fit/model",
    fixed = TRUE
  )
  expect_match(warning_message, "average timing regression: 10%", fixed = TRUE)
  expect_match(warning_message, "unweighted mean across 3 calls", fixed = TRUE)
})


test_that("large absolute scenario timing regressions are highlighted in red", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  old_colors <- getOption("cli.num_colors")
  on.exit(options(cli.num_colors = old_colors), add = TRUE)
  options(cli.num_colors = 8L)

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "large", 10)
  .scenario_register_timing("fit", "boundary", 8)
  expect_no_warning(.scenario_finalize_timing())

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "large", 12.1)
  .scenario_register_timing("fit", "boundary", 10)
  warning_message <- character()
  withCallingHandlers(
    .scenario_finalize_timing(),
    warning = function(warning) {

      warning_message <<- conditionMessage(warning)
      invokeRestart("muffleWarning")
    }
  )

  expect_match(
    warning_message,
    cli::col_red("+2.1 s (+21%; 10.0 s -> 12.1 s) fit/large"),
    fixed = TRUE
  )
  expect_false(grepl(
    cli::col_red("+2.0 s (+25%; 8.0 s -> 10.0 s) fit/boundary"),
    warning_message,
    fixed = TRUE
  ))
  expect_match(
    warning_message,
    "+2.0 s (+25%; 8.0 s -> 10.0 s) fit/boundary",
    fixed = TRUE
  )
})


test_that("scenario_text adds, compares, and regenerates tracked output", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root, create_missing = TRUE, width = 80L)
  expected_value <- data.frame(value = 1:2)
  value <- scenario_text("summary", expected_value)
  path <- file.path(root, "results", "unit", "summary.txt")

  expect_identical(value, expected_value)
  expect_true(file.exists(path))
  expect_equal(readLines(path, warn = FALSE), c("  value", "1     1", "2     2"))
  expect_equal(.scenario_current_timings()[c("type", "name")], data.frame(
    type = "text", name = "summary"
  ))

  explicit_value <- scenario_text("summary", print(expected_value))
  expect_identical(explicit_value, expected_value)
  expect_failure(
    scenario_text("summary", "changed"),
    "candidate is cached"
  )
  expect_true(file.exists(file.path(
    root, "results", "unit", "summary.new.txt"
  )))

  scenario_start(
    "unit",
    root           = root,
    regenerate     = TRUE,
    create_missing = TRUE
  )
  scenario_text("summary", "changed")
  expect_equal(readLines(path, warn = FALSE), "[1] \"changed\"")
  expect_false(file.exists(file.path(
    root, "results", "unit", "summary.new.txt"
  )))
})


test_that("scenario_text seeds stochastic snapshots independently", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root, create_missing = TRUE)
  first <- scenario_text("random", stats::runif(3L))
  stats::runif(20L)
  second <- scenario_text("random", stats::runif(3L))

  expect_identical(first, second)
})


test_that("scenario review combines table and figure decisions", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env   <- environment(scenario_text)
  helper_names <- c(".scenario_is_interactive", ".scenario_snapshot_review")
  old_helpers  <- mget(helper_names, envir = helper_env, inherits = FALSE)
  on.exit(list2env(old_helpers, envir = helper_env), add = TRUE)
  assign(".scenario_is_interactive", function() TRUE, envir = helper_env)
  state <- new.env(parent = emptyenv())
  state[["reviews"]] <- 0L
  assign(
    ".scenario_snapshot_review",
    function(path, files = NULL) {

      state[["reviews"]] <- state[["reviews"]] + 1L
      candidates <- list.files(
        file.path(path, "_snaps"),
        pattern    = "[.]new[.](txt|svg)$",
        recursive  = TRUE,
        full.names = TRUE
      )
      accept <- candidates[grepl("accept.*[.]new[.](txt|svg)$", candidates)]
      reject <- candidates[grepl("reject.*[.]new[.](txt|svg)$", candidates)]
      file.copy(
        accept,
        sub("[.]new[.]", ".", accept),
        overwrite = TRUE
      )
      unlink(accept)
      unlink(reject)
      return(TRUE)
    },
    envir = helper_env
  )

  scenario_start(
    "00-scenario-helpers",
    root           = root,
    update         = FALSE,
    show_output    = FALSE,
    create_missing = TRUE
  )
  scenario_text("accept", "old accept")
  scenario_text("reject", "old reject")
  plot_path <- file.path(
    root, "_snaps", "00-scenario-helpers", "accept-plot.svg"
  )
  .scenario_write_lines("old plot", plot_path)
  .scenario_write_lines(
    "new plot",
    .scenario_candidate_path(plot_path)
  )
  accept_path <- file.path(
    root, "results", "00-scenario-helpers", "accept.txt"
  )
  reject_path <- file.path(
    root, "results", "00-scenario-helpers", "reject.txt"
  )

  expect_failure(
    scenario_text("accept", "new accept"),
    "candidate is cached"
  )
  expect_failure(
    scenario_text("reject", "new reject"),
    "candidate is cached"
  )
  expect_identical(state[["reviews"]], 0L)
  expect_equal(readLines(accept_path, warn = FALSE), '[1] "old accept"')
  expect_equal(readLines(reject_path, warn = FALSE), '[1] "old reject"')

  changes <- expect_message(
    .scenario_review_snapshots(root, "00-scenario-helpers"),
    "Scenario snapshot mismatches"
  )
  expect_identical(state[["reviews"]], 1L)
  expect_equal(
    changes[c("type", "status")],
    data.frame(
      type   = c("table", "table", "figure"),
      status = c("accepted", "rejected", "accepted")
    )
  )
  expect_equal(readLines(accept_path, warn = FALSE), '[1] "new accept"')
  expect_equal(readLines(reject_path, warn = FALSE), '[1] "old reject"')
  expect_equal(readLines(plot_path, warn = FALSE), "new plot")
  expect_false(file.exists(.scenario_candidate_path(accept_path)))
  expect_false(file.exists(.scenario_candidate_path(reject_path)))
  expect_false(file.exists(.scenario_candidate_path(plot_path)))
})


test_that("scenario review tolerates selected scenarios without candidates", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  changes <- expect_message(
    .scenario_review_snapshots(root, c("empty-one", "empty-two")),
    "No cached scenario snapshots"
  )
  expect_equal(nrow(changes), 0L)
})


test_that("review_scenario_snapshots uses and filters the latest run", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env <- environment(review_scenario_snapshots)
  old_review <- get(
    ".scenario_review_snapshots",
    envir    = helper_env,
    inherits = FALSE
  )
  on.exit(assign(
    ".scenario_review_snapshots",
    old_review,
    envir = helper_env
  ), add = TRUE)
  had_last_run <- exists(
    "last_scenario_run", envir = .scenario_state, inherits = FALSE
  )
  if (had_last_run) {
    old_last_run <- get(
      "last_scenario_run", envir = .scenario_state, inherits = FALSE
    )
    on.exit(assign(
      "last_scenario_run",
      old_last_run,
      envir = .scenario_state
    ), add = TRUE)
  } else {
    on.exit(rm("last_scenario_run", envir = .scenario_state), add = TRUE)
  }

  state <- new.env(parent = emptyenv())
  assign(
    "last_scenario_run",
    list(root = root, scenarios = c("alpha", "beta")),
    envir = .scenario_state
  )
  assign(
    ".scenario_review_snapshots",
    function(root, scenarios, types = c("table", "figure")) {

      state[["root"]]      <- root
      state[["scenarios"]] <- scenarios
      state[["types"]]     <- types
      return(invisible(data.frame()))
    },
    envir = helper_env
  )

  review_scenario_snapshots(filter = "beta")
  expect_identical(state[["root"]], root)
  expect_identical(state[["scenarios"]], "beta")
  expect_identical(state[["types"]], c("table", "figure"))
})


test_that("review_test_snapshots forwards testthat snapshot selection", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env   <- environment(review_test_snapshots)
  helper_names <- c(".scenario_is_interactive", ".scenario_snapshot_review")
  old_helpers  <- mget(helper_names, envir = helper_env, inherits = FALSE)
  on.exit(list2env(old_helpers, envir = helper_env), add = TRUE)
  state <- new.env(parent = emptyenv())
  assign(".scenario_is_interactive", function() TRUE, envir = helper_env)
  assign(
    ".scenario_snapshot_review",
    function(path, files = NULL) {

      state[["path"]]  <- path
      state[["files"]] <- files
      return(TRUE)
    },
    envir = helper_env
  )

  reviewed <- review_test_snapshots(files = "test-plot/", root = root)
  expect_true(reviewed)
  expect_identical(
    state[["path"]],
    normalizePath(root, winslash = "/", mustWork = TRUE)
  )
  expect_identical(state[["files"]], "test-plot/")
})


test_that("review_test_references discovers and filters text candidates", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env <- environment(review_test_references)
  old_review <- get(
    ".scenario_review_snapshots",
    envir    = helper_env,
    inherits = FALSE
  )
  on.exit(assign(
    ".scenario_review_snapshots",
    old_review,
    envir = helper_env
  ), add = TRUE)
  .scenario_write_lines(
    "old",
    file.path(root, "results", "interpret", "output.txt")
  )
  .scenario_write_lines(
    "new",
    file.path(root, "results", "interpret", "output.new.txt")
  )
  .scenario_write_lines(
    "new",
    file.path(root, "results", "marginal_means", "table.new.txt")
  )

  state <- new.env(parent = emptyenv())
  assign(
    ".scenario_review_snapshots",
    function(root, scenarios, types = c("table", "figure"),
             subject = "Scenario") {

      state[["root"]]      <- root
      state[["scenarios"]] <- scenarios
      state[["types"]]     <- types
      state[["subject"]]   <- subject
      return(invisible(data.frame()))
    },
    envir = helper_env
  )

  review_test_references(filter = "interpret", root = root)
  expect_identical(
    state[["root"]],
    normalizePath(root, winslash = "/", mustWork = TRUE)
  )
  expect_identical(state[["scenarios"]], "interpret")
  expect_identical(state[["types"]], "table")
  expect_identical(state[["subject"]], "Test reference")
})


test_that("scenario_text rejects missing baselines when creation is disabled", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start("unit", root = root, create_missing = FALSE)
  expect_failure(
    scenario_text("missing", "output"),
    "Missing locked scenario text"
  )
  expect_false(file.exists(file.path(
    root, "results", "unit", "missing.txt"
  )))
})


test_that("scenario_text replays captured output when requested", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  scenario_start(
    "unit",
    root           = root,
    show_output    = TRUE,
    create_missing = TRUE
  )
  expect_output(
    scenario_text("visible", "interactive output"),
    '[1] "interactive output"',
    fixed = TRUE
  )
})


test_that("scenario_text suppresses evaluation messages", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  scenario_start("unit", root = root, create_missing = FALSE)
  .scenario_write_lines(
    '[1] "result"',
    file.path(root, "results", "unit", "quiet.txt")
  )

  expect_no_message(scenario_text("quiet", {
    message("JAGS fitting progress")
    "result"
  }))
  expect_warning(
    scenario_text("quiet", {
      warning("important warning")
      "result"
    }),
    "important warning"
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


test_that("scenario_plot locks interactive plots and seeds each rendering", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper_env   <- environment(scenario_plot)
  helper_names <- ".scenario_snapshot_context"
  old_helpers <- mget(helper_names, envir = helper_env, inherits = FALSE)
  on.exit(list2env(old_helpers, envir = helper_env), add = TRUE)
  assign(".scenario_snapshot_context", function() NULL, envir = helper_env)
  state <- new.env(parent = emptyenv())
  state[["draws"]] <- 0L
  state[["random"]] <- numeric()
  state[["preview_device"]] <- logical()

  scenario_start(
    "unit",
    root           = root,
    show_output    = TRUE,
    create_missing = TRUE
  )
  .with_temp_plot_device({
    preview_device <- grDevices::dev.cur()
    graphics::plot(1:3, 1:3)
    graphics::par(new = TRUE)
    scenario_plot("mu_BF_comparison", {
      state[["entry_new"]] <- graphics::par("new")
      state[["draws"]] <- state[["draws"]] + 1L
      state[["random"]] <- c(state[["random"]], stats::runif(1L))
      state[["preview_device"]] <- c(
        state[["preview_device"]],
        identical(grDevices::dev.cur(), preview_device)
      )
      graphics::plot(1:3, 1:3)
      graphics::par(new = TRUE)
    })
    state[["exit_new"]] <- graphics::par("new")
  })

  set.seed(1)
  expected_random <- stats::runif(1L)
  expect_identical(state[["draws"]], 2L)
  expect_equal(state[["random"]], rep(expected_random, 2L))
  expect_identical(state[["preview_device"]], c(FALSE, TRUE))
  expect_false(state[["entry_new"]])
  expect_false(state[["exit_new"]])
  expect_true(file.exists(file.path(
    root, "_snaps", "unit", "mu-bf-comparison.svg"
  )))

  .with_temp_plot_device({
    expect_failure(
      scenario_plot("mu_BF_comparison", {
        graphics::plot(1:4, 1:4)
      }),
      "Scenario plot.*changed"
    )
  })
  expect_true(file.exists(file.path(
    root, "_snaps", "unit", "mu-bf-comparison.new.svg"
  )))
  expect_equal(.scenario_current_timings()[c("type", "name")], data.frame(
    type = "plot", name = "mu_BF_comparison"
  ))
})


test_that("scenario_plot delegates plot comparison to vdiffr", {

  skip_if_not_installed("vdiffr")
  skip_if_no_vdiffr_snapshots()
  scenario_start(
    "00-scenario-helpers",
    root           = testthat::test_path(),
    create_missing = TRUE
  )
  .with_temp_plot_device({
    scenario_plot("simple-base-plot", {
      graphics::plot(1:3, 1:3, xlab = "x", ylab = "y")
    })
  })
})


test_that("scenario_agreement_plot validates and supports boundary cases", {

  expect_error(
    scenario_agreement_plot(1:2, 1),
    "same length"
  )
  expect_error(
    scenario_agreement_plot(c(NA_real_, Inf), c(1, 2)),
    "unavailable: no finite"
  )
  expect_warning(
    .with_temp_plot_device({
      scenario_agreement_plot(
        c(1, NA_real_), c(1, 2),
        reference_label = "baseline",
        estimate_label  = "candidate"
      )
    }),
    "1 non-finite reference-estimate pair"
  )
  expect_no_error(
    .with_temp_plot_device({
      scenario_agreement_plot(1, 1)
    })
  )
})


test_that("scenario estimate extractors select named metafor and RoBMA values", {

  local_mocked_s3_method(
    "summary", "scenario_metafor_fit",
    function(object, ...) structure(object, class = "scenario_metafor_summary")
  )
  local_mocked_s3_method(
    "coef", "scenario_metafor_summary",
    function(object, ...) {
      matrix(
        c(0.25, 0.10, 0.05, 0.45,
          0.40, 0.15, 0.10, 0.70),
        nrow = 2L, byrow = TRUE,
        dimnames = list(
          c("intrcpt", "groupB"),
          c("estimate", "se", "ci.lb", "ci.ub")
        )
      )
    }
  )
  metafor_fit <- structure(
    list(
      sigma2   = c(1, 4),
      s.names  = c("study", "study/esid"),
      tau2     = c(9, 16),
      g.levels.f = list(c("a", "b"), as.character(1:3)),
      rho      = c(0.3, 0.4)
    ),
    class = c("scenario_metafor_fit", "rma")
  )

  expect_equal(ex_m(metafor_fit, "mu"), 0.25)
  expect_equal(ex_m(metafor_fit, "intercept", "SE"), 0.10)
  expect_equal(ex_m(metafor_fit, "group[B]", "CI_0.975"), 0.70)
  expect_equal(ex_m(metafor_fit, "sigma[study]"), 1)
  expect_equal(ex_m(metafor_fit, "sigma[total]^2"), 5)
  expect_equal(ex_m(metafor_fit, "tau[b]^2"), 16)
  expect_equal(ex_m(metafor_fit, "rho[2]"), 0.4)
  expect_error(ex_m(metafor_fit, "sigma"), "ambiguous")
  expect_error(ex_m(metafor_fit, "rho"), "ambiguous")
  expect_true(is.na(ex_m(metafor_fit, "sigma[missing]")))
  expect_equal(
    ex_m(metafor_fit, c("mu", study_sd = "sigma[study]", absent = "gamma")),
    c(mu = 0.25, study_sd = 1, absent = NA_real_)
  )

  local_mocked_s3_method(
    "summary", "scenario_robma_fit",
    function(object, ...) data.frame(
      component = c("location", "random"),
      parameter = c("intercept", "sd_total"),
      Mean       = c(0.30, 0.50),
      SD         = c(0.08, 0.12),
      CI_0.025   = c(0.14, 0.28),
      CI_0.975   = c(0.46, 0.75)
    )
  )
  helper_env <- environment(ex_r)
  had_local_summary_heterogeneity <- exists(
    "summary_heterogeneity",
    envir    = helper_env,
    inherits = FALSE
  )
  if (had_local_summary_heterogeneity) {
    old_summary_heterogeneity <- get(
      "summary_heterogeneity",
      envir    = helper_env,
      inherits = FALSE
    )
  }
  on.exit({
    if (had_local_summary_heterogeneity) {
      assign(
        "summary_heterogeneity",
        old_summary_heterogeneity,
        envir = helper_env
      )
    } else {
      rm("summary_heterogeneity", envir = helper_env)
    }
  }, add = TRUE)
  assign(
    "summary_heterogeneity",
    function(fit, ...) {

      if (inherits(fit, "scenario_robma_single_fit")) {
        return(data.frame(
          component = c("study", "study"),
          parameter = c("sd", "var"),
          Mean       = c(0.25, 0.0625),
          Median     = c(0.24, 0.0576)
        ))
      }
      data.frame(
        component = c("study", "observation", "total"),
        parameter = c("sd", "sd", "var_total"),
        Mean       = c(0.20, 0.40, 0.20),
        Median     = c(0.18, 0.38, 0.19)
      )
    },
    envir = helper_env
  )
  robma_fit <- structure(list(), class = c("scenario_robma_fit", "brma"))

  expect_equal(ex_r(robma_fit, "mu"), 0.30)
  expect_equal(ex_r(robma_fit, "intercept", statistic = "SD"), 0.08)
  expect_equal(ex_r(robma_fit, "sd", "study"), 0.20)
  expect_equal(ex_r(robma_fit, "sd", "study", "Median"), 0.18)
  expect_equal(ex_r(robma_fit, "var_total"), 0.20)
  expect_error(ex_r(robma_fit, "sd"), "ambiguous")
  expect_true(is.na(ex_r(robma_fit, "missing")))
  expect_equal(
    ex_r(
      robma_fit,
      c(mu = "mu", study = "sd", observation = "sd", absent = "missing"),
      component = c(NA, "study", "observation", NA)
    ),
    c(mu = 0.30, study = 0.20, observation = 0.40, absent = NA_real_)
  )

  local_mocked_s3_method(
    "summary", "scenario_robma_single_fit",
    function(object, ...) data.frame(
      component = "location",
      parameter = "intercept",
      Mean       = 0.35
    )
  )
  robma_single_fit <- structure(list(), class = c("scenario_robma_single_fit", "brma"))
  expect_equal(
    ex_r(
      robma_single_fit,
      c(total = "sd", study = "sd", observation = "sd"),
      component = c(NA, "study", "observation")
    ),
    c(total = 0.25, study = 0.25, observation = NA_real_)
  )
  expect_equal(
    ex_r(
      robma_single_fit,
      c(study = "sd", observation = "sd", total = "sd"),
      component = c("study", "observation", "study")
    ),
    c(study = 0.25, observation = NA_real_, total = 0.25)
  )
  expect_error(ex_r(robma_single_fit, "sd", statistic = "SD"), "Statistic")

  local_mocked_s3_method(
    "summary", "scenario_robma_mu_fit",
    function(object, ...) data.frame(
      component = "common",
      parameter = "mu",
      Mean       = 0.35
    )
  )
  robma_mu_fit <- structure(list(), class = c("scenario_robma_mu_fit", "brma"))
  expect_equal(ex_r(robma_mu_fit, "intercept"), 0.35)

  expect_equal(
    ex(
      robma_fit,
      c(mu = "mu", study = "sd", observation = "sd", absent = "missing"),
      component = c(NA, "study", "observation", NA)
    ),
    c(mu = 0.30, study = 0.20, observation = 0.40, absent = NA_real_)
  )
  expect_equal(
    ex(list(metafor = metafor_fit, RoBMA = robma_fit), "mu"),
    c(metafor = 0.25, RoBMA = 0.30)
  )
  expect_equal(
    ex(
      list(metafor = metafor_fit, RoBMA = robma_fit),
      c(mu = "mu", absent = "missing")
    ),
    data.frame(
      mu = c(0.25, 0.30), absent = c(NA_real_, NA_real_),
      row.names = c("metafor", "RoBMA")
    )
  )
  expect_error(ex(robma_fit, "sd"), "ambiguous")
  expect_error(ex(robma_fit, "intercept", statistic = "missing"), "Statistic")
  expect_error(ex(1, "mu"), "metafor model")
})


test_that("ex_p selects stable pooled-effect comparison statistics", {

  local_mocked_s3_method(
    "pooled_effect", "scenario_pooled_fit",
    function(object, ...) data.frame(
      component = "location",
      parameter = "mu",
      Mean      = 0.30,
      Median    = 0.29,
      CI_0.025  = 0.14,
      CI_0.975  = 0.46,
      PI_0.025  = -0.20,
      PI_0.975  = 0.80
    )
  )
  fit <- structure(list(), class = "scenario_pooled_fit")

  expect_equal(
    ex_p(fit),
    data.frame(
      mu        = c(0.30, 0.14, 0.46, -0.20, 0.80),
      row.names = c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")
    )
  )
})


test_that("plot_marginal_diagnostics compares the shared marginal targets", {

  state <- new.env(parent = emptyenv())
  state[["plots"]] <- list()

  local_mocked_s3_method(
    "residuals", "scenario_diagnostic_fit",
    function(object, ...) {

      if (identical(object[["role"]], "estimate")) {
        state[["residual_args"]] <- list(...)
      }
      return(object[["offset"]] + c(1, 2))
    }
  )
  local_mocked_s3_method(
    "rstandard", "scenario_diagnostic_fit",
    function(model, ...) {

      if (identical(model[["role"]], "estimate")) {
        state[["rstandard_args"]] <- list(...)
        warning("expected diagnostic warning")
      }
      return(list(z = model[["offset"]] + c(3, 4)))
    }
  )
  local_mocked_s3_method(
    "hatvalues", "scenario_diagnostic_fit",
    function(model, ...) {

      if (identical(model[["role"]], "estimate")) {
        warning("expected diagnostic warning")
      }
      return(model[["offset"]] + c(5, 6))
    }
  )
  local_mocked_s3_method(
    "cooks.distance", "scenario_diagnostic_fit",
    function(model, ...) model[["offset"]] + c(7, 8)
  )
  local_mocked_s3_method(
    "dfbetas", "scenario_diagnostic_fit",
    function(model, ...) matrix(model[["offset"]] + c(9, 10), ncol = 1L)
  )

  helper_env <- environment(plot_marginal_diagnostics)
  old_agreement_plot <- get(
    "scenario_agreement_plot",
    envir    = helper_env,
    inherits = FALSE
  )
  on.exit(assign(
    "scenario_agreement_plot",
    old_agreement_plot,
    envir = helper_env
  ), add = TRUE)
  assign(
    "scenario_agreement_plot",
    function(reference, estimate, main = "", ...) {

      state[["plots"]][[main]] <- list(
        reference = reference,
        estimate  = estimate
      )
      return(invisible(NULL))
    },
    envir = helper_env
  )

  reference <- structure(
    list(role = "reference", offset = 0),
    class = "scenario_diagnostic_fit"
  )
  estimate <- structure(
    list(role = "estimate", offset = 10),
    class = "scenario_diagnostic_fit"
  )
  expect_no_warning(.with_temp_plot_device({
    value <- plot_marginal_diagnostics(reference, estimate)
    state[["mfrow"]] <- graphics::par("mfrow")
  }))

  expect_null(value)
  expect_identical(
    names(state[["plots"]]),
    c("Residuals", "Rstandard", "Hat values", "Cooks distance", "DFBETAS")
  )
  expected_reference <- list(
    Residuals        = c(1, 2),
    Rstandard        = c(3, 4),
    `Hat values`     = c(5, 6),
    `Cooks distance` = c(7, 8),
    DFBETAS          = matrix(c(9, 10), ncol = 1L)
  )
  expect_equal(
    lapply(state[["plots"]], `[[`, "reference"),
    expected_reference
  )
  expect_equal(
    lapply(state[["plots"]], `[[`, "estimate"),
    lapply(expected_reference, function(x) x + 10)
  )
  expect_identical(
    state[["residual_args"]],
    list(type = "outcome", conditioning_depth = "marginal")
  )
  expect_identical(
    state[["rstandard_args"]],
    list(conditioning_depth = "marginal")
  )
  expect_identical(state[["mfrow"]], c(3L, 2L))
})


test_that("scenario orphan audit reports only unreferenced locked output", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  scenario_file <- file.path(root, "test-example.R")
  .scenario_write_lines(c(
    "values <- matrix(1)",
    "values[, 1]",
    'scenario_text("kept-text", 1)',
    'scenario_plot("kept-plot", plot(1))'
  ), scenario_file)

  .scenario_write_lines(
    "kept",
    file.path(root, "results", "example", "kept-text.txt")
  )
  .scenario_write_lines(
    "candidate",
    file.path(root, "results", "example", "kept-text.new.txt")
  )
  .scenario_write_lines(
    "orphan",
    file.path(root, "results", "example", "orphan-text.txt")
  )
  .scenario_write_lines(
    "kept",
    file.path(root, "_snaps", "example", "kept-plot.svg")
  )
  .scenario_write_lines(
    "orphan",
    file.path(root, "_snaps", "example", "orphan-plot.svg")
  )

  orphans <- expect_message(
    .scenario_report_orphans(scenario_file),
    "orphan-plot[.]svg"
  )
  expect_equal(
    orphans,
    c(
      "_snaps/example/orphan-plot.svg",
      "results/example/orphan-text.txt"
    )
  )
})


test_that("test_scenario filters scenario names and forwards run controls", {

  root <- .scenario_test_root()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  scenario_root <- file.path(root, "tests", "scenarios")
  dir.create(scenario_root, recursive = TRUE)
  .scenario_write_lines("", file.path(scenario_root, "test-alpha.R"))
  .scenario_write_lines("", file.path(scenario_root, "test-beta.R"))

  helper_env <- environment(test_scenario)
  old_test_file <- get(
    ".scenario_test_file",
    envir    = helper_env,
    inherits = FALSE
  )
  old_review <- get(
    "review_scenario_snapshots",
    envir    = helper_env,
    inherits = FALSE
  )
  old_finalize <- get(
    ".scenario_finalize_timing",
    envir    = helper_env,
    inherits = FALSE
  )
  on.exit(assign(".scenario_test_file", old_test_file, envir = helper_env),
          add = TRUE)
  on.exit(assign(
    "review_scenario_snapshots",
    old_review,
    envir = helper_env
  ), add = TRUE)
  on.exit(assign(
    ".scenario_finalize_timing",
    old_finalize,
    envir = helper_env
  ), add = TRUE)
  state <- new.env(parent = emptyenv())
  state[["paths"]]  <- character()
  state[["refit"]]  <- character()
  state[["update"]] <- character()
  state[["update_timings"]] <- character()
  state[["runner"]] <- character()
  state[["reporter"]] <- character()
  state[["review"]] <- character()
  assign(
    ".scenario_test_file",
    function(path, reporter, stop_on_failure) {

      state[["paths"]] <- c(state[["paths"]], basename(path))
      state[["refit"]] <- c(
        state[["refit"]],
        Sys.getenv("ROBMA_SCENARIO_REFIT")
      )
      state[["update"]] <- c(
        state[["update"]],
        Sys.getenv("ROBMA_SCENARIO_UPDATE")
      )
      state[["update_timings"]] <- c(
        state[["update_timings"]],
        Sys.getenv("ROBMA_SCENARIO_UPDATE_TIMINGS")
      )
      state[["runner"]] <- c(
        state[["runner"]],
        Sys.getenv("ROBMA_SCENARIO_RUNNER")
      )
      state[["reporter"]] <- c(state[["reporter"]], reporter)
      return(invisible(list()))
    },
    envir = helper_env
  )
  assign(
    ".scenario_finalize_timing",
    function(scenario, allow_update) {

      state[["finalize"]] <- paste(scenario, allow_update, sep = ":")
      return(invisible(NULL))
    },
    envir = helper_env
  )
  assign(
    "review_scenario_snapshots",
    function(filter = NULL, root = NULL) {

      last_run <- get(
        "last_scenario_run", envir = .scenario_state, inherits = FALSE
      )
      state[["review"]] <- c(
        state[["review"]],
        paste(
          basename(last_run[["root"]]),
          paste(last_run[["scenarios"]], collapse = ","),
          sep = ":"
        )
      )
      return(invisible(data.frame()))
    },
    envir = helper_env
  )

  expect_message(
    test_scenario(
      filter         = "alpha",
      refit          = TRUE,
      update         = FALSE,
      update_timings = TRUE,
      load_package   = FALSE,
      root           = scenario_root
    ),
    "Running scenario 'alpha'"
  )
  expect_identical(state[["paths"]], "test-alpha.R")
  expect_identical(state[["refit"]], "TRUE")
  expect_identical(state[["update"]], "FALSE")
  expect_identical(state[["update_timings"]], "TRUE")
  expect_identical(state[["runner"]], "TRUE")
  expect_identical(state[["reporter"]], "progress")
  expect_identical(state[["finalize"]], "alpha:TRUE")
  expect_identical(state[["review"]], "scenarios:alpha")

  assign(
    ".scenario_test_file",
    function(path, reporter, stop_on_failure) {

      stop("scenario stopped", call. = FALSE)
    },
    envir = helper_env
  )
  expect_error(
    test_scenario(
      filter       = "beta",
      load_package = FALSE,
      root         = scenario_root
    ),
    "scenario stopped"
  )
  expect_identical(state[["finalize"]], "beta:FALSE")
  expect_identical(
    state[["review"]],
    c("scenarios:alpha", "scenarios:beta")
  )

  expect_error(
    test_scenario(
      filter       = "missing",
      load_package = FALSE,
      root         = scenario_root
    ),
    "No scenarios matched"
  )
})
