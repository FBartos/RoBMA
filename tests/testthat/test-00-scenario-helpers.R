context("Maintainer scenario helpers")

source(testthat::test_path("..", "scenarios", "helper-scenarios.R"))


.scenario_test_root <- function() {

  root <- tempfile("robma-scenario-test-")
  dir.create(root, recursive = TRUE)
  return(root)
}


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

  direct <- scenario_start("unit", root = root)
  expect_false(direct[["refit"]])
  expect_true(direct[["show_output"]])
  expect_true(direct[["update"]])
  expect_false(direct[["update_timings"]])
  expect_true(direct[["create_missing"]])
  direct_output <- capture.output({
    scenario_text("result", "old")
    scenario_text("result", "new")
  })
  expect_match(paste(direct_output, collapse = "\n"), '[1] "new"',
               fixed = TRUE)
  path <- file.path(root, "results", "unit", "result.txt")
  expect_equal(readLines(path, warn = FALSE), '[1] "new"')
  .with_temp_plot_device({
    scenario_plot("direct-plot", graphics::plot(1:3, 1:3))
    scenario_plot("direct-plot", graphics::plot(1:4, 1:4))
  })
  plot_path <- file.path(root, "_snaps", "unit", "direct-plot.svg")
  expect_true(file.exists(plot_path))
  expect_false(file.exists(.scenario_candidate_path(plot_path)))

  Sys.setenv(ROBMA_SCENARIO_RUNNER = "TRUE")
  managed <- scenario_start("unit", root = root)
  expect_false(managed[["refit"]])
  expect_false(managed[["show_output"]])
  expect_false(managed[["update"]])
  expect_false(managed[["update_timings"]])
  expect_false(managed[["create_missing"]])
  expect_failure(scenario_text("result", "managed"), "candidate is cached")
  expect_equal(readLines(path, warn = FALSE), '[1] "new"')
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
    "average timing regression: 20.0%"
  )
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 10)
  expect_equal(.scenario_read_timings(candidate_path)[["elapsed"]], 12)

  scenario_start("unit", root = root, update_timings = TRUE)
  expect_warning(
    .scenario_register_timing("text", "summary", 12),
    "average timing regression: 20.0%"
  )
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 12)
  expect_false(file.exists(candidate_path))
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


test_that("scenario_fit reuses cached production timing", {

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

  scenario_start("unit", root = root, create_missing = TRUE)
  second <- scenario_fit("timed", {
    counter <- counter + 1L
    counter
  })
  cached_timing <- .scenario_current_timings()
  expect_no_warning(.scenario_finalize_timing())

  expect_identical(first, 1L)
  expect_identical(second, 1L)
  expect_identical(counter, 1L)
  expect_equal(fit_metadata[["elapsed"]], 3)
  expect_equal(cached_timing[["elapsed"]], 3)
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], 3)
})


test_that("old fit caches require refitting before timing updates", {

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
  expect_warning(
    .scenario_finalize_timing(),
    "cached fit predates valid timing metadata"
  )

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
  expect_warning(.scenario_finalize_timing(), "fit/model: 12.100 s vs 10.000 s")
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
  expect_match(warning_message, "fit/model: 10.000 s vs 8.000 s", fixed = TRUE)
  expect_match(warning_message, "threshold 20%", fixed = TRUE)
  expect_match(warning_message, "average timing regression: 12.5%", fixed = TRUE)
  retained <- .scenario_read_timings(timing_path)
  expect_equal(retained[["elapsed"]], c(8, 8))
  expect_true(file.exists(file.path(root, "timings", "unit.new.tsv")))

  scenario_start("unit", root = root, update_timings = TRUE)
  .scenario_register_timing("fit", "model", 10)
  .scenario_register_timing("text", "summary", 8)
  expect_warning(
    .scenario_finalize_timing(),
    "average timing regression: 12.5%"
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
  expect_match(warning_message, "average timing regression: 10.0%", fixed = TRUE)
  expect_false(grepl("text/summary:", warning_message, fixed = TRUE))

  scenario_start("unit", root = root)
  .scenario_register_timing("fit", "model", 10.4)
  .scenario_register_timing("text", "summary", 8.4)
  expect_no_warning(.scenario_finalize_timing())
  expect_equal(.scenario_read_timings(timing_path)[["elapsed"]], c(10, 8))
  expect_true(file.exists(file.path(root, "timings", "unit.new.tsv")))
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

  scenario_start(
    "unit",
    root           = root,
    show_output    = TRUE,
    create_missing = TRUE
  )
  .with_temp_plot_device({
    graphics::plot(1:3, 1:3)
    graphics::par(new = TRUE)
    scenario_plot("mu_BF_comparison", {
      state[["entry_new"]] <- graphics::par("new")
      state[["draws"]] <- state[["draws"]] + 1L
      state[["random"]] <- c(state[["random"]], stats::runif(1L))
      graphics::plot(1:3, 1:3)
      graphics::par(new = TRUE)
    })
    state[["exit_new"]] <- graphics::par("new")
  })

  set.seed(1)
  expected_random <- stats::runif(1L)
  expect_identical(state[["draws"]], 2L)
  expect_equal(state[["random"]], rep(expected_random, 2L))
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
