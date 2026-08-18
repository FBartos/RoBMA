#!/usr/bin/env Rscript

args       <- commandArgs(trailingOnly = TRUE)
regenerate <- "--regenerate" %in% args
refit      <- regenerate || "--refit" %in% args
update     <- regenerate || "--update" %in% args
update_timings <- update || "--update-timings" %in% args
list_only  <- "--list" %in% args
requested  <- setdiff(
  args,
  c("--regenerate", "--refit", "--update", "--update-timings", "--list")
)

cmd      <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "test-scenario.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
scenario_dir <- file.path(project_root, "tests", "scenarios")
source(file.path(scenario_dir, "helper-scenarios.R"))

scenario_files <- .scenario_list_files(scenario_dir)
scenario_names <- names(scenario_files)

if (list_only) {
  if (length(scenario_names) == 0L) {
    message("No scenarios are defined.")
  } else {
    writeLines(scenario_names)
  }
  quit(save = "no", status = 0L)
}

if (length(requested) == 0L) {
  stop(
    "Specify a scenario name or 'all'. Use --list to list available scenarios.",
    call. = FALSE
  )
}
if (identical(requested, "all")) {
  requested <- scenario_names
}
if (length(requested) == 0L) {
  stop("No scenarios are defined.", call. = FALSE)
}

missing <- setdiff(requested, scenario_names)
if (length(missing) > 0L) {
  stop(
    "Unknown scenario: ", paste(missing, collapse = ", "),
    ". Available scenarios: ", paste(scenario_names, collapse = ", "),
    call. = FALSE
  )
}

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
Sys.setenv(
  AGENT                             = "1",
  NOT_CRAN                          = "true",
  ROBMA_SCENARIO_REGENERATE         = if (regenerate) "TRUE" else "FALSE",
  ROBMA_SCENARIO_REFIT              = if (refit) "TRUE" else "FALSE",
  ROBMA_SCENARIO_UPDATE             = if (update) "TRUE" else "FALSE",
  ROBMA_SCENARIO_UPDATE_TIMINGS     = if (update_timings) "TRUE" else "FALSE",
  ROBMA_SCENARIO_RUNNER             = "TRUE"
)

setwd(project_root)
devtools::load_all(quiet = TRUE)

for (name in requested) {
  path      <- scenario_files[match(name, scenario_names)]
  completed <- FALSE
  message("Running scenario '", name, "'.")
  tryCatch(
    {
      result <- testthat::test_file(
        path,
        reporter        = "llm",
        package         = "RoBMA",
        load_package    = "none",
        stop_on_failure = TRUE
      )
      completed <- TRUE
      result
    },
    finally = {
      .scenario_finalize_timing(name, allow_update = completed)
      .scenario_report_orphans(path)
    }
  )
}
