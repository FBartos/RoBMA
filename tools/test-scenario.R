#!/usr/bin/env Rscript

args       <- commandArgs(trailingOnly = TRUE)
regenerate <- "--regenerate" %in% args
list_only  <- "--list" %in% args
requested  <- setdiff(args, c("--regenerate", "--list"))

cmd      <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "test-scenario.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
scenario_dir <- file.path(project_root, "tests", "scenarios")

scenario_files <- list.files(
  scenario_dir,
  pattern    = "^test-[a-z0-9][a-z0-9_-]*\\.[Rr]$",
  full.names = TRUE
)
scenario_names <- sub(
  "^test-([a-z0-9][a-z0-9_-]*)\\.[Rr]$",
  "\\1",
  basename(scenario_files)
)

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

old_regenerate <- Sys.getenv("ROBMA_SCENARIO_REGENERATE", unset = NA_character_)
on.exit({
  if (is.na(old_regenerate)) {
    Sys.unsetenv("ROBMA_SCENARIO_REGENERATE")
  } else {
    Sys.setenv(ROBMA_SCENARIO_REGENERATE = old_regenerate)
  }
}, add = TRUE)
Sys.setenv(
  AGENT                       = "1",
  NOT_CRAN                    = "true",
  ROBMA_SCENARIO_REGENERATE   = if (regenerate) "TRUE" else "FALSE"
)

setwd(project_root)
devtools::load_all(quiet = TRUE)

for (name in requested) {
  path <- scenario_files[match(name, scenario_names)]
  message("Running scenario '", name, "'.")
  testthat::test_file(
    path,
    reporter        = "llm",
    package         = "RoBMA",
    load_package    = "none",
    stop_on_failure = TRUE
  )
}
