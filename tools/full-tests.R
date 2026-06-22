#!/usr/bin/env Rscript

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "full-tests.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

status <- system2(
  command = file.path(R.home("bin"), "Rscript"),
  args    = c(file.path(project_root, "tools", "test-profile.R"), "release")
)

quit(status = status)
