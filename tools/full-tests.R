#!/usr/bin/env Rscript

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "full-tests.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

# '--vanilla': the profile runs with this process's environment (its test
# controls such as ROBMA_TEST_FILES_DIR), not with values that user or site
# startup files would set again in the child.
status <- system2(
  command = file.path(R.home("bin"), "Rscript"),
  args    = c(
    "--vanilla",
    shQuote(file.path(project_root, "tools", "test-profile.R")),
    "release"
  )
)

quit(status = status)
