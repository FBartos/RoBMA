#!/usr/bin/env Rscript

args    <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) > 0L) args[[1L]] else "standard"
if (!profile %in% c("standard", "certification")) {
  stop("Cache profile must be 'standard' or 'certification'.", call. = FALSE)
}

cmd      <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "test-cache-key.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

Sys.setenv(ROBMA_TEST_PROFILE = profile)
source(file.path("tests", "testthat", "common-functions.R"))

catalog      <- active_fit_catalog()
source_files <- unique(catalog[["source_file"]])
source_hashes <- vapply(source_files, source_file_md5, character(1))
jags_version <- tryCatch(
  as.character(rjags::jags.version()),
  error = function(error) "unavailable"
)
key_data <- c(
  paste0("schema=", FIT_CACHE_SCHEMA_VERSION),
  paste0("profile=", profile),
  paste0(
    "RoBMA=",
    unname(read.dcf("DESCRIPTION", fields = "Version")[1L, 1L])
  ),
  paste0("BayesTools=", as.character(utils::packageVersion("BayesTools"))),
  paste0("R=", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("JAGS=", jags_version),
  paste0("package=", package_source_md5()),
  paste0("BayesTools-backend=", bayestools_backend_fingerprint()),
  paste0(source_files, "=", source_hashes)
)

key_file <- tempfile("robma-test-cache-key-")
on.exit(unlink(key_file), add = TRUE)
writeLines(key_data, key_file, useBytes = TRUE)
cat(profile, "-", unname(tools::md5sum(key_file)), "\n", sep = "")
