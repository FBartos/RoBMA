context("Density profiler controls")


test_that("density profiler controls resolve with their unchanged ordinate budgets", {

  path <- testthat::test_path("..", "..", "tools", "profile-iwmde-qcmde.R")
  skip_if_not(file.exists(path), "Developer profiler requires a source checkout.")
  expressions <- parse(path, keep.source = FALSE)
  definition <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name(".profile_iwmde_estimator"))
  }, as.list(expressions))[[1L]][[3L]]
  find_estimator <- function(expr) {
    if (!is.call(expr)) return(NULL)
    if (identical(expr[[1L]], as.name(".iwmde_estimate"))) return(expr)
    for (part in as.list(expr)) {
      found <- find_estimator(part)
      if (!is.null(found)) return(found)
    }
    return(NULL)
  }
  estimator <- find_estimator(definition)
  controls <- eval(estimator[["density_control"]],
                   list(n_points = 40L, max_samples = 40L,
                        normalization_points = 101L))
  resolved <- tryCatch(RoBMA:::.iwmde_density_control_resolve(
    controls, density_method = "qCMDE", purpose = "ordinate"
  ), error = identity)
  expect_false(inherits(resolved, "error"),
               info = if (inherits(resolved, "error")) conditionMessage(resolved))
  if (!inherits(resolved, "error")) {
    expect_identical(resolved[["samples"]], 40L)
    expect_identical(resolved[["n_points"]], 40L)
    expect_identical(resolved[["normalization_points"]], 101L)
    expect_identical(resolved[["normalization_prob"]], .9999)
    expect_identical(resolved[["display_grid"]], "ordinate")
  }
})


test_that("profiler loader resolves its own root and reports build and omitted cases", {

  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  skip_if_not(file.exists(file.path(package_root, "tools/profile-iwmde-qcmde.R")),
              "Developer profiler requires a source checkout.")
  root <- tempfile("robma-profiler-loader-")
  dir.create(file.path(root, "tools"), recursive = TRUE)
  dir.create(file.path(root, "R"))
  dir.create(file.path(root, "tests", "testthat"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (file in c("tools/profile-iwmde-qcmde.R", "tools/optimized-dll.R")) {
    file.copy(file.path(package_root, file), file.path(root, file))
  }
  writeLines(c("Package: RoBMA", "Version: 0.0.0.9000",
               "Title: Synthetic Profiler Loader Fixture",
               "Description: A minimal profiler loader fixture without native code.",
               "Author: RoBMA tests", "Maintainer: RoBMA tests <tests@example.org>",
               "License: MIT"), file.path(root, "DESCRIPTION"))
  writeLines("export(RoBMA.options,RoBMA.get_option)", file.path(root, "NAMESPACE"))
  writeLines(c('.state <- new.env(parent = emptyenv())', '.state$threads <- NA_integer_',
               'RoBMA.options <- function(native_threads) .state$threads <- native_threads',
               'RoBMA.get_option <- function(name) .state$threads'), file.path(root, "R/options.R"))
  writeLines('list_fits <- function(validate) character()',
             file.path(root, "tests/testthat/common-functions.R"))
  invocation <- tempfile("robma-profiler-invocation-")
  dir.create(invocation)
  on.exit(unlink(invocation, recursive = TRUE), add = TRUE)
  executable <- file.path(R.home("bin"),
                          if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  withr::local_envvar(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  for (debug in c("FALSE", "TRUE")) {
    withr::local_envvar(ROBMA_DEBUG_DLL = debug)
    run <- processx::run(executable,
                         c("--vanilla", file.path(root, "tools/profile-iwmde-qcmde.R"), "--list"),
                         wd = invocation, error_on_status = FALSE, timeout = 60000)
    cat("\nSynthetic profiler loader exit ", run$status, "\n", run$stdout, run$stderr, sep = "")
    expect_equal(run$status, 0L)
    expect_match(run$stdout, normalizePath(root, winslash = "/"), fixed = TRUE)
    expect_match(run$stdout, if (debug == "TRUE") "debug opt-out" else "release", fixed = TRUE)
    expect_match(run$stdout, "native_threads: 1", fixed = TRUE)
    expect_match(run$stdout, "Omitted profiler cases", fixed = TRUE)
  }
})
