context("Scenario benchmark preparation")


.scenario_benchmark_fixture <- function(block, fits = list(needed = list(value = 7))) {

  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  root <- tempfile("robma-bench-fixture-")
  dir.create(file.path(root, "tools"), recursive = TRUE)
  dir.create(file.path(root, "tests", "scenarios"), recursive = TRUE)
  dir.create(file.path(root, "tests", "testthat"), recursive = TRUE)
  for (file in c("tools/bench-scenario.R", "tools/optimized-dll.R",
                  "tests/scenarios/helper-scenarios.R",
                  "tests/testthat/helper-visual-writer.R")) {
    stopifnot(file.copy(file.path(package_root, file), file.path(root, file)))
  }
  writeLines(c(
    'library(RoBMA)',
    'RoBMA.options(max_cores = 2L, native_threads = 1L)',
    'for (package in c("RoBMA", "BayesTools")) {',
    '  path <- getNamespaceInfo(asNamespace(package), "path")',
    '  dll <- getLoadedDLLs()[[package]][["path"]]',
    '  stopifnot(identical(tolower(normalizePath(dirname(path))), tolower(normalizePath(.libPaths()[[1L]]))))',
    '  cat("CHILD_PRODUCT", package, path, dll, unname(tools::md5sum(dll)), "\\n")',
    '}',
    'sys.source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])), "bench-scenario.R"), envir = globalenv())'
  ), file.path(root, "tools", "bench-bootstrap.R"))
  writeLines(c('scenario_start("unit")', 'testthat::test_that("synthetic workload", {',
               block, '})'), file.path(root, "tests/scenarios/test-unit.R"))
  cache <- file.path(root, "cache", "unit")
  dir.create(cache, recursive = TRUE)
  for (name in names(fits)) saveRDS(fits[[name]], file.path(cache, paste0(name, ".rds")))
  return(root)
}


.scenario_benchmark_run <- function(root, args, target = root, installed = TRUE) {

  out <- tempfile("robma-bench-output-")
  executable <- file.path(R.home("bin"),
                          if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  withr::local_envvar(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  result <- processx::run(executable, c("--vanilla", file.path(root, "tools/bench-bootstrap.R"),
                          "unit", if (installed) "--installed", paste0("--root=", target),
                          paste0("--cache=", file.path(root, "cache")),
                          paste0("--out=", out), "--threads=1", args),
                          error_on_status = FALSE, timeout = 60000)
  cat("\nSynthetic benchmark exit ", result[["status"]], "\n",
      result[["stdout"]], result[["stderr"]], sep = "")
  result[["out"]] <- out
  return(result)
}


test_that("listing and filtering defer an irrelevant fit-list loop", {

  root <- .scenario_benchmark_fixture(c(
    'fit_needed <- scenario_fit("needed", stop("must never fit"))',
    'fit_unneeded <- list()',
    'for (key in c("missing")) {',
    '  fit_unneeded[[key]] <- scenario_fit(key, stop("must never fit"))',
    '}',
    'scenario_text("chosen", print(fit_needed$value))',
    'scenario_text("dependent", print(fit_unneeded[[1L]]$value))'
  ))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  listed <- .scenario_benchmark_run(root, "--list")
  on.exit(unlink(listed[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(listed[["status"]], 0L)
  listing <- file.path(listed[["out"]], "listing-unit.tsv")
  if (file.exists(listing)) expect_identical(read.delim(listing)$name, c("chosen", "dependent"))
  chosen <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(chosen[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(chosen[["status"]], 0L)
  text <- file.path(chosen[["out"]], "unit--chosen.txt")
  expect_true(file.exists(text))
  if (file.exists(text)) expect_identical(readLines(text), "[1] 7")
  dependent <- .scenario_benchmark_run(root, "--only=^dependent$")
  on.exit(unlink(dependent[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(dependent[["status"]], 1L)
  expect_match(dependent[["stderr"]], "No cached fit 'missing'", fixed = TRUE)
})


test_that("benchmark workload and helpers belong to the driver tree", {

  root <- .scenario_benchmark_fixture('scenario_text("chosen", print("MAIN WORKLOAD"))')
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  target <- tempfile("robma-helperless-target-")
  dir.create(file.path(target, "R"), recursive = TRUE)
  on.exit(unlink(target, recursive = TRUE), add = TRUE)
  writeLines(c("Package: RoBMA", "Version: 0.0.0.9000",
               "Title: Synthetic Benchmark Loader Fixture",
               "Description: A minimal API-compatible loader fixture without native code.",
               "Author: RoBMA tests", "Maintainer: RoBMA tests <tests@example.org>",
               "License: MIT"), file.path(target, "DESCRIPTION"))
  writeLines("export(RoBMA.options,RoBMA.get_option)", file.path(target, "NAMESPACE"))
  writeLines(c('.state <- new.env(parent = emptyenv())', '.state$threads <- 1L',
               'RoBMA.options <- function(native_threads) .state$threads <- native_threads',
               'RoBMA.get_option <- function(name) .state$threads'), file.path(target, "R", "options.R"))
  dir.create(file.path(target, "tests", "scenarios"), recursive = TRUE)
  writeLines('stop("TARGET WORKLOAD MUST NOT RUN")',
             file.path(target, "tests/scenarios/test-unit.R"))
  run <- .scenario_benchmark_run(root, "--only=^chosen$", target = target,
                                 installed = FALSE)
  on.exit(unlink(run[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(run[["status"]], 0L)
  text <- file.path(run[["out"]], "unit--chosen.txt")
  expect_true(file.exists(text))
  if (file.exists(text)) expect_identical(readLines(text), '[1] "MAIN WORKLOAD"')
  expect_false(dir.exists(file.path(target, "tools")))

  # A real incompatible target API must fail clearly, before artifact work.
  writeLines(c('RoBMA.options <- function(native_threads) stop("TARGET API INCOMPATIBLE")',
               'RoBMA.get_option <- function(name) 1L'),
             file.path(target, "R", "options.R"))
  incompatible <- .scenario_benchmark_run(root, "--only=^chosen$", target = target,
                                          installed = FALSE)
  on.exit(unlink(incompatible[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(incompatible[["status"]], 1L)
  expect_match(incompatible[["stderr"]], "TARGET API INCOMPATIBLE", fixed = TRUE)
})


test_that("benchmark plots use the canonical scenario SVG renderer", {

  skip_if_not_installed("vdiffr")
  root <- .scenario_benchmark_fixture('scenario_plot("chosen", graphics::plot(1:3))')
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  run <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(run[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(run[["status"]], 0L)
  env <- new.env(parent = globalenv())
  source(file.path(root, "tests/testthat/helper-visual-writer.R"), local = env)
  expected <- file.path(run[["out"]], "canonical-control.svg")
  env$.write_canonical_svg(function() graphics::plot(1:3), expected, "chosen")
  actual <- file.path(run[["out"]], "unit--chosen.svg")
  expect_identical(unname(tools::md5sum(actual)), unname(tools::md5sum(expected)))
})


test_that("benchmark memory reads MB maxima in ordinary and finite-heap layouts", {

  root <- .scenario_benchmark_fixture('scenario_text("chosen", print(7))')
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "tests/scenarios/helper-scenarios.R"), local = env)
  find_memory <- function(expr) {
    if (!is.call(expr)) return(NULL)
    if (identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("memory_gb"))) return(expr[[3L]])
    for (part in as.list(expr)) {
      found <- find_memory(part)
      if (!is.null(found)) return(found)
    }
    return(NULL)
  }
  source <- parse(file.path(root, "tools/bench-scenario.R"), keep.source = FALSE)
  memory_expr <- find_memory(as.call(c(list(as.name("{")), as.list(source))))
  ordinary <- matrix(c(10, 100, 1, 10, 20, 200, 2, 20, 40, 400, 4, 40), nrow = 2L,
                     dimnames = list(c("Ncells", "Vcells"),
                                     c("used", "(Mb)", "gc trigger", "(Mb)", "max used", "(Mb)")))
  finite <- cbind(ordinary[, 1:4], `limit (Mb)` = c(NA, 4096), ordinary[, 5:6])
  for (layout in list(ordinary, finite)) {
    env$gc <- function(...) layout
    evaluation <- list2env(list(gc = env$gc, bench_env = env), parent = globalenv())
    reference <- sum(layout[, match("max used", colnames(layout)) + 1L]) / 1024
    expect_equal(eval(memory_expr, evaluation), reference)
  }
  run <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(run[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(run[["status"]], 0L)
  result <- read.delim(file.path(run[["out"]], "results.tsv"))
  # This tiny real driver call cannot have thousands of GB of R heap maxima.
  expect_true(is.finite(result$memory_gb) && result$memory_gb >= 0 && result$memory_gb < 1)
})


test_that("Bem-shaped two-list preparation retains real named outputs and loop bindings", {

  root <- .scenario_benchmark_fixture(c(
    'mu0_seq <- c(-0.25, 0, 0.25)',
    'preparation_count <- 0L',
    'fit_3PSM_null_list <- list()',
    'fit_PET_null_list <- list()',
    'for (i in seq_along(mu0_seq)) {',
    '  fit_3PSM_null_list[[i]] <- scenario_fit(paste0("first_", mu0_seq[i]), stop("must never fit"))',
    '  fit_PET_null_list[[i]] <- scenario_fit(paste0("second_", mu0_seq[i]), stop("must never fit"))',
    '}',
    'names(fit_3PSM_null_list) <- { preparation_count <- preparation_count + 1L; paste0("mu=", mu0_seq) }',
    'names(fit_PET_null_list) <- { preparation_count <- preparation_count + 1L; paste0("mu=", mu0_seq) }',
    'scenario_text("both", {',
    '  print(list(first = fit_3PSM_null_list, second = fit_PET_null_list, i = i, preparations = preparation_count))',
    '})',
    'scenario_text("again", print(c(identical(fit_PET_null_list[[3L]]$value, 6), preparation_count == 2L)))'
  ), fits = stats::setNames(lapply(1:6, function(i) list(value = as.numeric(i))),
                            c(paste0("first_", c(-0.25, 0, 0.25)),
                              paste0("second_", c(-0.25, 0, 0.25)))))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  run <- .scenario_benchmark_run(root, "--only=both|again")
  on.exit(unlink(run[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(run[["status"]], 0L)
  path <- file.path(run[["out"]], "unit--both.txt")
  expect_true(file.exists(path))
  if (file.exists(path)) {
    expected <- capture.output(print(list(
      first = stats::setNames(lapply(1:3, function(i) list(value = as.numeric(i))),
                              paste0("mu=", c(-0.25, 0, 0.25))),
      second = stats::setNames(lapply(4:6, function(i) list(value = as.numeric(i))),
                               paste0("mu=", c(-0.25, 0, 0.25))), i = 3L,
      preparations = 2L)))
    expect_identical(readLines(path), expected)
  }
  again <- file.path(run[["out"]], "unit--again.txt")
  if (file.exists(again)) expect_identical(readLines(again), "[1] TRUE TRUE")
})


test_that("unsupported compound fit loading refuses before reading caches", {

  root <- .scenario_benchmark_fixture(c(
    'fit_needed <- scenario_fit("missing", stop("must never fit"))',
    'scenario_text("chosen", print(fit_needed$value))',
    'fits <- list()',
    'for (key in c("missing")) {',
    '  message("unsupported side effect")',
    '  fits[[key]] <- scenario_fit(key, stop("must never fit"))',
    '}'
  ))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  run <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(run[["out"]], recursive = TRUE), add = TRUE)
  expect_equal(run[["status"]], 1L)
  expect_match(run[["stderr"]], "Unsupported compound.*scenario_fit")
  expect_false(grepl("No cached fit|unsupported side effect", run[["stderr"]]))
})
