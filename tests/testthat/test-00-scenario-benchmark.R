context("Scenario benchmark preparation")


.scenario_benchmark_fixture <- function(block, fits = list(needed = list(value = 7))) {

  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  skip_if_not(file.exists(file.path(package_root, "tools/bench-scenario.R")),
              "Developer benchmark requires a source checkout.")
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
  for (name in names(fits)) {
    saveRDS(fits[[name]], file.path(cache, paste0(name, ".rds")))
    writeLines('stop("must never fit")', file.path(cache, paste0(name, ".call.txt")))
  }
  return(root)
}


.scenario_benchmark_run <- function(root, args, target = root, installed = TRUE, threads = 1L) {

  # Benchmark rows include HEAD: give disposable roots an identifiable HEAD.
  if (!dir.exists(file.path(target, ".git"))) {
    stopifnot(system2("git", c("-C", shQuote(target), "init", "--quiet")) == 0L)
    stopifnot(system2("git", c("-C", shQuote(target), "-c", "user.name=RoBMA-tests",
      "-c", "user.email=tests@example.org", "commit", "--quiet", "--allow-empty",
      "-m", "fixture")) == 0L)
  }
  out <- tempfile("robma-bench-output-")
  executable <- file.path(R.home("bin"),
                          if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  withr::local_envvar(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  result <- processx::run(executable, c("--vanilla", file.path(root, "tools/bench-bootstrap.R"),
                          "unit", if (installed) "--installed", paste0("--root=", target),
                          paste0("--cache=", file.path(root, "cache")),
                          paste0("--out=", out), if (!is.null(threads)) paste0("--threads=", threads), args),
                          error_on_status = FALSE, timeout = 60000)
  cat("\nSynthetic benchmark exit ", result[["status"]], "\n",
      result[["stdout"]], result[["stderr"]], sep = "")
  result[["out"]] <- out
  return(result)
}


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


test_that("Bem-shaped fit-list construction keeps its named outputs and loop binding", {

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


.scenario_marginal_fixture <- function(source = NULL) {

  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  if (is.null(source)) source <- file.path(package_root, "tests/scenarios/test-bcg.R")
  skip_if_not(file.exists(source), "Canonical BCG scenario requires a source checkout.")
  preparation <- NULL
  head <- function(x) if (is.call(x)) paste(deparse(x[[1L]]), collapse = "") else ""
  for (expr in parse(source, keep.source = FALSE)) {
    if (!grepl("test_that$", head(expr))) next
    statements <- as.list(expr[[3L]])[-1L]
    first <- which(vapply(statements, function(st) head(st) == "<-" &&
      is.symbol(st[[2L]]) && as.character(st[[2L]]) %in% c("marginal_draws", "marginal_inputs"), logical(1)))
    last <- which(vapply(statements, function(st) head(st) == "scenario_text" &&
      identical(st[[2L]], "fit_reg3_marginal_means_comparison"), logical(1)))
    if (length(first) && length(last)) preparation <- statements[seq.int(first[[1L]], last[[1L]])]
  }
  stopifnot(!is.null(preparation))
  draws <- list(a = c(1, 3), b = c(2, 4))
  attr(draws, "data") <- data.frame(alloc = factor(c("a", "b")), ablat = c(0, 1))
  primary <- list(inference = list(averaged = list(mu_alloc__xXx__ablat = draws)))
  reference <- structure(list(intercept = 0.5, beta = c(0.25, 0.1, 0.2)),
                         class = "scenario_benchmark_reference")
  root <- .scenario_benchmark_fixture(c(
    'dat <- data.frame(ablat = c(10, 20))',
    'fit_reg3 <- scenario_fit("primary", stop("must never fit"))',
    'fit_reg3_metafor <- scenario_fit("reference", stop("must never fit"))',
    'fit_reg3_emm <- scenario_time("fit_reg3_emm", fit_reg3)',
    unlist(lapply(preparation, deparse), use.names = FALSE)
  ), fits = list(primary = primary, reference = reference))
  # Explicit synthetic cached reference: this is not a fitted metafor object.
  scenario <- file.path(root, "tests/scenarios/test-unit.R")
  writeLines(c('predict.scenario_benchmark_reference <- function(object, newmods, ...) {',
               '  list(pred = unname(object$intercept + as.vector(newmods %*% object$beta)))',
               '}', readLines(scenario)), scenario)
  bootstrap <- file.path(root, "tools/bench-bootstrap.R")
  lines <- readLines(bootstrap)
  marker <- gsub("\\\\", "/", file.path(root, "cache-reads.log"))
  trace <- sprintf('trace("readRDS", where = baseenv(), print = FALSE, tracer = quote(if (is.character(file) && grepl("/cache/unit/", file, fixed = TRUE)) cat(file, "\\n", file = "%s", append = TRUE)))', marker)
  writeLines(append(lines, trace, after = length(lines) - 1L), bootstrap)
  attr(root, "preparation") <- preparation
  root
}


test_that("BCG marginal preparation agrees for full and selected synthetic cached inputs", {

  root <- .scenario_marginal_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  full <- .scenario_benchmark_run(root, "--only=.")
  on.exit(unlink(full$out, recursive = TRUE), add = TRUE)
  selected <- .scenario_benchmark_run(root, "--only=^fit_reg3_marginal_means_comparison$")
  on.exit(unlink(selected$out, recursive = TRUE), add = TRUE)
  expect_equal(full$status, 0L)
  expect_equal(selected$status, 0L)
  name <- "unit--fit_reg3_marginal_means_comparison.txt"
  expect_true(file.exists(file.path(full$out, name)))
  expect_true(file.exists(file.path(selected$out, name)))
  for (run in list(full, selected)) {
    if (file.exists(file.path(run$out, "results.tsv"))) {
      expect_true(all(read.delim(file.path(run$out, "results.tsv"))$status == "ok"))
    }
  }
  if (file.exists(file.path(full$out, name)) && file.exists(file.path(selected$out, name))) {
    expect_identical(readLines(file.path(full$out, name)), readLines(file.path(selected$out, name)))
  }
  marker <- file.path(root, "cache-reads.log")
  unlink(marker)
  listed <- .scenario_benchmark_run(root, "--list")
  on.exit(unlink(listed$out, recursive = TRUE), add = TRUE)
  expect_equal(listed$status, 0L)
  expect_false(file.exists(marker))
  expect_false(file.exists(file.path(listed$out, name)))
  listing <- file.path(listed$out, "listing-unit.tsv")
  expect_true(file.exists(listing))
  if (file.exists(listing)) expect_identical(read.delim(listing)$name,
    c("fit_reg3_emm", "fit_reg3_marginal_means_comparison"))
})


test_that("benchmark cache versions match literals without forcing expressions", {

  block <- c(
    'fit_needed <- scenario_fit("needed", stop("must never fit"), cache_version = 2L)',
    'scenario_text("chosen", print(fit_needed$value))'
  )
  root <- .scenario_benchmark_fixture(block)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  call_path <- file.path(root, "cache", "unit", "needed.call.txt")
  writeLines(c("cache_version: 2", 'stop("must never fit")'), call_path)
  matched <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(matched$out, recursive = TRUE), add = TRUE)
  expect_equal(matched$status, 0L)
  text <- file.path(matched$out, "unit--chosen.txt")
  expect_true(file.exists(text))
  if (file.exists(text)) expect_identical(readLines(text), "[1] 7")

  writeLines(c("cache_version: 1", 'stop("must never fit")'), call_path)
  stale <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(stale$out, recursive = TRUE), add = TRUE)
  expect_equal(stale$status, 1L)
  expect_match(stale$stderr, "Missing or changed fitting call", fixed = TRUE)

  # Count an actual version-expression side effect outside the failed child.
  # Both the old ignore behavior and the new refusal leave the count at zero;
  # the intended new behavior is refusal rather than silent version omission.
  counter <- file.path(root, "version-counter.txt")
  writeLines("0", counter)
  counter_literal <- paste(deparse(normalizePath(counter, winslash = "/")), collapse = "")
  opaque <- paste0('fit_needed <- scenario_fit("needed", stop("must never fit"), ',
    'cache_version = { writeLines(as.character(as.integer(readLines(', counter_literal,
    ')) + 1L), ', counter_literal, '); 2L })')
  writeLines(c('scenario_start("unit")', 'testthat::test_that("synthetic workload", {',
    opaque, block[[2L]], '})'), file.path(root, "tests", "scenarios", "test-unit.R"))
  refused <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(refused$out, recursive = TRUE), add = TRUE)
  expect_equal(refused$status, 1L)
  expect_match(refused$stderr, "expressions are not evaluated", fixed = TRUE)
  expect_identical(readLines(counter), "0")
})


test_that("benchmark cache loading validates genuine persisted fit contracts", {

  setup <- review_compile_object(brma(yi = c(.1, .2, -.1), sei = rep(.2, 3L),
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = prior("normal", list(0, 1)),
    prior_heterogeneity = prior("point", list(.1)), only_priors = TRUE, silent = TRUE))
  object <- review_attach_replay(setup,
    cbind(mu = seq(-.2, .3, length.out = 24L), tau = .1))
  required <- sub("_version$", "", setdiff(BayesTools::JAGS_fit_contract_schema()$field, "schema_version"))
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = required)
  root <- .scenario_benchmark_fixture(c(
    'fit_needed <- scenario_fit("needed", stop("must never fit"))',
    'scenario_text("chosen", print(nobs(fit_needed)))'
  ), fits = list(needed = object))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  helper <- new.env(parent = globalenv())
  source(file.path(root, "tests/scenarios/helper-scenarios.R"), local = helper)
  helper$scenario_start("unit", root = file.path(root, "tests/scenarios"), show_output = FALSE)
  expect_identical(helper$.scenario_check_cached_fit("needed", object), object)
  accepted <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(accepted$out, recursive = TRUE), add = TRUE)
  expect_equal(accepted$status, 0L)
  output <- file.path(accepted$out, "unit--chosen.txt")
  expect_true(file.exists(output))
  if (file.exists(output)) expect_identical(readLines(output), "[1] 3")
  stale <- object
  attr(stale$fit, "fit_contract") <- NULL
  # This is the actual finalized replay with its persisted contract removed,
  # not a synthetic accepted class or forged schema field.
  refusal <- tryCatch(helper$.scenario_check_cached_fit("needed", stale), error = identity)
  expect_s3_class(refusal, "RoBMA_scenario_stale_cache")
  expect_identical(refusal$fit, "needed")
  expect_identical(refusal$scenario, "unit")
  saveRDS(stale, file.path(root, "cache/unit/needed.rds"))
  rejected <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(rejected$out, recursive = TRUE), add = TRUE)
  expect_equal(rejected$status, 1L)
  expect_match(rejected$stderr, "Cached scenario fit 'needed' of scenario 'unit' is stale.", fixed = TRUE)
})


test_that("benchmark option defaults and missing call metadata remain explicit", {

  block <- function(version = "") c(
    paste0('fit_needed <- scenario_fit("needed", stop("must never fit")', version, ')'),
    'scenario_text("chosen", dput(list(value = fit_needed$value, threads = RoBMA.get_option("native_threads"))))')
  root <- .scenario_benchmark_fixture(block())
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  bootstrap <- file.path(root, "tools/bench-bootstrap.R")
  writeLines(sub("native_threads = 1L", "native_threads = 2L", readLines(bootstrap), fixed = TRUE), bootstrap)
  scenario <- file.path(root, "tests/scenarios/test-unit.R")
  write_block <- function(code) writeLines(c('scenario_start("unit")',
    'testthat::test_that("synthetic workload", {', code, '})'), scenario)
  for (version in c("", ", cache_version = NULL")) {
    write_block(block(version))
    for (threads in list(NULL, 1L)) {
      result <- .scenario_benchmark_run(root, "--only=^chosen$", threads = threads)
      on.exit(unlink(result$out, recursive = TRUE), add = TRUE)
      expect_equal(result$status, 0L)
      output <- file.path(result$out, "unit--chosen.txt")
      expect_true(file.exists(output))
      if (file.exists(output)) expect_identical(readLines(output), "list(value = 7, threads = 1L)")
    }
  }
  call_path <- file.path(root, "cache/unit/needed.call.txt")
  stopifnot(file.remove(call_path))
  missing <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(missing$out, recursive = TRUE), add = TRUE)
  expect_equal(missing$status, 1L)
  expect_match(missing$stderr, "Missing or changed fitting call", fixed = TRUE)
  writeLines('stop("must never fit")', call_path)
  for (version in c(", cache_version = 0", ", cache_version = 1.5")) {
    write_block(block(version))
    invalid <- .scenario_benchmark_run(root, "--only=^chosen$")
    on.exit(unlink(invalid$out, recursive = TRUE), add = TRUE)
    expect_equal(invalid$status, 1L)
    expect_match(invalid$stderr, "'cache_version' must be NULL or one positive integer.", fixed = TRUE)
  }
  write_block(block())
  for (threads in c(0, -1, 1.5, Inf)) {
    invalid <- .scenario_benchmark_run(root, "--only=^chosen$", threads = threads)
    on.exit(unlink(invalid$out, recursive = TRUE), add = TRUE)
    expect_equal(invalid$status, 1L)
    expect_match(invalid$stderr, "--threads must be a positive integer.", fixed = TRUE)
  }
})
