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
  for (name in names(fits)) saveRDS(fits[[name]], file.path(cache, paste0(name, ".rds")))
  return(root)
}


.scenario_benchmark_run <- function(root, args, target = root, installed = TRUE) {

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


.scenario_admission_cases <- function() {

  loop <- c('fits <- list()', 'for (key in keys) {',
            '  fits[[key]] <- scenario_fit(key, stop("must never fit"))', '}')
  chosen <- 'scenario_text("chosen", print(fits))'
  cached <- list(old = list(value = "OLD_CACHE"), new = list(value = "NEW_CACHE"))
  case <- function(block, fits = list(), args = "--only=^chosen$", refuse = TRUE,
                   top = character(), text = "[1] 7") {
    list(block = block, fits = fits, args = args, refuse = refuse, top = top, text = text)
  }
  dormant <- c('keys <- "old"', loop, 'keys <- "new"',
               'scenario_text("chosen", print(7))', 'scenario_text("dependent", print(fits))')
  lazy <- c('fit_other <- scenario_fit("irrelevant", stop("must never fit"))',
            'keys <- if (fit_other$value) "old" else "old"', loop, 'keys <- "new"',
            'scenario_text("chosen", print(7))', 'scenario_text("dependent", print(fits))')
  cases <- list(
    rebind_both = case(c('keys <- "old"', loop, 'keys <- "new"', chosen), cached),
    rebind_only_old = case(c('keys <- "old"', loop, 'keys <- "new"', chosen), cached["old"]),
    stable = case(c('keys <- "old"', loop, chosen), cached, refuse = FALSE, text = "OLD_CACHE"),
    names_rebind = case(c('keys <- "old"', 'labels <- "OLD_LABEL"', loop,
                         'names(fits) <- labels', 'labels <- "NEW_LABEL"', chosen), cached),
    unselected_list = case(dormant, args = "--list", refuse = FALSE),
    unselected_filter = case(dormant, refuse = FALSE),
    lazy_dependency_list = case(lazy, args = "--list", refuse = FALSE),
    lazy_dependency_filter = case(lazy, refuse = FALSE),
    derived_selected = case(c('keys <- "old"', loop, 'keys <- "new"',
                             'view <- fits[[1L]]$value', 'scenario_text("chosen", print(view))'), cached),
    derived_unused = case(c('keys <- "old"', loop, 'keys <- "new"', 'view <- fits[[1L]]$value',
                           'scenario_text("chosen", print(7))',
                           'scenario_text("dependent", print(view))'), refuse = FALSE),
    transitive_selector = case(c('selector <- "old"', 'keys <- selector', loop,
                                'selector <- "new"', chosen), cached),
    intervening_replacement = case(c('keys <- "old"', 'keys[1L] <- "new"', loop, chosen), cached),
    intervening_for = case(c('keys <- "old"', 'for (keys in "new") NULL', loop, chosen), cached),
    latest_plain_definition = case(c('keys <- "old"', 'keys[1L] <- "new"',
                                    'keys <- "old"', loop, chosen), cached,
                                   refuse = FALSE, text = "OLD_CACHE"),
    lazy_selector_selected = case(c('fit_other <- scenario_fit("irrelevant", stop("must never fit"))',
                                   'selector <- if (fit_other$value) "old" else "old"',
                                   'keys <- selector', loop, 'selector <- "new"', chosen)),
    unresolved_payload = case(c('keys <- "old"', loop, 'keys <- "new"',
                               'artifact_name <- "chosen"', 'scenario_text(artifact_name, print(fits))'), cached),
    name_root_list = case(c('keys <- "old"', loop, 'keys <- "new"',
                           'artifact_name <- paste0("artifact_", fits[[1L]]$value)',
                           'scenario_text(artifact_name, print(7))'), cached, args = "--list"),
    dynamic_selected = case(c('keys <- "old"', loop, 'keys <- "new"',
                             'scenario_text("chosen", print(get("fits")))'), cached),
    qualified_ordinary_payload = case(c('keys <- "old"', loop,
                                       'scenario_text("chosen", print(base::identity(fits)))'),
                                      cached, refuse = FALSE, text = "OLD_CACHE"),
    qualified_reflective = case(c('keys <- "old"', loop, 'keys <- "new"',
                                 'scenario_text("chosen", print(base::get("fits")))'), cached),
    local_callable = case(c('keys <- "old"', loop, 'keys <- "new"',
                           'view <- function() fits', 'scenario_text("chosen", print(view()))'), cached),
    empty_sequence = case(c('keys <- c()', loop, chosen)),
    unknown_sequence = case(c('keys <- character()', loop, chosen)),
    helper_shadow = case(c('keys <- "old"', 'paste0 <- function(...) "new"',
                           'fits <- list()', 'for (key in keys) {',
                           '  fits[[key]] <- scenario_fit(paste0(key), stop("must never fit"))', '}', chosen), cached),
    top_helper_shadow = case(c('keys <- "old"', 'fits <- list()', 'for (key in keys) {',
                               '  fits[[key]] <- scenario_fit(paste0(key), stop("must never fit"))', '}', chosen),
                             cached, top = 'paste0 <- function(...) "new"'),
    counter_only = case(c('keys <- "old"', 'preparation_count <- 0L', loop,
                         'names(fits) <- { preparation_count <- preparation_count + 1L; keys }',
                         'scenario_text("chosen", print(preparation_count))'), cached),
    counter_intermediate = case(c('keys <- "old"', 'preparation_count <- 0L', loop,
                                 'names(fits) <- { preparation_count <- preparation_count + 1L; keys }',
                                 'observed <- preparation_count',
                                 'scenario_text("chosen", print(list(fits, observed)))'), cached),
    invalid_counter_only = case(c('keys <- "old"', 'preparation_count <- 0L', loop,
                                 'names(fits) <- { preparation_count <- preparation_count + missing_increment; keys }',
                                 'scenario_text("chosen", print(preparation_count))'), cached),
    invalid_counter_dormant = case(c('keys <- "old"', 'preparation_count <- 0L', loop,
                                    'names(fits) <- { preparation_count <- preparation_count + missing_increment; keys }',
                                    'scenario_text("chosen", print(7))',
                                    'scenario_text("dependent", print(fits))'), refuse = FALSE),
    counter_indirect_name = case(c('keys <- "old"', 'preparation_count <- 0L', loop,
                                  'names(fits) <- { preparation_count <- preparation_count + 1L; keys }',
                                  'artifact_name <- paste0("artifact_", preparation_count, length(fits))',
                                  'scenario_text(artifact_name, print(7))'), cached, args = "--only=."),
    counter_direct_name = case(c('keys <- "old"', 'preparation_count <- 0L', loop,
                                'names(fits) <- { preparation_count <- preparation_count + 1L; keys }',
                                'scenario_text(paste0("artifact_", preparation_count, length(fits)), print(7))'),
                               cached, args = "--only=."),
    structural = case(c('fit_needed <- scenario_fit("missing", stop("must never fit"))',
                        'scenario_text("chosen", print(fit_needed))', 'fits <- list()',
                        'for (key in c("old")) {', '  message("unsupported side effect")',
                        '  fits[[key]] <- scenario_fit(key, stop("must never fit"))', '}'))
  )
  mutations <- c(arrow = '"new" -> keys', equal = 'keys = "new"',
                 replacement = 'keys[1L] <- "new"', nested_replacement = 'names(keys)[1L] <- "new"',
                 loop_target = 'for (keys in "new") NULL', nonlocal = 'keys <<- "new"',
                 reflective = 'assign("keys", "new")')
  for (name in names(mutations)) {
    cases[[paste0("mutation_", name)]] <- case(c('keys <- "old"', loop, mutations[[name]], chosen), cached)
  }
  for (shape in c("direct", "compound")) {
    timing <- 'scenario_time("excluded", print(fits))'
    if (shape == "compound") timing <- paste0("if (TRUE) ", timing)
    for (mode in c("filter", "list")) {
      cases[[paste0("time_", shape, "_", mode)]] <- case(
        c('keys <- "old"', loop, 'keys <- "new"', timing, 'scenario_text("chosen", print(7))'),
        args = if (mode == "list") "--list" else "--only=^chosen$",
        refuse = shape == "compound" && mode == "filter")
    }
  }
  cases
}


test_that("source admission refuses demanded binding changes and keeps dormant payloads lazy", {

  for (name in names(.scenario_admission_cases())) {
    fixture <- .scenario_admission_cases()[[name]]
    root <- .scenario_benchmark_fixture(fixture$block, fixture$fits)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    if (length(fixture$top)) {
      path <- file.path(root, "tests/scenarios/test-unit.R")
      writeLines(c(fixture$top, readLines(path)), path)
    }
    # Observe actual readRDS calls after package loading, without replacing the loader.
    bootstrap <- file.path(root, "tools/bench-bootstrap.R")
    lines <- readLines(bootstrap)
    marker <- gsub("\\\\", "/", file.path(root, "cache-reads.log"))
    trace <- sprintf('trace("readRDS", where = baseenv(), print = FALSE, tracer = quote(if (is.character(file) && grepl("/cache/unit/", file, fixed = TRUE)) cat(file, "\\n", file = "%s", append = TRUE)))', marker)
    writeLines(append(lines, trace, after = length(lines) - 1L), bootstrap)
    cache_paths <- list.files(file.path(root, "cache"), full.names = TRUE, recursive = TRUE)
    before <- tools::md5sum(cache_paths)
    run <- .scenario_benchmark_run(root, fixture$args)
    on.exit(unlink(run$out, recursive = TRUE), add = TRUE)
    expect_identical(tools::md5sum(cache_paths), before, info = name)
    if (fixture$refuse) {
      expect_equal(run$status, 1L, info = name)
      expect_match(run$stderr, "Unsupported (binding-sensitive|compound) scenario_fit", info = name)
      expect_false(grepl("No cached fit|unsupported side effect", run$stderr), info = name)
      expect_equal(length(list.files(run$out, recursive = TRUE)), 0L, info = name)
      expect_false(file.exists(marker), info = name)
    } else {
      expect_equal(run$status, 0L, info = name)
      if (fixture$args == "--list") {
        expect_true(file.exists(file.path(run$out, "listing-unit.tsv")), info = name)
      } else {
        path <- file.path(run$out, "unit--chosen.txt")
        expect_true(file.exists(path), info = name)
        if (file.exists(path)) expect_match(paste(readLines(path), collapse = "\n"), fixture$text,
                                           fixed = TRUE, info = name)
      }
      if (!name %in% c("stable", "latest_plain_definition", "qualified_ordinary_payload")) {
        expect_false(file.exists(marker), info = name)
      }
    }
  }
})


test_that("all scenario bodies are admitted before earlier user code or cache reads", {

  root <- .scenario_benchmark_fixture(c(
    'keys <- "old"', 'fits <- list()',
    'for (key in keys) fits[[key]] <- scenario_fit(key, stop("must never fit"))',
    'keys <- "new"', 'scenario_text("chosen", print(fits))'
  ), fits = list(old = list(value = "OLD_CACHE")))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "tests/scenarios/test-unit.R")
  marker <- gsub("\\\\", "/", file.path(root, "earlier-user-code"))
  writeLines(c('testthat::test_that("earlier effect", {',
               sprintf('  writeLines("must not execute", "%s")', marker),
               '  early <- scenario_fit("absent", stop("must never fit"))',
               '  scenario_text("chosen", print(early))', '})', readLines(path)), path)
  before <- tools::md5sum(list.files(file.path(root, "cache"), recursive = TRUE, full.names = TRUE))
  run <- .scenario_benchmark_run(root, "--only=^chosen$")
  on.exit(unlink(run$out, recursive = TRUE), add = TRUE)
  expect_equal(run$status, 1L)
  expect_match(run$stderr, "Unsupported binding-sensitive scenario_fit")
  expect_false(grepl("No cached fit", run$stderr))
  expect_false(file.exists(marker))
  expect_length(list.files(run$out, recursive = TRUE), 0L)
  expect_identical(tools::md5sum(names(before)), before)
})


.scenario_marginal_fixture <- function() {

  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  source <- file.path(package_root, "tests/scenarios/test-bcg.R")
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
