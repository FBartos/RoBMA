#!/usr/bin/env Rscript
# Benchmark selected scenario artifacts (scenario_plot / scenario_text /
# scenario_time) from a tests/scenarios/test-<scenario>.R file against the
# cached fits, without running the rest of the scenario.
#
#   Rscript tools/bench-scenario.R <scenario> --out=<dir> [--only=<regex>]
#       [--index=<from:to>] [--list] [--profile] [--limit=<seconds>]
#       [--threads=<n>] [--root=<RoBMA tree>] [--cache=<scenario cache root>]
#       [--installed] [--tag=<label>]
#
# Fits and every assignment that depends on a fit are bound lazily, so only what
# a selected artifact needs is loaded or computed; that preparation is timed
# separately (prep_elapsed) and excluded from the artifact's elapsed time.
#
# The run directory collects `results.tsv` (one row per evaluated artifact: wall
# time, peak memory, md5 fingerprint of the printed text or the rendered SVG),
# the artifact files themselves and, with --profile, an `Rprof` file plus a
# `.prof.txt` summary. Every per-artifact file is prefixed `<scenario>--`.
# This runner never writes timing baselines or testthat snapshots.

args <- commandArgs(trailingOnly = TRUE)
flag <- function(name) paste0("--", name) %in% args
opt  <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

cmd      <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[length(file_arg)]])
} else {
  file.path("tools", "bench-scenario.R")
}
# The package this script ships with owns the scenario definitions, the fit
# caches and the committed timings, also when --root points at another tree.
main_root <- normalizePath(file.path(dirname(script_path), ".."),
                           winslash = "/", mustWork = TRUE)

scenario  <- args[!grepl("^--", args)][1L]
if (is.na(scenario)) {
  stop("Specify a scenario name, for example 'bem2011'.", call. = FALSE)
}
only      <- opt("only", ".")
root      <- normalizePath(opt("root", main_root), winslash = "/", mustWork = TRUE)
cache     <- normalizePath(opt("cache", file.path(main_root, "tests", "scenarios", "cache")),
                           winslash = "/", mustWork = TRUE)
out_dir   <- opt("out")
if (is.null(out_dir)) {
  stop("Specify the run directory with --out=<dir>.", call. = FALSE)
}
out_dir   <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
# This runner must never be able to leave a file among the timings, snapshots,
# results or fit caches, whatever the caller passes, so a run directory inside
# either tree's tests/ is refused rather than created.
local({
  comparable <- function(path) {
    if (.Platform$OS.type == "windows") tolower(path) else path
  }
  for (tree in unique(c(main_root, root))) {
    tests_dir <- file.path(tree, "tests")
    if (startsWith(paste0(comparable(out_dir), "/"),
                   paste0(comparable(tests_dir), "/"))) {
      stop("--out must not point inside '", tests_dir,
           "'; use a run directory outside the package tree.", call. = FALSE)
    }
  }
})
limit     <- as.numeric(opt("limit", "Inf"))
threads   <- opt("threads")
tag       <- opt("tag", "")
profile   <- flag("profile")
list_only <- flag("list")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(NOT_CRAN = "true", AGENT = "1")
setwd(root)
if (flag("installed")) {
  suppressPackageStartupMessages(library(RoBMA))
  build <- paste0("installed ", as.character(utils::packageVersion("RoBMA")))
} else {
  # Timings describe the release build, not pkgbuild's default -O0 debug DLL.
  source(file.path(main_root, "tools", "optimized-dll.R"))
  ensure_optimized_dll(root, quiet = TRUE)
  pkgload::load_all(root, quiet = TRUE, compile = FALSE, helpers = FALSE,
                   debug = FALSE)
  build <- "load_all release"
}
if (!is.null(threads)) RoBMA.options(native_threads = as.integer(threads))
sha <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "--short", "HEAD"), stdout = TRUE),
                error = function(e) NA_character_)

bench_env <- new.env(parent = globalenv())
sys.source(file.path(main_root, "tests", "scenarios", "helper-scenarios.R"), envir = bench_env,
           chdir = TRUE)

state <- new.env()
state$lazy    <- character()
state$rows    <- list()
state$listing <- list()

index_range <- opt("index")   # "from:to" positions in listing-<scenario>.tsv (run --list first)
index_names <- NULL
if (!is.null(index_range)) {
  listing_path <- file.path(out_dir, paste0("listing-", scenario, ".tsv"))
  if (!file.exists(listing_path)) {
    stop("No artifact listing at ", listing_path,
         ". Run the same --out directory with --list first.", call. = FALSE)
  }
  listing     <- utils::read.delim(listing_path, stringsAsFactors = FALSE)
  bounds      <- as.integer(strsplit(index_range, ":", fixed = TRUE)[[1L]])
  index_names <- listing$name[seq.int(bounds[[1L]], min(bounds[[2L]], nrow(listing)))]
}
selected <- function(name) {
  grepl(only, name, perl = TRUE) && (is.null(index_names) || name %in% index_names)
}

make_lazy <- function(lhs, expr, env, loader = NULL) {
  force(lhs); force(expr); force(env); force(loader)
  if (is.null(loader)) {
    delayedAssign(lhs, eval(expr, env), assign.env = env)
  } else {
    delayedAssign(lhs, loader(), assign.env = env)
  }
  state$lazy <- union(state$lazy, lhs)
}

load_fit <- function(name) {
  path <- file.path(cache, scenario, paste0(name, ".rds"))
  if (!file.exists(path)) stop("No cached fit '", name, "' at ", path, call. = FALSE)
  readRDS(path)
}

bench_run <- function(type, name, expr, env) {
  if (list_only) {
    state$listing[[length(state$listing) + 1L]] <- data.frame(type = type, name = name)
    return(invisible(NULL))
  }
  if (!selected(name)) {
    if (type == "time") return(eval(expr, env))   # value may be consumed downstream
    return(invisible(NULL))
  }
  # force lazy inputs outside the timed region
  prep_started <- proc.time()[["elapsed"]]
  for (v in intersect(all.vars(expr), state$lazy)) {
    if (exists(v, envir = env, inherits = TRUE)) force(get(v, envir = env, inherits = TRUE))
  }
  prep_elapsed <- proc.time()[["elapsed"]] - prep_started

  file_stub <- file.path(out_dir, paste0(scenario, "--", gsub("[^A-Za-z0-9_.-]", "_", name)))
  prof_file <- paste0(file_stub, ".Rprof")
  status <- "ok"; value <- NULL; fingerprint <- NA_character_
  invisible(gc(reset = TRUE, full = TRUE))
  set.seed(1)
  if (type == "plot") {
    svg_file <- paste0(file_stub, ".svg")
  }
  if (profile) utils::Rprof(prof_file, interval = 0.02)
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    if (is.finite(limit)) setTimeLimit(elapsed = limit, transient = TRUE)
    if (type == "plot") {
      get(".write_canonical_svg", bench_env)(function() {
        set.seed(1)
        get(".scenario_evaluate_plot", bench_env)(expr, env, restore_par = FALSE)
      }, svg_file, name)
    } else if (type == "text") {
      old <- options(width = 150L)
      on.exit(options(old), add = TRUE)
      text <- utils::capture.output({
        visible <- suppressMessages(withVisible(eval(expr, env)))
        if (visible[["visible"]]) print(visible[["value"]])
      }, type = "output")
      writeLines(text, paste0(file_stub, ".txt"))
      visible[["value"]]
    } else {
      eval(expr, env)
    }
  }, error = function(e) {
    status <<- paste0("error: ", gsub("[\r\n\t]+", " ", conditionMessage(e)))
    NULL
  })
  elapsed <- proc.time()[["elapsed"]] - started
  setTimeLimit(elapsed = Inf)
  if (profile) utils::Rprof(NULL)
  if (type == "plot") {
    if (file.exists(svg_file)) fingerprint <- unname(tools::md5sum(svg_file))
  } else if (type == "text" && file.exists(paste0(file_stub, ".txt"))) {
    fingerprint <- unname(tools::md5sum(paste0(file_stub, ".txt")))
  } else if (type == "time" && status == "ok") {
    saveRDS(result, paste0(file_stub, ".rds"))
    txt <- utils::capture.output(print(result))
    writeLines(txt, paste0(file_stub, ".txt"))
    fingerprint <- unname(tools::md5sum(paste0(file_stub, ".txt")))
  }
  memory_gb <- sum(gc()[, 6L]) / 1024
  if (profile && file.exists(prof_file)) {
    summary <- tryCatch(utils::summaryRprof(prof_file), error = function(e) NULL)
    if (!is.null(summary)) {
      sink(paste0(file_stub, ".prof.txt"))
      cat("artifact:", name, " elapsed:", round(elapsed, 2), " status:", status, "\n\n== by.self (top 40)\n")
      print(utils::head(summary[["by.self"]], 40L))
      cat("\n== by.total (top 80)\n")
      print(utils::head(summary[["by.total"]], 80L))
      sink()
    }
  }
  row <- data.frame(
    scenario = scenario, type = type, name = name, status = status,
    elapsed = round(elapsed, 2), prep_elapsed = round(prep_elapsed, 2),
    memory_gb = round(memory_gb, 3), fingerprint = fingerprint,
    threads = as.character(RoBMA.get_option("native_threads")),
    build = build, sha = sha, tag = tag,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE
  )
  path <- file.path(out_dir, "results.tsv")
  utils::write.table(row, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = !file.exists(path), append = file.exists(path))
  message(sprintf("[bench] %-5s %-55s %9.2f s  prep %7.2f s  %5.2f GB  %s",
                  type, name, elapsed, prep_elapsed, memory_gb, status))
  if (type == "time") return(result)
  invisible(result)
}

assign("scenario_start", function(name, ...) invisible(NULL), envir = bench_env)
assign("scenario_fit", function(name, code, cache_version = NULL) load_fit(name), envir = bench_env)
assign("scenario_time", function(name, code) bench_run("time", name, substitute(code), parent.frame()),
       envir = bench_env)
assign("scenario_text", function(name, code) bench_run("text", name, substitute(code), parent.frame()),
       envir = bench_env)
assign("scenario_plot", function(name, code) bench_run("plot", name, substitute(code), parent.frame()),
       envir = bench_env)

call_head <- function(x) if (is.call(x)) paste(deparse(x[[1L]]), collapse = "") else ""
has_scenario_call <- function(x) any(c("scenario_time", "scenario_text", "scenario_plot") %in% all.names(x))

fit_list_blocks <- function(statements) {
  blocks <- list()
  index <- 1L
  initializer <- function(st) {
    call_head(st) %in% c("<-", "=") && is.symbol(st[[2L]]) &&
      call_head(st[[3L]]) == "list" && length(st[[3L]]) == 1L
  }
  while (index <= length(statements)) {
    first <- index
    outputs <- character()
    while (index <= length(statements) && initializer(statements[[index]])) {
      outputs <- c(outputs, as.character(statements[[index]][[2L]]))
      index <- index + 1L
    }
    if (length(outputs) > 0L && index <= length(statements) &&
        call_head(statements[[index]]) == "for" &&
        "scenario_fit" %in% all.names(statements[[index]])) {
      loop <- statements[[index]]
      body <- if (call_head(loop[[4L]]) == "{") as.list(loop[[4L]])[-1L] else list(loop[[4L]])
      valid_assignment <- function(st) {
        if (!call_head(st) %in% c("<-", "=") ||
            !call_head(st[[2L]]) %in% c("[", "[[") ||
            !is.symbol(st[[2L]][[2L]]) ||
            !as.character(st[[2L]][[2L]]) %in% outputs ||
            call_head(st[[3L]]) != "scenario_fit") return(FALSE)
        # The loop itself only constructs the declared lists. Cached fit code
        # stays unevaluated, exactly as for an ordinary scalar fit assignment.
        !"scenario_fit" %in% all.names(st[[2L]])
      }
      if (length(body) > 0L && all(vapply(body, valid_assignment, logical(1)))) {
        index <- index + 1L
        while (index <= length(statements)) {
          st <- statements[[index]]
          names_assignment <- call_head(st) %in% c("<-", "=") &&
            call_head(st[[2L]]) == "names" && length(st[[2L]]) == 2L &&
            is.symbol(st[[2L]][[2L]]) &&
            as.character(st[[2L]][[2L]]) %in% outputs &&
            !"scenario_fit" %in% all.names(st[[3L]])
          if (!names_assignment) break
          index <- index + 1L
        }
        blocks[[as.character(first)]] <- list(
          last = index - 1L,
          outputs = unique(c(outputs, as.character(loop[[2L]]))),
          expr = as.call(c(list(as.name("{")), statements[first:(index - 1L)]))
        )
        next
      }
    }
    # Preflight every statement before any selected artifact can read a fit.
    # This is an explicit supported construction block, not an R interpreter.
    index <- first
    st <- statements[[index]]
    scalar_fit <- call_head(st) %in% c("<-", "=") &&
      is.symbol(st[[2L]]) && call_head(st[[3L]]) == "scenario_fit"
    if ("scenario_fit" %in% all.names(st) && !scalar_fit) {
      stop("Unsupported compound scenario_fit construction; use scalar fit assignments ",
           "or list initializers followed by a fit-only for loop and names assignments.",
           call. = FALSE)
    }
    index <- index + 1L
  }
  blocks
}

make_lazy_fit_lists <- function(block, env) {
  force(block); force(env)
  result <- NULL
  evaluated <- FALSE
  load <- function(output) {
    if (!evaluated) {
      eval(block$expr, env)
      result <<- mget(block$outputs, envir = env, inherits = FALSE)
      evaluated <<- TRUE
    }
    result[[output]]
  }
  for (output in block$outputs) {
    local({ name <- output; make_lazy(name, NULL, env, loader = function() load(name)) })
  }
  invisible(NULL)
}

run_block <- function(body, env) {
  statements <- as.list(body)[-1L]
  blocks <- fit_list_blocks(statements)
  skip <- integer()
  for (index in seq_along(statements)) {
    if (index %in% skip) next
    block <- blocks[[as.character(index)]]
    if (!is.null(block)) {
      make_lazy_fit_lists(block, env)
      skip <- seq.int(index, block$last)
      next
    }
    st <- statements[[index]]
    head <- call_head(st)
    is_assign <- head %in% c("<-", "=") && is.symbol(st[[2L]])
    if (is_assign) {
      lhs <- as.character(st[[2L]]); rhs <- st[[3L]]; rhs_head <- call_head(rhs)
      if (rhs_head == "scenario_fit") {
        fit_name <- eval(rhs[[2L]], env)
        local({ n <- fit_name; make_lazy(lhs, NULL, env, loader = function() load_fit(n)) })
        next
      }
      if (rhs_head %in% c("scenario_time", "scenario_text", "scenario_plot")) {
        art_name <- eval(rhs[[2L]], env)
        if (list_only || selected(art_name)) {
          eval(st, env)
        } else {
          make_lazy(lhs, rhs[[3L]], env)
        }
        next
      }
      depends <- any(all.vars(rhs) %in% state$lazy) || grepl("^metafor::(rma|selmodel)", rhs_head)
      if (depends && !(lhs %in% all.vars(rhs))) {
        make_lazy(lhs, rhs, env)
        next
      }
      if (depends && list_only) next
      eval(st, env)
      next
    }
    if (head %in% c("scenario_time", "scenario_text", "scenario_plot")) {
      if (list_only || selected(eval(st[[2L]], env))) eval(st, env)
      next
    }
    if (has_scenario_call(st) || "scenario_fit" %in% all.names(st)) {
      eval(st, env)
      next
    }
    if (any(all.vars(st) %in% state$lazy)) next   # checks / prints on fits: not benchmarked
    if (grepl("expect_", head)) next
    eval(st, env)
  }
}

file  <- file.path(main_root, "tests", "scenarios", paste0("test-", scenario, ".R"))
exprs <- parse(file, keep.source = FALSE)
for (e in exprs) {
  head <- call_head(e)
  if (grepl("test_that$", head)) {
    run_block(e[[3L]], new.env(parent = bench_env))
  } else if (head == "if" && grepl("helper-scenarios", paste(deparse(e), collapse = ""))) {
    next
  } else {
    eval(e, bench_env)
  }
}

if (list_only) {
  listing <- do.call(rbind, state$listing)
  timing  <- file.path(main_root, "tests", "scenarios", "timings", paste0(scenario, c(".tsv", ".new.tsv")))
  for (i in seq_along(timing)) {
    if (file.exists(timing[[i]])) {
      t <- utils::read.delim(timing[[i]], stringsAsFactors = FALSE)
      listing[[c("best", "last")[[i]]]] <- t$elapsed[match(paste(listing$type, listing$name), paste(t$type, t$name))]
    }
  }
  utils::write.table(listing, file.path(out_dir, paste0("listing-", scenario, ".tsv")), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  print(listing, row.names = FALSE)
}
