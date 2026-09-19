#!/usr/bin/env Rscript
# Micro-benchmark of the Poisson GLMM adaptive Gauss-Hermite kernel
# (`RoBMA_glmm_pois_aghq`) as `add_loo()` drives it on the cached nielweise2008
# fits.
#
#   Rscript tools/bench-glmm-aghq.R [--root=<RoBMA tree>] [--installed]
#       [--cache=<scenario cache root>] [--fits=<regex>] [--reps=<n>]
#       [--threads=<n>] [--out=<tsv>]
#
# Why this exists: the kernel's speed depends on where the linker places the
# statically linked libm bodies it calls (`exp`, `log`, `lgamma`) relative to
# the kernel itself. An unlucky exact placement costs a factor of 2.2 with the
# kernel's own machine code byte-identical and at the same address. It is not
# an R-level or a numerical change - every repetition returns the same result -
# so no test and no scenario fingerprint can see it. Run this after any native
# change; `.agents/instructions/validation.md` records what to do with a slow
# reading.
#
# The tool takes no scratch inputs. It loads each cached fit, wraps the
# package's own `.glmm_pois_aghq()` to record the arguments `add_loo()` hands
# it, rebuilds the `.Call` argument list from those recorded arguments with the
# package's own coercion helpers, and then times the bare `.Call`. The first
# timed repetition's result is checked against the result the live `add_loo()`
# run produced, so a drift between this script's argument list and the
# package's fails loudly instead of timing something else.

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
  file.path("tools", "bench-glmm-aghq.R")
}
# The package this script ships with owns the fit caches, also when --root
# points at another tree.
main_root <- normalizePath(file.path(dirname(script_path), ".."),
                           winslash = "/", mustWork = TRUE)
root  <- normalizePath(opt("root", main_root), winslash = "/", mustWork = TRUE)
cache <- normalizePath(
  opt("cache", file.path(main_root, "tests", "scenarios", "cache")),
  winslash = "/", mustWork = TRUE
)
reps    <- as.integer(opt("reps", "5"))
threads <- as.integer(opt("threads", "1"))
fits_re <- opt("fits", ".")
out_tsv <- opt("out")
if (is.na(reps) || reps < 1L) {
  stop("--reps must be a positive integer.", call. = FALSE)
}

FIT_NAMES <- c("fit_BMA", "fit_brma", "fit_brma_null")
# The published reference state of this benchmark, measured on the development
# machine (AMD Ryzen 9 9950X, R 4.6.0 UCRT, Rtools45 g++ 14.2.0) with one
# thread: `fit_BMA` runs its 9 calls in 12.7 s in a fast layout and in 28 s in
# the slow one. Compare a reading against a run of the same session on a tree
# whose layout is known, not against these numbers on another machine.
REFERENCE_FIT  <- "fit_BMA"
REFERENCE_FAST <- 12.7
REFERENCE_SLOW <- 28.0

Sys.setenv(NOT_CRAN = "true", AGENT = "1")
setwd(root)
if (flag("installed")) {
  suppressPackageStartupMessages(library(RoBMA))
  build <- paste0("installed ", as.character(utils::packageVersion("RoBMA")))
} else {
  # Timings describe the release build, not pkgbuild's default -O0 debug DLL.
  source(file.path(root, "tools", "optimized-dll.R"))
  ensure_optimized_dll(root, quiet = TRUE)
  pkgload::load_all(root, quiet = TRUE, debug = FALSE)
  build <- "load_all release"
}
RoBMA.options(native_threads = threads)
sha <- tryCatch(
  system2("git", c("-C", shQuote(root), "rev-parse", "--short", "HEAD"),
          stdout = TRUE),
  error = function(e) NA_character_
)
dll <- tryCatch(getLoadedDLLs()[["RoBMA"]][["path"]], error = function(e) NA_character_)

cat("root:      ", root, "\n", sep = "")
cat("build:     ", build, " (", as.character(sha), ")\n", sep = "")
cat("dll:       ", as.character(dll), " md5 ",
    if (is.na(dll)) NA_character_ else unname(tools::md5sum(dll)), "\n", sep = "")
cat("cache:     ", cache, "\n", sep = "")
cat("threads:   ", threads, "   reps: ", reps, "\n", sep = "")


object_md5 <- function(value) {
  path <- tempfile("robma-aghq-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(value, path, compress = FALSE)
  unname(tools::md5sum(path))
}

# Rebuild the `.Call` argument list of `RoBMA:::.glmm_pois_aghq()` from the
# arguments it was handed. Keep this in the same order as that function; the
# md5 check below is what catches a drift.
call_arguments <- function(recorded) {
  control <- recorded[["control"]]
  list(
    symbol = if (isTRUE(recorded[["row_sum"]])) {
      "RoBMA_glmm_pois_aghq_row_sum"
    } else {
      "RoBMA_glmm_pois_aghq"
    },
    args = list(
      RoBMA:::.native_integer_vector(recorded[["x1i"]]),
      RoBMA:::.native_integer_vector(recorded[["x2i"]]),
      RoBMA:::.native_numeric_vector(recorded[["t1i"]]),
      RoBMA:::.native_numeric_vector(recorded[["t2i"]]),
      RoBMA:::.native_numeric_matrix(recorded[["mu_samples"]]),
      RoBMA:::.native_numeric_matrix(recorded[["tau_within"]]),
      if (is.null(recorded[["weights"]])) NULL else {
        RoBMA:::.native_numeric_vector(RoBMA:::.glmm_likelihood_weights(
          recorded[["weights"]], ncol(recorded[["mu_samples"]])
        ))
      },
      as.numeric(recorded[["prior_spec"]][["mean"]]),
      as.numeric(recorded[["prior_spec"]][["sd"]]),
      control[["rules"]][["nodes"]],
      control[["rules"]][["log_weights"]],
      control[["tolerance"]],
      control[["consecutive"]],
      control[["mode_tolerance"]]
    )
  )
}

capture_calls <- function(fit_name) {
  path <- file.path(cache, "nielweise2008", paste0(fit_name, ".rds"))
  if (!file.exists(path)) {
    stop("No cached fit at ", path,
         ". Run the nielweise2008 scenario cache first.", call. = FALSE)
  }
  object <- readRDS(path)

  original <- RoBMA:::.glmm_pois_aghq
  store    <- new.env(parent = emptyenv())
  store$recorded <- list()
  wrapper <- function(x1i, x2i, t1i, t2i, mu_samples, tau_within, weights,
                      prior_spec, row_sum = FALSE,
                      control = RoBMA:::.glmm_aghq_control()) {
    value <- original(x1i, x2i, t1i, t2i, mu_samples, tau_within, weights,
                      prior_spec, row_sum, control)
    store$recorded[[length(store$recorded) + 1L]] <- list(
      x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i, mu_samples = mu_samples,
      tau_within = tau_within, weights = weights, prior_spec = prior_spec,
      row_sum = row_sum, control = control, value = value
    )
    value
  }
  utils::assignInNamespace(".glmm_pois_aghq", wrapper, ns = "RoBMA")
  on.exit(
    utils::assignInNamespace(".glmm_pois_aghq", original, ns = "RoBMA"),
    add = TRUE
  )

  started <- proc.time()[["elapsed"]]
  invisible(add_loo(object, parallel = FALSE))
  elapsed <- proc.time()[["elapsed"]] - started

  recorded <- store$recorded
  if (length(recorded) == 0L) {
    stop("'", fit_name, "' made no Poisson AGHQ call; the cached fit or the ",
         "log-likelihood route has changed.", call. = FALSE)
  }
  list(
    calls    = lapply(recorded, call_arguments),
    expected = object_md5(lapply(recorded, `[[`, "value")),
    S        = nrow(recorded[[1L]][["mu_samples"]]),
    prep     = elapsed
  )
}

run_once <- function(calls) {
  values  <- vector("list", length(calls))
  seconds <- numeric(length(calls))
  for (index in seq_along(calls)) {
    call    <- calls[[index]]
    started <- proc.time()[["elapsed"]]
    values[[index]] <- do.call(
      ".Call", c(list(call[["symbol"]]), call[["args"]], list(PACKAGE = "RoBMA"))
    )
    seconds[[index]] <- proc.time()[["elapsed"]] - started
  }
  list(seconds = sum(seconds), per_call = seconds, md5 = object_md5(values))
}


selected <- FIT_NAMES[grepl(fits_re, FIT_NAMES, perl = TRUE)]
if (length(selected) == 0L) {
  stop("--fits selected none of: ", paste(FIT_NAMES, collapse = ", "),
       call. = FALSE)
}

captured <- list()
for (fit_name in selected) {
  captured[[fit_name]] <- capture_calls(fit_name)
  cat(sprintf(
    "captured:  %-14s %2d calls, S = %d (live add_loo() %.1f s)\n",
    fit_name, length(captured[[fit_name]][["calls"]]),
    captured[[fit_name]][["S"]], captured[[fit_name]][["prep"]]
  ))
}

timings <- stats::setNames(
  vector("list", length(selected)), selected
)
digests <- stats::setNames(character(length(selected)), selected)
for (rep in seq_len(reps)) {
  for (fit_name in selected) {
    run <- run_once(captured[[fit_name]][["calls"]])
    timings[[fit_name]] <- c(timings[[fit_name]], run[["seconds"]])
    if (rep == 1L) {
      digests[[fit_name]] <- run[["md5"]]
      if (!identical(run[["md5"]], captured[[fit_name]][["expected"]])) {
        stop("The replayed ", fit_name, " calls do not reproduce the result of ",
             "the live add_loo() run; this script's argument list and ",
             "RoBMA:::.glmm_pois_aghq() have drifted apart.", call. = FALSE)
      }
    } else if (!identical(run[["md5"]], digests[[fit_name]])) {
      stop("The ", fit_name, " kernel returned a different result between ",
           "repetitions.", call. = FALSE)
    }
    cat(sprintf("rep %d:     %-14s %7.2f s\n", rep, fit_name, run[["seconds"]]))
  }
}

rows <- do.call(rbind, lapply(selected, function(fit_name) {
  seconds <- timings[[fit_name]]
  data.frame(
    fit      = fit_name,
    calls    = length(captured[[fit_name]][["calls"]]),
    median_s = stats::median(seconds),
    min_s    = min(seconds),
    max_s    = max(seconds),
    reps_s   = paste(sprintf("%.2f", seconds), collapse = ","),
    md5      = digests[[fit_name]],
    stringsAsFactors = FALSE
  )
}))
total <- data.frame(
  fit      = "<total>",
  calls    = sum(rows[["calls"]]),
  median_s = stats::median(Reduce(`+`, timings)),
  min_s    = min(Reduce(`+`, timings)),
  max_s    = max(Reduce(`+`, timings)),
  reps_s   = paste(sprintf("%.2f", Reduce(`+`, timings)), collapse = ","),
  md5      = "",
  stringsAsFactors = FALSE
)
rows <- rbind(rows, total)

cat("\n")
cat(sprintf("%-14s %6s %10s %10s %10s  %s\n",
            "fit", "calls", "median s", "min s", "max s", "per-rep s"))
for (index in seq_len(nrow(rows))) {
  cat(sprintf("%-14s %6d %10.2f %10.2f %10.2f  %s\n",
              rows[["fit"]][[index]], rows[["calls"]][[index]],
              rows[["median_s"]][[index]], rows[["min_s"]][[index]],
              rows[["max_s"]][[index]], rows[["reps_s"]][[index]]))
}
cat("\n")
for (index in seq_len(nrow(rows) - 1L)) {
  cat(sprintf("md5 %-14s %s\n", rows[["fit"]][[index]], rows[["md5"]][[index]]))
}

reference <- rows[rows[["fit"]] == REFERENCE_FIT, "median_s"]
cat("\nReference (AMD Ryzen 9 9950X, R 4.6.0 UCRT, one thread): ",
    REFERENCE_FIT, " takes about ", sprintf("%.1f", REFERENCE_FAST),
    " s in a fast layout and about ", sprintf("%.0f", REFERENCE_SLOW),
    " s in the slow one.\n", sep = "")
if (length(reference) == 1L && is.finite(reference)) {
  cat("This run: ", sprintf("%.2f", reference), " s (",
      sprintf("%.2f", reference / REFERENCE_FAST),
      "x the fast reference).\n", sep = "")
  if (reference > 1.5 * REFERENCE_FAST) {
    cat("SLOW: this build is in the collided layout. See ",
        "'.agents/instructions/validation.md', section 'Performance ",
        "investigation'.\n", sep = "")
  }
}
cat("A reading is only comparable within one machine and session; measure the ",
    "tree you are comparing against in the same session.\n", sep = "")

if (!is.null(out_tsv)) {
  rows[["root"]]    <- root
  rows[["sha"]]     <- as.character(sha)
  rows[["threads"]] <- threads
  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(rows, out_tsv, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  cat("wrote: ", out_tsv, "\n", sep = "")
}
