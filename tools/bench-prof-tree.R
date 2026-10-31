#!/usr/bin/env Rscript
# Text call tree from a raw Rprof file written by tools/bench-scenario.R.
#
#   Rscript tools/bench-prof-tree.R <file.Rprof> [min_pct=2] [max_depth=40] [root_regex]
args      <- commandArgs(trailingOnly = TRUE)
file      <- args[[1L]]
min_pct   <- if (length(args) > 1L) as.numeric(args[[2L]]) else 2
max_depth <- if (length(args) > 2L) as.integer(args[[3L]]) else 40L
root_rx   <- if (length(args) > 3L) args[[4L]] else NULL
lines <- readLines(file, warn = FALSE)[-1L]
lines <- lines[nzchar(lines)]
stacks <- lapply(strsplit(lines, " ", fixed = TRUE), function(x) rev(gsub("\"", "", x[nzchar(x)])))
skip <- c("tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch", "eval", "withCallingHandlers",
          "suppressMessages", "suppressWarnings", "withVisible", "bench_run", "run_block",
          "scenario_plot", "scenario_text", "scenario_time", "<Anonymous>", "do.call", "FUN", "lapply",
          "vapply", "sapply", "force", "try", "withRestarts", "withOneRestart", "doWithOneRestart",
          "standardGeneric", "NextMethod", "mapply", "Map", "simplify2array", "unlist")
stacks <- lapply(stacks, function(s) s[!s %in% skip])
if (!is.null(root_rx)) {
  stacks <- lapply(stacks, function(s) {
    hit <- grep(root_rx, s)
    if (length(hit) == 0L) return(NULL)
    s[hit[[1L]]:length(s)]
  })
  stacks <- stacks[!vapply(stacks, is.null, logical(1))]
}
total <- length(stacks)
cat("samples:", total, "\n")
node <- function(stacks, depth, prefix) {
  if (depth > max_depth || length(stacks) == 0L) return(invisible())
  heads <- vapply(stacks, function(s) if (length(s) >= depth) s[[depth]] else NA_character_, character(1))
  tab   <- sort(table(heads), decreasing = TRUE)
  for (name in names(tab)) {
    pct <- 100 * tab[[name]] / total
    if (pct < min_pct) next
    sub  <- stacks[!is.na(heads) & heads == name]
    self <- sum(vapply(sub, length, integer(1)) == depth)
    cat(sprintf("%s%-s  %.1f%%  (self %.1f%%)\n", prefix, name, pct, 100 * self / total))
    node(sub, depth + 1L, paste0(prefix, "  "))
  }
}
node(stacks, 1L, "")
