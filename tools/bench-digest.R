#!/usr/bin/env Rscript
# Digest a bench-scenario run directory: per artifact elapsed, native share, top
# self-time functions. Writes `digest.txt` next to the results it reads.
#
#   Rscript tools/bench-digest.R <run dir> [min seconds]
args    <- commandArgs(trailingOnly = TRUE)
run_dir <- args[[1L]]
min_s   <- if (length(args) > 1L) as.numeric(args[[2L]]) else 2
results <- utils::read.delim(file.path(run_dir, "results.tsv"), stringsAsFactors = FALSE, quote = "")
results <- results[results$elapsed >= min_s, , drop = FALSE]
results <- results[order(-results$elapsed), , drop = FALSE]
native_names <- c("\".Call\"", "\".External\"", "\".C\"")
out <- character()
for (i in seq_len(nrow(results))) {
  stub <- file.path(run_dir, paste0(results$scenario[[i]], "--", gsub("[^A-Za-z0-9_.-]", "_", results$name[[i]])))
  if (!file.exists(paste0(stub, ".Rprof"))) stub <- file.path(run_dir, gsub("[^A-Za-z0-9_.-]", "_", results$name[[i]]))
  prof <- paste0(stub, ".Rprof")
  line <- sprintf("%-18s %-5s %-46s %8.1fs %5.2fGB %s", results$scenario[[i]], results$type[[i]],
                  results$name[[i]], results$elapsed[[i]], results$memory_gb[[i]],
                  if (results$status[[i]] == "ok") "" else substr(results$status[[i]], 1, 80))
  if (file.exists(prof)) {
    s <- tryCatch(utils::summaryRprof(prof), error = function(e) NULL)
    if (!is.null(s)) {
      self   <- s$by.self
      total  <- s$sampling.time
      native <- sum(self$self.time[rownames(self) %in% native_names])
      gcs    <- sum(self$self.time[rownames(self) %in% c("\"gc\"", "\"<GC>\"")])
      top    <- utils::head(self[!rownames(self) %in% native_names, , drop = FALSE], 6L)
      line <- paste0(line, sprintf("  sampled %6.1fs native %4.0f%%  | ", total, 100 * native / total),
                     paste0(gsub("\"", "", rownames(top)), " ", round(100 * top$self.time / total), "%",
                            collapse = ", "))
    }
  }
  out <- c(out, line)
}
writeLines(out, file.path(run_dir, "digest.txt"))
cat(out, sep = "\n")
