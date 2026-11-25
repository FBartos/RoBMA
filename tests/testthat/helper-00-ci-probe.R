# TEMPORARY CI PROBE (scratch branch only, never merge): records the session state at the top of every test
# file, so that the CI run shows where it changes.
.ci_probe_env <- new.env()
.ci_probe_dir <- function() {
  dir <- file.path(Sys.getenv("GITHUB_WORKSPACE", unset = tempdir()), "ci-diagnostic")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}
.ci_probe_line <- function(...) {
  cat(paste0(format(Sys.time(), "%H:%M:%S"), " ", paste0(..., collapse = ""), "\n"),
      file = file.path(.ci_probe_dir(), "probe.txt"), append = TRUE)
}
probe_state <- function(tag) {
  tryCatch({
    factories <- tryCatch(rjags::list.factories("sampler"), error = function(e) NULL)
    factories <- if (is.null(factories)) "NA" else paste0(factories$factory, ifelse(factories$status, "", "(off)"), collapse = " ")
    options_now <- tryCatch(RoBMA.options(), error = function(e) list(error = conditionMessage(e)))
    state <- list(
      rng       = paste(RNGkind(), collapse = "/"),
      seed      = exists(".Random.seed", envir = globalenv(), inherits = FALSE),
      options   = paste(utils::capture.output(utils::str(options_now, give.attr = FALSE)), collapse = " ; "),
      factories = factories,
      modules   = paste(tryCatch(rjags::list.modules(), error = function(e) "NA"), collapse = ","),
      r_options = paste(getOption("digits"), getOption("scipen"), getOption("OutDec"), getOption("warn"),
                        getOption("mc.cores", "NULL"), getOption("contrasts")[1], sep = "/"),
      env       = paste(Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "ROBMA_TEST_PROFILE", "NOT_CRAN")), collapse = "/"),
      locale    = Sys.getlocale("LC_COLLATE")
    )
    changed <- names(state)[vapply(names(state), function(n) !identical(state[[n]], .ci_probe_env$last[[n]]), logical(1))]
    .ci_probe_line(tag, " | changed: ", if (length(changed)) paste(changed, collapse = ",") else "-")
    for (n in changed) .ci_probe_line("    ", n, " = ", state[[n]])
    .ci_probe_env$last <- state
  }, error = function(e) .ci_probe_line(tag, " | probe error: ", conditionMessage(e)))
  invisible()
}
probe_fit <- function(tag, fit) {
  tryCatch({
    chains <- lapply(coda::as.mcmc.list(fit$fit), as.matrix)
    m <- do.call(rbind, chains)
    .ci_probe_line(tag, " | fit hash ", rlang::hash(unname(m)), " | dim ", paste(dim(m), collapse = "x"), " | samplers ",
                   paste(names(table(fit$fit$samplers$Sampler)), table(fit$fit$samplers$Sampler), collapse = "; "))
    for (col in grep("indicator", colnames(m), value = TRUE)) {
      .ci_probe_line("    ", col, " ", paste(vapply(chains, function(x) sprintf("mean %.4f switches %d", mean(x[, col]), sum(diff(x[, col]) != 0)), ""), collapse = " | "))
    }
    .ci_probe_line("    first draw: ", paste(signif(m[1, seq_len(min(8, ncol(m)))], 10), collapse = " "))
  }, error = function(e) .ci_probe_line(tag, " | fit probe error: ", conditionMessage(e)))
  invisible()
}
