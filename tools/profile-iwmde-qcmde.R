# Profile qCMDE/IWMDE ordinate evaluation against validated cached fits.

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1L]]) else
  file.path("tools", "profile-iwmde-qcmde.R")
script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."),
                              winslash = "/", mustWork = TRUE)
source(file.path(project_root, "tools", "optimized-dll.R"))
preparation <- ensure_optimized_dll(project_root, quiet = TRUE)
if (identical(preparation, "unavailable")) {
  stop("Release preparation requires the 'pkgbuild' package.", call. = FALSE)
}
pkgload::load_all(project_root, quiet = TRUE, compile = FALSE,
                 helpers = FALSE, debug = FALSE)
RoBMA.options(native_threads = 1L)
source(file.path(project_root, "tests", "testthat", "common-functions.R"))

build <- if (identical(preparation, "debug")) "load_all debug opt-out" else "load_all release"
package_path <- getNamespaceInfo(asNamespace("RoBMA"), "path")
dll <- getLoadedDLLs()[["RoBMA"]]
dll_path <- if (is.null(dll)) NA_character_ else dll[["path"]]
git_root <- suppressWarnings(tryCatch(system2("git",
  c("-C", shQuote(project_root), "rev-parse", "--show-toplevel"),
  stdout = TRUE, stderr = FALSE), error = function(e) character()))
source_head <- if (length(git_root) == 1L &&
  identical(tolower(normalizePath(git_root, winslash = "/")), tolower(project_root))) {
  system2("git", c("-C", shQuote(project_root), "rev-parse", "HEAD"), stdout = TRUE)
} else NA_character_
source_status <- if (!is.na(source_head)) system2("git",
  c("-C", shQuote(project_root), "status", "--porcelain=v1"), stdout = TRUE) else NA_character_
identity <- list(source_root = project_root, source_head = source_head,
  source_status = source_status,
  package_source_md5 = if (exists("package_source_md5", mode = "function")) package_source_md5() else NA_character_,
  source_script_md5 = unname(tools::md5sum(script_path)),
  package_version = as.character(getNamespaceVersion("RoBMA")),
  package_path = package_path, dll_path = dll_path,
  dll_md5 = if (is.na(dll_path)) NA_character_ else unname(tools::md5sum(dll_path)),
  build = build, preparation = preparation,
  native_threads = RoBMA.get_option("native_threads"))
cat("source_root: ", project_root, "\nsource_head: ", source_head,
    "\nsource_status: ", paste(source_status, collapse = "; "),
    "\npackage_source_md5: ", identity$package_source_md5,
    "\nsource_script_md5: ", identity$source_script_md5,
    "\npackage: ", identity$package_version, " ", package_path,
    "\ndll: ", dll_path, " md5 ", identity$dll_md5,
    "\nbuild: ", build, " (", preparation, ")",
    "\nnative_threads: ", identity$native_threads,
    "\nElapsed values include CPU and allocation profiling.\n", sep = "")

.profile_iwmde_estimator <- function(name, likelihood_family, parameter, value,
                                     max_samples, density_method = "qCMDE",
                                     n_points = 40L,
                                     normalization_points = 101L) {

  fit     <- load_fit(name, validate = TRUE)
  context <- .iwmde_context(fit)
  memory  <- tempfile(fileext = ".out")
  cpu     <- tempfile(fileext = ".out")

  gc()
  Rprofmem(memory)
  utils::Rprof(cpu, interval = .005)
  on.exit({
    Rprofmem(NULL)
    utils::Rprof(NULL)
    unlink(c(memory, cpu))
  }, add = TRUE)

  elapsed <- system.time({
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = list(
        n_points             = n_points,
        samples              = max_samples,
        normalization_points = normalization_points,
        normalization_prob   = .9999,
        display_grid         = "ordinate"
      ),
      outputs          = "ordinate",
      values           = value,
      parameter_spec   = list(
        type             = "primitive",
        conditional      = NULL,
        conditional_rule = "AND"
      ),
      metadata         = list(parameter = parameter),
      cache            = .iwmde_estimate_cache()
    )
  })[["elapsed"]]
  Rprofmem(NULL)
  utils::Rprof(NULL)

  allocations <- readLines(memory, warn = FALSE)
  sizes       <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
  finite      <- sizes[is.finite(sizes)]
  profile_summary <- utils::summaryRprof(cpu)
  profile         <- list(
    by_total = utils::head(profile_summary[["by.total"]], 20L),
    by_self  = utils::head(profile_summary[["by.self"]], 20L)
  )
  ordinate <- estimate[["posterior_ordinate"]]
  status   <- "accepted"
  if (is.null(ordinate)) {
    ordinate <- estimate[["rejected_posterior_ordinate"]]
    status   <- "rejected"
  }
  if (is.null(ordinate)) {
    stop("Profiled ordinate was unavailable for '", name, "'.", call. = FALSE)
  }
  diagnostics      <- ordinate[["diagnostics"]]
  evaluation_value <- diagnostics[["evaluation_value"]]

  result <- data.frame(
    fit                   = name,
    likelihood_family     = likelihood_family,
    parameter             = parameter,
    density_method        = density_method,
    status                = status,
    requested_value       = value,
    evaluation_value      = evaluation_value,
    max_samples           = max_samples,
    n_points              = n_points,
    normalization_points  = normalization_points,
    profiled_elapsed_seconds = elapsed,
    native_threads        = RoBMA.get_option("native_threads"),
    allocated_megabytes   = sum(finite) / 1024^2,
    largest_allocation_mb = if (length(finite) > 0L) {
      max(finite) / 1024^2
    } else {
      NA_real_
    },
    posterior_ordinate    = ordinate[["ordinate"]],
    stringsAsFactors      = FALSE
  )
  attr(result, "profile") <- profile

  return(result)
}

cases <- list(
  list(
    name                 = "brma.mv_block_mvn_fixed_random_null",
    likelihood_family    = "known-V normal",
    parameter            = "mu",
    value                = 0,
    max_samples          = 40L,
    normalization_points = 51L
  ),
  list(
    name                 = "brma.mv_block_mvn_fixed_random_null",
    likelihood_family    = "known-V normal",
    parameter            = "mu",
    value                = 0,
    max_samples          = 240L,
    normalization_points = 51L
  ),
  list(
    name                 = "brma.mv_block_mvn_fixed_random_null",
    likelihood_family    = "known-V normal",
    parameter            = "mu",
    value                = 0,
    max_samples          = 240L,
    n_points             = 80L,
    normalization_points = 201L
  ),
  list(
    name                 = "brma.mv_block_mvn",
    likelihood_family    = "known-V normal",
    parameter            = "tau",
    value                = 0,
    max_samples          = 240L,
    n_points             = 60L,
    normalization_points = 120L
  ),
  list(
    name                 = "brma.mv_block_mvn",
    likelihood_family    = "known-V normal",
    parameter            = "tau",
    value                = 0,
    max_samples          = 240L,
    density_method       = "IWMDE",
    n_points             = 60L,
    normalization_points = 120L
  ),
  list(
    name              = "bcg_meta-analysis",
    likelihood_family = "ordinary normal",
    parameter         = "mu",
    value             = 0,
    max_samples       = 500L
  ),
  list(
    name              = "nielweise2008_glmm",
    likelihood_family = "Poisson GLMM",
    parameter         = "mu",
    value             = 0,
    max_samples       = 500L
  ),
  list(
    name              = "dat.lehmann2018-3PSM",
    likelihood_family = "selection normal",
    parameter         = "mu",
    value             = 0,
    max_samples       = 500L
  )
)

available <- list_fits(validate = TRUE)
available_cases <- vapply(cases, function(case) case[["name"]] %in% available,
                          logical(1))
omitted <- which(!available_cases)
if (length(omitted) > 0L) {
  cat("Omitted profiler cases (required cached fit unavailable):\n")
  for (index in omitted) cat("- case ", index, ": ", cases[[index]][["name"]], "\n", sep = "")
}
if ("--list" %in% args) {
  print(data.frame(case = seq_along(cases),
                    fit = vapply(cases, `[[`, character(1), "name"),
                    cache_available = available_cases), row.names = FALSE)
  quit(save = "no", status = 0L)
}
cases <- cases[available_cases]
if (length(cases) == 0L) {
  stop(
    "Required cached fits are unavailable; run the test-01 cache suite first.",
    call. = FALSE
  )
}

results <- lapply(cases, function(case) {

  do.call(.profile_iwmde_estimator, case)
})
for (result in results) {
  print(result, row.names = FALSE, digits = 6)
  cat("Top functions by total CPU time:\n")
  print(attr(result, "profile")[["by_total"]], digits = 4)
  cat("Top functions by self CPU time:\n")
  print(attr(result, "profile")[["by_self"]], digits = 4)
}
