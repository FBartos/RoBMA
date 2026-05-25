# ============================================================================ #
# CONFIGURATION: Set to TRUE to regenerate reference files, FALSE to run tests
# ============================================================================ #
if (!exists("GENERATE_REFERENCE_FILES")) {
  GENERATE_REFERENCE_FILES <- FALSE
}
if (!exists("FIT_CACHE_VERSION")) {
  FIT_CACHE_VERSION <- 3L
}

.common_functions_dir <- function() {

  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    ofile <- frame[["ofile"]]
    if (is.null(ofile)) "" else ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]

  if (length(ofiles) > 0L) {
    return(dirname(normalizePath(ofiles[length(ofiles)], winslash = "/", mustWork = FALSE)))
  }

  fallback <- file.path("tests", "testthat", "common-functions.R")
  if (file.exists(fallback)) {
    return(dirname(normalizePath(fallback, winslash = "/", mustWork = TRUE)))
  }

  return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

# Get the directory where prefitted models are stored. Local development uses
# a persistent ignored folder; CRAN checks fall back to tempdir().
test_files_dir <- Sys.getenv("ROBMA_TEST_FILES_DIR")
if (test_files_dir == "") {
  on_cran <- get("on_cran", envir = asNamespace("testthat"), inherits = FALSE)
  test_files_dir <- if (on_cran()) {
    file.path(tempdir(), "RoBMA_test_files")
  } else {
    file.path(.common_functions_dir(), "test_files")
  }
}
test_files_dir <- normalizePath(test_files_dir, winslash = "/", mustWork = FALSE)

# Setup directory for saving fitted models
temp_fits_dir     <- file.path(test_files_dir, "fits")
temp_info_dir     <- file.path(test_files_dir, "info")
temp_metadata_dir <- file.path(test_files_dir, "metadata")
temp_temp_dir     <- file.path(test_files_dir, "temp")
if (!dir.exists(test_files_dir)) dir.create(test_files_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_fits_dir)) dir.create(temp_fits_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_info_dir)) dir.create(temp_info_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_metadata_dir)) dir.create(temp_metadata_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_temp_dir)) dir.create(temp_temp_dir, showWarnings = FALSE, recursive = TRUE)

# Set environment variable so other test files can locate pre-fitted models
Sys.setenv(ROBMA_TEST_FILES_DIR = test_files_dir)

# Use skip_if_no_fits() for tests that need pre-fitted models.

# ============================================================================ #
# HELPER FUNCTIONS: Reference File Testing
# ============================================================================ #

# Process reference file: save if GENERATE_REFERENCE_FILES=TRUE, test otherwise
test_reference_table <- function(table, filename, info_msg = NULL,
                                 print_dir = REFERENCE_DIR) {
  if (GENERATE_REFERENCE_FILES) {
    # Save mode
    if (!dir.exists(print_dir)) {
      dir.create(print_dir, recursive = TRUE)
    }
    writeLines(
      capture_output_lines(table, print = TRUE, width = 150),
      file.path(print_dir, filename)
    )
  } else {
    # Test mode
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      actual_output <- capture_output_lines(table, print = TRUE, width = 150)
      expect_equal(actual_output, expected_output, info = info_msg)
    } else {
      skip(paste("Reference file", filename, "not found."))
    }
  }
}

test_reference_text <- function(text, filename, info_msg = NULL,
                                print_dir = REFERENCE_DIR) {
  if (GENERATE_REFERENCE_FILES) {
    # Save mode
    if (!dir.exists(print_dir)) {
      dir.create(print_dir, recursive = TRUE)
    }
    writeLines(text, file.path(print_dir, filename))
  } else {
    # Test mode
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      expected_output <- paste0(expected_output, collapse = "\n")
      expect_equal(text, expected_output, info = info_msg)
    } else {
      skip(paste("Reference file", filename, "not found."))
    }
  }
}

vdiffr_snapshots_available <- function() {

  snap_dir <- testthat::test_path("_snaps")

  return(
    dir.exists(snap_dir) &&
      length(list.files(snap_dir, pattern = "\\.svg$", recursive = TRUE)) > 0L
  )
}

skip_if_no_vdiffr_snapshots <- function() {

  if (!vdiffr_snapshots_available() &&
      !is_true_env("ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS")) {
    testthat::skip(paste(
      "vdiffr snapshots are not available.",
      "`tests/testthat/_snaps` is intentionally excluded from source builds.",
      "Set ROBMA_TEST_ALLOW_MISSING_SNAPSHOTS=TRUE only when regenerating snapshots."
    ))
  }

  return(invisible(FALSE))
}

# ============================================================================ #
# HELPER FUNCTIONS: Cached Fit Catalog
# ============================================================================ #

fit_catalog <- function() {

  catalog <- data.frame(
    name = c(
      "bcg_meta-analysis",
      "bcg_meta-regression",
      "bcg_meta-regression2",
      "bcg_meta-regression2b",
      "bcg_meta-regression3",
      "bcg_meta-regression3b",
      "bcg_meta-regression4",
      "bcg_meta-regression4b",
      "bangertdrowns2004_location-scale",
      "konstantopoulos2011_3lvl",
      "konstantopoulos2011_3lvl2",
      "bcg_glmm",
      "bcg_glmm_reg",
      "nielweise2008_glmm",
      "bcg_glmm_3lvl_scale",
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESE_neg",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-3PSM",
      "dat.lehmann2018-4PSM",
      "dat.lehmann2018-3PSM_neg",
      "dat.lehmann2018-3PSMreg",
      "dat.lehmann2018_BMA.norm",
      "dat.lehmann2018_BMA.norm_custom",
      "dat.lehmann2018_BMA.norm_mods",
      "dat.lehmann2018_BMA.norm_scale",
      "bcg_BMA.glmm",
      "nielweise2008_BMA.glmm",
      "bcg_BMA.glmm_custom",
      "bcg_BMA.glmm_3lvl_location_scale",
      "dat.lehmann2018_RoBMA",
      "dat.lehmann2018_RoBMA_custom",
      "dat.lehmann2018_RoBMA_mods",
      "dat.lehmann2018_RoBMA_mods2",
      "dat.lehmann2018_RoBMA_3lvl_mods_scale"
    ),
    class = c(
      rep("brma.norm", 11),
      rep("brma.glmm", 4),
      rep("bPET", 3),
      rep("bPEESE", 3),
      rep("bselmodel", 4),
      rep("BMA.norm", 4),
      rep("BMA.glmm", 4),
      rep("RoBMA", 5)
    ),
    family = c(
      rep("norm", 11),
      rep("glmm", 4),
      rep("norm", 3),
      rep("norm", 3),
      rep("norm", 4),
      rep("norm", 4),
      rep("glmm", 4),
      rep("norm", 5)
    ),
    source_file = c(
      rep("test-01-brma.norm.R", 11),
      rep("test-01-brma.glmm.R", 4),
      rep("test-01-bPET.R", 3),
      rep("test-01-bPEESE.R", 3),
      rep("test-01-bselmodel.R", 4),
      rep("test-01-BMA.norm.R", 4),
      rep("test-01-BMA.glmm.R", 4),
      rep("test-01-RoBMA.R", 5)
    ),
    has_metafor = c(
      TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, FALSE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE,
      rep(FALSE, 13)
    ),
    has_loo = TRUE,
    has_waic = FALSE,
    has_marglik = c(
      rep(TRUE, 25),
      rep(FALSE, 13)
    ),
    tier = c(
      "core", "core", "extended", "extended", "core", "extended", "extended", "extended",
      "core", "core", "extended",
      "core", "core", "core", "extended",
      "core", "extended", "extended",
      "core", "extended", "core",
      "core", "extended", "extended", "extended",
      "core", "extended", "core", "extended",
      "core", "core", "extended", "extended",
      "core", "extended", "core", "extended", "extended"
    ),
    stringsAsFactors = FALSE
  )

  catalog[["features"]] <- I(list(
    c("metafor", "normal", "simple"),
    c("metafor", "normal", "mods"),
    c("metafor", "normal", "factor_mods"),
    c("normal", "factor_mods", "meandif"),
    c("metafor", "normal", "interaction"),
    c("normal", "interaction", "meandif"),
    c("metafor", "normal", "interaction"),
    c("normal", "interaction", "meandif"),
    c("metafor", "normal", "mods", "scale"),
    c("metafor", "normal", "multilevel"),
    c("metafor", "normal", "multilevel", "mods"),
    c("metafor", "glmm", "binomial"),
    c("metafor", "glmm", "binomial", "mods"),
    c("metafor", "glmm", "poisson"),
    c("glmm", "binomial", "multilevel", "scale"),
    c("metafor", "PET"),
    c("metafor", "PET", "negative"),
    c("metafor", "PET", "mods"),
    c("metafor", "PEESE"),
    c("metafor", "PEESE", "negative"),
    c("metafor", "PEESE", "mods"),
    c("metafor", "selection"),
    c("metafor", "selection", "steps"),
    c("metafor", "selection", "negative"),
    c("metafor", "selection", "mods"),
    c("BMA.norm", "normal", "simple"),
    c("BMA.norm", "normal", "custom_priors"),
    c("BMA.norm", "normal", "mods"),
    c("BMA.norm", "normal", "scale"),
    c("BMA.glmm", "glmm", "binomial"),
    c("BMA.glmm", "glmm", "poisson"),
    c("BMA.glmm", "glmm", "binomial", "custom_priors"),
    c("BMA.glmm", "glmm", "binomial", "multilevel", "mods", "scale"),
    c("RoBMA", "normal", "simple", "selection", "PET"),
    c("RoBMA", "normal", "custom_priors"),
    c("RoBMA", "normal", "mods"),
    c("RoBMA", "normal", "interaction"),
    c("RoBMA", "normal", "multilevel", "mods", "scale")
  ))

  return(catalog)
}

fit_catalog_entry <- function(name) {

  catalog <- fit_catalog()
  index   <- match(name, catalog[["name"]])

  if (is.na(index)) {
    return(NULL)
  }

  return(catalog[index, , drop = FALSE])
}

catalog_group_fits <- function(name) {

  catalog <- fit_catalog()

  if (name %in% catalog[["name"]]) {
    return(name)
  }

  if (name %in% catalog[["class"]]) {
    return(catalog[catalog[["class"]] == name, "name"])
  }

  source_file <- paste0("test-01-", name, ".R")
  if (source_file %in% catalog[["source_file"]]) {
    return(catalog[catalog[["source_file"]] == source_file, "name"])
  }

  return(character())
}

catalog_fits <- function(feature, class, family, has_metafor, has_loo,
                         has_waic, has_marglik, tier) {

  catalog <- fit_catalog()

  if (!missing(feature) && !is.null(feature)) {
    catalog <- catalog[vapply(catalog[["features"]], function(x) all(feature %in% x), TRUE), , drop = FALSE]
  }
  if (!missing(class) && !is.null(class)) {
    catalog <- catalog[catalog[["class"]] %in% class, , drop = FALSE]
  }
  if (!missing(family) && !is.null(family)) {
    catalog <- catalog[catalog[["family"]] %in% family, , drop = FALSE]
  }
  if (!missing(has_metafor) && !is.null(has_metafor)) {
    catalog <- catalog[catalog[["has_metafor"]] %in% has_metafor, , drop = FALSE]
  }
  if (!missing(has_loo) && !is.null(has_loo)) {
    catalog <- catalog[catalog[["has_loo"]] %in% has_loo, , drop = FALSE]
  }
  if (!missing(has_waic) && !is.null(has_waic)) {
    catalog <- catalog[catalog[["has_waic"]] %in% has_waic, , drop = FALSE]
  }
  if (!missing(has_marglik) && !is.null(has_marglik)) {
    catalog <- catalog[catalog[["has_marglik"]] %in% has_marglik, , drop = FALSE]
  }
  if (!missing(tier) && !is.null(tier)) {
    catalog <- catalog[catalog[["tier"]] %in% tier, , drop = FALSE]
  }

  return(catalog[["name"]])
}

.source_hash_call_name <- function(expr) {

  if (!is.call(expr)) {
    return("")
  }

  fun <- expr[[1]]
  if (is.symbol(fun)) {
    return(as.character(fun))
  }
  if (is.call(fun) && identical(as.character(fun[[1]]), "::")) {
    return(as.character(fun[[3]]))
  }

  return("")
}

.source_hash_skip_call <- function(expr) {

  call_name <- .source_hash_call_name(expr)

  return(
    call_name %in% c("context", "source", "skip_on_cran",
                     "skip_if_not_installed", "skip_refit_if_cached") ||
      grepl("^expect_", call_name)
  )
}

.source_hash_normalize_expr <- function(expr) {

  if (is.call(expr)) {
    call_name <- .source_hash_call_name(expr)

    if (.source_hash_skip_call(expr)) {
      return(character())
    }

    if (identical(call_name, "test_that")) {
      body <- if (length(expr) >= 3L) expr[[3]] else NULL
      return(.source_hash_normalize_expr(body))
    }

    if (identical(call_name, "{")) {
      normalized <- unlist(
        lapply(as.list(expr)[-1], .source_hash_normalize_expr),
        use.names = FALSE
      )
      return(normalized[nzchar(normalized)])
    }
  }

  return(paste(deparse(expr, width.cutoff = 500L), collapse = "\n"))
}

source_file_md5 <- function(source_file) {

  if (is.null(source_file) || is.na(source_file) || !nzchar(source_file)) {
    return(NA_character_)
  }

  path <- testthat::test_path(source_file)
  if (!file.exists(path)) {
    return(NA_character_)
  }

  parsed <- try(parse(path), silent = TRUE)
  if (inherits(parsed, "try-error")) {
    lines <- readLines(path, warn = FALSE)
    lines <- sub("#.*$", "", lines)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
  } else {
    lines <- unlist(lapply(parsed, .source_hash_normalize_expr), use.names = FALSE)
    lines <- as.character(lines)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
  }

  normalized <- tempfile("robma-fit-source-", fileext = ".R")
  writeLines(lines, normalized, useBytes = TRUE)
  on.exit(unlink(normalized), add = TRUE)

  return(unname(tools::md5sum(normalized)))
}

.package_source_md5_cache <- new.env(parent = emptyenv())

.fit_cache_required_source_files <- function() {

  # Keep this scoped to code that can change saved fit objects or cached fit
  # extensions. Post-fit methods are tested against the cached objects.
  r_files   <- c(
    "R/BMA.glmm.R",
    "R/BMA.norm.R",
    "R/RoBMA.R",
    "R/bPEESE.R",
    "R/bPET.R",
    "R/brma.glmm.R",
    "R/brma.norm.R",
    "R/bselmodel.R",
    "R/fit.R",
    "R/input-data.R",
    "R/input-object.R",
    "R/input-priors.R",
    "R/input-priors-assignment.R",
    "R/input-priors-check-list.R",
    "R/input-priors-documentation.R",
    "R/input-priors-formula.R",
    "R/input-priors-heterogeneity-allocation.R",
    "R/loo.R",
    "R/marglik.R",
    "R/log-lik-cluster.R",
    "R/log-lik-cluster-glmm.R",
    "R/log-lik-cluster-normal.R",
    "R/log-lik.R",
    "R/pdf-outcome-glmm.R",
    "R/pdf-outcome-normal.R",
    "R/pdf-utils.R",
    "R/priors.R",
    "R/selection-mapping.R",
    "R/utilities.R",
    "R/zzz.R"
  )
  src_files <- c(
    "src/RoBMA.cc",
    "src/init.c",
    "src/distributions/DSELNORMKERNEL.cc",
    "src/distributions/DSELNORMKERNEL.h",
    "src/distributions/DSELNORMSTEP.cc",
    "src/distributions/DSELNORMSTEP.h",
    "src/distributions/DSELNORMSTEPSWITCH.cc",
    "src/distributions/DSELNORMSTEPSWITCH.h",
    "src/distributions/DWB.cc",
    "src/distributions/DWB.h",
    "src/distributions/DWN.cc",
    "src/distributions/DWN.h",
    "src/distributions/DWP.cc",
    "src/distributions/DWP.h",
    "src/r-glmm.cc",
    "src/r-selnorm.cc",
    "src/r-selnorm-common.cc.inc",
    "src/r-selnorm-funnel-zcurve.cc.inc",
    "src/r-selnorm-kernel.cc.inc",
    "src/r-selnorm-loglik.cc.inc",
    "src/selnorm/selnorm.cc",
    "src/selnorm/selnorm-api.cc.inc",
    "src/selnorm/selnorm-boundary.cc.inc",
    "src/selnorm/selnorm-phack.cc.inc",
    "src/selnorm/selnorm-probability.cc.inc",
    "src/selnorm/selnorm-step.cc.inc",
    "src/selnorm/selnorm.h"
  )

  return(sort(c(r_files, src_files)))
}

.path_ancestors <- function(path) {

  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  out  <- character()

  repeat {
    out    <- c(out, path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      break
    }
    path <- parent
  }

  return(unique(out))
}

.fit_cache_source_root_candidates <- function() {

  package_name <- "RoBMA"
  roots        <- c(
    testthat::test_path("..", ".."),
    getwd(),
    system.file(package = package_name)
  )
  roots <- roots[nzchar(roots)]
  roots <- normalizePath(roots, winslash = "/", mustWork = FALSE)

  bases      <- unique(unlist(lapply(roots, .path_ancestors), use.names = FALSE))
  candidates <- unique(c(
    roots,
    file.path(bases, package_name),
    file.path(bases, "00_pkg_src", package_name)
  ))

  return(normalizePath(candidates, winslash = "/", mustWork = FALSE))
}

.fit_cache_source_tree_available <- function(package_root) {

  package_root <- normalizePath(package_root, winslash = "/", mustWork = FALSE)
  source_files <- .fit_cache_required_source_files()

  return(all(file.exists(file.path(package_root, source_files))))
}

.fit_cache_source_root <- function() {

  if (exists("source_root", envir = .package_source_md5_cache, inherits = FALSE)) {
    return(get("source_root", envir = .package_source_md5_cache, inherits = FALSE))
  }

  candidates <- .fit_cache_source_root_candidates()
  for (candidate in candidates) {
    if (.fit_cache_source_tree_available(candidate)) {
      assign("source_root", candidate, envir = .package_source_md5_cache)
      return(candidate)
    }
  }

  assign("source_root", NA_character_, envir = .package_source_md5_cache)
  return(NA_character_)
}

.fit_cache_source_files <- function(package_root = NULL, relative = FALSE) {

  if (is.null(package_root)) {
    package_root <- .fit_cache_source_root()
    if (is.na(package_root)) {
      return(character())
    }
  } else {
    package_root <- normalizePath(package_root, winslash = "/", mustWork = FALSE)
  }

  source_files <- .fit_cache_required_source_files()
  source_paths <- file.path(package_root, source_files)
  if (!all(file.exists(source_paths))) {
    return(character())
  }

  if (relative) {
    return(source_files)
  }

  return(sort(normalizePath(
    file.path(package_root, source_files),
    winslash = "/",
    mustWork = TRUE
  )))
}

.fit_cache_source_file_md5 <- function(path) {

  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "r")) {
    parsed <- try(parse(file = path, keep.source = FALSE), silent = TRUE)
    if (!inherits(parsed, "try-error")) {
      lines <- unlist(lapply(parsed, function(expr) {
        paste(deparse(expr, width.cutoff = 500L), collapse = "\n")
      }), use.names = FALSE)
      lines <- as.character(lines)

      normalized <- tempfile("robma-package-source-", fileext = ".R")
      writeLines(lines, normalized, useBytes = TRUE)
      on.exit(unlink(normalized), add = TRUE)

      return(unname(tools::md5sum(normalized)))
    }
  }

  return(unname(tools::md5sum(path)))
}

package_source_md5 <- function() {

  if (exists("value", envir = .package_source_md5_cache, inherits = FALSE)) {
    return(get("value", envir = .package_source_md5_cache, inherits = FALSE))
  }

  package_root <- .fit_cache_source_root()
  if (is.na(package_root)) {
    assign("value", NA_character_, envir = .package_source_md5_cache)
    return(NA_character_)
  }

  source_files <- .fit_cache_source_files(package_root = package_root)
  if (length(source_files) == 0L) {
    assign("value", NA_character_, envir = .package_source_md5_cache)
    return(NA_character_)
  }

  file_hashes    <- vapply(source_files, .fit_cache_source_file_md5, character(1))
  relative_files <- substring(source_files, nchar(package_root) + 2L)
  hash_input     <- paste(
    paste0(relative_files, ":", unname(file_hashes)),
    collapse = "\n"
  )

  normalized <- tempfile("robma-package-source-", fileext = ".txt")
  writeLines(hash_input, normalized, useBytes = TRUE)
  on.exit(unlink(normalized), add = TRUE)

  value <- unname(tools::md5sum(normalized))
  assign("value", value, envir = .package_source_md5_cache)

  return(value)
}

fit_cache_metadata <- function(name, fit, info = NULL) {

  entry       <- fit_catalog_entry(name)
  source_file <- if (is.null(entry)) NA_character_ else entry[["source_file"]]

  metadata <- list(
    version            = FIT_CACHE_VERSION,
    name               = name,
    saved_at           = format(Sys.time(), usetz = TRUE),
    fit_class          = class(fit),
    source_file        = source_file,
    source_file_md5    = source_file_md5(source_file),
    package_source_md5 = package_source_md5(),
    has_loo            = !is.null(fit[["loo"]]),
    has_waic           = !is.null(fit[["waic"]]),
    has_marglik        = !is.null(fit[["marglik"]]),
    has_metafor_info   = !is.null(info) && "metafor" %in% names(info) && !is.null(info[["metafor"]])
  )

  return(metadata)
}

fit_cache_paths <- function(name) {

  return(list(
    fit      = file.path(temp_fits_dir, paste0(name, ".RDS")),
    info     = file.path(temp_info_dir, paste0(name, ".RDS")),
    metadata = file.path(temp_metadata_dir, paste0(name, ".RDS")),
    marker   = file.path(temp_temp_dir, paste0(name, ".txt"))
  ))
}

load_fit_metadata <- function(name) {

  path <- fit_cache_paths(name)[["metadata"]]
  metadata <- suppressWarnings(try(readRDS(path), silent = TRUE))

  if (inherits(metadata, "try-error")) {
    return(NULL)
  }

  return(metadata)
}

read_cached_fit <- function(name) {

  path <- fit_cache_paths(name)[["fit"]]
  fit  <- suppressWarnings(try(readRDS(path), silent = TRUE))

  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  return(fit)
}

read_cached_info <- function(name) {

  path <- fit_cache_paths(name)[["info"]]
  info <- suppressWarnings(try(readRDS(path), silent = TRUE))

  if (inherits(info, "try-error")) {
    return(list())
  }

  return(info)
}

.fit_object_cache  <- new.env(parent = emptyenv())
.info_object_cache <- new.env(parent = emptyenv())

.clear_fit_object_cache <- function(names = NULL) {

  if (is.null(names)) {
    names <- ls(envir = .fit_object_cache)
  } else {
    names <- intersect(names, ls(envir = .fit_object_cache))
  }
  if (length(names) > 0L) {
    rm(list = names, envir = .fit_object_cache)
  }

  return(invisible(TRUE))
}

.clear_info_object_cache <- function(names = NULL) {

  if (is.null(names)) {
    names <- ls(envir = .info_object_cache)
  } else {
    names <- intersect(names, ls(envir = .info_object_cache))
  }
  if (length(names) > 0L) {
    rm(list = names, envir = .info_object_cache)
  }

  return(invisible(TRUE))
}

is_true_env <- function(name) {

  value <- Sys.getenv(name)
  value <- tolower(value)

  return(value %in% c("1", "true", "yes", "y"))
}

is_false_env <- function(name) {

  value <- Sys.getenv(name)
  value <- tolower(value)

  return(value %in% c("0", "false", "no", "n"))
}

validate_cached_fit <- function(name, fit = NULL, info = NULL,
                                metadata = NULL, check_source = TRUE,
                                check_files = TRUE, deep = FALSE) {

  messages <- character()
  paths    <- fit_cache_paths(name)
  entry    <- fit_catalog_entry(name)

  if (check_files && !file.exists(paths[["fit"]]) && is.null(fit)) {
    messages <- c(messages, "fit file is missing")
    return(messages)
  }
  if (check_files && !file.exists(paths[["fit"]])) {
    messages <- c(messages, "fit file is missing")
  }

  if (is.null(entry)) {
    messages <- c(messages, "fit is missing from fit_catalog()")
  }

  if (is.null(metadata)) {
    metadata <- load_fit_metadata(name)
  }
  if (is.null(metadata)) {
    messages <- c(messages, "metadata file is missing")
  }

  if (!is.null(entry) && !is.null(metadata)) {
    if (is.null(metadata[["version"]]) || !identical(metadata[["version"]], FIT_CACHE_VERSION)) {
      messages <- c(messages, "cache version changed")
    }
    if (!is.null(metadata[["name"]]) && !identical(metadata[["name"]], name)) {
      messages <- c(messages, "metadata name mismatch")
    }
    if (is.null(metadata[["fit_class"]]) || !entry[["class"]] %in% metadata[["fit_class"]]) {
      messages <- c(messages, paste0("metadata fit class does not include '", entry[["class"]], "'"))
    }

    if (isTRUE(entry[["has_loo"]]) && !isTRUE(metadata[["has_loo"]])) {
      messages <- c(messages, "metadata reports missing LOO")
    }
    if (isFALSE(entry[["has_loo"]]) && isTRUE(metadata[["has_loo"]])) {
      messages <- c(messages, "metadata reports unexpected LOO")
    }

    if (isTRUE(entry[["has_waic"]]) && !isTRUE(metadata[["has_waic"]])) {
      messages <- c(messages, "metadata reports missing WAIC")
    }
    if (isFALSE(entry[["has_waic"]]) && isTRUE(metadata[["has_waic"]])) {
      messages <- c(messages, "metadata reports unexpected WAIC")
    }

    if (isTRUE(entry[["has_marglik"]]) && !isTRUE(metadata[["has_marglik"]])) {
      messages <- c(messages, "metadata reports missing marginal likelihood")
    }
    if (isFALSE(entry[["has_marglik"]]) && isTRUE(metadata[["has_marglik"]])) {
      messages <- c(messages, "metadata reports unexpected marginal likelihood")
    }

    if (isTRUE(entry[["has_metafor"]])) {
      if (!isTRUE(metadata[["has_metafor_info"]])) {
        messages <- c(messages, "metadata reports missing metafor reference")
      }
      if (check_files && !file.exists(paths[["info"]]) && is.null(info)) {
        messages <- c(messages, "metafor info file is missing")
      }
    }

    expected_md5 <- source_file_md5(entry[["source_file"]])
    if (check_source && !is.na(expected_md5) &&
        (is.null(metadata[["source_file_md5"]]) || !identical(metadata[["source_file_md5"]], expected_md5))) {
      messages <- c(messages, "source file hash changed")
    }

    expected_package_md5 <- package_source_md5()
    if (check_source && !is.na(expected_package_md5) &&
        (is.null(metadata[["package_source_md5"]]) ||
         !identical(metadata[["package_source_md5"]], expected_package_md5))) {
      messages <- c(messages, "cache source hash changed")
    }
  }

  if (!deep && is.null(fit) && is.null(info)) {
    return(messages)
  }

  if (is.null(fit)) {
    fit <- read_cached_fit(name)
  }
  if (is.null(fit)) {
    messages <- c(messages, "fit file cannot be read")
    return(messages)
  }

  if (is.null(info)) {
    info <- read_cached_info(name)
  }

  if (!is.null(entry) && !inherits(fit, entry[["class"]])) {
    messages <- c(messages, paste0("fit does not inherit from '", entry[["class"]], "'"))
  }

  if (!is.null(entry) && isTRUE(entry[["has_loo"]]) && is.null(fit[["loo"]])) {
    messages <- c(messages, "LOO is missing")
  }
  if (!is.null(entry) && isFALSE(entry[["has_loo"]]) && !is.null(fit[["loo"]])) {
    messages <- c(messages, "LOO is unexpectedly present")
  }

  if (!is.null(entry) && isTRUE(entry[["has_waic"]]) && is.null(fit[["waic"]])) {
    messages <- c(messages, "WAIC is missing")
  }
  if (!is.null(entry) && isFALSE(entry[["has_waic"]]) && !is.null(fit[["waic"]])) {
    messages <- c(messages, "WAIC is unexpectedly present")
  }

  if (!is.null(entry) && isTRUE(entry[["has_marglik"]]) && is.null(fit[["marglik"]])) {
    messages <- c(messages, "marginal likelihood is missing")
  }
  if (!is.null(entry) && isFALSE(entry[["has_marglik"]]) && !is.null(fit[["marglik"]])) {
    messages <- c(messages, "marginal likelihood is unexpectedly present")
  }

  if (!is.null(entry) && isTRUE(entry[["has_metafor"]])) {
    if (is.null(info) || !"metafor" %in% names(info) || is.null(info[["metafor"]])) {
      messages <- c(messages, "metafor reference is missing")
    }
  }

  return(messages)
}

is_cached_fit_valid <- function(name, check_source = TRUE, deep = FALSE) {

  return(length(validate_cached_fit(name, check_source = check_source, deep = deep)) == 0)
}

# Skip if pre-fitted models are not available
skip_if_no_fits <- function() {

  if (length(list_fits(validate = TRUE)) == 0) {
    skip("No valid cached fits available. Run `devtools::test(filter = '01-')` first.")
  }
}

skip_if_missing_fits <- function(names) {

  missing <- setdiff(names, list_fits(validate = TRUE))
  if (length(missing) > 0) {
    skip(paste0("Required pre-fitted models missing or stale: ", paste(missing, collapse = ", ")))
  }
}

skip_if_not_full_diagnostics <- function(reason) {

  if (!is_true_env("ROBMA_TEST_FULL_DIAGNOSTICS")) {
    skip(paste(
      "Skipping extended diagnostic redundancy check by default.",
      reason,
      "Set ROBMA_TEST_FULL_DIAGNOSTICS=TRUE to run it."
    ))
  }
}

# Skip model fitting if a valid cached fit exists. Refit by setting
# ROBMA_TEST_FORCE_REFIT=TRUE or ROBMA_TEST_SKIP_REFIT=FALSE.
skip_refit_if_cached <- function(name) {

  if (is_true_env("ROBMA_TEST_FORCE_REFIT")) {
    return(invisible(FALSE))
  }

  skip_refit_env <- Sys.getenv("ROBMA_TEST_SKIP_REFIT")
  skip_refit     <- if (skip_refit_env == "") TRUE else !is_false_env("ROBMA_TEST_SKIP_REFIT")
  fit_names      <- catalog_group_fits(name)

  if (skip_refit && length(fit_names) > 0 &&
      all(vapply(fit_names, is_cached_fit_valid, TRUE))) {
    skip("Skipping model refit: valid cached fit exists.")
  }

  return(invisible(FALSE))
}

# ============================================================================ #
# HELPER FUNCTIONS: Model Fit Saving / Loading
# ============================================================================ #

save_fit <- function(name, fit, info = NULL) {
  metadata <- fit_cache_metadata(name, fit, info)
  problems <- validate_cached_fit(
    name         = name,
    fit          = fit,
    info         = info,
    metadata     = metadata,
    check_source = FALSE
  )
  problems <- problems[problems != "fit file is missing"]

  if (length(problems) > 0) {
    stop(
      paste0("Cached fit '", name, "' is invalid: ", paste(problems, collapse = "; ")),
      call. = FALSE
    )
  }

  paths <- fit_cache_paths(name)

  # Save model fit (marglik and loo are stored in the fit object)
  saveRDS(fit, file = paths[["fit"]])

  # Save info if provided
  if (!is.null(info)) {
    saveRDS(info, file = paths[["info"]])
  } else if (file.exists(paths[["info"]])) {
    file.remove(paths[["info"]])
  }

  # Save cache metadata last so interrupted fits do not validate later.
  saveRDS(metadata, file = paths[["metadata"]])

  assign(name, fit, envir = .fit_object_cache)
  if (!is.null(info)) {
    assign(name, info, envir = .info_object_cache)
  } else {
    .clear_info_object_cache(name)
  }

  return(invisible(TRUE))
}

load_fit <- function(name, validate = TRUE) {

  if (validate) {
    problems <- validate_cached_fit(name)
    if (length(problems) > 0) {
      skip(paste0(
        "Cached fit '", name, "' is missing, stale, or invalid: ",
        paste(problems, collapse = "; "),
        ". Run `test(filter = 'test-01')`."
      ))
    }
  }

  if (exists(name, envir = .fit_object_cache, inherits = FALSE)) {
    return(get(name, envir = .fit_object_cache, inherits = FALSE))
  }

  fit <- read_cached_fit(name)
  if (is.null(fit)) {
    stop(paste0("Cached fit '", name, "' could not be read."), call. = FALSE)
  }

  assign(name, fit, envir = .fit_object_cache)

  return(fit)
}

load_info <- function(name, validate = TRUE) {

  if (validate) {
    problems <- validate_cached_fit(name)
    if (length(problems) > 0) {
      skip(paste0(
        "Cached fit '", name, "' is missing, stale, or invalid: ",
        paste(problems, collapse = "; "),
        ". Run `test(filter = 'test-01')`."
      ))
    }
  }

  if (exists(name, envir = .info_object_cache, inherits = FALSE)) {
    return(get(name, envir = .info_object_cache, inherits = FALSE))
  }

  info <- read_cached_info(name)

  assign(name, info, envir = .info_object_cache)

  return(info)
}

list_fits <- function(name, feature, class, family, has_metafor, has_loo,
                      has_waic, has_marglik, tier, validate = TRUE,
                      deep = FALSE) {

  files <- suppressWarnings(list.files(temp_fits_dir, pattern = "\\.RDS$"))
  files <- sub("\\.RDS$", "", files)

  if (!missing(name)) {
    files <- intersect(files, name)
  }

  if (!missing(feature) || !missing(class) || !missing(family) ||
      !missing(has_metafor) || !missing(has_loo) || !missing(has_waic) ||
      !missing(has_marglik) || !missing(tier)) {
    selected <- catalog_fits(
      feature     = if (!missing(feature)) feature else NULL,
      class       = if (!missing(class)) class else NULL,
      family      = if (!missing(family)) family else NULL,
      has_metafor = if (!missing(has_metafor)) has_metafor else NULL,
      has_loo     = if (!missing(has_loo)) has_loo else NULL,
      has_waic    = if (!missing(has_waic)) has_waic else NULL,
      has_marglik = if (!missing(has_marglik)) has_marglik else NULL,
      tier        = if (!missing(tier)) tier else NULL
    )
    files <- intersect(files, selected)
  }

  if (validate) {
    files <- files[vapply(files, is_cached_fit_valid, TRUE, deep = deep)]
  }

  return(files)
}

lazy_cached_objects <- function(names, loader) {

  names <- unique(names)

  return(structure(
    list(
      .names  = names,
      .loader = loader,
      .cache  = new.env(parent = emptyenv())
    ),
    class = "lazy_cached_objects"
  ))
}

lazy_fits <- function(names = list_fits(), validate = TRUE) {

  if (validate) {
    skip_if_missing_fits(names)
  }

  return(lazy_cached_objects(
    names,
    function(name) load_fit(name, validate = FALSE)
  ))
}

lazy_infos <- function(names = list_fits(), validate = TRUE) {

  if (validate) {
    skip_if_missing_fits(names)
  }

  return(lazy_cached_objects(
    names,
    function(name) load_info(name, validate = FALSE)
  ))
}

load_fits <- function(names = list_fits(), validate = TRUE) {

  out <- lapply(names, load_fit, validate = validate)
  names(out) <- names

  return(out)
}

load_infos <- function(names = list_fits(), validate = TRUE) {

  out <- lapply(names, load_info, validate = validate)
  names(out) <- names

  return(out)
}

names.lazy_cached_objects <- function(x) {

  return(unclass(x)[[".names"]])
}

length.lazy_cached_objects <- function(x) {

  return(length(names(x)))
}

`[[.lazy_cached_objects` <- function(x, i, exact = TRUE) {

  object_names <- names(x)

  if (is.numeric(i)) {
    i <- object_names[[i]]
  }
  if (!is.character(i) || length(i) != 1L || !i %in% object_names) {
    return(NULL)
  }

  cache <- unclass(x)[[".cache"]]
  if (!exists(i, envir = cache, inherits = FALSE)) {
    assign(i, unclass(x)[[".loader"]](i), envir = cache)
  }

  return(get(i, envir = cache, inherits = FALSE))
}

`[.lazy_cached_objects` <- function(x, i) {

  object_names <- names(x)
  selected     <- object_names[i]

  return(lazy_cached_objects(
    selected,
    unclass(x)[[".loader"]]
  ))
}

as.list.lazy_cached_objects <- function(x, ...) {

  out <- lapply(names(x), function(name) x[[name]])
  names(out) <- names(x)

  return(out)
}

clean_cached_fits <- function(name) {

  if (!missing(name)) {
    fit_names <- catalog_group_fits(name)
    if (length(fit_names) == 0) {
      fit_names <- name
    }
    for (fit_name in fit_names) {
      paths <- unlist(fit_cache_paths(fit_name), use.names = FALSE)
      paths <- paths[file.exists(paths)]
      if (length(paths) > 0L) {
        file.remove(paths)
      }
    }
    .clear_fit_object_cache(fit_names)
    .clear_info_object_cache(fit_names)
  } else {
    # Remove all cached files from test directories
    unlink(temp_fits_dir, recursive = TRUE)
    unlink(temp_info_dir, recursive = TRUE)
    unlink(temp_metadata_dir, recursive = TRUE)
    unlink(temp_temp_dir, recursive = TRUE)

    # Recreate empty directories
    dir.create(temp_fits_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_info_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_metadata_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_temp_dir, showWarnings = FALSE, recursive = TRUE)

    .clear_fit_object_cache()
    .clear_info_object_cache()
  }

  message("Cleaned cached fits in: ", test_files_dir)

  return(invisible(TRUE))
}


# Construct conditional standardized residuals for multilevel metafor::rma.mv fits.
# metafor::rstandard() only provides marginal residuals for rma.mv, so the
# conditional oracle must be built from the marginal fit, BLUP random effects,
# and the GLS residual variance matrix.
metafor_rstandard_conditional_mv <- function(model) {
  if (!inherits(model, "rma.mv")) {
    stop("'model' must inherit from 'rma.mv'.", call. = FALSE)
  }
  if (is.null(model[["yi"]]) || is.null(model[["X"]]) ||
      is.null(model[["M"]])  || is.null(model[["V"]])) {
    stop("Model object does not contain the information needed for conditional residuals.", call. = FALSE)
  }

  fitted_values <- as.vector(stats::fitted(model))

  ranef_components <- metafor::ranef(model, expand = TRUE)
  if (length(ranef_components) == 0L) {
    stop("No random-effect BLUPs available for this 'rma.mv' object.", call. = FALSE)
  }

  random_effects <- Reduce(
    `+`,
    lapply(ranef_components, function(component) component[["intrcpt"]])
  )

  residuals <- as.vector(model[["yi"]] - (fitted_values + random_effects))

  hat_matrix <- stats::hatvalues(model, type = "matrix")
  marginal_v <- model[["M"]]
  sampling_v <- model[["V"]]
  weights    <- chol2inv(chol(marginal_v))
  k          <- model[["k"]]
  identity   <- diag(k)
  imh        <- identity - hat_matrix

  residual_v <- sampling_v %*% weights %*% imh %*% marginal_v %*%
    t(imh) %*% weights %*% sampling_v
  se <- sqrt(diag(residual_v))

  out <- data.frame(
    resid = residuals,
    se    = as.vector(se),
    z     = residuals / se
  )
  rownames(out) <- NULL

  return(out)
}
