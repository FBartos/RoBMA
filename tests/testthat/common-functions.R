FIT_CACHE_SCHEMA_VERSION <- 7L
FIT_CACHE_METADATA_FIELDS <- c(
  "schema_version",
  "name",
  "saved_at",
  "robma_version",
  "bayestools_version",
  "r_version",
  "jags_version",
  "fit_class",
  "source_file",
  "source_file_md5",
  "package_source_md5",
  "bayestools_backend_fingerprint",
  "has_loo",
  "has_waic",
  "has_marglik",
  "has_metafor_info"
)

TEST_PROFILES <- c("standard", "certification")

test_profile <- function() {

  profile <- Sys.getenv("ROBMA_TEST_PROFILE", unset = "standard")
  profile <- tolower(profile)

  if (!profile %in% TEST_PROFILES) {
    stop(
      "'ROBMA_TEST_PROFILE' must be one of: ",
      paste(TEST_PROFILES, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  return(profile)
}

is_certification_profile <- function() {

  return(identical(test_profile(), "certification"))
}

.quiet_llm_reporter_class <- R6::R6Class(
  classname = "RoBMAQuietLlmReporter",
  inherit   = testthat::LlmReporter,
  public    = list(
    add_result = function(context, test, result) {

      if (self$is_full()) {
        return(invisible())
      }
      if (inherits(result, "expectation_skip")) {
        self$n_skip <- self$n_skip + 1L
        return(invisible())
      }

      super$add_result(context, test, result)
    }
  )
)

quiet_llm_reporter <- function(...) {

  return(.quiet_llm_reporter_class$new(...))
}

test_glmm_fit_settings <- function() {

  if (is_certification_profile()) {
    return(list(
      chains = 3L,
      sample = 5000L,
      burnin = 2000L,
      adapt  = 500L
    ))
  }

  return(list(
    chains = 3L,
    sample = 1500L,
    burnin = 500L,
    adapt  = 500L
  ))
}

.test_profile_cache_dir <- function(root, profile = test_profile()) {

  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  return(file.path(root, profile))
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

# Get the directory where prefitted models are stored. An explicit directory
# is a shared cache root; each profile always receives its own subdirectory.
test_files_root <- Sys.getenv("ROBMA_TEST_FILES_DIR")
if (test_files_root == "") {
  test_files_root <- file.path(.common_functions_dir(), "test_files")
}
test_files_dir <- .test_profile_cache_dir(test_files_root)

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

# Use skip_if_no_fits() for tests that need pre-fitted models.

# ============================================================================ #
# HELPER FUNCTIONS: Reference File Testing
# ============================================================================ #

.reference_candidate_path <- function(path) {

  if (grepl("[.]txt$", path)) {
    return(sub("[.]txt$", ".new.txt", path))
  }

  return(paste0(path, ".new.txt"))
}


.write_reference_candidate <- function(lines, path) {

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(lines, connection, useBytes = TRUE)

  return(invisible(path))
}


.test_reference_lines <- function(actual_output, filename, info_msg,
                                  print_dir) {

  ref_file      <- file.path(print_dir, filename)
  candidate     <- .reference_candidate_path(ref_file)
  candidate_msg <- paste0(
    "Reference candidate cached at ", candidate,
    ". Review it with 'review_test_snapshots()'."
  )

  if (!file.exists(ref_file)) {
    skip(paste("Reference file", filename, "not found."))
    return(invisible(actual_output))
  }

  expected_output <- readLines(ref_file, warn = FALSE, encoding = "UTF-8")
  if (identical(actual_output, expected_output)) {
    unlink(candidate)
  } else {
    .write_reference_candidate(actual_output, candidate)
  }

  expect_equal(
    actual_output,
    expected_output,
    info = paste(c(info_msg, candidate_msg), collapse = "\n")
  )
  return(invisible(actual_output))
}


test_reference_table <- function(table, filename, info_msg = NULL,
                                 print_dir = REFERENCE_DIR) {

  actual_output <- capture_output_lines(table, print = TRUE, width = 150)
  return(.test_reference_lines(
    actual_output = actual_output,
    filename      = filename,
    info_msg      = info_msg,
    print_dir     = print_dir
  ))
}


test_reference_text <- function(text, filename, info_msg = NULL,
                                print_dir = REFERENCE_DIR) {

  actual_output <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  return(.test_reference_lines(
    actual_output = actual_output,
    filename      = filename,
    info_msg      = info_msg,
    print_dir     = print_dir
  ))
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

hypothesis_plan_fit_names <- function() {

  return(c(
    treatment = "hypothesis_plan_treatment",
    meandif = "hypothesis_plan_meandif",
    orthonormal = "hypothesis_plan_orthonormal",
    independent = "hypothesis_plan_independent",
    ordered_2 = "hypothesis_plan_ordered_2",
    ordered_3 = "hypothesis_plan_ordered_3",
    ordered_4 = "hypothesis_plan_ordered_4",
    scaled_interaction = "hypothesis_plan_scaled_interaction",
    location_scale = "hypothesis_plan_location_scale",
    model_averaged = "hypothesis_plan_model_averaged",
    robma_mixture = "hypothesis_plan_robma_mixture",
    two_scale = "hypothesis_plan_two_scale",
    meandif_mt = "hypothesis_plan_meandif_mt"
  ))
}


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
      "dat.lehmann2018_RoBMA_3lvl_mods_scale",
      "brma.mv_latent",
      "brma.mv_whitened",
      "brma.mv_block_mvn",
      "brma.mv_block_mvn_fixed_random_null",
      "brma.mv_block_mvn_random",
      "brma.mv_block_mvn_random_sampled",
      "brma.mv_block_mvn_known_R",
      "brma.mv_latent_estimate_scale",
      "brma.mv_block_mvn_estimate_scale",
      "brma.mv_block_mvn_random_scale",
      "brma.mv_block_mvn_3lvl_scale_total",
      "brma.mv_block_mvn_3lvl_scale_top",
      "brma.mv_block_mvn_3lvl_scale_bottom",
      "brma.mv_block_mvn_mods",
      "brma.mv_block_mvn_random_mods_scale",
      "brma.mv_v14_konstantopoulos2011_cs",
      "brma.mv_v14_assink2016_nested",
      "brma.mv_v14_ishak2007_har",
      "brma.mv_v14_begg1989_study_treatment",
      "brma.mv_singular_regularized_whitened",
      "brma.mv_singular_regularized_block_mvn"
    ),
    class = c(
      rep("brma.norm", 11),
      rep("brma.glmm", 4),
      rep("bPET", 3),
      rep("bPEESE", 3),
      rep("bselmodel", 4),
      rep("BMA.norm", 4),
      rep("BMA.glmm", 4),
      rep("RoBMA", 5),
      rep("brma.mv", 21)
    ),
    family = c(
      rep("norm", 11),
      rep("glmm", 4),
      rep("norm", 3),
      rep("norm", 3),
      rep("norm", 4),
      rep("norm", 4),
      rep("glmm", 4),
      rep("norm", 5),
      rep("norm", 21)
    ),
    source_file = c(
      rep("test-01-brma.norm.R", 11),
      rep("test-01-brma.glmm.R", 4),
      rep("test-01-bPET.R", 3),
      rep("test-01-bPEESE.R", 3),
      rep("test-01-bselmodel.R", 4),
      rep("test-01-BMA.norm.R", 4),
      rep("test-01-BMA.glmm.R", 4),
      rep("test-01-RoBMA.R", 5),
      rep("test-01-brma.mv.R", 21)
    ),
    has_metafor = c(
      TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, FALSE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE,
      rep(FALSE, 16),
      TRUE,
      FALSE,
      FALSE,
      TRUE,
      rep(FALSE, 8),
      rep(TRUE, 4),
      rep(FALSE, 2)
    ),
    has_waic = FALSE,
    tier = c(
      "core", "core", "extended", "extended", "core", "extended", "extended", "extended",
      "core", "core", "extended",
      "core", "core", "core", "extended",
      "core", "extended", "extended",
      "core", "extended", "core",
      "core", "extended", "extended", "extended",
      "core", "extended", "core", "extended",
      "core", "core", "extended", "extended",
      "core", "extended", "core", "extended", "extended",
      "core", "core", "core", "core", "core", "core", "core", "core", "core", "core",
      "core", "core", "core", "core", "core", "extended", "extended", "extended", "extended", "core",
      "core"
    ),
    stringsAsFactors = FALSE
  )
  catalog[["has_loo"]] <- !catalog[["name"]] %in% c(
    "brma.mv_block_mvn_random_mods_scale"
  )
  catalog[["has_marglik"]] <- !catalog[["class"]] %in%
    c("BMA.norm", "BMA.glmm", "RoBMA", "brma.mv")
  catalog[["has_marglik"]][catalog[["name"]] %in% c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment"
  )] <- TRUE

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
    c("RoBMA", "normal", "multilevel", "mods", "scale"),
    c("brma.mv", "normal", "known_v", "latent"),
    c("brma.mv", "normal", "known_v", "whitened"),
    c("brma.mv", "normal", "known_v", "block_mvn"),
    c("brma.mv", "metafor", "normal", "known_v", "block_mvn", "fixed", "random_null"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "sampled"),
    c("brma.mv", "metafor", "normal", "known_v", "known_R", "block_mvn", "random", "mods"),
    c("brma.mv", "normal", "known_v", "latent", "random", "scale"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "scale"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "scale", "allocation"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "scale", "3lvl", "total_sd"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "scale", "3lvl", "top_sd"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "scale", "3lvl", "bottom_sd"),
    c("brma.mv", "normal", "known_v", "block_mvn", "mods"),
    c("brma.mv", "normal", "known_v", "block_mvn", "random", "mods", "scale"),
    c("brma.mv", "metafor", "normal", "known_v", "v14", "multilevel", "cs"),
    c("brma.mv", "metafor", "normal", "known_v", "v14", "multilevel", "nested"),
    c("brma.mv", "metafor", "normal", "known_v", "v14", "multilevel", "mods", "har"),
    c("brma.mv", "metafor", "normal", "known_v", "v14", "multilevel", "mods", "cs"),
    c("brma.mv", "normal", "known_v", "singular", "whitened"),
    c("brma.mv", "normal", "known_v", "singular", "block_mvn")
  ))

  bma_mv_catalog <- data.frame(
    name          = "BMA.mv_random_components",
    class         = "BMA.mv",
    family        = "norm",
    source_file   = "test-01-BMA.mv.R",
    has_metafor   = FALSE,
    has_waic      = FALSE,
    tier          = "core",
    has_loo       = TRUE,
    has_marglik   = FALSE,
    features      = I(list(c(
      "BMA.mv", "normal", "known_v", "random", "mods",
      "product_space", "random_inclusion"
    ))),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, bma_mv_catalog)

  robma_mv_catalog <- data.frame(
    name          = "RoBMA.mv_marg_product_space",
    class         = "RoBMA.mv",
    family        = "norm",
    source_file   = "test-01-RoBMA.mv.R",
    has_metafor   = FALSE,
    has_waic      = FALSE,
    tier          = "core",
    has_loo       = TRUE,
    has_marglik   = FALSE,
    features      = I(list(c(
      "RoBMA.mv", "normal", "known_v", "random", "mods",
      "product_space", "random_inclusion", "selection", "marg", "PET",
      "PEESE"
    ))),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, robma_mv_catalog)

  bselmodel_mv_catalog <- data.frame(
    name          = "bselmodel.mv_marg_random",
    class         = "bselmodel.mv",
    family        = "norm",
    source_file   = "test-01-bselmodel.mv.R",
    has_metafor   = FALSE,
    has_waic      = FALSE,
    tier          = "core",
    has_loo       = TRUE,
    has_marglik   = TRUE,
    features      = I(list(c(
      "bselmodel.mv", "normal", "known_v", "random", "mods",
      "selection", "marg"
    ))),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, bselmodel_mv_catalog)

  pet_peese_mv_catalog <- data.frame(
    name = c("bPET.mv_random", "bPEESE.mv_random"),
    class = c("bPET.mv", "bPEESE.mv"),
    family = c("norm", "norm"),
    source_file = rep("test-01-PET-PEESE.mv.R", 2L),
    has_metafor = FALSE,
    has_waic = FALSE,
    tier = "core",
    has_loo = TRUE,
    has_marglik = TRUE,
    features = I(list(
      c("bPET.mv", "normal", "known_v", "random", "mods", "PET"),
      c("bPEESE.mv", "normal", "known_v", "random", "mods", "PEESE")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, pet_peese_mv_catalog)

  parity_catalog <- data.frame(
    name = c(
      "vif_parity_brma",
      "vif_parity_brma_mv"
    ),
    class = c("brma.norm", "brma.mv"),
    family = c("norm", "norm"),
    source_file = rep("test-01-vif-parity.R", 2),
    has_metafor = FALSE,
    has_waic = FALSE,
    tier = "core",
    has_loo = TRUE,
    has_marglik = c(TRUE, FALSE),
    features = I(list(
      c("normal", "mods", "vif", "parity"),
      c("brma.mv", "normal", "known_v", "mods", "vif", "parity")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, parity_catalog)

  iwmde_oracle_catalog <- data.frame(
    name = c(
      "nielweise2008_glmm_effect_null",
      "dat.lehmann2018-3PSM_effect_null",
      "iwmde_known_v_tau_full",
      "iwmde_known_v_tau_null"
    ),
    class = c("brma.glmm", "bselmodel", "brma.mv", "brma.mv"),
    family = c("glmm", "norm", "norm", "norm"),
    source_file = rep("test-01-iwmde-oracle-nested.R", 4),
    has_metafor = FALSE,
    has_waic = FALSE,
    tier = "extended",
    has_loo = c(FALSE, FALSE, TRUE, TRUE),
    has_marglik = TRUE,
    features = I(list(
      c("glmm", "poisson", "iwmde", "nested_null"),
      c("normal", "selection", "iwmde", "nested_null"),
      c("normal", "known_v", "iwmde", "nested_full"),
      c("normal", "known_v", "iwmde", "nested_null")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, iwmde_oracle_catalog)

  iwmde_multilevel_catalog <- data.frame(
    name          = "dat.lehmann2018-3PSM_3lvl_mods_marginal",
    class         = "bselmodel",
    family        = "norm",
    source_file   = "test-01-iwmde-multilevel-selection.R",
    has_metafor   = FALSE,
    has_waic      = FALSE,
    tier          = "extended",
    has_loo       = TRUE,
    has_marglik   = TRUE,
    features      = I(list(c(
      "normal", "selection", "multilevel", "mods", "iwmde",
      "marginal_selection"
    ))),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, iwmde_multilevel_catalog)

  fixed_model_catalog <- data.frame(
    name = c(
      "fixed_null_brma",
      "fixed_nonzero_brma",
      "fixed_nonzero_brma_mv",
      "fixed_null_brma_glmm",
      "fixed_null_brma_glmm_pois",
      "fixed_null_bPET",
      "fixed_null_bPEESE",
      "fixed_null_bselmodel"
    ),
    class = c(
      "brma.norm",
      "brma.norm",
      "brma.mv",
      "brma.glmm",
      "brma.glmm",
      "bPET",
      "bPEESE",
      "bselmodel"
    ),
    family = c("norm", "norm", "norm", "glmm", "glmm", "norm", "norm", "norm"),
    source_file = rep("test-01-fixed-null-models.R", 8),
    has_metafor = FALSE,
    has_waic = FALSE,
    tier = "core",
    has_loo = TRUE,
    has_marglik = TRUE,
    features = I(list(
      c("normal", "fixed", "null", "exact_marglik"),
      c("normal", "fixed", "nonzero_tau", "exact_marglik"),
      c("brma.mv", "normal", "known_v", "fixed", "nonzero_tau", "exact_marglik"),
      c("glmm", "binomial", "fixed", "null", "point_nuisance"),
      c("glmm", "poisson", "fixed", "null", "point_nuisance"),
      c("normal", "PET", "fixed", "null", "exact_marglik"),
      c("normal", "PEESE", "fixed", "null", "exact_marglik"),
      c("normal", "selection", "fixed", "null", "exact_marglik")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, fixed_model_catalog)

  catalog[["profile"]] <- "certification"
  catalog[["profile"]][catalog[["name"]] %in% c(
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
    "dat.lehmann2018-PET",
    "dat.lehmann2018-PEESE",
    "dat.lehmann2018-3PSM",
    "dat.lehmann2018_BMA.norm",
    "dat.lehmann2018_BMA.norm_mods",
    "bcg_BMA.glmm",
    "dat.lehmann2018_RoBMA",
    "dat.lehmann2018_RoBMA_mods",
    "dat.lehmann2018_RoBMA_mods2",
    "dat.lehmann2018_RoBMA_3lvl_mods_scale",
    "BMA.mv_random_components",
    "RoBMA.mv_marg_product_space",
    "bselmodel.mv_marg_random",
    "bPET.mv_random",
    "bPEESE.mv_random",
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_fixed_random_null",
    "brma.mv_block_mvn_random",
    "vif_parity_brma",
    "vif_parity_brma_mv",
    "fixed_null_brma",
    "fixed_nonzero_brma",
    "fixed_nonzero_brma_mv",
    "fixed_null_brma_glmm",
    "fixed_null_brma_glmm_pois",
    "fixed_null_bPET",
    "fixed_null_bPEESE",
    "fixed_null_bselmodel"
  )] <- "standard"

  hypothesis_plan_catalog <- data.frame(
    name = unname(hypothesis_plan_fit_names()),
    class = c(
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "brma.norm",
      "BMA.norm",
      "RoBMA",
      "brma.mv",
      "brma.norm"
    ),
    family = "norm",
    source_file = "test-01-hypothesis-plan.R",
    has_metafor = FALSE,
    has_loo = FALSE,
    has_waic = FALSE,
    has_marglik = FALSE,
    tier = "core",
    features = I(list(
      c("hypothesis_plan", "normal", "factor_mods"),
      c("hypothesis_plan", "normal", "factor_mods"),
      c("hypothesis_plan", "normal", "factor_mods"),
      c("hypothesis_plan", "normal", "factor_mods"),
      c("hypothesis_plan", "normal", "factor_mods", "ordered"),
      c("hypothesis_plan", "normal", "factor_mods", "ordered"),
      c("hypothesis_plan", "normal", "factor_mods", "ordered"),
      c("hypothesis_plan", "normal", "factor_mods", "interaction"),
      c("hypothesis_plan", "normal", "factor_mods", "scale"),
      c("hypothesis_plan", "normal", "factor_mods", "product_space"),
      c("hypothesis_plan", "normal", "factor_mods", "product_space", "selection"),
      c("hypothesis_plan", "normal", "factor_mods", "known_v", "multilevel", "scale"),
      c("hypothesis_plan", "normal", "factor_mods", "multivariate_t")
    )),
    profile = "standard",
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, hypothesis_plan_catalog)
  expression_catalog <- hypothesis_plan_catalog[1L, , drop = FALSE]
  expression_catalog$name <- "marginal_means_expression"
  expression_catalog$features <- I(list(c("normal", "factor_mods", "ordered", "expression_total", "marginal_means")))
  catalog <- rbind(catalog, expression_catalog)
  zero_ordinate_catalog <- hypothesis_plan_catalog[1L, , drop = FALSE]
  zero_ordinate_catalog$name <- "marginal_means_zero_ordinate"
  zero_ordinate_catalog$source_file <- "test-01-marginal-means-zero-ordinate.R"
  zero_ordinate_catalog$features <- I(list(c("normal", "continuous_mods", "truncated_effect", "marginal_means")))
  catalog <- rbind(catalog, zero_ordinate_catalog)

  ordinary_lane_catalog <- data.frame(
    name = c(
      "lane_carrier_treatment",
      "lane_carrier_meandif",
      "lane_carrier_nonzero",
      "lane_carrier_averaged",
      "lane_carrier_ordinary",
      "lane_carrier_gated",
      "lane_contribution_constant2",
      "lane_contribution_constant1",
      "lane_contribution_tau",
      "lane_contribution_point",
      "lane_contribution_discovery",
      "lane_exp_affine",
      "lane_factor_averaged_reference",
      "lane_factor_ordered_labels",
      "lane_factor_numeric_labels",
      "lane_scalar_standardized",
      "lane_scalar_raw",
      "lane_scalar_convolution",
      "lane_scalar_averaged",
      "lane_scalar_averaged_default",
      "lane_known_v_factor_scale",
      "lane_random_correlation",
      "lane_random_correlation_direct",
      "lane_random_correlation_gated",
      "lane_nested_total",
      "lane_selection_reference",
      "lane_glmm_point_baserate",
      "lane_glmm_boundary_beta",
      "lane_glmm_truncated_beta",
      "lane_glmm_truncated_lograte"
    ),
    class = c("brma.norm", "brma.norm", "brma.norm", "BMA.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "BMA.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "brma.norm", "BMA.norm", "BMA.norm", "brma.mv", "brma.mv", "brma.mv", "BMA.mv", "BMA.mv", "bselmodel", "brma.glmm", "brma.glmm", "brma.glmm", "brma.glmm"),
    family = c("norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm", "glmm", "glmm", "glmm", "glmm"),
    source_file = "test-01-ordinary-lane-fixtures.R",
    has_metafor = FALSE, has_loo = FALSE, has_waic = FALSE, has_marglik = FALSE,
    tier = "core", profile = "standard",
    features = I(list(
      c("ordinary_lane", "factor_mods", "carrier"),
      c("ordinary_lane", "factor_mods", "carrier"),
      c("ordinary_lane", "factor_mods", "carrier"),
      c("ordinary_lane", "factor_mods", "carrier"),
      c("ordinary_lane", "factor_mods", "carrier"),
      c("ordinary_lane", "factor_mods", "carrier"),
      c("ordinary_lane", "factor_mods", "formula_contribution"),
      c("ordinary_lane", "factor_mods", "formula_contribution"),
      c("ordinary_lane", "factor_mods", "formula_contribution"),
      c("ordinary_lane", "formula_contribution", "point"),
      c("ordinary_lane", "formula_contribution", "interaction"),
      c("ordinary_lane", "mods", "scale"),
      c("ordinary_lane", "factor_mods", "product_space"),
      c("ordinary_lane", "factor_mods", "ordered"),
      c("ordinary_lane", "factor_mods", "ordered"),
      c("ordinary_lane", "scalar_eligibility", "mods"),
      c("ordinary_lane", "scalar_eligibility", "mods"),
      c("ordinary_lane", "scalar_eligibility", "mods"),
      c("ordinary_lane", "scalar_eligibility", "mods"),
      c("ordinary_lane", "scalar_eligibility", "mods"),
      c("ordinary_lane", "known_v", "scale", "factor_mods"),
      c("ordinary_lane", "known_v", "random_correlation"),
      c("ordinary_lane", "known_v", "random_correlation", "direct_sd"),
      c("ordinary_lane", "known_v", "random_correlation", "product_space"),
      c("ordinary_lane", "known_v", "nested", "product_space"),
      c("ordinary_lane", "selection", "structural_monitor"),
      c("ordinary_lane", "binomial", "point_baserate"),
      c("ordinary_lane", "binomial", "boundary_beta"),
      c("ordinary_lane", "binomial", "truncated_nuisance"),
      c("ordinary_lane", "Poisson", "truncated_nuisance")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, ordinary_lane_catalog)

  # ensembles whose members differ in their bias terms (PET / PEESE designs)
  # and the fits of the random-effect prior plots
  bias_branch_catalog <- data.frame(
    name = c(
      "bias_branch_PP_mods",
      "bias_branch_PP_intercept",
      "bias_branch_bPET_mods",
      "bias_branch_known_v_ensemble"
    ),
    class = c("RoBMA", "RoBMA", "bPET", "RoBMA.mv"),
    family = "norm",
    source_file = "test-01-bias-branch-design.R",
    has_metafor = FALSE, has_loo = c(TRUE, FALSE, FALSE, FALSE),
    has_waic = FALSE, has_marglik = FALSE,
    tier = "core", profile = "standard",
    features = I(list(
      c("RoBMA", "normal", "mods", "PET", "PEESE", "bias_branch"),
      c("RoBMA", "normal", "PET", "PEESE", "bias_branch"),
      c("bPET", "normal", "mods", "PET", "bias_branch"),
      c("RoBMA.mv", "normal", "known_v", "mods", "PET", "PEESE", "bias_branch")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, bias_branch_catalog)

  plot_prior_random_catalog <- data.frame(
    name = c("plot_prior_random_allocated", "plot_prior_random_scale"),
    class = "brma.mv",
    family = "norm",
    source_file = "test-01-plot-prior-random.R",
    has_metafor = FALSE, has_loo = FALSE, has_waic = FALSE, has_marglik = FALSE,
    tier = "core", profile = "standard",
    features = I(list(
      c("brma.mv", "normal", "known_v", "random", "allocation", "plot_prior"),
      c("brma.mv", "normal", "known_v", "random", "scale", "plot_prior")
    )),
    stringsAsFactors = FALSE
  )
  catalog <- rbind(catalog, plot_prior_random_catalog)

  catalog[["dependencies"]] <- rep(list(data.frame(
    name = character(), artifact = character()
  )), nrow(catalog))
  for (name in c("iwmde_known_v_tau_full", "iwmde_known_v_tau_null")) {
    catalog[["dependencies"]][[match(name, catalog[["name"]])]] <- data.frame(
      name = c("brma.mv_block_mvn_random", "brma.mv_block_mvn_fixed_random_null"),
      artifact = c("info", "fit"), stringsAsFactors = FALSE
    )
  }

  return(catalog)
}

# `catalog` lets a caller that looks up many fits build the catalog once.
fit_catalog_entry <- function(name, catalog = fit_catalog()) {

  index <- match(name, catalog[["name"]])

  if (is.na(index)) {
    return(NULL)
  }

  return(catalog[index, , drop = FALSE])
}

active_fit_catalog <- function(catalog = fit_catalog()) {

  if (!is_certification_profile()) {
    catalog <- catalog[catalog[["profile"]] == "standard", , drop = FALSE]
  }

  active_fits <- Sys.getenv("ROBMA_TEST_ACTIVE_FITS", unset = "")
  if (nzchar(active_fits)) {
    requested <- if (identical(active_fits, "__none__")) {
      character()
    } else {
      strsplit(active_fits, ",", fixed = TRUE)[[1L]]
    }
    unknown <- setdiff(requested, catalog[["name"]])
    if (length(unknown) > 0L) {
      stop(
        "Unknown active cached fit: ", paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    catalog <- catalog[catalog[["name"]] %in% requested, , drop = FALSE]
  }

  return(catalog)
}

is_fit_active <- function(name) {

  return(name %in% active_fit_catalog()[["name"]])
}

skip_if_fit_not_active <- function(name) {

  if (!is_fit_active(name)) {
    testthat::skip(paste0(
      "Cached fit '", name, "' belongs to the certification profile."
    ))
  }

  return(invisible(FALSE))
}

catalog_group_fits <- function(name, active_only = TRUE) {

  catalog <- if (active_only) active_fit_catalog() else fit_catalog()

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
                         has_waic, has_marglik, tier,
                         catalog = fit_catalog()) {

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
                     "skip_if_not_installed", "skip_refit_if_cached",
                     "skip_if_fit_not_active", "skip_if_not_certification")
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

# The shared GLMM fit settings (chains, draws, burn-in, adaptation) are read at
# run time from test_glmm_fit_settings(); a generating file that references the
# helper hashes the value it returns for the active profile, because the call
# itself reads the same after a change of the settings.
.source_hash_fit_settings <- function(references) {

  if (!isTRUE(references)) {
    return(character())
  }

  return(paste0(
    "test_glmm_fit_settings() = ",
    paste(deparse(test_glmm_fit_settings(), width.cutoff = 500L), collapse = " ")
  ))
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
    references <- any(grepl("test_glmm_fit_settings", lines, fixed = TRUE))
    lines <- sub("#.*$", "", lines)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
  } else {
    references <- "test_glmm_fit_settings" %in% all.names(parsed)
    lines <- unlist(lapply(parsed, .source_hash_normalize_expr), use.names = FALSE)
    lines <- as.character(lines)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
  }
  lines <- c(lines, .source_hash_fit_settings(references))

  normalized <- tempfile("robma-fit-source-", fileext = ".R")
  writeLines(lines, normalized, useBytes = TRUE)
  on.exit(unlink(normalized), add = TRUE)

  return(unname(tools::md5sum(normalized)))
}

.package_source_md5_cache <- getOption("RoBMA.test.package_source_md5_cache")
if (!is.environment(.package_source_md5_cache)) {
  .package_source_md5_cache <- new.env(parent = emptyenv())
  options(RoBMA.test.package_source_md5_cache = .package_source_md5_cache)
}

.clear_package_source_md5_cache <- function() {

  rm(list = ls(envir = .package_source_md5_cache),
     envir = .package_source_md5_cache)

  return(invisible(TRUE))
}

.fit_cache_required_source_files <- function() {

  # Keep this scoped to code that can change saved fit objects or cached fit
  # extensions: a file is listed when it executes while a model is constructed
  # (constructors, data and prior preparation, JAGS model text, native code of
  # the JAGS module and of the R-native routines they call), while the sampler
  # runs, or while the LOO, WAIC, and marginal-likelihood results stored with
  # a cached fit are computed. Post-fit methods (summaries, plots, predictions,
  # hypotheses and densities, diagnostics) are tested against the cached
  # objects and stay unlisted.
  r_files   <- c(
    "R/BMA.glmm.R",
    "R/BMA.mv.R",
    "R/BMA.norm.R",
    "R/RoBMA.R",
    "R/RoBMA.mv.R",
    "R/bPET-PEESE.mv.R",
    "R/bPEESE.R",
    "R/bPET.R",
    "R/block-covariance.R",
    "R/brma.glmm.R",
    "R/brma.mv.R",
    "R/brma.norm.R",
    "R/brma-mv-targets.R",
    "R/bselmodel.R",
    "R/bselmodel.mv.R",
    "R/brma-mv-known-r.R",
    "R/covariance-factorization.R",
    "R/distributions.R",
    "R/evaluate.R",
    "R/fit.R",
    "R/formula-design.R",
    "R/glmm-aghq.R",
    "R/input-data.R",
    "R/input-data-mv.R",
    "R/input-object.R",
    "R/input-priors.R",
    "R/input-priors-assignment.R",
    "R/input-priors-check-list.R",
    "R/input-priors-documentation.R",
    "R/input-priors-formula.R",
    "R/input-priors-heterogeneity-allocation.R",
    "R/jags-formula-args.R",
    "R/known-v-gls.R",
    "R/known-v-preflight.R",
    "R/known-v-representation.R",
    "R/loo.R",
    "R/marglik.R",
    "R/outcome-helpers.R",
    "R/log-lik-cluster.R",
    "R/log-lik-cluster-glmm.R",
    "R/log-lik-cluster-normal.R",
    "R/log-lik.R",
    "R/log-lik-known-v.R",
    "R/pdf-outcome-glmm.R",
    "R/pdf-outcome-normal.R",
    "R/pdf-utils.R",
    "R/pdf.R",
    "R/predict-brma-mv.R",
    "R/priors.R",
    "R/random-effects-compile.R",
    "R/random-inclusion.R",
    "R/random-marginal-vcov.R",
    "R/rng.R",
    "R/selection-cache.R",
    "R/selection-likelihood.R",
    "R/selection-mapping.R",
    "R/selection-normalizer-grid.R",
    "R/unit_level.R",
    "R/utilities.R",
    "R/zzz.R"
  )
  src_files <- c(
    "src/Makevars.in",
    "src/Makevars.ucrt",
    "src/Makevars.win",
    "src/RoBMA.cc",
    "src/init.c",
    "src/distributions/DSELNORMCLUSTERSTEP.cc",
    "src/distributions/DSELNORMCLUSTERSTEP.h",
    "src/distributions/DSELNORMFACTORSTEP.cc",
    "src/distributions/DSELNORMFACTORSTEP.h",
    "src/distributions/DSELNORMKERNEL.cc",
    "src/distributions/DSELNORMKERNEL.h",
    "src/distributions/DSELNORMMVSTEP.cc",
    "src/distributions/DSELNORMMVSTEP.h",
    "src/distributions/DSELNORMSAMPLINGCONDITIONED.cc",
    "src/distributions/DSELNORMSAMPLINGCONDITIONED.h",
    "src/distributions/DSELNORMSTEP.cc",
    "src/distributions/DSELNORMSTEP.h",
    "src/distributions/DSELNORMSTEPSWITCH.cc",
    "src/distributions/DSELNORMSTEPSWITCH.h",
    "src/distributions/selnorm-jags-bounds.h",
    "src/distributions/DKNOWNVMNORM.cc",
    "src/distributions/DKNOWNVMNORM.h",
    "src/distributions/DWB.cc",
    "src/distributions/DWB.h",
    "src/distributions/DWN.cc",
    "src/distributions/DWN.h",
    "src/distributions/DWP.cc",
    "src/distributions/DWP.h",
    "src/glmm-aghq.cc",
    "src/glmm-aghq.h",
    "src/glmm-binomial-loglik.h",
    "src/r-native-api.h",
    "src/r-native-boundary.h",
    "src/r-glmm.cc",
    "src/r-glmm-common.cc.inc",
    "src/r-glmm-marginal.cc.inc",
    "src/r-known-v.cc",
    "src/r-known-v-common.cc.inc",
    "src/r-known-v-evaluate.cc.inc",
    "src/r-known-v-low-rank.cc.inc",
    "src/r-known-v-plan.cc.inc",
    "src/r-known-v-structured.cc.inc",
    "src/r-selnorm.cc",
    "src/r-selnorm-cluster.cc.inc",
    "src/r-selnorm-common.cc.inc",
    "src/r-selnorm-funnel-common.cc.inc",
    "src/r-selnorm-funnel-zcurve.cc.inc",
    "src/r-selnorm-kernel.cc.inc",
    "src/r-selnorm-loglik.cc.inc",
    "src/r-selnorm-mv.cc.inc",
    "src/r-selnorm-sampling-conditioned.cc.inc",
    "src/samplers/CoarseCorrectedSlice.cc",
    "src/samplers/CoarseCorrectedSlice.h",
    "src/selnorm/selnorm.cc",
    "src/selnorm/selnorm-api.cc.inc",
    "src/selnorm/selnorm-boundary.cc.inc",
    "src/selnorm/selnorm-cdf-table-data.h",
    "src/selnorm/selnorm-event.cc.inc",
    "src/selnorm/selnorm-fma.h",
    "src/selnorm/selnorm-mixture-compression.h",
    "src/selnorm/selnorm-mv.cc",
    "src/selnorm/selnorm-mv.h",
    "src/selnorm/selnorm-parallel.h",
    "src/selnorm/selnorm-phack.cc.inc",
    "src/selnorm/selnorm-probability.cc.inc",
    "src/selnorm/selnorm-sampling-conditioned.cc.inc",
    "src/selnorm/selnorm-step.cc.inc",
    "src/selnorm/selnorm-tail-kernel.cc.inc",
    "src/selnorm/selnorm-tail.h",
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

  text_extensions <- c(
    "c", "cc", "cpp", "cxx", "f", "f77", "f90", "f95",
    "h", "hh", "hpp", "hxx", "inc", "in"
  )
  text_basenames <- c(
    "NAMESPACE", "Makevars", "Makevars.in", "Makevars.ucrt", "Makevars.win",
    "configure", "configure.win", "cleanup", "cleanup.win"
  )
  if (extension %in% text_extensions || basename(path) %in% text_basenames) {
    lines      <- readLines(path, warn = FALSE)
    normalized <- tempfile("robma-package-source-", fileext = ".txt")
    writeLines(lines, normalized, useBytes = TRUE)
    on.exit(unlink(normalized), add = TRUE)

    return(unname(tools::md5sum(normalized)))
  }

  return(unname(tools::md5sum(path)))
}

.fit_cache_source_fingerprint <- function(package_root, source_files) {

  # Only size and modification time enter the fingerprint; skipping the extra
  # columns avoids reading every file's header on Windows.
  info <- file.info(source_files, extra_cols = FALSE)
  paste(
    package_root,
    paste(source_files, info[["size"]], info[["mtime"]], sep = ":", collapse = "\n"),
    sep = "\n"
  )
}

package_source_md5 <- function() {

  package_root <- .fit_cache_source_root()
  if (is.na(package_root)) {
    assign("value", NA_character_, envir = .package_source_md5_cache)
    assign("fingerprint", NA_character_, envir = .package_source_md5_cache)
    return(NA_character_)
  }

  source_files <- .fit_cache_source_files(package_root = package_root)
  if (length(source_files) == 0L) {
    assign("value", NA_character_, envir = .package_source_md5_cache)
    assign("fingerprint", NA_character_, envir = .package_source_md5_cache)
    return(NA_character_)
  }

  fingerprint <- .fit_cache_source_fingerprint(
    package_root = package_root,
    source_files = source_files
  )
  if (exists("value", envir = .package_source_md5_cache, inherits = FALSE) &&
      exists("fingerprint", envir = .package_source_md5_cache, inherits = FALSE) &&
      identical(
        get("fingerprint", envir = .package_source_md5_cache, inherits = FALSE),
        fingerprint
      )) {
    return(get("value", envir = .package_source_md5_cache, inherits = FALSE))
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
  assign("fingerprint", fingerprint, envir = .package_source_md5_cache)

  return(value)
}

bayestools_backend_fingerprint <- function(
    required = TRUE,
    value    = BayesTools::fit_backend_fingerprint()) {

  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[[:xdigit:]]{32}$", value)
  if (!valid && required) {
    stop(
      "BayesTools backend fingerprint is unavailable or invalid.",
      call. = FALSE
    )
  }
  if (!valid) {
    return(NA_character_)
  }

  return(value)
}


.fit_cache_versions_cache <- new.env(parent = emptyenv())


fit_cache_versions <- function() {

  if (exists("versions", envir = .fit_cache_versions_cache, inherits = FALSE)) {
    return(get("versions", envir = .fit_cache_versions_cache, inherits = FALSE))
  }
  source_root <- .fit_cache_source_root()
  description <- if (is.na(source_root)) {
    ""
  } else {
    file.path(source_root, "DESCRIPTION")
  }
  robma_version <- if (file.exists(description)) {
    unname(read.dcf(description, fields = "Version")[1L, 1L])
  } else {
    as.character(utils::packageVersion("RoBMA"))
  }
  jags_version <- tryCatch(
    as.character(rjags::jags.version()),
    error = function(error) NA_character_
  )

  versions <- list(
    robma_version      = robma_version,
    bayestools_version = as.character(utils::packageVersion("BayesTools")),
    r_version          = paste(R.version$major, R.version$minor, sep = "."),
    jags_version       = jags_version
  )
  assign("versions", versions, envir = .fit_cache_versions_cache)

  return(versions)
}


# Bind discovery to the helper's actual checkout, including when a different
# source root is tested in the same R process. Validation memoization stays live.
.fit_cache_helper_root <- dirname(dirname(.common_functions_dir()))
if (.fit_cache_source_tree_available(.fit_cache_helper_root) &&
    !identical(.fit_cache_source_root(), .fit_cache_helper_root)) {
  .clear_package_source_md5_cache()
  assign("source_root", .fit_cache_helper_root, .package_source_md5_cache)
}

.fit_cache_generations <- getOption("RoBMA.test.producer_generations")
if (!is.environment(.fit_cache_generations)) {
  .fit_cache_generations <- new.env(parent = emptyenv())
  options(RoBMA.test.producer_generations = .fit_cache_generations)
}

.fit_cache_loaded_product <- function() {

  packages <- c("RoBMA", "BayesTools")
  return(lapply(setNames(packages, packages), function(package) {
    namespace <- asNamespace(package)
    dll <- getLoadedDLLs()[[package]][["path"]]
    list(namespace = namespace, dll = dll, dll_md5 = unname(tools::md5sum(dll)))
  }))
}

# Call with force = TRUE only at an explicit runner/product-load lifecycle.
# Re-sourcing a fit file's helpers must retain the original producer stamp.
.fit_cache_generation_begin <- function(force = FALSE) {

  root <- .fit_cache_source_root()
  key <- if (is.na(root)) "<source unavailable>" else root
  product <- .fit_cache_loaded_product()
  previous <- if (exists(key, .fit_cache_generations, inherits = FALSE)) {
    get(key, .fit_cache_generations, inherits = FALSE)
  } else NULL
  same_namespaces <- !is.null(previous) && all(vapply(names(product), function(package) {
    identical(product[[package]][["namespace"]], previous$product[[package]][["namespace"]])
  }, TRUE))
  if (!force && same_namespaces) return(previous)

  rm(list = ls(.fit_cache_versions_cache), envir = .fit_cache_versions_cache)
  generation <- new.env(parent = emptyenv())
  generation$product <- product
  generation$versions <- fit_cache_versions()
  generation$package_md5 <- package_source_md5()
  generation$backend <- bayestools_backend_fingerprint(required = FALSE)
  files <- unique(fit_catalog()[["source_file"]])
  generation$builders <- setNames(vapply(files, source_file_md5, character(1)), files)
  generation$consumed <- list()
  generation$conflicts <- character()
  assign(key, generation, .fit_cache_generations)
  return(generation)
}

.fit_cache_assert_producer <- function(name) {

  generation <- .fit_cache_generation_begin()
  entry <- fit_catalog_entry(name)
  source_file <- if (is.null(entry)) NA_character_ else entry[["source_file"]]
  source_md5 <- generation$builders[source_file]
  valid <- length(source_md5) == 1L && !is.na(source_md5) &&
    !is.na(generation$package_md5) && !is.na(generation$backend) &&
    identical(unname(source_md5), source_file_md5(source_file)) &&
    identical(generation$package_md5, package_source_md5()) &&
    identical(generation$backend, bayestools_backend_fingerprint(required = FALSE)) &&
    identical(generation$product, .fit_cache_loaded_product())
  if (!valid) {
    stop("Cached fit producer source or build changed. Reload the product and rerun the fit generation.",
         call. = FALSE)
  }
  return(generation)
}

# Devtools has loaded the product before test helpers run. A profile runner may
# source helpers earlier: let it load its own product before beginning a run.
if ("RoBMA" %in% loadedNamespaces()) .fit_cache_generation_begin()

fit_cache_metadata <- function(name, fit, info = NULL) {

  entry       <- fit_catalog_entry(name)
  source_file <- if (is.null(entry)) NA_character_ else entry[["source_file"]]
  generation  <- .fit_cache_generation_begin()
  versions    <- generation$versions

  metadata <- list(
    schema_version                 = FIT_CACHE_SCHEMA_VERSION,
    name                           = name,
    saved_at                       = format(Sys.time(), usetz = TRUE),
    robma_version                  = versions[["robma_version"]],
    bayestools_version             = versions[["bayestools_version"]],
    r_version                      = versions[["r_version"]],
    jags_version                   = versions[["jags_version"]],
    fit_class                      = class(fit),
    source_file                    = source_file,
    source_file_md5                = unname(generation$builders[source_file]),
    package_source_md5             = generation$package_md5,
    bayestools_backend_fingerprint = generation$backend,
    has_loo                        = !is.null(fit[["loo"]]),
    has_waic                       = !is.null(fit[["waic"]]),
    has_marglik                    = !is.null(fit[["marglik"]]),
    has_metafor_info               = !is.null(info) &&
      "metafor" %in% names(info) && !is.null(info[["metafor"]])
  )
  if (!is.null(entry) && nrow(entry[["dependencies"]][[1L]]) > 0L) {
    metadata[["dependencies"]] <- .fit_cache_dependencies_for_save(entry)
  }

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

.fit_cache_artifact_receipt <- function(name, artifact) {

  path <- fit_cache_paths(name)[[artifact]]
  if (!file.exists(path)) {
    stop("Required dependency artifact is missing: ", artifact, " ", name, call. = FALSE)
  }
  metadata <- load_fit_metadata(name)
  list(name = name, artifact = artifact, content_md5 = unname(tools::md5sum(path)),
       provenance = metadata[setdiff(names(metadata), "saved_at")])
}

.fit_cache_receipt_key <- function(name, artifact) {
  normalizePath(fit_cache_paths(name)[[artifact]], winslash = "/", mustWork = FALSE)
}

.fit_cache_record_consumption <- function(receipt) {

  generation <- .fit_cache_generation_begin()
  key <- .fit_cache_receipt_key(receipt[["name"]], receipt[["artifact"]])
  first <- generation$consumed[[key]]
  if (is.null(first)) {
    generation$consumed[[key]] <- receipt
  } else if (!identical(first, receipt)) {
    generation$conflicts <- union(generation$conflicts, key)
  }
  return(invisible(TRUE))
}

.fit_cache_dependencies_for_save <- function(entry) {

  dependencies <- entry[["dependencies"]][[1L]]
  generation <- .fit_cache_generation_begin()
  receipts <- vector("list", nrow(dependencies))
  for (i in seq_len(nrow(dependencies))) {
    name <- dependencies[["name"]][[i]]
    artifact <- dependencies[["artifact"]][[i]]
    key <- .fit_cache_receipt_key(name, artifact)
    if (key %in% generation$conflicts) {
      stop("Cached fit has conflicting dependency consumption: ", artifact, " ", name,
           ". Reload the product and rerun generation.", call. = FALSE)
    }
    receipt <- generation$consumed[[key]]
    if (is.null(receipt)) {
      stop("Cached fit requires dependency consumption in this generation: ", artifact, " ", name,
           ". Load the declared input before building the fit.", call. = FALSE)
    }
    problems <- validate_cached_fit(name)
    if (length(problems) > 0L) {
      stop("Cached fit dependency changed or is invalid: ", name, ": ",
           paste(problems, collapse = "; "), call. = FALSE)
    }
    current <- .fit_cache_artifact_receipt(name, artifact)
    if (!identical(receipt, current)) {
      stop("Cached fit dependency changed after consumption: ", artifact, " ", name,
           ". Rerun generation with the current inputs.", call. = FALSE)
    }
    receipts[[i]] <- receipt
  }
  return(receipts)
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
                                check_files = TRUE, deep = FALSE,
                                bayestools_fingerprint =
                                  bayestools_backend_fingerprint(
                                    required = FALSE
                                  ),
                                package_md5 = NULL, catalog = NULL,
                                .active = character()) {

  messages <- character()
  paths    <- fit_cache_paths(name)
  if (is.null(catalog)) {
    catalog <- fit_catalog()
  }
  entry    <- fit_catalog_entry(name, catalog)
  if (name %in% .active) {
    return(paste0("cache dependency cycle: ", paste(c(.active, name), collapse = " -> ")))
  }

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
  } else if (!is.list(metadata)) {
    messages <- c(messages, "metadata must be a list")
    metadata <- NULL
  }

  if (!is.null(entry) && !is.null(metadata)) {
    dependencies <- entry[["dependencies"]][[1L]]
    fields <- c(FIT_CACHE_METADATA_FIELDS, if (nrow(dependencies) > 0L) "dependencies")
    if (!identical(names(metadata), fields)) {
      messages <- c(messages, "cache metadata fields changed")
    }
    if (nrow(dependencies) > 0L) {
      saved_dependencies <- metadata[["dependencies"]]
      if (!is.list(saved_dependencies) || length(saved_dependencies) != nrow(dependencies)) {
        messages <- c(messages, "cache dependency fingerprints missing or changed")
      }
      for (i in seq_len(nrow(dependencies))) {
        parent <- dependencies[["name"]][[i]]
        artifact <- dependencies[["artifact"]][[i]]
        problems <- validate_cached_fit(parent, deep = deep, catalog = catalog,
                                         package_md5 = package_md5,
                                         bayestools_fingerprint = bayestools_fingerprint,
                                         .active = c(.active, name))
        if (length(problems) > 0L) {
          messages <- c(messages, paste0("cache dependency ", parent, ": ", problems))
        }
        current <- tryCatch(.fit_cache_artifact_receipt(parent, artifact), error = identity)
        if (inherits(current, "error")) {
          messages <- c(messages, conditionMessage(current))
        } else if (length(saved_dependencies) < i || !identical(saved_dependencies[[i]], current)) {
          messages <- c(messages, paste0("cache dependency artifact changed: ", artifact, " ", parent))
        }
      }
    }
    if (is.null(metadata[["schema_version"]]) ||
        !identical(
          metadata[["schema_version"]],
          FIT_CACHE_SCHEMA_VERSION
        )) {
      messages <- c(messages, "cache schema changed")
    }
    if (is.null(metadata[["name"]]) || !identical(metadata[["name"]], name)) {
      messages <- c(messages, "metadata name mismatch")
    }
    versions <- fit_cache_versions()
    version_messages <- c(
      robma_version      = "RoBMA version changed",
      bayestools_version = "BayesTools version changed",
      r_version          = "R version changed",
      jags_version       = "JAGS version changed"
    )
    for (field in names(version_messages)) {
      if (is.null(metadata[[field]]) ||
          !identical(metadata[[field]], versions[[field]])) {
        messages <- c(messages, unname(version_messages[[field]]))
      }
    }
    if (is.null(metadata[["source_file"]]) ||
        !identical(metadata[["source_file"]], entry[["source_file"]])) {
      messages <- c(messages, "metadata source file mismatch")
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

    if (check_source) {
      expected_md5 <- source_file_md5(entry[["source_file"]])
      if (is.na(expected_md5)) {
        messages <- c(messages, "source file fingerprint unavailable")
      } else if (is.null(metadata[["source_file_md5"]]) ||
                 !identical(metadata[["source_file_md5"]], expected_md5)) {
        messages <- c(messages, "source file hash changed")
      }

      if (is.null(package_md5)) {
        package_md5 <- package_source_md5()
      }
      if (is.na(package_md5)) {
        messages <- c(messages, "cache source fingerprint unavailable")
      } else if (is.null(metadata[["package_source_md5"]]) ||
                 !identical(metadata[["package_source_md5"]], package_md5)) {
        messages <- c(messages, "cache source hash changed")
      }

      bayestools_fingerprint <- bayestools_backend_fingerprint(
        required = FALSE,
        value    = bayestools_fingerprint
      )
      if (is.na(bayestools_fingerprint)) {
        messages <- c(messages, "BayesTools backend fingerprint unavailable")
      } else if (is.null(metadata[["bayestools_backend_fingerprint"]]) ||
                 !identical(
                   metadata[["bayestools_backend_fingerprint"]],
                   bayestools_fingerprint
                 )) {
        messages <- c(messages, "BayesTools backend changed")
      }
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
    if (exists(".announce_existing_visual_snapshots", mode = "function")) {
      .announce_existing_visual_snapshots()
    }
    skip("No valid cached fits available. Run `devtools::test(filter = '01-')` first.")
  }
}

skip_if_missing_fits <- function(names, active_only = TRUE) {

  # Validate only the requested fits; every other cached fit is irrelevant to
  # this skip decision and validating all of them dominated cached test time.
  missing <- setdiff(
    names,
    list_fits(name = names, validate = TRUE, active_only = active_only)
  )
  if (length(missing) > 0) {
    if (exists(".announce_existing_visual_snapshots", mode = "function")) {
      .announce_existing_visual_snapshots()
    }
    skip(paste0("Required pre-fitted models missing or stale: ", paste(missing, collapse = ", ")))
  }
}

skip_if_not_full_diagnostics <- function(reason) {

  if (!is_certification_profile()) {
    skip(paste(
      "Skipping extended diagnostic redundancy check by default.",
      reason,
      "Run the certification profile to include it."
    ))
  }
}

skip_if_not_certification <- function(reason = NULL) {

  if (!is_certification_profile()) {
    detail <- if (is.null(reason)) "" else paste0(" ", reason)
    skip(paste0(
      "Skipping numerical certification case in the standard profile.",
      detail
    ))
  }

  return(invisible(FALSE))
}

# Skip model fitting if a valid cached fit exists. Refit by setting
# ROBMA_TEST_FORCE_REFIT=TRUE or ROBMA_TEST_SKIP_REFIT=FALSE.
skip_refit_if_cached <- function(name) {

  if (is_true_env("ROBMA_TEST_FORCE_REFIT")) {
    return(invisible(FALSE))
  }

  skip_refit_env <- Sys.getenv("ROBMA_TEST_SKIP_REFIT")
  skip_refit     <- if (skip_refit_env == "") TRUE else !is_false_env("ROBMA_TEST_SKIP_REFIT")

  catalog     <- active_fit_catalog()
  source_file <- paste0("test-01-", name, ".R")
  fit_names   <- if (name %in% catalog[["name"]]) {
    name
  } else if (source_file %in% catalog[["source_file"]]) {
    catalog[catalog[["source_file"]] == source_file, "name"]
  } else {
    catalog_group_fits(name)
  }

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
  .fit_cache_assert_producer(name)
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

  assign(name, list(object = fit, receipt = .fit_cache_artifact_receipt(name, "fit")),
         envir = .fit_object_cache)
  if (!is.null(info)) {
    assign(name, list(object = info, receipt = .fit_cache_artifact_receipt(name, "info")),
           envir = .info_object_cache)
  } else {
    .clear_info_object_cache(name)
  }

  return(invisible(TRUE))
}

.load_cache_artifact <- function(name, artifact, validate) {

  path <- fit_cache_paths(name)[[artifact]]
  before <- if (file.exists(path)) .fit_cache_artifact_receipt(name, artifact) else NULL

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

  cache <- if (artifact == "fit") .fit_object_cache else .info_object_cache
  cached <- if (exists(name, cache, inherits = FALSE)) get(name, cache, inherits = FALSE) else NULL
  if (!is.null(cached) && identical(cached[["receipt"]], before)) {
    object <- cached[["object"]]
  } else {
    object <- if (artifact == "fit") read_cached_fit(name) else read_cached_info(name)
  }
  if (artifact == "fit" && is.null(object)) {
    stop(paste0("Cached fit '", name, "' could not be read."), call. = FALSE)
  }
  after <- if (file.exists(path)) .fit_cache_artifact_receipt(name, artifact) else NULL
  if (!identical(before, after)) {
    stop("Cached artifact changed during read: ", artifact, " ", name, call. = FALSE)
  }
  assign(name, list(object = object, receipt = after), cache)
  if (validate && !is.null(after)) .fit_cache_record_consumption(after)
  return(object)
}

load_fit <- function(name, validate = TRUE) {
  return(.load_cache_artifact(name, "fit", validate))
}

load_info <- function(name, validate = TRUE) {
  return(.load_cache_artifact(name, "info", validate))
}

list_fits <- function(name, feature, class, family, has_metafor, has_loo,
                      has_waic, has_marglik, tier, validate = TRUE,
                      deep = FALSE, active_only = TRUE) {

  files <- suppressWarnings(list.files(temp_fits_dir, pattern = "\\.RDS$"))
  files <- sub("\\.RDS$", "", files)

  # One catalog per call: validating a fit looks its entry up in the catalog.
  catalog <- fit_catalog()

  if (active_only) {
    files <- intersect(files, active_fit_catalog(catalog)[["name"]])
  }

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
      tier        = if (!missing(tier)) tier else NULL,
      catalog     = catalog
    )
    files <- intersect(files, selected)
  }

  if (validate) {
    package_md5 <- package_source_md5()
    files <- files[vapply(
      files,
      function(name) {
        length(validate_cached_fit(
          name        = name,
          deep        = deep,
          package_md5 = package_md5,
          catalog     = catalog
        )) == 0L
      },
      logical(1)
    )]
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
  names <- intersect(
    names,
    list_fits(validate = FALSE, active_only = TRUE)
  )

  return(lazy_cached_objects(
    names,
    function(name) load_fit(name, validate = FALSE)
  ))
}

lazy_infos <- function(names = list_fits(), validate = TRUE) {

  if (validate) {
    skip_if_missing_fits(names)
  }
  names <- intersect(
    names,
    list_fits(validate = FALSE, active_only = TRUE)
  )

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
  if (!is.character(i) || length(i) != 1L) {
    return(NULL)
  }
  if (!i %in% object_names) {
    skip(paste0(
      "Required pre-fitted model missing or stale: ", i,
      ". Run `devtools::test(filter = '01-', reporter = 'llm')` first."
    ))
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
    fit_names <- unique(unlist(lapply(name, function(value) {
      selected <- catalog_group_fits(value, active_only = FALSE)
      if (length(selected) == 0L) value else selected
    }), use.names = FALSE))
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
  .clear_package_source_md5_cache()

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
