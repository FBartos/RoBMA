# ============================================================================ #
# CONFIGURATION: Set to TRUE to regenerate reference files, FALSE to run tests
# ============================================================================ #
if (!exists("GENERATE_REFERENCE_FILES")) {
  GENERATE_REFERENCE_FILES <- FALSE
}

# Get the directory where prefitted models are stored
test_files_dir <- Sys.getenv("ROBMA_TEST_FILES_DIR")
if (test_files_dir == "" || !dir.exists(test_files_dir)) {
  test_files_dir <- file.path(tempdir(), "RoBMA_test_files")
}

# Setup directory for saving fitted models
temp_fits_dir    <- file.path(test_files_dir, "fits")
temp_marglik_dir <- file.path(test_files_dir, "margliks")
temp_info_dir    <- file.path(test_files_dir, "info")
temp_temp_dir    <- file.path(test_files_dir, "temp")
if (!dir.exists(temp_fits_dir))    dir.create(temp_fits_dir,    showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_marglik_dir)) dir.create(temp_marglik_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_info_dir))    dir.create(temp_info_dir,    showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(temp_temp_dir))    dir.create(temp_temp_dir,    showWarnings = FALSE, recursive = TRUE)

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
    writeLines(capture_output_lines(table, print = TRUE, width = 150),
               file.path(print_dir, filename))
  } else {
    # Test mode
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      actual_output   <- capture_output_lines(table, print = TRUE, width = 150)
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

# Skip if pre-fitted models are not available
skip_if_no_fits <- function() {
  if (length(list.files(temp_fits_dir)) == 0) {
    skip("Pre-fitted models not found. Run `test(filter = 'test-01')` first.")
  }
}

# Skip model fitting if cached fits exist and ROBMA_TEST_SKIP_REFIT is TRUE
skip_refit_if_cached <- function(name) {
  # refitting settings
  skip_refit <- Sys.getenv("ROBMA_TEST_SKIP_REFIT")
  skip_refit <- skip_refit != "" && as.logical(skip_refit)

  # fitted indicator
  fitted_indicator <- file.exists(file.path(temp_temp_dir, paste0(name, ".txt")))

  if (skip_refit && fitted_indicator) {
    skip("Skipping model refitting: cached fits exist and ROBMA_TEST_SKIP_REFIT=TRUE.")
  }

  # tests are not going to be skipped -- add fits done indicator into `temp_temp_dir`
  file.create(file.path(temp_temp_dir, paste0(name, ".txt")))
}

# ============================================================================ #
# HELPER FUNCTIONS: Model Fit Saving / Loading
# ============================================================================ #

save_fit     <- function(name, fit, marglik = NULL, info = NULL) {

  # Save model fit
  saveRDS(fit, file = file.path(temp_fits_dir, paste0(name, ".RDS")))

  # Save marglik if provided
  if (!is.null(marglik)) {
    saveRDS(marglik, file = file.path(temp_marglik_dir, paste0(name, ".RDS")))
  }

  # Save info if provided
  if (!is.null(info)) {
    saveRDS(info, file = file.path(temp_info_dir, paste0(name, ".RDS")))
  }

  return(invisible(TRUE))
}
load_fit     <- function(name) {

  # load model fit
  fit <- readRDS(file = file.path(temp_fits_dir, paste0(name, ".RDS")))
  return(fit)
}
load_marglik <- function(name) {

  # load model marglik
  marglik <- try(readRDS(file = file.path(temp_marglik_dir, paste0(name, ".RDS"))), silent = TRUE)
  if (inherits(marglik, "try-error")) {
    return(list())
  } else {
    return(marglik)
  }
}
load_info    <- function(name) {

  # load model info
  info <- suppressWarnings(try(readRDS(file = file.path(temp_info_dir, paste0(name, ".RDS"))), silent = TRUE))
  if (inherits(info, "try-error")) {
    return(list())
  } else {
    return(info)
  }
}

list_fits    <- function(name) {

  files <- suppressWarnings(list.files(temp_fits_dir))
  files <- gsub(".RDS", "", files)

  return(files)
}

clean_cached_fits <- function(name) {

  if (!missing(name)) {
    # remove only the specific `name`` fitted indicator files side-effects from `temp_temp_dir`
    file.remove(file.path(temp_temp_dir, paste0(name, ".txt")))
  } else {
    # Remove all cached files from test directories
    unlink(temp_fits_dir,    recursive = TRUE)
    unlink(temp_marglik_dir, recursive = TRUE)
    unlink(temp_info_dir,    recursive = TRUE)
    unlink(temp_temp_dir,    recursive = TRUE)

    # Recreate empty directories
    dir.create(temp_fits_dir,    showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_marglik_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_info_dir,    showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_temp_dir,    showWarnings = FALSE, recursive = TRUE)
  }

  return(invisible(TRUE))
}
