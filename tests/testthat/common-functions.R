# ============================================================================ #
# CONFIGURATION: Set to TRUE to regenerate reference files, FALSE to run tests
# ============================================================================ #
if (!exists("GENERATE_REFERENCE_FILES")) {
  GENERATE_REFERENCE_FILES <- FALSE
}

# Get the directory where prefitted models are stored
temp_fits_dir <- Sys.getenv("ROBMA_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  temp_fits_dir <- file.path(tempdir(), "RoBMA_test_fits")
}

# Setup directory for saving fitted models
temp_fits_dir <- file.path(tempdir(), "RoBMA_test_fits")
if (!dir.exists(temp_fits_dir)) {
  dir.create(temp_fits_dir, showWarnings = FALSE, recursive = TRUE)
}

# Set environment variable so other test files can locate pre-fitted models
Sys.setenv(ROBMA_TEST_FITS_DIR = temp_fits_dir)

# NOTE: File-level skip_on_cran() was removed intentionally.
# Each test file should manage its own skip conditions appropriately.
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
  model_registry_file <- file.path(temp_fits_dir, "model_registry.RDS")
  if (!file.exists(model_registry_file)) {
    skip("Pre-fitted models not found. Run `test(filter = 'test-01')` first.")
  }
}

# ============================================================================ #
# HELPER FUNCTIONS: Model Fit Saving / Loading
# ============================================================================ #

save_fit     <- function(fit, name, marglik = NULL) {

  # Save model fit
  saveRDS(fit, file = file.path(temp_fits_dir, paste0(name, ".RDS")))

  # Save marglik if provided
  if (!is.null(marglik)) {
    saveRDS(marglik, file = file.path(temp_fits_dir, paste0(name, "_marglik.RDS")))
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
  marglik <- readRDS(file = file.path(temp_fits_dir, paste0(name, "_marglik.RDS")))
  return(marglik)
}
list_fits    <- function(name) {

  files <- list.files(temp_fits_dir)
  files <- files[!grepl("_marglik.RDS", files)]
  files <- gsub(".RDS", "", files)

  return(files)
}
