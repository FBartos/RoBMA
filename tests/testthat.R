library(testthat)
library(RoBMA)

on_cran <- get("on_cran", envir = asNamespace("testthat"), inherits = FALSE)

if (on_cran()) {
  if (!nzchar(Sys.getenv("ROBMA_TEST_FILES_DIR"))) {
    Sys.setenv(ROBMA_TEST_FILES_DIR = file.path(tempdir(), "robma-test-files"))
  }
  test_check("RoBMA", filter = "cran-smoke")
} else {
  test_check("RoBMA")
}
