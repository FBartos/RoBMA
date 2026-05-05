library(testthat)
library(RoBMA)

on_cran <- get("on_cran", envir = asNamespace("testthat"), inherits = FALSE)

if (on_cran()) {
  test_check("RoBMA", filter = "cran-smoke")
} else {
  test_check("RoBMA")
}
