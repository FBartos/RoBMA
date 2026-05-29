Sys.setenv(AGENT = "1")
source(file.path("tests", "testthat", "common-functions.R"))
clean_cached_fits()
Sys.setenv(ROBMA_TEST_EXTENDED = "TRUE")
Sys.setenv(ROBMA_TEST_FULL_VISUALS = "TRUE")
Sys.setenv(ROBMA_TEST_FULL_DIAGNOSTICS = "TRUE")

devtools::test(reporter = "llm", stop_on_failure = TRUE)
