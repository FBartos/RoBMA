clean_cached_fits()
Sys.unsetenv("AGENT")
Sys.setenv(ROBMA_TEST_FULL_VISUALS = TRUE)
Sys.setenv(ROBMA_TEST_FULL_DIAGNOSTICS = TRUE)

devtools::test()
