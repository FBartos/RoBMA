context("Funnel plots")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()


test_that("Funnel plots for brma models", {

  # basic funnel plots
  for (name in names(fits)) {
    # only test normal outcome models
    if (RoBMA:::.outcome_type(fits[[name]]) == "norm") {
      vdiffr::expect_doppelganger(paste0("base_funnel_", name), function() funnel(fits[[name]]))
      vdiffr::expect_doppelganger(paste0("ggplot_funnel_", name), funnel(fits[[name]], plot_type = "ggplot"))

      # test without publication bias
      vdiffr::expect_doppelganger(paste0("base_funnel_no_pb_", name), function() funnel(fits[[name]], bias_adjusted = FALSE))

      # test without heterogeneity
      vdiffr::expect_doppelganger(paste0("base_funnel_no_het_", name), function() funnel(fits[[name]], heterogeneity_adjusted = FALSE))
    }
  }

  # customization
  name <- "bcg_meta-analysis"
  vdiffr::expect_doppelganger("base_funnel_custom", function() funnel(fits[[name]], shape = 24, size = 1.5))
  vdiffr::expect_doppelganger("ggplot_funnel_custom", funnel(fits[[name]], plot_type = "ggplot", shape = 24, size = 3))

})
