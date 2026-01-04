
# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()


test_that("Prior and posterior distributions for model parameters", {

  # effect
  for (name in names(fits)) {
    vdiffr::expect_doppelganger(paste0("baseplot_pp_mu_", name), plot(fits[[name]], "mu", prior = TRUE))
    vdiffr::expect_doppelganger(paste0("ggplot_pp_mu_", name), plot(fits[[name]], "mu", prior = TRUE, plot_type = "ggplot"))
  }

  # heterogeneity
  for (name in names(fits)) {
    vdiffr::expect_doppelganger(paste0("base_pp_tau_", name), plot(fits[[name]], "tau", prior = TRUE))
    vdiffr::expect_doppelganger(paste0("ggplot_pp_tau_", name), plot(fits[[name]], "tau", prior = TRUE, plot_type = "ggplot"))
  }

  # variance allocation
  for (name in names(fits)) {
    if (RoBMA:::.is_multilevel(fits)) {
      vdiffr::expect_doppelganger(paste0("base_pp_rho_", name), plot(fits[[name]], "rho", prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_rho_", name), plot(fits[[name]], "rho", prior = TRUE, plot_type = "ggplot"))
    }
  }

  # effect regression
  for (name in names(fits)) {
    if (RoBMA:::.is_mods(fits)) {
      for (pars in info[[name]][["mods"]]) {
        vdiffr::expect_doppelganger(paste0("base_pp_mods_", pars, "_", name), plot(fits[[name]], parameter_mods = pars, prior = TRUE))
        vdiffr::expect_doppelganger(paste0("ggplot_pp_scale_", pars, "_",, name), plot(fits[[name]], parameter_mods = pars, prior = TRUE, plot_type = "ggplot"))
      }
    }
  }

  # scale regression
  for (name in names(fits)) {
    if (RoBMA:::.is_scale(fits)) {
      for (pars in info[[name]][["scale"]]) {
        vdiffr::expect_doppelganger(paste0("base_pp_scale_", pars, "_", name), plot(fits[[name]], parameter_mods = pars, prior = TRUE))
        vdiffr::expect_doppelganger(paste0("ggplot_pp_scale_", pars, "_",, name), plot(fits[[name]], parameter_mods = pars, prior = TRUE, plot_type = "ggplot"))
      }
    }
  }

  # test plot customization


})
