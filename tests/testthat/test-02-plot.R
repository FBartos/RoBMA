
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
    vdiffr::expect_doppelganger(paste0("baseplot_pp_mu_", name), function() plot(fits[[name]], "mu", prior = TRUE))
    vdiffr::expect_doppelganger(paste0("ggplot_pp_mu_", name), plot(fits[[name]], "mu", prior = TRUE, plot_type = "ggplot"))
  }

  # heterogeneity
  for (name in names(fits)) {
    vdiffr::expect_doppelganger(paste0("base_pp_tau_", name), function() plot(fits[[name]], "tau", prior = TRUE))
    vdiffr::expect_doppelganger(paste0("ggplot_pp_tau_", name), plot(fits[[name]], "tau", prior = TRUE, plot_type = "ggplot"))
  }

  # variance allocation
  for (name in names(fits)) {
    if (RoBMA:::.is_multilevel(fits[[name]])) {
      vdiffr::expect_doppelganger(paste0("base_pp_rho_", name), function() plot(fits[[name]], "rho", prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_rho_", name), plot(fits[[name]], "rho", prior = TRUE, plot_type = "ggplot"))
    }
  }

  # effect regression
  for (name in names(fits)) {
    if (RoBMA:::.is_mods(fits[[name]])) {
      for (pars in info[[name]][["mods"]]) {
        vdiffr::expect_doppelganger(paste0("base_pp_mods_", pars, "_", name), function() plot(fits[[name]], parameter_mods = pars, prior = TRUE))
        vdiffr::expect_doppelganger(paste0("ggplot_pp_scale_", pars, "_",, name), plot(fits[[name]], parameter_mods = pars, prior = TRUE, plot_type = "ggplot"))
      }
    }
  }

  # scale regression
  for (name in names(fits)) {
    if (RoBMA:::.is_scale(fits[[name]])) {
      for (pars in info[[name]][["scale"]]) {
        vdiffr::expect_doppelganger(paste0("base_pp_scale_", pars, "_", name), function() plot(fits[[name]], parameter_mods = pars, prior = TRUE))
        vdiffr::expect_doppelganger(paste0("ggplot_pp_scale_", pars, "_",, name), plot(fits[[name]], parameter_mods = pars, prior = TRUE, plot_type = "ggplot"))
      }
    }
  }

  # publication bias parameters
  for (name in names(fits)) {
    if (RoBMA:::.is_PET(fits[[name]])) {
      vdiffr::expect_doppelganger(paste0("base_pp_PET_", name), function() plot(fits[[name]], parameter = "PET", prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_PET_", name), plot(fits[[name]], parameter = "PET", prior = TRUE, plot_type = "ggplot"))
      # test complete PET-PEESE regression
      vdiffr::expect_doppelganger(paste0("base_pp_PETPEESE_", name), function() plot_PETPEESE(fits[[name]], prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_PETPEESE_", name), plot_PETPEESE(fits[[name]], prior = TRUE, plot_type = "ggplot"))
    }
    if (RoBMA:::.is_PEESE(fits[[name]])) {
      vdiffr::expect_doppelganger(paste0("base_pp_PEESE_", name), function() plot(fits[[name]], parameter = "PEESE", prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_PEESE_", name), plot(fits[[name]], parameter = "PEESE", prior = TRUE, plot_type = "ggplot"))
      # test complete PET-PEESE regression
      vdiffr::expect_doppelganger(paste0("base_pp_PETPEESE_", name), function() plot_PETPEESE(fits[[name]], prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_PETPEESE_", name), plot_PETPEESE(fits[[name]], prior = TRUE, plot_type = "ggplot"))
    }
    if (RoBMA:::.is_weightfunction(fits[[name]])) {
      vdiffr::expect_doppelganger(paste0("base_pp_omega_", name), function() {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(mar = oldpar[["mar"]]))
        par(mar = c(4, 4, 1, 4), mfrow = c(1, 4))
        plot(fits[[name]], parameter = "omega", prior = TRUE)
      })
      vdiffr::expect_doppelganger(paste0("ggplot_pp_omega_", name), function() {
        temp <- plot(fits[[name]], parameter = "omega", prior = TRUE, plot_type = "ggplot")
        gridExtra::grid.arrange(grobs = temp, nrow = 1)
      })
      # test complete weightfunction
      vdiffr::expect_doppelganger(paste0("base_pp_weightfunction_", name), function() plot_weightfunction(fits[[name]], parameter = "omega", prior = TRUE))
      vdiffr::expect_doppelganger(paste0("ggplot_pp_weightfunction_", name), plot_weightfunction(fits[[name]], parameter = "omega", prior = TRUE, plot_type = "ggplot"))
    }
  }

})


# test plot customization
