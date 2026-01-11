
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
    set.seed(1)
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

test_that("Customization of prior and posterior distributions for model parameters", {

  ### simple parameter
  name <- "bcg_meta-analysis"

  # no prior
  vdiffr::expect_doppelganger("baseplot_pp_no_prior", function() plot(fits[[name]], "mu", prior = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_no_prior", plot(fits[[name]], "mu", prior = FALSE, plot_type = "ggplot"))

  # change range
  vdiffr::expect_doppelganger("baseplot_pp_range", function() plot(fits[[name]], "mu", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 5)))
  vdiffr::expect_doppelganger("ggplot_pp_range", plot(fits[[name]], "mu", prior = TRUE, plot_type = "ggplot", xlim = c(-1, 1), ylim = c(0, 5)))

  # change aesthetics
  vdiffr::expect_doppelganger("baseplot_pp_aesthetics", function() plot(fits[[name]], "mu", prior = TRUE, lwd = 3, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))
  vdiffr::expect_doppelganger("ggplot_pp_aesthetics", plot(fits[[name]], "mu", prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))


  ### PET-PEESE
  name <- "dat.lehmann2018-PET"
  set.seed(1)

  # no prior
  vdiffr::expect_doppelganger("baseplot_pp_PETPEESE_no_prior", function() plot_PETPEESE(fits[[name]], prior = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_PETPEESE_no_prior", plot_PETPEESE(fits[[name]], prior = FALSE, plot_type = "ggplot"))

  # change range
  vdiffr::expect_doppelganger("baseplot_pp_PETPEESE_range", function() plot_PETPEESE(fits[[name]], prior = TRUE, xlim = c(0, 0.5), ylim = c(-3, 3)))
  vdiffr::expect_doppelganger("ggplot_pp_PETPEESE_range", plot_PETPEESE(fits[[name]], prior = TRUE, plot_type = "ggplot", xlim = c(0, 0.5), ylim = c(-3, 3)))

  # change aesthetics
  vdiffr::expect_doppelganger("baseplot_pp_PETPEESE_aesthetics", function() plot_PETPEESE(fits[[name]], prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                               dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  vdiffr::expect_doppelganger("ggplot_pp_PETPEESE_aesthetics", plot_PETPEESE(fits[[name]], prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                  dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))


  ### weight function
  name <- "dat.lehmann2018-3PSM"
  set.seed(1)

  # no prior
  vdiffr::expect_doppelganger("baseplot_pp_weightfunction_no_prior", function() plot_weightfunction(fits[[name]], prior = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_weightfunction_no_prior", plot_weightfunction(fits[[name]], prior = FALSE, plot_type = "ggplot"))

  # change range
  vdiffr::expect_doppelganger("baseplot_pp_weightfunction_range", function() plot_weightfunction(fits[[name]], prior = TRUE, rescale_p_values = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_weightfunction_range", plot_weightfunction(fits[[name]], prior = TRUE, plot_type = "ggplot", rescale_p_values = FALSE))

  # change aesthetics
  vdiffr::expect_doppelganger("baseplot_pp_weightfunction_aesthetics", function() plot_weightfunction(fits[[name]], prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                        dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  vdiffr::expect_doppelganger("ggplot_pp_weightfunction_aesthetics", plot_weightfunction(fits[[name]], prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                           dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

