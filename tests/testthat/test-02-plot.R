context("Prior and posterior plots")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)


test_that("Prior and posterior distributions for brma.norm models", {

  ### simple meta-analysis ----
  name <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  ### effect size
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, "mu", prior = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot"))

  # change range
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_range"), function() plot(temp_fit, "mu", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 5)))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_range"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot", xlim = c(-1, 1), ylim = c(0, 5)))

  # change aesthetics
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_aesthetics"), function() plot(temp_fit, "mu", prior = TRUE, lwd = 3, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_aesthetics"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))

  ### heterogeneity
  vdiffr::expect_doppelganger(paste0(name, "-tau_baseplot_tau_pp_prior"), function() plot(temp_fit, "tau", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-tau_ggplot_pp_prior"), plot(temp_fit, "tau", prior = TRUE, plot_type = "ggplot"))

  ### meta-regression (continuous)
  name <- "bcg_meta-regression"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_orig"), function() plot(temp_fit, "mu", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_std"), function() plot(temp_fit, "mu", prior = TRUE, standardized_coefficients = TRUE))

  vdiffr::expect_doppelganger(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression (categorical: dummy)
  name <- "bcg_meta-regression2"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression (categorical: meandif)
  name <- "bcg_meta-regression2b"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression with interactions (categorical: dummy)
  name <- "bcg_meta-regression3"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-reg0_baseplot_orig"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg0_baseplot_std"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg1_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg1_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg2_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg2_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg3_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg3_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression with interactions (categorical: meandif)
  name <- "bcg_meta-regression3b"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-reg0_baseplot_orig"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg0_baseplot_std"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg1_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg1_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg2_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg2_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg3_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-reg3_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE, standardized_coefficients = TRUE))

  ### location-scale model
  name <- "bangertdrowns2004_location-scale"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-scale0_baseplot_orig"), function() plot(temp_fit, parameter_scale = "intercept", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-scale0_baseplot_std"), function() plot(temp_fit, parameter_scale = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-scale1_baseplot_orig"), function() plot(temp_fit, parameter_scale = "ni100", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-scale1_baseplot_std"), function() plot(temp_fit, parameter_scale = "ni100", prior = TRUE, standardized_coefficients = TRUE))

  ### between-study heterogeneity and multilevel models
  name <- "konstantopoulos2011_3lvl"
  temp_fit <- fits[[name]]

  vdiffr::expect_doppelganger(paste0(name, "-rho_baseplot"), function() plot(temp_fit, parameter = "rho", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-tau_baseplot"), function() plot(temp_fit, parameter = "tau", prior = TRUE))
})

test_that("Prior and posterior distributions for bPET / bPEESE objects", {

  ### PET model
  name <- "dat.lehmann2018-PET"
  temp_fit <- fits[[name]]
  set.seed(1)

  # no prior
  vdiffr::expect_doppelganger("baseplot_pp_PETPEESE_no_prior", function() plot_PETPEESE(temp_fit, prior = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_PETPEESE_no_prior", plot_PETPEESE(temp_fit, prior = FALSE, plot_type = "ggplot"))

  # change range
  vdiffr::expect_doppelganger("baseplot_pp_PETPEESE_range", function() plot_PETPEESE(temp_fit, prior = TRUE, xlim = c(0, 0.5), ylim = c(-3, 3)))
  vdiffr::expect_doppelganger("ggplot_pp_PETPEESE_range", plot_PETPEESE(temp_fit, prior = TRUE, plot_type = "ggplot", xlim = c(0, 0.5), ylim = c(-3, 3)))

  # change aesthetics
  vdiffr::expect_doppelganger("baseplot_pp_PETPEESE_aesthetics", function() plot_PETPEESE(temp_fit, prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                          dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  vdiffr::expect_doppelganger("ggplot_pp_PETPEESE_aesthetics", plot_PETPEESE(temp_fit, prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                             dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

test_that("Prior and posterior distributions for bselmod objects", {


  ### weight function
  name <- "dat.lehmann2018-3PSM"
  temp_fit <- fits[[name]]
  set.seed(1)

  # no prior
  vdiffr::expect_doppelganger("baseplot_pp_weightfunction_no_prior", function() plot_weightfunction(temp_fit, prior = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_weightfunction_no_prior", plot_weightfunction(temp_fit, prior = FALSE, plot_type = "ggplot"))

  # change range
  vdiffr::expect_doppelganger("baseplot_pp_weightfunction_range", function() plot_weightfunction(temp_fit, prior = TRUE, rescale_p_values = FALSE))
  vdiffr::expect_doppelganger("ggplot_pp_weightfunction_range", plot_weightfunction(temp_fit, prior = TRUE, plot_type = "ggplot", rescale_p_values = FALSE))

  # change aesthetics
  vdiffr::expect_doppelganger("baseplot_pp_weightfunction_aesthetics", function() plot_weightfunction(temp_fit, prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                                      dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  vdiffr::expect_doppelganger("ggplot_pp_weightfunction_aesthetics", plot_weightfunction(temp_fit, prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                         dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

test_that("Prior and posterior distributions for BMA.norm objects", {

  name <- "dat.lehmann2018_BMA.norm_mods"
  temp_fit <- fits[[name]]
  set.seed(1)

  # effect
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, "mu", prior = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot"))

  # moderation
  vdiffr::expect_doppelganger(paste0(name, "-mods_ggplot_pp_no_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0(name, "-mods_ggplot_pp_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = TRUE, plot_type = "ggplot"))

  # heterogeneity
  vdiffr::expect_doppelganger(paste0(name, "-tau_baseplot_pp_no_prior"), function() plot(temp_fit, "tau", prior = FALSE))
  vdiffr::expect_doppelganger(paste0(name, "-tau_baseplot_pp_prior"), function() plot(temp_fit, "tau", prior = TRUE))

})

test_that("Prior and posterior distributions for RoBMA objects", {

  name <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  temp_fit <- fits[[name]]
  set.seed(1)

  # effect
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  vdiffr::expect_doppelganger(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))

  # moderation
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = FALSE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = TRUE, plot_type = "ggplot"))

  # heterogeneity
  vdiffr::expect_doppelganger(paste0(name, "-tau_baseplot_pp_no_prior"), function() plot(temp_fit, "tau", prior = TRUE))
  vdiffr::expect_doppelganger(paste0(name, "-rho_baseplot_pp_no_prior"), function() plot(temp_fit, "rho", prior = TRUE))

})
