context("Prior and posterior plots")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)


test_that("plot.brma uses parameter x-axis labels by default", {

  temp_fit <- fits[["bcg_meta-analysis"]]

  x_label <- function(plot) {

    return(plot$scales$get_scales("x")$name)
  }

  expect_identical(
    x_label(plot(temp_fit, "mu", plot_type = "ggplot")),
    "Effect Size"
  )
  expect_identical(
    x_label(plot(temp_fit, "tau", plot_type = "ggplot")),
    "Heterogeneity"
  )
  expect_identical(
    x_label(plot(temp_fit, "mu", plot_type = "ggplot", xlab = "Custom Label")),
    "Custom Label"
  )
})


test_that("Prior and posterior distributions for brma.norm models", {

  ### simple meta-analysis ----
  name <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  ### effect size
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, "mu", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot"))

  # change range
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_range"), function() plot(temp_fit, "mu", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 5)))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_range"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot", xlim = c(-1, 1), ylim = c(0, 5)))

  # change aesthetics
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_aesthetics"), function() plot(temp_fit, "mu", prior = TRUE, lwd = 3, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_aesthetics"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))

  ### heterogeneity
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_tau_pp_prior"), function() plot(temp_fit, "tau", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-tau_ggplot_pp_prior"), plot(temp_fit, "tau", prior = TRUE, plot_type = "ggplot"))

  ### meta-regression (continuous)
  name <- "bcg_meta-regression"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_orig"), function() plot(temp_fit, "mu", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_std"), function() plot(temp_fit, "mu", prior = TRUE, standardized_coefficients = TRUE))

  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression (categorical: dummy)
  name <- "bcg_meta-regression2"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression (categorical: meandif)
  name <- "bcg_meta-regression2b"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression with interactions (categorical: dummy)
  name <- "bcg_meta-regression3"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_orig"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_std"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression with interactions (categorical: meandif)
  name <- "bcg_meta-regression3b"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_orig"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_std"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE, standardized_coefficients = TRUE))

  ### location-scale model
  name <- "bangertdrowns2004_location-scale"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-scale0_baseplot_orig"), function() plot(temp_fit, parameter_scale = "intercept", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-scale0_baseplot_std"), function() plot(temp_fit, parameter_scale = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-scale1_baseplot_orig"), function() plot(temp_fit, parameter_scale = "ni100", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-scale1_baseplot_std"), function() plot(temp_fit, parameter_scale = "ni100", prior = TRUE, standardized_coefficients = TRUE))

  ### between-study heterogeneity and multilevel models
  name <- "konstantopoulos2011_3lvl"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-rho_baseplot"), function() plot(temp_fit, parameter = "rho", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot"), function() plot(temp_fit, parameter = "tau", prior = TRUE))
})

test_that("Transformed effect-size plots render", {

  name     <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  .with_temp_plot_device(
    expect_silent(plot(temp_fit, "mu", transform = "EXP", plot_type = "base"))
  )
  expect_true(.is_ggplot(
    plot(temp_fit, "mu", transform = "EXP", plot_type = "ggplot")
  ))
})

test_that("Prior and posterior plots transform effect-size axis", {

  skip_if_not_full_visuals("Effect-size transform snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-mu_transform_exp_comparison"), function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    plot(temp_fit, "mu", plot_type = "base", main = "log RR")
    plot(temp_fit, "mu", transform = "EXP", plot_type = "base", main = "RR")
  })
})

test_that("Prior and posterior distributions for bPET / bPEESE objects", {

  skip_if_not_full_visuals("PET/PEESE prior-posterior snapshots are visual-gallery coverage.")

  ### PET model
  name <- "dat.lehmann2018-PET"
  temp_fit <- fits[[name]]
  set.seed(1)

  # no prior
  expect_vdiffr_snapshot("baseplot_pp_PETPEESE_no_prior", function() plot_pet_peese(temp_fit, prior = FALSE))
  expect_vdiffr_snapshot("ggplot_pp_PETPEESE_no_prior", plot_pet_peese(temp_fit, prior = FALSE, plot_type = "ggplot"))

  # change range
  expect_vdiffr_snapshot("baseplot_pp_PETPEESE_range", function() plot_pet_peese(temp_fit, prior = TRUE, xlim = c(0, 0.5), ylim = c(-3, 3)))
  expect_vdiffr_snapshot("ggplot_pp_PETPEESE_range", plot_pet_peese(temp_fit, prior = TRUE, plot_type = "ggplot", xlim = c(0, 0.5), ylim = c(-3, 3)))

  # change aesthetics
  expect_vdiffr_snapshot("baseplot_pp_PETPEESE_aesthetics", function() plot_pet_peese(temp_fit, prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                          dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  expect_vdiffr_snapshot("ggplot_pp_PETPEESE_aesthetics", plot_pet_peese(temp_fit, prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                             dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

test_that("Prior and posterior distributions for bselmodel objects", {

  skip_if_not_full_visuals("Selection prior-posterior snapshots are visual-gallery coverage.")


  ### weight function
  name <- "dat.lehmann2018-3PSM"
  temp_fit <- fits[[name]]
  set.seed(1)

  # no prior
  expect_vdiffr_snapshot("baseplot_pp_weightfunction_no_prior", function() plot_weightfunction(temp_fit, prior = FALSE))
  expect_vdiffr_snapshot("ggplot_pp_weightfunction_no_prior", plot_weightfunction(temp_fit, prior = FALSE, plot_type = "ggplot"))

  # change range
  expect_vdiffr_snapshot("baseplot_pp_weightfunction_range", function() plot_weightfunction(temp_fit, prior = TRUE, rescale_p_values = FALSE))
  expect_vdiffr_snapshot("ggplot_pp_weightfunction_range", plot_weightfunction(temp_fit, prior = TRUE, plot_type = "ggplot", rescale_p_values = FALSE))

  # change aesthetics
  expect_vdiffr_snapshot("baseplot_pp_weightfunction_aesthetics", function() plot_weightfunction(temp_fit, prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                                      dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  expect_vdiffr_snapshot("ggplot_pp_weightfunction_aesthetics", plot_weightfunction(temp_fit, prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                         dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

test_that("Weightfunction plot supports observed p-value rug", {

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  .with_temp_plot_device(
    expect_silent(plot_weightfunction(fits[[name]], show_data = TRUE))
  )
  expect_true(.is_ggplot(
    plot_weightfunction(
      fits[[name]],
      show_data = TRUE,
      plot_type = "ggplot",
      dots_data = list(color = "red", alpha = .5, linewidth = .4, rug_side = "top")
    )
  ))
  expect_silent(
    plot_weightfunction(
      fits[[name]],
      show_data        = TRUE,
      rescale_p_values = FALSE,
      dots_data        = list(col = "blue", lwd = 1, side = "bottom")
    )
  )

  rug_plot <- plot_weightfunction(
    fits[[name]],
    show_data        = TRUE,
    plot_type        = "ggplot",
    rescale_p_values = FALSE
  )
  rug_data <- ggplot_geom_layer_data(rug_plot, "GeomRug")

  expect_equal(
    rug_data[["p"]],
    .weightfunction_observed_p_values(fits[[name]]),
    tolerance = sqrt(.Machine$double.eps)
  )
})

test_that("Prior and posterior distributions for BMA.norm objects", {

  skip_if_not_full_visuals("BMA prior-posterior snapshots are visual-gallery coverage.")

  name <- "dat.lehmann2018_BMA.norm_mods"
  temp_fit <- fits[[name]]
  set.seed(1)

  # effect
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, "mu", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot"))

  # moderation
  expect_vdiffr_snapshot(paste0(name, "-mods_ggplot_pp_no_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mods_ggplot_pp_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = TRUE, plot_type = "ggplot"))

  # heterogeneity
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_pp_no_prior"), function() plot(temp_fit, "tau", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_pp_prior"), function() plot(temp_fit, "tau", prior = TRUE))

})

test_that("Prior and posterior distributions for RoBMA objects", {

  skip_if_not_full_visuals("RoBMA prior-posterior snapshots are visual-gallery coverage.")

  name <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  temp_fit <- fits[[name]]
  set.seed(1)

  # effect
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))

  # moderation
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = TRUE, plot_type = "ggplot"))

  # heterogeneity
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_pp_no_prior"), function() plot(temp_fit, "tau", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-rho_baseplot_pp_no_prior"), function() plot(temp_fit, "rho", prior = TRUE))

})
