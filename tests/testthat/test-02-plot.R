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


test_that("plot.brma component disambiguates shared location-scale terms", {

  skip_if_missing_fits("dat.lehmann2018_RoBMA_3lvl_mods_scale")

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  fit <- fits[["dat.lehmann2018_RoBMA_3lvl_mods_scale"]]

  expect_error(
    plot(fit, parameter = "Preregistered", plot_type = "ggplot"),
    "ambiguous"
  )

  out <- plot(
    fit,
    parameter = "Preregistered",
    component = "scale",
    plot_type = "ggplot"
  )
  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "log_tau_Preregistered")

  out <- plot(
    fit,
    parameter = "Preregistered",
    component = "mods",
    plot_type = "ggplot"
  )
  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu_Preregistered")

  out <- plot(
    fit,
    parameter = "Preregistered",
    component = "location",
    plot_type = "ggplot"
  )
  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu_Preregistered")
})


test_that("plot.brma uses KDE by default", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter = "mu",
    plot_type = "ggplot"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["dots"]][["density_method"]], "KDE")
  expect_null(attr(captured[["samples"]][["mu"]], "posterior_density"))
})


test_that("plot.brma forwards attached qCMDE posterior density", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter          = "mu",
    plot_type          = "ggplot",
    density_method     = "qCMDE",
    density_control    = list(n_points = 20, max_samples = 20)
  )

  posterior_density <- attr(captured[["samples"]][["mu"]], "posterior_density")

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu")
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_equal(posterior_density[["status"]], "ok")
  expect_equal(posterior_density[["density_method"]], "qCMDE")
  expect_true(all(is.finite(posterior_density[["x"]])))
  expect_true(all(is.finite(posterior_density[["y"]])))
})


test_that("plot.brma forwards attached IWMDE posterior density", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter          = "mu",
    plot_type          = "ggplot",
    density_method     = "IWMDE",
    density_control    = list(n_points = 20, max_samples = 50)
  )

  posterior_density <- attr(captured[["samples"]][["mu"]], "posterior_density")

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu")
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_equal(posterior_density[["status"]], "ok")
  expect_equal(posterior_density[["density_method"]], "IWMDE")
  expect_equal(posterior_density[["diagnostics"]][["estimator"]], "iwmde")
})


test_that("plot.brma aligns qCMDE density to plotted coefficient scale", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  plot(
    fits[["bcg_meta-regression"]],
    parameter_mods    = "year",
    plot_type         = "ggplot",
    density_method    = "qCMDE",
    density_control   = list(n_points = 20, max_samples = 20)
  )

  plotted_samples   <- captured[["samples"]][[captured[["parameter"]]]]
  posterior_density <- attr(plotted_samples, "posterior_density")

  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_false(is.null(posterior_density[["diagnostics"]][["plot_scale_transform"]]))
  expect_lte(
    max(abs(posterior_density[["x"]])),
    max(abs(as.numeric(plotted_samples))) + .05
  )
})


test_that("plot.brma drops qCMDE density when scale alignment is not affine", {

  posterior_density <- list(
    x            = seq(-1, 1, length.out = 20),
    y            = stats::dnorm(seq(-1, 1, length.out = 20)),
    diagnostics  = list(),
    point_masses = data.frame(x = numeric(), mass = numeric())
  )

  out <- .plot_brma_align_iwmde_density(
    posterior_density = posterior_density,
    raw_samples       = seq(-1, 1, length.out = 20),
    plotted_samples   = seq(-1, 1, length.out = 20)^2
  )

  expect_null(out)
})


test_that("plot.brma forwards attached qCMDE density for PET and PEESE parameters", {

  .skip_if_missing_raw_fits(c("dat.lehmann2018-PET", "dat.lehmann2018-PEESE"))

  for (case in list(
    list(fit = "dat.lehmann2018-PET", parameter = "PET"),
    list(fit = "dat.lehmann2018-PEESE", parameter = "PEESE")
  )) {
    captured <- NULL
    testthat::local_mocked_bindings(
      plot_posterior = function(samples, parameter, ...) {
        captured <<- list(samples = samples, parameter = parameter, dots = list(...))
        return(structure(list(), class = "mock_plot"))
      },
      .package = "BayesTools"
    )

    out <- plot(
      load_fit(case[["fit"]], validate = FALSE),
      parameter         = case[["parameter"]],
      plot_type         = "ggplot",
      density_method    = "qCMDE",
      density_control   = list(n_points = 20, max_samples = 20)
    )

    density_source <- if (!is.null(captured[["samples"]][[case[["parameter"]]]])) {
      case[["parameter"]]
    } else {
      "bias"
    }
    posterior_density <- attr(captured[["samples"]][[density_source]], "posterior_density")

    expect_s3_class(out, "mock_plot")
    expect_equal(captured[["dots"]][["density_method"]], "precomputed")
    expect_equal(posterior_density[["status"]], "ok")
    expect_equal(posterior_density[["density_method"]], "qCMDE")
  }
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
