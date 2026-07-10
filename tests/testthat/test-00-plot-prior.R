context("Prior plotting")

test_data <- data.frame(
  effect     = c(0.10, 0.25, 0.15, 0.30, 0.05),
  std_err    = sqrt(c(0.04, 0.06, 0.05, 0.08, 0.03)),
  mod_cont   = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_factor = factor(c("A", "B", "A", "B", "A")),
  scale_var  = c(0.5, 1.0, 0.8, 1.2, 0.6),
  stringsAsFactors = FALSE
)

.with_temp_plot_device <- function(expr) {

  file <- tempfile(fileext = ".png")
  grDevices::png(filename = file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)

  force(expr)
}

.is_ggplot <- function(x) {

  inherits(x, "ggplot")
}

.ggplot_x_range <- function(x) {

  build_data <- ggplot2::ggplot_build(x)[["data"]]
  x_values   <- unlist(lapply(build_data, function(layer) layer[["x"]]))

  range(x_values, finite = TRUE)
}

test_that("plot_prior plots outcome priors from only_priors objects", {

  skip_on_cran()

  priors <- BMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_priors = TRUE
  )

  expect_true(.is_ggplot(plot_prior(priors, parameter = "mu",  plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(priors, parameter = "tau", plot_type = "ggplot")))

  plots <- plot_prior(priors, parameter = c("mu", "tau"), plot_type = "ggplot")
  expect_named(plots, c("mu", "tau"))
  expect_true(all(vapply(plots, .is_ggplot, logical(1))))

  .with_temp_plot_device(
    expect_silent(plot_prior(priors, parameter = "mu"))
  )
})

test_that("plot_prior selects moderator and scale priors", {

  skip_on_cran()

  reg_priors <- BMA(
    yi = effect, sei = std_err,
    mods = ~ mod_cont + mod_factor,
    data = test_data, measure = "SMD", only_priors = TRUE
  )

  expect_true(.is_ggplot(plot_prior(reg_priors, parameter = "mu", plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(reg_priors, parameter_mods = "mod_cont", plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(reg_priors, parameter = "mod_factor", plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(reg_priors, parameter = "mod_factor", standardized_coefficients = FALSE, plot_type = "ggplot")))

  brma_factor_priors <- brma(
    yi = effect, sei = std_err,
    mods = ~ mod_factor,
    data = test_data, measure = "SMD",
    set_contrast_factor_predictors = "meandif",
    only_priors = TRUE
  )

  .with_temp_plot_device(
    expect_silent(plot_prior(brma_factor_priors, parameter_mods = "mod_factor"))
  )

  expect_error(
    plot_prior(reg_priors, parameter_mods = "missing", plot_type = "ggplot"),
    regexp = "not available"
  )
  expect_error(
    plot_prior(reg_priors, parameter_mods = "mod_cont", standardized_coefficients = "no", plot_type = "ggplot"),
    regexp = "standardized_coefficients"
  )

  default_plot      <- plot_prior(reg_priors, parameter_mods = "mod_cont", plot_type = "ggplot")
  standardized_plot <- plot_prior(reg_priors, parameter_mods = "mod_cont", standardized_coefficients = TRUE,  plot_type = "ggplot")
  raw_plot          <- plot_prior(reg_priors, parameter_mods = "mod_cont", standardized_coefficients = FALSE, plot_type = "ggplot")

  expect_equal(.ggplot_x_range(default_plot), .ggplot_x_range(standardized_plot))
  expect_gt(diff(.ggplot_x_range(raw_plot)), diff(.ggplot_x_range(standardized_plot)))

  scale_priors <- suppressWarnings(BMA(
    yi = effect, sei = std_err,
    scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  ))

  expect_true(.is_ggplot(plot_prior(scale_priors, parameter = "tau", plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(scale_priors, component = "scale", plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(scale_priors, parameter_scale = "scale_var", plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(scale_priors, parameter_scale = "scale_var", standardized_coefficients = FALSE, plot_type = "ggplot")))

  both_priors <- suppressWarnings(BMA(
    yi = effect, sei = std_err,
    mods = ~ scale_var, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  ))

  expect_error(
    plot_prior(both_priors, parameter = "scale_var", plot_type = "ggplot"),
    regexp = "component"
  )
  expect_true(.is_ggplot(plot_prior(
    both_priors, parameter = "scale_var", component = "mods",
    plot_type = "ggplot"
  )))
  expect_true(.is_ggplot(plot_prior(
    both_priors, parameter = "scale_var", component = "location",
    plot_type = "ggplot"
  )))
  expect_true(.is_ggplot(plot_prior(
    both_priors, parameter = "scale_var", component = "scale",
    plot_type = "ggplot"
  )))
})

test_that("plot_prior selects brma.mv location and component-scale priors", {

  skip_on_cran()

  mv_data <- data.frame(
    effect    = test_data[["effect"]],
    std_err   = test_data[["std_err"]],
    mod_cont  = test_data[["mod_cont"]],
    scale_var = test_data[["scale_var"]],
    study     = c("s1", "s1", "s2", "s2", "s3"),
    site      = c("a", "b", "a", "b", "a"),
    stringsAsFactors = FALSE
  )

  mv_priors <- brma.mv(
    yi                        = effect,
    V                         = diag(std_err^2),
    mods                      = ~ mod_cont,
    random                    = list(study = ~ 1 | study, site = ~ 1 | site),
    scale                     = list(study = ~ scale_var, site = ~ scale_var),
    data                      = mv_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  catalog <- .brma_parameter_catalog(mv_priors)
  expect_equal(
    nrow(.brma_random_parameter_bundle(mv_priors)[["specs"]]),
    0L
  )
  expect_true(any(
    catalog[["parameter"]] == "mu_mod_cont" &
      catalog[["source"]] == "location"
  ))
  expect_true(any(
    catalog[["parameter"]] == "log_tau_study_scale_var" &
      catalog[["formula_parameter"]] == "log_tau_study"
  ))

  expect_true(BayesTools::is.prior(print_prior(
    mv_priors, parameter = "mod_cont", component = "location", silent = TRUE
  )))
  expect_true(BayesTools::is.prior(print_prior(
    mv_priors,
    parameter = "log_tau_study_scale_var",
    component = "scale",
    silent    = TRUE
  )))
  expect_error(
    print_prior(mv_priors, parameter = "scale_var", component = "scale", silent = TRUE),
    "ambiguous"
  )
  expect_error(
    print_prior(mv_priors, component = "scale", silent = TRUE),
    "component-specific scale"
  )

  selected <- print_prior(mv_priors, silent = TRUE)
  expect_true("mu_mod_cont" %in% names(selected))
  expect_true("log_tau_study_scale_var" %in% names(selected))
  expect_true("log_tau_site_scale_var" %in% names(selected))
  expect_s3_class(selected[["random"]], "prior_random")
  expect_s3_class(
    print_prior(mv_priors, parameter = "random", silent = TRUE),
    "prior_random"
  )
  selected_random <- print_prior(
    mv_priors,
    parameter = c("mod_cont", "random"),
    silent    = TRUE
  )
  expect_true(BayesTools::is.prior(selected_random[["mod_cont"]]))
  expect_s3_class(selected_random[["random"]], "prior_random")
  expect_error(
    print_prior(mv_priors, parameter = "random", component = "scale", silent = TRUE),
    "not available"
  )
  expect_error(
    print_prior(mv_priors, parameter = "random", component = "nope", silent = TRUE),
    "should be one of"
  )

  expect_true(.is_ggplot(plot_prior(
    mv_priors,
    parameter                 = "mod_cont",
    component                 = "location",
    standardized_coefficients = FALSE,
    plot_type                 = "ggplot"
  )))
  expect_true(.is_ggplot(plot_prior(
    mv_priors,
    parameter                 = "log_tau_study_scale_var",
    component                 = "scale",
    standardized_coefficients = FALSE,
    plot_type                 = "ggplot"
  )))
})

test_that("print_prior respects regular terms named random", {

  skip_on_cran()

  mv_data <- data.frame(
    effect  = test_data[["effect"]],
    std_err = test_data[["std_err"]],
    random  = test_data[["mod_cont"]],
    study   = c("s1", "s1", "s2", "s2", "s3"),
    stringsAsFactors = FALSE
  )

  mv_priors <- brma.mv(
    yi                        = effect,
    V                         = diag(std_err^2),
    mods                      = ~ random,
    random                    = ~ 1 | study,
    data                      = mv_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  auto_selected <- print_prior(mv_priors, parameter = "random", silent = TRUE)
  mods_selected <- print_prior(
    mv_priors,
    parameter = "random", component = "mods",
    silent    = TRUE
  )

  expect_true(BayesTools::is.prior(auto_selected))
  expect_true(BayesTools::is.prior(mods_selected))
  expect_false(inherits(auto_selected, "prior_random"))
  expect_false(inherits(mods_selected, "prior_random"))
  expect_s3_class(
    print_prior(mv_priors, parameter = "mu_random", silent = TRUE),
    "prior"
  )
  expect_error(
    print_prior(mv_priors, parameter = "random", component = "scale", silent = TRUE),
    "not available"
  )
  expect_true(.is_ggplot(plot_prior(
    mv_priors,
    parameter = "random", component = "mods",
    plot_type = "ggplot"
  )))
})

test_that("plot_prior omits absent scalar scale priors for random brma.mv", {

  skip_on_cran()

  mv_data <- data.frame(
    effect  = test_data[["effect"]],
    std_err = test_data[["std_err"]],
    study   = c("s1", "s1", "s2", "s2", "s3"),
    site    = c("a", "b", "a", "b", "a"),
    stringsAsFactors = FALSE
  )

  mv_priors <- brma.mv(
    yi                        = effect,
    V                         = diag(std_err^2),
    random                    = list(study = ~ 1 | study, site = ~ 1 | site),
    data                      = mv_data,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  catalog <- .brma_parameter_catalog(mv_priors)
  expect_false("tau" %in% catalog[["parameter"]])

  selected <- print_prior(mv_priors, silent = TRUE)
  expect_true("mu_intercept" %in% names(selected))
  expect_false("tau" %in% names(selected))
  expect_s3_class(selected[["random"]], "prior_random")
  full_printed <- capture.output(print_prior(mv_priors))
  expect_true("random:" %in% full_printed)
  expect_error(
    print_prior(mv_priors, component = "scale", silent = TRUE),
    "does not contain scale priors"
  )
})

test_that("plot_prior supports direct prior objects", {

  skip_on_cran()

  prior_object <- prior("normal", parameters = list(mean = 0, sd = 1))

  expect_true(.is_ggplot(plot_prior(prior_object, plot_type = "ggplot")))
})

test_that("print_prior prints selected priors", {

  skip_on_cran()

  priors <- BMA(
    yi = effect, sei = std_err,
    mods = ~ mod_cont + mod_factor,
    data = test_data, measure = "SMD", only_priors = TRUE
  )

  expect_true(BayesTools::is.prior(print_prior(priors, parameter = "mu", silent = TRUE)))
  expect_true(BayesTools::is.prior(print_prior(priors, parameter_mods = "mod_cont", silent = TRUE)))
  expect_true(BayesTools::is.prior(print_prior(priors, parameter = "mod_factor", silent = TRUE)))

  both_priors <- suppressWarnings(BMA(
    yi = effect, sei = std_err,
    mods = ~ scale_var, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  ))
  expect_true(BayesTools::is.prior(print_prior(
    both_priors, parameter = "scale_var", component = "scale",
    silent = TRUE
  )))
  expect_true(BayesTools::is.prior(print_prior(
    both_priors, component = "scale", silent = TRUE
  )))

  full_output <- capture.output(full_selected <- print_prior(priors, silent = TRUE))
  expect_identical(full_output, character(0))
  expect_named(full_selected, c("mu_intercept", "mu_mod_cont", "mu_mod_factor", "tau"))
  expect_true(all(vapply(full_selected, BayesTools::is.prior, logical(1))))

  full_printed   <- capture.output(print_prior(priors))
  parameter_rows <- grep("^(mu_|tau:)", full_printed, value = TRUE)
  expect_equal(sub(":.*$", "", parameter_rows), names(full_selected))

  printed <- capture.output(print_prior(priors, parameter = "tau"))
  expect_equal(printed[1:3], c(
    "tau:",
    "  alternative:",
    "    (1/2) * Normal(0, 0.35)[0, Inf]"
  ))
  expect_true("  null:" %in% printed)

  selected <- print_prior(priors, parameter = c("mu", "tau"), silent = TRUE)
  expect_named(selected, c("mu", "tau"))
  expect_true(all(vapply(selected, BayesTools::is.prior, logical(1))))

  prior_object <- prior("normal", parameters = list(mean = 0, sd = 1))
  expect_true(BayesTools::is.prior(print_prior(prior_object, silent = TRUE)))
})

test_that("plot_prior handles publication-bias prior components", {

  skip_on_cran()

  priors <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_priors = TRUE
  )

  expect_error(
    plot_prior(priors, parameter = "bias", plot_type = "ggplot"),
    regexp = "mixes weight-function and PET-PEESE"
  )
  omega_output <- capture.output(
    omega_plot <- plot_prior(priors, parameter = "omega", plot_type = "ggplot")
  )
  expect_identical(omega_output, character(0))
  expect_true(.is_ggplot(omega_plot))
  expect_true(.is_ggplot(plot_prior(priors, parameter = "PET",   plot_type = "ggplot")))
  expect_true(.is_ggplot(plot_prior(priors, parameter = "PEESE", plot_type = "ggplot")))

  expect_true(BayesTools::is.prior(print_prior(priors, parameter = "bias",  silent = TRUE)))
  expect_true(BayesTools::is.prior(print_prior(priors, parameter = "omega", silent = TRUE)))
  expect_true(BayesTools::is.prior(print_prior(priors, parameter = "PET",   silent = TRUE)))
})

test_that("only_priors objects print and plot via prior methods", {

  skip_on_cran()

  priors <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", model_type = "PP", only_priors = TRUE
  )

  expect_s3_class(priors, "only_priors.brma")

  printed <- capture.output(print(priors))
  expect_true(any(grepl("^mu:", printed)))
  expect_true(any(grepl("^bias:", printed)))

  .with_temp_plot_device(
    expect_silent(plot(priors))
  )
  expect_true(.is_ggplot(plot(priors, plot_type = "ggplot")))
})
