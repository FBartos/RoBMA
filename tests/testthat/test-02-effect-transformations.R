context("Effect-size transformations")

source(testthat::test_path("common-functions.R"))

skip_if_not_installed("metafor")

expect_effect_transform_matches_metafor <- function(input_measure, output_measure,
                                                    values, expected) {

  info   <- .effect_output_setup_measure(
    input_measure  = input_measure,
    output_measure = output_measure
  )
  actual <- .transform_effect_vector(values, info)

  expect_equal(
    actual,
    expected,
    tolerance = sqrt(.Machine$double.eps),
    info      = paste(input_measure, "to", output_measure, "matches metafor")
  )
}

expect_effect_jacobian_matches_numeric <- function(input_measure, values,
                                                   output_measure = NULL,
                                                   transform = NULL) {

  info           <- .effect_output_setup_measure(
    input_measure  = input_measure,
    output_measure = output_measure,
    transform      = transform
  )
  transformation <- info[["transformation"]]
  step           <- 1e-6
  expected       <- (
    transformation[["fun"]](values + step) -
    transformation[["fun"]](values - step)
  ) / (2 * step)

  expect_equal(
    transformation[["jac"]](values),
    expected,
    tolerance = 1e-5,
    info      = paste(input_measure, "to", info[["output_measure"]], "uses forward Jacobian")
  )
}

test_that("Fisher transformation preserves and validates its exact domain", {

  expect_identical(
    .cor_to_z(c(-1, 0, 1)),
    c(-Inf, 0, Inf)
  )
  expect_error(
    .cor_to_z(c(0, 1 + .Machine$double.eps)),
    "numeric, finite, and within \\[-1, 1\\]"
  )
  expect_error(
    .cor_to_z(c(0, NA_real_)),
    "numeric, finite, and within"
  )
})


test_that("effect-size measure transformations match metafor", {

  d_values <- c(-1, 0, 1)
  r_values <- c(-0.4, 0, 0.4)
  z_values <- metafor::transf.rtoz(r_values)
  o_values <- c(-1, 0, 1)

  expect_effect_transform_matches_metafor(
    input_measure  = "SMD",
    output_measure = "COR",
    values         = d_values,
    expected       = metafor::transf.dtorpb(d_values)
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "SMD",
    output_measure = "ZCOR",
    values         = d_values,
    expected       = metafor::transf.rtoz(metafor::transf.dtorpb(d_values))
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "SMD",
    output_measure = "OR",
    values         = d_values,
    expected       = metafor::transf.dtolnor.logis(d_values)
  )

  expect_effect_transform_matches_metafor(
    input_measure  = "COR",
    output_measure = "SMD",
    values         = r_values,
    expected       = metafor::transf.rpbtod(r_values)
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "COR",
    output_measure = "ZCOR",
    values         = r_values,
    expected       = metafor::transf.rtoz(r_values)
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "COR",
    output_measure = "OR",
    values         = r_values,
    expected       = metafor::transf.dtolnor.logis(metafor::transf.rpbtod(r_values))
  )

  expect_effect_transform_matches_metafor(
    input_measure  = "ZCOR",
    output_measure = "SMD",
    values         = z_values,
    expected       = metafor::transf.rpbtod(metafor::transf.ztor(z_values))
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "ZCOR",
    output_measure = "COR",
    values         = z_values,
    expected       = metafor::transf.ztor(z_values)
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "ZCOR",
    output_measure = "OR",
    values         = z_values,
    expected       = metafor::transf.dtolnor.logis(
      metafor::transf.rpbtod(metafor::transf.ztor(z_values))
    )
  )

  expect_effect_transform_matches_metafor(
    input_measure  = "OR",
    output_measure = "SMD",
    values         = o_values,
    expected       = metafor::transf.lnortod.logis(o_values)
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "OR",
    output_measure = "COR",
    values         = o_values,
    expected       = metafor::transf.dtorpb(metafor::transf.lnortod.logis(o_values))
  )
  expect_effect_transform_matches_metafor(
    input_measure  = "OR",
    output_measure = "ZCOR",
    values         = o_values,
    expected       = metafor::transf.rtoz(
      metafor::transf.dtorpb(metafor::transf.lnortod.logis(o_values))
    )
  )
})

test_that("effect-size transformations use forward Jacobians", {

  expect_effect_jacobian_matches_numeric(
    input_measure  = "SMD",
    output_measure = "COR",
    values         = c(-1, 0, 1)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "COR",
    output_measure = "SMD",
    values         = c(-0.4, 0, 0.4)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "COR",
    output_measure = "ZCOR",
    values         = c(-0.4, 0, 0.4)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "ZCOR",
    output_measure = "COR",
    values         = c(-0.5, 0, 0.5)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "SMD",
    output_measure = "OR",
    values         = c(-1, 0, 1)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "OR",
    output_measure = "SMD",
    values         = c(-1, 0, 1)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure = "RR",
    transform     = "EXP",
    values        = c(log(0.5), 0, log(2))
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "SMD",
    output_measure = "OR",
    transform      = "EXP",
    values         = c(-1, 0, 1)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "SMD",
    output_measure = "ZCOR",
    values         = c(-1, 0, 1)
  )
  expect_effect_jacobian_matches_numeric(
    input_measure  = "OR",
    output_measure = "COR",
    values         = c(-1, 0, 1)
  )
})

test_that("EXP is explicit for log-scale ratio output", {

  log_values <- c(log(0.5), 0, log(2))

  info <- .effect_output_setup_measure(
    input_measure  = "OR",
    output_measure = "OR",
    transform      = "EXP"
  )

  expect_equal(
    .transform_effect_vector(log_values, info),
    metafor::transf.exp.int(log_values, targs = list(tau2 = 0))
  )
  expect_equal(info[["label"]], "odds ratio")

  expect_error(
    .effect_output_setup_measure(
      input_measure = "SMD",
      transform     = "EXP"
    ),
    "EXP"
  )
  expect_error(
    .effect_output_setup_measure(
      input_measure = "RR",
      transform     = "LOG"
    ),
    "Unknown 'transform'"
  )
})

test_that("non-core measures are not converted across measures", {

  expect_error(
    .effect_output_setup_measure(
      input_measure  = "RR",
      output_measure = "SMD"
    ),
    "not available"
  )

  info <- .effect_output_setup_measure(
    input_measure = "RR",
    transform     = "EXP"
  )
  expect_equal(
    .transform_effect_vector(c(0, log(2)), info),
    metafor::transf.exp.int(c(0, log(2)), targs = list(tau2 = 0))
  )
})

test_that("scale-regression plot coefficients can use EXP", {

  data <- data.frame()
  attr(data, "measure") <- "SMD"
  object <- list(data = data)
  entry  <- list(
    component         = "scale",
    term              = "ni100",
    formula_parameter = "log_tau",
    role              = "fixed_coefficient"
  )

  info           <- .plot_output_setup(
    object          = object,
    parameter       = "log_tau_ni100",
    parameter_entry = entry,
    transform       = "EXP"
  )
  transformation <- info[["transformation"]]

  expect_true(info[["active"]])
  expect_identical(info[["label"]], "multiplicative scale")
  expect_equal(transformation[["fun"]](c(0, log(2))), c(1, 2))
  expect_equal(transformation[["inv"]](c(1, 2)), c(0, log(2)))
  expect_equal(transformation[["jac"]](c(0, log(2))), c(1, 2))

  entry[["term"]] <- "intercept"
  expect_error(
    .plot_output_setup(
      object          = object,
      parameter       = "log_tau_intercept",
      parameter_entry = entry,
      transform       = "EXP"
    ),
    "scale-regression coefficients"
  )
})

test_that("plot transformations use BayesTools forward Jacobian convention", {

  original_x    <- c(-1, 0, 1)
  transformed_x <- exp(original_x)
  info          <- .effect_output_setup_measure(
    input_measure = "RR",
    transform     = "EXP"
  )
  transformation <- .effect_plot_transformation(info)

  expect_equal(
    transformation[["jac"]](original_x),
    transformed_x,
    tolerance = sqrt(.Machine$double.eps)
  )
  expect_warning(
    density_y <- BayesTools:::.density.prior_transformation_y(
      x              = transformed_x,
      y              = rep(1, length(transformed_x)),
      transformation = transformation
    ),
    NA
  )
  expect_equal(
    density_y,
    1 / transformed_x,
    tolerance = sqrt(.Machine$double.eps)
  )

  positive_x         <- c(.5, 1, 2)
  log_transformation <- .log_plot_transformation()
  expect_equal(
    log_transformation[["jac"]](positive_x),
    1 / positive_x,
    tolerance = sqrt(.Machine$double.eps)
  )
  expect_warning(
    density_y <- BayesTools:::.density.prior_transformation_y(
      x              = log(positive_x),
      y              = rep(1, length(positive_x)),
      transformation = log_transformation
    ),
    NA
  )
  expect_equal(
    density_y,
    positive_x,
    tolerance = sqrt(.Machine$double.eps)
  )
})

test_that("precomputed posterior densities use EXP plot Jacobian", {

  raw_x          <- log(c(.5, 1, 2, 4))
  raw_y          <- c(.1, .4, .3, .05)
  transformed_x  <- exp(raw_x)
  sample         <- c(log(.6), log(1.2), log(2.5), log(3))
  point_mass_x   <- log(2.25)
  transformation <- .effect_plot_transformation(.effect_output_setup_measure(
    input_measure = "OR",
    transform     = "EXP"
  ))
  plot_data_samples_simple <- get(
    ".plot_data_samples.simple",
    envir    = asNamespace("BayesTools"),
    inherits = FALSE
  )

  density_methods <- list(
    qCMDE = list(method = "q_grid_cmde"),
    IWMDE = list(method = "iwmde")
  )

  for (density_method in names(density_methods)) {
    sample_with_density <- sample
    attr(sample_with_density, "prior_list") <- list(BayesTools::prior(
      "normal",
      parameters = list(mean = 0, sd = 1)
    ))
    attr(sample_with_density, "models_ind") <- rep(1L, length(sample_with_density))
    attr(sample_with_density, "posterior_density") <-
      BayesTools::posterior_density_attribute(
        x              = raw_x,
        y              = raw_y,
        method         = density_methods[[density_method]][["method"]],
        density_method = density_method,
        diagnostics    = list(estimator = density_methods[[density_method]][["method"]]),
        point_masses   = data.frame(x = point_mass_x, mass = .2)
      )

    plot_data <- plot_data_samples_simple(
      samples                  = list(mu = sample_with_density),
      parameter                = "mu",
      n_points                 = 20,
      transformation           = transformation,
      transformation_arguments = NULL,
      transformation_settings  = FALSE,
      density_method           = "precomputed"
    )

    expect_equal(plot_data[["density"]][["x"]], transformed_x)
    expect_equal(plot_data[["density"]][["y"]], raw_y / transformed_x)
    expect_equal(plot_data[["density"]][["samples"]], exp(sample))
    expect_equal(
      attr(plot_data[["density"]], "posterior_density_method"),
      density_methods[[density_method]][["method"]]
    )
    expect_equal(
      attr(plot_data[["density"]], "posterior_density_diagnostics")[["estimator"]],
      density_methods[[density_method]][["method"]]
    )
    expect_equal(plot_data[["points1"]][["x"]], exp(point_mass_x))
    expect_equal(plot_data[["points1"]][["y"]], .2)
  }
})

test_that("precomputed posterior densities use LOG plot Jacobian", {

  raw_x          <- c(.5, 1, 2, 4)
  raw_y          <- c(.1, .4, .3, .05)
  sample         <- c(.6, 1.2, 2.5, 3)
  point_mass_x   <- 2.25
  transformation <- .log_plot_transformation()
  plot_data_samples_simple <- get(
    ".plot_data_samples.simple",
    envir    = asNamespace("BayesTools"),
    inherits = FALSE
  )

  sample_with_density <- sample
  attr(sample_with_density, "prior_list") <- list(BayesTools::prior(
    "normal",
    parameters = list(mean = 0, sd = 1)
  ))
  attr(sample_with_density, "models_ind") <- rep(1L, length(sample_with_density))
  attr(sample_with_density, "posterior_density") <-
    BayesTools::posterior_density_attribute(
      x              = raw_x,
      y              = raw_y,
      method         = "q_grid_cmde",
      density_method = "qCMDE",
      point_masses   = data.frame(x = point_mass_x, mass = .2)
    )

  plot_data <- plot_data_samples_simple(
    samples                  = list(log_tau_intercept = sample_with_density),
    parameter                = "log_tau_intercept",
    n_points                 = 20,
    transformation           = transformation,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  expect_equal(plot_data[["density"]][["x"]], log(raw_x))
  expect_equal(plot_data[["density"]][["y"]], raw_y * raw_x)
  expect_equal(plot_data[["density"]][["samples"]], log(sample))
  expect_equal(plot_data[["points1"]][["x"]], log(point_mass_x))
  expect_equal(plot_data[["points1"]][["y"]], .2)
})

test_that("transformed brma_samples preserve posterior draw integration", {

  samples <- matrix(c(0, log(2), log(3)), ncol = 1)
  colnames(samples) <- "mu"

  info <- .effect_output_setup_measure(
    input_measure = "OR",
    transform     = "EXP"
  )

  out <- .new_effect_brma_samples(
    samples          = samples,
    n_chains         = 1,
    n_iter           = 3,
    title            = "Pooled Effect Size",
    effect_transform = info
  )

  expect_s3_class(out, "brma_samples")
  expect_equal(as.vector(as.matrix(out)), c(1, 2, 3))
  expect_equal(attr(out, "nchains"), 1)
  expect_equal(attr(out, "niter"), 3)

  skip_if_not_installed("posterior")
  draws <- posterior::as_draws_matrix(out)
  expect_equal(posterior::variables(draws), "mu")
  expect_equal(as.numeric(draws[, "mu"]), c(1, 2, 3))
})

test_that("marginal posterior samples are transformed without losing metadata", {

  marginal_sample <- list(intercept = c(0, log(2)))
  class(marginal_sample) <- c("marginal_posterior.formula", "marginal_posterior")
  attr(marginal_sample, "formula_parameter") <- "mu"

  samples <- list(mu_intercept = marginal_sample)
  info    <- .effect_output_setup_measure(
    input_measure = "OR",
    transform     = "EXP"
  )

  out <- .transform_marginal_samples_effect(
    samples          = samples,
    effect_transform = info
  )

  expect_s3_class(out[["mu_intercept"]], "marginal_posterior.formula")
  expect_equal(attr(out[["mu_intercept"]], "formula_parameter"), "mu")
  expect_equal(out[["mu_intercept"]][["intercept"]], c(1, 2))
})
