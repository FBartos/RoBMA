source(testthat::test_path("common-functions.R"))

# plot_prior() and print_prior() of the random-effect selectors of fitted
# objects. A quantity that is exactly one declared prior under the identity
# transform shows that prior; any other quantity shows the prior law
# BayesTools derives from the declared priors (parameter_prior_density(), the
# law plot(prior = TRUE) overlays) without a posterior. The expectations are
# analytic identities of the fits' own declared priors.

test_that("plot_prior and print_prior show the declared prior of a random SD", {

  skip_on_cran()
  fit   <- load_fit("brma.mv_block_mvn_random")
  prior <- attr(fit[["fit"]], "prior_list")[["mu__xREx__estimate_intercept"]]
  expect_true(BayesTools::is.prior(prior))

  plot  <- plot_prior(fit, "tau", plot_type = "ggplot")
  curve <- .random_prior_curve(plot)
  expect_gt(nrow(curve), 500L)
  .expect_random_prior_curve(
    curve, exp(BayesTools::lpdf(prior, curve[["x"]])), 1e-8
  )
  expect_equal(
    .random_prior_curve(plot_prior(
      fit, "tau", plot_type = "ggplot", standardized_coefficients = FALSE
    )),
    curve
  )
  expect_equal(
    .random_prior_curve(plot_prior(
      fit, "tau", plot_type = "ggplot", component = "random"
    )),
    curve
  )

  # base plots run and return nothing
  .random_prior_device(expect_null(plot_prior(fit, "tau")))

  printed <- capture.output(returned <- print_prior(fit, "tau"))
  expect_identical(returned, prior)
  expect_identical(printed[[1L]], "tau:")
  expect_true(all(
    paste0("  ", print(prior, silent = TRUE)) %in% printed
  ))
  expect_silent(print_prior(fit, "tau", silent = TRUE))

  # several selectors, and the whole specification of print_prior(fit, "random")
  plots <- plot_prior(fit, c("tau", "tau2"), plot_type = "ggplot")
  expect_named(plots, c("tau", "tau2"))
  expect_s3_class(
    print_prior(fit, "random", silent = TRUE), "prior_random"
  )
})

test_that("plot_prior draws the prior of a variance as the square of the SD prior", {

  skip_on_cran()
  fit   <- load_fit("brma.mv_block_mvn_random")
  prior <- attr(fit[["fit"]], "prior_list")[["mu__xREx__estimate_intercept"]]

  # the density of v = sigma^2 is f(sqrt(v)) / (2 sqrt(v))
  plot  <- plot_prior(fit, "tau2", plot_type = "ggplot")
  curve <- .random_prior_curve(plot)
  expect_gt(nrow(curve), 500L)
  .expect_random_prior_curve(
    curve,
    exp(BayesTools::lpdf(prior, sqrt(curve[["x"]]))) / (2 * sqrt(curve[["x"]])),
    1e-6
  )
  .random_prior_device(expect_null(plot_prior(fit, "tau2")))

  # the derivation line and the declared source prior
  printed <- capture.output(returned <- print_prior(fit, "tau2"))
  expect_identical(printed[[1L]], "tau2 = tau^2")
  expect_identical(printed[[2L]], "  tau:")
  expect_true(all(
    paste0("    ", print(prior, silent = TRUE)) %in% printed
  ))
  expect_s3_class(returned, "prior_linear_density")
  expect_silent(print_prior(fit, "tau2", silent = TRUE))

  # xlim, n_points and par_name select the plotted range, grid and label
  limited <- plot_prior(
    fit, "tau2", plot_type = "ggplot", xlim = c(0.05, 0.4), n_points = 40,
    par_name = "variance"
  )
  curve <- .random_prior_curve(limited)
  expect_gte(min(curve[["x"]]), 0.05)
  expect_lte(max(curve[["x"]]), 0.4)
  expect_lte(nrow(curve), 45L)
  .expect_random_prior_curve(
    curve,
    exp(BayesTools::lpdf(prior, sqrt(curve[["x"]]))) / (2 * sqrt(curve[["x"]])),
    1e-6, minimum = 30L
  )
  expect_match(
    paste(format(ggplot2::get_labs(limited)[["x"]]), collapse = ""),
    "variance", fixed = TRUE
  )
  .random_prior_device(expect_null(
    plot_prior(fit, "tau2", xlim = c(0.05, 0.4), n_points = 40, col = "red")
  ))

  # the evaluation grid of a prior law is not selected by x_seq
  expect_error(plot_prior(fit, "tau2", x_seq = c(0.1, 0.2)), "x_seq")
})

test_that("random selectors of gated allocations show prior laws with their point masses", {

  skip_on_cran()
  fit <- load_fit("BMA.mv_random_components")
  prior_list <- attr(fit[["fit"]], "prior_list")
  total   <- prior_list[["mu__xRE_ALLOCx_heterogeneity__allocation_sd"]]
  weights <- prior_list[["mu__xRE_ALLOCx_heterogeneity__weight"]]
  selection <- function(selector) {
    .brma_parameter_select_entry(fit, selector)[["selection"]]
  }
  # two components with independent inclusion gates of probability 0.5

  # the total SD is 0 when no component is included (probability 0.25) and
  # the scale prior otherwise: the mixed law has the atom and is drawn with
  # an arrow at 0
  law <- BayesTools::parameter_prior_density(fit[["fit"]], selection("tau_total"))
  expect_equal(law[["points"]][["x"]], 0)
  expect_equal(law[["points"]][["p"]], 0.25, tolerance = 1e-12)
  layers <- .random_prior_layers(plot_prior(fit, "tau_total", plot_type = "ggplot"))
  expect_length(layers[["arrows"]], 1L)
  expect_equal(layers[["arrows"]][[1L]][["x"]], 0)
  expect_gt(layers[["arrows"]][[1L]][["yend"]], 0)
  printed <- capture.output(print_prior(fit, "tau_total"))
  expect_match(
    printed[[1L]], "tau_total = sigma_total * sqrt(sum of I_c * w[c] over the components c)",
    fixed = TRUE
  )
  expect_true(any(grepl("^  w:$", printed)))
  expect_true(any(grepl("^  p_study:$", printed)))
  expect_true(any(grepl("^  p_observation:$", printed)))
  expect_true(any(grepl("Spike(0.5)", printed, fixed = TRUE)))
  expect_true(all(paste0("    ", print(weights, silent = TRUE)) %in% printed))

  # the share of a component is a mixture of the atoms at 0 and 1 and the
  # Beta(1, 1) marginal, given an active component: each of the three
  # configurations (only this component, only the other one, both) has
  # probability 1/3
  law <- BayesTools::parameter_prior_density(fit[["fit"]], selection("tau2_prop(study)"))
  expect_equal(law[["points"]][["x"]], c(0, 1))
  expect_equal(law[["points"]][["p"]], rep(1 / 3, 2), tolerance = 1e-12)
  expect_equal(law[["density"]][["mass"]], 1 / 3, tolerance = 1e-12)
  printed <- capture.output(print_prior(fit, "tau2_prop(study)"))
  expect_identical(
    printed[[1L]],
    paste0(
      "tau2_prop(study) = I_study * w[1] / sum of I_c * w[c] over the ",
      "components c, the share of the total variance of allocation ",
      "'heterogeneity' of the active components"
    )
  )
  expect_identical(
    grep(":$", printed, value = TRUE), c("  w:", "  p_study:", "  p_observation:")
  )
  layers <- .random_prior_layers(plot_prior(fit, "tau2_prop(study)", plot_type = "ggplot"))
  expect_length(layers[["arrows"]], 2L)
  expect_equal(
    sort(vapply(layers[["arrows"]], function(arrow) arrow[["x"]], numeric(1))),
    c(0, 1)
  )
  curve <- .random_prior_curve(plot_prior(fit, "tau2_prop(study)", plot_type = "ggplot"))
  expect_equal(curve[["y"]], rep(1 / 3, nrow(curve)), tolerance = 1e-6)

  # a component SD is drawn the same way, with the gate priors in its derivation
  layers <- .random_prior_layers(plot_prior(fit, "study: tau", plot_type = "ggplot"))
  expect_length(layers[["arrows"]], 1L)
  expect_equal(layers[["arrows"]][[1L]][["x"]], 0)
  printed <- capture.output(print_prior(fit, "study: tau"))
  expect_match(
    printed[[1L]], "study: tau = sigma_total * I_study * sqrt(w[1])", fixed = TRUE
  )
  expect_true(any(grepl("^  p_study:$", printed)))
  .random_prior_device(expect_null(plot_prior(fit, "tau2_prop(study)")))
  .random_prior_device(expect_null(plot_prior(fit, "tau_total", xlim = c(0, 2))))
})

test_that("random selectors of a nested gated allocation derive from its allocation chain", {

  skip_on_cran()
  fit <- load_fit("lane_nested_total")
  prior_list <- attr(fit[["fit"]], "prior_list")
  total <- prior_list[["mu__xRE_ALLOCx_heterogeneity__allocation_sd"]]
  included <- 0.5

  # the total is the scale prior times the gate: the atom at 0 and the scaled
  # scale prior, whose derivation does not involve the split weights
  law <- BayesTools::parameter_prior_density(
    fit[["fit"]], .brma_parameter_select_entry(fit, "tau_total")[["selection"]]
  )
  expect_equal(law[["points"]][["x"]], 0)
  expect_equal(law[["points"]][["p"]], 1 - included, tolerance = 1e-12)
  plot <- plot_prior(fit, "tau_total", plot_type = "ggplot")
  curve <- .random_prior_curve(plot)
  .expect_random_prior_curve(
    curve[curve[["x"]] > 0, ], included * exp(BayesTools::lpdf(total, curve[["x"]][curve[["x"]] > 0])),
    1e-6
  )
  expect_length(.random_prior_layers(plot)[["arrows"]], 1L)
  printed <- capture.output(print_prior(fit, "tau_total"))
  expect_identical(printed[[1L]], "tau_total = sigma_total * I_component_1")
  expect_identical(grep(":$", printed, value = TRUE), c("  sigma_total:", "  p_component_1:"))

  # a component SD follows the gate of the root and the split weights
  printed <- capture.output(print_prior(fit, "esid_study: tau"))
  expect_identical(
    printed[[1L]], "esid_study: tau = sigma_total * I_component_1 * sqrt(w[1])"
  )
  expect_identical(
    grep(":$", printed, value = TRUE), c("  sigma_total:", "  p_component_1:", "  w:")
  )
  .random_prior_device(expect_null(plot_prior(fit, "esid_study: tau")))

  # a share of the split depends on the split weights alone
  printed <- capture.output(print_prior(fit, "tau2_prop(study)"))
  expect_identical(
    printed[[1L]],
    "tau2_prop(study) = w[2], the share of the total variance of allocation 'component_1_split'"
  )
  expect_identical(grep(":$", printed, value = TRUE), "  w:")
})

test_that("random selectors without a prior density are refused with the target classes", {

  skip_on_cran()
  fit <- load_fit("lane_random_correlation_direct")

  expect_error(
    plot_prior(fit, "tau(x)"),
    "prior of the random-effect quantity 'tau(x)' is unavailable",
    fixed = TRUE, class = "RoBMA_hypothesis_target"
  )
  for (call in list(
    quote(plot_prior(fit, "tau(x)")),
    quote(print_prior(fit, "tau(x)")),
    quote(plot_prior(fit, "rho(intercept,x)", standardized_coefficients = FALSE))
  )) {
    error <- tryCatch(eval(call), error = identity)
    expect_s3_class(error, "RoBMA_hypothesis_target")
    expect_s3_class(error, "RoBMA_hypothesis_unavailable")
  }

  # the correlation has its exact prior on the fitted scale
  expect_error(
    plot_prior(fit, "rho(intercept,x)", standardized_coefficients = FALSE),
    "standardized_coefficients = TRUE", fixed = TRUE
  )
  .random_prior_device(expect_null(plot_prior(fit, "rho(intercept,x)")))
  expect_gt(
    nrow(.random_prior_curve(plot_prior(fit, "rho(intercept,x)", plot_type = "ggplot"))),
    100L
  )
})

test_that("only_priors objects point to the random-effect specification", {

  skip_on_cran()
  data <- data.frame(
    yi    = c(0.10, 0.25, 0.15, 0.30, 0.05),
    sei   = c(0.2, 0.25, 0.22, 0.28, 0.17),
    study = c("s1", "s1", "s2", "s2", "s3")
  )
  random <- brma.mv(
    yi = yi, V = diag(data[["sei"]]^2), random = ~ 1 | study, data = data,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE
  )

  for (call in list(
    quote(plot_prior(random, "tau")),
    quote(print_prior(random, "tau", silent = TRUE))
  )) {
    error <- tryCatch(eval(call), error = identity)
    expect_identical(class(error), c("simpleError", "error", "condition"))
    expect_match(
      conditionMessage(error), "The specified parameter 'tau' is not available."
    )
    expect_match(
      conditionMessage(error),
      "print_prior(object, parameter = \"random\")", fixed = TRUE
    )
  }
  expect_s3_class(
    print_prior(random, parameter = "random", silent = TRUE), "prior_random"
  )

  # objects without random-effect priors give the original message only
  plain <- BMA(
    yi = yi, sei = sei, data = data, measure = "SMD", only_priors = TRUE
  )
  error <- tryCatch(print_prior(plain, "tau2", silent = TRUE), error = identity)
  expect_match(conditionMessage(error), "is not available.")
  expect_false(grepl("random", conditionMessage(error), fixed = TRUE))
})
