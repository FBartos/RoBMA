context("Prior input handling for brma.norm")

# Test data for prior specification tests
test_data <- data.frame(
  effect     = c(0.10, 0.25, 0.15, 0.30, 0.05),
  std_err    = sqrt(c(0.04, 0.06, 0.05, 0.08, 0.03)),
  n          = c(50L, 75L, 60L, 40L, 90L),
  mod_cont   = c(1.5, 2.3, 1.8, 3.1, 0.9),
  mod_factor = factor(c("A", "B", "A", "B", "A")),
  scale_var  = c(0.5, 1.0, 0.8, 1.2, 0.6),
  cluster    = c("g1", "g1", "g2", "g2", "g3"),
  stringsAsFactors = FALSE
)


# ============================================================================
# Tests for measure specification and UISD
# ============================================================================

test_that("Default priors scale correctly with known UISD for each measure", {

  skip_on_cran()

  # Test all known UISD measures
  measures_uisd <- list(SMD = sqrt(2), ZCOR = 1, OR = 2, RR = 2, HR = 2)

  for (measure in names(measures_uisd)) {
    result <- brma.norm(
      yi = effect, sei = std_err, data = test_data,
      measure = measure, only_priors = TRUE
    )[["priors"]]

    expected_effect_sd <- measures_uisd[[measure]] * RoBMA.get_option("default_UISD.effect")
    expected_het_sd    <- measures_uisd[[measure]] * RoBMA.get_option("default_UISD.heterogeneity")

    expect_true(BayesTools::is.prior(result$outcome$mu))
    expect_true(BayesTools::is.prior(result$outcome$tau))
    expect_equal(result$outcome$mu$parameters$sd, expected_effect_sd)
    expect_equal(result$outcome$tau$parameters$sd, expected_het_sd)
    expect_equal(result$outcome$tau$truncation$lower, 0)
  }
})


test_that("UISD is estimated from sample sizes for GEN measure", {

  skip_on_cran()

  result <- brma.norm(
    yi = effect, sei = std_err, ni = n, data = test_data,
    measure = "GEN", only_priors = TRUE
  )[["priors"]]

  expected_uisd <- estimate_unit_information_sd(sei = test_data$std_err, ni = test_data$n)
  expected_effect_sd <- expected_uisd * RoBMA.get_option("default_UISD.effect")

  expect_equal(result$outcome$mu$parameters$sd, expected_effect_sd, tolerance = 1e-6)
})


test_that("GEN measure without ni throws error", {

  skip_on_cran()

  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, measure = "GEN", only_priors = TRUE),
    regexp = "ni|unit_information_sd|UISD"
  )
})


test_that("Invalid measure throws error", {

  skip_on_cran()

  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, measure = "INVALID", only_priors = TRUE),
    regexp = "measure"
  )
})


# ============================================================================
# Tests for prior_unit_information_sd and estimate_unit_information_sd
# ============================================================================

test_that("Manual prior_unit_information_sd is used and overrides known UISD", {

  skip_on_cran()

  custom_uisd <- 3.0

  # Should work for GEN measure
  result_gen <- brma.norm(
    yi = effect, sei = std_err, ni = n, data = test_data,
    measure = "GEN", prior_unit_information_sd = custom_uisd, only_priors = TRUE
  )[["priors"]]

  # Should override SMD's known UISD
  result_smd <- brma.norm(
    yi = effect, sei = std_err, ni = n, data = test_data,
    measure = "SMD", prior_unit_information_sd = custom_uisd, only_priors = TRUE
  )[["priors"]]

  expected_sd <- custom_uisd * RoBMA.get_option("default_UISD.effect")
  expect_equal(result_gen$outcome$mu$parameters$sd, expected_sd)
  expect_equal(result_smd$outcome$mu$parameters$sd, expected_sd)
})


test_that("estimate_unit_information_sd computes correctly and validates inputs", {

  skip_on_cran()

  sei <- c(0.2, 0.3, 0.25, 0.15)
  ni  <- c(50, 30, 40, 80)

  # Correct computation
  result <- estimate_unit_information_sd(sei = sei, ni = ni)
  expected <- sqrt(sum(ni) / sum(sei^(-2)))
  expect_equal(result, expected, tolerance = 1e-10)

  # Invalid inputs
  expect_error(estimate_unit_information_sd(sei = c(-0.2, 0.3), ni = c(50, 30)), regexp = "sei")
  expect_error(estimate_unit_information_sd(sei = c(0, 0.3), ni = c(50, 30)), regexp = "sei")
  expect_error(estimate_unit_information_sd(sei = c(0.2, 0.3), ni = c(-50, 30)), regexp = "ni")
  expect_error(estimate_unit_information_sd(sei = c(0.2, 0.3, 0.4), ni = c(50, 30)), regexp = "length|ni")
  expect_error(estimate_unit_information_sd(sei = c(0.2, NA), ni = c(50, 30)), regexp = "sei")
})


test_that("Invalid prior_unit_information_sd throws error", {

  skip_on_cran()

  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, prior_unit_information_sd = -1, only_priors = TRUE),
    regexp = "prior_unit_information_sd"
  )
  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, prior_unit_information_sd = 0, only_priors = TRUE),
    regexp = "prior_unit_information_sd"
  )
})


# ============================================================================
# Tests for rescale_priors
# ============================================================================

test_that("rescale_priors correctly scales prior distributions", {

  skip_on_cran()

  result_base <- brma.norm(yi = effect, sei = std_err, data = test_data, measure = "SMD", only_priors = TRUE)[["priors"]]

  for (scale_factor in c(0.5, 2)) {
    result_rescaled <- brma.norm(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD", rescale_priors = scale_factor, only_priors = TRUE
    )[["priors"]]
    expect_equal(result_rescaled$outcome$mu$parameters$sd, result_base$outcome$mu$parameters$sd * scale_factor)
    expect_equal(result_rescaled$outcome$tau$parameters$sd, result_base$outcome$tau$parameters$sd * scale_factor)
  }

  # Invalid values
  expect_error(brma.norm(yi = effect, sei = std_err, data = test_data, rescale_priors = -1, only_priors = TRUE))
  expect_error(brma.norm(yi = effect, sei = std_err, data = test_data, rescale_priors = 0, only_priors = TRUE))
})


# ============================================================================
# Tests for custom priors
# ============================================================================

test_that("Custom prior_effect and prior_heterogeneity work with rescaling", {

  skip_on_cran()

  custom_effect <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5))
  custom_het    <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.3), truncation = list(0, Inf))

  result <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    prior_effect = custom_effect, prior_heterogeneity = custom_het,
    only_priors = TRUE
  )[["priors"]]

  expect_equal(result$outcome$mu$parameters$sd, 0.5)
  expect_equal(result$outcome$tau$parameters$sd, 0.3)

  # With rescaling
  result_rescaled <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    prior_effect = custom_effect, prior_heterogeneity = custom_het, rescale_priors = 2, only_priors = TRUE
  )[["priors"]]
  expect_equal(result_rescaled$outcome$mu$parameters$sd,  0.5 * 2)
  expect_equal(result_rescaled$outcome$tau$parameters$sd, 0.3 * 2)
})


test_that("Different prior distribution types work", {

  skip_on_cran()

  # Spike prior
  spike <- BayesTools::prior("spike", parameters = list(location = 0))
  result_spike <- brma.norm(yi = effect, sei = std_err, data = test_data, prior_effect = spike, measure = "OR", only_priors = TRUE)[["priors"]]
  expect_true(BayesTools::is.prior.point(result_spike$outcome$mu))

  # t/Cauchy prior (rescalable via scale parameter)
  cauchy <- BayesTools::prior("cauchy", parameters = list(location = 0, scale = 0.707))
  result_cauchy <- brma.norm(yi = effect, sei = std_err, data = test_data, prior_effect = cauchy, measure = "RR", only_priors = TRUE)[["priors"]]
  expect_equal(result_cauchy$outcome$mu$distribution, "t")

  result_cauchy_rescaled <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    prior_effect = cauchy, rescale_priors = 2, measure = "HR", only_priors = TRUE
  )[["priors"]]
  expect_equal(result_cauchy_rescaled$outcome$mu$parameters$scale, 0.707 * 2)
})


# ============================================================================
# Tests for informed priors
# ============================================================================

test_that("Informed priors work for medicine field with rescaling", {

  skip_on_cran()

  # Test general priors
  result <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "OR", prior_informed_field = "medicine", only_priors = TRUE
  )[["priors"]]
  expect_equal(result$outcome$mu, BayesTools::prior_informed(
    name = "cochrane", parameter = "effect", type = "logOR"))
  expect_equal(result$outcome$tau, BayesTools::prior_informed(
    name = "cochrane", parameter = "heterogeneity", type = "logOR"))

  # With subfield
  result_subfield <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_informed_field = "medicine", prior_informed_subfield = "neonatal",
    only_priors = TRUE
  )[["priors"]]
  expect_equal(result_subfield$outcome$mu, BayesTools::prior_informed(
    name = "neonatal", parameter = "effect", type = "smd"))
  expect_equal(result_subfield$outcome$tau, BayesTools::prior_informed(
    name = "neonatal", parameter = "heterogeneity", type = "smd"))

  # With rescaling
  result_base <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_informed_field = "medicine", only_priors = TRUE
  )[["priors"]]
  result_rescaled <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_informed_field = "medicine", rescale_priors = 2, only_priors = TRUE
  )[["priors"]]

  scale_param <- if ("sd" %in% names(result_base$outcome$mu$parameters)) "sd" else "scale"
  expect_equal(
    result_rescaled$outcome$mu$parameters[[scale_param]],
    result_base$outcome$mu$parameters[[scale_param]] * 2
  )
})


test_that("Conflicting prior specifications throw errors", {

  skip_on_cran()

  # Invalid field
  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, prior_informed_field = "invalid", only_priors = TRUE),
    regexp = "prior_informed_field"
  )

  # Both UISD and informed priors
  expect_error(
    brma.norm(
      yi = effect, sei = std_err, data = test_data,
      prior_unit_information_sd = 1.5, prior_informed_field = "medicine", only_priors = TRUE
    ),
    regexp = "prior_unit_information_sd|prior_informed_field|mutually exclusive"
  )
})


# ============================================================================
# Tests for heterogeneity allocation prior (study_ids)
# ============================================================================

test_that("Heterogeneity allocation prior works correctly", {

  skip_on_cran()

  # Default Beta(1,1) when study_ids specified
  result <- brma.norm(
    yi = effect, sei = std_err, study_ids = cluster, data = test_data,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]
  expect_true("rho" %in% names(result$outcome))
  expect_equal(result$outcome$rho$distribution, "beta")
  expect_equal(result$outcome$rho$parameters$alpha, 1)
  expect_equal(result$outcome$rho$parameters$beta, 1)

  # Custom prior
  custom_prior <- BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2))
  result_custom <- brma.norm(
    yi = effect, sei = std_err, study_ids = cluster, data = test_data,
    prior_heterogeneity_allocation = custom_prior, measure = "OR", only_priors = TRUE
  )[["priors"]]
  expect_equal(result_custom$outcome$rho$parameters$alpha, 2)

  # No rho when study_ids not specified
  result_no_ids <- brma.norm(yi = effect, sei = std_err, data = test_data, measure = "SMD", only_priors = TRUE)
  expect_false("rho" %in% names(result_no_ids$outcome))
})


# ============================================================================
# Tests for moderator and scale priors
# ============================================================================

test_that("Moderator priors are assigned correctly", {

  skip_on_cran()

  result <- brma.norm(
    yi = effect, sei = std_err, mods = ~ mod_cont + mod_factor,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(!is.null(result$mods))
  expect_true("intercept" %in% names(result$mods))
  expect_false("mu" %in% names(result$outcome))  # mu forwarded to intercept

  # Intercept inherits from custom effect prior
  custom_prior <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.8))
  result_custom <- brma.norm(
    yi = effect, sei = std_err, mods = ~ mod_cont, data = test_data,
    prior_effect = custom_prior, measure = "SMD", only_priors = TRUE
  )[["priors"]]
  expect_equal(result_custom$mods$intercept$parameters$sd, 0.8)

  # Moderator priors scale with UISD
  expected_mod_sd <- sqrt(2) * RoBMA.get_option("default_UISD.mods")
  expect_equal(result$mods$mod_cont$parameters$sd, expected_mod_sd)

  # Custom prior_mods
  custom_mod <- list(mod_cont = BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.3)))
  result_mod <- brma.norm(
    yi = effect, sei = std_err, mods = ~ mod_cont, data = test_data,
    prior_mods = custom_mod, measure = "SMD", only_priors = TRUE
  )[["priors"]]
  expect_equal(result_mod$mods$mod_cont$parameters$sd, 0.3)
})


test_that("Scale priors are assigned correctly", {

  skip_on_cran()

  result <- brma.norm(
    yi = effect, sei = std_err, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(!is.null(result$scale))
  expect_true("intercept" %in% names(result$scale))
  expect_false("tau" %in% names(result$outcome))  # tau forwarded to intercept

  # Intercept inherits from custom heterogeneity prior
  custom_prior <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5), truncation = list(0, Inf))
  result_custom <- brma.norm(
    yi = effect, sei = std_err, scale = ~ scale_var, data = test_data,
    prior_heterogeneity = custom_prior, prior_unit_information_sd = 10, only_priors = TRUE
  )[["priors"]]
  expect_equal(result_custom$scale$intercept$parameters$sd, 0.5)

  # Scale priors use fixed UISD (not measure-dependent)
  result_smd  <- brma.norm(yi = effect, sei = std_err, scale = ~ scale_var, data = test_data, measure = "SMD",  only_priors = TRUE)[["priors"]]
  result_zcor <- brma.norm(yi = effect, sei = std_err, scale = ~ scale_var, data = test_data, measure = "ZCOR", only_priors = TRUE)[["priors"]]

  if (!is.null(result_smd$scale$scale_var) && !is.null(result_zcor$scale$scale_var)) {
    expect_equal(result_smd$scale$scale_var$parameters$sd, result_zcor$scale$scale_var$parameters$sd)
  }
})


test_that("Both mods and scale priors work together", {

  skip_on_cran()

  result <- brma.norm(
    yi = effect, sei = std_err, mods = ~ mod_cont + mod_factor, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(!is.null(result$mods))
  expect_true(!is.null(result$scale))
  expect_true("intercept" %in% names(result$mods))
  expect_true("intercept" %in% names(result$scale))
  expect_false("mu" %in% names(result$outcome))
  expect_false("tau" %in% names(result$outcome))
})


# ============================================================================
# Tests for contrast
# ============================================================================

test_that("set_contrast_factor_predictors options work", {

  skip_on_cran()

  for (contrast in c("treatment", "meandif", "orthonormal")) {
    result <- brma.norm(
      yi = effect, sei = std_err, mods = ~ mod_factor, data = test_data,
      measure = "SMD", set_contrast_factor_predictors = contrast, only_priors = TRUE
    )[["priors"]]
    expect_true(any(grepl(contrast, class(result$mods[["mod_factor"]]))))
  }

  expect_error(
    brma.norm(
      yi = effect, sei = std_err, mods = ~ mod_factor, data = test_data,
      set_contrast_factor_predictors = "invalid", only_priors = TRUE
    ),
    regexp = "set_contrast_factor_predictors"
  )
})
