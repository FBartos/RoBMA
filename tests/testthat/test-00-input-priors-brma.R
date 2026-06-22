context("Prior input handling for brma.norm")

source(testthat::test_path("helper-contracts.R"))

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

test_that("Default priors scale with known UISD for each measure", {

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


test_that("GEN measure without ni is rejected", {

  skip_on_cran()

  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, measure = "GEN", only_priors = TRUE),
    regexp = "ni|unit_information_sd|UISD"
  )
})

test_that("GEN measure accepts explicit default moderator prior without UISD", {

  skip_on_cran()

  explicit_effect <- BayesTools::prior(
    distribution = "normal",
    parameters   = list(mean = 0, sd = 10)
  )
  explicit_heterogeneity <- BayesTools::prior(
    distribution = "normal",
    parameters   = list(mean = 0, sd = 1),
    truncation   = list(0, Inf)
  )
  explicit_mods <- BayesTools::prior_factor(
    distribution = "normal",
    parameters   = list(mean = 0, sd = 2),
    contrast     = "independent"
  )

  result <- brma.norm(
    yi = effect, sei = std_err, mods = ~ 0 + mod_factor,
    prior_effect        = explicit_effect,
    prior_heterogeneity = explicit_heterogeneity,
    prior_mods          = explicit_mods,
    data = test_data, measure = "GEN", only_priors = TRUE
  )[["priors"]]

  expect_equal(result[["mods"]][["mod_factor"]][["parameters"]][["sd"]], 2)
  expect_equal(result[["outcome"]][["tau"]][["parameters"]][["sd"]], 1)
})

test_that("GEN measure rejects incomplete ni for UISD defaults", {

  skip_on_cran()

  expect_error(
    brma.norm(
      yi          = effect,
      sei         = std_err,
      ni          = c(50L, NA, 60L, 40L, 90L),
      data        = test_data,
      measure     = "GEN",
      only_priors = TRUE
    ),
    regexp = "ni.*complete|unit information"
  )
})


test_that("Invalid measure is rejected", {

  skip_on_cran()

  expect_error(
    brma.norm(yi = effect, sei = std_err, data = test_data, measure = "INVALID", only_priors = TRUE),
    regexp = "measure"
  )
})


test_that("Normal constructors require explicit measure", {

  skip_on_cran()

  constructors <- list(
    brma      = brma,
    BMA       = BMA,
    RoBMA     = RoBMA,
    bPET      = bPET,
    bPEESE    = bPEESE,
    bselmodel = bselmodel
  )

  for (constructor in constructors) {
    expect_error(
      do.call(constructor, list(
        yi                        = test_data$effect,
        sei                       = test_data$std_err,
        prior_unit_information_sd = 1,
        only_priors               = TRUE
      )),
      regexp = "requires explicit 'measure'"
    )
  }
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


test_that("estimate_unit_information_sd computes UISD and validates inputs", {

  skip_on_cran()

  sei <- c(0.2, 0.3, 0.25, 0.15)
  ni  <- c(50, 30, 40, 80)

  # Correct computation
  result <- estimate_unit_information_sd(sei = sei, ni = ni)
  expected <- sqrt(sum(ni) / sum(sei^(-2)))
  expect_equal(result, expected, tolerance = 1e-10)

  expect_error_cases(list(
    list(
      label  = "negative sei",
      expr   = quote(estimate_unit_information_sd(sei = c(-0.2, 0.3), ni = c(50, 30))),
      regexp = "sei"
    ),
    list(
      label  = "zero sei",
      expr   = quote(estimate_unit_information_sd(sei = c(0, 0.3), ni = c(50, 30))),
      regexp = "sei"
    ),
    list(
      label  = "negative ni",
      expr   = quote(estimate_unit_information_sd(sei = c(0.2, 0.3), ni = c(-50, 30))),
      regexp = "ni"
    ),
    list(
      label  = "length mismatch",
      expr   = quote(estimate_unit_information_sd(sei = c(0.2, 0.3, 0.4), ni = c(50, 30))),
      regexp = "length|ni"
    ),
    list(
      label  = "missing sei",
      expr   = quote(estimate_unit_information_sd(sei = c(0.2, NA), ni = c(50, 30))),
      regexp = "sei"
    )
  ))
})


test_that("Invalid prior_unit_information_sd is rejected", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "negative UISD",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "GEN",
        prior_unit_information_sd = -1, only_priors = TRUE
      )),
      regexp = "prior_unit_information_sd"
    ),
    list(
      label  = "zero UISD",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "GEN",
        prior_unit_information_sd = 0, only_priors = TRUE
      )),
      regexp = "prior_unit_information_sd"
    )
  ))
})


# ============================================================================
# Tests for constructor aliases and bias prior defaults
# ============================================================================

test_that("brma.norm alias is not overwritten by bselmodel", {

  skip_on_cran()

  result_alias <- brma(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_data = TRUE
  )
  result       <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_data = TRUE
  )

  expect_identical(class(result_alias), class(result))
  expect_s3_class(result, "brma.norm")
  expect_false("bselmodel" %in% class(result))
  expect_false("prior_bias" %in% names(formals(brma.norm)))
})


test_that("Explicit NULL prior_bias uses defaults for single bias models", {

  skip_on_cran()

  result_sel <- bselmodel(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_bias = NULL, only_priors = TRUE
  )[["priors"]]
  expect_true(BayesTools::is.prior.weightfunction(result_sel$outcome$bias))

  result_pet <- bPET(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_bias = NULL, only_priors = TRUE
  )[["priors"]]
  expect_true(BayesTools::is.prior.PET(result_pet$outcome$bias))
})

test_that("bselmodel keeps only steps shortcut for default weightfunction", {

  skip_on_cran()

  result <- bselmodel(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", steps = c(0.05, 0.10), only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.weightfunction(result$outcome$bias))
  expect_equal(result$outcome$bias$side, "one-sided")
  expect_equal(result$outcome$bias$steps, c(0.05, 0.10))
  expect_equal(result$outcome$bias$weights$type, "cumulative")
  expect_equal(result$outcome$bias$weights$alpha, rep(RoBMA.get_option("default_bias_weightfunction.alpha"), 3))

  prior_bias <- prior_weightfunction("two-sided", c(0.05), wf_fixed(c(1, 0.5)))
  result <- bselmodel(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_bias = prior_bias, only_priors = TRUE
  )[["priors"]]

  expect_equal(result$outcome$bias$side, "two-sided")
  expect_equal(result$outcome$bias$weights$type, "fixed")
  expect_equal(result$outcome$bias$weights$omega, c(1, 0.5))
})

test_that("RoBMA accepts fixed zero weightfunctions at prior setup", {

  skip_on_cran()

  zero_prior <- BayesTools::prior_weightfunction(
    "one-sided", c(.05), BayesTools::wf_fixed(c(1, 0))
  )
  valid_prior <- BayesTools::prior_weightfunction(
    "one-sided", c(.05), BayesTools::wf_fixed(c(1, .5))
  )

  expect_equal(zero_prior$weights$omega, c(1, 0))
  expect_silent(
    bselmodel(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD", prior_bias = zero_prior, only_priors = TRUE
    )
  )
  expect_silent(
    bselmodel(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD",
      prior_bias = BayesTools::prior_bias(selection = zero_prior),
      only_priors = TRUE
    )
  )
  expect_silent(
    RoBMA(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD",
      prior_bias = list(valid_prior, zero_prior),
      prior_unit_information_sd = 1,
      only_priors = TRUE,
      silent = TRUE
    )
  )
})

test_that("selection backend consumes BayesTools omega p-order directly", {

  skip_on_cran()

  prior_bias <- prior_weightfunction(
    "one-sided", c(0.025, 0.05), wf_fixed(c(1, 0.5, 0.1))
  )
  object <- bselmodel(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_bias = prior_bias, only_priors = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true(all(c(
    "sel_z_lower", "sel_z_upper", "sel_obs_bin", "sel_sign"
  ) %in% names(fit_data)))
  expect_false(any(grepl("sel_phack|phack_z|sel_segment|sel_kernel_mode", names(fit_data))))
  expect_match(syntax, "dselnorm_step", fixed = TRUE)
  expect_false(grepl("dselnorm_kernel|sel_phack|phack_z|sel_segment|sel_kernel_mode", syntax))

  omega_p <- matrix(c(1, 0.5, 0.1), nrow = 1)
  sei     <- object[["data"]][["outcome"]][["sei"]][1]
  yi_eval <- stats::qnorm(.04, lower.tail = FALSE) * sei

  selection_spec <- .selection_spec(
    priors           = object[["priors"]],
    yi               = yi_eval,
    sei              = sei,
    effect_direction = .effect_direction(object)
  )
  log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi_eval,
    mu_num         = matrix(0, nrow = 1, ncol = 1),
    sigma_num      = matrix(sei, nrow = 1, ncol = 1),
    sei            = sei,
    omega          = omega_p,
    selection_spec = selection_spec
  )
  prob <- stats::pnorm(
    selection_spec[["z_upper"]],
    mean = 0, sd = sei / sei
  ) - stats::pnorm(
    selection_spec[["z_lower"]],
    mean = 0, sd = sei / sei
  )
  ref_log_lik <- stats::dnorm(yi_eval, mean = 0, sd = sei, log = TRUE) +
    log(omega_p[selection_spec[["obs_bin"]]]) -
    log(sum(omega_p * prob))

  expect_equal(
    as.numeric(log_lik), as.numeric(ref_log_lik),
    tolerance = 1e-12
  )
})


test_that("bPEESE default respects supplied prior_unit_information_sd", {

  skip_on_cran()

  result <- bPEESE(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_unit_information_sd = 10,
    only_priors = TRUE
  )[["priors"]]

  expect_equal(
    result$outcome$bias$parameters$scale,
    RoBMA.get_option("default_bias_PEESE.scale") * sqrt(2) / 10
  )
})


# ============================================================================
# Tests for rescale_priors
# ============================================================================

test_that("rescale_priors scales prior distributions", {

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

  expect_error_cases(list(
    list(
      label  = "negative rescale_priors",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "SMD",
        rescale_priors = -1, only_priors = TRUE
      )),
      regexp = "rescale_priors"
    ),
    list(
      label  = "zero rescale_priors",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "SMD",
        rescale_priors = 0, only_priors = TRUE
      )),
      regexp = "rescale_priors"
    )
  ))
})

test_that("rescale_priors reports unsupported prior distributions once", {

  skip_on_cran()

  unsupported <- BayesTools::prior("gamma", parameters = list(shape = 1, rate = 1))

  message <- tryCatch(
    brma.norm(
      yi             = effect,
      sei            = std_err,
      data           = test_data,
      measure        = "SMD",
      prior_effect   = unsupported,
      rescale_priors = 2,
      only_priors    = TRUE
    ),
    error = conditionMessage
  )

  expect_match(message, "The 'gamma' prior distribution cannot be rescaled", fixed = TRUE)
  expect_match(message, "'normal', 'mnormal', 'cauchy', 'mcauchy', 't', 'mt', 'invgamma'", fixed = TRUE)
  expect_equal(length(gregexpr("prior distribution cannot be rescaled", message, fixed = TRUE)[[1]]), 1L)
})


# ============================================================================
# Tests for custom priors
# ============================================================================

test_that("Custom prior_effect and prior_heterogeneity rescale", {

  skip_on_cran()

  custom_effect <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5))
  custom_het    <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.3), truncation = list(0, Inf))

  result <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD",
    prior_effect = custom_effect, prior_heterogeneity = custom_het,
    only_priors = TRUE
  )[["priors"]]

  expect_equal(result$outcome$mu$parameters$sd, 0.5)
  expect_equal(result$outcome$tau$parameters$sd, 0.3)

  # With rescaling
  result_rescaled <- brma.norm(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD",
    prior_effect = custom_effect, prior_heterogeneity = custom_het, rescale_priors = 2, only_priors = TRUE
  )[["priors"]]
  expect_equal(result_rescaled$outcome$mu$parameters$sd,  0.5 * 2)
  expect_equal(result_rescaled$outcome$tau$parameters$sd, 0.3 * 2)
})


test_that("Different prior distribution types are accepted", {

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

test_that("Informed priors are accepted for medicine field with rescaling", {

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


test_that("Conflicting prior specifications are rejected", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "invalid informed-prior field",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "SMD",
        prior_informed_field = "invalid", only_priors = TRUE
      )),
      regexp = "prior_informed_field"
    ),
    list(
      label  = "informed-prior subfield without field",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "SMD",
        prior_informed_subfield = "neonatal", only_priors = TRUE
      )),
      regexp = "prior_informed_subfield.*prior_informed_field"
    ),
    list(
      label  = "UISD conflicts with informed priors",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "SMD",
        prior_unit_information_sd = 1.5, prior_informed_field = "medicine",
        only_priors = TRUE
      )),
      regexp = "prior_unit_information_sd|prior_informed_field|mutually exclusive"
    )
  ))
})


test_that("Informed scale priors are assigned", {

  skip_on_cran()

  result <- brma.norm(
    yi = effect, sei = std_err, scale = ~ scale_var,
    data = test_data, measure = "SMD",
    prior_informed_field = "medicine", only_priors = TRUE
  )[["priors"]]

  expect_equal(
    result$scale$scale_var$parameters$sd,
    RoBMA.get_option("default_informed_priors.scale")
  )
})


# ============================================================================
# Tests for heterogeneity allocation prior (cluster)
# ============================================================================

test_that("Heterogeneity allocation priors are assigned", {

  skip_on_cran()

  # Default Beta(1,1) when cluster specified
  result <- brma.norm(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]
  expect_true("rho" %in% names(result$outcome))
  expect_equal(result$outcome$rho$distribution, "beta")
  expect_equal(result$outcome$rho$parameters$alpha, 1)
  expect_equal(result$outcome$rho$parameters$beta, 1)

  # Custom prior
  custom_prior <- BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2))
  result_custom <- brma.norm(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    prior_heterogeneity_allocation = custom_prior, measure = "OR", only_priors = TRUE
  )[["priors"]]
  expect_equal(result_custom$outcome$rho$parameters$alpha, 2)

  # No rho when cluster not specified
  result_no_ids <- brma.norm(yi = effect, sei = std_err, data = test_data, measure = "SMD", only_priors = TRUE)
  expect_false("rho" %in% names(result_no_ids$outcome))
})


test_that("Constrained point priors outside support are rejected", {

  skip_on_cran()

  expect_error_cases(list(
    list(
      label  = "negative heterogeneity spike",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, data = test_data,
        measure = "SMD",
        prior_heterogeneity = BayesTools::prior("spike", parameters = list(location = -0.1)),
        only_priors = TRUE
      )),
      regexp = "prior_heterogeneity"
    ),
    list(
      label  = "heterogeneity-allocation spike above one",
      expr   = quote(brma.norm(
        yi = effect, sei = std_err, cluster = cluster, data = test_data,
        measure = "SMD",
        prior_heterogeneity_allocation = BayesTools::prior("spike", parameters = list(location = 1.1)),
        only_priors = TRUE
      )),
      regexp = "heterogeneity_allocation"
    )
  ))
})


# ============================================================================
# Tests for moderator and scale priors
# ============================================================================

test_that("Moderator priors are assigned", {

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

  # Single prior_mods applies to all terms
  result_mod_default <- brma.norm(
    yi = effect, sei = std_err, mods = ~ mod_cont,
    data = test_data, prior_mods = custom_mod[[1]],
    measure = "SMD", only_priors = TRUE
  )[["priors"]]
  expect_equal(result_mod_default$mods$mod_cont$parameters$sd, 0.3)
})


test_that("No-intercept moderator formulas use fixed-zero intercept", {

  object <- brma.norm(
    yi = effect, sei = std_err, mods = ~ 0 + mod_cont,
    data = test_data, measure = "SMD", only_priors = TRUE
  )
  result <- object[["priors"]]

  expect_null(result$outcome$mu)
  expect_true(BayesTools::is.prior.point(result$mods$intercept))
  expect_equal(mean(result$mods$intercept), 0)
  expect_false(BayesTools::is.prior.point(result$mods$mod_cont))

  formula_output <- BayesTools::JAGS_formula(
    parameter  = "mu",
    data       = object[["data"]][["mods"]],
    formula    = .create_fit_formula_list(data = object[["data"]], parameter = "mods"),
    prior_list = .create_fit_formula_prior_list(priors = object[["priors"]], parameter = "mods")
  )

  expect_true(BayesTools::is.prior.point(formula_output[["prior_list"]][["mu_intercept"]]))
  expect_equal(mean(formula_output[["prior_list"]][["mu_intercept"]]), 0)
  expect_match(formula_output[["formula_syntax"]], "mu_intercept", fixed = TRUE)
})

test_that("Scale formulas repair no-intercept terms through BayesTools public API", {

  env <- new.env(parent = emptyenv())

  data <- list(
    scale = test_data
  )
  formula_minus <- ~ scale_var - 1
  environment(formula_minus) <- env
  attr(data[["scale"]], "formula") <- formula_minus

  expect_warning(
    repaired <- .create_fit_formula_list(data = data, parameter = "scale"),
    "Intercept cannot be omitted"
  )
  expect_equal(attr(stats::terms(repaired), "intercept"), 1)
  expect_equal(attr(stats::terms(repaired), "term.labels"), "scale_var")
  expect_true(identical(environment(repaired), env))
  expect_true(isTRUE(attr(repaired, "log(intercept)")))

  formula_zero <- ~ 0 + scale_var
  environment(formula_zero) <- env
  attr(data[["scale"]], "formula") <- formula_zero
  expect_warning(
    repaired_zero <- .create_fit_formula_list(data = data, parameter = "scale"),
    "Intercept cannot be omitted"
  )
  expect_equal(attr(stats::terms(repaired_zero), "intercept"), 1)
  expect_equal(attr(stats::terms(repaired_zero), "term.labels"), "scale_var")

  formula_special <- ~ I(scale_var - 1) + offset(n - 1)
  environment(formula_special) <- env
  attr(data[["scale"]], "formula") <- formula_special
  expect_silent(
    repaired_special <- .create_fit_formula_list(data = data, parameter = "scale")
  )
  expect_equal(
    attr(stats::terms(repaired_special), "term.labels"),
    "I(scale_var - 1)"
  )
  expect_equal(attr(stats::terms(repaired_special), "offset"), 2)
})


test_that("Scale priors are assigned", {

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
    measure = "SMD",
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


test_that("Both mods and scale priors are assigned together", {

  skip_on_cran()

  object <- brma.norm(
    yi = effect, sei = std_err, mods = ~ mod_cont + mod_factor, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  )
  result <- object[["priors"]]

  expect_true(!is.null(result$mods))
  expect_true(!is.null(result$scale))
  expect_true("intercept" %in% names(result$mods))
  expect_true("intercept" %in% names(result$scale))
  expect_false("mu" %in% names(result$outcome))
  expect_false("tau" %in% names(result$outcome))

  formula_args <- .create_jags_formula_args(
    data   = object[["data"]],
    priors = object[["priors"]]
  )

  expect_named(formula_args[["formula_list"]],       c("mu", "log_tau"))
  expect_named(formula_args[["formula_data_list"]],  c("mu", "log_tau"))
  expect_named(formula_args[["formula_prior_list"]], c("mu", "log_tau"))
  expect_named(formula_args[["formula_scale_list"]], c("mu", "log_tau"))
  expect_length(formula_args[["formula_random_prior_list"]], 0L)
  expect_length(formula_args[["formula_random_effects_compile_list"]], 0L)
  expect_length(formula_args[["add_parameters"]], 0L)
  expect_true(isTRUE(attr(formula_args[["formula_list"]][["log_tau"]], "log(intercept)")))
})


# ============================================================================
# Tests for contrast
# ============================================================================

test_that("set_contrast_factor_predictors options are applied", {

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
      measure = "SMD",
      set_contrast_factor_predictors = "invalid", only_priors = TRUE
    ),
    regexp = "set_contrast_factor_predictors"
  )
})
