context("Prior input handling for RoBMA (mixture priors)")

skip_on_cran()

# Helper function to get is_null logical vector from mixture prior
# BayesTools stores this as "components" attribute with values "null" and "alternative"
.get_is_null <- function(prior) {
  components <- attr(prior, "components")
  if (is.null(components)) return(NULL)
  return(components == "null")
}

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
# Tests for mixture prior creation for effect and heterogeneity
# ============================================================================

test_that("Default priors create mixture with spike(0) null and normal alternative for effect", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # Check that mu is a mixture prior

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))

  # Check mixture components
  expect_equal(length(result$outcome$mu), 2)  # null + alternative
  expect_true(.get_is_null(result$outcome$mu)[1])
  expect_false(.get_is_null(result$outcome$mu)[2])

  # Check null is spike(0)
  expect_true(BayesTools::is.prior.point(result$outcome$mu[[1]]))
  expect_equal(mean(result$outcome$mu[[1]]), 0)

  # Check alternative scales with UISD
  expected_effect_sd <- sqrt(2) * RoBMA.get_option("default_UISD.effect")
  expect_equal(result$outcome$mu[[2]]$parameters$sd, expected_effect_sd)
})


test_that("Default priors create mixture with spike(0) null and truncated normal alternative for heterogeneity", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # Check that tau is a mixture prior
  expect_true(BayesTools::is.prior.mixture(result$outcome$tau))

  # Check mixture components
  expect_equal(length(result$outcome$tau), 2)  # null + alternative
  expect_true(.get_is_null(result$outcome$tau)[1])
  expect_false(.get_is_null(result$outcome$tau)[2])

  # Check null is spike(0)
  expect_true(BayesTools::is.prior.point(result$outcome$tau[[1]]))
  expect_equal(mean(result$outcome$tau[[1]]), 0)

  # Check alternative is truncated at 0
  expect_equal(result$outcome$tau[[2]]$truncation$lower, 0)

  # Check alternative scales with UISD
  expected_het_sd <- sqrt(2) * RoBMA.get_option("default_UISD.heterogeneity")
  expect_equal(result$outcome$tau[[2]]$parameters$sd, expected_het_sd)
})


test_that("Mixture priors scale correctly with different measures", {

  measures_uisd <- list(SMD = sqrt(2), ZCOR = 1, OR = 2, RR = 2, HR = 2)

  for (measure in names(measures_uisd)) {
    result <- RoBMA(
      yi = effect, sei = std_err, data = test_data,
      measure = measure, only_priors = TRUE
    )[["priors"]]

    expected_effect_sd <- measures_uisd[[measure]] * RoBMA.get_option("default_UISD.effect")
    expected_het_sd    <- measures_uisd[[measure]] * RoBMA.get_option("default_UISD.heterogeneity")

    # Check alternative component (second in mixture)
    expect_equal(result$outcome$mu[[2]]$parameters$sd, expected_effect_sd)
    expect_equal(result$outcome$tau[[2]]$parameters$sd, expected_het_sd)
  }
})


test_that("rescale_priors applies to alternative components in mixture", {

  result_base <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  for (scale_factor in c(0.5, 2)) {
    result_rescaled <- RoBMA(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD", rescale_priors = scale_factor, only_priors = TRUE
    )[["priors"]]

    # Alternative component (index 2) should be rescaled
    expect_equal(
      result_rescaled$outcome$mu[[2]]$parameters$sd,
      result_base$outcome$mu[[2]]$parameters$sd * scale_factor
    )
    expect_equal(
      result_rescaled$outcome$tau[[2]]$parameters$sd,
      result_base$outcome$tau[[2]]$parameters$sd * scale_factor
    )
  }
})


# ============================================================================
# Tests for custom null/alternative prior lists
# ============================================================================

test_that("Custom prior_effect replaces alternative component", {

  custom_effect <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5))

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect = custom_effect, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_equal(result$outcome$mu[[2]]$parameters$sd, 0.5)
})


test_that("Custom prior_effect_null replaces null component", {

  custom_null <- BayesTools::prior("spike", parameters = list(location = 0.1))

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect_null = custom_null, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_true(BayesTools::is.prior.point(result$outcome$mu[[1]]))
  expect_equal(mean(result$outcome$mu[[1]]), 0.1)
})


test_that("List of priors creates mixture with multiple components", {

  alt_priors <- list(
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5)),
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1.0))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect = alt_priors, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_equal(length(result$outcome$mu), 3)  # 1 null + 2 alternatives

  # Check is_null attribute
  is_null <- .get_is_null(result$outcome$mu)
  expect_equal(sum(is_null), 1)   # 1 null
  expect_equal(sum(!is_null), 2)  # 2 alternatives
})


test_that("NULL/FALSE omits null hypothesis component", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect_null = NULL, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_equal(length(result$outcome$mu), 1)  # only alternative

  # All components should be alternative (not null)
  is_null <- .get_is_null(result$outcome$mu)
  expect_equal(sum(is_null), 0)
})


test_that("NULL/FALSE omits alternative hypothesis component", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect = NULL, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_equal(length(result$outcome$mu), 1)  # only null

  # All components should be null
  is_null <- .get_is_null(result$outcome$mu)
  expect_equal(sum(is_null), 1)
})


test_that("FALSE works same as NULL for omitting components", {

  result_null <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect_null = NULL, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  result_false <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect_null = FALSE, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_equal(length(result_null$outcome$mu), length(result_false$outcome$mu))
  expect_equal(.get_is_null(result_null$outcome$mu), .get_is_null(result_false$outcome$mu))
})


# ============================================================================
# Tests for publication bias model types
# ============================================================================

test_that("model_type '2w' creates correct bias mixture", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", model_type = "2w", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$bias))

  # Check for prior_none as null
  is_null <- .get_is_null(result$outcome$bias)
  expect_true(any(is_null))

  # Check that weight functions are present (not PET/PEESE)
  has_wf <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.weightfunction(p))
  expect_true(any(has_wf))

  # No PET/PEESE in 2w
  has_pet <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PET(p))
  has_peese <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PEESE(p))
  expect_false(any(has_pet))
  expect_false(any(has_peese))
})


test_that("model_type '6w' creates correct bias mixture", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", model_type = "6w", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$bias))

  # Count weight functions (should be 6)
  has_wf <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.weightfunction(p))
  expect_equal(sum(has_wf), 6)

  # No PET/PEESE
  has_pet <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PET(p))
  has_peese <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PEESE(p))
  expect_false(any(has_pet))
  expect_false(any(has_peese))
})


test_that("model_type 'PP' creates correct bias mixture", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", model_type = "PP", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$bias))

  # Should have PET and PEESE
  has_pet <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PET(p))
  has_peese <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PEESE(p))
  expect_true(any(has_pet))
  expect_true(any(has_peese))

  # No weight functions
  has_wf <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.weightfunction(p))
  expect_false(any(has_wf))
})


test_that("model_type 'PSMA' creates correct bias mixture", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", model_type = "PSMA", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$bias))

  # Should have all: weight functions, PET, and PEESE
  has_wf <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.weightfunction(p))
  has_pet <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PET(p))
  has_peese <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.PEESE(p))

  expect_true(any(has_wf))
  expect_true(any(has_pet))
  expect_true(any(has_peese))
  expect_equal(sum(has_wf), 6)  # 6 weight functions in PSMA
})


test_that("Custom prior_bias and prior_bias_null override model_type", {

  custom_bias <- BayesTools::prior_weightfunction(
    distribution = "two.sided",
    parameters = list(alpha = c(1, 1), steps = c(0.05))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_bias = custom_bias, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$bias))

  # Should have the custom weight function as alternative
  has_wf <- sapply(result$outcome$bias, function(p) BayesTools::is.prior.weightfunction(p))
  expect_equal(sum(has_wf), 1)
})


test_that("prior_bias_null = NULL omits null bias hypothesis", {

  custom_bias <- BayesTools::prior_weightfunction(
    distribution = "two.sided",
    parameters = list(alpha = c(1, 1), steps = c(0.05))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_bias = custom_bias, prior_bias_null = NULL,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  is_null <- .get_is_null(result$outcome$bias)
  expect_equal(sum(is_null), 0)  # no null component
})


test_that("Bias presets only require UISD when PEESE is present", {

  custom_effect <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5))
  custom_het    <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.3), truncation = list(0, Inf))

  expect_no_error(
    RoBMA(
      yi = effect, sei = std_err, data = test_data,
      measure = "GEN", model_type = "6w",
      prior_effect = custom_effect, prior_heterogeneity = custom_het,
      only_priors = TRUE
    )
  )
})


test_that("Partial bias customization respects model_type alternatives", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", model_type = "PP",
    prior_bias_null = NULL, only_priors = TRUE
  )[["priors"]]

  has_wf    <- sapply(result$outcome$bias, BayesTools::is.prior.weightfunction)
  has_pet   <- sapply(result$outcome$bias, BayesTools::is.prior.PET)
  has_peese <- sapply(result$outcome$bias, BayesTools::is.prior.PEESE)

  expect_equal(sum(has_wf), 0)
  expect_equal(sum(has_pet), 1)
  expect_equal(sum(has_peese), 1)
})


test_that("Empty top-level mixtures throw clear errors", {

  expect_error(
    RoBMA(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD",
      prior_effect = NULL, prior_effect_null = NULL,
      only_priors = TRUE
    ),
    regexp = "effect"
  )

  expect_error(
    RoBMA(
      yi = effect, sei = std_err, data = test_data,
      measure = "SMD",
      prior_bias = NULL, prior_bias_null = NULL,
      only_priors = TRUE
    ),
    regexp = "bias"
  )
})


# ============================================================================
# Tests for heterogeneity allocation mixture (multilevel models)
# ============================================================================

test_that("Heterogeneity allocation creates mixture with cluster", {

  result <- RoBMA(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true("rho" %in% names(result$outcome))
  expect_true(BayesTools::is.prior.mixture(result$outcome$rho))

  # Default: no null hypothesis for rho (only alternative Beta(1,1))
  is_null <- .get_is_null(result$outcome$rho)
  expect_equal(sum(is_null), 0)

  # Check alternative is Beta(1,1)
  expect_equal(result$outcome$rho[[1]]$distribution, "beta")
  expect_equal(result$outcome$rho[[1]]$parameters$alpha, 1)
  expect_equal(result$outcome$rho[[1]]$parameters$beta, 1)
})


test_that("Custom prior_heterogeneity_allocation works", {

  custom_prior <- BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2))

  result <- RoBMA(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    prior_heterogeneity_allocation = custom_prior,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_equal(result$outcome$rho[[1]]$parameters$alpha, 2)
  expect_equal(result$outcome$rho[[1]]$parameters$beta, 2)
})


test_that("prior_heterogeneity_allocation_null adds null hypothesis", {

  null_prior <- BayesTools::prior("spike", parameters = list(location = 0.5))

  result <- RoBMA(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    prior_heterogeneity_allocation_null = null_prior,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  is_null <- .get_is_null(result$outcome$rho)
  expect_equal(sum(is_null), 1)  # now has null component
  expect_true(BayesTools::is.prior.point(result$outcome$rho[[1]]))
})


# ============================================================================
# Tests for moderator and scale priors with mixtures
# ============================================================================

test_that("Moderator priors create mixture with intercept from effect prior", {

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(!is.null(result$mods))
  expect_true("intercept" %in% names(result$mods))
  expect_null(result$outcome$mu)  # mu forwarded to intercept

  # Intercept should be a mixture prior
  expect_true(BayesTools::is.prior.mixture(result$mods$intercept))
})


test_that("Moderator term priors include null and alternative", {

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # mod_cont should have a mixture prior
  expect_true(BayesTools::is.prior.mixture(result$mods$mod_cont))

  is_null <- .get_is_null(result$mods$mod_cont)
  expect_true(any(is_null))   # has null component
  expect_true(any(!is_null))  # has alternative component
})


test_that("Custom prior_mods replaces alternative for terms", {

  custom_mod <- list(mod_cont = BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.3)))

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont,
    prior_mods = custom_mod, data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # Alternative component should have custom prior
  is_null <- .get_is_null(result$mods$mod_cont)
  alt_idx <- which(!is_null)[1]
  expect_equal(result$mods$mod_cont[[alt_idx]]$parameters$sd, 0.3)
})


test_that("Custom prior_mods_null replaces null for terms", {

  custom_null <- list(mod_cont = BayesTools::prior("spike", parameters = list(location = 0.1)))

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont,
    prior_mods_null = custom_null, data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # Null component should have custom prior
  is_null <- .get_is_null(result$mods$mod_cont)
  null_idx <- which(is_null)[1]
  expect_true(BayesTools::is.prior.point(result$mods$mod_cont[[null_idx]]))
  expect_equal(mean(result$mods$mod_cont[[null_idx]]), 0.1)
})


# ============================================================================
# Tests for informed priors with mixtures
# ============================================================================

test_that("Informed priors work with mixture structure", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "OR", prior_informed_field = "medicine", only_priors = TRUE
  )[["priors"]]

  # Check mixture structure preserved
  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(result$outcome$tau))

  # Alternative component should be informed prior - compare key properties
  expected_mu  <- BayesTools::prior_informed(name = "cochrane", parameter = "effect", type = "logOR")
  expected_tau <- BayesTools::prior_informed(name = "cochrane", parameter = "heterogeneity", type = "logOR")

  expect_equal(result$outcome$mu[[2]]$distribution, expected_mu$distribution)
  expect_equal(result$outcome$mu[[2]]$parameters, expected_mu$parameters)
  expect_equal(result$outcome$tau[[2]]$distribution, expected_tau$distribution)
  expect_equal(result$outcome$tau[[2]]$parameters, expected_tau$parameters)
})


test_that("Informed priors with rescaling apply to alternative", {

  result_base <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_informed_field = "medicine", only_priors = TRUE
  )[["priors"]]

  result_rescaled <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    measure = "SMD", prior_informed_field = "medicine",
    rescale_priors = 2, only_priors = TRUE
  )[["priors"]]

  scale_param <- if ("sd" %in% names(result_base$outcome$mu[[2]]$parameters)) "sd" else "scale"
  expect_equal(
    result_rescaled$outcome$mu[[2]]$parameters[[scale_param]],
    result_base$outcome$mu[[2]]$parameters[[scale_param]] * 2
  )
})


# ============================================================================
# Tests for error conditions
# ============================================================================

test_that("Invalid model_type throws error", {

  expect_error(
    RoBMA(yi = effect, sei = std_err, data = test_data, model_type = "invalid", only_priors = TRUE),
    regexp = "model_type"
  )
})


test_that("Conflicting prior specifications throw errors", {

  # Both UISD and informed priors
  expect_error(
    RoBMA(
      yi = effect, sei = std_err, data = test_data,
      prior_unit_information_sd = 1.5, prior_informed_field = "medicine", only_priors = TRUE
    ),
    regexp = "prior_unit_information_sd|prior_informed_field|mutually exclusive"
  )
})


test_that("Invalid prior types throw appropriate errors", {

  # Non-prior object
  expect_error(
    RoBMA(yi = effect, sei = std_err, data = test_data, prior_effect = "not a prior", only_priors = TRUE),
    regexp = "prior_effect"
  )
})


test_that("GEN measure requires ni or prior_unit_information_sd for RoBMA", {

  expect_error(
    RoBMA(yi = effect, sei = std_err, data = test_data, measure = "GEN", only_priors = TRUE),
    regexp = "ni|unit_information_sd|UISD|Sample size"
  )

  # Should work with ni
  expect_no_error(
    RoBMA(yi = effect, sei = std_err, ni = n, data = test_data, measure = "GEN", only_priors = TRUE)
  )

  # Should work with prior_unit_information_sd
  expect_no_error(
    RoBMA(yi = effect, sei = std_err, data = test_data, measure = "GEN", prior_unit_information_sd = 1.5, only_priors = TRUE)
  )
})


# ============================================================================
# Tests for scale priors (simple priors, not mixtures for terms)
# ============================================================================

test_that("Scale priors create mixture with intercept from heterogeneity prior", {

  result <- suppressWarnings(RoBMA(
    yi = effect, sei = std_err, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]])

  expect_true(!is.null(result$scale))
  expect_true("intercept" %in% names(result$scale))
  expect_null(result$outcome$tau)  # tau forwarded to intercept

  # Intercept should be a mixture prior
  expect_true(BayesTools::is.prior.mixture(result$scale$intercept))

  # Note: spike(0) null is removed for scale since log(0) is undefined
  # Check that remaining components are truncated (positive heterogeneity)
  for (i in seq_along(result$scale$intercept)) {
    if (!BayesTools::is.prior.point(result$scale$intercept[[i]])) {
      expect_equal(result$scale$intercept[[i]]$truncation$lower, 0)
    }
  }
})


test_that("Scale term priors are mixture priors (same as mods)", {

  # Scale terms use .assign_prior_list.terms_mixture() like mods
  # returning mixture priors for model-averaging

  result <- suppressWarnings(RoBMA(
    yi = effect, sei = std_err, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]])

  # scale_var should be a mixture prior (like mods terms)
  expect_true(BayesTools::is.prior.mixture(result$scale$scale_var))
  # Check it has both null and alternative components
  is_null <- .get_is_null(result$scale$scale_var)
  expect_true(any(is_null))   # has null component(s)
  expect_true(any(!is_null))  # has alternative component(s)
})


test_that("Custom prior_scale replaces alternative for scale terms", {

  custom_scale <- list(scale_var = BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.4)))

  result <- suppressWarnings( RoBMA(
    yi = effect, sei = std_err, scale = ~ scale_var,
    prior_scale = custom_scale, data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]])

  # Result should be a mixture
  expect_true(BayesTools::is.prior.mixture(result$scale$scale_var))
  # Custom prior should be in the alternative component
  is_null <- .get_is_null(result$scale$scale_var)
  alt_indices <- which(!is_null)
  expect_equal(result$scale$scale_var[[alt_indices[1]]]$parameters$sd, 0.4)
})


# ============================================================================
# Tests for custom heterogeneity prior lists
# ============================================================================

test_that("Custom prior_heterogeneity replaces alternative component", {

  custom_het <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.3), truncation = list(0, Inf))

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_heterogeneity = custom_het, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$tau))
  expect_equal(result$outcome$tau[[2]]$parameters$sd, 0.3)
})


test_that("Custom prior_heterogeneity_null replaces null component", {

  custom_null <- BayesTools::prior("spike", parameters = list(location = 0.1))

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_heterogeneity_null = custom_null, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$tau))
  expect_true(BayesTools::is.prior.point(result$outcome$tau[[1]]))
  expect_equal(mean(result$outcome$tau[[1]]), 0.1)
})


test_that("List of heterogeneity priors creates mixture with multiple alternatives", {

  alt_priors <- list(
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.2), truncation = list(0, Inf)),
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.5), truncation = list(0, Inf))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_heterogeneity = alt_priors, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$tau))
  expect_equal(length(result$outcome$tau), 3)  # 1 null + 2 alternatives

  is_null <- .get_is_null(result$outcome$tau)
  expect_equal(sum(is_null), 1)
  expect_equal(sum(!is_null), 2)
})


test_that("prior_heterogeneity = NULL creates null-only mixture", {

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_heterogeneity = NULL, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$tau))
  expect_equal(length(result$outcome$tau), 1)

  is_null <- .get_is_null(result$outcome$tau)
  expect_equal(sum(is_null), 1)
})


# ============================================================================
# Tests for factor moderator priors
# ============================================================================

test_that("Factor moderator creates mixture prior", {

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_factor,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # Factor moderator is stored as mod_factor (not individual levels)
  expect_true("mod_factor" %in% names(result$mods))
  expect_true(BayesTools::is.prior.mixture(result$mods$mod_factor))

  is_null <- .get_is_null(result$mods$mod_factor)
  expect_true(any(is_null))
  expect_true(any(!is_null))
})


test_that("Multiple moderators each get mixture priors", {

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont + mod_factor,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # Both moderators should have mixture priors
  expect_true(BayesTools::is.prior.mixture(result$mods$mod_cont))
  expect_true(BayesTools::is.prior.mixture(result$mods$mod_factor))
})


# ============================================================================
# Tests for mixed moderator specifications (some with null, some without)
# ============================================================================

test_that("Different moderators can have different null specifications", {

  # mod_cont: omit null (alternative only)
  # mod_factor: keep default (null + alternative)
  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont + mod_factor,
    prior_mods_null = list(mod_cont = NULL),  # omit null for mod_cont only
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # mod_cont should have no null component
  is_null_cont <- .get_is_null(result$mods$mod_cont)
  expect_equal(sum(is_null_cont), 0)

  # mod_factor should still have null component
  is_null_factor <- .get_is_null(result$mods$mod_factor)
  expect_equal(sum(is_null_factor), 1)
})


test_that("Different moderators can have different alternative specifications", {

  # Custom alternative for mod_cont, default for mod_factor
  custom_mod <- list(mod_cont = BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.2)))

  result <- RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont + mod_factor,
    prior_mods = custom_mod, data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  # mod_cont should have custom alternative
  is_null_cont <- .get_is_null(result$mods$mod_cont)
  alt_idx <- which(!is_null_cont)[1]
  expect_equal(result$mods$mod_cont[[alt_idx]]$parameters$sd, 0.2)

  # mod_factor should have default alternative
  is_null_factor <- .get_is_null(result$mods$mod_factor)
  expect_true(any(!is_null_factor))  # has alternative component
})


# ============================================================================
# Tests for list of bias priors
# ============================================================================

test_that("List of bias priors creates mixture with multiple alternatives", {

  bias_priors <- list(
    BayesTools::prior_weightfunction("two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05))),
    BayesTools::prior_PET("normal", parameters = list(mean = 0, sd = 1))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_bias = bias_priors, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$bias))

  # Should have 1 null (prior_none) + 2 custom alternatives
  expect_equal(length(result$outcome$bias), 3)

  has_wf  <- sapply(result$outcome$bias, BayesTools::is.prior.weightfunction)
  has_pet <- sapply(result$outcome$bias, BayesTools::is.prior.PET)
  expect_equal(sum(has_wf), 1)
  expect_equal(sum(has_pet), 1)
})


# ============================================================================
# Tests for heterogeneity allocation list priors
# ============================================================================

test_that("List of heterogeneity allocation priors creates mixture", {

  rho_priors <- list(
    BayesTools::prior("beta", parameters = list(alpha = 1, beta = 1)),
    BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    prior_heterogeneity_allocation = rho_priors,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$rho))
  expect_equal(length(result$outcome$rho), 2)  # 2 alternatives (no null by default)

  # Both should be beta distributions with different parameters
  expect_equal(result$outcome$rho[[1]]$parameters$alpha, 1)
  expect_equal(result$outcome$rho[[2]]$parameters$alpha, 2)
})


test_that("Heterogeneity allocation with both null and alternative lists", {

  alt_priors  <- list(BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2)))
  null_priors <- list(BayesTools::prior("spike", parameters = list(location = 0.5)))

  result <- RoBMA(
    yi = effect, sei = std_err, cluster = cluster, data = test_data,
    prior_heterogeneity_allocation = alt_priors,
    prior_heterogeneity_allocation_null = null_priors,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$rho))
  expect_equal(length(result$outcome$rho), 2)  # 1 null + 1 alternative

  is_null <- .get_is_null(result$outcome$rho)
  expect_equal(sum(is_null), 1)
  expect_equal(sum(!is_null), 1)
})


# ============================================================================
# Tests for multiple null priors (list of nulls)
# ============================================================================

test_that("List of null priors creates mixture with multiple null components", {

  null_priors <- list(
    BayesTools::prior("spike", parameters = list(location = 0)),
    BayesTools::prior("spike", parameters = list(location = 0.1))
  )

  result <- RoBMA(
    yi = effect, sei = std_err, data = test_data,
    prior_effect_null = null_priors, measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_true(BayesTools::is.prior.mixture(result$outcome$mu))
  expect_equal(length(result$outcome$mu), 3)  # 2 nulls + 1 default alternative

  is_null <- .get_is_null(result$outcome$mu)
  expect_equal(sum(is_null), 2)
  expect_equal(sum(!is_null), 1)
})


# ============================================================================
# Tests for interaction between effect/heterogeneity and mods/scale
# ============================================================================

test_that("Effect prior mixture is correctly forwarded to mods intercept", {

  custom_effect <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 0.8))
  custom_null   <- BayesTools::prior("spike", parameters = list(location = 0.1))

  result <- suppressWarnings(RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont,
    prior_effect = custom_effect, prior_effect_null = custom_null,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]])

  # Intercept should have the custom priors
  expect_true(BayesTools::is.prior.mixture(result$mods$intercept))

  is_null  <- .get_is_null(result$mods$intercept)
  null_idx <- which(is_null)[1]
  alt_idx  <- which(!is_null)[1]

  expect_equal(mean(result$mods$intercept[[null_idx]]), 0.1)
  expect_equal(result$mods$intercept[[alt_idx]]$parameters$sd, 0.8)
})


test_that("Both mods and scale work together with mixtures", {

  result <- suppressWarnings(RoBMA(
    yi = effect, sei = std_err, mods = ~ mod_cont, scale = ~ scale_var,
    data = test_data, measure = "SMD", only_priors = TRUE
  )[["priors"]])

  expect_true(!is.null(result$mods))
  expect_true(!is.null(result$scale))
  expect_true("intercept" %in% names(result$mods))
  expect_true("intercept" %in% names(result$scale))
  expect_null(result$outcome$mu)
  expect_null(result$outcome$tau)

  # mods intercept should be mixture
  expect_true(BayesTools::is.prior.mixture(result$mods$intercept))
})
