context("Selection kernel")
skip_on_cran()

test_that("selection model fit data and syntax use only the selected-normal kernel", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- bselmodel(
    yi                        = c(.1, .2, .3),
    sei                       = c(.1, .1, .1),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true(all(c(
    "sel_z_lower", "sel_z_upper", "sel_obs_bin", "sel_sign"
  ) %in% names(fit_data)))
  expect_false(any(grepl("sel_phack|phack_z|sel_segment|sel_kernel_mode", names(fit_data))))
  expect_match(syntax, "dselnorm_step", fixed = TRUE)
  expect_match(syntax, "sel_total_sd\\[i\\]")
  expect_false(grepl("dselnorm_kernel|sel_phack|phack_z|sel_segment|sel_kernel_mode", syntax))
})

test_that("mixed normal-step bias syntax uses scalar step switch", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- RoBMA(
    yi                        = c(.1, .2, .3),
    sei                       = c(.1, .1, .1),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_bias_null           = BayesTools::prior_none(),
    prior_effect              = BayesTools::prior("normal", parameters = list(0, 1)),
    prior_effect_null         = NULL,
    prior_heterogeneity       = BayesTools::prior("invgamma", parameters = list(1, .15)),
    prior_heterogeneity_null  = NULL,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  syntax   <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_true(all(c(
    "sel_z_lower", "sel_z_upper", "sel_obs_bin", "sel_sign"
  ) %in% names(fit_data)))
  expect_false(any(grepl("sel_phack|phack_z|sel_segment", names(fit_data))))
  expect_match(syntax, "sel_kernel_mode_active", fixed = TRUE)
  expect_match(syntax, "dselnorm_step_switch", fixed = TRUE)
  expect_false(grepl("dselnorm_kernel", syntax, fixed = TRUE))
})

test_that("selection spec sets probability telescoping flag once", {

  yi  <- c(.1, .2, .3)
  sei <- c(.1, .1, .1)
  safe_prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = BayesTools::wf_fixed(c(1, .5))
  )
  safe_spec <- .selection_spec(
    priors           = list(outcome = list(bias = safe_prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive",
    signed_data      = TRUE
  )

  expect_true(safe_spec[["telescope_probabilities"]])
  expect_identical(safe_spec[["jags_data"]][["sel_telescope_probabilities"]], 1L)

  wide_prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = BayesTools::wf_fixed(c(1, 200))
  )
  expect_warning(
    wide_spec <- .selection_spec(
      priors           = list(outcome = list(bias = wide_prior)),
      yi               = yi,
      sei              = sei,
      effect_direction = "positive",
      signed_data      = TRUE
    ),
    "probability telescoping disabled"
  )

  expect_false(wide_spec[["telescope_probabilities"]])
  expect_identical(wide_spec[["jags_data"]][["sel_telescope_probabilities"]], 0L)
})

test_that("branch kernel modes follow BayesTools phack and combined labels", {

  selection <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025),
    weights = BayesTools::wf_fixed(c(1, .5))
  )
  phacking <- BayesTools::prior_phacking(form = "linear")
  bias <- BayesTools::prior_mixture(list(
    BayesTools::prior_none(),
    selection,
    phacking,
    BayesTools::prior_bias(selection, phacking)
  ))
  backend <- BayesTools::selection_backend_spec(bias)

  expect_equal(
    backend[["branch_type"]],
    c("none", "weightfunction", "phack", "combined")
  )
  expect_equal(
    .selection_branch_kernel_modes(backend),
    c(
      SELKERNEL_NORMAL,
      SELKERNEL_STEP,
      SELKERNEL_PHACK_POWER,
      SELKERNEL_STEP_PHACK_POWER
    )
  )
})

test_that("single two-sided bselmodel uses active full-grid omega in JAGS", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "two-sided",
    steps   = c(.05, .10),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- bselmodel(
    yi                        = c(.24, .31, -.18, .05),
    sei                       = c(.10, .12, .09, .20),
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  priors   <- .create_fit_priors(object[["data"]], object[["priors"]])
  syntax   <- BayesTools::JAGS_add_priors(
    .create_model_syntax(object[["data"]], object[["priors"]]),
    priors
  )

  expect_equal(length(fit_data[["sel_z_lower"]]), 5L)
  expect_equal(length(fit_data[["sel_z_upper"]]), 5L)
  expect_match(syntax, "omega_local\\[2\\] <- 0.5")
  expect_match(syntax, "omega\\[4\\] <- omega_local\\[2\\]")
  expect_match(syntax, "omega\\[5\\] <- omega_local\\[1\\]")
  expect_false(grepl("omega\\[6\\]", syntax))
  expect_true("omega" %in% BayesTools::JAGS_to_monitor(priors))

})

test_that("JAGS permits fixed zero weights for empty p-value bins", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, 0, .5))
  )
  sei <- c(.10, .10)
  yi  <- stats::qnorm(c(.01, .20), lower.tail = FALSE) * sei

  object <- bselmodel(
    yi                        = yi,
    sei                       = sei,
    measure                   = "SMD",
    prior_bias                = prior_bias,
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  priors   <- .create_fit_priors(object[["data"]], object[["priors"]])
  syntax   <- BayesTools::JAGS_add_priors(
    .create_model_syntax(object[["data"]], object[["priors"]]),
    priors
  )

  expect_equal(fit_data[["sel_obs_bin"]], c(1L, 3L))
  expect_match(syntax, "omega[2] <- 0", fixed = TRUE)

  expect_false(grepl("omega[2] ~", syntax, fixed = TRUE))
})

test_that("selection p-value bins are upper-closed at step boundaries", {

  p_cuts <- c(0, .025, .50, 1)
  sei    <- rep(1, 5)
  p_val  <- c(.025, .0250001, .50, .5000001, .0249999)
  yi     <- stats::qnorm(p_val, lower.tail = FALSE) * sei

  expect_equal(
    .selection_obs_bin(yi, sei, p_cuts, sign = 1L),
    c(1L, 2L, 2L, 3L, 1L)
  )
  expect_equal(.selection_step_bin_from_z(0, p_cuts), 2L)
  expect_equal(.selection_step_bin_from_z(stats::qnorm(.025, lower.tail = FALSE), p_cuts), 1L)
})

test_that("selection p-value cut breaks include endpoints and are sorted", {

  expect_error(
    .selection_p_bin(c(.01, .20), c(.025, .05, 1)),
    regexp = "endpoints 0 and 1"
  )
  expect_error(
    .selection_p_bin(c(.01, .20), c(0, .05, .025, 1)),
    regexp = "strictly increasing"
  )
})

test_that("selection omega extraction orders indexed posterior columns numerically", {

  posterior_samples <- matrix(
    seq_len(30),
    nrow = 3,
    dimnames = list(NULL, c("omega[10]", "mu", "omega[2]", "omega[1]", "omega[3]",
                            "omega[4]", "omega[5]", "omega[6]", "omega[7]", "omega[8]"))
  )
  posterior_samples <- cbind(posterior_samples, "omega[9]" = 31:33)

  selection_spec <- list(jags_omega = "omega", n_bins = 10L)
  omega          <- .extract_selection_omega_samples(posterior_samples, selection_spec)

  expect_equal(colnames(omega), paste0("omega[", 1:10, "]"))
  expect_equal(omega[, 1], posterior_samples[, "omega[1]"])
  expect_equal(omega[, 10], posterior_samples[, "omega[10]"])

  custom_samples            <- posterior_samples
  colnames(custom_samples) <- sub("^omega", "custom_omega", colnames(custom_samples))
  selection_spec            <- list(jags_omega = "custom_omega", n_bins = 10L)
  custom_omega              <- .extract_selection_omega_samples(custom_samples, selection_spec)

  expect_equal(colnames(custom_omega), paste0("custom_omega[", 1:10, "]"))
  expect_equal(custom_omega[, 2], custom_samples[, "custom_omega[2]"])
})
