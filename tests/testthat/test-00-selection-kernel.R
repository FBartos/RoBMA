context("Selection kernel")
skip_on_cran()

.test_step_spec <- function(yi, sei, effect_direction = "positive") {

  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05, .50),
    weights = BayesTools::wf_fixed(c(1, .7, .35, .2))
  )

  .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = effect_direction,
    signed_data      = FALSE
  )
}

.test_two_sided_step_spec <- function(yi, sei, effect_direction = "positive") {

  prior <- BayesTools::prior_weightfunction(
    side    = "two-sided",
    steps   = c(.05, .10),
    weights = BayesTools::wf_fixed(c(1, .7, .35))
  )

  .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = effect_direction,
    signed_data      = FALSE
  )
}

test_that("selection model fit data and syntax use only the selected-normal kernel", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )

  object <- bselmodel(
    yi                        = c(.1, .2, .3),
    sei                       = c(.1, .1, .1),
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

  skip_if_not_installed("rjags")

  fit <- suppressWarnings(BayesTools::JAGS_fit(
    model_syntax      = .create_model_syntax(object[["data"]], object[["priors"]]),
    data              = fit_data,
    prior_list        = priors,
    chains            = 1,
    adapt             = 50,
    burnin            = 50,
    sample            = 100,
    required_packages = "RoBMA",
    silent            = TRUE,
    seed              = 91
  ))

  expect_false(inherits(fit, "error"))
  posterior <- as.matrix(coda::as.mcmc(fit))
  expect_true(all(paste0("omega[", 1:5, "]") %in% colnames(posterior)))
  expect_equal(unname(posterior[, "omega[2]"]), rep(.5, nrow(posterior)))
  expect_equal(unname(posterior[, "omega[4]"]), rep(.5, nrow(posterior)))
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

  skip_if_not_installed("rjags")

  fit <- suppressWarnings(BayesTools::JAGS_fit(
    model_syntax      = .create_model_syntax(object[["data"]], object[["priors"]]),
    data              = fit_data,
    prior_list        = priors,
    chains            = 1,
    adapt             = 50,
    burnin            = 50,
    sample            = 100,
    required_packages = "RoBMA",
    silent            = TRUE,
    seed              = 92
  ))

  posterior <- as.matrix(coda::as.mcmc(fit))
  expect_equal(unname(posterior[, "omega[2]"]), rep(0, nrow(posterior)))
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

test_that("step selected-normal kernel matches an independent p-bin reference", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- c(-.15, 0, .196, .35, .9)
  sei   <- c(.10, .12, .10, .25, .50)
  mu    <- matrix(c(
    -.05, .02, .10, .20, .25,
     .01, .04, .12, .18, .30,
     .08, .06, .14, .24, .35
  ), nrow = 3, byrow = TRUE)
  sigma <- matrix(c(
    .11, .18, .20, .30, .55,
    .15, .22, .25, .35, .60,
    .18, .25, .30, .40, .70
  ), nrow = 3, byrow = TRUE)
  omega <- matrix(c(
    1, .70, .35, .20,
    1, .20, .80, .05,
    1, .95, .40, .10
  ), nrow = 3, byrow = TRUE)
  spec <- .test_step_spec(yi, sei)

  out <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )

  expect_equal(out, .test_step_reference(yi, mu, sigma, sei, omega, spec),
               tolerance = 1e-12)
})

test_that("unit weights reduce one- and two-sided step kernels to normal", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- c(-.20, .00, .15, .55)
  sei   <- c(.08, .10, .17, .26)
  mu    <- matrix(c(
    -.05, .02, .10, .30,
     .04, .01, .12, .22
  ), nrow = 2, byrow = TRUE)
  sigma <- matrix(c(
    .12, .18, .28, .34,
    .16, .20, .25, .40
  ), nrow = 2, byrow = TRUE)

  for (spec in list(
    .test_step_spec(yi, sei),
    .test_two_sided_step_spec(yi, sei)
  )) {
    omega <- matrix(1, nrow = nrow(mu), ncol = spec[["n_bins"]])

    log_pdf <- .selnorm_kernel_loglik_matrix(
      yi             = yi,
      mu_num         = mu,
      sigma_num      = sigma,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )
    cdf <- .selnorm_kernel_cdf_matrix(
      q              = yi,
      mean           = mu,
      sd             = sigma,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )
    moments <- .selnorm_kernel_moments_matrix(
      mean           = mu,
      sd             = sigma,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )

    expect_equal(log_pdf, stats::dnorm(
      matrix(yi, nrow = nrow(mu), ncol = length(yi), byrow = TRUE),
      mean = mu,
      sd   = sigma,
      log  = TRUE
    ), tolerance = 1e-12)
    expect_equal(cdf, stats::pnorm(
      matrix(yi, nrow = nrow(mu), ncol = length(yi), byrow = TRUE),
      mean = mu,
      sd   = sigma
    ), tolerance = 1e-12)
    expect_equal(moments[["mean"]], mu, tolerance = 1e-12)
    expect_equal(moments[["second"]], sigma^2 + mu^2, tolerance = 1e-12)
  }
})

test_that("step normalizer and CDF use stable tail interval probabilities", {

  skip_if_not(.has_native_selnorm_kernel())

  sei  <- 1
  yi   <- 9.5
  spec <- .test_step_spec(yi, sei)
  spec[["z_lower"]] <- c(9, 10, -Inf)
  spec[["z_upper"]] <- c(10, Inf, 9)
  spec[["obs_bin"]] <- 1L
  spec[["n_bins"]]  <- 3L

  omega <- matrix(c(1, 1e-300, 1e-300), nrow = 1)
  mean  <- matrix(0, nrow = 1, ncol = 1)
  sd    <- matrix(1, nrow = 1, ncol = 1)

  log_mass <- log(omega[1, ]) + vapply(
    seq_len(spec[["n_bins"]]),
    function(b) .test_interval_log_prob(
      spec[["z_lower"]][b],
      spec[["z_upper"]][b],
      0,
      1
    ),
    numeric(1)
  )
  expected_log_norm <- .test_logsumexp(log_mass)

  expected_log_cdf <- .test_logsumexp(log(omega[1, ]) + c(
    .test_interval_log_prob(9, yi, 0, 1),
    .test_interval_log_prob(10, yi, 0, 1),
    .test_interval_log_prob(-Inf, 9, 0, 1)
  ))
  expected_cdf <- exp(expected_log_cdf - expected_log_norm)

  log_norm <- .selnorm_kernel_log_norm_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )
  log_pdf <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mean,
    sigma_num      = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )
  cdf <- .selnorm_kernel_cdf_matrix(
    q              = yi,
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )

  expect_true(is.finite(log_norm[1, 1]))
  expect_equal(log_norm[1, 1], expected_log_norm, tolerance = 1e-10)
  expect_equal(
    log_pdf[1, 1],
    stats::dnorm(yi, mean = 0, sd = 1, log = TRUE) - expected_log_norm,
    tolerance = 1e-10
  )
  expect_equal(cdf[1, 1], expected_cdf, tolerance = 1e-10)
})

test_that("negative effect direction uses the same signed kernel", {

  skip_if_not(.has_native_selnorm_kernel())

  yi_pos <- c(.08, .20, .45)
  yi_neg <- -yi_pos
  sei    <- c(.10, .16, .22)
  mu_pos <- matrix(c(.04, .12, .30, .02, .10, .25), nrow = 2, byrow = TRUE)
  sigma  <- matrix(c(.16, .22, .28, .18, .24, .30), nrow = 2, byrow = TRUE)
  omega  <- matrix(c(1, .6, .3, .2, 1, .8, .4, .1), nrow = 2, byrow = TRUE)

  spec_pos <- .test_step_spec(yi_pos, sei, effect_direction = "positive")
  spec_neg <- .test_step_spec(yi_neg, sei, effect_direction = "negative")

  out_pos <- .selnorm_kernel_loglik_matrix(
    yi_pos, mu_pos, sigma, sei = sei, omega = omega,
    selection_spec = spec_pos
  )
  out_neg <- .selnorm_kernel_loglik_matrix(
    yi_neg, -mu_pos, sigma, sei = sei, omega = omega,
    selection_spec = spec_neg
  )

  expect_equal(out_neg, out_pos, tolerance = 1e-12)
})

test_that("inactive kernel rows reduce to ordinary normal likelihood", {

  skip_if_not(.has_native_selnorm_kernel())

  yi     <- c(.10, .40)
  sei    <- c(.10, .20)
  mu     <- matrix(c(.02, .20, .05, .22), nrow = 2, byrow = TRUE)
  sigma  <- matrix(c(.15, .30, .18, .35), nrow = 2, byrow = TRUE)
  omega  <- matrix(c(.01, .02, .03, .04, 1, .6, .3, .1), nrow = 2, byrow = TRUE)
  weights <- c(1.5, .25)
  spec   <- .test_step_spec(yi, sei)

  out <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega,
    selection_spec = spec,
    kernel_mode    = c(SELKERNEL_NORMAL, SELKERNEL_STEP),
    weights        = weights
  )
  ref <- .test_step_reference(yi, mu, sigma, sei, omega, spec, weights = weights)
  ref[1, ] <- weights * stats::dnorm(yi, mean = mu[1, ], sd = sigma[1, ], log = TRUE)

  expect_equal(out, ref, tolerance = 1e-12)
})

test_that("step selected-normal RNG follows exact bin masses and moments", {

  skip_if_not(.has_native_selnorm_kernel())

  set.seed(101)
  yi    <- 0
  sei   <- .11
  spec  <- .test_step_spec(yi, sei)
  S     <- 30000L
  mu    <- matrix(.05, nrow = S, ncol = 1)
  sigma <- matrix(.23, nrow = S, ncol = 1)
  omega <- matrix(
    rep(c(1, .15, 2.5, .05), each = S),
    nrow = S,
    ncol = spec[["n_bins"]]
  )

  context <- spec
  context[["omega"]]       <- omega
  context[["alpha"]]       <- rep(0, S)
  context[["phack_kind"]]  <- rep(0L, S)
  context[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)

  out  <- .selection_step_rng_matrix(mu, sigma, sei, context)
  bins <- .selection_step_bin_from_z(out[, 1] / sei, spec[["p_cuts"]])

  mass <- numeric(spec[["n_bins"]])
  for (b in seq_len(spec[["n_bins"]])) {
    mass[b] <- omega[1, b] * .test_interval_prob(
      spec[["z_lower"]][b],
      spec[["z_upper"]][b],
      spec[["sign"]] * mu[1, 1] / sei,
      sigma[1, 1] / sei
    )
  }
  expected_prob <- mass / sum(mass)
  observed_prob <- as.numeric(prop.table(table(
    factor(bins, levels = seq_len(spec[["n_bins"]]))
  )))

  moments <- .test_step_moments_reference(
    mu    = matrix(mu[1, 1], nrow = 1),
    sigma = matrix(sigma[1, 1], nrow = 1),
    sei   = sei,
    omega = matrix(omega[1, ], nrow = 1),
    spec  = spec
  )

  expect_equal(observed_prob, expected_prob, tolerance = .01)
  expect_equal(mean(out[, 1]), moments[["mean"]][1, 1], tolerance = .006)
  expect_equal(mean(out[, 1]^2), moments[["second"]][1, 1], tolerance = .006)
})

test_that("step selected-normal RNG handles signed extreme tail truncation", {

  skip_if_not(.has_native_selnorm_kernel())

  set.seed(102)
  sei   <- .10
  spec  <- .test_step_spec(-.10, sei, effect_direction = "negative")
  S     <- 128L
  omega <- matrix(
    rep(c(1, 1e-12, 1e-12, 1e-12), each = S),
    nrow = S,
    ncol = spec[["n_bins"]]
  )

  context <- spec
  context[["omega"]]       <- omega
  context[["alpha"]]       <- rep(0, S)
  context[["phack_kind"]]  <- rep(0L, S)
  context[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)

  out <- .selection_step_rng_matrix(
    mean              = matrix(-.35, nrow = S, ncol = 1),
    sd                = matrix(.05, nrow = S, ncol = 1),
    sei               = sei,
    selection_context = context
  )
  signed_z <- spec[["sign"]] * out[, 1] / sei

  expect_true(all(is.finite(out)))
  expect_true(all(signed_z >= spec[["z_lower"]][1] - 1e-10))
  expect_true(all(out[, 1] < 0))
})

test_that("selected-normal native boundary handles structural zero omega", {

  skip_if_not(.has_native_selnorm_kernel())

  yi     <- c(.10, .20)
  sei    <- c(.10, .10)
  spec   <- .test_step_spec(yi, sei)
  mu     <- matrix(0, nrow = 1, ncol = 2)
  sigma  <- matrix(.20, nrow = 1, ncol = 2)
  omega_empty_bin <- matrix(c(1, 0, .5, .25), nrow = 1)
  omega_obs_bin   <- matrix(c(1, .5, 0, .25), nrow = 1)

  empty_bin_log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega_empty_bin,
    selection_spec = spec
  )
  obs_bin_log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu,
    sigma_num      = sigma,
    sei            = sei,
    omega          = omega_obs_bin,
    selection_spec = spec
  )

  expect_true(all(is.finite(empty_bin_log_lik)))
  expect_equal(obs_bin_log_lik[1, 1], -Inf)
  expect_true(is.finite(obs_bin_log_lik[1, 2]))
})

test_that("active p-hacking kernels error on unsupported non-likelihood paths", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- .10
  sei   <- .10
  spec  <- .test_step_spec(yi, sei)
  mean  <- matrix(c(0, .02), nrow = 2, ncol = 1)
  sd    <- matrix(.20, nrow = 2, ncol = 1)
  omega <- matrix(1, nrow = 2, ncol = spec[["n_bins"]])
  alpha <- c(.50, 0)
  phack <- c(1L, 0L)
  mode  <- c(SELKERNEL_STEP_PHACK_POWER, SELKERNEL_STEP)

  log_lik <- .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mean,
    sigma_num      = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec,
    alpha          = alpha,
    phack_kind     = phack,
    kernel_mode    = mode
  )
  expect_true(all(is.finite(log_lik)))

  expect_error(
    .selnorm_kernel_cdf_matrix(
      q              = yi,
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )
  expect_error(
    .selnorm_kernel_moments_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )
  expect_error(
    .selnorm_kernel_rng_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )
  expect_error(
    .selnorm_kernel_weighted_summary(
      yi             = yi,
      mean           = mean,
      sd             = sd,
      sei            = sei,
      psis_weights   = matrix(c(.5, .5), nrow = 2, ncol = 1),
      omega          = omega,
      selection_spec = spec,
      alpha          = alpha,
      phack_kind     = phack,
      kernel_mode    = mode
    ),
    "active p-hacking"
  )

  context <- spec
  context[["omega"]]          <- omega
  context[["alpha"]]          <- alpha
  context[["phack_kind"]]     <- phack
  context[["kernel_mode"]]    <- mode
  context[["bias_indicator"]] <- c(1L, 1L)
  context[["use_normal"]]     <- c(FALSE, FALSE)

  setup <- list(
    mu        = c(0, .02),
    selection = context
  )

  expect_error(
    .funnel_selected_cdf(
      q                = yi,
      rows             = 1:2,
      se               = 0,
      total_sd         = c(.20, .20),
      setup            = setup,
      effect_direction = "positive"
    ),
    "active p-hacking"
  )
})

test_that("funnel selected-normal CDF is finite at zero standard error", {

  skip_if_not(.has_native_selnorm_kernel())

  spec <- .test_step_spec(yi = c(.10, .25), sei = c(.10, .15))

  selection_context <- spec
  selection_context[["omega"]]          <- matrix(c(1, .70, .35, .20), nrow = 1)
  selection_context[["alpha"]]          <- 0
  selection_context[["phack_kind"]]     <- 0L
  selection_context[["kernel_mode"]]    <- SELKERNEL_STEP
  selection_context[["bias_indicator"]] <- 1L
  selection_context[["use_normal"]]     <- FALSE

  setup <- list(
    mu                = 0,
    tau               = .20,
    PET               = 0,
    PEESE             = 0,
    is_weightfunction = TRUE,
    selection         = selection_context
  )

  out <- .funnel_model_averaged_cdf(
    q                = 0,
    se               = 0,
    setup            = setup,
    effect_direction = "positive"
  )

  expected <- (.20 * stats::pnorm(0, mean = 0, sd = .20)) /
    (stats::pnorm(0, mean = 0, sd = .20, lower.tail = FALSE) +
      .20 * stats::pnorm(0, mean = 0, sd = .20))

  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("selection-only prior_bias wrapper is delegated to BayesTools priors", {

  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  object <- bselmodel(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    prior_bias                = BayesTools::prior_bias(selection = selection),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  priors <- .create_fit_priors(object[["data"]], object[["priors"]])

  expect_true("bias" %in% names(priors))
  expect_true(BayesTools::is_prior_bias(priors[["bias"]]))
  expect_equal(BayesTools::JAGS_to_monitor(priors), c("mu", "tau", "omega"))
})

test_that("selection-only prior_bias wrapper is visible to zplot thresholds", {

  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  object <- bselmodel(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    prior_bias                = BayesTools::prior_bias(selection = selection),
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )

  expect_equal(.zplot_threshold(object[["priors"]]), stats::qnorm(1 - .025))
})

test_that("native zplot density matches selected-normal R reference", {

  skip_if_not(.has_native_zplot_density(selection = TRUE))

  S        <- 4L
  K        <- 3L
  sei      <- c(.07, .13, .29)
  mu       <- matrix(c(
    -.20, .05, .30,
     .10, .12, .18,
     .40, .00, -.15,
     .05, -.25, .20
  ), nrow = S, byrow = TRUE)
  tau      <- matrix(c(
    .05, .07, .10,
    .08, .06, .11,
    .03, .09, .12,
    .04, .10, .15
  ), nrow = S, byrow = TRUE)
  z_values <- c(
    -3,
    0,
    stats::qnorm(.05, lower.tail = FALSE),
    stats::qnorm(.025, lower.tail = FALSE),
    3
  )
  spec     <- .test_step_spec(yi = c(.10, .25, .40), sei = sei)

  selection <- spec
  selection[["omega"]]       <- matrix(c(
    1, .7, .2, .1,
    .2, .4, .8, 1.5,
    1.5, .2, 2.0, .05,
    1, 1, 1, 1
  ), nrow = S, byrow = TRUE)
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- c(
    SELKERNEL_STEP,
    SELKERNEL_NORMAL,
    SELKERNEL_STEP,
    SELKERNEL_NORMAL
  )
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  reference_density <- function(extrapolate) {

    sei_mat     <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
    total_sd    <- sqrt(tau^2 + sei_mat^2)
    densities   <- matrix(NA_real_, nrow = S, ncol = length(z_values))
    const_rows  <- if (extrapolate) which(!selection[["use_normal"]]) else integer(0)
    weight_rows <- if (extrapolate) integer(0) else which(!selection[["use_normal"]])

    log_norm <- NULL
    if (length(const_rows) > 0L) {
      log_norm <- .selection_step_log_norm_matrix(
        mean              = mu[const_rows, , drop = FALSE],
        sd                = total_sd[const_rows, , drop = FALSE],
        sei               = sei,
        selection_context = .selection_context_subset_rows(selection, const_rows)
      )
    }

    for (ell in seq_along(z_values)) {
      y_seq <- z_values[ell] * sei
      dens  <- stats::dnorm(
        matrix(y_seq, nrow = S, ncol = K, byrow = TRUE),
        mean = mu,
        sd   = total_sd
      ) * sei_mat

      if (length(const_rows) > 0L) {
        dens[const_rows, ] <- dens[const_rows, , drop = FALSE] / exp(log_norm)
      }
      if (length(weight_rows) > 0L) {
        selection_weight <- .selection_context_subset_rows(selection, weight_rows)
        log_pdf <- .selection_step_logpdf_matrix(
          y                 = y_seq,
          mean              = mu[weight_rows, , drop = FALSE],
          sd                = total_sd[weight_rows, , drop = FALSE],
          sei               = sei,
          selection_context = selection_weight
        )
        dens[weight_rows, ] <- exp(log_pdf) * sei_mat[weight_rows, , drop = FALSE]
      }

      densities[, ell] <- rowMeans(dens)
    }

    return(densities)
  }

  for (extrapolate in c(FALSE, TRUE)) {
    expect_equal(
      .zplot_density_vectorized(
        z_sequence       = z_values,
        mu_samples       = mu,
        tau_within       = tau,
        sei              = sei,
        selection        = selection,
        extrapolate      = extrapolate,
        effect_direction = "positive"
      ),
      reference_density(extrapolate),
      tolerance = 1e-12
    )
  }
})


test_that("zplot threshold separates fitted, extrapolated, and missing-study targets", {

  skip_if_not(.has_native_selnorm_kernel())

  S         <- 4L
  K         <- 3L
  sei       <- c(.07, .13, .29)
  mu        <- matrix(c(
    -.20, .05, .30,
     .10, .12, .18,
     .40, .00, -.15,
     .05, -.25, .20
  ), nrow = S, byrow = TRUE)
  tau       <- matrix(c(
    .05, .07, .10,
    .08, .06, .11,
    .03, .09, .12,
    .04, .10, .15
  ), nrow = S, byrow = TRUE)
  threshold <- stats::qnorm(.025, lower.tail = FALSE)
  spec      <- .test_step_spec(yi = c(.10, .25, .40), sei = sei)

  selection <- spec
  selection[["omega"]]       <- matrix(c(
    1, .7, .2, .1,
    .2, .4, .8, 1.5,
    1.5, .2, 2.0, .05,
    1, 1, 1, 1
  ), nrow = S, byrow = TRUE)
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- c(
    SELKERNEL_STEP,
    SELKERNEL_NORMAL,
    SELKERNEL_STEP,
    SELKERNEL_NORMAL
  )
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau^2 + sei_mat^2)
  q_upper  <- threshold * sei
  q_lower  <- -threshold * sei

  unselected_thresholds <- stats::pnorm(
    matrix(q_upper, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu,
    sd         = total_sd,
    lower.tail = FALSE
  ) + stats::pnorm(
    matrix(q_lower, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu,
    sd         = total_sd,
    lower.tail = TRUE
  )

  step_rows      <- which(!selection[["use_normal"]])
  selection_step <- .selection_context_subset_rows(selection, step_rows)
  log_norm       <- .selection_step_log_norm_matrix(
    mean              = mu[step_rows, , drop = FALSE],
    sd                = total_sd[step_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_step
  )
  expected_weights <- matrix(1, nrow = S, ncol = K)
  expected_weights[step_rows, ] <- exp(-log_norm)

  selected_thresholds <- unselected_thresholds
  selected_thresholds[step_rows, ] <- .selection_step_cdf_matrix(
    q                 = q_upper,
    mean              = mu[step_rows, , drop = FALSE],
    sd                = total_sd[step_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_step,
    lower.tail        = FALSE
  ) + .selection_step_cdf_matrix(
    q                 = q_lower,
    mean              = mu[step_rows, , drop = FALSE],
    sd                = total_sd[step_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_step,
    lower.tail        = TRUE
  )

  extrapolated <- .zplot_threshold_vectorized(
    z_threshold      = threshold,
    mu_samples       = mu,
    tau_within       = tau,
    sei              = sei,
    selection        = selection,
    extrapolate      = TRUE,
    effect_direction = "positive"
  )
  fitted <- .zplot_threshold_vectorized(
    z_threshold      = threshold,
    mu_samples       = mu,
    tau_within       = tau,
    sei              = sei,
    selection        = selection,
    extrapolate      = FALSE,
    effect_direction = "positive"
  )

  expect_equal(
    extrapolated[["EDR"]],
    rowSums(unselected_thresholds * expected_weights) / rowSums(expected_weights),
    tolerance = 1e-12
  )
  expect_equal(extrapolated[["weights"]], rowMeans(expected_weights), tolerance = 1e-12)
  expect_equal(fitted[["EDR"]], rowMeans(selected_thresholds), tolerance = 1e-12)
  expect_equal(fitted[["weights"]], rowMeans(expected_weights), tolerance = 1e-12)
})


test_that("native regplot selection intervals match R reference", {

  skip_if_not(.has_native_regplot_selection_mixture())

  set.seed(131)
  S            <- 15L
  K            <- 4L
  se           <- .14
  probs        <- c(.025, .975)
  mean_samples <- matrix(rnorm(S * K, 0, .35), nrow = S, ncol = K)
  sd_samples   <- matrix(runif(S * K, .02, .30), nrow = S, ncol = K)
  sd_samples[c(3, 11)] <- sqrt(.Machine$double.eps) / 10
  spec         <- .test_step_spec(yi = seq_len(K) / 10, sei = rep(se, K))

  selection <- spec
  selection[["omega"]]       <- matrix(
    runif(S * spec[["n_bins"]], .02, 3),
    nrow = S,
    ncol = spec[["n_bins"]]
  )
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection[["kernel_mode"]][seq(2, S, by = 4)] <- SELKERNEL_NORMAL
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  setup <- list(
    is_weightfunction = !selection[["use_normal"]],
    selection         = selection
  )

  native <- .regplot_selnorm_mixture_interval_quantiles(
    mean_samples      = mean_samples,
    sd_samples        = sd_samples,
    se                = se,
    probs             = probs,
    selection_context = selection
  )
  ref <- .regplot_selection_mixture_interval_quantiles_r(
    mean_samples     = mean_samples,
    sd_samples       = sd_samples,
    se               = se,
    probs            = probs,
    setup            = setup,
    effect_direction = "positive"
  )

  expect_equal(native[["lower"]], ref[["lower"]], tolerance = 1e-5)
  expect_equal(native[["upper"]], ref[["upper"]], tolerance = 1e-5)
})

test_that("native weighted selected-normal summary matches matrix reductions", {

  skip_if_not(.has_native_selnorm_kernel())

  set.seed(141)
  S    <- 17L
  K    <- 5L
  yi   <- -seq_len(K) / 20
  sei  <- seq(.06, .21, length.out = K)
  spec <- .test_step_spec(yi = yi, sei = sei, effect_direction = "negative")
  mean <- matrix(rnorm(S * K, 0, .25), nrow = S, ncol = K)
  sd   <- matrix(runif(S * K, .04, .35), nrow = S, ncol = K)
  psis <- matrix(rexp(S * K), nrow = S, ncol = K)
  psis <- sweep(psis, 2, colSums(psis), "/")

  selection <- spec
  selection[["omega"]]       <- matrix(
    runif(S * spec[["n_bins"]], .03, 2),
    nrow = S,
    ncol = spec[["n_bins"]]
  )
  selection[["alpha"]]       <- rep(0, S)
  selection[["phack_kind"]]  <- rep(0L, S)
  selection[["kernel_mode"]] <- rep(SELKERNEL_STEP, S)
  selection[["kernel_mode"]][seq(3, S, by = 5)] <- SELKERNEL_NORMAL
  selection[["use_normal"]] <- selection[["kernel_mode"]] == SELKERNEL_NORMAL

  native <- .selection_step_weighted_summary(
    yi                = yi,
    mean              = mean,
    sd                = sd,
    sei               = sei,
    psis_weights      = psis,
    selection_context = selection
  )
  cdf     <- .selection_step_cdf_matrix(yi, mean, sd, sei, selection)
  moments <- .selection_step_moments_matrix(mean, sd, sei, selection)

  expect_equal(native[["cdf"]], colSums(psis * cdf), tolerance = 1e-12)
  expect_equal(native[["mean"]], colSums(psis * moments[["mean"]]), tolerance = 1e-12)
  expect_equal(native[["second"]], colSums(psis * moments[["second"]]), tolerance = 1e-12)
})

test_that("p-hacking priors are deferred before RoBMA model construction", {

  phacking <- BayesTools::prior_phacking(form = "linear")
  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  expect_false("prior_phacking" %in% getNamespaceExports("RoBMA"))
  expect_false("prior_bias" %in% getNamespaceExports("RoBMA"))
  expect_error(
    .use_plot_prior_list_dispatch(phacking),
    "not enabled"
  )
  expect_error(
    bselmodel(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      prior_bias                = phacking,
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    bselmodel(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      prior_bias                = BayesTools::prior_bias(
        selection = selection,
        phacking  = phacking
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    RoBMA(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      prior_bias                = phacking,
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    RoBMA(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      prior_bias                = BayesTools::prior_bias(
        selection = selection,
        phacking  = phacking
      ),
      prior_unit_information_sd = 1,
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "not enabled"
  )
  expect_error(
    .selection_spec(
      priors           = list(outcome = list(bias = phacking)),
      yi               = c(.1, .2),
      sei              = c(.1, .1),
      effect_direction = "positive"
    ),
    "not enabled"
  )
  expect_error(
    .selection_spec(
      priors = list(outcome = list(
        bias = BayesTools::prior_bias(selection = selection, phacking = phacking)
      )),
      yi               = c(.1, .2),
      sei              = c(.1, .1),
      effect_direction = "positive"
    ),
    "not enabled"
  )
})
