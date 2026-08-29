context("Selection kernel")
skip_on_cran()

test_that("selection likelihood exposes only exact and approximate targets", {

  expect_error(
    bselmodel(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection_likelihood      = "conditional",
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "one of.*exact.*approximate"
  )
  expect_error(
    RoBMA(
      yi                        = c(.1, .2),
      sei                       = c(.1, .1),
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection_likelihood      = "conditional",
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "one of.*exact.*approximate"
  )
})

test_that("selection model titles omit the likelihood implementation", {

  args <- list(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact       <- do.call(bselmodel, c(args, selection_likelihood = "exact"))
  approximate <- do.call(
    bselmodel,
    c(args, selection_likelihood = "approximate")
  )

  expect_identical(
    .summary.brma_model_names(exact),
    "Bayesian Random-Effects Selection Model (k = 2)"
  )
  expect_identical(
    .summary.brma_model_names(approximate),
    "Bayesian Random-Effects Selection Model (k = 2)"
  )
})

test_that("selection reference weights are structural convergence parameters", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_cumulative(c(1, 1, 1))
  )

  object <- bselmodel(
    yi          = c(.1, .2, .3),
    sei         = c(.1, .1, .1),
    measure     = "SMD",
    prior_bias  = prior_bias,
    only_priors = TRUE,
    silent      = TRUE
  )

  expect_identical(
    .convergence_structural_parameters(object[["priors"]]),
    c("omega[1]", "omega[0,0.025]")
  )

  object[["priors"]][["outcome"]][["bias"]] <- NULL
  expect_length(.convergence_structural_parameters(object[["priors"]]), 0L)
})

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
    selection_likelihood      = "approximate",
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
    selection_likelihood      = "approximate",
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


test_that("RoBMA defaults to exact finite-vector selection", {

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = .025,
    weights = BayesTools::wf_fixed(c(1, .5))
  )
  args <- list(
    yi                        = c(.10, .20, .05, .15),
    sei                       = rep(.10, 4L),
    cluster                   = c("a", "a", "b", "b"),
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
  exact <- do.call(RoBMA, args)

  expect_true(.is_data_exact_selection(exact[["data"]]))
  expect_identical(exact[["selection_likelihood"]][["type"]], "exact")
  expect_identical(
    exact[["selection_likelihood"]][["target"]],
    "finite_vector_product_selection"
  )
  expect_identical(
    .data_exact_selection_setup(exact[["data"]])[["row_blocks"]],
    list(1:2, 3:4)
  )
  exact_priors <- .create_fit_priors(exact[["data"]], exact[["priors"]])
  exact_syntax <- .create_model_syntax(exact[["data"]], exact[["priors"]])
  expect_null(exact_priors[["gamma"]])
  expect_match(exact_syntax, "dselnorm_mnorm_step", fixed = TRUE)
  expect_match(exact_syntax, "sel_kernel_mode_active", fixed = TRUE)
  expect_false(grepl("gamma[", exact_syntax, fixed = TRUE))
  expect_false(grepl("dselnorm_step_switch", exact_syntax, fixed = TRUE))

  approximate <- do.call(
    RoBMA,
    c(args, selection_likelihood = "approximate")
  )
  approximate_syntax <- .create_model_syntax(
    approximate[["data"]],
    approximate[["priors"]]
  )
  expect_false(.is_data_exact_selection(approximate[["data"]]))
  expect_identical(
    approximate[["selection_likelihood"]][["target"]],
    "row_selected_normal_conditional_on_random_effects"
  )
  expect_match(approximate_syntax, "dselnorm_step_switch", fixed = TRUE)
  expect_match(approximate_syntax, "gamma[cluster[i]]", fixed = TRUE)
})


test_that("RoBMA exact selection preserves explicit target boundaries", {

  weighted_args <- list(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    weights                   = c(1, .5, 1),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  expect_error(
    do.call(RoBMA, weighted_args),
    "'weights' are unavailable with 'selection_likelihood = \"exact\"'",
    fixed = TRUE
  )
  approximate <- do.call(
    RoBMA,
    c(weighted_args, selection_likelihood = "approximate")
  )
  expect_identical(
    approximate[["selection_likelihood"]][["type"]],
    "approximate"
  )

  no_selection <- RoBMA(
    yi                        = c(.10, .20, .05),
    sei                       = rep(.10, 3L),
    measure                   = "SMD",
    model_type                = "PP",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  expect_false(.is_priors_weightfunction(no_selection[["priors"]]))
  expect_null(no_selection[["selection_likelihood"]])
  expect_null(attr(
    no_selection[["data"]],
    "selection_likelihood",
    exact = TRUE
  ))
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

  wrapped_safe_prior <- BayesTools::prior_bias(selection = safe_prior)
  expect_silent(
    wrapped_safe_spec <- .selection_spec(
      priors           = list(outcome = list(bias = wrapped_safe_prior)),
      yi               = yi,
      sei              = sei,
      effect_direction = "positive",
      signed_data      = TRUE
    )
  )
  expect_true(wrapped_safe_spec[["telescope_probabilities"]])
  expect_identical(wrapped_safe_spec[["jags_data"]][["sel_telescope_probabilities"]], 1L)

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
    backend[["branch_kernel_mode"]],
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
    selection_likelihood      = "approximate",
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
    selection_likelihood      = "approximate",
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

  custom_samples <- posterior_samples
  colnames(custom_samples) <- sub("^omega", "custom.omega+beta", colnames(custom_samples))
  custom_samples <- cbind(
    custom_samples,
    "omega[1]" = 101:103,
    "omega[2]" = 201:203
  )
  selection_spec <- list(jags_omega = "custom.omega+beta", n_bins = 10L)
  custom_omega   <- .extract_selection_omega_samples(custom_samples, selection_spec)

  expect_equal(colnames(custom_omega), paste0("custom.omega+beta[", 1:10, "]"))
  expect_equal(custom_omega[, 2], custom_samples[, "custom.omega+beta[2]"])
  expect_false(any(custom_omega[, 1] == custom_samples[, "omega[1]"]))

  missing_custom <- custom_samples[, grepl("^omega\\[", colnames(custom_samples)), drop = FALSE]
  expect_error(
    .extract_selection_omega_samples(missing_custom, selection_spec),
    "custom.omega\\+beta"
  )
})

test_that("exact selection constructors marginalize Gaussian dependence", {

  cluster_object <- bselmodel(
    yi                        = c(.10, .20, .05, .15),
    sei                       = rep(.10, 4L),
    cluster                   = c("a", "a", "b", "b"),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  cluster_fit_priors <- .create_fit_priors(
    cluster_object[["data"]],
    cluster_object[["priors"]]
  )
  cluster_syntax <- .create_model_syntax(
    cluster_object[["data"]],
    cluster_object[["priors"]]
  )

  expect_true(.is_data_exact_selection(cluster_object[["data"]]))
  expect_identical(cluster_object[["selection_likelihood"]][["type"]], "exact")
  expect_identical(
    .data_exact_selection_setup(cluster_object[["data"]])[["row_blocks"]],
    list(1:2, 3:4)
  )
  expect_null(cluster_fit_priors[["gamma"]])
  expect_match(cluster_syntax, "dselnorm_mnorm_step", fixed = TRUE)
  expect_false(grepl("gamma[", cluster_syntax, fixed = TRUE))

  data <- data.frame(study = factor(c("a", "a", "b")))
  V <- matrix(
    c(.010, .003, 0, .003, .014, 0, 0, 0, .012),
    nrow = 3L,
    byrow = TRUE
  )
  mv_object <- bselmodel.mv(
    yi                        = c(.10, .20, .05),
    V                         = V,
    random                    = ~ diag(1 | study),
    data                      = data,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  random_terms <- mv_object[["formula_design"]][["mu"]][["random_effects"]]

  expect_s3_class(mv_object, "bselmodel.mv")
  expect_true(.is_data_exact_selection(mv_object[["data"]]))
  expect_true(all(vapply(
    random_terms,
    function(term) identical(term[["compile_mode"]], "marginalized"),
    logical(1L)
  )))
  expect_match(
    .create_model_syntax(mv_object[["data"]], mv_object[["priors"]]),
    "sel_exact_random_block_1_lower",
    fixed = TRUE
  )

  expect_error(
    bselmodel.mv(
      yi                        = c(.10, .20, .05),
      V                         = V,
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection_likelihood      = "approximate",
      known_v_parameterization  = "whitened",
      only_priors               = TRUE,
      silent                    = TRUE
    ),
    "requires.*latent"
  )

  approximate_auto <- bselmodel.mv(
    yi                        = c(.10, .20, .05),
    V                         = V,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "approximate",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  known_V <- .data_known_v_data(approximate_auto[["data"]])
  expect_identical(.known_v_effective_backend(known_V), "latent")
  expect_identical(.known_v_requested_parameterization(known_V), "auto")
})


test_that("exact selection prediction partitions covariance without known-V metadata", {

  prediction_data <- list(outcome = list(sei = c(.10, .15, .20)))
  random_terms <- list(list(
    block_name = "study",
    group_map  = c(1L, 1L, 2L)
  ))

  expect_identical(
    .selection_exact_dependency_blocks(prediction_data, random_terms),
    list(1:2, 3L)
  )
})


test_that("exact selection covariance batches preserve multilevel algebra", {

  object <- bselmodel(
    yi                        = c(.10, .20),
    sei                       = c(.10, .15),
    cluster                   = c("a", "a"),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    only_priors               = TRUE,
    silent                    = TRUE
  )
  tau_within <- rbind(c(.20, .30), c(.25, .35))
  tau_between <- rbind(c(.40, .50), c(.45, .55))
  setup <- list(
    data          = object[["data"]],
    S             = 2L,
    tau_within    = tau_within,
    tau_between   = tau_between,
    is_multilevel = TRUE
  )

  observed <- .selection_exact_covariance_lower(setup, rows = 1:2)
  expected_covariance <- array(NA_real_, dim = c(2L, 2L, 2L))
  for (draw in seq_len(2L)) {
    expected_covariance[draw, , ] <- diag(c(.10, .15)^2) +
      diag(tau_within[draw, ]^2) +
      outer(tau_between[draw, ], tau_between[draw, ])
  }
  expected <- t(vapply(seq_len(2L), function(draw) {
    expected_covariance[draw, , ][
      cbind(c(1L, 2L, 2L), c(1L, 1L, 2L))
    ]
  }, numeric(3L)))

  expect_equal(observed, expected)

  response_setup <- .predict_exact_selection_response_setup(
    context        = list(
      object             = object,
      posterior_samples  = matrix(0, nrow = 2L, ncol = 1L),
      conditioning_depth = "marginal",
      new_data            = object[["data"]],
      known_V_new         = NULL,
      K                   = 2L,
      is_known_v          = FALSE,
      outcome_data        = object[["data"]][["outcome"]],
      random_mv           = FALSE,
      is_multilevel       = TRUE,
      same_data           = TRUE
    ),
    location_state = list(fixed_mu = matrix(0, nrow = 2L, ncol = 2L)),
    scale_state    = list(within = tau_within, between = tau_between)
  )
  expect_equal(response_setup[["covariance"]], expected_covariance)
})

test_that("exact bivariate selection kernel matches rectangle integration", {

  y     <- c(.30, .50)
  mu    <- c(.10, .15)
  sigma <- matrix(c(.09, .03, .03, .16), 2L, 2L)
  sei   <- c(.20, .25)
  z     <- stats::qnorm(.025, lower.tail = FALSE)
  omega <- c(.4, 1)
  plan  <- BayesTools::selection_likelihood_plan(
    block_sizes         = 2L,
    points_per_scramble = 16384L,
    scrambles           = 8L,
    seed                 = 5L,
    relative_tolerance   = .01
  )
  lower <- matrix(
    sigma[cbind(c(1L, 2L, 2L), c(1L, 1L, 2L))],
    nrow = 1L
  )
  actual <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    y, matrix(mu, nrow = 1L), lower, sei, matrix(omega, nrow = 1L),
    c(z, -Inf), c(Inf, z), c(2L, 1L), 1L, TRUE, 1L,
    as.double(plan[["designs"]][["2"]]), 16384L, 8L, .01,
    PACKAGE = "RoBMA"
  )
  reflected <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    -y, matrix(-mu, nrow = 1L), lower, sei, matrix(omega, nrow = 1L),
    c(z, -Inf), c(Inf, z), c(2L, 1L), -1L, TRUE, 1L,
    as.double(plan[["designs"]][["2"]]), 16384L, 8L, .01,
    PACKAGE = "RoBMA"
  )

  threshold <- z * sei
  rectangles <- list(
    list(lower = threshold, upper = c(Inf, Inf), weight = .4^2),
    list(lower = c(threshold[[1L]], -Inf),
         upper = c(Inf, threshold[[2L]]), weight = .4),
    list(lower = c(-Inf, threshold[[2L]]),
         upper = c(threshold[[1L]], Inf), weight = .4),
    list(lower = c(-Inf, -Inf), upper = threshold, weight = 1)
  )
  normalizer <- sum(vapply(rectangles, function(rectangle) {
    rectangle[["weight"]] * as.numeric(mvtnorm::pmvnorm(
      lower = rectangle[["lower"]],
      upper = rectangle[["upper"]],
      mean  = mu,
      sigma = sigma
    ))
  }, numeric(1L)))
  expected <- mvtnorm::dmvnorm(y, mean = mu, sigma = sigma, log = TRUE) +
    log(prod(omega)) - log(normalizer)

  expect_equal(actual[["log_density"]], expected, tolerance = 5e-4)
  expect_equal(reflected[["log_density"]], actual[["log_density"]],
               tolerance = 1e-12)
  expect_lt(actual[["relative_mcse"]], 5e-4)
})


test_that("exact diagonal selection kernel reduces to analytic row factors", {

  y     <- c(.30, .50)
  mu    <- c(.10, .15)
  sd    <- c(.30, .40)
  sei   <- c(.20, .25)
  z     <- stats::qnorm(.025, lower.tail = FALSE)
  omega <- c(.4, 1)
  plan  <- BayesTools::selection_likelihood_plan(
    block_sizes         = 2L,
    points_per_scramble = 8L,
    scrambles           = 2L,
    seed                 = 9L,
    relative_tolerance   = .01
  )
  actual <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    y, matrix(mu, nrow = 1L), matrix(c(sd[[1L]]^2, 0, sd[[2L]]^2), nrow = 1L),
    sei, matrix(omega, nrow = 1L), c(z, -Inf), c(Inf, z), c(2L, 1L),
    1L, TRUE, 1L, as.double(plan[["designs"]][["2"]]), 8L, 2L, .01,
    PACKAGE = "RoBMA"
  )

  threshold       <- z * sei
  significant     <- y >= threshold
  observed_weight <- ifelse(significant, omega[[1L]], omega[[2L]])
  normalizer      <- omega[[1L]] * stats::pnorm(
    threshold,
    mean       = mu,
    sd         = sd,
    lower.tail = FALSE
  ) + omega[[2L]] * stats::pnorm(threshold, mean = mu, sd = sd)
  expected <- sum(
    stats::dnorm(y, mean = mu, sd = sd, log = TRUE) +
      log(observed_weight) - log(normalizer)
  )

  expect_equal(actual[["log_density"]], expected, tolerance = 1e-13)
  expect_identical(actual[["relative_mcse"]], 0)
})


test_that("exact selected multivariate response RNG matches region masses", {

  set.seed(41)
  S     <- 6000L
  mu    <- c(.10, .15)
  sigma <- matrix(c(.09, .03, .03, .16), 2L, 2L)
  sei   <- c(.20, .25)
  omega <- c(40, 100)
  context <- list(
    omega       = matrix(omega, nrow = S, ncol = 2L, byrow = TRUE),
    kernel_mode = rep(SELKERNEL_STEP, S),
    use_normal  = rep(FALSE, S),
    p_cuts      = c(0, .025, 1),
    sign        = 1L
  )
  covariance <- array(
    rep(sigma, each = S),
    dim = c(S, 2L, 2L)
  )
  draws <- .outcome_rng.selnorm_mvn(
    mu_samples         = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
    covariance_samples = covariance,
    sei                = sei,
    selection_context  = context,
    dependency_blocks  = list(1:2)
  )

  threshold <- stats::qnorm(.025, lower.tail = FALSE) * sei
  region <- 1L + (draws[, 1L] >= threshold[[1L]]) +
    2L * (draws[, 2L] >= threshold[[2L]])
  rectangles <- list(
    list(lower = c(-Inf, -Inf), upper = threshold, weight = 1),
    list(lower = c(threshold[[1L]], -Inf),
         upper = c(Inf, threshold[[2L]]), weight = .4),
    list(lower = c(-Inf, threshold[[2L]]),
         upper = c(threshold[[1L]], Inf), weight = .4),
    list(lower = threshold, upper = c(Inf, Inf), weight = .4^2)
  )
  expected <- vapply(rectangles, function(rectangle) {
    rectangle[["weight"]] * as.numeric(mvtnorm::pmvnorm(
      lower = rectangle[["lower"]],
      upper = rectangle[["upper"]],
      mean  = mu,
      sigma = sigma
    ))
  }, numeric(1L))
  expected <- expected / sum(expected)

  expect_equal(as.numeric(prop.table(table(factor(region, levels = 1:4)))),
               expected, tolerance = .025)
})


test_that("exact selected response RNG preserves singular covariance support", {

  S <- 4L
  loading <- c(.2, -.3)
  covariance_matrix <- tcrossprod(loading)
  covariance <- array(
    rep(covariance_matrix, each = S),
    dim = c(S, 2L, 2L)
  )
  context <- list(
    omega       = matrix(1, nrow = S, ncol = 1L),
    kernel_mode = rep(SELKERNEL_NORMAL, S),
    use_normal  = rep(TRUE, S),
    p_cuts      = c(0, 1),
    sign        = 1L
  )
  mean <- matrix(c(.1, -.2), nrow = S, ncol = 2L, byrow = TRUE)

  set.seed(842)
  actual <- .outcome_rng.selnorm_mvn(
    mu_samples         = mean,
    covariance_samples = covariance,
    sei                = c(.1, .1),
    selection_context  = context,
    dependency_blocks  = list(1:2)
  )
  set.seed(842)
  factor <- .covariance_sampling_factor(
    .covariance_factorization(covariance_matrix)
  )
  expected <- matrix(NA_real_, nrow = S, ncol = 2L)
  for (draw in seq_len(S)) {
    expected[draw, ] <- mean[draw, ] +
      as.vector(stats::rnorm(2L) %*% factor)
  }

  expect_equal(actual, expected, tolerance = 1e-15)

  indefinite <- array(
    rep(matrix(c(1, 2, 2, 1), nrow = 2L), each = S),
    dim = c(S, 2L, 2L)
  )
  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = mean,
      covariance_samples = indefinite,
      sei                = c(.1, .1),
      selection_context  = context,
      dependency_blocks  = list(1:2)
    ),
    "Selected response covariance is not positive semidefinite\\."
  )
})


test_that("exact selected response RNG rejects asymmetric covariance", {

  context <- list(
    omega       = matrix(1, nrow = 1L, ncol = 1L),
    kernel_mode = SELKERNEL_NORMAL,
    use_normal  = TRUE,
    p_cuts      = c(0, 1),
    sign        = 1L
  )
  covariance <- array(c(1, .1, .2, 1), dim = c(1L, 2L, 2L))

  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = matrix(0, nrow = 1L, ncol = 2L),
      covariance_samples = covariance,
      sei                = c(1, 1),
      selection_context  = context,
      dependency_blocks  = list(1:2)
    ),
    "Selected response covariance must be symmetric.",
    fixed = TRUE
  )
})


test_that("exact selected response RNG reports exhausted rejection sampling", {

  context <- list(
    omega       = matrix(c(0, 1), nrow = 1L),
    kernel_mode = SELKERNEL_STEP,
    use_normal  = FALSE,
    p_cuts      = c(0, .025, 1),
    sign        = 1L
  )
  expect_error(
    .outcome_rng.selnorm_mvn(
      mu_samples         = matrix(1, nrow = 1L),
      covariance_samples = array(0, dim = c(1L, 1L, 1L)),
      sei                = .1,
      selection_context  = context,
      dependency_blocks  = list(1L),
      max_attempts       = 3L
    ),
    paste0(
      "Exact selected response RNG was rejected by diagnostics: no ",
      "proposal was accepted in 3 attempts for a dependency block of size ",
      "1. Use 'bias_adjusted = TRUE' to draw responses before selection."
    ),
    fixed = TRUE
  )
})


test_that("fixed selection mixtures require exact branch indicators", {

  selection_spec <- list(
    fixed_omega = matrix(c(1, 0.5), ncol = 1L),
    n_bins      = 1L,
    jags_omega  = "omega"
  )
  samples <- matrix(
    c(1, 2),
    ncol     = 1L,
    dimnames = list(NULL, "bias_indicator")
  )

  expect_equal(
    .extract_selection_fixed_omega_samples(samples, selection_spec),
    matrix(c(1, 0.5), ncol = 1L,
           dimnames = list(NULL, "omega[1]"))
  )

  samples[1L, 1L] <- 1 + .Machine$double.eps
  expect_error(
    .extract_selection_fixed_omega_samples(samples, selection_spec),
    "integer-valued"
  )
})
