context("Selection kernel")
skip_on_cran()

test_that("selection-only prior_bias wrapper is delegated to BayesTools priors", {

  selection <- BayesTools::prior_weightfunction(
    "one-sided", c(.025), BayesTools::wf_fixed(c(1, .5))
  )

  object <- bselmodel(
    yi                        = c(.1, .2),
    sei                       = c(.1, .1),
    measure                   = "SMD",
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
    measure                   = "SMD",
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
        selection_context = BayesTools::selection_context_subset_rows(selection, const_rows)
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
        selection_weight <- BayesTools::selection_context_subset_rows(selection, weight_rows)
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
  selection_step <- BayesTools::selection_context_subset_rows(selection, step_rows)
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
      measure                   = "SMD",
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
      measure                   = "SMD",
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
      measure                   = "SMD",
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
      measure                   = "SMD",
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
