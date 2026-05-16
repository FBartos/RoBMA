context("Selection kernel")
skip_on_cran()

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

test_that("default RoBMA step normalizer matches independent Bem2011 reference", {

  skip_if_not(.has_native_selnorm_kernel())

  object <- RoBMA(
    yi           = Bem2011[["d"]],
    sei          = Bem2011[["se"]],
    measure      = "SMD",
    only_priors  = TRUE,
    silent       = TRUE
  )
  spec <- .selection_spec(
    priors           = object[["priors"]],
    yi               = object[["data"]][["outcome"]][["yi"]],
    sei              = object[["data"]][["outcome"]][["sei"]],
    effect_direction = attr(object[["data"]], "effect_direction")
  )

  set.seed(107)
  S     <- 4L
  K     <- length(Bem2011[["d"]])
  mean  <- matrix(stats::rnorm(S * K, mean = .12, sd = .08), nrow = S)
  sd    <- matrix(stats::runif(S * K, min = .09, max = .42), nrow = S)
  omega <- matrix(stats::runif(S * spec[["n_bins"]], min = .05, max = 1.50), nrow = S)

  log_norm <- .selnorm_kernel_log_norm_matrix(
    mean           = mean,
    sd             = sd,
    sei            = Bem2011[["se"]],
    omega          = omega,
    selection_spec = spec
  )

  expect_equal(
    log_norm,
    .test_step_log_norm_reference(mean, sd, Bem2011[["se"]], omega, spec),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(log_norm)))
})

test_that("trusted step normalizer handles nonmonotone weights and tail underflow", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- Bem2011[["d"]][1:3]
  sei   <- Bem2011[["se"]][1:3]
  spec  <- .test_step_spec(yi, sei)
  mean  <- matrix(c(.07, .15, .23, -.04, .05, .17), nrow = 2, byrow = TRUE)
  sd    <- matrix(c(.12, .19, .31, .10, .21, .28), nrow = 2, byrow = TRUE)
  omega <- matrix(c(.2, 1.4, .3, 1.1, 1.3, .1, 1.2, .2), nrow = 2, byrow = TRUE)

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    ),
    .test_step_log_norm_reference(mean, sd, sei, omega, spec),
    tolerance = 1e-12
  )

  tail_spec              <- .test_step_spec(0, 1)
  tail_spec[["n_bins"]]  <- 3L
  tail_spec[["z_lower"]] <- c(40, 39, -Inf)
  tail_spec[["z_upper"]] <- c(Inf, 40, 39)
  tail_spec[["obs_bin"]] <- 1L
  tail_spec[["telescope_probabilities"]] <- FALSE
  tail_omega             <- matrix(c(1e300, 0, 1e-300), nrow = 1)
  tail_mean              <- matrix(0, nrow = 1, ncol = 1)
  tail_sd                <- matrix(1, nrow = 1, ncol = 1)

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = tail_mean,
      sd             = tail_sd,
      sei            = 1,
      omega          = tail_omega,
      selection_spec = tail_spec
    )[1, 1],
    .test_step_log_norm_reference(tail_mean, tail_sd, 1, tail_omega, tail_spec)[1, 1],
    tolerance = 1e-10
  )

  cancellation_spec              <- .test_step_spec(0, 1)
  cancellation_spec[["n_bins"]]  <- 3L
  cancellation_spec[["z_lower"]] <- c(1e-8, 0, -Inf)
  cancellation_spec[["z_upper"]] <- c(Inf, 1e-8, 0)
  cancellation_spec[["obs_bin"]] <- 2L
  cancellation_omega             <- matrix(c(1, 1e16, 1), nrow = 1)

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = tail_mean,
      sd             = tail_sd,
      sei            = 1,
      omega          = cancellation_omega,
      selection_spec = cancellation_spec
    )[1, 1],
    .test_step_log_norm_reference(
      tail_mean, tail_sd, 1, cancellation_omega, cancellation_spec
    )[1, 1],
    tolerance = 1e-10
  )
})

test_that("step normalizer handles normal mode and full-support unit weights", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- Bem2011[["d"]][1:4]
  sei   <- Bem2011[["se"]][1:4]
  spec  <- .test_step_spec(yi, sei)
  mean  <- matrix(c(.05, .12, .18, .24, .03, .10, .20, .30), nrow = 2, byrow = TRUE)
  sd    <- matrix(c(.12, .16, .20, .28, .10, .18, .22, .35), nrow = 2, byrow = TRUE)
  omega <- matrix(c(1, .6, .6, .2, 1, 1.4, .7, .7), nrow = 2, byrow = TRUE)

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      kernel_mode    = SELKERNEL_NORMAL
    ),
    matrix(0, nrow = nrow(mean), ncol = ncol(mean)),
    tolerance = 1e-12
  )
  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = matrix(1, nrow = nrow(mean), ncol = spec[["n_bins"]]),
      selection_spec = spec
    ),
    matrix(0, nrow = nrow(mean), ncol = ncol(mean)),
    tolerance = 1e-12
  )
  expect_false(all(abs(.selnorm_kernel_log_norm_matrix(
    mean           = mean,
    sd             = sd,
    sei            = sei,
    omega          = omega,
    selection_spec = spec
  )) < 1e-12))
})

test_that("step normalizer is unchanged by sentinels and adjacent repeated bins", {

  skip_if_not(.has_native_selnorm_kernel())

  yi      <- Bem2011[["d"]][1:3]
  sei     <- Bem2011[["se"]][1:3]
  spec    <- .test_step_spec(yi, sei)
  mean    <- matrix(c(.06, .12, .18, -.03, .04, .10), nrow = 2, byrow = TRUE)
  sd      <- matrix(c(.11, .18, .25, .14, .20, .30), nrow = 2, byrow = TRUE)
  omega   <- matrix(c(1, .5, .5, .2, 1, .8, .8, .3), nrow = 2, byrow = TRUE)
  sentinel_spec <- spec

  sentinel_spec[["z_lower"]] <- .selection_jags_bounds(sentinel_spec[["z_lower"]])
  sentinel_spec[["z_upper"]] <- .selection_jags_bounds(sentinel_spec[["z_upper"]])

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = sentinel_spec
    ),
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    ),
    tolerance = 1e-12
  )

  collapsed_spec              <- spec
  collapsed_spec[["n_bins"]]  <- 3L
  collapsed_spec[["z_lower"]] <- c(spec[["z_lower"]][1L], spec[["z_lower"]][3L], spec[["z_lower"]][4L])
  collapsed_spec[["z_upper"]] <- c(spec[["z_upper"]][1L], spec[["z_upper"]][2L], spec[["z_upper"]][4L])
  collapsed_omega             <- omega[, c(1L, 2L, 4L), drop = FALSE]

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = collapsed_omega,
      selection_spec = collapsed_spec
    ),
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    ),
    tolerance = 1e-12
  )
})

test_that("step normalizer sign handling is symmetric for negative effects", {

  skip_if_not(.has_native_selnorm_kernel())

  yi_pos <- Bem2011[["d"]][1:4]
  yi_neg <- -yi_pos
  sei    <- Bem2011[["se"]][1:4]
  mean   <- matrix(c(.04, .10, .16, .22, .02, .08, .14, .20), nrow = 2, byrow = TRUE)
  sd     <- matrix(c(.10, .14, .18, .22, .12, .16, .20, .24), nrow = 2, byrow = TRUE)
  omega  <- matrix(c(1, .7, .4, .2, 1, .9, .5, .3), nrow = 2, byrow = TRUE)

  spec_pos <- .test_step_spec(yi_pos, sei, effect_direction = "positive")
  spec_neg <- .test_step_spec(yi_neg, sei, effect_direction = "negative")

  expect_equal(
    .selnorm_kernel_log_norm_matrix(
      mean           = -mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec_neg
    ),
    .selnorm_kernel_log_norm_matrix(
      mean           = mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec_pos
    ),
    tolerance = 1e-12
  )
})

test_that("trusted step CDF, moments, and RNG match exact fallback paths", {

  skip_if_not(.has_native_selnorm_kernel())

  yi    <- c(.04, .18, .31, .47)
  sei   <- c(.08, .12, .17, .24)
  mean  <- matrix(c(
    .02, .10, .18, .26,
    .05, .08, .20, .33,
    .01, .14, .22, .29
  ), nrow = 3, byrow = TRUE)
  sd    <- matrix(c(
    .11, .19, .25, .31,
    .13, .22, .28, .36,
    .16, .20, .30, .42
  ), nrow = 3, byrow = TRUE)
  omega <- matrix(c(
    1, .35, 1.60, .20,
    .80, .80, .25, 1.10,
    1.40, .10, .95, .30
  ), nrow = 3, byrow = TRUE)

  for (effect_direction in c("positive", "negative")) {
    sign <- if (effect_direction == "positive") 1 else -1
    spec <- .test_step_spec(sign * yi, sei, effect_direction = effect_direction)
    spec_fallback <- spec
    spec_fallback[["telescope_probabilities"]] <- FALSE

    cdf_lower <- .selnorm_kernel_cdf_matrix(
      q              = sign * yi,
      mean           = sign * mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )
    cdf_upper <- .selnorm_kernel_cdf_matrix(
      q              = sign * yi,
      mean           = sign * mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec,
      lower.tail     = FALSE
    )
    moments <- .selnorm_kernel_moments_matrix(
      mean           = sign * mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )

    expect_equal(
      cdf_lower,
      .test_step_cdf_reference(sign * yi, sign * mean, sd, sei, omega, spec),
      tolerance = 1e-12
    )
    expect_equal(
      cdf_upper,
      .test_step_cdf_reference(
        sign * yi, sign * mean, sd, sei, omega, spec, lower.tail = FALSE
      ),
      tolerance = 1e-12
    )
    expect_equal(
      moments,
      .test_step_moments_reference(sign * mean, sd, sei, omega, spec),
      tolerance = 1e-12
    )

    set.seed(511)
    rng_fast <- .selnorm_kernel_rng_matrix(
      mean           = sign * mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec
    )
    set.seed(511)
    rng_fallback <- .selnorm_kernel_rng_matrix(
      mean           = sign * mean,
      sd             = sd,
      sei            = sei,
      omega          = omega,
      selection_spec = spec_fallback
    )
    expect_equal(rng_fast, rng_fallback, tolerance = 1e-12)
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
