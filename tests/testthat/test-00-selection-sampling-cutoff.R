test_that("sampling-conditioned deletion keeps an observed cutoff only while that outcome is fixed", {

  p_cut <- .49
  cut   <- stats::qnorm(p_cut, lower.tail = FALSE)
  mu    <- 2
  # Adding the fixed sampling residual back to the mean misses the cutoff.
  expect_lt(mu + (cut - mu), cut)
  V     <- matrix(c(1, .5, .5, 1), 2L)
  total <- V + diag(c(0, .25))
  y     <- c(cut, -.0625)
  e2    <- drop(V %*% solve(total, y - mu))[2L]
  mean1 <- mu + .5 * e2
  for (sign in c(1, -1)) for (rule in c("product", "best")) {
    dat <- data.frame(yi = sign * y, x = c(0, 1), study = "a")
    model <- bselmodel.mv(
      yi = yi, V = V, random = ~ 0 + x | study, data = dat, measure = "GEN",
      effect_direction = if (sign == 1) "positive" else "negative",
      standardize_continuous_predictors = FALSE,
      prior_effect = prior("point", list(location = sign * mu)),
      prior_heterogeneity = prior("point", list(location = .5)),
      prior_bias = prior_weightfunction("one-sided", p_cut, wf_fixed(c(1, .2)),
        model = selection_model(known_sampling_variance = "condition",
          other_random_effects = "integrate", weight_rule = rule, group = study)),
      only_priors = TRUE, silent = TRUE
    )
    fit <- structure(list(), formula_design = model$formula_design,
      prior_list = c(model$formula_design$mu$prior_list, .create_fit_priors(model$data, model$priors)))
    draw <- c(mu_intercept = sign * mu, mu__xREx__study_x = .5, "mu__xREx__study_xRE_Zx[1,1]" = 0,
      "sampling_z[1]" = 0, "sampling_z[2]" = 0, "omega[1]" = 1, "omega[2]" = .2)
    setup <- .log_lik_posterior_setup(fit, matrix(draw, 1L, dimnames = list(NULL, names(draw))),
      model$data, model$priors, "estimate", NULL)

    block <- .selection_conditioned_sampling_deletion_loglik(setup, list(1L, 2L, 1:2))
    if (rule == "product") {
      # The fixed row's own weight cancels at every candidate sampling error.
      score1 <- stats::dnorm(cut, mean1, sqrt(.75), log = TRUE)
      targets <- .selection_conditioned_sampling_estimate_targets(setup, "log_density")
      expect_equal(targets$log_density[1L, 1L], score1, tolerance = 1e-12)
      expect_equal(block[1L, 1L], score1, tolerance = 1e-12)
      next
    }

    # Signed coordinates below point in the favored direction. With row 1
    # retained at the leading cutoff, every candidate for row 2 has weight
    # one. Conditioning e2 on e1 leaves variance .75, plus .25 for the slope.
    mean2 <- mu + .5 * (cut - mu)
    score2 <- stats::dnorm(y[2], mean2, 1, log = TRUE)
    targets <- .selection_conditioned_sampling_estimate_targets(setup,
      c("log_density", "cdf", "mean", "variance"))
    expect_equal(targets$log_density[1L, 2L], score2, tolerance = 5e-5)
    expect_equal(targets$cdf[1L, 2L], stats::pnorm(y[2], mean2, 1, lower.tail = sign == 1),
                 tolerance = 5e-5)
    expect_equal(targets$mean[1L, 2L], sign * mean2, tolerance = 5e-5)
    expect_equal(targets$variance[1L, 2L], 1, tolerance = 5e-5)

    # When row 1 is deleted, e1 varies. Below the cutoff the observed row
    # 2 is in the low-weight bin, while the new slope can cross the cutoff.
    # At and above the cutoff, the ratio of observed to expected weight is one.
    below <- stats::pnorm(cut, mean1, sqrt(.75))
    mass2 <- .2 + .8 * stats::pnorm(cut, mu + e2, .5, lower.tail = FALSE)
    score1 <- stats::dnorm(cut, mean1, sqrt(.75), log = TRUE) -
      log((1 - below) + below * .2 / mass2)
    expect_equal(targets$log_density[1L, 1L], score1, tolerance = 5e-5)

    # The block integrator exercises both partial and complete deletions.
    # In a complete deletion, conditioning on y fixes the first error
    # exactly. Its leading bin makes the joint density the ordinary Gaussian.
    expect_equal(block[1L, 1L], score1, tolerance = .001)
    expect_equal(block[1L, 2L], score2, tolerance = 1e-12)
    expect_equal(block[1L, 3L], mvtnorm::dmvnorm(y, rep(mu, 2L), total, log = TRUE),
                 tolerance = 1e-12)
  }
})
