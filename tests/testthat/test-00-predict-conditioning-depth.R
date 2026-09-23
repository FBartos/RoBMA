context("Prediction conditioning depth")


test_that("standardized estimate effects preserve original scale and row identities", {

  scales <- matrix(c(.3, .5), 2L, 1L)
  samples <- matrix(c(-.4, .8), 2L, dimnames = list(NULL, "theta"))
  expect_equal(unname(.evaluate.brma.estimate_effects(NULL, scales, TRUE, 1L, samples)),
    matrix(c(-.12, .4), 2L))
  set.seed(811)
  expected <- matrix(stats::rnorm(2L), 2L) * scales
  expected_rng <- .Random.seed
  set.seed(811)
  actual <- .evaluate.brma.estimate_effects(NULL, scales, FALSE, 1L)
  expect_identical(actual, expected)
  expect_identical(.Random.seed, expected_rng)
})


test_that("selection partitions specialized estimate and cluster sources independently", {

  dat <- data.frame(yi = c(.8, -.1), vi = c(.04, .09), paper = c("a", "a"))
  theta <- matrix(c(.5, -.75, .2, .9), 2L, byrow = TRUE)
  gamma <- c(-.3, .6)
  samples <- cbind(mu = c(.1, -.2), tau = .5, rho = .64,
    "theta[2]" = theta[, 2L], "gamma[1]" = gamma, "theta[1]" = theta[, 1L])
  fixed <- matrix(samples[, "mu"], 2L, 2L)
  within <- matrix(.3, 2L, 2L)
  between <- matrix(.4, 2L, 2L)
  for (estimate in c("integrate", "condition")) {
    for (other in c("integrate", "condition")) {
      object <- bselmodel(
        yi = yi, vi = vi, cluster = paper, data = dat, measure = "GEN",
        prior_unit_information_sd = 1, only_priors = TRUE,
        prior_bias = BayesTools::prior_weightfunction(
          "one-sided", .05, BayesTools::wf_fixed(c(1, .4)),
          model = BayesTools::selection_model(estimate_random_effects = estimate,
            other_random_effects = other, known_sampling_variance = "integrate", group = "paper")
        )
      )
      set.seed(812)
      before <- .Random.seed
      parts <- .predict_joint_selection_gaussian_parts(object, object$data, samples,
        fixed, within, between, fitted_context = TRUE)
      expect_identical(.Random.seed, before)
      expected_estimate <- if (estimate == "condition") theta * within else matrix(0, 2L, 2L)
      expected_other <- if (other == "condition") matrix(gamma * .4, 2L, 2L) else matrix(0, 2L, 2L)
      expected_mean <- fixed + expected_estimate + expected_other
      integrated <- (if (estimate == "integrate") diag(.3^2, 2L) else matrix(0, 2L, 2L)) +
        (if (other == "integrate") matrix(.4^2, 2L, 2L) else matrix(0, 2L, 2L))
      retained <- diag(.3^2, 2L) + matrix(.4^2, 2L, 2L) - integrated
      expect_equal(unname(parts$estimate_mean), expected_estimate, tolerance = 0)
      expect_equal(unname(parts$latent_means), unname(expected_mean), tolerance = 1e-15)
      expected_posterior <- expected_mean
      for (s in 1:2) {
        expect_equal(.block_covariance_dense(parts$random_covariance, s),
                     integrated, tolerance = 1e-15)
        expected_posterior[s, ] <- expected_mean[s, ] + integrated %*%
          solve(integrated + diag(dat$vi), dat$yi - expected_mean[s, ])
      }
      expect_equal(unname(.predict_joint_selection_source_posterior(parts, dat$yi, FALSE)),
        unname(expected_posterior), tolerance = 1e-12)
      prediction_context <- .predict_brma_context(object, newdata = NULL, V_new = NULL,
        type = "cluster", conditioning_depth = "cluster", conditioning_depth_specified = TRUE,
        as_measure = TRUE, output_measure = NULL, transform = NULL, probs = c(.025, .975),
        bias_adjusted = FALSE, quiet = TRUE, conditional = FALSE,
        dots = list(.posterior_samples = samples))
      prediction_context$type <- "location"
      location <- .predict_brma_location_state(prediction_context,
        list(within = within, between = between))
      expected_cluster <- fixed + expected_other
      if (other == "integrate") {
        for (s in 1:2) expected_cluster[s, ] <- fixed[s, ] + matrix(.4^2, 2L, 2L) %*%
          solve(integrated + diag(dat$vi), dat$yi - expected_mean[s, ])
      }
      expect_equal(unname(location$mu), unname(expected_cluster), tolerance = 1e-12)
      mixture <- .predict_joint_selection_gaussian_parts(object, object$data, samples,
        fixed, within, between, draw_context = FALSE)
      expect_identical(.Random.seed, before)
      expect_equal(mixture$means, fixed, tolerance = 0)
      for (s in 1:2) {
        expect_equal(.block_covariance_dense(mixture$context_covariance, s),
                     retained, tolerance = 1e-15)
      }
    }
  }
  # A new estimate within a fitted cluster draws fresh estimate context, never fitted theta.
  cluster_mean <- fixed + matrix(gamma * .4, 2L, 2L)
  context <- list(object = object, new_data = object$data, posterior_samples = samples,
    K = 2L, conditioning_depth = "cluster", outcome_data = object$data$outcome,
    same_data = TRUE, known_V_new = NULL)
  set.seed(813)
  fresh_estimate <- matrix(stats::rnorm(4L), 2L, 2L) * within
  set.seed(813)
  response <- .predict_joint_selection_response_setup(context,
    list(mu = cluster_mean, fixed_mu = fixed), list(within = within, between = between))
  expect_equal(unname(response$latent_means), unname(cluster_mean + fresh_estimate), tolerance = 1e-15)
  # The normal outcome kernel uses the fitted SE; its squared value can differ
  # from the original vi literal by a rounding bit after sqrt(vi).
  for (s in 1:2) expect_equal(.block_covariance_dense(response$covariance, s),
    diag(object$data$outcome$sei^2), tolerance = 0)
})

test_that("marginal prediction is independent of implicit versus explicit design", {

  dat <- data.frame(
    yi = c(0.2, 0.5),
    vi = c(0.04, 0.09)
  )
  object <- brma(
    yi                        = yi,
    vi                        = vi,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(
      0.10, 0.30,
      0.20, 0.40,
      0.15, 0.35
    ),
    nrow     = 3L,
    byrow    = TRUE,
    dimnames = list(NULL, c("mu", "tau"))
  )

  set.seed(804)
  implicit <- predict(
    object,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  set.seed(804)
  explicit <- predict(
    object,
    newdata            = dat,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(implicit, explicit, tolerance = 0)
  expect_error(
    predict(
      object,
      newdata            = dat,
      type               = "estimate",
      conditioning_depth = "estimate",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "fitted observation identities"
  )
  expect_error(
    predict(
      object,
      type               = "terms",
      conditioning_depth = "estimate",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "only available for type = 'estimate' or type = 'response'",
    fixed = TRUE
  )
})


test_that("estimate depth includes conditional latent uncertainty", {

  S   <- 30000L
  mu  <- 0.2
  tau <- 0.4
  yi  <- 0.8
  sei <- 0.3
  tau2 <- tau^2
  vi   <- sei^2
  expected_mean <- mu + tau2 / (tau2 + vi) * (yi - mu)
  expected_var  <- tau2 * vi / (tau2 + vi)

  set.seed(805)
  observed <- .evaluate.brma.true_effects_posterior.norm(
    mu_samples = matrix(mu, nrow = S, ncol = 1L),
    tau_within = matrix(tau, nrow = S, ncol = 1L),
    yi         = yi,
    sei        = sei
  )

  expect_equal(mean(observed), expected_mean, tolerance = 0.004)
  expect_equal(as.numeric(stats::var(observed)), expected_var,
               tolerance = 0.002)
})


test_that("known-V estimate depth preserves the joint conditional posterior", {

  S   <- 30000L
  mu  <- c(0.1, -0.2)
  tau <- c(0.4, 0.3)
  yi  <- c(0.5, 0.1)
  V   <- matrix(c(0.09, 0.03, 0.03, 0.16), nrow = 2L)
  D   <- diag(tau^2)
  conditional_mean <- mu + D %*% solve(D + V, yi - mu)
  conditional_v    <- D - D %*% solve(D + V, D)
  known_V <- .known_v_newdata_prepare(V, k = 2L)

  set.seed(806)
  observed <- .evaluate.brma.known_v_posterior.norm(
    mu_samples = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
    tau_within = matrix(tau, nrow = S, ncol = 2L, byrow = TRUE),
    yi         = yi,
    known_V    = known_V
  )

  expect_equal(colMeans(observed), as.vector(conditional_mean), tolerance = 0.004)
  expect_equal(stats::cov(observed), conditional_v, tolerance = 0.003)
})


test_that("specialized multilevel estimate depth preserves joint latent uncertainty", {

  S       <- 15000L
  mu      <- c(0.1, -0.2)
  tau_b   <- 0.3
  tau_w   <- c(0.4, 0.25)
  yi      <- c(0.5, 0.1)
  vi      <- c(0.09, 0.16)
  latent_v <- tau_b^2 * matrix(1, nrow = 2L, ncol = 2L) +
    diag(tau_w^2)
  marginal_v     <- latent_v + diag(vi)
  conditional_mean <- mu + latent_v %*% solve(marginal_v, yi - mu)
  conditional_v    <- latent_v -
    latent_v %*% solve(marginal_v, latent_v)

  set.seed(807)
  observed <- matrix(mu, nrow = S, ncol = 2L, byrow = TRUE) +
    .evaluate.brma.multilevel_posterior.norm(
      mu_samples = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
      tau_within = matrix(tau_w, nrow = S, ncol = 2L, byrow = TRUE),
      tau_between = matrix(tau_b, nrow = S, ncol = 2L),
      yi          = yi,
      vi          = vi,
      cluster     = c(1L, 1L)
    )

  expect_equal(colMeans(observed), as.vector(conditional_mean), tolerance = 0.006)
  expect_equal(stats::cov(observed), conditional_v, tolerance = 0.005)
})


test_that("specialized multilevel cluster depth retains posterior uncertainty", {

  S       <- 20000L
  mu      <- c(0.1, -0.2)
  tau_b   <- 0.3
  tau_w   <- c(0.4, 0.25)
  yi      <- c(0.5, 0.1)
  vi      <- c(0.09, 0.16)
  cluster_v  <- tau_b^2 * matrix(1, nrow = 2L, ncol = 2L)
  marginal_v <- cluster_v + diag(tau_w^2 + vi)
  expected_mean <- cluster_v %*% solve(marginal_v, yi - mu)
  expected_v    <- cluster_v -
    cluster_v %*% solve(marginal_v, cluster_v)

  set.seed(808)
  observed <- .evaluate.brma.multilevel_posterior.norm(
    mu_samples  = matrix(mu, nrow = S, ncol = 2L, byrow = TRUE),
    tau_within  = matrix(tau_w, nrow = S, ncol = 2L, byrow = TRUE),
    tau_between = matrix(tau_b, nrow = S, ncol = 2L),
    yi          = yi,
    vi          = vi,
    cluster     = c(1L, 1L),
    component   = "cluster"
  )

  expect_equal(colMeans(observed), as.vector(expected_mean), tolerance = 0.004)
  expect_equal(stats::cov(observed), expected_v, tolerance = 0.002)
  expect_gt(stats::var(observed[, 1L]), 0)
})


test_that("conditional random selection integrates the complete sampling covariance", {

  dat <- data.frame(
    yi     = c(0.8, -0.1),
    vi     = c(0.04, 0.09),
    study  = c("s1", "s1"),
    effect = c("e1", "e2")
  )
  sampled_study <- matrix(c(0.4, -0.2), nrow = 2L, ncol = 2L)
  testthat::local_mocked_bindings(
    .evaluate.brma.random_effects = function(...) sampled_study,
    .package = "RoBMA"
  )

  for (with_sampling_factor in c(FALSE, TRUE)) {
    loading <- c(0.1, 0.2)
    covariance_args <- if (with_sampling_factor) {
      list(V = known_v_factor(dat$vi, matrix(loading, ncol = 1L)))
    } else {
      list(vi = dat$vi)
    }
    object <- do.call(bselmodel.mv, c(
      list(
        yi                        = dat$yi,
        random                    = list(study = ~ 1 | study, effect = ~ 1 | effect),
        data                      = dat,
        measure                   = "GEN",
        prior_unit_information_sd = 1,
        prior_bias = BayesTools::prior_weightfunction(
          "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
          model = BayesTools::selection_model(group = "study", known_sampling_variance = "integrate")
        ),
        only_priors               = TRUE
      ),
      covariance_args
    ))
    random_terms <- .fitted_formula_design(object, "mu", required = TRUE)[["random_effects"]]
    sd_names <- vapply(random_terms, function(term) term[["sd_parameter_names"]][[1L]], character(1))
    names(sd_names) <- vapply(random_terms, `[[`, character(1), "block_name")
    posterior_samples <- matrix(
      c(0.1, 0.5, 0.4, 0.3, 0.2, 0.8, 0.3, -0.4),
      nrow     = 2L,
      byrow    = TRUE,
      dimnames = list(NULL, c("mu", sd_names[["study"]], sd_names[["effect"]], "sampling_z[1]"))
    )

    # Conditional normal identity: shrink the integrated local effect using
    # the full V, irrespective of its computational factor representation.
    conditional_location <- posterior_samples[, "mu"] + sampled_study
    within_variance <- matrix(c(0.4, 0.3)^2, nrow = 2L, ncol = 2L)
    observed_y <- matrix(dat$yi, nrow = 2L, ncol = 2L, byrow = TRUE)
    V <- diag(dat$vi) + if (with_sampling_factor) tcrossprod(loading) else 0
    expected <- conditional_location
    shrinkage <- vector("list", 2L)
    for (draw in 1:2) {
      Q <- diag(within_variance[draw, ])
      shrinkage[[draw]] <- Q %*% solve(Q + V)
      expected[draw, ] <- expected[draw, ] +
        as.vector(shrinkage[[draw]] %*% (observed_y[draw, ] - conditional_location[draw, ]))
    }

    set.seed(809)
    random_state <- .Random.seed
    actual <- blup(object, .posterior_samples = posterior_samples)
    components <- ranef(object, expand = TRUE, simplify = FALSE,
                        .posterior_samples = posterior_samples)
    expect_identical(.Random.seed, random_state)
    expect_equal(unname(as.matrix(actual)), expected, tolerance = 1e-12)
    expect_equal(unname(as.matrix(components[["study"]])), sampled_study,
                 tolerance = 1e-12)
    expect_equal(unname(as.matrix(components[["effect"]])),
                 expected - conditional_location, tolerance = 1e-12)

    set.seed(810)
    random_noise   <- matrix(stats::rnorm(4L), nrow = 2L, byrow = TRUE) * sqrt(within_variance)
    sampling_noise <- matrix(0, 2L, 2L)
    for (draw in 1:2) sampling_noise[draw, ] <- stats::rnorm(2L) %*% chol(V)
    expected_random_state <- .Random.seed
    expected_draws <- expected + random_noise
    for (draw in 1:2) {
      Q <- diag(within_variance[draw, ])
      H <- shrinkage[[draw]]
      expect_equal((diag(2L) - H) %*% Q %*% t(diag(2L) - H) + H %*% V %*% t(H),
        Q - Q %*% solve(Q + V) %*% Q, tolerance = 1e-14)
      expected_draws[draw, ] <- expected_draws[draw, ] -
        as.vector(H %*% (random_noise[draw, ] + sampling_noise[draw, ]))
    }
    set.seed(810)
    observed_draws <- predict(
      object,
      type               = "estimate",
      conditioning_depth = "estimate",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    )
    expect_equal(unname(as.matrix(observed_draws)), expected_draws,
                 tolerance = 1e-12)
    expect_identical(.Random.seed, expected_random_state)
  }

})


test_that("one integrated random source retains correlated and zero-variance updates", {

  dat <- data.frame(yi = c(.8, -.1, .3), study = c("a", "a", "b"))
  V <- matrix(c(.10, .02, .01, .02, .20, .03, .01, .03, .15), 3L)
  object <- bselmodel.mv(
    yi = yi, V = V, random = list(study = ~ 1 | study), data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    selection = selection_model(other_random_effects = "integrate")
  )
  term <- .fitted_formula_design(object, "mu", required = TRUE)[["random_effects"]][[1L]]
  samples <- cbind(mu = c(.1, -.2), c(.4, 0))
  colnames(samples)[[2L]] <- term[["sd_parameter_names"]][[1L]]
  expected <- matrix(0, 2L, 3L)
  Q <- .4^2 * outer(dat$study, dat$study, `==`)
  expected[1L, ] <- Q %*% solve(V + Q, dat$yi - samples[1L, "mu"])

  withr::local_seed(813)
  initial_seed <- .Random.seed
  actual <- ranef(object, expand = TRUE, simplify = FALSE,
                  .posterior_samples = samples)
  expect_identical(.Random.seed, initial_seed)
  expect_equal(unname(as.matrix(actual[["study"]])), expected, tolerance = 1e-12)

  # Return the source update directly; subtracting two large fitted locations
  # would erase this small, identifiable contribution.
  shifted <- object
  shifted[["data"]][["outcome"]][["yi"]] <- 1e16 + c(2, -2, 4)
  shifted_samples <- samples[1L, , drop = FALSE]
  shifted_samples[1L, ] <- c(1e16, .04)
  Q <- .04^2 * outer(dat$study, dat$study, `==`)
  expected_shifted <- Q %*% solve(V + Q, c(2, -2, 4))
  actual_shifted <- ranef(shifted, expand = TRUE, simplify = FALSE,
                          .posterior_samples = shifted_samples)
  expect_equal(as.numeric(as.matrix(actual_shifted[["study"]])),
               as.vector(expected_shifted), tolerance = 1e-12)
})
