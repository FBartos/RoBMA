.pr82_sampling_memory_setup <- function(S = 4L, K = 6L) {

  object <- bselmodel(
    yi = seq(-.3, .5, length.out = K), sei = seq(.7, 1, length.out = K),
    cluster = rep(seq_len(K / 2L), each = 2L), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, .4)), model = BayesTools::selection_model(
        other_random_effects = "condition", known_sampling_variance = "condition")),
    only_priors = TRUE
  )
  samples <- cbind(mu = seq(-.1, .2, length.out = S),
                   tau = seq(.2, .5, length.out = S), rho = .3)
  add <- function(stem, n) {
    out <- matrix(seq(-.4, .4, length.out = S * n), S, n)
    colnames(out) <- paste0(stem, "[", seq_len(n), "]")
    out
  }
  samples <- cbind(samples, add("theta", K), add("gamma", K / 2L), add("sampling_z", K))
  .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
}

test_that("conditioned sampling blocks preserve dense Gaussian source reconstruction", {

  setup <- .pr82_sampling_memory_setup()
  withr::local_seed(729L)
  seed <- .Random.seed
  state <- .selection_conditioned_sampling_state(setup)
  expect_identical(.Random.seed, seed)
  expect_s3_class(state[["integrated_covariance"]], "RoBMA_block_covariance")
  expect_s3_class(state[["total_covariance"]], "RoBMA_block_covariance")
  S <- setup[["S"]]
  K <- setup[["K"]]
  clusters <- setup[["data"]][["outcome"]][["cluster"]]
  samples <- setup[["posterior_samples"]]
  V <- diag(setup[["data"]][["outcome"]][["sei"]]^2)
  y <- setup[["data"]][["outcome"]][["yi"]]
  sources <- .selection_random_source_posterior(setup, state)
  source_names <- .data_selection_model(setup[["data"]])[["sources"]][["random"]]
  estimate_name <- Filter(function(source) source[["role"]] == "estimate", source_names)[[1L]][["name"]]
  other_name <- Filter(function(source) source[["role"]] == "other", source_names)[[1L]][["name"]]
  for (draw in seq_len(S)) {
    within <- setup[["tau_within"]][draw, ]
    between <- setup[["tau_between"]][draw, ]
    estimate <- samples[draw, paste0("theta[", seq_len(K), "]")] * within
    other <- samples[draw, paste0("gamma[", match(clusters, unique(clusters)), "]")] * between
    e0 <- samples[draw, paste0("sampling_z[", seq_len(K), "]")] * sqrt(diag(V))
    Q <- diag(within^2)
    R <- outer(between, between) * outer(clusters, clusters, `==`)
    total <- V + R + Q
    delta <- as.vector(solve(total, y - setup[["mu"]][draw, ] - e0 - other - estimate))
    expected_e <- as.vector(e0 + V %*% delta)
    expected_other <- as.vector(other + R %*% delta)
    expected_estimate <- as.vector(estimate + Q %*% delta)
    baseline <- setup[["mu"]][draw, ] + expected_other
    expect_equal(state[["e"]][draw, ], expected_e, tolerance = 1e-13)
    expect_equal(state[["baseline_mu"]][draw, ], baseline, tolerance = 1e-13)
    expect_equal(as.numeric(sources[[other_name]][draw, ]), expected_other, tolerance = 1e-13)
    expect_equal(as.numeric(sources[[estimate_name]][draw, ]), expected_estimate, tolerance = 1e-13)
    expect_identical(.selection_covariance_draw(state[["integrated_covariance"]], draw), Q)
    expect_identical(.selection_covariance_draw(state[["total_covariance"]], draw), total)
    root <- chol(total)
    z <- forwardsolve(t(root), y - setup[["mu"]][draw, ])
    log_gaussian <- -.5 * K * log(2 * pi) - sum(log(diag(root))) - sum(z^2) / 2
    log_mass <- sum(log(.4 + .6 * pnorm((baseline + expected_e) / within)))
    expected <- log_gaussian + sum(log(ifelse(y >= 0, 1, .4))) - log_mass
    expect_equal(state[["log_lik"]][draw], expected, tolerance = 1e-12)
  }
})

test_that("conditioned source storage follows retained rather than candidate dependencies", {

  setup <- .pr82_sampling_memory_setup(S = 12L, K = 10L)
  state <- .selection_conditioned_sampling_state(setup)
  sources <- .data_selection_model(setup[["data"]])[["sources"]][["random"]]
  retained <- Filter(function(source) source[["retained"]], sources)[[1L]]
  covariance <- .selection_random_source_covariance(setup, retained)
  expect_equal(length(covariance[["blocks"]]), 5L)
  expect_true(all(lengths(covariance[["blocks"]]) == 2L))
  dense <- .selection_covariance_draw(covariance, 1L)
  expect_gt(dense[1L, 2L], 0)
  expect_equal(dense[1L, 3L], 0)
  integrated <- .selection_covariance_draw(state[["integrated_covariance"]], 1L)
  expect_equal(integrated[1L, 2L], 0)
  stored <- sum(vapply(state[["total_covariance"]][["values"]], length, integer(1L)))
  expect_equal(stored, setup[["S"]] * 5L * 2L^2)
})

test_that("formula source covariance chunks retain declared cross-row dependencies", {

  dat <- data.frame(study = factor(rep(1:2, each = 2L)))
  object <- bselmodel.mv(yi = c(.1, -.2, .3, .4), V = diag(4), random = ~ 1 | study,
    data = dat, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1,
    prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior("point", list(.3))),
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, .4)), model = BayesTools::selection_model(
        known_sampling_variance = "condition")), only_priors = TRUE)
  fit <- structure(list(), formula_design = object[["formula_design"]],
    prior_list = c(object[["formula_design"]][["mu"]][["prior_list"]],
                   .create_fit_priors(object[["data"]], object[["priors"]])))
  samples <- matrix(seq(-.1, .1, length.out = 5L), 5L, 1L,
                     dimnames = list(NULL, "mu_intercept"))
  setup <- list(object = list(fit = fit, data = object[["data"]], priors = object[["priors"]]),
    fit = fit, data = object[["data"]], priors = object[["priors"]],
    posterior_samples = samples, S = 5L, K = 4L)
  source <- .data_selection_model(object[["data"]])[["sources"]][["random"]][[1L]]
  withr::local_options(RoBMA.known_v_covariance_max_bytes = .known_v_covariance_peak_bytes(1L, 4L))
  actual <- .selection_random_source_covariance(setup, source)
  expected <- matrix(0, 4L, 4L)
  expected[1:2, 1:2] <- .3^2
  expected[3:4, 3:4] <- .3^2
  for (draw in seq_len(5L)) {
    expect_equal(.selection_covariance_draw(actual, draw), expected, tolerance = 0)
  }
})
