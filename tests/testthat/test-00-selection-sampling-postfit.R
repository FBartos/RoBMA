test_that("sampling endpoints retain the whole covariance independently of its representation", {

  variances <- c(.4, .9, 1.2)
  loading <- cbind(c(.2, .3, 0), c(0, -.2, .4))
  V <- diag(variances) + tcrossprod(loading)
  inputs <- list(NULL, diag(V), V, known_v_factor(variances, loading))
  within <- matrix(c(.2, .3, .4), 2L, 3L, byrow = TRUE)
  samples <- matrix(0, 2L, 0L)
  for (input in inputs) for (sampling_mode in c("condition", "integrate")) {
    data <- list(outcome = data.frame(yi = c(.2, -.1, .4), sei = sqrt(diag(V))))
    if (!is.null(input)) {
      attr(data, "known_V") <- TRUE
      attr(data, "known_V_data") <- .known_v_prepare(
        V = input, keep_rows = rep(TRUE, 3L),
        known_v_parameterization = "block_mvn", warn_singular = FALSE
      )
    }
    attr(data, "selection_model") <- structure(list(
      schema_version = 3L, estimate_random_effects = "integrate",
      other_random_effects = "condition", known_sampling_variance = sampling_mode,
      sources = list(random = list(list(name = "estimate", role = "estimate", retained = FALSE))),
      applicability = list(estimate_random_effects = TRUE, other_random_effects = FALSE,
                           known_sampling_variance = TRUE),
      groups = list(row_index = 1:3)
    ), class = c("RoBMA_selection_model", "list"))
    parts <- .predict_joint_selection_gaussian_parts(
      object = list(data = data), data = data, posterior_samples = samples,
      fixed_mu = matrix(.1, 2L, 3L), within = within, between = matrix(0, 2L, 3L),
      draw_context = FALSE
    )
    expected_V <- if (is.null(input) || is.null(dim(input)) && is.numeric(input)) diag(diag(V)) else V
    random <- diag(within[1L, ]^2)
    expect_equal(.block_covariance_dense(parts$covariance, 1L),
      random + if (sampling_mode == "integrate") expected_V else 0, tolerance = 1e-14)
    expect_equal(.block_covariance_dense(parts$context_covariance, 1L),
      if (sampling_mode == "condition") expected_V else matrix(0, 3L, 3L), tolerance = 1e-14)
    expect_identical(.selection_postfit_target_metadata(data)$sampling_structure$source,
      "whole_sampling_error")
  }
})


test_that("fitted source reconstruction assigns each covariance correction to its owner", {

  data <- list(outcome = data.frame(yi = c(.7, -.2), sei = c(.5, .8), cluster = c(1L, 1L)))
  attr(data, "cluster") <- TRUE
  attr(data, "effect_direction") <- "positive"
  samples <- cbind("theta[1]" = c(.3, -.2), "theta[2]" = c(-.5, .4),
                   "gamma[1]" = c(.2, -.4))
  setup <- list(data = data, fit = NULL, priors = list(), posterior_samples = samples,
    S = 2L, K = 2L, tau_within = matrix(c(.3, .4), 2L, 2L, byrow = TRUE),
    tau_between = matrix(.2, 2L, 2L))
  V <- matrix(c(.25, .07, .07, .64), 2L)
  e0 <- rbind(c(.1, -.3), c(-.2, .15))
  for (retain_estimate in c(FALSE, TRUE)) for (retain_cluster in c(FALSE, TRUE)) {
    sources <- list(list(name = "estimate", role = "estimate", retained = retain_estimate),
                    list(name = "cluster", role = "other", retained = retain_cluster))
    attr(setup$data, "selection_model") <- structure(list(
      schema_version = 3L, sources = list(random = sources)
    ), class = c("RoBMA_selection_model", "list"))
    estimate <- samples[, 1:2, drop = FALSE] * setup$tau_within
    cluster <- matrix(samples[, 3L] * .2, 2L, 2L)
    Q_estimate <- diag(c(.3, .4)^2)
    Q_cluster <- matrix(.2^2, 2L, 2L)
    # Conditioning Gaussian auxiliaries on the observed vector reconstructs
    # every source. Selection roles alter normalization, not this identity.
    Q <- Q_estimate + Q_cluster
    delta <- t(solve(V + Q, t(matrix(data$outcome$yi, 2L, 2L, byrow = TRUE) -
      .1 - e0 - estimate - cluster)))
    expected_estimate <- estimate + delta %*% Q_estimate
    expected_cluster <- cluster + delta %*% Q_cluster
    reconstructed <- .selection_random_source_posterior(setup, list(correction = delta))
    expect_equal(unname(reconstructed$estimate), unname(expected_estimate), tolerance = 1e-13)
    expect_equal(unname(reconstructed$cluster), unname(expected_cluster), tolerance = 1e-13)
    expect_equal(unname(.1 + reconstructed$estimate + reconstructed$cluster + e0 + delta %*% V),
      matrix(data$outcome$yi, 2L, 2L, byrow = TRUE), tolerance = 1e-13)
    integrated_Q <- (if (retain_estimate) 0 else Q_estimate) +
      (if (retain_cluster) 0 else Q_cluster)
    baseline <- matrix(.1, 2L, 2L) +
      (if (retain_estimate) reconstructed$estimate else 0) +
      (if (retain_cluster) reconstructed$cluster else 0)
    state <- list(e = e0 + delta %*% V, baseline_mu = baseline,
      integrated_covariance = array(rep(integrated_Q, each = 2L), c(2L, 2L, 2L)))
    means <- .selection_random_source_conditional_means(setup, state, reconstructed)
    if (!retain_estimate && !retain_cluster) {
      total <- reconstructed$estimate + reconstructed$cluster
      expect_equal(unname(means$estimate), unname(total %*% solve(Q) %*% Q_estimate), tolerance = 1e-13)
      expect_equal(unname(means$cluster), unname(total %*% solve(Q) %*% Q_cluster), tolerance = 1e-13)
    } else {
      expect_identical(means, reconstructed)
    }
  }
})


test_that("a retained complete sampling error identifies total truth with singular random variation", {

  Q <- matrix(c(1, 2, 2, 4), 2L)
  parts <- list(
    means = matrix(c(.5, -.1), 1L), latent_means = matrix(c(.2, .3), 1L),
    random_covariance = .as_block_covariance(array(Q, c(1L, 2L, 2L))),
    sampling_covariance = .block_covariance_zero(1L, 2L),
    covariance = .as_block_covariance(array(Q, c(1L, 2L, 2L)))
  )
  y <- c(.7, .3)
  withr::local_seed(291)
  before <- .Random.seed
  for (draw in c(FALSE, TRUE)) {
    expect_equal(as.vector(.predict_joint_selection_source_posterior(parts, y, draw)), c(.4, .7))
  }
  expect_identical(.Random.seed, before)
})


test_that("ranef allocation means retain a singular integrated Gaussian law", {

  data <- list(outcome = data.frame(yi = rep(.4, 3L)))
  sources <- list(list(name = "first", retained = FALSE), list(name = "second", retained = FALSE))
  attr(data, "selection_model") <- structure(list(schema_version = 3L,
    sources = list(random = sources)), class = c("RoBMA_selection_model", "list"))
  setup <- list(data = data, S = 1L, K = 3L)
  state <- list(e = matrix(0, 1L, 3L), baseline_mu = matrix(0, 1L, 3L),
    integrated_covariance = array(5, c(1L, 3L, 3L)))
  testthat::local_mocked_bindings(
    .selection_random_source_covariance = function(setup, source) {
      array(if (source$name == "first") 1 else 4, c(1L, 3L, 3L))
    }, .package = "RoBMA"
  )
  # Independent scalar identity: two shared effects with variances 1 and 4
  # receive one fifth and four fifths of their observed sum.
  for (first in c(-1, .2)) {
    contributions <- list(first = matrix(first, 1L, 3L), second = matrix(.4 - first, 1L, 3L))
    means <- .selection_random_source_conditional_means(setup, state, contributions)
    expect_equal(means$first, matrix(.08, 1L, 3L), tolerance = 1e-14)
    expect_equal(means$second, matrix(.32, 1L, 3L), tolerance = 1e-14)
  }
})


test_that("singleton JAGS sampling and cluster auxiliaries retain their row identities", {

  skip_if_not_installed("rjags")
  connection <- textConnection(paste0(
    "model { for (i in 1:1) { ",
    "sampling_z[i] ~ dnorm(0, 1); gamma[i] ~ dnorm(0, 1) } }"
  ))
  on.exit(close(connection), add = TRUE)
  model <- rjags::jags.model(
    connection, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
    inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 731L)
  )
  samples <- as.matrix(rjags::coda.samples(
    model, c("sampling_z", "gamma"), n.iter = 3L, progress.bar = "none"
  ))
  expect_true(all(c("sampling_z", "gamma") %in% colnames(samples)))
  data <- list(outcome = data.frame(yi = c(.2, -.1), sei = c(.5, 1)))
  attr(data, "known_V") <- TRUE
  attr(data, "effect_direction") <- "positive"
  attr(data, "known_V_data") <- .known_v_resolve_selection_structure(
    .known_v_prepare(tcrossprod(c(.5, 1)), c(TRUE, TRUE), "block_mvn",
                     warn_singular = FALSE)
  )
  sampling <- .selection_sampling_auxiliary(NULL, data, samples)
  expect_equal(unname(sampling), outer(samples[, "sampling_z"], c(.5, 1)),
               tolerance = 1e-14)
  cluster <- .evaluate.brma.cluster_effects(
    fit = NULL, tau_between = matrix(.4, 3L, 2L), cluster = c(1L, 1L),
    same_data = TRUE, effect_direction = "positive", posterior_samples = samples
  )
  expect_equal(unname(cluster), matrix(samples[, "gamma"] * .4, 3L, 2L),
               tolerance = 0)
  indexed <- samples
  colnames(indexed) <- paste0(colnames(indexed), "[1]")
  expect_equal(.selection_sampling_auxiliary(NULL, data, indexed), sampling,
               tolerance = 0)
  for (parameter in c("sampling_z", "gamma")) {
    wrong_index <- matrix(samples[, parameter], ncol = 1L,
                          dimnames = list(NULL, paste0(parameter, "[2]")))
    expect_error(.extract_indexed_parameter_samples(
      wrong_index, parameter, n_expected = 1L
    ), paste0("Missing posterior column(s): ", parameter, "[1]."), fixed = TRUE)
    expect_error(.extract_indexed_parameter_samples(
      samples, parameter, n_expected = 2L
    ), paste0("Expected 2 posterior ", parameter, " column(s), found 1."),
    fixed = TRUE)
  }
})
