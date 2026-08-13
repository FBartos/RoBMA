test_that("shared Gaussian bridge covariance matches exact known-V likelihoods", {

  dat <- data.frame(
    yi = c(0.10, -0.20, 0.15, 0.05),
    study = factor(c("s1", "s1", "s2", "s3"))
  )
  V <- matrix(
    c(
      0.04, 0.01, 0, 0,
      0.01, 0.09, 0, 0,
      0, 0, 0.05, 0.012,
      0, 0, 0.012, 0.06
    ),
    nrow = 4L,
    byrow = TRUE
  )
  random_covariance <- 0.3^2 * outer(dat$study, dat$study, "==")
  bridge_context <- structure(
    list(
      nodes = numeric(),
      marginalized_random = list(
        mu = list(covariance = random_covariance)
      )
    ),
    class = c(
      "BayesTools_bridge_marginal_context",
      "BayesTools_bridge_context",
      "list"
    )
  )

  for (backend in c("block_mvn", "whitened")) {
    object <- brma.mv(
      yi = yi,
      V = V,
      random = ~ 1 | study,
      data = dat,
      measure = "GEN",
      known_v_parameterization = backend,
      marginalize_estimate_level = FALSE,
      prior_unit_information_sd = 1,
      only_priors = TRUE
    )
    fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
    fit_data <- .marglik_add_random_covariance_bridge_data(
      fit_data = fit_data,
      model_data = object[["data"]],
      marginalizing = TRUE,
      dependency_blocks = list(1:2, 3:4)
    )
    actual <- .log_posterior(
      parameters = list(mu = 0),
      data = fit_data,
      is_mods = FALSE,
      is_scale = FALSE,
      is_random = TRUE,
      is_multilevel = FALSE,
      is_weights = FALSE,
      is_known_v = TRUE,
      is_PET = FALSE,
      is_PEESE = FALSE,
      is_weightfunction = FALSE,
      effect_direction = "positive",
      outcome_type = "norm",
      model_data = object[["data"]],
      bridge_context = bridge_context
    )
    expected <- .marglik_mvn_log_density(
      y = dat$yi,
      mean = rep(0, nrow(dat)),
      covariance = V + random_covariance
    )

    expect_equal(actual, expected, tolerance = 1e-12, info = backend)
  }
})

test_that("latent known-V bridge conditions on sampling effects and integrates random effects", {

  dat <- data.frame(
    yi = c(0.10, -0.20, 0.15),
    study = factor(c("s1", "s1", "s2"))
  )
  V <- matrix(
    c(
      0.04, 0.01, 0,
      0.01, 0.09, 0.012,
      0, 0.012, 0.05
    ),
    nrow = 3L,
    byrow = TRUE
  )
  object <- brma.mv(
    yi = yi,
    V = V,
    random = ~ 1 | study,
    data = dat,
    measure = "GEN",
    known_v_parameterization = "latent",
    known_v_residual_fraction = 0.5,
    marginalize_estimate_level = FALSE,
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  known_V <- .data_known_v_data(object[["data"]])
  sampling_z <- seq_len(.known_v_rank(known_V)) / 10
  parameters <- list(mu = 0, sampling_z = sampling_z)
  random_covariance <- 0.25^2 * outer(dat$study, dat$study, "==")
  bridge_context <- structure(
    list(
      nodes = numeric(),
      marginalized_random = list(
        mu = list(covariance = random_covariance)
      )
    ),
    class = c(
      "BayesTools_bridge_marginal_context",
      "BayesTools_bridge_context",
      "list"
    )
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  fit_data <- .marglik_add_random_covariance_bridge_data(
    fit_data = fit_data,
    model_data = object[["data"]],
    marginalizing = TRUE,
    dependency_blocks = list(1:2, 3L)
  )
  conditional_mean <- as.numeric(.marglik_get_sampling_dependency(
    parameters = parameters,
    model_data = object[["data"]],
    effect_direction = "positive",
    K = nrow(dat)
  ))
  conditional_covariance <- diag(.known_v_residual_variance(known_V)) +
    random_covariance

  actual <- .log_posterior(
    parameters = parameters,
    data = fit_data,
    is_mods = FALSE,
    is_scale = FALSE,
    is_random = TRUE,
    is_multilevel = FALSE,
    is_weights = FALSE,
    is_known_v = TRUE,
    is_PET = FALSE,
    is_PEESE = FALSE,
    is_weightfunction = FALSE,
    effect_direction = "positive",
    outcome_type = "norm",
    model_data = object[["data"]],
    bridge_context = bridge_context
  )
  expected <- .marglik_mvn_log_density(
    y = dat$yi,
    mean = conditional_mean,
    covariance = conditional_covariance
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("latent known-V sampling coordinates are integrated exactly", {

  dat <- data.frame(yi = c(0.10, -0.20, 0.15))
  V <- matrix(
    c(
      0.04, 0.01, 0,
      0.01, 0.09, 0.012,
      0, 0.012, 0.05
    ),
    nrow = 3L,
    byrow = TRUE
  )
  object <- brma.mv(
    yi = yi,
    V = V,
    random = NULL,
    data = dat,
    measure = "GEN",
    known_v_parameterization = "latent",
    marginalize_estimate_level = FALSE,
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
  setup <- .marglik_sampling_latent_setup(
    data       = object[["data"]],
    priors     = object[["priors"]],
    fit_priors = fit_priors
  )

  expect_true(setup[["marginalized"]])
  expect_false("sampling_z" %in% names(setup[["fit_priors"]]))
  expect_identical(
    setup[["diagnostics"]][["included"]],
    paste0("sampling_z[", seq_len(nrow(V)), "]")
  )

  dependency_blocks <- .marglik_random_dependency_blocks(
    model_data                   = object[["data"]],
    formula_design               = NULL,
    blocks                       = character(),
    sampling_latent_marginalized = TRUE
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  fit_data <- .marglik_add_random_covariance_bridge_data(
    fit_data                    = fit_data,
    model_data                  = object[["data"]],
    marginalizing              = TRUE,
    sampling_latent_marginalized = TRUE,
    dependency_blocks          = dependency_blocks
  )
  tau <- 0.12
  actual <- .log_posterior(
    parameters = list(mu = 0, tau = tau),
    data = fit_data,
    is_mods = FALSE,
    is_scale = FALSE,
    is_random = FALSE,
    is_multilevel = FALSE,
    is_weights = FALSE,
    is_known_v = TRUE,
    is_PET = FALSE,
    is_PEESE = FALSE,
    is_weightfunction = FALSE,
    effect_direction = "positive",
    outcome_type = "norm",
    model_data = object[["data"]],
    sampling_latent_marginalized = TRUE
  )
  expected <- .marglik_mvn_log_density(
    y = dat$yi,
    mean = rep(0, nrow(dat)),
    covariance = V + diag(tau^2, nrow(dat))
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("compiled marginalized variance plan matches generic evaluation", {

  dat <- data.frame(
    yi    = c(0.10, -0.20, 0.15, 0.05),
    study = factor(c("s1", "s1", "s2", "s3")),
    esid  = factor(seq_len(4L))
  )
  object <- brma.mv(
    yi = yi,
    V = diag(rep(0.04, nrow(dat))),
    random = ~ 1 | study / esid,
    data = dat,
    measure = "GEN",
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  plan <- .marglik_marginalized_variance_plan(object[["data"]])
  expect_length(plan$terms, 1L)

  sd_name <- plan$terms[[1L]]$parameter
  nodes <- stats::setNames(0.35, sd_name)
  context <- structure(
    list(nodes = nodes),
    class = c("BayesTools_bridge_nodes_context", "BayesTools_bridge_context", "list")
  )
  compiled <- .marglik_evaluate_marginalized_variance_plan(
    plan           = plan,
    bridge_context = context,
    K              = nrow(dat)
  )
  generic <- .evaluate_marginalized_random_variance(
    data = object[["data"]],
    posterior_samples = matrix(
      nodes,
      nrow = 1L,
      dimnames = list(NULL, names(nodes))
    )
  )

  expect_equal(compiled, generic, tolerance = 0)
  expect_error(
    .marglik_evaluate_marginalized_variance_plan(
      plan           = plan,
      bridge_context = structure(
        list(nodes = numeric()),
        class = c("BayesTools_bridge_nodes_context", "BayesTools_bridge_context", "list")
      ),
      K              = nrow(dat)
    ),
    "missing marginalized random-effect SD node",
    fixed = TRUE
  )
})

test_that("native factor likelihood equals independently materialized ZGZ'", {

  dat <- data.frame(
    yi = c(0.10, -0.20, 0.15, 0.05, -0.10),
    study = factor(c("s1", "s1", "s1", "s2", "s2"))
  )
  V <- matrix(0, nrow = 5L, ncol = 5L)
  V[1:3, 1:3] <- matrix(
    c(0.04, 0.01, 0, 0.01, 0.09, 0.02, 0, 0.02, 0.06),
    nrow = 3L,
    byrow = TRUE
  )
  V[4:5, 4:5] <- matrix(c(0.05, 0.01, 0.01, 0.07), nrow = 2L)
  object <- brma.mv(
    yi = yi,
    V = V,
    random = ~ 1 | study,
    data = dat,
    measure = "GEN",
    known_v_parameterization = "block_mvn",
    marginalize_estimate_level = FALSE,
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  Z <- cbind(1, c(-1, 0, 1, 0.5, -0.5))
  G <- matrix(c(0.16, 0.03, 0.03, 0.09), nrow = 2L)
  G_factor <- t(chol(G))
  G <- tcrossprod(G_factor)
  group_map <- c(1L, 1L, 1L, 2L, 2L)
  random_covariance <- (outer(group_map, group_map, "==") * 1) *
    tcrossprod(Z %*% G, Z)
  bridge_context <- structure(
    list(
      nodes = numeric(),
      marginalized_random = list(
        mu = list(
          representation = "factor",
          row_blocks = list(1:3, 4:5),
          factors = list(list(
            type = "group",
            model_matrix = Z,
            group_map = group_map,
            coefficient_covariance = G,
            coefficient_factor = G_factor
          ))
        )
      )
    ),
    class = c(
      "BayesTools_bridge_marginal_context",
      "BayesTools_bridge_context",
      "list"
    )
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  fit_data <- .marglik_add_random_covariance_bridge_data(
    fit_data = fit_data,
    model_data = object[["data"]],
    marginalizing = TRUE,
    dependency_blocks = list(1:3, 4:5)
  )

  actual <- .log_posterior(
    parameters = list(mu = 0),
    data = fit_data,
    is_mods = FALSE,
    is_scale = FALSE,
    is_random = TRUE,
    is_multilevel = FALSE,
    is_weights = FALSE,
    is_known_v = TRUE,
    is_PET = FALSE,
    is_PEESE = FALSE,
    is_weightfunction = FALSE,
    effect_direction = "positive",
    outcome_type = "norm",
    model_data = object[["data"]],
    bridge_context = bridge_context
  )
  expected <- .marglik_mvn_log_density(
    y = dat$yi,
    mean = rep(0, nrow(dat)),
    covariance = V + random_covariance
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("native known-group factor equals independently materialized covariance", {

  y <- c(0.1, -0.2, 0.3, 0.05)
  mean <- c(0.02, 0.02, -0.01, -0.01)
  sampling_covariance <- diag(c(0.04, 0.05, 0.06, 0.07))
  group_map <- c(2L, 1L, 3L, 2L)
  group_covariance <- matrix(
    c(1, 0.3, 0.1, 0.3, 1.5, 0.2, 0.1, 0.2, 0.8),
    nrow = 3L,
    byrow = TRUE
  )
  Z <- matrix(c(1, 0.5, 1.5, -0.5), ncol = 1L)
  sd <- 0.4
  coefficient_factor <- matrix(sd, nrow = 1L, ncol = 1L)
  coefficient_covariance <- tcrossprod(coefficient_factor)
  factor <- list(
    type = "known_group",
    model_matrix = Z,
    group_map = group_map,
    group_covariance = group_covariance,
    coefficient_covariance = coefficient_covariance,
    coefficient_factor = coefficient_factor
  )
  actual <- .Call(
    "RoBMA_known_v_block_mvn_loglik",
    as.double(y),
    as.double(mean),
    sampling_covariance,
    list(factor),
    list(seq_along(y)),
    double(length(y)),
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    group_covariance[group_map, group_map, drop = FALSE] *
    tcrossprod(as.numeric(Z) * sd)
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance
  )

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("native covariance plan reuses exact low-rank group geometry", {

  y <- c(0.1, -0.2, 0.3, 0.05, -0.1, 0.2, -0.05, 0.15)
  mean <- seq(-0.03, 0.04, length.out = length(y))
  sampling_covariance <- diag(seq(0.04, 0.075, length.out = length(y)))
  Z <- cbind(1, seq(-1, 1, length.out = length(y)))
  group_map <- rep(1:2, each = 4L)
  coefficient_factor <- matrix(c(0.4, 0.1, 0, 0.3), nrow = 2L)
  coefficient_covariance <- tcrossprod(coefficient_factor)
  factor <- list(
    type = "group",
    model_matrix = Z,
    group_map = group_map,
    coefficient_covariance = coefficient_covariance,
    coefficient_factor = coefficient_factor
  )
  blocks <- list(1:4, 5:8)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    blocks,
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    double(length(y)),
    PACKAGE = "RoBMA"
  )
  same_group <- outer(group_map, group_map, "==") * 1
  random_covariance <- same_group * tcrossprod(Z %*% coefficient_covariance, Z)
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance
  )

  expect_identical(attr(plan, "low_rank_blocks"), 2L)
  expect_identical(attr(plan, "dense_blocks"), 0L)
  expect_equal(actual, expected, tolerance = 1e-12)

  updated_factor <- factor
  updated_factor$coefficient_factor <- matrix(
    c(0.25, -0.04, 0, 0.2),
    nrow = 2L
  )
  updated_factor$coefficient_covariance <- tcrossprod(
    updated_factor$coefficient_factor
  )
  updated_actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(updated_factor)),
    double(length(y)),
    PACKAGE = "RoBMA"
  )
  updated_random_covariance <- same_group * tcrossprod(
    Z %*% updated_factor$coefficient_covariance,
    Z
  )
  updated_expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + updated_random_covariance
  )

  expect_equal(updated_actual, updated_expected, tolerance = 1e-12)
})

test_that("native diagonal coefficient plans preserve the full covariance", {

  set.seed(701L)
  n <- 12L
  q <- 5L
  y <- stats::rnorm(n)
  mean <- seq(-0.1, 0.15, length.out = n)
  sampling_covariance <- diag(seq(0.06, 0.12, length.out = n))
  Z <- matrix(stats::rnorm(n * q), nrow = n)
  group_map <- rep(1:2, each = n / 2L)
  row_scale <- c(0, seq(0.5, 1.5, length.out = n - 1L))
  scale <- c(0, 0.1, 0.25, 0.4, 0.7)
  coefficient_factor <- diag(scale, q)
  factor <- list(
    type = "row_group",
    model_matrix = Z,
    group_map = group_map,
    row_scale = row_scale,
    coefficient_structure = "diagonal",
    coefficient_factor = coefficient_factor
  )
  blocks <- list(seq_len(n / 2L), (n / 2L + 1L):n)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    blocks,
    PACKAGE = "RoBMA"
  )
  state <- list(.marglik_covariance_factor_state(factor))
  actual_joint <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    state,
    double(n),
    PACKAGE = "RoBMA"
  )
  actual_conditional <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    state,
    double(n),
    PACKAGE = "RoBMA"
  )
  effective_Z <- Z * row_scale
  covariance <- sampling_covariance +
    outer(group_map, group_map, "==") *
      tcrossprod(effective_Z %*% coefficient_factor)
  precision <- solve(covariance)
  precision_residual <- drop(precision %*% (y - mean))
  expected_conditional <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )

  expect_equal(
    actual_joint,
    .marglik_mvn_log_density(y, mean, covariance),
    tolerance = 1e-13
  )
  expect_equal(actual_conditional, expected_conditional, tolerance = 1e-13)

  invalid_factor <- factor
  invalid_factor$coefficient_factor[2L, 1L] <- 0.01
  expect_error(
    .Call(
      "RoBMA_known_v_covariance_plan_create",
      as.double(y),
      sampling_covariance,
      list(invalid_factor),
      blocks,
      PACKAGE = "RoBMA"
    ),
    "must be exactly diagonal",
    fixed = TRUE
  )
  invalid_state <- state
  invalid_state[[1L]]$coefficient_factor[2L, 1L] <- 0.01
  expect_error(
    .Call(
      "RoBMA_known_v_covariance_plan_loglik",
      plan,
      as.double(mean),
      invalid_state,
      double(n),
      PACKAGE = "RoBMA"
    ),
    "must be exactly diagonal",
    fixed = TRUE
  )
})

test_that("native Markov plans match dense joint and conditional likelihoods", {

  transition <- c(-0.45, -0.45, -0.45, 0.35, 0.6, 0.25, 0.7)
  scale <- c(0, 0.2, 0.3, 0.15, 0.4, 0.25, 0.1, 0.35)
  level <- c(8L, 1L, 4L, 4L, 6L, 2L)
  loading <- c(1, -0.5, 1.2, 0.7, -1.1, 0.4)
  row_scale <- c(1, 0.5, 1.3, 0, 0.8, 1.1)
  n <- length(level)
  q <- length(scale)
  Z <- matrix(0, nrow = n, ncol = q)
  Z[cbind(seq_len(n), level)] <- loading

  correlation <- diag(q)
  for (left in seq_len(q - 1L)) {
    correlation_value <- 1
    for (right in (left + 1L):q) {
      correlation_value <- correlation_value * transition[[right - 1L]]
      correlation[left, right] <- correlation_value
      correlation[right, left] <- correlation_value
    }
  }
  coefficient_factor <- t(chol(correlation)) * scale
  sampling_covariance <- diag(seq(0.15, 0.3, length.out = n))
  extra_variance <- seq(0.01, 0.04, length.out = n)
  y <- c(0.2, -0.7, 0.4, 1.1, -0.3, 0.8)
  mean <- c(-0.1, 0.2, 0, 0.15, -0.05, 0.1)
  factor <- list(
    type = "row_group",
    model_matrix = Z,
    group_map = rep.int(1L, n),
    coefficient_structure = "markov",
    coefficient_factor = coefficient_factor,
    coefficient_scale = scale,
    markov_transition = transition,
    markov_innovation_variance = 1 - transition^2,
    row_scale = row_scale
  )
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(seq_len(n)),
    PACKAGE = "RoBMA"
  )
  state <- .marglik_covariance_factor_state(factor)
  actual_joint <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(state),
    as.double(extra_variance),
    PACKAGE = "RoBMA"
  )
  actual_conditional <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    list(state),
    as.double(extra_variance),
    PACKAGE = "RoBMA"
  )
  effective_Z <- Z * row_scale
  covariance <- sampling_covariance + diag(extra_variance) +
    tcrossprod(effective_Z %*% coefficient_factor)
  root <- chol(covariance)
  precision <- chol2inv(root)
  residual <- y - mean
  precision_residual <- as.vector(precision %*% residual)
  expected_conditional <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )

  expect_identical(attr(plan, "markov_blocks"), 1L)
  expect_equal(
    actual_joint,
    .marglik_mvn_log_density(y, mean, covariance),
    tolerance = 1e-13
  )
  expect_equal(actual_conditional, expected_conditional, tolerance = 1e-13)

  correlated_sampling <- sampling_covariance
  correlated_sampling[1L, 2L] <- 0.01
  correlated_sampling[2L, 1L] <- 0.01
  fallback_plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    correlated_sampling,
    list(factor),
    list(seq_len(n)),
    PACKAGE = "RoBMA"
  )
  fallback_covariance <- correlated_sampling + diag(extra_variance) +
    tcrossprod(effective_Z %*% coefficient_factor)
  fallback_joint <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    fallback_plan,
    as.double(mean),
    list(state),
    as.double(extra_variance),
    PACKAGE = "RoBMA"
  )

  expect_identical(attr(fallback_plan, "markov_blocks"), 0L)
  expect_equal(
    fallback_joint,
    .marglik_mvn_log_density(y, mean, fallback_covariance),
    tolerance = 1e-13
  )

  dense_contribution <- diag(seq(0.005, 0.01, length.out = n))
  multi_factor_plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    c(list(factor), list(list(
      type = "dense",
      covariance = dense_contribution
    ))),
    list(seq_len(n)),
    PACKAGE = "RoBMA"
  )
  multi_factor_joint <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    multi_factor_plan,
    as.double(mean),
    list(state, list(covariance = dense_contribution)),
    as.double(extra_variance),
    PACKAGE = "RoBMA"
  )

  expect_identical(attr(multi_factor_plan, "markov_blocks"), 0L)
  expect_equal(
    multi_factor_joint,
    .marglik_mvn_log_density(
      y,
      mean,
      covariance + dense_contribution
    ),
    tolerance = 1e-13
  )
})

test_that("native covariance plan assembles sparse nested latent geometry", {

  K <- 12L
  y <- seq(-0.2, 0.35, length.out = K)
  mean <- seq(0.03, -0.02, length.out = K)
  sampling_variance <- seq(0.04, 0.095, length.out = K)
  parent_map <- rep(1:2, each = 6L)
  child_map <- rep(1:4, each = 3L)
  parent_factor <- list(
    type = "group",
    model_matrix = matrix(1, nrow = K, ncol = 1L),
    group_map = parent_map,
    coefficient_covariance = matrix(0.3^2),
    coefficient_factor = matrix(0.3)
  )
  child_factor <- list(
    type = "group",
    model_matrix = matrix(1, nrow = K, ncol = 1L),
    group_map = child_map,
    coefficient_covariance = matrix(0.2^2),
    coefficient_factor = matrix(0.2)
  )
  factors <- list(parent_factor, child_factor)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    diag(sampling_variance),
    factors,
    list(seq_len(K)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    lapply(factors, .marglik_covariance_factor_state),
    double(K),
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    0.3^2 * outer(parent_map, parent_map, "==") +
    0.2^2 * outer(child_map, child_map, "==")
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = diag(sampling_variance) + random_covariance
  )
  precision <- solve(diag(sampling_variance) + random_covariance)
  precision_residual <- drop(precision %*% (y - mean))
  conditional_expected <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )
  conditional_actual <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    lapply(factors, .marglik_covariance_factor_state),
    double(K),
    PACKAGE = "RoBMA"
  )

  expect_identical(attr(plan, "low_rank_blocks"), 1L)
  expect_identical(attr(plan, "sparse_assembly_blocks"), 1L)
  expect_identical(attr(plan, "sparse_factor_blocks"), 1L)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(conditional_actual, conditional_expected, tolerance = 1e-12)

  updated_factors <- factors
  updated_factors[[1L]]$coefficient_factor <- matrix(0.22)
  updated_factors[[1L]]$coefficient_covariance <- matrix(0.22^2)
  updated_factors[[2L]]$coefficient_factor <- matrix(0.14)
  updated_factors[[2L]]$coefficient_covariance <- matrix(0.14^2)
  updated_actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    lapply(updated_factors, .marglik_covariance_factor_state),
    double(K),
    PACKAGE = "RoBMA"
  )
  updated_random_covariance <-
    0.22^2 * outer(parent_map, parent_map, "==") +
    0.14^2 * outer(child_map, child_map, "==")
  updated_expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = diag(sampling_variance) + updated_random_covariance
  )

  expect_equal(updated_actual, updated_expected, tolerance = 1e-12)
  expect_equal(
    .Call(
      "RoBMA_known_v_covariance_plan_loglik",
      plan,
      as.double(mean),
      lapply(factors, .marglik_covariance_factor_state),
      double(K),
      PACKAGE = "RoBMA"
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("native sparse factor supports more latent groups than observations", {

  K <- 6L
  y <- seq(-0.15, 0.25, length.out = K)
  mean <- seq(0.02, -0.03, length.out = K)
  sampling_variance <- seq(0.04, 0.065, length.out = K)
  parent_map <- rep(1:2, each = 3L)
  child_map <- seq_len(K)
  parent_factor <- list(
    type = "group",
    model_matrix = matrix(1, nrow = K, ncol = 1L),
    group_map = parent_map,
    coefficient_covariance = matrix(0.3^2),
    coefficient_factor = matrix(0.3)
  )
  child_factor <- list(
    type = "group",
    model_matrix = matrix(1, nrow = K, ncol = 1L),
    group_map = child_map,
    coefficient_covariance = matrix(0.2^2),
    coefficient_factor = matrix(0.2)
  )
  factors <- list(parent_factor, child_factor)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    diag(sampling_variance),
    factors,
    list(seq_len(K)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    lapply(factors, .marglik_covariance_factor_state),
    double(K),
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    0.3^2 * outer(parent_map, parent_map, "==") +
    0.2^2 * diag(K)
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = diag(sampling_variance) + random_covariance
  )
  precision <- solve(diag(sampling_variance) + random_covariance)
  precision_residual <- drop(precision %*% (y - mean))
  conditional_expected <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )
  conditional_actual <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    lapply(factors, .marglik_covariance_factor_state),
    double(K),
    PACKAGE = "RoBMA"
  )

  expect_gt(
    length(unique(parent_map)) + length(unique(child_map)),
    K
  )
  expect_identical(attr(plan, "sparse_factor_blocks"), 1L)
  expect_identical(attr(plan, "root_dense_blocks"), 0L)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(conditional_actual, conditional_expected, tolerance = 1e-12)

})

test_that("native Woodbury quadratic remains stable for explained residuals", {

  K <- 12L
  direction <- seq(-1, 1, length.out = K)
  direction <- direction / sqrt(sum(direction ^ 2))
  orthogonal <- rep(c(-1, 1), length.out = K)
  orthogonal <- orthogonal -
    sum(orthogonal * direction) * direction
  orthogonal <- orthogonal / sqrt(sum(orthogonal ^ 2))
  random_root <- 1e6 * direction
  y <- 0.7 * random_root + 0.25 * orthogonal
  factor <- list(
    type                   = "group",
    model_matrix           = matrix(random_root, ncol = 1L),
    group_map              = rep(1L, K),
    coefficient_covariance = matrix(1, nrow = 1L, ncol = 1L),
    coefficient_factor     = matrix(1, nrow = 1L, ncol = 1L)
  )
  sampling_covariance <- diag(1, nrow = K)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(seq_len(K)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    double(K),
    list(.marglik_covariance_factor_state(factor)),
    double(K),
    PACKAGE = "RoBMA"
  )
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = double(K),
    covariance = sampling_covariance + tcrossprod(random_root)
  )

  expect_identical(attr(plan, "low_rank_blocks"), 1L)
  expect_equal(actual, expected, tolerance = 1e-5)
})

test_that("native covariance plan uses exact spectral Woodbury blocks", {

  y <- c(0.1, -0.2, 0.3, 0.05, -0.1, 0.2, -0.05, 0.15)
  mean <- seq(-0.03, 0.04, length.out = length(y))
  base_block <- matrix(
    c(
      0.08, 0.012, 0.006, 0.004,
      0.012, 0.07, 0.009, 0.005,
      0.006, 0.009, 0.06, 0.008,
      0.004, 0.005, 0.008, 0.05
    ),
    nrow = 4L,
    byrow = TRUE
  )
  sampling_covariance <- matrix(0, nrow = length(y), ncol = length(y))
  sampling_covariance[1:4, 1:4] <- base_block
  sampling_covariance[5:8, 5:8] <- 1.2 * base_block
  Z <- matrix(1, nrow = length(y), ncol = 1L)
  group_map <- rep(1:2, each = 4L)
  coefficient_factor <- matrix(0.3, nrow = 1L, ncol = 1L)
  factor <- list(
    type = "group",
    model_matrix = Z,
    group_map = group_map,
    coefficient_covariance = tcrossprod(coefficient_factor),
    coefficient_factor = coefficient_factor
  )
  blocks <- list(1:4, 5:8)
  extra_variance <- rep(c(0.01, 0.02), each = 4L)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    blocks,
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    extra_variance,
    PACKAGE = "RoBMA"
  )
  same_group <- outer(group_map, group_map, "==") * 1
  random_covariance <- same_group * tcrossprod(
    Z %*% factor$coefficient_covariance,
    Z
  )
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance +
      diag(extra_variance)
  )

  expect_identical(attr(plan, "low_rank_blocks"), 0L)
  expect_identical(attr(plan, "spectral_blocks"), 2L)
  expect_identical(attr(plan, "dense_blocks"), 0L)
  expect_equal(actual, expected, tolerance = 1e-12)

  varying_extra <- extra_variance
  varying_extra[[2L]] <- 0.011
  fallback <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    varying_extra,
    PACKAGE = "RoBMA"
  )
  fallback_expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance +
      diag(varying_extra)
  )
  expect_equal(fallback, fallback_expected, tolerance = 1e-12)
})

test_that("native covariance plan uses exact block-base Woodbury algebra", {

  K <- 12L
  y <- sin(seq(0.2, 2.4, length.out = K))
  mean <- seq(-0.05, 0.06, length.out = K)
  sampling_covariance <- matrix(0, nrow = K, ncol = K)
  sampling_block <- matrix(
    c(
      0.08, 0.012, 0.006, 0.004,
      0.012, 0.07, 0.009, 0.005,
      0.006, 0.009, 0.06, 0.008,
      0.004, 0.005, 0.008, 0.05
    ),
    nrow = 4L,
    byrow = TRUE
  )
  for (block_i in seq_len(3L)) {
    index <- 4L * (block_i - 1L) + seq_len(4L)
    sampling_covariance[index, index] <- block_i * sampling_block
  }
  group_map <- rep(1:2, length.out = K)
  coefficient_factor <- matrix(0.25, nrow = 1L, ncol = 1L)
  factor <- list(
    type = "group",
    model_matrix = matrix(1, nrow = K, ncol = 1L),
    group_map = group_map,
    coefficient_covariance = tcrossprod(coefficient_factor),
    coefficient_factor = coefficient_factor
  )
  extra_variance <- seq(0.005, 0.016, length.out = K)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(seq_len(K)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    extra_variance,
    PACKAGE = "RoBMA"
  )
  random_covariance <- outer(group_map, group_map, "==") * 0.25^2
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance +
      diag(extra_variance)
  )
  covariance <- sampling_covariance + random_covariance +
    diag(extra_variance)
  precision <- solve(covariance)
  precision_residual <- drop(precision %*% (y - mean))
  conditional_expected <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )
  conditional_actual <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    extra_variance,
    PACKAGE = "RoBMA"
  )

  expect_identical(attr(plan, "low_rank_blocks"), 0L)
  expect_identical(attr(plan, "block_base_blocks"), 1L)
  expect_identical(attr(plan, "spectral_blocks"), 0L)
  expect_identical(attr(plan, "dense_blocks"), 0L)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(conditional_actual, conditional_expected, tolerance = 1e-12)

  zero_extra <- double(K)
  zero_covariance <- sampling_covariance + random_covariance
  zero_precision <- solve(zero_covariance)
  zero_precision_residual <- drop(zero_precision %*% (y - mean))
  zero_expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = zero_covariance
  )
  zero_conditional_expected <- 0.5 * (
    log(diag(zero_precision)) - log(2 * pi) -
      zero_precision_residual^2 / diag(zero_precision)
  )
  zero_actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    zero_extra,
    PACKAGE = "RoBMA"
  )
  zero_conditional_actual <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    zero_extra,
    PACKAGE = "RoBMA"
  )
  repeated_actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    extra_variance,
    PACKAGE = "RoBMA"
  )

  expect_equal(zero_actual, zero_expected, tolerance = 1e-12)
  expect_equal(
    zero_conditional_actual,
    zero_conditional_expected,
    tolerance = 1e-12
  )
  expect_identical(repeated_actual, actual)
})

test_that("native sparse factor supports correlated sampling blocks", {

  K <- 8L
  y <- sin(seq(0.15, 1.2, length.out = K))
  mean <- seq(-0.04, 0.05, length.out = K)
  sampling_covariance <- matrix(0, nrow = K, ncol = K)
  for (block_i in seq_len(K / 2L)) {
    index <- 2L * (block_i - 1L) + 1:2
    sampling_covariance[index, index] <- matrix(
      c(0.06, 0.012, 0.012, 0.08),
      nrow = 2L
    )
  }
  make_factor <- function(group_map, standard_deviation) {

    list(
      type = "group",
      model_matrix = matrix(1, nrow = K, ncol = 1L),
      group_map = group_map,
      coefficient_covariance = matrix(standard_deviation^2),
      coefficient_factor = matrix(standard_deviation)
    )
  }
  parent_map <- rep(1:2, each = K / 2L)
  child_map <- seq_len(K)
  crossed_map <- rep(1:2, length.out = K)
  factors <- list(
    make_factor(parent_map, 0.30),
    make_factor(child_map, 0.20),
    make_factor(crossed_map, 0.15)
  )
  extra_variance <- seq(0.004, 0.011, length.out = K)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    factors,
    list(seq_len(K)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    lapply(factors, .marglik_covariance_factor_state),
    extra_variance,
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    0.30^2 * outer(parent_map, parent_map, "==") +
    0.20^2 * diag(K) +
    0.15^2 * outer(crossed_map, crossed_map, "==")
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance +
      diag(extra_variance)
  )
  covariance <- sampling_covariance + random_covariance +
    diag(extra_variance)
  precision <- solve(covariance)
  precision_residual <- drop(precision %*% (y - mean))
  conditional_expected <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )
  conditional_actual <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    lapply(factors, .marglik_covariance_factor_state),
    extra_variance,
    PACKAGE = "RoBMA"
  )

  expect_gt(
    length(unique(parent_map)) + length(unique(child_map)) +
      length(unique(crossed_map)),
    K
  )
  expect_identical(attr(plan, "block_base_blocks"), 1L)
  expect_identical(attr(plan, "sparse_factor_blocks"), 1L)
  expect_identical(attr(plan, "root_dense_blocks"), 0L)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(conditional_actual, conditional_expected, tolerance = 1e-12)
})

test_that("block-base Woodbury falls back exactly for singular base blocks", {

  K <- 4L
  y <- c(0.1, -0.2, 0.3, -0.1)
  sampling_covariance <- matrix(0, nrow = K, ncol = K)
  sampling_covariance[1:2, 1:2] <- 0.05
  sampling_covariance[3:4, 3:4] <- 0.08
  Z <- cbind(
    1,
    c(1, -1, 0, 0),
    c(0, 0, 1, -1)
  )
  coefficient_factor <- diag(c(0.2, 0.3, 0.25))
  factor <- list(
    type = "group",
    model_matrix = Z,
    group_map = rep(1L, K),
    coefficient_covariance = tcrossprod(coefficient_factor),
    coefficient_factor = coefficient_factor
  )
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(seq_len(K)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    double(K),
    list(.marglik_covariance_factor_state(factor)),
    double(K),
    PACKAGE = "RoBMA"
  )
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = double(K),
    covariance = sampling_covariance + tcrossprod(Z %*% coefficient_factor)
  )

  expect_identical(attr(plan, "block_base_blocks"), 1L)
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("compiled factor validation retains changing-value checks", {

  K <- 4L
  coefficient_factor <- matrix(c(0.4, 0.1, 0, 0.3), nrow = 2L)
  factor <- list(
    type = "group",
    model_matrix = cbind(1, c(-1, 0, 1, 2)),
    group_map = c(1L, 1L, 2L, 2L),
    coefficient_covariance = tcrossprod(coefficient_factor),
    coefficient_factor = coefficient_factor
  )
  context <- list(marginalized_random = list(mu = list(
    representation = "factor",
    row_blocks = list(1:2, 3:4),
    factors = list(factor)
  )))
  cache <- new.env(parent = emptyenv())

  first <- .marglik_bridge_random_covariance(
    bridge_context   = context,
    K                = K,
    validation_cache = cache
  )
  updated_factor <- factor
  updated_factor$coefficient_factor <- matrix(
    c(0.25, -0.04, 0, 0.2),
    nrow = 2L
  )
  updated_factor$coefficient_covariance <- tcrossprod(
    updated_factor$coefficient_factor
  )
  updated <- context
  updated$marginalized_random$mu$factors[[1L]] <- updated_factor
  second <- .marglik_bridge_random_covariance(
    bridge_context   = updated,
    K                = K,
    validation_cache = cache
  )

  expect_equal(
    second$factors[[1L]]$coefficient_covariance,
    updated_factor$coefficient_covariance,
    tolerance = 0
  )
  expect_identical(second$row_blocks, first$row_blocks)

  asymmetric <- updated
  asymmetric_factor <- asymmetric$marginalized_random$mu$factors[[1L]]
  asymmetric_factor$coefficient_covariance[1L, 2L] <- 0.123
  asymmetric$marginalized_random$mu$factors[[1L]] <- asymmetric_factor
  expect_error(
    .marglik_bridge_random_covariance(
      bridge_context   = asymmetric,
      K                = K,
      validation_cache = cache
    ),
    "must be symmetric",
    fixed = TRUE
  )

  changed_blocks <- updated
  changed_blocks$marginalized_random$mu$row_blocks <- list(1:3, 4L)
  expect_error(
    .marglik_bridge_random_covariance(
      bridge_context   = changed_blocks,
      K                = K,
      validation_cache = cache
    ),
    "row blocks changed",
    fixed = TRUE
  )
})

test_that("compact bridge factor states retain the exact covariance contract", {

  K <- 4L
  coefficient_factor <- matrix(c(0.4, 0.1, 0, 0.3), nrow = 2L)
  factor_plan <- list(
    type = "row_group",
    model_matrix = cbind(1, c(-1, 0, 1, 2)),
    group_map = c(1L, 1L, 2L, 2L)
  )
  contract_id <- new.env(parent = emptyenv())
  context <- list(marginalized_random = list(mu = list(
    representation = "factor_state",
    contract_id = contract_id,
    row_blocks = list(1:2, 3:4),
    factor_plans = list(factor_plan),
    factor_states = list(list(
      coefficient_factor = coefficient_factor,
      row_scale = c(0.8, 1.1, 0.9, 1.2)
    ))
  )))
  cache <- new.env(parent = emptyenv())

  first <- .marglik_bridge_random_covariance(
    bridge_context = context,
    K = K,
    validation_cache = cache
  )
  updated <- context
  updated_factor <- matrix(c(0.25, -0.04, 0, 0.2), nrow = 2L)
  updated$marginalized_random$mu$factor_states[[1L]]$coefficient_factor <-
    updated_factor
  second <- .marglik_bridge_random_covariance(
    bridge_context = updated,
    K = K,
    validation_cache = cache
  )

  expect_identical(first$representation, "factor")
  expect_identical(second$factors[[1L]]$model_matrix, factor_plan$model_matrix)
  expect_identical(second$factors[[1L]]$coefficient_factor, updated_factor)
  expect_null(second$factors[[1L]]$coefficient_covariance)

  invalid <- updated
  invalid$marginalized_random$mu$factor_states[[1L]]$row_scale[[2L]] <- -1
  expect_error(
    .marglik_bridge_random_covariance(
      bridge_context = invalid,
      K = K,
      validation_cache = cache
    ),
    "row scale is invalid",
    fixed = TRUE
  )

  changed_contract <- updated
  changed_contract$marginalized_random$mu$contract_id <-
    new.env(parent = emptyenv())
  expect_error(
    .marglik_bridge_random_covariance(
      bridge_context = changed_contract,
      K = K,
      validation_cache = cache
    ),
    "contract changed",
    fixed = TRUE
  )
})

test_that("native low-rank plan supports row-specific external scales", {

  y <- c(0.1, -0.1, 0.2, -0.2, 0.05, 0.15)
  mean <- rep(0.02, length(y))
  sampling_covariance <- diag(seq(0.03, 0.055, length.out = length(y)))
  Z <- cbind(1, seq(-0.5, 0.5, length.out = length(y)))
  group_map <- rep(1:2, each = 3L)
  row_scale <- c(0.8, 1.1, 0.9, 1.2, 0.7, 1.05)
  coefficient_factor <- matrix(c(0.3, 0.05, 0, 0.2), nrow = 2L)
  coefficient_covariance <- tcrossprod(coefficient_factor)
  factor <- list(
    type = "row_group",
    model_matrix = Z,
    group_map = group_map,
    row_scale = row_scale,
    coefficient_covariance = coefficient_covariance,
    coefficient_factor = coefficient_factor
  )
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(1:3, 4:6),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    double(length(y)),
    PACKAGE = "RoBMA"
  )
  scaled_Z <- Z * row_scale
  random_covariance <- (outer(group_map, group_map, "==") * 1) *
    tcrossprod(scaled_Z %*% coefficient_covariance, scaled_Z)
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance
  )

  expect_identical(attr(plan, "low_rank_blocks"), 2L)
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("native dense plan supports known group covariance and random slopes", {

  y <- c(0.1, -0.2, 0.3, 0.05, -0.1, 0.15)
  mean <- rep(c(0.02, -0.01), each = 3L)
  sampling_covariance <- diag(seq(0.04, 0.065, length.out = length(y)))
  Z <- cbind(1, seq(-1, 1, length.out = length(y)))
  group_map <- c(1L, 1L, 2L, 2L, 3L, 3L)
  group_covariance <- matrix(
    c(1, 0.3, 0.1, 0.3, 1.2, 0.2, 0.1, 0.2, 0.8),
    nrow = 3L,
    byrow = TRUE
  )
  coefficient_factor <- matrix(c(0.35, 0.08, 0, 0.25), nrow = 2L)
  coefficient_covariance <- tcrossprod(coefficient_factor)
  factor <- list(
    type = "known_group",
    model_matrix = Z,
    group_map = group_map,
    group_covariance = group_covariance,
    coefficient_covariance = coefficient_covariance,
    coefficient_factor = coefficient_factor
  )
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(seq_along(y)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    double(length(y)),
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    group_covariance[group_map, group_map, drop = FALSE] *
    tcrossprod(Z %*% coefficient_covariance, Z)
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance
  )

  expect_identical(attr(plan, "low_rank_blocks"), 0L)
  expect_identical(attr(plan, "root_dense_blocks"), 1L)
  expect_identical(attr(plan, "dense_blocks"), 1L)
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("native plan factors repeated-observation known group covariance", {

  y <- c(0.1, -0.2, 0.3, 0.05, -0.1, 0.15, -0.05, 0.2)
  mean <- seq(-0.02, 0.05, length.out = length(y))
  sampling_covariance <- diag(seq(0.04, 0.075, length.out = length(y)))
  Z <- cbind(1, seq(-1, 1, length.out = length(y)))
  group_map <- rep(1:2, each = 4L)
  group_covariance <- matrix(c(1, 0.3, 0.3, 1.2), nrow = 2L)
  coefficient_factor <- matrix(c(0.35, 0.08, 0, 0.25), nrow = 2L)
  coefficient_covariance <- tcrossprod(coefficient_factor)
  factor <- list(
    type = "known_group",
    model_matrix = Z,
    group_map = group_map,
    group_covariance = group_covariance,
    coefficient_covariance = coefficient_covariance,
    coefficient_factor = coefficient_factor
  )
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(factor),
    list(seq_along(y)),
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    list(.marglik_covariance_factor_state(factor)),
    double(length(y)),
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    group_covariance[group_map, group_map, drop = FALSE] *
    tcrossprod(Z %*% coefficient_covariance, Z)
  expected <- .marglik_mvn_log_density(
    y = y,
    mean = mean,
    covariance = sampling_covariance + random_covariance
  )

  expect_identical(attr(plan, "low_rank_blocks"), 1L)
  expect_identical(attr(plan, "dense_blocks"), 0L)
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("native known group cache preserves diagonal coefficient models", {

  set.seed(811L)
  n <- 12L
  y <- stats::rnorm(n)
  mean <- seq(-0.2, 0.15, length.out = n)
  sampling_variance <- seq(0.04, 0.09, length.out = n)
  Z <- cbind(1, seq(-1, 1, length.out = n))
  group_map <- rep(1:3, each = 4L)
  group_covariance <- matrix(
    c(1.2, 0.3, 0.1, 0.3, 0.9, -0.2, 0.1, -0.2, 1.1),
    nrow = 3L,
    byrow = TRUE
  )
  coefficient_factor <- diag(c(0, 0.35))
  factor <- list(
    type = "known_group",
    model_matrix = Z,
    group_map = group_map,
    group_covariance = group_covariance,
    coefficient_structure = "diagonal",
    coefficient_factor = coefficient_factor
  )
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    diag(sampling_variance),
    list(factor),
    list(seq_len(n)),
    PACKAGE = "RoBMA"
  )
  state <- list(list(coefficient_factor = coefficient_factor))
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    state,
    double(n),
    PACKAGE = "RoBMA"
  )
  actual_conditional <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    state,
    double(n),
    PACKAGE = "RoBMA"
  )
  random_covariance <-
    group_covariance[group_map, group_map, drop = FALSE] *
    tcrossprod(Z %*% coefficient_factor)
  covariance <- diag(sampling_variance) + random_covariance
  precision <- solve(covariance)
  precision_residual <- drop(precision %*% (y - mean))
  expected_conditional <- 0.5 * (
    log(diag(precision)) - log(2 * pi) -
      precision_residual^2 / diag(precision)
  )

  expect_identical(attr(plan, "fixed_known_group_blocks"), 1L)
  expect_equal(
    actual,
    .marglik_mvn_log_density(y, mean, covariance),
    tolerance = 1e-12
  )
  expect_equal(actual_conditional, expected_conditional, tolerance = 1e-12)

  extra_variance <- seq(0.001, 0.012, length.out = n)
  fallback <- .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    as.double(mean),
    state,
    extra_variance,
    PACKAGE = "RoBMA"
  )
  expect_equal(
    fallback,
    .marglik_mvn_log_density(
      y,
      mean,
      covariance + diag(extra_variance)
    ),
    tolerance = 1e-12
  )

  correlated_sampling <- diag(sampling_variance)
  correlated_sampling[1L, 2L] <- 0.01
  correlated_sampling[2L, 1L] <- 0.01
  fallback_plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    correlated_sampling,
    list(factor),
    list(seq_len(n)),
    PACKAGE = "RoBMA"
  )
  expect_identical(attr(fallback_plan, "fixed_known_group_blocks"), 0L)
  expect_equal(
    .Call(
      "RoBMA_known_v_covariance_plan_loglik",
      fallback_plan,
      as.double(mean),
      state,
      double(n),
      PACKAGE = "RoBMA"
    ),
    .marglik_mvn_log_density(
      y,
      mean,
      correlated_sampling + random_covariance
    ),
    tolerance = 1e-12
  )
})

test_that("marglik requests every fitted sampled random block", {

  dat <- data.frame(
    yi = c(0.10, -0.20, 0.15, 0.05),
    study = factor(c("s1", "s1", "s2", "s3")),
    esid = factor(seq_len(4L))
  )
  object <- brma.mv(
    yi = yi,
    V = diag(rep(0.04, nrow(dat))),
    random = ~ 1 | study / esid,
    data = dat,
    measure = "GEN",
    prior_unit_information_sd = 1,
    only_priors = TRUE
  )
  formula_args <- .create_jags_formula_args(
    data = object[["data"]],
    priors = object[["priors"]]
  )
  mu_design <- BayesTools::JAGS_formula(
    formula = formula_args$formula_list$mu,
    parameter = "mu",
    data = formula_args$formula_data_list$mu,
    prior_list = formula_args$formula_prior_list$mu,
    formula_scale = formula_args$formula_scale_list$mu,
    prior_random = formula_args$formula_random_prior_list$mu,
    random_effects_compile = formula_args$formula_random_effects_compile_list$mu
  )$formula_design
  fit <- list()
  attr(fit, "formula_design") <- list(mu = mu_design)
  setup <- .marglik_bridge_random_marginalization(
    object = object,
    fit = fit,
    fixed_zero_random = FALSE
  )
  expected <- .formula_design_sampled_random_effect_blocks(
    mu_design
  )

  expect_length(expected, 1L)
  expect_identical(setup$blocks$mu, expected)
  expect_true(setup$request$mu$factor_state)
  expect_identical(setup$dependency_blocks, list(1:2, 3L, 4L))
  expect_identical(setup$diagnostics$included, expected)
  expect_true(setup$diagnostics$exact)
  expect_equal(
    setup$diagnostics$skipped$reason,
    "already marginalized by the fitted likelihood"
  )
})

test_that("marglik accepts every brma.mv random covariance structure", {

  dat <- data.frame(
    yi = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2"),
    x = c(0, 1, 0, 1),
    out = c("a", "b", "a", "b"),
    time = c(1, 2, 1, 2)
  )
  formulas <- list(
    id   = ~ id(1 | study),
    diag = ~ diag(1 + x | study),
    us   = ~ us(1 + x | study),
    cs   = ~ cs(out | study),
    hcs  = ~ hcs(out | study),
    ar1  = ~ ar1(time | study),
    har  = ~ har(time | study),
    car  = ~ car(time | study)
  )

  for (structure in names(formulas)) {
    object <- brma.mv(
      yi = yi,
      V = diag(rep(0.04, nrow(dat))),
      random = formulas[[structure]],
      data = dat,
      measure = "GEN",
      prior_unit_information_sd = 1,
      marginalize_estimate_level = FALSE,
      only_priors = TRUE
    )
    formula_args <- .create_jags_formula_args(
      data = object[["data"]],
      priors = object[["priors"]]
    )
    mu_design <- BayesTools::JAGS_formula(
      formula = formula_args$formula_list$mu,
      parameter = "mu",
      data = formula_args$formula_data_list$mu,
      prior_list = formula_args$formula_prior_list$mu,
      formula_scale = formula_args$formula_scale_list$mu,
      prior_random = formula_args$formula_random_prior_list$mu,
      random_effects_compile =
        formula_args$formula_random_effects_compile_list$mu
    )$formula_design
    fit <- list()
    attr(fit, "formula_design") <- list(mu = mu_design)
    setup <- .marglik_bridge_random_marginalization(
      object = object,
      fit = fit,
      fixed_zero_random = FALSE
    )

    expect_identical(
      setup$blocks$mu,
      .formula_design_sampled_random_effect_blocks(mu_design),
      info = structure
    )
    expect_true(setup$diagnostics$exact, info = structure)
  }
})

test_that("marglik accepts sampled random intercepts with known group covariance", {

  dat <- data.frame(
    yi = c(0.10, 0.20, 0.30, 0.40),
    study = c("s1", "s1", "s2", "s2")
  )
  R <- matrix(
    c(1, 0.35, 0.35, 1),
    nrow = 2L,
    dimnames = list(c("s1", "s2"), c("s1", "s2"))
  )
  object <- brma.mv(
    yi = yi,
    V = diag(rep(0.04, nrow(dat))),
    random = ~ 1 | study,
    R = list(study = R),
    data = dat,
    measure = "GEN",
    prior_unit_information_sd = 1,
    marginalize_estimate_level = FALSE,
    only_priors = TRUE
  )
  formula_args <- .create_jags_formula_args(
    data = object[["data"]],
    priors = object[["priors"]]
  )
  mu_design <- BayesTools::JAGS_formula(
    formula = formula_args$formula_list$mu,
    parameter = "mu",
    data = formula_args$formula_data_list$mu,
    prior_list = formula_args$formula_prior_list$mu,
    formula_scale = formula_args$formula_scale_list$mu,
    prior_random = formula_args$formula_random_prior_list$mu,
    random_effects_compile =
      formula_args$formula_random_effects_compile_list$mu
  )$formula_design
  fit <- list()
  attr(fit, "formula_design") <- list(mu = mu_design)
  setup <- .marglik_bridge_random_marginalization(
    object = object,
    fit = fit,
    fixed_zero_random = FALSE
  )

  expect_identical(
    setup$blocks$mu,
    .formula_design_sampled_random_effect_blocks(mu_design)
  )
  expect_identical(setup$dependency_blocks, list(seq_len(nrow(dat))))
  expect_true(setup$diagnostics$exact)
})
