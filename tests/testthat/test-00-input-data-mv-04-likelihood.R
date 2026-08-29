context("brma.mv known-V likelihood contracts")

test_that("brma.mv known-V models pass bridge marginal likelihood availability", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = diag(c(0.04, 0.09)),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_marglik_available(object, "test()")
  )
})


test_that("bridge SD-source spec expands row sources and attaches bounds", {

  posterior <- matrix(
    1:6,
    nrow = 2L,
    dimnames = list(NULL, c("tau[1]", "tau[2]", "sigma"))
  )

  spec <- .marglik_bridge_sd_source_spec(
    add_parameters = c("tau", "sigma"),
    fit            = posterior,
    K              = 2L
  )

  expect_equal(spec[["parameters"]], c("tau[1]", "tau[2]", "sigma"))
  expect_equal(names(spec[["bounds"]][["lb"]]), spec[["parameters"]])
  expect_equal(names(spec[["bounds"]][["ub"]]), spec[["parameters"]])
  expect_equal(unname(spec[["bounds"]][["lb"]]), c(0, 0, 0))
  expect_true(all(is.infinite(spec[["bounds"]][["ub"]])))

  empty <- .marglik_bridge_sd_source_spec(
    add_parameters = character(),
    fit            = posterior,
    K              = 2L
  )
  expect_equal(empty[["parameters"]], character())
  expect_null(empty[["bounds"]])
})


test_that("brma.mv known-V bridge log posterior matches exact normal targets", {

  V <- matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2)
  parameters <- list(mu = 0, tau = 0.10)

  block_mvn <- brma.mv(
    yi                       = c(0.10, 0.20),
    V                        = V,
    known_v_parameterization = "block_mvn",
    measure                  = "GEN",
    prior_unit_information_sd = 1,
    only_priors              = TRUE
  )
  block_fit_data <- .create_fit_data(block_mvn[["data"]], block_mvn[["priors"]])
  block_expected <- .marglik_mvn_log_density(
    y          = c(0.10, 0.20),
    mean       = c(0, 0),
    covariance = V + diag(parameters[["tau"]]^2, nrow = 2)
  )
  expect_equal(
    .log_posterior(
      parameters                 = parameters,
      data                       = block_fit_data,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = block_mvn[["data"]]
    ),
    block_expected,
    tolerance = 1e-10
  )
  expect_error(
    .log_posterior(
      parameters                 = parameters,
      data                       = block_fit_data,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = TRUE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = block_mvn[["data"]]
    ),
    "Known-V bridge likelihoods"
  )
  block_setup <- list(
    fit               = NULL,
    priors            = block_mvn[["priors"]],
    data              = block_mvn[["data"]],
    yi                = block_mvn[["data"]][["outcome"]][["yi"]],
    sei               = block_mvn[["data"]][["outcome"]][["sei"]],
    selection_sei     = block_mvn[["data"]][["outcome"]][["sei"]],
    K                 = 2L,
    S                 = 1L,
    mu                = matrix(c(0, 0), nrow = 1),
    tau_within        = matrix(c(parameters[["tau"]], parameters[["tau"]]),
                               nrow = 1),
    tau_between       = matrix(0, nrow = 1, ncol = 2),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(block_mvn),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )
  expect_equal(
    .log_lik_estimate_sum_from_setup(block_setup),
    block_expected,
    tolerance = 1e-8
  )

  whitened <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  whitened_fit_data <- .create_fit_data(whitened[["data"]], whitened[["priors"]])
  expect_equal(
    .log_posterior(
      parameters                 = parameters,
      data                       = whitened_fit_data,
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = whitened[["data"]]
    ),
    block_expected,
    tolerance = 1e-10
  )

  rank_one_sei        <- c(0.20, 0.30, 0.40)
  rank_one_V          <- tcrossprod(rank_one_sei)
  rank_one_parameters <- list(mu = 0, tau = 0.10, sampling_z = 0)
  rank_one_expected   <- sum(stats::dnorm(
    c(0.10, 0.20, -0.05),
    mean = 0,
    sd   = rank_one_parameters[["tau"]],
    log  = TRUE
  ))

  expect_warning(
    rank_one_block <- brma.mv(
      yi                        = c(0.10, 0.20, -0.05),
      V                         = rank_one_V,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    .log_posterior(
      parameters                 = rank_one_parameters,
      data                       = .create_fit_data(
        rank_one_block[["data"]],
        rank_one_block[["priors"]]
      ),
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = rank_one_block[["data"]]
    ),
    rank_one_expected,
    tolerance = 1e-10
  )

  expect_warning(
    rank_one_whitened <- brma.mv(
      yi                        = c(0.10, 0.20, -0.05),
      V                         = rank_one_V,
      known_v_parameterization  = "whitened",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    .log_posterior(
      parameters                 = rank_one_parameters,
      data                       = .create_fit_data(
        rank_one_whitened[["data"]],
        rank_one_whitened[["priors"]]
      ),
      is_scale                   = FALSE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = rank_one_whitened[["data"]]
    ),
    rank_one_expected,
    tolerance = 1e-10
  )

  dat_scale <- data.frame(
    yi = c(0.10, 0.20, -0.05),
    x  = c(-1, 0, 1)
  )
  V_scale <- matrix(
    c(
      0.05, 0.015, 0.010,
      0.015, 0.08, 0.020,
      0.010, 0.020, 0.07
    ),
    nrow  = 3,
    byrow = TRUE
  )
  scale_tau <- c(0.04, 0.10, 0.18)
  scale_parameters <- list(
    mu      = 0,
    log_tau = log(scale_tau)
  )
  scale_block <- brma.mv(
    yi                        = yi,
    V                         = V_scale,
    scale                     = ~ x,
    data                      = dat_scale,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  scale_expected <- .marglik_mvn_log_density(
    y          = dat_scale[["yi"]],
    mean       = rep(0, 3),
    covariance = V_scale + diag(scale_tau^2, nrow = 3)
  )
  expect_equal(
    .log_posterior(
      parameters                 = scale_parameters,
      data                       = .create_fit_data(scale_block[["data"]], scale_block[["priors"]]),
      is_scale                   = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = scale_block[["data"]]
    ),
    scale_expected,
    tolerance = 1e-10
  )

  dat_random <- data.frame(
    yi    = c(0.10, 0.20),
    study = c("s1", "s2")
  )
  random_sd <- 0.12
  random_parameters <- list(mu = c(0, 0))
  random_parameters[["mu__xREx__study_intercept"]] <- random_sd
  random_expected <- .marglik_mvn_log_density(
    y          = dat_random[["yi"]],
    mean       = c(0, 0),
    covariance = V + diag(random_sd^2, nrow = 2)
  )

  random_block <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat_random,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(
    .log_posterior(
      parameters                 = random_parameters,
      data                       = .create_fit_data(random_block[["data"]], random_block[["priors"]]),
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = random_block[["data"]]
    ),
    random_expected,
    tolerance = 1e-10
  )
  expect_true(.marglik_needs_bridge_context(random_block[["data"]]))

  random_whitened <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat_random,
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  row_parameters <- list(mu = c(0, 0))
  attr(row_parameters, "posterior_samples") <- matrix(
    random_sd,
    nrow     = 1L,
    dimnames = list(NULL, "mu__xREx__study_intercept")
  )
  expect_equal(
    unname(.parameters_as_sample_matrix(row_parameters)[1L, "mu__xREx__study_intercept"]),
    random_sd
  )
  expect_equal(
    .log_posterior(
      parameters                 = row_parameters,
      data                       = .create_fit_data(random_whitened[["data"]], random_whitened[["priors"]]),
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = random_whitened[["data"]]
    ),
    random_expected,
    tolerance = 1e-10
  )

  allocation_dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    study    = c("s1", "s1", "s2", "s2"),
    estimate = paste0("e", 1:4),
    x        = c(0, 1, 0, 1)
  )
  allocation_V <- diag(rep(0.04, 4))
  allocation_tau <- c(0.10, 0.20, 0.30, 0.40)
  allocation_block <- brma.mv(
    yi                        = yi,
    V                         = allocation_V,
    scale                     = ~ x,
    random                    = ~ 1 | study / estimate,
    data                      = allocation_dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  allocation_term <- .data_marginalized_random_effects(
    allocation_block[["data"]]
  )[[1L]]
  allocation_factor <- .marginalized_random_effect_allocation_factors(
    allocation_term
  )[[1L]]
  allocation_weight_name <- allocation_factor[["weight_name"]]
  allocation_bridge_context <- structure(
    list(nodes = c(
      stats::setNames(allocation_tau, paste0("tau[", seq_along(allocation_tau), "]")),
      stats::setNames(
        c(0.25, 0.75),
        paste0(allocation_weight_name, "[", 1:2, "]")
      )
    )),
    class = c("BayesTools_bridge_context", "list")
  )
  allocation_parameters <- list(mu = 0)
  allocation_posterior <- .marglik_bridge_posterior_samples(
    parameters     = allocation_parameters,
    bridge_context = allocation_bridge_context
  )
  expect_equal(
    unname(allocation_posterior[1L, paste0(allocation_weight_name, "[1]")]),
    0.25
  )
  allocation_expected <- .marglik_mvn_log_density(
    y          = allocation_dat[["yi"]],
    mean       = rep(0, 4),
    covariance = allocation_V + diag(allocation_tau^2 * 0.25, nrow = 4)
  )
  expect_true(.marglik_needs_bridge_context(allocation_block[["data"]]))
  expect_equal(
    .log_posterior(
      parameters                 = allocation_parameters,
      data                       = .create_fit_data(allocation_block[["data"]], allocation_block[["priors"]]),
      is_scale                   = TRUE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = allocation_block[["data"]],
      bridge_context             = allocation_bridge_context
    ),
    allocation_expected,
    tolerance = 1e-10
  )

  sd_component_dat <- transform(
    allocation_dat,
    criterion = c("c1", "c2", "c1", "c2")
  )
  sd_component_block <- brma.mv(
    yi                        = yi,
    V                         = allocation_V,
    random                    = list(~ 1 | study / estimate, ~ 1 | criterion),
    data                      = sd_component_dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  sd_component_terms <- .data_marginalized_random_effects(
    sd_component_block[["data"]]
  )
  sd_component_term       <- sd_component_terms[[1L]]
  sd_component_binding    <- sd_component_term[["sd_binding"]]
  sd_component_allocation <- sd_component_binding[["allocations"]][[1L]]
  sd_component_allocation[["target"]]               <- "sd_component"
  sd_component_allocation[["parent_factors"]]       <- sd_component_allocation[["factors"]]
  sd_component_allocation[["leaf_index_by_column"]] <- 1L
  sd_component_allocation[["weight_name"]]          <- "mu__xRE_ALLOCx_leaf__weight"
  sd_component_allocation[["scale"]]                <- "mean_variance"
  sd_component_allocation[["n_targets"]]            <- 3L
  sd_component_binding[["allocations"]][[1L]]       <- sd_component_allocation
  sd_component_term[["sd_parameter_names"]]         <- NA_character_
  sd_component_term[["sd_binding"]]                 <- sd_component_binding
  sd_component_terms[[1L]]      <- sd_component_term
  sd_component_location         <- sd_component_block[["data"]][["location"]]
  attr(sd_component_location, "marginalized_random_effects") <- sd_component_terms
  sd_component_block[["data"]][["location"]] <- sd_component_location

  parent_factors <- sd_component_allocation[["parent_factors"]]
  parent_nodes   <- stats::setNames(
    c(0.25, 0.40),
    vapply(parent_factors, function(factor) {
      paste0(factor[["weight_name"]], "[", factor[["index"]], "]")
    }, character(1))
  )
  sd_component_context <- structure(
    list(nodes = c(
      stats::setNames(
        0.50,
        sd_component_allocation[["source"]][["name"]]
      ),
      parent_nodes,
      stats::setNames(
        0.20,
        paste0(sd_component_allocation[["weight_name"]], "[1]")
      )
    )),
    class = c("BayesTools_bridge_context", "list")
  )
  sd_component_parameters <- list(mu = 0)
  sd_component_posterior  <- .marglik_bridge_posterior_samples(
    parameters     = sd_component_parameters,
    bridge_context = sd_component_context
  )
  sd_component_variance   <- 0.50^2 * 0.25 * 0.40 * (3 * 0.20)
  expect_equal(
    .evaluate_marginalized_random_variance(
      data              = sd_component_block[["data"]],
      posterior_samples = sd_component_posterior
    ),
    matrix(sd_component_variance, nrow = 1L, ncol = 4L)
  )
  sd_component_expected   <- .marglik_mvn_log_density(
    y          = sd_component_dat[["yi"]],
    mean       = rep(0, 4),
    covariance = allocation_V + diag(sd_component_variance, nrow = 4)
  )
  expect_equal(
    .log_posterior(
      parameters                 = sd_component_parameters,
      data                       = .create_fit_data(sd_component_block[["data"]], sd_component_block[["priors"]]),
      is_scale                   = FALSE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = sd_component_block[["data"]],
      bridge_context             = sd_component_context
    ),
    sd_component_expected,
    tolerance = 1e-10
  )
})


test_that("brma.mv preserves tiny nonzero known-V covariance", {

  V <- matrix(c(0.04, 1e-12, 1e-12, 0.09), nrow = 2)

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    known_v_parameterization  = "latent",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_data                 = TRUE
  )
  known_V <- attr(object[["data"]], "known_V_data")

  expect_equal(
    lapply(.known_v_blocks(known_V), `[[`, "index"),
    list(1:2)
  )
  expect_gt(known_V[["rank"]], 0L)
  reconstructed <- diag(known_V[["residual_variance"]], nrow = 2)
  for (block in known_V[["latent_blocks"]]) {
    index <- block[["index"]]
    reconstructed[index, index] <- reconstructed[index, index] +
      tcrossprod(block[["B"]])
  }
  expect_equal(
    reconstructed,
    V,
    tolerance = 1e-10
  )
})


test_that("brma.mv known-V estimate log-likelihood uses Schur conditional target", {

  V <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)
  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  mu <- matrix(
    c(
      0.00, 0.10,
      0.05, 0.15
    ),
    nrow  = 2,
    byrow = TRUE
  )
  tau <- matrix(
    c(
      0.10, 0.20,
      0.05, 0.15
    ),
    nrow  = 2,
    byrow = TRUE
  )
  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = object[["data"]][["outcome"]][["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 2L,
    S                 = 2L,
    mu                = mu,
    tau_within        = tau,
    tau_between       = matrix(0, nrow = 2, ncol = 2),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 2, ncol = 0L)
  )

  log_lik <- .log_lik_estimate_from_setup(setup)
  expected <- matrix(NA_real_, nrow = 2, ncol = 2)
  for (s in seq_len(2)) {
    covariance <- V + diag(tau[s, ]^2, nrow = 2)
    for (i in seq_len(2)) {
      j <- setdiff(seq_len(2), i)
      conditional_mean <- mu[s, i] +
        covariance[i, j] / covariance[j, j] * (setup[["yi"]][j] - mu[s, j])
      conditional_variance <- covariance[i, i] -
        covariance[i, j]^2 / covariance[j, j]
      expected[s, i] <- stats::dnorm(
        setup[["yi"]][i],
        mean = conditional_mean,
        sd   = sqrt(conditional_variance),
        log  = TRUE
      )
    }
  }

  expect_equal(log_lik, expected, tolerance = 1e-12)
  expected_joint <- vapply(seq_len(2), function(s) {
    .marglik_mvn_log_density(
      y          = setup[["yi"]],
      mean       = mu[s, ],
      covariance = V + diag(tau[s, ]^2, nrow = 2)
    )
  }, numeric(1))
  expect_equal(
    .log_lik_estimate_sum_from_setup(setup),
    expected_joint,
    tolerance = 1e-12
  )

  setup_negative <- setup
  setup_negative[["effect_direction"]] <- "negative"
  log_lik_negative <- .log_lik_estimate_from_setup(setup_negative)
  expected_negative <- matrix(NA_real_, nrow = 2, ncol = 2)
  for (s in seq_len(2)) {
    covariance <- V + diag(tau[s, ]^2, nrow = 2)
    for (i in seq_len(2)) {
      j <- setdiff(seq_len(2), i)
      conditional_mean <- -mu[s, i] +
        covariance[i, j] / covariance[j, j] * (-setup[["yi"]][j] + mu[s, j])
      conditional_variance <- covariance[i, i] -
        covariance[i, j]^2 / covariance[j, j]
      expected_negative[s, i] <- stats::dnorm(
        -setup[["yi"]][i],
        mean = conditional_mean,
        sd   = sqrt(conditional_variance),
        log  = TRUE
      )
    }
  }
  expect_equal(log_lik_negative, expected_negative, tolerance = 1e-12)

  target <- .estimate_log_lik_target_metadata(
    setup             = setup,
    data_hash         = .get_outcome_hash(object),
    dependency_blocks = .known_v_dependency_blocks(
      object[["data"]],
      setup[["K"]]
    )
  )
  expect_true(target[["known_v_estimate_backend"]])
  expect_true(target[["dependency_conditioning"]])
  expect_equal(target[["target"]], "estimate_log_score")
})


test_that("known-V dependency blocks use canonical block metadata", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20, -0.05),
    V                         = matrix(
      c(
        0.04, 0.01, 0.00,
        0.01, 0.09, 0.00,
        0.00, 0.00, 0.16
      ),
      nrow = 3,
      byrow = TRUE
    ),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  blocks <- .known_v_dependency_blocks(object[["data"]], K = 3L)
  expect_equal(blocks, list(1:2, 3L))

  data_bad_blocks <- object[["data"]]
  known_V         <- .data_known_v_data(data_bad_blocks)
  known_V[["blocks"]] <- list(list(
    index      = c(1L, 1L),
    covariance = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2L)
  ))
  attr(data_bad_blocks, "known_V_data") <- known_V
  expect_error(
    .known_v_dependency_blocks(data_bad_blocks, K = 3L),
    "partition"
  )
})


test_that("singular PSD known-V Cholesky targets fail with targeted messages", {

  sei <- c(0.20, 0.30, 0.40)
  V   <- tcrossprod(sei)

  expect_warning(
    object <- brma.mv(
      yi                        = c(0.10, 0.20, -0.05),
      V                         = V,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )

  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = object[["data"]][["outcome"]][["yi"]],
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = matrix(c(0, 0, 0), nrow = 1),
    tau_within        = matrix(0, nrow = 1, ncol = 3),
    tau_between       = matrix(0, nrow = 1, ncol = 3),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )

  expect_error(
    .log_lik_estimate_from_setup(setup),
    "not positive definite"
  )
})


test_that("brma.mv singular-V preflight requires structural regularization", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    study    = rep("s1", 3),
    estimate = paste0("e", 1:3),
    x        = c(-1, 0, 1)
  )
  sei <- c(0.20, 0.30, 0.40)
  V   <- tcrossprod(sei)
  A         <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, ncol = 2)
  V_general <- tcrossprod(A)

  prior_zero <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 0)
  )
  prior_positive <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 0.10)
  )
  prior_tiny_positive <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 1e-6)
  )

  expect_false(.known_v_nullspace_is_regularized(V, rep(FALSE, 3)))
  expect_false(.known_v_nullspace_is_regularized(V, c(TRUE, FALSE, FALSE)))
  expect_true(.known_v_nullspace_is_regularized(V, c(TRUE, TRUE, FALSE)))

  expect_warning(
    general_positive <- brma.mv(
      yi                        = yi,
      V                         = V_general,
      data                      = dat,
      prior_heterogeneity       = prior_positive,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "rank-deficient correlation structure"
  )
  expect_s3_class(general_positive, "brma.mv")
  expect_warning(
    scale_positive <- brma.mv(
      yi                        = yi,
      V                         = V_general,
      scale                     = ~ x,
      data                      = dat,
      prior_scale               = list(
        intercept = prior_positive,
        x         = BayesTools::prior(
          distribution = "spike",
          parameters   = list(location = log(2))
        )
      ),
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "rank-deficient correlation structure"
  )
  expect_equal(
    .brma_mv_fixed_integrated_variance(scale_positive, K = 3L),
    (0.10 * exp(log(2) * dat[["x"]]))^2,
    tolerance = 1e-14
  )
  V_tolerance <- matrix(c(1, 1 + 1e-9, 1 + 1e-9, 1), nrow = 2)
  expect_error(
    brma.mv(
      yi                        = yi,
      V                         = V_tolerance,
      data                      = dat[1:2, , drop = FALSE],
      prior_heterogeneity       = prior_tiny_positive,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  V_tight <- matrix(
    c(1, 1 + .Machine$double.eps, 1 + .Machine$double.eps, 1),
    nrow = 2
  )
  expect_error(
    .known_v_as_matrix(V_tight),
    "positive semidefinite"
  )
  V_rank_one <- tcrossprod(c(1, 1))
  prior_too_small <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 1e-9)
  )
  prior_scale_too_small <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = 1e-9)
  )
  prior_scale_invalid <- BayesTools::prior(
    distribution = "spike",
    parameters   = list(location = log(1e-9))
  )
  expect_warning(
    tiny_positive <- brma.mv(
      yi                        = yi,
      V                         = V_rank_one,
      data                      = dat[1:2, , drop = FALSE],
      prior_heterogeneity       = prior_too_small,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(
    .data_known_v_effective_backend(tiny_positive[["data"]]),
    "latent"
  )
  expect_error(
    suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V_rank_one,
      scale                     = ~ 1,
      data                      = dat[1:2, , drop = FALSE],
      prior_scale               = list(intercept = prior_scale_invalid),
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )),
    "strictly positive"
  )
  expect_warning(
    tiny_scale <- brma.mv(
      yi                        = yi,
      V                         = V_rank_one,
      scale                     = ~ 1,
      data                      = dat[1:2, , drop = FALSE],
      prior_scale               = list(intercept = prior_scale_too_small),
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(.data_known_v_effective_backend(tiny_scale[["data"]]), "latent")
  expect_error(
    suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V_general,
      data                      = dat,
      prior_heterogeneity       = prior_zero,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )),
    "not structurally regularized"
  )

  for (parameterization in c("whitened", "block_mvn")) {
    expect_error(
      suppressWarnings(brma.mv(
        yi                        = yi,
        V                         = V,
        data                      = dat,
        prior_heterogeneity       = prior_zero,
        known_v_parameterization  = parameterization,
        measure                   = "GEN",
        prior_unit_information_sd = 1,
        only_priors               = TRUE
      )),
      "not structurally regularized"
    )
    expect_warning(
      positive <- brma.mv(
        yi                        = yi,
        V                         = V,
        data                      = dat,
        prior_heterogeneity       = prior_positive,
        known_v_parameterization  = parameterization,
        measure                   = "GEN",
        prior_unit_information_sd = 1,
        only_priors               = TRUE
      ),
      "rank-deficient correlation structure"
    )
    expect_s3_class(positive, "brma.mv")
  }

  expect_error(
    suppressWarnings(brma.mv(
      yi                         = yi,
      V                          = V,
      random                     = ~ 1 | study,
      data                       = dat,
      known_v_parameterization   = "block_mvn",
      measure                    = "GEN",
      prior_unit_information_sd  = 1,
      only_priors                = TRUE
    )),
    "Sampled random effects change the conditional mean"
  )

  expect_warning(
    marginalized <- brma.mv(
      yi                         = yi,
      V                          = V,
      random                     = ~ 1 | estimate,
      data                       = dat,
      known_v_parameterization   = "block_mvn",
      measure                    = "GEN",
      prior_unit_information_sd  = 1,
      only_priors                = TRUE
    ),
    "rank-deficient correlation structure"
  )
  expect_true(.data_has_marginalized_random_effects(marginalized[["data"]]))

  expect_warning(
    allocated <- brma.mv(
      yi                         = yi,
      V                          = V,
      scale                      = ~ x,
      random                     = ~ 1 | study / estimate,
      data                       = dat,
      known_v_parameterization   = "block_mvn",
      measure                    = "GEN",
      prior_unit_information_sd  = 1,
      only_priors                = TRUE
    ),
    "rank-deficient correlation structure"
  )
  expect_true(.data_has_marginalized_random_effects(allocated[["data"]]))

  expect_error(
    suppressWarnings(brma.mv(
      yi                        = yi,
      V                         = V,
      random                    = ~ 1 | estimate,
      data                      = dat,
      prior_heterogeneity       = BayesTools::prior_random(
        estimate = BayesTools::random_block(sd = prior_zero)
      ),
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    )),
    "not structurally regularized"
  )
  expect_warning(
    tiny_random <- brma.mv(
      yi                        = yi,
      V                         = V,
      random                    = ~ 1 | estimate,
      data                      = dat,
      prior_heterogeneity       = BayesTools::prior_random(
        estimate = BayesTools::random_block(sd = prior_too_small)
      ),
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_unit_information_sd = 1,
      only_priors               = TRUE
    ),
    "positive semidefinite"
  )
  expect_equal(.data_known_v_effective_backend(tiny_random[["data"]]), "latent")
})


test_that("brma.mv Schur estimate target matches partitioned MVN conditionals", {

  V <- matrix(c(
    0.10, 0.03, 0.02,
    0.03, 0.12, 0.04,
    0.02, 0.04, 0.11
  ), nrow = 3)
  yi <- c(0.10, -0.05, 0.20)
  mu <- matrix(c(0.00, -0.02, 0.10), nrow = 1)
  tau <- matrix(c(0.05, 0.10, 0.08), nrow = 1)

  object <- brma.mv(
    yi                        = yi,
    V                         = V,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  setup <- list(
    fit               = NULL,
    priors            = object[["priors"]],
    data              = object[["data"]],
    yi                = yi,
    sei               = object[["data"]][["outcome"]][["sei"]],
    selection_sei     = object[["data"]][["outcome"]][["sei"]],
    K                 = 3L,
    S                 = 1L,
    mu                = mu,
    tau_within        = tau,
    tau_between       = matrix(0, nrow = 1, ncol = 3),
    cluster           = NULL,
    weights           = NULL,
    data_hash         = .get_outcome_hash(object),
    is_weightfunction = FALSE,
    outcome_type      = "norm",
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1, ncol = 0L)
  )

  covariance <- V + diag(as.numeric(tau)^2, nrow = 3)
  expected   <- matrix(NA_real_, nrow = 1, ncol = 3)
  for (i in seq_len(3)) {
    j <- setdiff(seq_len(3), i)
    conditional_mean <- mu[1, i] +
      covariance[i, j, drop = FALSE] %*%
      solve(covariance[j, j, drop = FALSE], yi[j] - mu[1, j])
    conditional_variance <- covariance[i, i] -
      covariance[i, j, drop = FALSE] %*%
      solve(covariance[j, j, drop = FALSE],
            covariance[j, i, drop = FALSE])
    expected[1, i] <- stats::dnorm(
      yi[i],
      mean = conditional_mean,
      sd   = sqrt(conditional_variance),
      log  = TRUE
    )
  }

  expect_true(any(lengths(.known_v_dependency_blocks(
    object[["data"]], K = 3L
  )) > 1L))
  expect_equal(.log_lik_estimate_from_setup(setup), expected,
               tolerance = 1e-12)
})


test_that("brma.mv known-V target hash includes covariance structure", {

  object_1 <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.00, 0.00, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  object_2 <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_false(identical(
    .get_outcome_hash(object_1),
    .get_outcome_hash(object_2)
  ))
})


test_that("brma.mv known-V target hash normalizes signed zero", {

  V_positive_zero <- matrix(
    c(
      0.04, 0.01, 0,
      0.01, 0.09, 0.02,
      0, 0.02, 0.16
    ),
    nrow = 3L,
    byrow = TRUE
  )
  V_negative_zero       <- V_positive_zero
  V_negative_zero[1, 3] <- -0
  V_negative_zero[3, 1] <- -0
  object_1 <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V_positive_zero,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  object_2 <- brma.mv(
    yi                        = c(0.10, 0.20, 0.30),
    V                         = V_negative_zero,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_identical(.get_outcome_hash(object_1), .get_outcome_hash(object_2))
})


test_that("brma.mv allows estimate log-likelihood diagnostics for correlated V", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    .check_log_lik_target_available(object, "cluster", "test()"),
    "unit = 'cluster'"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
})


test_that("brma.mv allows estimate log-likelihood diagnostics for correlated whitened V", {

  object <- brma.mv(
    yi                        = c(0.10, 0.20),
    V                         = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
    known_v_parameterization  = "whitened",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_silent(
    .check_log_lik_target_available(object, "estimate", "test()")
  )
  expect_error(
    .check_log_lik_target_available(object, "cluster", "test()"),
    "unit = 'cluster'"
  )
  expect_silent(
    .check_log_lik_target_available(object, "estimate", "add_waic()")
  )
})


test_that("brma.mv regplot PI/SI no longer use the old feature guard", {

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    mods                      = ~ x,
    data                      = data.frame(
      yi = c(0.10, 0.20, 0.30, 0.40),
      x  = c(0, 1, 0, 1)
    ),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_error(
    regplot(object, pi = TRUE, as_data = TRUE),
    "Posterior samples are required"
  )
  expect_error(
    regplot(object, si = TRUE, as_data = TRUE),
    "Posterior samples are required"
  )
})
