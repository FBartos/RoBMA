context("Exact selection certified factor routing")
skip_on_cran()

test_that("known_v_factor preserves exact provenance and structural routing", {

  K <- 6L
  dat <- data.frame(
    yi    = seq(-.3, .7, length.out = K),
    study = factor(rep(c("a", "b"), each = 3L))
  )
  loading <- cbind(
    c(.12, .08, .04, .10, .06, .03),
    c(-.03, .07, .11, .04, .09, .13)
  )
  V_factor <- known_v_factor(rep(.02, K), loading)
  covariance <- diag(.02, K) + tcrossprod(loading)

  factor_object <- bselmodel.mv(
    yi                        = yi,
    V                         = V_factor,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  factor_setup <- .data_exact_selection_setup(factor_object[["data"]])
  factor_syntax <- .create_model_syntax(
    factor_object[["data"]], factor_object[["priors"]]
  )

  expect_identical(
    .known_v_storage(.data_known_v_data(factor_object[["data"]])),
    "factor"
  )
  expect_equal(
    .selection_exact_sampling_block(factor_setup[["sampling"]], seq_len(K)),
    covariance
  )
  expect_identical(factor_setup[["schema_version"]], 1L)
  expect_identical(factor_setup[["exactness"]], "EF")
  expect_identical(
    factor_setup[["factor_ranks"]],
    2L
  )
  expect_match(factor_syntax, "dselnorm_factor_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", factor_syntax, fixed = TRUE))

  dense_object <- bselmodel.mv(
    yi                        = yi,
    V                         = covariance,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  dense_setup <- .data_exact_selection_setup(dense_object[["data"]])
  dense_syntax <- .create_model_syntax(
    dense_object[["data"]], dense_object[["priors"]]
  )
  expect_identical(dense_setup[["exactness"]], "E2")
  expect_match(dense_syntax, "dselnorm_mnorm_step", fixed = TRUE)
  expect_false(grepl("dselnorm_factor_step", dense_syntax, fixed = TRUE))
})


test_that("factor routing fails closed outside the certified kernel contract", {

  K <- 6L
  dat <- data.frame(yi = seq(-.3, .7, length.out = K))
  high_rank <- matrix(seq_len(K * 5L) / 100, nrow = K, ncol = 5L)
  zero_residual <- matrix(c(
    .2, .1,
    .1, .2,
    .2, .1,
    .1, .2,
    .2, .1,
    .1, .2
  ), nrow = K, byrow = TRUE)

  make_object <- function(V) {
    bselmodel.mv(
      yi                        = yi,
      V                         = V,
      data                      = dat,
      measure                   = "SMD",
      prior_unit_information_sd = 1,
      selection_likelihood      = "exact",
      only_priors               = TRUE,
      silent                    = TRUE
    )
  }
  high_rank_object <- make_object(known_v_factor(rep(.02, K), high_rank))
  expect_warning(
    zero_residual_object <- make_object(known_v_factor(
      rep(0, K), zero_residual
    )),
    "positive semidefinite"
  )

  for (object in list(high_rank_object, zero_residual_object)) {
    setup <- .data_exact_selection_setup(object[["data"]])
    syntax <- .create_model_syntax(object[["data"]], object[["priors"]])
    expect_identical(setup[["exactness"]], "E2")
    expect_match(syntax, "dselnorm_mnorm_step", fixed = TRUE)
    expect_false(grepl("dselnorm_factor_step", syntax, fixed = TRUE))
  }

  random_dat <- data.frame(
    yi    = dat[["yi"]],
    study = factor(rep("a", K))
  )
  random_high_rank_object <- bselmodel.mv(
    yi                        = yi,
    V                         = known_v_factor(rep(.02, K), high_rank),
    random                    = ~ 1 | study,
    data                      = random_dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  random_setup    <- .data_exact_selection_setup(
    random_high_rank_object[["data"]]
  )
  random_fit_data <- .create_fit_data(
    random_high_rank_object[["data"]],
    random_high_rank_object[["priors"]]
  )
  random_syntax   <- .create_model_syntax(
    random_high_rank_object[["data"]],
    random_high_rank_object[["priors"]]
  )
  expect_identical(
    random_setup[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_identical(
    random_setup[["block_methods"]],
    "dense"
  )
  expect_true("sel_exact_block_1_diagonal" %in% names(random_fit_data))
  expect_match(random_syntax, "sel_exact_block_1_diagonal", fixed = TRUE)

  expect_error(known_v_factor(c(.1, -.1), matrix(0, 2L, 1L)),
               "non-negative")
  expect_error(known_v_factor(c(.1, .1), matrix(0, 3L, 1L)),
               "one row")
  expect_error(known_v_factor(c(.1, .1), matrix(c(0, NA), 2L, 1L)),
               "finite")
})


test_that("higher-rank factor controls retain budget semantics", {

  default_control <- set_selection_likelihood_control()
  odd_control <- set_selection_likelihood_control(
    points_per_scramble     = 9L,
    max_points_per_scramble = 9L,
    scrambles               = 2L
  )
  sampling <- structure(list(
    representation = "diagonal_factor",
    diagonal       = rep(1, 4L),
    loading        = matrix(0, nrow = 4L, ncol = 0L)
  ), class = c("RoBMA_selection_sampling_plan", "list"))
  factor_plan <- .selection_exact_execution_plan(
    row_blocks       = list(1:4),
    block_methods     = "factor",
    factor_ranks      = 3L,
    selection_control = odd_control,
    sampling          = sampling,
    sampling_factor_blocks = NULL,
    random_covariance = NULL
  )
  dense_plan <- .selection_exact_execution_plan(
    row_blocks       = list(1:4),
    block_methods     = "dense",
    factor_ranks      = NA_integer_,
    selection_control = odd_control,
    sampling          = sampling,
    sampling_factor_blocks = NULL,
    random_covariance = NULL
  )

  expect_identical(default_control[["points_per_scramble"]], 512L)
  expect_identical(default_control[["max_points_per_scramble"]], 8192L)
  expect_identical(default_control[["scrambles"]], 8L)
  expect_identical(default_control[["relative_tolerance"]], .005)
  expect_identical(default_control[["seed"]], 1L)
  expect_identical(factor_plan[["points_per_scramble"]], 9L)
  expect_identical(factor_plan[["max_points_per_scramble"]], 9L)
  expect_identical(factor_plan[["factor_points_per_proposal"]], 5L)
  expect_identical(factor_plan[["factor_max_points_per_proposal"]], 5L)
  expect_identical(
    dim(factor_plan[["designs"]][["factor_3"]]),
    c(2L, 5L, 6L)
  )
  expect_identical(
    dim(dense_plan[["designs"]][["4"]]),
    c(2L, 9L, 8L)
  )
})


test_that("post-fit likelihood retains certified sampling factors", {

  K <- 4L
  loading <- cbind(
    c(.12, .08, .04, .10),
    c(-.03, .07, .11, .04)
  )
  object <- bselmodel.mv(
    yi                        = seq(-.3, .6, length.out = K),
    V                         = known_v_factor(rep(.02, K), loading),
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  factor_call <- NULL
  testthat::local_mocked_bindings(
    .estimate_normal_covariance_target_location_from_setup = function(setup) {
      list(
        y     = seq(-.3, .6, length.out = K),
        means = matrix(0, nrow = 2L, ncol = K)
      )
    },
    .selection_exact_signed_context = function(setup, signed_yi) {
      list(obs_bin = rep(1L, K))
    },
    .selection_exact_random_covariance_samples = function(setup) NULL,
    .selection_exact_covariance_lower = function(...) {
      stop("Certified post-fit factors must not use dense covariance.")
    },
    .selection_exact_factor_loglik_block = function(
        yi, means, residual_sd, loading, sei, selection_context,
        execution_plan, block_index) {
      factor_call <<- list(
        residual_sd = residual_sd,
        loading     = loading,
        block_index = block_index
      )
      c(11, 12)
    },
    .package = "RoBMA"
  )

  observed <- .selection_exact_block_loglik_from_setup(list(
    data          = object[["data"]],
    S             = 2L,
    selection_sei = sqrt(.known_v_diagonal(
      .data_known_v_data(object[["data"]])
    )),
    tau_within    = matrix(.15, nrow = 2L, ncol = K),
    tau_between   = matrix(0, nrow = 2L, ncol = K),
    is_multilevel = FALSE
  ))

  expect_equal(observed, matrix(c(11, 12), nrow = 2L))
  expect_identical(factor_call[["block_index"]], 1L)
  expect_equal(
    factor_call[["residual_sd"]],
    matrix(sqrt(.02 + .15^2), nrow = 2L, ncol = K),
    tolerance = 0
  )
  expect_equal(
    matrix(factor_call[["loading"]], nrow = 2L),
    matrix(array(
      rep(loading, each = 2L),
      dim = c(2L, K, ncol(loading))
    ), nrow = 2L),
    tolerance = 0
  )
})


test_that("bridge factor states retain the certified covariance exactly", {

  dat <- data.frame(
    yi    = c(-.20, .05, .35),
    vi    = c(.012, .018, .025),
    study = factor(rep("a", 3L)),
    time  = 1:3
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ har(time | study),
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact_setup <- .data_exact_selection_setup(object[["data"]])
  random_term <- object[["formula_design"]][["mu"]][["random_effects"]][[1L]]
  coefficient_factor <- matrix(
    c(.25, 0, 0, .08, .22, 0, -.03, .07, .18),
    nrow = 3L,
    byrow = TRUE
  )
  bridge_factor <- list(
    representation = "factor_state",
    row_blocks     = exact_setup[["row_blocks"]],
    factor_plans   = list(study = list(
      type                  = "group",
      model_matrix          = random_term[["model_matrix"]],
      group_map             = random_term[["group_map"]],
      coefficient_structure = "markov"
    )),
    factor_states  = list(study = list(
      coefficient_factor = coefficient_factor
    ))
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  bridge_data <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = object[["priors"]],
    effect_direction = .data_effect_direction(object[["data"]]),
    model_data       = object[["data"]]
  )
  factor_call <- NULL
  testthat::local_mocked_bindings(
    .marglik_bridge_random_covariance = function(...) bridge_factor,
    .marglik_random_covariance_dense = function(...) {
      stop("Certified bridge factors must not be materialized as dense.")
    },
    .marglik_selection_context = function(parameters, data) {
      list(obs_bin = rep(1L, 3L))
    },
    .selection_exact_factor_loglik_block = function(
        yi, means, residual_sd, loading, sei, selection_context,
        execution_plan, block_index) {
      factor_call <<- list(
        residual_sd = residual_sd,
        loading     = loading,
        block_index = block_index
      )
      17
    },
    .package = "RoBMA"
  )

  observed <- .marglik_exact_selection_log_lik(
    parameters            = list(),
    data                  = bridge_data,
    model_data            = object[["data"]],
    bridge_context        = list(),
    covariance_plan_cache = NULL,
    mu_samples            = matrix(0, nrow = 1L, ncol = 3L),
    tau_within_samples    = matrix(0, nrow = 1L, ncol = 3L),
    tau_between_samples   = matrix(0, nrow = 1L, ncol = 3L),
    is_random             = TRUE,
    is_multilevel         = FALSE,
    fixed_zero_random     = FALSE,
    K                     = 3L
  )

  rank <- exact_setup[["factor_ranks"]][[1L]]
  loading <- matrix(factor_call[["loading"]], nrow = 3L, ncol = rank)
  reconstructed <- diag(as.numeric(factor_call[["residual_sd"]])^2) +
    tcrossprod(loading)
  basis <- random_term[["model_matrix"]] %*% coefficient_factor
  expected <- diag(dat$vi) + tcrossprod(basis)

  expect_identical(observed, 17)
  expect_identical(factor_call[["block_index"]], 1L)
  expect_identical(rank, 2L)
  expect_equal(unname(reconstructed), unname(expected), tolerance = 1e-14)
})


test_that("bridge dense fallback reconstructs the exact required block", {

  dat <- data.frame(
    yi    = c(-.20, .05, .35),
    study = factor(rep("a", 3L))
  )
  sampling <- diag(c(.02, .03, .04)) + tcrossprod(c(.03, .02, .01))
  object <- bselmodel.mv(
    yi                        = yi,
    V                         = sampling,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    selection_likelihood      = "exact",
    only_priors               = TRUE,
    silent                    = TRUE
  )
  exact_setup <- .data_exact_selection_setup(object[["data"]])
  random_term <- object[["formula_design"]][["mu"]][["random_effects"]][[1L]]
  bridge_factor <- list(
    representation = "factor_state",
    row_blocks     = exact_setup[["row_blocks"]],
    factor_plans   = list(study = list(
      type                  = "group",
      model_matrix          = random_term[["model_matrix"]],
      group_map             = random_term[["group_map"]],
      coefficient_structure = "diagonal"
    )),
    factor_states  = list(study = list(
      coefficient_factor = matrix(.2, 1L, 1L)
    ))
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  bridge_data <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = object[["priors"]],
    effect_direction = .data_effect_direction(object[["data"]]),
    model_data       = object[["data"]]
  )
  observed_lower <- NULL
  testthat::local_mocked_bindings(
    .marglik_bridge_random_covariance = function(...) bridge_factor,
    .marglik_selection_context = function(parameters, data) {
      list(obs_bin = rep(1L, 3L))
    },
    .selection_exact_joint_loglik_block = function(
        yi, means, covariance_lower, sei, selection_context,
        execution_plan, block_size) {
      observed_lower <<- covariance_lower
      23
    },
    .package = "RoBMA"
  )

  observed <- .marglik_exact_selection_log_lik(
    parameters            = list(),
    data                  = bridge_data,
    model_data            = object[["data"]],
    bridge_context        = list(),
    covariance_plan_cache = NULL,
    mu_samples            = matrix(0, nrow = 1L, ncol = 3L),
    tau_within_samples    = matrix(0, nrow = 1L, ncol = 3L),
    tau_between_samples   = matrix(0, nrow = 1L, ncol = 3L),
    is_random             = TRUE,
    is_multilevel         = FALSE,
    fixed_zero_random     = FALSE,
    K                     = 3L
  )

  expected <- sampling + matrix(.04, 3L, 3L)
  pairs <- .selection_exact_lower_pairs(
    exact_setup,
    1:3
  )
  expected_lower <- matrix(
    expected[cbind(pairs[["row_1"]], pairs[["row_2"]])],
    nrow = 1L
  )
  expect_identical(exact_setup[["exactness"]], "E2")
  expect_identical(observed, 23)
  expect_equal(observed_lower, expected_lower, tolerance = 1e-14)
})


test_that("native factor likelihood agrees with an independent MVN oracle", {

  reference_loglik <- function(y, mu, covariance, sei, omega, z_lower,
                               z_upper, obs_bin, sign) {
    assignments <- expand.grid(rep(list(seq_along(omega)), length(y)))
    terms <- apply(assignments, 1L, function(bins) {
      lower <- if (sign == 1L) {
        z_lower[bins] * sei
      } else {
        -z_upper[bins] * sei
      }
      upper <- if (sign == 1L) {
        z_upper[bins] * sei
      } else {
        -z_lower[bins] * sei
      }
      probability <- suppressWarnings(as.numeric(mvtnorm::pmvnorm(
        lower     = lower,
        upper     = upper,
        mean      = mu,
        sigma     = covariance,
        algorithm = mvtnorm::Miwa(steps = 4096L)
      )))
      prod(omega[bins]) * probability
    })
    mvtnorm::dmvnorm(y, mu, covariance, log = TRUE) +
      sum(log(omega[obs_bin])) - log(sum(terms))
  }

  set.seed(927)
  errors <- diagnostics <- normal_errors <- numeric(3L)
  for (rank in 2:4) {
    K <- 4L
    residual_sd <- seq(.28, .48, length.out = K)
    loading <- matrix(stats::rnorm(K * rank, sd = .22), K, rank)
    loading[, 1L] <- loading[, 1L] + seq(.15, .35, length.out = K)
    mu  <- seq(-.35, .45, length.out = K)
    y   <- c(-.42, -.05, .35, .91)
    sei <- seq(.20, .36, length.out = K)
    sign     <- if (rank == 3L) -1L else 1L
    omega    <- c(.18, .45, 1)
    z_lower  <- stats::qnorm(c(.025, .10, 1), lower.tail = FALSE)
    z_upper  <- stats::qnorm(c(0, .025, .10), lower.tail = FALSE)
    score    <- sign * y / sei
    obs_bin  <- vapply(score, function(value) {
      which(value >= z_lower)[[1L]]
    }, integer(1L))
    qmc <- BayesTools::selection_qmc_design(
      dimensions = 2L * rank,
      points     = 4096L,
      scrambles  = 16L,
      seed       = 417L
    )
    quadrature <- .selection_exact_cluster_quadrature_rules(
      SELNORM_FACTOR_QUADRATURE_ORDERS[[as.character(rank)]]
    )
    observed <- .Call(
      "RoBMA_selnorm_factor_step_loglik_batch",
      as.double(y),
      matrix(mu, nrow = 1L),
      matrix(residual_sd, nrow = 1L),
      matrix(as.double(loading), nrow = 1L),
      as.double(sei),
      matrix(omega, nrow = 1L),
      as.double(z_lower),
      as.double(z_upper),
      as.integer(obs_bin),
      sign,
      TRUE,
      1L,
      quadrature[["nodes"]],
      quadrature[["log_weights"]],
      as.double(quadrature[["orders"]]),
      as.double(qmc),
      4096L,
      4096L,
      16L,
      .005,
      PACKAGE = "RoBMA"
    )
    covariance <- diag(residual_sd^2) + tcrossprod(loading)
    normal_observed <- .Call(
      "RoBMA_selnorm_factor_step_loglik_batch",
      as.double(y), matrix(mu, nrow = 1L),
      matrix(residual_sd, nrow = 1L),
      matrix(as.double(loading), nrow = 1L),
      as.double(sei), matrix(omega, nrow = 1L),
      as.double(z_lower), as.double(z_upper), as.integer(obs_bin),
      sign, TRUE, 0L, quadrature[["nodes"]],
      quadrature[["log_weights"]], as.double(quadrature[["orders"]]),
      as.double(qmc), 4096L, 4096L, 16L, .005,
      PACKAGE = "RoBMA"
    )
    reference <- reference_loglik(
      y, mu, covariance, sei, omega, z_lower, z_upper, obs_bin, sign
    )
    errors[[rank - 1L]] <- observed[["log_density"]] - reference
    normal_errors[[rank - 1L]] <- normal_observed[["log_density"]] -
      mvtnorm::dmvnorm(y, mu, covariance, log = TRUE)
    diagnostics[[rank - 1L]] <- max(
      observed[["relative_mcse"]], observed[["relative_change"]]
    )
  }

  expect_true(all(abs(errors) < .005))
  expect_equal(normal_errors, rep(0, 3L), tolerance = 1e-11)
  expect_true(all(diagnostics < .005))

  K <- 4L
  residual_sd <- c(.28, .32, .36, .40)
  loading <- cbind(
    c(.25, .20, .15, .10),
    c(.12, .08, .10, .05),
    c(0, .09, .07, .06),
    c(0, 0, 0, .08)
  )
  mu      <- c(-.20, -.05, .15, .35)
  y       <- c(-.31, .02, .29, .61)
  sei     <- c(.20, .24, .28, .32)
  omega   <- c(.35, 1)
  z_lower <- c(stats::qnorm(.025, lower.tail = FALSE), -Inf)
  z_upper <- c(Inf, stats::qnorm(.025, lower.tail = FALSE))
  obs_bin <- ifelse(y / sei >= z_lower[[1L]], 1L, 2L)
  quadrature <- .selection_exact_cluster_quadrature_rules(
    SELNORM_FACTOR_QUADRATURE_ORDERS[["4"]]
  )
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 8L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 419L
  )
  observed <- .Call(
    "RoBMA_selnorm_factor_step_loglik_batch",
    y, matrix(mu, nrow = 1L), matrix(residual_sd, nrow = 1L),
    matrix(as.double(loading), nrow = 1L), sei, matrix(omega, nrow = 1L),
    z_lower, z_upper, as.integer(obs_bin), 1L, TRUE, 1L,
    quadrature[["nodes"]], quadrature[["log_weights"]],
    as.double(quadrature[["orders"]]), as.double(qmc),
    8L, 8L, 2L, .005,
    PACKAGE = "RoBMA"
  )
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  reference <- reference_loglik(
    y, mu, covariance, sei, omega, z_lower, z_upper, obs_bin, 1L
  )

  expect_equal(observed[["relative_mcse"]], 0, tolerance = 0)
  expect_lte(observed[["relative_change"]], .005)
  expect_lt(abs(observed[["log_density"]] - reference), .005)
})


test_that("factor quadrature resolves an Assink state before QMC fallback", {

  y <- c(
    2.2844, 2.1771, 1.7777, 1.5480, 1.4855,
    1.4836, 1.2777, 1.0311, .9409, .6263
  )
  vi <- c(
    .3325, .3073, .2697, .4533, .1167,
    .1706, .1538, .3132, .1487, .2139
  )
  mean <- rep(-.462179, length(y))
  residual_sd <- sqrt(.3 * vi + .315042^2)
  loading <- cbind(sqrt(.7 * vi), rep(.245830, length(y)))
  sei <- sqrt(vi)
  omega <- c(1, .562858)
  z_lower <- c(stats::qnorm(.025, lower.tail = FALSE), -Inf)
  z_upper <- c(Inf, stats::qnorm(.025, lower.tail = FALSE))
  obs_bin <- ifelse(y / sei >= z_lower[[1L]], 1L, 2L)
  quadrature <- .selection_exact_cluster_quadrature_rules(
    SELNORM_FACTOR_QUADRATURE_ORDERS[["2"]]
  )
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 4L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 1L
  )
  observed <- .Call(
    "RoBMA_selnorm_factor_step_loglik_batch",
    y, matrix(mean, nrow = 1L), matrix(residual_sd, nrow = 1L),
    matrix(as.double(loading), nrow = 1L), sei, matrix(omega, nrow = 1L),
    z_lower, z_upper, as.integer(obs_bin), 1L, TRUE, 1L,
    quadrature[["nodes"]], quadrature[["log_weights"]],
    as.double(quadrature[["orders"]]), as.double(qmc), 8L, 8L, 2L, .005,
    PACKAGE = "RoBMA"
  )
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  normal_log_density <- mvtnorm::dmvnorm(
    y,
    mean       = mean,
    sigma      = covariance,
    log        = TRUE
  )
  observed_log_normalizer <- normal_log_density +
    sum(log(omega[obs_bin])) - observed[["log_density"]]

  expect_equal(observed[["relative_mcse"]], 0, tolerance = 0)
  expect_lte(observed[["relative_change"]], .005)
  expect_equal(
    observed_log_normalizer,
    -5.5163303949665297,
    tolerance = 5e-5
  )
})


test_that("JAGS instantiates the certified factor distribution", {

  skip_if_not_installed("rjags")
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 4L,
    points     = 8L,
    scrambles  = 2L,
    seed       = 19L
  )
  quadrature <- .selection_exact_cluster_quadrature_rules(
    SELNORM_FACTOR_QUADRATURE_ORDERS[["2"]]
  )
  model_text <- paste0(
    "model{\n",
    "  y[1:3] ~ dselnorm_factor_step(",
    "mu[1:3],residual_sd[1:3],loading[1:3,1:2],",
    "sei[1:3],omega,z_lower,z_upper,obs_bin[1:3],",
    "1,1,1,nodes,log_weights,orders,",
    "qmc[1:2,1:8,1:4],8,8,2,0.005)\n",
    "}\n"
  )
  connection <- textConnection(model_text)
  on.exit(close(connection), add = TRUE)

  model <- rjags::jags.model(
    file = connection,
    data = list(
      y           = c(-.2, .1, .4),
      mu          = c(0, .05, .1),
      residual_sd = c(.2, .25, .3),
      loading     = matrix(c(.1, .05, .03, -.04, .02, .08), 3L, 2L),
      sei         = c(.15, .2, .25),
      omega       = 1,
      z_lower     = -1e300,
      z_upper     = 1e300,
      obs_bin     = rep(1L, 3L),
      nodes       = quadrature[["nodes"]],
      log_weights = quadrature[["log_weights"]],
      orders      = quadrature[["orders"]],
      qmc         = qmc
    ),
    n.chains = 1L,
    n.adapt  = 0L,
    quiet    = TRUE
  )

  expect_s3_class(model, "jags")
})
